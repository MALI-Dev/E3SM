! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_velocity_gnn
!
!> \brief MPAS land-ice GNN velocity emulator
!> \author Trevor Hillebrand
!> \date   June 2026
!> \details
!>  This module couples a PyTorch Geometric GNN velocity emulator
!>  to MALI via FTorch. The GNN predicts depth-averaged x/y velocities
!>  from cell-based features and graph connectivity.
!>  Data is gathered to rank 0 for global inference and scattered back.
!
!-----------------------------------------------------------------------

module li_velocity_gnn

#ifdef USE_FTORCH
   use ftorch, only: torch_model, torch_tensor, torch_kCPU, torch_kCUDA, &
                     torch_tensor_from_array, torch_model_load, &
                     torch_model_forward, torch_delete
#endif

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use mpi

   use li_setup
   use li_mask
   use li_constants
   use ieee_arithmetic
   use ieee_exceptions

   implicit none
   private

   public :: li_velocity_gnn_init, &
             li_velocity_gnn_solve, &
             li_velocity_gnn_finalize

   !--------------------------------------------------------------------
   ! Private module variables (persistent state)
   !--------------------------------------------------------------------
#ifdef USE_FTORCH
   type(torch_model), save :: gnn_model
#endif

   ! Global undirected edge pairs (0-based cell indices, built once in init)
   integer, dimension(:,:), allocatable, save :: global_edge_pairs
   integer, save :: nUndirectedEdgesGlobal  ! number of undirected edges

   ! MPI gather/scatter bookkeeping
   integer, save :: nCellsGlobal          ! total number of cells across all ranks
   integer, save :: nLocalOwned           ! nCellsSolve on this rank
   integer, allocatable, save :: recvcounts(:)  ! number of owned cells per rank
   integer, allocatable, save :: displs(:)      ! displacements for gatherv/scatterv

   ! Device type for FTorch
   integer, save :: ftorch_device_type

   ! Number of input features. Features 1-12 are assembled per-rank; features
   ! 13-14 (surface-gradient dsurf/dx, dsurf/dy) are computed on rank 0 after
   ! the gather, because the least-squares gradient needs each cell's neighbors
   ! and rank 0 holds the full global mesh (features + global_edge_pairs).
   integer, parameter :: nFeatures = 14

   logical, save :: gnn_initialized = .false.

contains

!***********************************************************************
!
!  routine li_velocity_gnn_init
!
!> \brief   Initialize GNN velocity solver
!> \author  Trevor Hillebrand
!> \date    June 2026
!> \details
!>  Load TorchScript model, build global edge_index, set up MPI
!>  gather/scatter infrastructure.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_gnn_init(domain, err)

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      ! local variables
      character (len=StrKIND), pointer :: config_gnn_model_path
      character (len=StrKIND), pointer :: config_gnn_device_type
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      integer, pointer :: nCellsSolve, nEdgesSolve
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:), pointer :: indexToCellID, indexToEdgeID
      integer :: myRank, nProcs, iEdge, idx
      integer :: nEdgesLocal, nValidEdges
      ! Temporary arrays for gathering edge connectivity
      integer, allocatable :: local_edges(:,:)  ! (2, nValidEdges) with global cell IDs
      integer, allocatable :: edgecounts(:), edgedispls(:)
      integer, allocatable :: all_edges_flat(:)
      integer :: nEdgesAll, i

      err = 0

#ifdef USE_FTORCH
      call mpas_pool_get_config(liConfigs, 'config_gnn_model_path', config_gnn_model_path)
      call mpas_pool_get_config(liConfigs, 'config_gnn_device_type', config_gnn_device_type)

      ! Set device type
      if (trim(config_gnn_device_type) == 'cuda') then
         ftorch_device_type = torch_kCUDA
      else
         ftorch_device_type = torch_kCPU
      endif

      myRank = domain % dminfo % my_proc_id
      nProcs = domain % dminfo % nProcs

      ! --- Determine local owned cells and set up gather/scatter ---
      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      nLocalOwned = nCellsSolve
      nEdgesLocal = nEdgesSolve

      ! Gather nCellsSolve from each rank to build displacements
      allocate(recvcounts(nProcs))
      allocate(displs(nProcs))
      call MPI_Allgather(nLocalOwned, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, &
                         domain % dminfo % comm, err)
      nCellsGlobal = sum(recvcounts)
      displs(1) = 0
      do i = 2, nProcs
         displs(i) = displs(i-1) + recvcounts(i-1)
      end do

      ! --- Build global edge_index ---
      ! Collect local edges as global cell ID pairs, skipping boundary
      ! edges where cellsOnEdge contains 0 (garbage cells)
      nValidEdges = 0
      do iEdge = 1, nEdgesLocal
         if (cellsOnEdge(1, iEdge) > 0 .and. cellsOnEdge(2, iEdge) > 0) then
            nValidEdges = nValidEdges + 1
         endif
      end do

      allocate(local_edges(2, nValidEdges))
      idx = 0
      do iEdge = 1, nEdgesLocal
         if (cellsOnEdge(1, iEdge) > 0 .and. cellsOnEdge(2, iEdge) > 0) then
            idx = idx + 1
            local_edges(1, idx) = indexToCellID(cellsOnEdge(1, iEdge))
            local_edges(2, idx) = indexToCellID(cellsOnEdge(2, iEdge))
         endif
      end do

      ! Gather edge counts from all ranks
      allocate(edgecounts(nProcs))
      allocate(edgedispls(nProcs))
      call MPI_Allgather(nValidEdges * 2, 1, MPI_INTEGER, edgecounts, 1, MPI_INTEGER, &
                         domain % dminfo % comm, err)
      nEdgesAll = sum(edgecounts) / 2  ! total undirected edges
      edgedispls(1) = 0
      do i = 2, nProcs
         edgedispls(i) = edgedispls(i-1) + edgecounts(i-1)
      end do

      ! Gather all edges to rank 0
      if (myRank == 0) then
         allocate(all_edges_flat(sum(edgecounts)))
      else
         allocate(all_edges_flat(1))  ! placeholder
      endif

      call MPI_Gatherv(local_edges, nValidEdges * 2, MPI_INTEGER, &
                       all_edges_flat, edgecounts, edgedispls, MPI_INTEGER, &
                       0, domain % dminfo % comm, err)

      ! On rank 0, store undirected edge pairs with 0-based indexing
      if (myRank == 0) then
         nUndirectedEdgesGlobal = nEdgesAll
         allocate(global_edge_pairs(2, nUndirectedEdgesGlobal))
         do i = 1, nUndirectedEdgesGlobal
            global_edge_pairs(1, i) = all_edges_flat((i-1)*2 + 1) - 1  ! 0-based
            global_edge_pairs(2, i) = all_edges_flat((i-1)*2 + 2) - 1  ! 0-based
         end do
      else
         nUndirectedEdgesGlobal = 0
         allocate(global_edge_pairs(2, 1))  ! placeholder
      endif

      deallocate(local_edges, edgecounts, edgedispls, all_edges_flat)

      ! Broadcast nUndirectedEdgesGlobal to all ranks (needed for consistency)
      call MPI_Bcast(nUndirectedEdgesGlobal, 1, MPI_INTEGER, 0, domain % dminfo % comm, err)

      ! --- Load the TorchScript model on rank 0 ---
      if (myRank == 0) then
         call torch_model_load(gnn_model, trim(config_gnn_model_path), ftorch_device_type)
         call mpas_log_write('GNN velocity model loaded from: ' // trim(config_gnn_model_path))
      endif

      gnn_initialized = .true.
      call mpas_log_write('GNN velocity solver initialized.')
      call mpas_log_write('  Global cells: $i, Undirected edges: $i', &
           intArgs=(/nCellsGlobal, nUndirectedEdgesGlobal/))
      if (myRank == 0) then
         call mpas_log_write('nCellsGlobal: $i', intArgs=(/nCellsGlobal/))
         call mpas_log_write('edge_pairs min: $i', intArgs=(/minval(global_edge_pairs)/))
         call mpas_log_write('edge_pairs max: $i', intArgs=(/maxval(global_edge_pairs)/))
      endif

#else
      call mpas_log_write('GNN velocity solver requires compilation with USE_FTORCH.', MPAS_LOG_CRIT)
      err = 1
#endif

   end subroutine li_velocity_gnn_init


!***********************************************************************
!
!  routine li_velocity_gnn_solve
!
!> \brief   Solve velocity using GNN emulator
!> \author  Trevor Hillebrand
!> \date    June 2026
!> \details
!>  Gathers input features to rank 0, runs GNN inference, scatters
!>  predicted velocities back, and fills uReconstructX/Y and
!>  normalVelocity fields.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_gnn_solve(domain, err)

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: err

      ! local variables
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool, geometryPool, thermalPool, velocityPool
      integer, pointer :: nCellsSolve, nEdges, nVertInterfaces
      integer, dimension(:), pointer :: cellMask, indexToCellID
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell, thickness, bedTopography
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell
      real (kind=RKIND), dimension(:), pointer :: muFriction, stiffnessFactor
      real (kind=RKIND), dimension(:,:), pointer :: temperature
      real (kind=RKIND), dimension(:), pointer :: basalTemperature
      real (kind=RKIND), dimension(:), pointer :: layerThicknessFractions
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY
      real (kind=RKIND), dimension(:,:), pointer :: edgeNormalVectors
      real (kind=RKIND), pointer :: config_sea_level, config_ice_density, config_ocean_density

      ! Local feature assembly
      real (kind=RKIND), allocatable :: local_features(:,:)  ! (nLocalOwned, nFeatures)
      real (kind=RKIND), allocatable :: local_velocity(:,:)  ! (nLocalOwned, 2)

      ! Global arrays on rank 0 (indexed by global cell ID)
      real (kind=RKIND), allocatable :: global_features(:,:)   ! (nCellsGlobal, nFeatures)
      real (kind=RKIND), allocatable :: global_velocity(:,:)   ! (nCellsGlobal, 2)
      real (kind=RKIND), allocatable :: gathered_features(:,:) ! (nCellsGlobal, nFeatures) in gather order
      real (kind=RKIND), allocatable :: gathered_velocity(:,:) ! (nCellsGlobal, 2) in gather order
      integer, allocatable :: global_cellIDs(:)                ! (nCellsGlobal) maps gathered pos -> global cell ID

      ! Ice-only subgraph arrays (rank 0 only)
      integer, allocatable :: ice_map(:)                                ! maps global cell ID to compact ice index
      integer :: nIceCells, nIceEdgesUndirected, nIceEdgesDirected
      real (kind=4), allocatable, target :: ice_features(:,:)           ! (nIceCells, nFeatures) single precision
      real (kind=4), allocatable, target :: ice_velocity(:,:)           ! (nIceCells, 2) single precision
      integer(kind=8), allocatable, target :: ice_edge_index(:,:)       ! (2, nIceEdgesDirected)
      integer(kind=8), allocatable, target :: ice_batch(:)              ! (nIceCells)

#ifdef USE_FTORCH
      type(torch_tensor), dimension(3) :: in_tensors
      type(torch_tensor), dimension(1) :: out_tensors
#endif

      integer :: iCell, iEdge, i, idx, k, myRank, cell1, cell2, c1, c2
      integer :: nVertInterfacesLocal
      real (kind=RKIND) :: xVelEdge, yVelEdge
      real (kind=RKIND) :: thicknessAboveFloatation, depthAvgTemp
      integer, pointer :: nVertLevels
      logical :: fpe_invalid, fpe_divzero, fpe_overflow  ! saved FPE trap state

      err = 0

#ifdef USE_FTORCH
      if (.not. gnn_initialized) then
         call mpas_log_write('GNN velocity solver not initialized.', MPAS_LOG_ERR)
         err = 1
         return
      endif

      myRank = domain % dminfo % my_proc_id

      ! --- Get pools from block ---
      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
      call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      nVertInterfacesLocal = nVertInterfaces

      ! Get mesh arrays
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormalVectors)
      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      ! Get geometry arrays
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      ! Get thermal arrays
      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'basalTemperature', basalTemperature)

      ! Get velocity arrays
      call mpas_pool_get_array(velocityPool, 'muFriction', muFriction)
      call mpas_pool_get_array(velocityPool, 'stiffnessFactor', stiffnessFactor)
      call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
      call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)

      ! Get configs
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)

      ! --- Assemble local node features ---
      allocate(local_features(nLocalOwned, nFeatures))
      do iCell = 1, nLocalOwned
         ! Feature 1: xCell (km)
         local_features(iCell, 1) = xCell(iCell) / 1000.0_RKIND

         ! Feature 2: yCell (km)
         local_features(iCell, 2) = yCell(iCell) / 1000.0_RKIND

         ! Feature 3: areaCell (km^2)
         local_features(iCell, 3) = areaCell(iCell) / 1.0e6_RKIND

         ! Feature 4: thickness (km)
         local_features(iCell, 4) = thickness(iCell) / 1000.0_RKIND

         ! Feature 5: bedTopography (km)
         local_features(iCell, 5) = bedTopography(iCell) / 1000.0_RKIND

         ! Feature 6: log10(muFriction)
         if (muFriction(iCell) > 0.0_RKIND) then
            local_features(iCell, 6) = log10(muFriction(iCell))
         else
            local_features(iCell, 6) = 0.0_RKIND
         endif

         ! Feature 7: basalTemperature (°C)
         local_features(iCell, 7) = basalTemperature(iCell) - kelvin_to_celsius

         ! Feature 8: depth-averaged temperature (°C)
         depthAvgTemp = 0.0_RKIND
         do k = 1, nVertLevels
            depthAvgTemp = depthAvgTemp + temperature(k, iCell) * layerThicknessFractions(k)
         end do
         local_features(iCell, 8) = depthAvgTemp - kelvin_to_celsius

         ! Feature 9: stiffnessFactor
         local_features(iCell, 9) = stiffnessFactor(iCell)

         ! Feature 10: thicknessAboveFloatation (km)
         thicknessAboveFloatation = thickness(iCell) - &
            max(0.0_RKIND, (config_sea_level - bedTopography(iCell)) * &
                config_ocean_density / config_ice_density)
         local_features(iCell, 10) = thicknessAboveFloatation / 1000.0_RKIND

         ! Feature 11: floatMask (1.0 if floating, 0.0 otherwise)
         if (li_mask_is_floating_ice(cellMask(iCell))) then
            local_features(iCell, 11) = 1.0_RKIND
         else
            local_features(iCell, 11) = 0.0_RKIND
         endif

         ! Feature 12: margin mask (1.0 if ice margin, 0.0 otherwise)
         if (li_mask_is_margin(cellMask(iCell))) then
            local_features(iCell, 12) = 1.0_RKIND
         else
            local_features(iCell, 12) = 0.0_RKIND
         endif

         ! Features 13-14 (surface-gradient dsurf/dx, dsurf/dy) cannot be
         ! computed per-rank without halo neighbor values; they are filled on
         ! rank 0 after the gather (see compute_surface_gradient_rank0). Zero
         ! them here so the gather never propagates uninitialized values.
         local_features(iCell, 13) = 0.0_RKIND
         local_features(iCell, 14) = 0.0_RKIND
      end do

      ! --- MPI Gather features and cell IDs to rank 0 ---
      if (myRank == 0) then
         allocate(gathered_features(nCellsGlobal, nFeatures))
         allocate(global_features(nCellsGlobal, nFeatures))
         allocate(global_velocity(nCellsGlobal, 2))
         allocate(gathered_velocity(nCellsGlobal, 2))
         allocate(global_cellIDs(nCellsGlobal))
         global_velocity(:,:) = 0.0_RKIND  ! zero velocity for non-ice cells
      else
         allocate(gathered_features(1, 1))  ! placeholder
         allocate(global_features(1, 1))  ! placeholder
         allocate(global_velocity(1, 1))  ! placeholder
         allocate(gathered_velocity(1, 1))  ! placeholder
         allocate(global_cellIDs(1))  ! placeholder
      endif

      ! Gather cell IDs so we can reorder from gathered position to global cell ID order
      call MPI_Gatherv(indexToCellID, nLocalOwned, MPI_INTEGER, &
                       global_cellIDs, recvcounts, displs, MPI_INTEGER, &
                       0, domain % dminfo % comm, err)

      ! Gather each feature column separately to handle Fortran column-major layout
      call gather_features_to_rank0(local_features, gathered_features, &
                                    nLocalOwned, nCellsGlobal, nFeatures, &
                                    recvcounts, displs, domain % dminfo % comm, myRank)

      ! Reorder features from gathered order to global cell ID order
      ! so that features and edge indices use the same index space
      if (myRank == 0) then
         do i = 1, nCellsGlobal
            global_features(global_cellIDs(i), :) = gathered_features(i, :)
         end do

         ! Fill features 13-14 (surface-gradient) now that the full global mesh
         ! is assembled on rank 0. Uses global_features columns and the global
         ! edge list, so no halo exchange is needed and it stays consistent with
         ! the other (already-gathered) features.
         call compute_surface_gradient_rank0(global_features, nCellsGlobal, &
                                             nFeatures)
      endif

      ! --- Filter to ice-only cells and run inference on rank 0 ---
      if (myRank == 0) then
         ! Identify ice cells (thickness > 0, feature 4)
         allocate(ice_map(nCellsGlobal))
         nIceCells = 0
         do iCell = 1, nCellsGlobal
            if (global_features(iCell, 4) > 0.0_RKIND) then  ! thickness in km
               nIceCells = nIceCells + 1
               ice_map(iCell) = nIceCells - 1  ! 0-based compact index
            else
               ice_map(iCell) = -1
            endif
         end do

         ! Build ice-only features
         allocate(ice_features(nIceCells, nFeatures))
         idx = 0
         do iCell = 1, nCellsGlobal
            if (ice_map(iCell) >= 0) then
               idx = idx + 1
               ice_features(idx, :) = global_features(iCell, :)
            endif
         end do

         ! Debug: print first ice cell features from global (double) and ice (single)
         call mpas_log_write('First global_features (double): x=$r y=$r thick=$r', &
              realArgs=(/global_features(1, 1), global_features(1, 2), global_features(1, 4)/))
         call mpas_log_write('First ice_features (single): x=$r y=$r thick=$r', &
              realArgs=(/dble(ice_features(1, 1)), dble(ice_features(1, 2)), dble(ice_features(1, 4))/))
         call mpas_log_write('nIceCells=$i, nIceEdgesDirected will be computed next', &
              intArgs=(/nIceCells/))

         ! Build ice-only bidirectional edge_index
         ! First count edges where both endpoints have ice
         nIceEdgesUndirected = 0
         do i = 1, nUndirectedEdgesGlobal
            c1 = global_edge_pairs(1, i) + 1  ! 1-based for ice_map lookup
            c2 = global_edge_pairs(2, i) + 1
            if (c1 >= 1 .and. c1 <= nCellsGlobal .and. &
                c2 >= 1 .and. c2 <= nCellsGlobal) then
               if (ice_map(c1) >= 0 .and. ice_map(c2) >= 0) then
                  nIceEdgesUndirected = nIceEdgesUndirected + 1
               endif
            endif
         end do

         nIceEdgesDirected = nIceEdgesUndirected * 2
         allocate(ice_edge_index(2, nIceEdgesDirected))
         idx = 0
         do i = 1, nUndirectedEdgesGlobal
            c1 = global_edge_pairs(1, i) + 1
            c2 = global_edge_pairs(2, i) + 1
            if (c1 >= 1 .and. c1 <= nCellsGlobal .and. &
                c2 >= 1 .and. c2 <= nCellsGlobal) then
               if (ice_map(c1) >= 0 .and. ice_map(c2) >= 0) then
                  idx = idx + 1
                  ice_edge_index(1, idx) = ice_map(c1)
                  ice_edge_index(2, idx) = ice_map(c2)
                  idx = idx + 1
                  ice_edge_index(1, idx) = ice_map(c2)
                  ice_edge_index(2, idx) = ice_map(c1)
               endif
            endif
         end do

         ! Build batch tensor
         allocate(ice_batch(nIceCells))
         ice_batch(:) = 0

         ! Build tensors and run inference
         allocate(ice_velocity(nIceCells, 2))

         ! Disable FPE traps around LibTorch call (LibTorch triggers them internally)
         call ieee_get_halting_mode(ieee_invalid, fpe_invalid)
         call ieee_get_halting_mode(ieee_divide_by_zero, fpe_divzero)
         call ieee_get_halting_mode(ieee_overflow, fpe_overflow)
         call ieee_set_halting_mode(ieee_invalid, .false.)
         call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
         call ieee_set_halting_mode(ieee_overflow, .false.)

         call torch_tensor_from_array(in_tensors(1), ice_features, [1, 2], ftorch_device_type)
         call torch_tensor_from_array(in_tensors(2), ice_edge_index, [1, 2], ftorch_device_type)
         call torch_tensor_from_array(in_tensors(3), ice_batch, [1], ftorch_device_type)
         call torch_tensor_from_array(out_tensors(1), ice_velocity, [1, 2], torch_kCPU)

         call torch_model_forward(gnn_model, in_tensors, out_tensors)

         call torch_delete(in_tensors)
         call torch_delete(out_tensors)

         ! Restore FPE traps
         call ieee_set_halting_mode(ieee_invalid, fpe_invalid)
         call ieee_set_halting_mode(ieee_divide_by_zero, fpe_divzero)
         call ieee_set_halting_mode(ieee_overflow, fpe_overflow)

         ! Map ice velocities back to global array
         do iCell = 1, nCellsGlobal
            if (ice_map(iCell) >= 0) then
               global_velocity(iCell, 1) = ice_velocity(ice_map(iCell) + 1, 1)
               global_velocity(iCell, 2) = ice_velocity(ice_map(iCell) + 1, 2)
            endif
         end do

         call mpas_log_write('GNN inference: $i ice cells, $i directed ice edges', &
              intArgs=(/nIceCells, nIceEdgesDirected/))
         call mpas_log_write('GNN velocity stats (km/yr): vx min=$r max=$r, vy min=$r max=$r', &
              realArgs=(/dble(minval(ice_velocity(:,1))), dble(maxval(ice_velocity(:,1))), &
                         dble(minval(ice_velocity(:,2))), dble(maxval(ice_velocity(:,2)))/))
         call mpas_log_write('GNN velocity cell 1: vx=$r vy=$r (km/yr)', &
              realArgs=(/dble(ice_velocity(1,1)), dble(ice_velocity(1,2))/))

         ! Dump GNN inputs/outputs as text for debugging (first solve call only)
         call dump_gnn_debug(ice_features, ice_edge_index, ice_velocity, &
                             nIceCells, nFeatures, nIceEdgesDirected)

         deallocate(ice_map, ice_features, ice_velocity, ice_edge_index, ice_batch)
      endif

      ! Reorder velocity from global cell ID order back to gathered order for scattering
      if (myRank == 0) then
         do i = 1, nCellsGlobal
            gathered_velocity(i, :) = global_velocity(global_cellIDs(i), :)
         end do
      endif

      ! --- Scatter predicted velocities back to all ranks ---
      allocate(local_velocity(nLocalOwned, 2))
      call scatter_velocity_from_rank0(gathered_velocity, local_velocity, &
                                       nLocalOwned, nCellsGlobal, &
                                       recvcounts, displs, domain % dminfo % comm, myRank)

      ! --- Fill uReconstructX/Y with uniform vertical profile ---
      ! GNN predicts in km/yr; convert to m/s
      do iCell = 1, nLocalOwned
         do k = 1, nVertInterfacesLocal
            uReconstructX(k, iCell) = local_velocity(iCell, 1) * 1000.0_RKIND / scyr
            uReconstructY(k, iCell) = local_velocity(iCell, 2) * 1000.0_RKIND / scyr
         end do
      end do

      ! Halo exchange uReconstructX/Y so that halo cells have current values
      ! before computing normalVelocity on owned edges with halo-cell neighbors
      call mpas_dmpar_field_halo_exch(domain, 'uReconstructX')
      call mpas_dmpar_field_halo_exch(domain, 'uReconstructY')

      ! --- Compute normalVelocity from cell-centered velocities ---
      ! At margin edges (ice cell <-> non-ice cell), use one-sided velocity
      ! from the ice cell rather than averaging with zero. This allows
      ! advection to transport ice into empty cells, analogous to how Albany
      ! provides velocity on its boundary edges.
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         do k = 1, nVertInterfacesLocal
            if (thickness(cell1) > 0.0_RKIND .and. thickness(cell2) <= 0.0_RKIND) then
               ! Margin edge: only cell1 has ice, use its velocity
               xVelEdge = uReconstructX(k, cell1)
               yVelEdge = uReconstructY(k, cell1)
            else if (thickness(cell2) > 0.0_RKIND .and. thickness(cell1) <= 0.0_RKIND) then
               ! Margin edge: only cell2 has ice, use its velocity
               xVelEdge = uReconstructX(k, cell2)
               yVelEdge = uReconstructY(k, cell2)
            else
               ! Interior edge or both cells empty: average
               xVelEdge = 0.5_RKIND * (uReconstructX(k, cell1) + uReconstructX(k, cell2))
               yVelEdge = 0.5_RKIND * (uReconstructY(k, cell1) + uReconstructY(k, cell2))
            endif
            normalVelocity(k, iEdge) = xVelEdge * edgeNormalVectors(1, iEdge) + &
                                       yVelEdge * edgeNormalVectors(2, iEdge)
         end do
      end do

      ! Clean up
      deallocate(local_features, local_velocity, global_features, global_velocity)
      deallocate(gathered_features, gathered_velocity, global_cellIDs)

#else
      call mpas_log_write('GNN velocity solver requires compilation with USE_FTORCH.', MPAS_LOG_CRIT)
      err = 1
#endif

   end subroutine li_velocity_gnn_solve


!***********************************************************************
!
!  routine li_velocity_gnn_finalize
!
!> \brief   Finalize GNN velocity solver
!> \author  Trevor Hillebrand
!> \date    June 2026
!> \details
!>  Clean up FTorch model and deallocate persistent arrays.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_gnn_finalize(err)

      integer, intent(out) :: err

      err = 0

#ifdef USE_FTORCH
      if (gnn_initialized) then
         call torch_delete(gnn_model)
         if (allocated(global_edge_pairs)) deallocate(global_edge_pairs)
         if (allocated(recvcounts)) deallocate(recvcounts)
         if (allocated(displs)) deallocate(displs)
         gnn_initialized = .false.
      endif
#endif

   end subroutine li_velocity_gnn_finalize


!***********************************************************************
! Private helper routines
!***********************************************************************

!-----------------------------------------------------------------------
!
!  routine gather_features_to_rank0
!
!> \brief   Gather node features from all ranks to rank 0
!
!-----------------------------------------------------------------------

   subroutine gather_features_to_rank0(local_features, global_features, &
                                       nLocal, nGlobal, nFeat, &
                                       counts, disps, comm, myRank)

      real (kind=RKIND), intent(in) :: local_features(:,:)
      real (kind=RKIND), intent(inout) :: global_features(:,:)
      integer, intent(in) :: nLocal, nGlobal, nFeat
      integer, intent(in) :: counts(:), disps(:)
      integer, intent(in) :: comm, myRank

      real (kind=RKIND), allocatable :: local_col(:), global_col(:)
      integer :: iFeat, ierr

      allocate(local_col(nLocal))
      if (myRank == 0) then
         allocate(global_col(nGlobal))
      else
         allocate(global_col(1))
      endif

      do iFeat = 1, nFeat
         local_col(1:nLocal) = local_features(1:nLocal, iFeat)
         call MPI_Gatherv(local_col, nLocal, MPI_DOUBLE_PRECISION, &
                          global_col, counts, disps, MPI_DOUBLE_PRECISION, &
                          0, comm, ierr)
         if (myRank == 0) then
            global_features(1:nGlobal, iFeat) = global_col(1:nGlobal)
         endif
      end do

      deallocate(local_col, global_col)

   end subroutine gather_features_to_rank0


!-----------------------------------------------------------------------
!
!  routine compute_surface_gradient_rank0
!
!> \brief   Fill features 13-14 with the surface gradient (dsurf/dx, dsurf/dy)
!> \author  Trevor Hillebrand
!> \details
!>  Rank 0 only. Computes a per-cell least-squares surface gradient over ice
!>  neighbors, matching utils/load_data.py:surface_gradient in the training
!>  pipeline. MUST stay consistent with that Python routine: same surface
!>  definition (bed+thickness grounded; thickness*(1-rhoi/rhosw) floating),
!>  same neighbor set (mesh cell adjacency), same closed-form 2x2 normal
!>  equations, same degenerate handling (<2 neighbors or near-singular -> 0),
!>  and all arithmetic in METERS (Python works in meters; the threshold 1e-12
!>  is unit-dependent). global_features stores thickness/bed in km, so they are
!>  scaled to meters here. Densities are hard-coded to match the Python
!>  training code (rhoi=910, rhosw=1028), NOT the run configs, so that the
!>  feature the model sees is identical to what it was trained on.
!
!-----------------------------------------------------------------------

   subroutine compute_surface_gradient_rank0(global_features, nGlobal, nFeat)

      real (kind=RKIND), intent(inout) :: global_features(:,:)
      integer, intent(in) :: nGlobal, nFeat

      ! Must match the Python training constants (utils/load_data.py), not the
      ! run configs, so inputs are identical to those seen during training.
      real (kind=RKIND), parameter :: rhoi_train = 910.0_RKIND
      real (kind=RKIND), parameter :: rhosw_train = 1028.0_RKIND

      real (kind=RKIND), allocatable :: surface(:), xm(:), ym(:)
      logical, allocatable :: is_ice(:)
      ! Per-cell accumulators for the 2x2 normal equations.
      real (kind=RKIND), allocatable :: sxx(:), sxy(:), syy(:), bx(:), by(:)
      integer, allocatable :: nnbr(:)
      integer :: iCell, iEdge, a, b
      real (kind=RKIND) :: dx, dy, db, det

      allocate(surface(nGlobal), xm(nGlobal), ym(nGlobal), is_ice(nGlobal))
      allocate(sxx(nGlobal), sxy(nGlobal), syy(nGlobal), bx(nGlobal), by(nGlobal))
      allocate(nnbr(nGlobal))

      ! Reconstruct surface (m), coordinates (m), and ice mask from the gathered
      ! feature columns: 1=x[km], 2=y[km], 4=thickness[km], 5=bed[km],
      ! 11=floatMask. (These columns are already consistent with the model's
      ! other inputs, so surface is consistent too.)
      do iCell = 1, nGlobal
         xm(iCell) = global_features(iCell, 1) * 1000.0_RKIND
         ym(iCell) = global_features(iCell, 2) * 1000.0_RKIND
         is_ice(iCell) = (global_features(iCell, 4) > 0.0_RKIND)
         if (global_features(iCell, 11) > 0.5_RKIND) then
            ! Floating: flotation surface.
            surface(iCell) = global_features(iCell, 4) * 1000.0_RKIND * &
                             (1.0_RKIND - rhoi_train / rhosw_train)
         else
            ! Grounded: bed + thickness.
            surface(iCell) = (global_features(iCell, 5) + &
                              global_features(iCell, 4)) * 1000.0_RKIND
         endif
      end do

      sxx(:) = 0.0_RKIND; sxy(:) = 0.0_RKIND; syy(:) = 0.0_RKIND
      bx(:) = 0.0_RKIND;  by(:) = 0.0_RKIND;  nnbr(:) = 0

      ! Accumulate normal-equation sums by iterating the undirected edge list.
      ! For each ice-ice edge (a,b), cell a gains neighbor b and vice versa --
      ! the same neighbor set Python builds from cellsOnCell (both encode cell
      ! adjacency across shared edges). Non-ice neighbors are excluded.
      do iEdge = 1, nUndirectedEdgesGlobal
         a = global_edge_pairs(1, iEdge) + 1  ! 0-based -> 1-based
         b = global_edge_pairs(2, iEdge) + 1
         if (a < 1 .or. a > nGlobal .or. b < 1 .or. b > nGlobal) cycle
         if (.not. is_ice(a) .or. .not. is_ice(b)) cycle

         ! Contribution to cell a (neighbor b).
         dx = xm(b) - xm(a); dy = ym(b) - ym(a); db = surface(b) - surface(a)
         sxx(a) = sxx(a) + dx*dx; sxy(a) = sxy(a) + dx*dy; syy(a) = syy(a) + dy*dy
         bx(a) = bx(a) + dx*db;   by(a) = by(a) + dy*db;   nnbr(a) = nnbr(a) + 1

         ! Contribution to cell b (neighbor a).
         dx = xm(a) - xm(b); dy = ym(a) - ym(b); db = surface(a) - surface(b)
         sxx(b) = sxx(b) + dx*dx; sxy(b) = sxy(b) + dx*dy; syy(b) = syy(b) + dy*dy
         bx(b) = bx(b) + dx*db;   by(b) = by(b) + dy*db;   nnbr(b) = nnbr(b) + 1
      end do

      ! Solve the 2x2 system per cell; degenerate cells keep gradient 0.
      do iCell = 1, nGlobal
         global_features(iCell, 13) = 0.0_RKIND
         global_features(iCell, 14) = 0.0_RKIND
         if (nnbr(iCell) < 2) cycle
         det = sxx(iCell) * syy(iCell) - sxy(iCell) * sxy(iCell)
         if (abs(det) < 1.0e-12_RKIND) cycle
         global_features(iCell, 13) = &
            (syy(iCell) * bx(iCell) - sxy(iCell) * by(iCell)) / det
         global_features(iCell, 14) = &
            (sxx(iCell) * by(iCell) - sxy(iCell) * bx(iCell)) / det
      end do

      deallocate(surface, xm, ym, is_ice, sxx, sxy, syy, bx, by, nnbr)

   end subroutine compute_surface_gradient_rank0


!-----------------------------------------------------------------------
!
!  routine scatter_velocity_from_rank0
!
!> \brief   Scatter predicted velocities from rank 0 to all ranks
!
!-----------------------------------------------------------------------

   subroutine scatter_velocity_from_rank0(global_velocity, local_velocity, &
                                          nLocal, nGlobal, &
                                          counts, disps, comm, myRank)

      real (kind=RKIND), intent(in) :: global_velocity(:,:)
      real (kind=RKIND), intent(out) :: local_velocity(:,:)
      integer, intent(in) :: nLocal, nGlobal
      integer, intent(in) :: counts(:), disps(:)
      integer, intent(in) :: comm, myRank

      real (kind=RKIND), allocatable :: local_col(:), global_col(:)
      integer :: iComp, ierr

      allocate(local_col(nLocal))
      if (myRank == 0) then
         allocate(global_col(nGlobal))
      else
         allocate(global_col(1))
      endif

      ! Scatter x-velocity (column 1) and y-velocity (column 2)
      do iComp = 1, 2
         if (myRank == 0) then
            global_col(1:nGlobal) = global_velocity(1:nGlobal, iComp)
         endif
         call MPI_Scatterv(global_col, counts, disps, MPI_DOUBLE_PRECISION, &
                           local_col, nLocal, MPI_DOUBLE_PRECISION, &
                           0, comm, ierr)
         local_velocity(1:nLocal, iComp) = local_col(1:nLocal)
      end do

      deallocate(local_col, global_col)

   end subroutine scatter_velocity_from_rank0


!-----------------------------------------------------------------------
!
!  routine dump_gnn_debug
!
!> \brief   Write GNN inputs/outputs to binary files for debugging
!> \details Only writes once (first call). Files can be read in Python:
!>    features = np.fromfile('gnn_debug_features.bin', dtype=np.float32).reshape(nIceCells, nFeatures)
!>    edge_index = np.fromfile('gnn_debug_edge_index.bin', dtype=np.int64).reshape(2, nEdges)
!>    velocity = np.fromfile('gnn_debug_velocity.bin', dtype=np.float32).reshape(nIceCells, 2)
!
!-----------------------------------------------------------------------

   subroutine dump_gnn_debug(features, edge_index, velocity, nCells, nFeat, nEdges)

      real (kind=4), intent(in) :: features(:,:)
      integer(kind=8), intent(in) :: edge_index(:,:)
      real (kind=4), intent(in) :: velocity(:,:)
      integer, intent(in) :: nCells, nFeat, nEdges

      logical, save :: dumped = .false.
      integer :: iu, ios
      integer :: i

      if (dumped) return
      dumped = .true.

      call mpas_log_write('Dumping GNN debug text files: nCells=$i, nFeat=$i, nEdges=$i', &
           intArgs=(/nCells, nFeat, nEdges/))

      ! Write all features as text (formatted I/O known to work)
      open(newunit=iu, file='gnn_debug_features.txt', status='replace', iostat=ios)
      write(iu, '(A,I8,A,I3)') '# nCells=', nCells, ' nFeat=', nFeat
      do i = 1, nCells
         write(iu, '(*(ES15.7,1X))') features(i, 1:nFeat)
      end do
      close(iu)

      ! Write all edge pairs as text
      open(newunit=iu, file='gnn_debug_edges.txt', status='replace', iostat=ios)
      write(iu, '(A,I10)') '# nEdges=', nEdges
      do i = 1, nEdges
         write(iu, '(2(I10,1X))') edge_index(1, i), edge_index(2, i)
      end do
      close(iu)

      ! Write all velocities as text
      open(newunit=iu, file='gnn_debug_velocity.txt', status='replace', iostat=ios)
      write(iu, '(A,I8)') '# nCells=', nCells
      do i = 1, nCells
         write(iu, '(2(ES15.7,1X))') velocity(i, 1), velocity(i, 2)
      end do
      close(iu)

      call mpas_log_write('GNN debug text files written.')

   end subroutine dump_gnn_debug


end module li_velocity_gnn
