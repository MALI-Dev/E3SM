! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!***********************************************************************
!
!  mpas_plane_slope_operators
!
!> \brief   MPAS plane-fit slope operators
!> \author  MPAS Development Team
!> \date    06/09/26
!> \details
!> This module computes per-cell unsigned slope magnitude from a linear
!> plane fit over a contiguous ring stencil on planar meshes.
!
!-----------------------------------------------------------------------
module mpas_plane_slope_operators

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   implicit none
   private
   save

   public :: mpas_max_plane_slope_ring

contains

!***********************************************************************
!
!  routine mpas_max_plane_slope_ring
!
!> \brief   Compute unsigned plane-fit slope on configurable ring stencils
!> \author  MPAS Development Team
!> \date    06/09/26
!> \details
!> For each owned cell, this routine fits a linear plane to neighboring
!> cell-center values over a contiguous ring stencil and returns the
!> slope magnitude sqrt(a^2+b^2) where z = a*x + b*y + c in local
!> cell-centered coordinates. The routine is planar-mesh only.
!
!-----------------------------------------------------------------------
   subroutine mpas_max_plane_slope_ring(fieldCell, meshPool, configPool, ringCount, maxSlope, validMask)!{{{

      real(kind=RKIND), dimension(:), intent(in) :: fieldCell  ! Cell-centered scalar field to fit
      type(mpas_pool_type), intent(in) :: meshPool             ! Mesh geometry/connectivity pool
      type(mpas_pool_type), intent(in) :: configPool           ! Config pool (used for halo-depth guard)
      integer, intent(in) :: ringCount                         ! Requested contiguous ring depth
      real(kind=RKIND), dimension(:), intent(out) :: maxSlope  ! Output unsigned slope magnitude per solve cell
      logical, dimension(:), intent(in), optional :: validMask ! Optional mask of cells eligible for stencils

      integer, pointer :: nCells, nCellsSolve
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      real(kind=RKIND), dimension(:), pointer :: xCell, yCell
      logical, pointer :: on_a_sphere
      integer, pointer :: config_num_halos

      integer, allocatable, dimension(:) :: visited, ringByCell, frontier, nextFrontier, stencil
      integer :: iCell, i, iRing, iFront, iNbr, iCellNbr
      integer :: frontierCount, nextFrontierCount, stencilCount
      integer :: maxBuiltRing, effectiveRing
      integer :: visitTag
      real(kind=RKIND) :: aCoeff, bCoeff
      logical :: fitOK

      ! Gather mesh/config metadata and connectivity needed to build local ring stencils.
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(configPool, 'config_num_halos', config_num_halos)

      ! Input and configuration validation.
      if (ringCount < 1) then
         call mpas_log_write('mpas_max_plane_slope_ring: ringCount must be >= 1.', MPAS_LOG_CRIT)
      end if

      if (ringCount > config_num_halos) then
         call mpas_log_write('mpas_max_plane_slope_ring: requested ring count $i exceeds config_num_halos $i.', &
            MPAS_LOG_CRIT, intArgs=(/ringCount, config_num_halos/))
      end if

      if (on_a_sphere) then
         call mpas_log_write('mpas_max_plane_slope_ring currently supports planar meshes only.', MPAS_LOG_CRIT)
      end if

      if (size(maxSlope) < nCellsSolve) then
         call mpas_log_write('mpas_max_plane_slope_ring: output array maxSlope is too small.', MPAS_LOG_CRIT)
      end if

      if (size(fieldCell) < nCells) then
         call mpas_log_write('mpas_max_plane_slope_ring: input array fieldCell is too small.', MPAS_LOG_CRIT)
      end if

      if (present(validMask)) then
         if (size(validMask) < nCells) then
            call mpas_log_write('mpas_max_plane_slope_ring: optional validMask is too small.', MPAS_LOG_CRIT)
         end if
      end if

      maxSlope = 0.0_RKIND

      ! Work arrays are sized once to avoid repeated allocation in the per-cell loop.
      allocate(visited(nCells))
      allocate(ringByCell(nCells))
      allocate(frontier(nCells))
      allocate(nextFrontier(nCells))
      allocate(stencil(nCells))

      visited = 0
      ringByCell = -1
      visitTag = 0

      do iCell = 1, nCellsSolve

         ! If the center cell is masked out, return the requested fallback immediately.
         if (present(validMask)) then
            if (.not. validMask(iCell)) then
               maxSlope(iCell) = 0.0_RKIND
               cycle
            end if
         end if

         visitTag = visitTag + 1
         if (visitTag == huge(visitTag)) then
            visited = 0
            visitTag = 1
         end if

         frontierCount = 1
         frontier(1) = iCell
         visited(iCell) = visitTag
         ringByCell(iCell) = 0

         stencilCount = 1
         stencil(1) = iCell
         maxBuiltRing = 0

         ! Breadth-first expansion over cellsOnCell builds a contiguous stencil up to ringCount.
         do iRing = 1, ringCount
            nextFrontierCount = 0

            do iFront = 1, frontierCount
               do iNbr = 1, nEdgesOnCell(frontier(iFront))
                  iCellNbr = cellsOnCell(iNbr, frontier(iFront))

                  if (iCellNbr < 1 .or. iCellNbr > nCells) cycle
                  if (visited(iCellNbr) == visitTag) cycle

                  if (present(validMask)) then
                     if (.not. validMask(iCellNbr)) cycle
                  end if

                  visited(iCellNbr) = visitTag
                  ringByCell(iCellNbr) = iRing

                  nextFrontierCount = nextFrontierCount + 1
                  nextFrontier(nextFrontierCount) = iCellNbr

                  stencilCount = stencilCount + 1
                  stencil(stencilCount) = iCellNbr
               end do
            end do

            if (nextFrontierCount == 0) exit

            frontier(1:nextFrontierCount) = nextFrontier(1:nextFrontierCount)
            frontierCount = nextFrontierCount
            maxBuiltRing = iRing
         end do

         ! Reduce outer rings until a stable linear fit is available.
         effectiveRing = maxBuiltRing
         fitOK = .false.

         do while (effectiveRing >= 1 .and. .not. fitOK)
            call mpas_try_plane_fit(iCell, effectiveRing, stencil, stencilCount, ringByCell, &
                                    xCell, yCell, fieldCell, aCoeff, bCoeff, fitOK)
            if (.not. fitOK) effectiveRing = effectiveRing - 1
         end do

         if (fitOK) then
            maxSlope(iCell) = sqrt(aCoeff*aCoeff + bCoeff*bCoeff)
         else
            maxSlope(iCell) = 0.0_RKIND
         end if

      end do

      ! Explicit cleanup of workspace arrays.
      deallocate(visited)
      deallocate(ringByCell)
      deallocate(frontier)
      deallocate(nextFrontier)
      deallocate(stencil)

   end subroutine mpas_max_plane_slope_ring!}}}


!***********************************************************************
!
!  routine mpas_try_plane_fit
!
!> \brief   Attempt linear plane fit over a reduced ring stencil
!> \author  MPAS Development Team
!> \date    06/09/26
!> \details
!> Fits dz = a*dx + b*dy + c and returns (a,b) when matrix is well posed.
!
!-----------------------------------------------------------------------
   subroutine mpas_try_plane_fit(iCenter, maxRing, stencil, stencilCount, ringByCell, xCell, yCell, fieldCell, aCoeff, bCoeff, fitOK)!{{{

      integer, intent(in) :: iCenter, maxRing, stencilCount          ! Center cell index and reduced-ring depth
      integer, dimension(:), intent(in) :: stencil, ringByCell       ! Candidate stencil cells and ring assignment
      real(kind=RKIND), dimension(:), intent(in) :: xCell, yCell, fieldCell
      real(kind=RKIND), intent(out) :: aCoeff, bCoeff                ! Plane slopes dz/dx and dz/dy
      logical, intent(out) :: fitOK                                  ! .true. when solve is well-posed

      integer :: i, iCell
      integer :: nPoints
      real(kind=RKIND) :: x0, y0, z0, dx, dy, dz
      real(kind=RKIND), dimension(3,3) :: mat
      real(kind=RKIND), dimension(3) :: rhs, soln
      real(kind=RKIND) :: detScale

      fitOK = .false.
      aCoeff = 0.0_RKIND
      bCoeff = 0.0_RKIND

      x0 = xCell(iCenter)
      y0 = yCell(iCenter)
      z0 = fieldCell(iCenter)

      mat = 0.0_RKIND
      rhs = 0.0_RKIND
      nPoints = 0

      ! Build normal equations for dz = a*dx + b*dy + c in center-relative coordinates.
      do i = 1, stencilCount
         iCell = stencil(i)
         if (ringByCell(iCell) > maxRing) cycle

         dx = xCell(iCell) - x0
         dy = yCell(iCell) - y0
         dz = fieldCell(iCell) - z0

         mat(1,1) = mat(1,1) + dx*dx
         mat(1,2) = mat(1,2) + dx*dy
         mat(1,3) = mat(1,3) + dx

         mat(2,2) = mat(2,2) + dy*dy
         mat(2,3) = mat(2,3) + dy

         mat(3,3) = mat(3,3) + 1.0_RKIND

         rhs(1) = rhs(1) + dx*dz
         rhs(2) = rhs(2) + dy*dz
         rhs(3) = rhs(3) + dz

         nPoints = nPoints + 1
      end do

      ! Fewer than 3 points cannot constrain a linear plane robustly.
      if (nPoints < 3) return

      mat(2,1) = mat(1,2)
      mat(3,1) = mat(1,3)
      mat(3,2) = mat(2,3)

      ! Use matrix scale to set a relative singularity tolerance.
      detScale = max(1.0_RKIND, maxval(abs(mat)))
      call mpas_solve_3x3(mat, rhs, detScale*1.0e-12_RKIND, soln, fitOK)
      if (.not. fitOK) return

      aCoeff = soln(1)
      bCoeff = soln(2)

   end subroutine mpas_try_plane_fit!}}}


!***********************************************************************
!
!  routine mpas_solve_3x3
!
!> \brief   Solve a 3x3 linear system with Cramer's rule
!> \author  MPAS Development Team
!> \date    06/09/26
!
!-----------------------------------------------------------------------
   subroutine mpas_solve_3x3(a, b, detTol, x, ok)!{{{

      real(kind=RKIND), dimension(3,3), intent(in) :: a     ! Coefficient matrix
      real(kind=RKIND), dimension(3), intent(in) :: b       ! Right-hand side
      real(kind=RKIND), intent(in) :: detTol                ! Determinant tolerance for singularity
      real(kind=RKIND), dimension(3), intent(out) :: x      ! Solution vector
      logical, intent(out) :: ok                            ! .true. when matrix is invertible

      real(kind=RKIND) :: detA, detX1, detX2, detX3
      real(kind=RKIND), dimension(3,3) :: work

      ! Cramer's rule is sufficient here because the system is fixed-size (3x3).
      detA = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) - &
             a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1)) + &
             a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))

      if (abs(detA) <= detTol) then
         ok = .false.
         x = 0.0_RKIND
         return
      end if

      work = a
      work(:,1) = b
      detX1 = work(1,1)*(work(2,2)*work(3,3) - work(2,3)*work(3,2)) - &
              work(1,2)*(work(2,1)*work(3,3) - work(2,3)*work(3,1)) + &
              work(1,3)*(work(2,1)*work(3,2) - work(2,2)*work(3,1))

      work = a
      work(:,2) = b
      detX2 = work(1,1)*(work(2,2)*work(3,3) - work(2,3)*work(3,2)) - &
              work(1,2)*(work(2,1)*work(3,3) - work(2,3)*work(3,1)) + &
              work(1,3)*(work(2,1)*work(3,2) - work(2,2)*work(3,1))

      work = a
      work(:,3) = b
      detX3 = work(1,1)*(work(2,2)*work(3,3) - work(2,3)*work(3,2)) - &
              work(1,2)*(work(2,1)*work(3,3) - work(2,3)*work(3,1)) + &
              work(1,3)*(work(2,1)*work(3,2) - work(2,2)*work(3,1))

      x(1) = detX1 / detA
      x(2) = detX2 / detA
      x(3) = detX3 / detA
      ok = .true.

   end subroutine mpas_solve_3x3!}}}

end module mpas_plane_slope_operators
