module pcd_const
    use shr_kind_mod, only : R8 => shr_kind_r8

    implicit none

    real(R8), parameter :: glacier_ice_density_reference = 910.0_R8
    real(R8), parameter :: freshwater_kinematic_viscosity_at_freezing_reference = 1.793e-6_R8
    real(R8), parameter :: ocean_freezing_temperature_pressure_depth_dependence_reference = -7.28617688e-4_R8
end module pcd_const
