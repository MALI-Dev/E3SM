module pcd_const
    use shr_kind_mod, only : R8 => shr_kind_r8

    implicit none

    real(R8), parameter :: glacier_ice_density_reference = 0.910e3_R8 ! [kg m^-3]
    real(R8), parameter :: freshwater_kinematic_viscosity_at_freezing_reference = 1.793e-6_R8 ! [m^2 s^-1]
    real(R8), parameter :: ocean_freezing_temperature_depth_dependence_reference = -7.28617688e-4_R8 ! [C m^-1]
end module pcd_const
