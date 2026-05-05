module pcd
    implicit none

    integer, parameter :: dp = selected_real_kind(12)

    real(dp), parameter :: glacier_ice_density_reference = 910.0_dp
    real(dp), parameter :: freshwater_kinematic_viscosity_at_freezing_reference = 1.793e-6_dp
    real(dp), parameter :: ocean_freezing_temperature_pressure_depth_dependence_reference = -7.28617688e-4_dp
end module pcd
