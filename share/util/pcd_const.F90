module pcd_const
    use shr_kind_mod, only : R8 => shr_kind_r8

    implicit none
    private

    public :: glacier_ice_density_reference
    public :: freshwater_kinematic_viscosity_at_freezing_reference
    public :: ocean_freezing_temperature_depth_dependence_reference

    real(R8), parameter :: glacier_ice_density_reference = 0.910e3_R8 ! reference glacier-ice density ~ kg/m^3
    real(R8), parameter :: freshwater_kinematic_viscosity_at_freezing_reference = 1.793e-6_R8 ! kinematic viscosity of freshwater at 0 C ~ m^2/s
    real(R8), parameter :: ocean_freezing_temperature_depth_dependence_reference = -7.28617688e-4_R8 ! pressure-driven ocean freezing temperature depth coefficient ~ C/m
end module pcd_const
