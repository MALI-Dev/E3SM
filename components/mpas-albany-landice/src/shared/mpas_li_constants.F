! Copyright (c) 2015-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_constants
!
!> \brief MPAS land ice specific constants
!> \author Matthew Hoffman
!> \date   17 Feb. 2015
!> \details
!>  This module contains constants specific to the land ice model.
!
!-----------------------------------------------------------------------

module li_constants

   use mpas_derived_types
   use mpas_kind_types

#ifdef MPAS_ESM_SHR_CONST
    use shr_const_mod, only: &
           cp_ice => SHR_CONST_CPICE,&
           latent_heat_ice => SHR_CONST_LATICE,&
           triple_point => SHR_CONST_TKTRIP,&
           rho_water => SHR_CONST_RHOFW, &
           cp_freshwater => SHR_CONST_CPFW, &
           cp_seawater => SHR_CONST_CPSW, &
           pii => SHR_CONST_PI, &
           gravity => SHR_CONST_G, &
           shr_gas_constant => SHR_CONST_RGAS, &
           iceConductivity => SHR_CONST_CONDICE, &
           kelvin_to_celsius => SHR_CONST_TKFRZ, &
           seconds_per_day => SHR_CONST_CDAY
    implicit none
    save

    real (kind=RKIND), parameter, public :: idealGasConstant = shr_gas_constant * 1.0e-3_RKIND  !< ideal gas constant (J mol^-1 K^-1), converted from SHR_CONST_RGAS (J kmol^-1 K^-1)

#else

   implicit none
   save

   ! physical constants
   real (kind=RKIND), parameter, public :: cp_ice = 2.11727e3_RKIND        !< heat capacity of ice (J/kg/K), aligned with SHR_CONST_CPICE
   real (kind=RKIND), parameter, public :: latent_heat_ice = 3.337e5_RKIND !< Latent heat of melting of ice (J/kg), aligned with SHR_CONST_LATICE
   real (kind=RKIND), parameter, public :: triple_point = 273.16_RKIND  !< Triple point of water (K)
   real (kind=RKIND), parameter, public :: rho_water = 1000.0_RKIND     !< Density of fresh water (kg m^-3)
   real (kind=RKIND), parameter, public :: cp_freshwater = 4.188e3_RKIND        !< heat capacity of freshwater (J/kg/K)
   real (kind=RKIND), parameter, public :: cp_seawater = 3.996e3_RKIND        !< heat capacity of seawater (J/kg/K), aligned with SHR_CONST_CPSW
   real (kind=RKIND), parameter, public :: pii = 3.14159265358979323846_RKIND !< Constant: Pi
   real (kind=RKIND), parameter, public :: gravity = 9.80616_RKIND       !< Constant: Acceleration due to gravity [m s-2]
   real (kind=RKIND), parameter, private :: shr_gas_constant = 8314.467591_RKIND
   real (kind=RKIND), parameter, public :: idealGasConstant = shr_gas_constant * 1.0e-3_RKIND  !< ideal gas constant (J mol^-1 K^-1), converted from SHR_CONST_RGAS
   real (kind=RKIND), parameter, public :: iceConductivity = 2.1_RKIND  !< thermal conductivity of ice (W m^-1 K^-1)
   real (kind=RKIND), parameter, public :: kelvin_to_celsius = 273.15_RKIND    !< factor to convert Kelvin to Celsius
   real (kind=RKIND), parameter, public :: seconds_per_day = 86400.0_RKIND  !< Conversion from days to seconds

#endif

   real (kind=RKIND), parameter, public :: glacier_ice_density = 910.0_RKIND !< Density of glacier ice (kg m^-3)
   real (kind=RKIND), parameter, private :: oceanFreezingTempPressureDependence = -7.43e-8_RKIND
   real (kind=RKIND), parameter, public :: kin_visc_water = 1.793e-6_RKIND  !< kinematic viscosity of freshwater at 0 deg C (m^2 s^-1)
   real (kind=RKIND), parameter, public :: oceanFreezingTempDepthDependence = oceanFreezingTempPressureDependence * rho_water * gravity  !< Rate of change of freezing temperature with depth (deg m^-1)

   real (kind=RKIND), parameter, public :: &
        oceanFreezingTempSurface = -1.92_RKIND  !< Freezing temperature of seawater (deg C) at surface pressure,
                                                 !< given S = 35 PSU

   real (kind=RKIND), parameter, public :: iceMeltingPointPressureDependence = 7.9e-8_RKIND
      !< Clausius-Clapeyron constant: Dependence of ice melting point on pressure (K Pa^-1)

   !  conversion factors
   real (kind=RKIND), parameter, public :: scyr = seconds_per_day * 365.0_RKIND       !< seconds in a 365-day year; used for diagnostics

!***********************************************************************


!***********************************************************************

end module li_constants
