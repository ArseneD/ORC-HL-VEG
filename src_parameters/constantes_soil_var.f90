! =================================================================================================================================
! MODULE 	: constantes_soil_var
!
! CONTACT       : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "constantes_soil_var" module contains the parameters related to soil and hydrology.
!!
!!\n DESCRIPTION : The non saturated hydraulic properties are defined from the  
!!                 formulations of van Genuchten (1980) and Mualem (1976), combined as  
!!                 explained in d'Orgeval (2006). \n
!!                 The related parameters for three main  
!!                 soil textures (coarse, medium and fine) come from Carsel and Parrish  
!!                 (1988).
!!
!! RECENT CHANGE(S): Sonke Zaehle changed hcrit_litter value according to Shilong Piao from 0.03 to 0.08, 080806
!!
!! REFERENCE(S)	:
!!- Roger A.Pielke, (2002), Mesoscale meteorological modeling, Academic Press Inc. 
!!- Polcher, J., Laval, K., Dümenil, L., Lean, J., et Rowntree, P. R. (1996).
!! Comparing three land surface schemes used in general circulation models. Journal of Hydrology, 180(1-4), 373--394.
!!- Ducharne, A., Laval, K., et Polcher, J. (1998). Sensitivity of the hydrological cycle
!! to the parametrization of soil hydrology in a GCM. Climate Dynamics, 14, 307--327. 
!!- Rosnay, P. de et Polcher, J. (1999). Modelling root water uptake in a complex land surface
!! scheme coupled to a GCM. Hydrol. Earth Syst. Sci., 2(2/3), 239--255.
!!- d'Orgeval, T. et Polcher, J. (2008). Impacts of precipitation events and land-use changes
!! on West African river discharges during the years 1951--2000. Climate Dynamics, 31(2), 249--262. 
!!- Carsel, R. and Parrish, R.: Developing joint probability distributions of soil water
!! retention characteristics, Water Resour. Res.,24, 755–769, 1988.
!!- Mualem Y (1976) A new model for predicting the hydraulic conductivity  
!! of unsaturated porous media. Water Resources Research 12(3):513-522
!!- Van Genuchten M (1980) A closed-form equation for predicting the  
!! hydraulic conductivity of unsaturated soils. Soil Sci Soc Am J, 44(5):892-898
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================

MODULE constantes_soil_var

  USE defprec

  IMPLICIT NONE


  LOGICAL, SAVE             :: check_waterbal=.TRUE.    !! The check the water balance (true/false)
!$OMP THREADPRIVATE(check_waterbal)

  !! Dimensioning parameters

  INTEGER(i_std), SAVE      :: ngrnd                    !! Number of soil level for thermo (unitless)
  !$OMP THREADPRIVATE(ngrnd)
  INTEGER(i_std),PARAMETER  :: ndeep=32
  REAL(r_std), PARAMETER :: zalph=1.18
  REAL(r_std), PARAMETER          :: z_deepsoil = 2.    !!depth below which soil humidity is set to fixed values
  INTEGER(i_std), PARAMETER :: nbdl=11                  !! Number of diagnostic levels in the soil
                                                        !! To compare hydrologic variables with tag 1.6 and lower,
                                                        !! set nbdl to 6 : INTEGER(i_std), PARAMETER :: nbdl = 6
                                                        !! (unitless)
  INTEGER(i_std), PARAMETER :: nslm=11                  !! Number of levels in CWRR (unitless)

  REAL(r_std), SAVE         :: dpu_max                  !! Maximum depth of soil reservoir (m). Default value is set 
                                                        !! in intsurf_config depending on Choisnel(4m) or CWRR(2m)
!$OMP THREADPRIVATE(dpu_max)

  !! Number of soil classes

  INTEGER(i_std), PARAMETER :: ntext=3                  !! Number of soil textures (Silt, Sand, Clay)
  INTEGER(i_std), PARAMETER :: nstm=3                   !! Number of soil tiles (unitless)
  CHARACTER(LEN=30)         :: soil_classif             !! Type of classification used for the map of soil types.
                                                        !! It must be consistent with soil file given by 
                                                        !! SOILCLASS_FILE parameter.
  INTEGER(i_std), PARAMETER :: nscm_fao=3               !! For FAO Classification (unitless)
  INTEGER(i_std), PARAMETER :: nscm_usda=12             !! For USDA Classification (unitless)
  INTEGER(i_std), SAVE      :: nscm=nscm_fao            !! Default value for nscm
!$OMP THREADPRIVATE(nscm)

  !! Parameters for soil thermodynamics

  REAL(r_std), SAVE :: so_capa_dry = 1.80e+6            !! Dry soil Heat capacity of soils 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex 
!$OMP THREADPRIVATE(so_capa_dry)
  REAL(r_std), SAVE :: so_cond_dry = 0.40               !! Dry soil Thermal Conductivity of soils
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(so_cond_dry)
  REAL(r_std), SAVE :: so_capa_wet = 3.03e+6            !! Wet soil Heat capacity of soils 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(so_capa_wet)
  REAL(r_std), SAVE :: so_cond_wet = 1.89               !! Wet soil Thermal Conductivity of soils 
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex 
!$OMP THREADPRIVATE(so_cond_wet)
  REAL(r_std), SAVE :: sn_cond = 0.3                    !! Thermal Conductivity of snow 
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex  
!$OMP THREADPRIVATE(sn_cond)
  REAL(r_std), SAVE :: sn_dens = 330.0                  !! Snow density for the soil thermodynamics
                                                        !! (kg/m3)
!$OMP THREADPRIVATE(sn_dens)
  REAL(r_std), SAVE :: sn_capa                          !! Heat capacity for snow 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(sn_capa)
  REAL(r_std), PARAMETER     :: poros_org = 0.92 !! for now just a number from dmitry's code
  REAL(r_std), PARAMETER                 :: cond_solid_org = 0.25 !! W/m/K from Farouki via Lawrence and Slater
  REAL(r_std), PARAMETER :: so_cond_dry_org = 0.25          !! W/m/K from Farouki via Lawrence and Slater
  REAL(r_std), PARAMETER :: so_capa_dry_org = 2.5e6         !! J/K/m^3 from Farouki via Lawrence and Slater

!! Arsene 18-01-2016 - START - Add new thermal capacity and conductivity for mosses
  REAL(r_std), SAVE :: so_capa_dry_moss = 0.29e+6       !! Dry Moss Heat capacity From Soudzilovskaia et al., 2013
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex 
!$OMP THREADPRIVATE(so_capa_dry_moss)
  REAL(r_std), SAVE :: so_capa_wet_moss = 4.29e+6       !! Wet moss Heat capacity From Soudzilovskaia et al., 2013
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(so_capa_wet_moss)
  REAL(r_std), SAVE :: so_capa_ice_moss = 3.26e+6       !! Ice moss Heat capacity Deduce from wet and ice ratio (without moss)
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(so_capa_ice_moss)
  REAL(r_std), SAVE :: so_cond_dry_moss = 0.092         !! Dry moss Thermal Conductivity of soils From Soudzilovskaia et al., 2013
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(so_cond_dry_moss)
  REAL(r_std), SAVE :: cond_solid_moss = 0.754          !! Wet moss Thermal Conductivity From Soudzilovskaia et al., 2013
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(cond_solid_moss)
  REAL(r_std), SAVE :: cond_ice_moss = 0.715            !! Ice moss Thermal Conductivity from wet and ice ratio (without moss)
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(cond_ice_moss)
!! Arsene 18-01-2016 - END - Add new thermal capacity and conductivity for mosses

  !REAL(r_std),PARAMETER :: sn_capa = 2100.0_r_std*sn_dens   !! Heat capacity
  !for snow @tex $(J.m^{-3}.K^{-1})$ @endtex
  REAL(r_std), PARAMETER   :: soilc_max =  130000.          !! g/m^3 from lawrence and slater



  !! Specific parameters for the Choisnel hydrology

  REAL(r_std), SAVE :: min_drain = 0.001                !! Diffusion constant for the slow regime
                                                        !! (This is for the diffusion between reservoirs)
                                                        !! @tex $(kg.m^{-2}.dt^{-1})$ @endtex
!$OMP THREADPRIVATE(min_drain)
  REAL(r_std), SAVE :: max_drain = 0.1                  !! Diffusion constant for the fast regime 
                                                        !! @tex $(kg.m^{-2}.dt^{-1})$ @endtex
!$OMP THREADPRIVATE(max_drain)
  REAL(r_std), SAVE :: exp_drain = 1.5                  !! The exponential in the diffusion law (unitless)
!$OMP THREADPRIVATE(exp_drain)
  REAL(r_std), SAVE :: qsintcst = 0.1                   !! Transforms leaf area index into size of interception reservoir
                                                        !! (unitless)
!$OMP THREADPRIVATE(qsintcst)
  REAL(r_std), SAVE :: qsintcst_moss_coef = 5.          !! Coefficient to compute value of qsintcst for moss: qsintcst_moss = qsintcst_moss_coef * qsintcst !! Arsene 20-01-2016 - Add
                                                        !! (unitless)
!$OMP THREADPRIVATE(qsintcst_moss_coef)
  REAL(r_std), SAVE :: mx_eau_nobio = 150.              !! Volumetric available soil water capacity in nobio fractions
                                                        !! @tex $(kg.m^{-3} of soil)$ @endtex
!$OMP THREADPRIVATE(mx_eau_nobio)
  REAL(r_std), SAVE :: rsol_cste = 33.E3                !! Constant in the computation of resistance for bare soil evaporation
                                                        !! @tex $(s.m^{-2})$ @endtex
!$OMP THREADPRIVATE(rsol_cste)
  REAL(r_std), SAVE :: hcrit_litter=0.08_r_std          !! Scaling depth for litter humidity (m)
!$OMP THREADPRIVATE(hcrit_litter)

  LOGICAL, SAVE :: moss_water_leack_ok = .false.         ! Activate or not the leave interception increase and leak for mosses     !! Arsene 23-02-2016 - Add
!$OMP THREADPRIVATE(reinf_slope_moss_ok)
  REAL(r_std), SAVE :: moss_water_leack = 14.           !! Number of day to loss the leave interception water by leack to the soil (via precisol), for moss (day) !! Arsene 20-01-2016 - Add
!$OMP THREADPRIVATE(moss_water_leack)
  LOGICAL, SAVE :: reinf_slope_moss_ok = .true.         !! Activate or not the different reinf_slope for moss (impact the runoff)  !! Arsene 23-02-2016 - Add
!$OMP THREADPRIVATE(reinf_slope_moss_ok)
  REAL(r_std), SAVE :: reinf_slope_moss = 0.6          !! Value for reinf_slope (not topography dependent), [0-1] (unitless)       !! Arsene 09-02-2016 - Add
!$OMP THREADPRIVATE(reinf_slope_moss)

!  REAL(r_std), SAVE :: reinf_slope_moss_coef = 2.       !! Coeficient to increase reinf_slope. 0 --> no impact. 1 --> slope x2. +1 --> slope +100% (unitless) !! Arsene 20-01-2016 - Add
!!$OMP THREADPRIVATE(reinf_slope_moss_coef)
!  REAL(r_std), SAVE :: moss_water_runoff = 1.           !! Number of day to loss the water from water2infil (day)             !! Arsene 09-02-2016 - Add
!!$OMP THREADPRIVATE(moss_water_runoff)

  LOGICAL, SAVE :: thermal_property_moss_ok = .true.     !! Activate or not the specific thermosoil properties for mosses (from Soudzilovskaia et al., 2013) !! Arsene 24-02-2016 - Add
!$OMP THREADPRIVATE(thermal_property_moss_ok)


  !! Parameters specific for the CWRR hydrology.

  !!  1. Parameters for FAO Classification

  !! Parameters for soil type distribution

  REAL(r_std),DIMENSION(nscm_fao),SAVE :: soilclass_default_fao = &   !! Default soil texture distribution for fao :
 & (/ 0.28, 0.52, 0.20 /)                                             !! in the following order : COARSE, MEDIUM, FINE (unitless)
!$OMP THREADPRIVATE(soilclass_default_fao)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: nvan_fao = &            !! Van genuchten coefficient n
 & (/ 1.89_r_std, 1.56_r_std, 1.31_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: avan_fao = &            !! Van genuchten coefficient a (mm^{-1}) 
  & (/ 0.0075_r_std, 0.0036_r_std, 0.0019_r_std /) 

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcr_fao = &             !! Residual soil water content
 & (/ 0.065_r_std, 0.078_r_std, 0.095_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcs_fao = &             !! Saturated soil water content
 & (/ 0.41_r_std, 0.43_r_std, 0.41_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: ks_fao = &              !! Hydraulic conductivity Saturation (mm/d)
 & (/ 1060.8_r_std, 249.6_r_std, 62.4_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: pcent_fao = &           !! Fraction of saturated volumetric soil moisture 
 & (/ 0.5_r_std, 0.5_r_std, 0.5_r_std /)                               !! above which transpir is max (0-1, unitless)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: free_drain_max_fao = &  !! Max value of the permeability coeff at 
 & (/ 1.0_r_std, 1.0_r_std, 1.0_r_std /)                               !! the bottom of the soil

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcf_fao = &             !! Volumetric water content field capacity
 & (/ 0.32_r_std, 0.32_r_std, 0.32_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mcw_fao = &             !! Volumetric water content Wilting pt
 & (/ 0.065_r_std, 0.078_r_std, 0.095_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mc_awet_fao = &         !! Vol. wat. cont. above which albedo is cst
 & (/ 0.25_r_std, 0.25_r_std, 0.25_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_fao) :: mc_adry_fao = &         !! Vol. wat. cont. below which albedo is cst
 & (/ 0.1_r_std, 0.1_r_std, 0.1_r_std /)


  !!  2. Parameters for USDA Classification

  !! Parameters for soil type distribution :
  !! Sand, Loamy Sand, Sandy Loam, Silt Loam, Silt, Loam, Sandy Clay Loam, Silty Clay Loam, Clay Loam, Sandy Clay, Silty Clay, Clay

  REAL(r_std),DIMENSION(nscm_usda),SAVE :: soilclass_default_usda = &    !! Default soil texture distribution in the following order :
 & (/ 0.28, 0.52, 0.20, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)   !!    sand, loam and clay ??? OR COARSE, MEDIUM, FINE???
!$OMP THREADPRIVATE(soilclass_default_usda)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: nvan_usda = &            !! Van genuchten coefficient n
 & (/ 2.68_r_std, 2.28_r_std, 1.89_r_std, 1.41_r_std, &
 &    1.37_r_std, 1.56_r_std, 1.48_r_std, 1.23_r_std, &
 &    1.31_r_std, 1.23_r_std, 1.09_r_std, 1.09_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: avan_usda = &            !! Van genuchten coefficient a (mm^{-1})
 & (/ 0.0145_r_std, 0.0124_r_std, 0.0075_r_std, 0.0020_r_std, &
 &    0.0016_r_std, 0.0036_r_std, 0.0059_r_std, 0.0010_r_std, &
 &    0.0019_r_std, 0.0027_r_std, 0.0005_r_std, 0.0008_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcr_usda = &             !! Residual soil water content
 & (/ 0.045_r_std, 0.057_r_std, 0.065_r_std, 0.067_r_std, &
 &    0.034_r_std, 0.078_r_std, 0.100_r_std, 0.089_r_std, &
 &    0.095_r_std, 0.100_r_std, 0.070_r_std, 0.068_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcs_usda = &             !! Saturated soil water content
 & (/ 0.43_r_std, 0.41_r_std, 0.41_r_std, 0.45_r_std, &
 &    0.46_r_std, 0.43_r_std, 0.39_r_std, 0.43_r_std, &
 &    0.41_r_std, 0.38_r_std, 0.36_r_std, 0.38_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: ks_usda = &              !! Hydraulic conductivity Saturation (mm/d)
 & (/ 7128.0_r_std, 3501.6_r_std, 1060.8_r_std, 108.0_r_std, &
 &    60.0_r_std, 249.6_r_std, 314.4_r_std, 16.8_r_std, &
 &    62.4_r_std, 28.8_r_std, 4.8_r_std, 48.0_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: pcent_usda = &           !! Soil moisture above which transpir is max
 & (/ 0.5_r_std, 0.5_r_std, 0.5_r_std, 0.5_r_std, &
 &    0.5_r_std, 0.5_r_std, 0.5_r_std, 0.5_r_std, &
 &    0.5_r_std, 0.5_r_std, 0.5_r_std, 0.5_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: free_drain_max_usda = &  !! Max value of the permeability coeff at 
 & (/ 1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std, &                      !! the bottom of the soil
 &    1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std, &
 &    1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcf_usda = &             !! Volumetric water content field capacity
 & (/ 0.32_r_std, 0.32_r_std, 0.32_r_std, 0.32_r_std, &
 &    0.32_r_std, 0.32_r_std, 0.32_r_std, 0.32_r_std, &
 &    0.32_r_std, 0.32_r_std, 0.32_r_std, 0.32_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcw_usda = &             !! Volumetric water content Wilting pt
 & (/ 0.10_r_std, 0.10_r_std, 0.10_r_std, 0.10_r_std, &
 &    0.10_r_std, 0.10_r_std, 0.10_r_std, 0.10_r_std, &
 &    0.10_r_std, 0.10_r_std, 0.10_r_std, 0.10_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mc_awet_usda = &         !! Vol. wat. cont. above which albedo is cst
 & (/ 0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &
 &    0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &
 &    0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mc_adry_usda = &         !! Vol. wat. cont. below which albedo is cst
 & (/ 0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &
 &    0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &
 &    0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std /)


  !! Parameters for the numerical scheme used by CWRR

  INTEGER(i_std), PARAMETER :: imin = 1                                 !! CWRR linearisation (unitless)
  INTEGER(i_std), PARAMETER :: nbint = 50                               !! Number of interval for CWRR (unitless)
  INTEGER(i_std), PARAMETER :: imax = nbint+1                           !! Number of points for CWRR (unitless)
  REAL(r_std), PARAMETER    :: w_time = 1.0_r_std                       !! Time weighting for discretisation (unitless)


  !! Diagnostic variables

  REAL(r_std),DIMENSION(nbdl),SAVE :: diaglev                           !! The lower limit of the layer on which soil moisture
                                                                        !! (relative) and temperature are going to be diagnosed.
                                                                        !! These variables are made for transfering the information
                                                                        !! to the biogeophyical processes modelled in STOMATE. 
                                                                        !! (unitless)
!$OMP THREADPRIVATE(diaglev)

  !! Variables related to soil freezing, in thermosoil : 
  LOGICAL, SAVE        :: ok_Ecorr                    !! Flag for energy conservation correction
  LOGICAL, SAVE        :: ok_freeze_thermix           !! Flag to activate thermal part of the soil freezing scheme
  LOGICAL, SAVE        :: read_reftemp                !! Flag to initialize soil temperature using climatological temperature
  LOGICAL, SAVE        :: read_permafrost_map         !! Read information about ice content, overburden and permafrost type from IPA map
  REAL(r_std), SAVE    :: poros                       !! Soil porosity (from USDA classification, mean value)(-)
  REAL(r_std), SAVE    :: fr_dT                       !! Freezing window (K)

  !! Variables related to soil freezing, in diffuco : 
  LOGICAL, SAVE        ::  ok_snowfact                !! Activate snow smoothering

  !! Variables related to soil freezing, in hydrol : 
  LOGICAL, SAVE        :: ok_freeze_cwrr              !! CWRR freezing scheme by I. Gouttevin
  LOGICAL, SAVE        :: ok_thermodynamical_freezing !! Calculate frozen fraction thermodynamically

  
END MODULE constantes_soil_var
