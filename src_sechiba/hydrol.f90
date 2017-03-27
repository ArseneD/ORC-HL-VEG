! =================================================================================================================================
! MODULE        : hydrol
!
! CONTACT       : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module computes the soil moisture processes on continental points. 
!!
!!\n DESCRIPTION : contains hydrol_init, hydrol_var_init, hydrol_waterbal, hydrol_alma,
!!                 hydrol_snow, hydrol_vegupd, hydrol_canop, hydrol_flood, hydrol_soil.
!!                 The assumption in this module is that very high vertical resolution is
!!                 needed in order to properly resolve the vertical diffusion of water in
!!                 the soils. Furthermore we have taken into account the sub-grid variability
!!                 of soil properties and vegetation cover by allowing the co-existence of
!!                 different soil moisture columns in the same grid box.
!!                 This routine was originaly developed by Patricia deRosnay.
!! 
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_sechiba/hydrol.f90 $
!! $Date: 2014-07-08 08:02:07 +0200 (Tue, 08 Jul 2014) $
!! $Revision: 2224 $
!! \n
!_ ================================================================================================================================

MODULE hydrol

  USE ioipsl
  USE xios_orchidee
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE sechiba_io
  USE slowproc
  USE grid
  USE explicitsnow

!pss:+USE TOPMODEL routines
!  USE extrac_cti
  USE ioipsl_para
  USE init_top
  USE hydro_subgrid
!pss:-


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: hydrol_main,hydrol_clear 

  !
  ! variables used inside hydrol module : declaration and initialisation
  !
  LOGICAL, SAVE                                   :: l_first_hydrol=.TRUE.   !! Initialisation has to be done one time
!$OMP THREADPRIVATE(l_first_hydrol)
  LOGICAL, SAVE                                   :: l_second_hydrol=.TRUE.  !! Initialisation has to be done one time
!$OMP THREADPRIVATE(l_second_hydrol)
  !
  LOGICAL, SAVE                                   :: check_cwrr=.FALSE.      !! The check the water balance
!$OMP THREADPRIVATE(check_cwrr)
  LOGICAL, SAVE                                   :: doponds=.FALSE.         !! Reinfiltration param.
!$OMP THREADPRIVATE(doponds)

!pss:+
  LOGICAL,SAVE                                    :: TOPM_calcul             !! flag of TOPMODEL usage
!$OMP THREADPRIVATE(TOPM_calcul)
!pss:-

  !
  CHARACTER(LEN=80) , SAVE                        :: var_name                !! To store variables names for I/O
!$OMP THREADPRIVATE(var_name)
  !
  REAL(r_std), PARAMETER                          :: allowed_err =  2.0E-8_r_std
  REAL(r_std), PARAMETER                          :: EPS1 = EPSILON(un)      !! A small number
  ! one dimension array allocated, computed, saved and got in hydrol module
  ! Values per soil type
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: nvan                !! Van Genuchten coeficients n
!$OMP THREADPRIVATE(nvan)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: avan                !! Van Genuchten coeficients a
!$OMP THREADPRIVATE(avan)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcr                 !! Residual humidity [m3/m3]
!$OMP THREADPRIVATE(mcr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcs                 !! Saturation humidity [m3/m3]
!$OMP THREADPRIVATE(mcs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: ks                  !! Hydraulic conductivity Saturation for each soil type 
                                                                         !! [mm/day] 
!$OMP THREADPRIVATE(ks)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: pcent               !! Fraction of saturated volumetric soil moisture above 
                                                                         !! which transpir is max 
!$OMP THREADPRIVATE(pcent)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: free_drain_max      !! Max value of the permeability coeff at the bottom of 
                                                                         !! the soil [mm/day] 
!$OMP THREADPRIVATE(free_drain_max)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcf                 !! Field capacity [m3/m3]
!$OMP THREADPRIVATE(mcf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mcw                 !! Wilting point [m3/m3]
!$OMP THREADPRIVATE(mcw)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mc_awet             !! Vol. wat. cont. above which albedo is cst [m3/m3]
!$OMP THREADPRIVATE(mc_awet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: mc_adry             !! Vol. wat. cont. below which albedo is cst [m3/m3]
!$OMP THREADPRIVATE(mc_adry)

  ! Values per grid point
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_water_beg    !! Total amount of water at start of time step
!$OMP THREADPRIVATE(tot_water_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_water_end    !! Total amount of water at end of time step
!$OMP THREADPRIVATE(tot_water_end)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_flux         !! Total water flux
!$OMP THREADPRIVATE(tot_flux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watveg_beg   !! Total amount of water on vegetation at start of time 
                                                                         !! step 
!$OMP THREADPRIVATE(tot_watveg_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watveg_end   !! Total amount of water on vegetation at end of time step
!$OMP THREADPRIVATE(tot_watveg_end)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watsoil_beg  !! Total amount of water in the soil at start of time step
!$OMP THREADPRIVATE(tot_watsoil_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_watsoil_end  !! Total amount of water in the soil at end of time step
!$OMP THREADPRIVATE(tot_watsoil_end)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snow_beg         !! Total amount of snow at start of time step
!$OMP THREADPRIVATE(snow_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snow_end         !! Total amount of snow at end of time step
!$OMP THREADPRIVATE(snow_end)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delsoilmoist     !! Change in soil moisture
!$OMP THREADPRIVATE(delsoilmoist)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delintercept     !! Change in interception storage
!$OMP THREADPRIVATE(delintercept)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: delswe           !! Change in SWE
!$OMP THREADPRIVATE(delswe)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: swi              !! Soil Wetness Index
!$OMP THREADPRIVATE(swi)

!pss+
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: mcs_grid              !! Saturation dim kjpindex
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: mcw_grid              !! Wilting point dim kjpindex  
!pss-

  ! array allocated, computed, saved and got in hydrol module
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mask_veget       !! zero/one when veget fraction is zero/higher
!$OMP THREADPRIVATE(mask_veget)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mask_soiltile    !! zero/one where soil tile is zero/higher 
!$OMP THREADPRIVATE(mask_soiltile)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: humrelv          !! humrel for each soil type
!$OMP THREADPRIVATE(humrelv)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: vegstressv       !! vegstress for each soil type  
!$OMP THREADPRIVATE(vegstressv)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:,:):: us               !! relative humidity 
!$OMP THREADPRIVATE(us)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: precisol         !! Eau tombee sur le sol
!$OMP THREADPRIVATE(precisol)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: precisol_ns      !! Eau tombee sur le sol par type de sol
!$OMP THREADPRIVATE(precisol_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ae_ns            !! Evaporation du sol nu par type de sol
!$OMP THREADPRIVATE(ae_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: evap_bare_lim_ns !! limitation of bare soil evaporation on each soil column 
                                                                         !! (used to deconvoluate vevapnu) 
!$OMP THREADPRIVATE(evap_bare_lim_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: free_drain_coef  !! Coefficient for free drainage at bottom
!$OMP THREADPRIVATE(free_drain_coef)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_bare_ns     !! evaporating bare soil fraction per tile
!$OMP THREADPRIVATE(frac_bare_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: rootsink         !! stress racinaire par niveau et type de sol
!$OMP THREADPRIVATE(rootsink)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: subsnowveg       !! Sublimation of snow on vegetation
!$OMP THREADPRIVATE(subsnowveg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: subsnownobio     !! Sublimation of snow on other surface types (ice, lakes, 
!$OMP THREADPRIVATE(subsnownobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: snowmelt          !! Quantite de glace fondue
!$OMP THREADPRIVATE(snowmelt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: icemelt          !! Quantite de glace fondue
!$OMP THREADPRIVATE(icemelt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: subsinksoil      !! Excess of sublimation as a sink for the soil
!$OMP THREADPRIVATE(subsinksoil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: vegtot           !! Total  vegetation (veget_max)
!$OMP THREADPRIVATE(vegtot)
  ! The last vegetation map which was used to distribute the reservoirs
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: resdist          !! Distribution of reservoirs
!$OMP THREADPRIVATE(resdist)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: mx_eau_var       !!
!$OMP THREADPRIVATE(mx_eau_var)

  ! arrays used by cwrr scheme
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: nroot            !! nvm * nstm * nslm
!$OMP THREADPRIVATE(nroot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: kfact_root       !! kjpindex * nslm * nstm
!$OMP THREADPRIVATE(kfact_root)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: kfact            !! nslm * nscm
!$OMP THREADPRIVATE(kfact)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: zz               !! nslm+1 * nstm
!$OMP THREADPRIVATE(zz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: dz               !! nslm+1 * nstm
!$OMP THREADPRIVATE(dz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: mc_lin           !! imin:imax * nscm
!$OMP THREADPRIVATE(mc_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: k_lin            !! imin:imax * nslm * nscm
!$OMP THREADPRIVATE(k_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: d_lin            !! imin:imax * nslm * nscm
!$OMP THREADPRIVATE(d_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: a_lin            !! imin:imax * nslm * nscm
!$OMP THREADPRIVATE(a_lin)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: b_lin            !! imin:imax * nslm * nscm
!$OMP THREADPRIVATE(b_lin)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: humtot           !! (:) Total Soil Moisture, Kg/m2
!$OMP THREADPRIVATE(humtot)
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION (:)          :: resolv           !! (:)
!$OMP THREADPRIVATE(resolv)

!! linarization coefficients of hydraulic conductivity K
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: k                !! (:,nslm)
!$OMP THREADPRIVATE(k)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: kk_moy           !! Mean hydraulic conductivity over soiltiles (mm/d)
!$OMP THREADPRIVATE(kk_moy)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: kk               !! Hydraulic conductivity for each soiltiles (mm/d)
!$OMP THREADPRIVATE(kk)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: a                !! (:,nslm)
!$OMP THREADPRIVATE(a)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: b                !!
!$OMP THREADPRIVATE(b)
!! linarization coefficients of hydraulic diffusivity D
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: d                !!
!$OMP THREADPRIVATE(d)
!! matrix coefficients
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: e                !!
!$OMP THREADPRIVATE(e)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: f                !!
!$OMP THREADPRIVATE(f)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: g1               !!
!$OMP THREADPRIVATE(g1)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ep               !!
!$OMP THREADPRIVATE(ep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: fp               !!
!$OMP THREADPRIVATE(fp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: gp               !!
!$OMP THREADPRIVATE(gp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: rhs              !!
!$OMP THREADPRIVATE(rhs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: srhs             !!
!$OMP THREADPRIVATE(srhs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: water2infilt     !! Water to be infiltrated
!$OMP THREADPRIVATE(water2infilt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc              !! (:,nstm) Total moisture content (mm)
!$OMP THREADPRIVATE(tmc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcr             !! (nstm) Total moisture constent at residual (mm)
!$OMP THREADPRIVATE(tmcr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmcs             !! (nstm) Total moisture constent at saturation (mm)
!$OMP THREADPRIVATE(tmcs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter       !! (:,nstm) Total moisture in the litter by soil type
!$OMP THREADPRIVATE(tmc_litter)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_mea     !! Total moisture in the litter over the grid
!$OMP THREADPRIVATE(tmc_litt_mea)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_wilt  !! (:,nstm) Moisture of litter at wilt pt
!$OMP THREADPRIVATE(tmc_litter_wilt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_field !! (:,nstm) Moisture of litter at field cap.
!$OMP THREADPRIVATE(tmc_litter_field)
!!! A CHANGER DANS TOUT HYDROL: tmc_litter_res et sat ne devraient pas dependre de ji - tdo
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_res   !! (:,nstm) Moisture of litter at residual moisture.
!$OMP THREADPRIVATE(tmc_litter_res)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_sat   !! (:,nstm) Moisture of litter at saturatiion
!$OMP THREADPRIVATE(tmc_litter_sat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_awet  !! (:,nstm) Moisture of litter at mc_awet
!$OMP THREADPRIVATE(tmc_litter_awet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tmc_litter_adry  !! (:,nstm) Moisture of litter at mc_dry
!$OMP THREADPRIVATE(tmc_litter_adry)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_wet_mea !! Total moisture in the litter over the grid below which 
                                                                         !! albedo is fixed 
!$OMP THREADPRIVATE(tmc_litt_wet_mea)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tmc_litt_dry_mea !! Total moisture in the litter over the grid above which 
                                                                         !! albedo is fixed 
!$OMP THREADPRIVATE(tmc_litt_dry_mea)
  LOGICAL, SAVE                                      :: tmc_init_updated = .FALSE. !! Flag allowing to determine if tmc is initialized.
!$OMP THREADPRIVATE(tmc_init_updated)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: v1               !! (:)
!$OMP THREADPRIVATE(v1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: qflux00          !! flux at the top of the soil column
!$OMP THREADPRIVATE(qflux00)

  !! par type de sol :
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ru_ns            !! ruissellement
!$OMP THREADPRIVATE(ru_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: dr_ns            !! drainage
!$OMP THREADPRIVATE(dr_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: tr_ns            !! transpiration 
!$OMP THREADPRIVATE(tr_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: cvs_over_veg     !! (:,nvm,nstm) old value of corr_veg_soil/veget_max kept 
                                                                         !! from diag to next split 
!$OMP THREADPRIVATE(cvs_over_veg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: corr_veg_soil    !! (:,nvm,nstm) percentage of each veg. type on each soil 
                                                                         !! of each grid point 
!$OMP THREADPRIVATE(corr_veg_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: mc               !! Total moisture Content (eg : liquid + frozen)
!$OMP THREADPRIVATE(mc)
   REAL(r_std),ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: mcl              !! Liquid moisture content
!$OMP THREADPRIVATE(mcl)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soilmoist        !! (:,nslm)
!$OMP THREADPRIVATE(soilmoist)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: soil_wet         !! soil wetness
!$OMP THREADPRIVATE(soil_wet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: soil_wet_litter  !! soil wetness of the litter
!$OMP THREADPRIVATE(soil_wet_litter)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: qflux            !! fluxes between the soil layers
!$OMP THREADPRIVATE(qflux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: tmat             !!
!$OMP THREADPRIVATE(tmat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: stmat            !!
!$OMP THREADPRIVATE(stmat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_hydro_diag       !! 
!$OMP THREADPRIVATE(frac_hydro_diag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: profil_froz_hydro     !! Frozen fraction for each hydrological soil layer
!$OMP THREADPRIVATE(profil_froz_hydro)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: profil_froz_hydro_ns  !! As  profil_froz_hydro per soiltile
!$OMP THREADPRIVATE(profil_froz_hydro_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: temp_hydro            !! Temp profile on hydrological levels
!$OMP THREADPRIVATE(temp_hydro)

!pss:+ TOPMODEL
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fsat             !! field capacity fraction
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwet             !! wetland fraction with WTD = 0 cm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwt1             !! wetland fraction with WTD entre 0 et -3cm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwt2             !! wetland fraction with WTD entre -3cm et -6cm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwt3             !! wetland fraction with WTD entre ... et ...
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: fwt4             !! etc.
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)      :: drunoff        !! runoff de Dunne
!pss:-

!pss:+ TOPMODEL parameter
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZMEAN            !! statistiques de la fonction de distribution des indices topo au sein de chaque maille
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZSTDT
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZSKEW
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZMIN
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZMAX
!!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: NB_PIXE
!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZWWILT           !! wilting point
!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZWSAT            !! saturation point
!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZWFC             !! field capacity
!  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZD_TOP           !! profondeur de sol pour TOPMODEL 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZM               !! parametre TOPMODEL
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: ZZPAS            !! pas des veceturs d indice topo au sein de chaque maille
! vecteurs calculees par TOPMODEL pour chaque maille (contenu = f(indice seuil); fsat = f(indice seuil); etc.)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ZTAB_FSAT        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ZTAB_WTOP
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ZTAB_FWET
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: ZTAB_WTOP_WET
!pss:-


CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_main
!!
!>\BRIEF         
!!
!! DESCRIPTION :
!! - called only one time for initialisation
!! - called every time step
!! - called one more time at last time step for writing _restart_ file
!!
!! - Choose between initialisation/Restart time/Time Step
!! - 1 Do initialisation ==> hydrol_init
!! - X if check_waterbal ==> hydrol_waterbal
!! - 2 prepares restart file for the next simulation
!! - 3 shared time step
!! - 3.1 computes snow  ==> hydrol_snow
!! - 3.2 computes vegetations reservoirs  ==> hydrol_vegupd
!! - 3.3 computes canopy  ==> hydrol_canop
!! - 3.4 computes surface reservoir  ==> hydrol_flood
!! - 3.5 computes soil hydrologie ==> hydrol_soil
!! - X if check_waterbal ==> hydrol_waterbal
!! - 4 write out file  ==> hydrol_alma/histwrite(*)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_main (kjit, kjpindex, & !pss+
       & lalo, resolution, & !pss-
       & dtradia, ldrestart_read, ldrestart_write, &
       & index, indexveg, indexsoil, indexlayer, indexnbdl, control_in, &
       & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max, njsc, &
       & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,  &
       & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, &
       & humrel, vegstress, drysoil_frac, evapot, evapot_penm, evap_bare_lim, flood_frac, flood_res, &
       & shumdiag,shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, reinf_slope, rest_id, hist_id, hist2_id,&
       & stempdiag, &
       & temp_air, pb, u, v, pgflux, &
       & snowrho,snowtemp,soiltemp,snowgrain,snowdz,snowheat,snowliq,&
       & grndflux,gtemp,gthick,gpkappa,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil, &
       & soilflxresid,snowflx,snowcap,pkappa_snow,lambda_snow,cgrnd_snow,dgrnd_snow,temp_sol_add, & !pss:+
       & fwet_out) !pss- !! Arsene 28-01-2016 - REMOVE drunoff_tot because never user and bug in sechiba_output.f90

!!!pss---note:here, 
! lalo and resolution are needed to extra_cti 
! fwet_out will be used if you wish to change the balance of energy on fractions wetland / saturated soil

    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
  
    ! input scalar 
    INTEGER(i_std), INTENT(in)                         :: kjit             !! Time step number 
    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
!pss:EXTRACT CTI resolution
!pss:+
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo       !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution !! size in x an y of the grid (m)
!pss:-
    INTEGER(i_std),INTENT (in)                         :: rest_id,hist_id  !! _Restart_ file and _history_ file identifier
    INTEGER(i_std),INTENT (in)                         :: hist2_id         !! _history_ file 2 identifier
    REAL(r_std), INTENT (in)                           :: dtradia          !! Time step in seconds
    LOGICAL, INTENT(in)                                :: ldrestart_read   !! Logical for _restart_ file to read
    LOGICAL, INTENT(in)                                :: ldrestart_write  !! Logical for _restart_ file to write
    ! input fields
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index            !! Indeces of the points on the map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg        !! Indeces of the points on the 3D map for veg
    INTEGER(i_std),DIMENSION (kjpindex*nstm), INTENT (in):: indexsoil      !! Indeces of the points on the 3D map for soil
    INTEGER(i_std),DIMENSION (kjpindex*nslm), INTENT (in):: indexlayer     !! Indeces of the points on the 3D map for soil layers
    INTEGER(i_std),DIMENSION (kjpindex*nslm), INTENT (in):: indexnbdl      !! Indeces of the points on the 3D map for of diagnostic soil layers
    TYPE(control_type), INTENT (in)                    :: control_in       !! Flags that (de)activate parts of the model
    !
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain      !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_snow      !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: returnflow       !! Routed water which comes back into the soil (from the 
                                                                           !! bottom) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: reinfiltration   !! Routed water which comes back into the soil (at the 
                                                                           !! top) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: irrigation       !! Water from irrigation returning to soil moisture  
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol_new     !! New soil temperature

    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio    !! Total fraction of ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: soilcap          !! Soil capacity
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile         !! Fraction of each soil tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vevapwet         !! Interception loss
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget            !! Fraction of vegetation type           
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max        !! Max. fraction of vegetation type (LAI -> infty)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintmax         !! Maximum water on vegetation for interception
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: transpir         !! Transpiration
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: reinf_slope      !! Slope coef
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot           !! Soil Potential Evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: evapot_penm      !! Soil Potential Evaporation Correction
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: flood_frac       !! flood fraction
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in) :: stempdiag        !! Diagnostic temp profile from thermosoil
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: temp_air         !! Air temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: u,v              !! Horizontal wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: pb               !! Surface pressure
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)       :: pgflux           !! Net energy into snowpack
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)       :: gtemp            !! First soil layer temperature
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)       :: gthick           !! First soil layer thickness
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)       :: gpkappa          !! First soil layer thermal conductivity
    REAL(r_std), DIMENSION (kjpindex),INTENT(inout)       :: soilflxresid     !! Energy flux to the snowpack
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(inout) :: pkappa_snow
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: humrel           !! Relative humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: vegstress        !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: drysoil_frac     !! function of litter wetness
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out):: shumdiag         !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out):: shumdiag_perma   !! Percent of porosity filled with water (mc/mcs) used for the thermal computations 
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: k_litt           !! litter approximate conductivity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: litterhumdiag    !! litter humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: tot_melt         !! Total melt    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out) :: qsintveg         !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)      :: evap_bare_lim    !! Limitation factor for bare soil evaporation
    
!pss:+
!!    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: drunoff_tot    !! Dunne runoff  !! Arsene 28-01-2016 - REMOVE because never user and bug in sechiba_output.f90
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: fwet_out   !! fwet: to change the balance of energy according to wetland fraction
!pss:-

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapnu          !! Bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapsno         !! Snow evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: vevapflo         !! Floodplain evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: flood_res        !! flood reservoir estimate
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: snow             !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)   :: snow_age         !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout) :: snow_nobio  !! Water balance on ice, lakes, .. [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout) :: snow_nobio_age !! Snow age on ice, lakes, ...
    !! We consider that any water on the ice is snow and we only peforme a water balance to have consistency.
    !! The water balance is limite to + or - 10^6 so that accumulation is not endless
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)       :: floodout     !! Flux out of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)       :: runoff       !! Complete runoff
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)       :: drainage     !! Drainage
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowrho      !! Snow density
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowtemp     !! Snow temperature
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(inout) :: soiltemp     !! Soil temperature
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowgrain    !! Snow grainsize
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowdz       !! Snow layer thickness
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowheat     !! Snow heat content
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(inout) :: snowliq      !! Snow liquid content (m)
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: grndflux     !! Net flux into soil W/m2
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)         :: zdz1_soil
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)         :: zdz2_soil
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)         :: cgrnd_soil
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)         :: dgrnd_soil
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout)   :: cgrnd_snow
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout)   :: dgrnd_snow
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)         :: snowflx
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)         :: snowcap
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)         :: lambda_snow
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)         :: temp_sol_add
    !! 0.4 Local variables

    INTEGER(i_std)                                     :: jst, jsl
    REAL(r_std),DIMENSION (kjpindex)                   :: soilwet          !! A temporary diagnostic of soil wetness
    REAL(r_std),DIMENSION (kjpindex)                   :: snowdepth        !! Depth of snow layer
    REAL(r_std),DIMENSION (kjpindex)                   :: njsc_tmp         !! Temporary REAL value for njsc to write it
    INTEGER(i_std), DIMENSION(kjpindex*imax)           :: mc_lin_axis_index
!pss:+
    logical                                           :: filealive, TOPMODEL_CTI
    INTEGER(i_std)                                    :: ind_spe, iet
!pss:-
!pss:+
    CHARACTER(LEN=80) :: filename   !! To store file names for I/O
    INTEGER(i_std) :: il, ils, ip, ix, iy, imin, jmin
    REAL(r_std) :: dlon, dlonmin, dlat, dlatmin
    INTEGER(i_std) :: iml, jml, lml, tml, fid
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: Zminf, Zmaxf, Zmeanf, Zstdf, Zskewf
    REAL(r_std),ALLOCATABLE,DIMENSION(:) :: lon_temp, lat_temp
    REAL(r_std) :: lev(1), pssdate, pssdt
    INTEGER(i_std) :: pssitau(1)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_rel, lon_rel
    INTEGER                  :: ALLOC_ERR

!pss-

   
!_ ================================================================================================================================

    !!_hydrol_main
    !
    !! Choose between initialisation/Restart time/Time Step
    !! 1 Do initialisation ==>hydrol_init
    !
    IF (l_first_hydrol) THEN

       IF (long_print) WRITE (numout,*) ' l_first_hydrol : call hydrol_init '

       CALL hydrol_init (kjit, ldrestart_read, kjpindex, index, rest_id, veget_max, soiltile, humrel,&
            & vegstress, snow, snow_age, snow_nobio, snow_nobio_age, qsintveg, &
            & grndflux,soiltemp,&
            & snowdz,snowgrain,snowrho,snowtemp,snowliq,snowheat,&
            & zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil,cgrnd_snow,dgrnd_snow,lambda_snow,snowflx,snowcap, & !pss+
            & fwet_out) !pss-
                 
       CALL hydrol_var_init (kjpindex, veget, veget_max, &
            & soiltile, njsc, mx_eau_var, shumdiag, shumdiag_perma, k_litt, &
            & litterhumdiag, drysoil_frac, evap_bare_lim,qsintveg) 

!pss+       
       TOPM_calcul  = .FALSE.
       CALL getin_p('TOPM_CALCUL', TOPM_calcul)

       IF (TOPM_calcul) THEN

          TOPMODEL_CTI = .TRUE.
          IF (TOPMODEL_CTI) THEN
            !  Needs to be a configurable variable
            !
            !
            !Config Key   = TOPMODEL_PARAMETERS_FILE
            !Config Desc  = Name of file from which TOPMODEL parameters file are read
            !Config Def   = TOPMODEL_param_1deg.nc
            !Config If    = NOT(IMPOSE_VEG)
            !Config Help  = The name of the file to be opened to read the TOPMODEL parameters. 
            !Config         
            !Config Units = [FILE]
            !
            filename = 'TOPMODEL_param_1deg.nc'
            CALL getin_p('TOPMODEL_PARAMETERS_FILE',filename)
            !
            IF (is_root_prc) THEN
               CALL flininfo(filename,iml, jml, lml, tml, fid)
               CALL flinclo(fid)
            ENDIF
            CALL bcast(iml)
            CALL bcast(jml)
            CALL bcast(lml)
            CALL bcast(tml)
            !
            ! soils_param.nc file is 1∞ soit texture file.
            !
            ALLOC_ERR=-1
            ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
            IF (ALLOC_ERR/=0) THEN
              WRITE(numout,*) "ERROR IN ALLOCATION of lat_rel : ",ALLOC_ERR
              STOP 
            ENDIF
            ALLOC_ERR=-1
            ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
            IF (ALLOC_ERR/=0) THEN
              WRITE(numout,*) "ERROR IN ALLOCATION of lon_rel : ",ALLOC_ERR
              STOP 
            ENDIF

            ALLOC_ERR=-1
            ALLOCATE(Zminf(iml,jml), STAT=ALLOC_ERR)
            IF (ALLOC_ERR/=0) THEN
              WRITE(numout,*) "ERROR IN ALLOCATION of ZMINf : ",ALLOC_ERR
              STOP 
            ENDIF
            ALLOC_ERR=-1
            ALLOCATE(Zmaxf(iml,jml), STAT=ALLOC_ERR)
            IF (ALLOC_ERR/=0) THEN
              WRITE(numout,*) "ERROR IN ALLOCATION of ZMAXf : ",ALLOC_ERR
              STOP 
            ENDIF
            ALLOC_ERR=-1
            ALLOCATE(Zmeanf(iml,jml), STAT=ALLOC_ERR)
            IF (ALLOC_ERR/=0) THEN
              WRITE(numout,*) "ERROR IN ALLOCATION of ZMEANf : ",ALLOC_ERR
              STOP 
            ENDIF
            ALLOC_ERR=-1
            ALLOCATE(Zstdf(iml,jml), STAT=ALLOC_ERR)
            IF (ALLOC_ERR/=0) THEN
              WRITE(numout,*) "ERROR IN ALLOCATION of ZSTDTf : ",ALLOC_ERR
              STOP 
            ENDIF
            ALLOC_ERR=-1
            ALLOCATE(Zskewf(iml,jml), STAT=ALLOC_ERR)
            IF (ALLOC_ERR/=0) THEN
              WRITE(numout,*) "ERROR IN ALLOCATION of ZSKEWf : ",ALLOC_ERR
              STOP 
            ENDIF

            ALLOC_ERR=-1
            ALLOCATE (lon_temp(iml),lat_temp(jml), STAT=ALLOC_ERR)
            IF (ALLOC_ERR/=0) THEN
              WRITE(numout,*) "ERROR IN ALLOCATION of lon_temp,lat_temp : ",ALLOC_ERR
              STOP 
            ENDIF
            !
            IF (is_root_prc) CALL flinopen(filename, .FALSE., iml, jml, lml, lon_rel,&
                                       & lat_rel, lev, tml, pssitau, pssdate, pssdt, fid)
            CALL bcast(lon_rel)
            CALL bcast(lat_rel)
            CALL bcast(pssitau)
            CALL bcast(pssdate)
            CALL bcast(pssdt)

            !
            IF (is_root_prc) CALL flinget(fid, 'Zmin', iml, jml, lml, tml, 1, 1, Zminf)
            IF (is_root_prc) CALL flinget(fid, 'Zmax', iml, jml, lml, tml, 1, 1, Zmaxf)
            IF (is_root_prc) CALL flinget(fid, 'Zmean', iml, jml, lml, tml, 1, 1, Zmeanf)
            IF (is_root_prc) CALL flinget(fid, 'Zstdev', iml, jml, lml, tml, 1, 1, Zstdf)
            IF (is_root_prc) CALL flinget(fid, 'Zskew', iml, jml, lml, tml, 1, 1, Zskewf)

            CALL bcast(Zminf)
            CALL bcast(Zmaxf)
            CALL bcast(Zmeanf)
            CALL bcast(Zstdf)
            CALL bcast(Zskewf)
            !
            IF (is_root_prc) CALL flinclo(fid)

        !!!! TOPMODEL parameters 2D into 1D 
               lon_temp(:) = lon_rel(:,1)
               lat_temp(:) = lat_rel(1,:)

               DO ip = 1, kjpindex
                  dlonmin = HUGE(1.)
                  DO ix = 1,iml
                     dlon = MIN( ABS(lalo(ip,2)-lon_temp(ix)), ABS(lalo(ip,2)+360.-lon_temp(ix)),&
                                          & ABS(lalo(ip,2)-360.-lon_temp(ix)) )
                     IF ( dlon .LT. dlonmin ) THEN
                        imin = ix
                        dlonmin = dlon
                     ENDIF
                  ENDDO
                  dlatmin = HUGE(1.)
                  DO iy = 1,jml
                     dlat = ABS(lalo(ip,1)-lat_temp(iy))
                     IF ( dlat .LT. dlatmin ) THEN
                        jmin = iy
                        dlatmin = dlat
                     ENDIF
                  ENDDO
                  ZMIN(ip) = Zminf(imin,jmin)
                  ZMAX(ip) = Zmaxf(imin,jmin)
                  ZMEAN(ip) = Zmeanf(imin,jmin)
                  ZSTDT(ip) = Zstdf(imin,jmin)
                  ZSKEW(ip) = Zskewf(imin,jmin)
               ENDDO

               DEALLOCATE (lon_temp)
               DEALLOCATE (lat_temp)
               DEALLOCATE (Zminf)
               DEALLOCATE (Zmaxf)
               DEALLOCATE (Zmeanf)
               DEALLOCATE (Zstdf)
               DEALLOCATE (Zskewf)
             
             TOPMODEL_CTI = .FALSE.
             write (numout,*) 'STATS CTI OK num1!'
             write (numout,*) 'psstest2'
          ELSE
             write (*,*) 'topmodel data already calculate!'
             write (numout,*) 'psstest3'
          ENDIF
       ELSE
          
          ZMIN(:)=0.
          ZMAX(:)=0.
          ZMEAN(:)=0.
          ZSTDT(:)=0.
          ZSKEW(:)=0.

       ENDIF

        !le deficit utilise pour TOPMODEL va etre calcule par rapport a la saturation
        !ZM(:)=(ZWFC(:)-ZWWILT(:))*ZD_TOP(:)/4.

        !ZM(:) = (mcs du grid_cell - mcw du grid_cell)*dpu_max/4.
        mcs_grid(:) = mcs(1)*soiltile(:,1)+mcs(2)*soiltile(:,2)+mcs(3)*soiltile(:,3)
        mcw_grid(:) = mcw(1)*soiltile(:,1)+mcw(2)*soiltile(:,2)+mcw(3)*soiltile(:,3)
        ZM(:) = ( mcs_grid(:) -  mcw_grid(:) )*dpu_max/4.


       IF(TOPM_calcul) THEN
  
          !2 obtention des differentes fonctions necessaires a TOPMODEL en chaque grid-cell  
          CALL init_top_main(kjpindex, lalo, veget, mcw_grid,mcs_grid,dpu_max, ZM,ZMIN, ZMAX, &
               & ZMEAN, ZSTDT, ZSKEW, ZTAB_FSAT, ZTAB_WTOP, ZTAB_FWET, ZTAB_WTOP_WET, ZZPAS)
          
       ELSE

          ZTAB_FSAT=0
          ZTAB_WTOP=0
          ZTAB_FWET=0
          ZTAB_WTOP_WET=0
          ZZPAS=0

       ENDIF

!pss:-

       ! If we check the water balance we first save the total amount of water
    !! X if check_waterbal ==> hydrol_waterbal
       IF (check_waterbal) THEN
          CALL hydrol_waterbal(kjpindex, index, .TRUE., dtradia, veget_max, &
               & totfrac_nobio, qsintveg, snow, snow_nobio,&
               & precip_rain, precip_snow, returnflow, reinfiltration, irrigation, tot_melt, &
               & vevapwet, transpir, vevapnu, vevapsno, vevapflo, floodout, runoff,drainage)
       ENDIF
       !
       IF (almaoutput) THEN
          CALL hydrol_alma(kjpindex, index, .TRUE., qsintveg, snow, snow_nobio, soilwet)
       ENDIF

       RETURN

    ENDIF

    !
    !! 2. prepares restart file for the next simulation
    !
    IF (ldrestart_write) THEN

       IF (long_print) WRITE (numout,*) ' we have to complete restart file with HYDROLOGIC variables '

       DO jst=1,nstm
          ! var_name= "mc_1" ... "mc_3"
          WRITE (var_name,"('moistc_',i1)") jst
          CALL restput_p(rest_id, var_name, nbp_glo,  nslm, 1, kjit, mc(:,:,jst), 'scatter',  nbp_glo, index_g)
       END DO
       !
       DO jst=1,nstm
          DO jsl=1,nslm
             ! var_name= "us_1_01" ... "us_3_11"
             WRITE (var_name,"('us_',i1,'_',i2.2)") jst,jsl
             CALL restput_p(rest_id, var_name, nbp_glo,nvm, 1,kjit,us(:,:,jst,jsl),'scatter',nbp_glo,index_g)
          END DO
       END DO
       !
       var_name= 'free_drain_coef'  
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  free_drain_coef, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'water2infilt'  
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  water2infilt, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'ae_ns'  
       CALL restput_p(rest_id, var_name, nbp_glo,   nstm, 1, kjit,  ae_ns, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'vegstress'
       CALL restput_p(rest_id, var_name, nbp_glo,   nvm, 1, kjit,  vegstress, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'snow'    
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  snow, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'snow_age'
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  snow_age, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'snow_nobio'    
       CALL restput_p(rest_id, var_name, nbp_glo,   nnobio, 1, kjit,  snow_nobio, 'scatter', nbp_glo, index_g)
       !
       var_name= 'snow_nobio_age'
       CALL restput_p(rest_id, var_name, nbp_glo,   nnobio, 1, kjit,  snow_nobio_age, 'scatter', nbp_glo, index_g)
       !
       var_name= 'qsintveg'
       CALL restput_p(rest_id, var_name, nbp_glo, nvm, 1, kjit,  qsintveg, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'resdist'
       CALL restput_p(rest_id, var_name, nbp_glo, nstm, 1, kjit,  resdist, 'scatter',  nbp_glo, index_g)       
       !
!pss:+
       !
        var_name= 'fwet_out'
        CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  fwet_out, 'scatter',  nbp_glo, index_g)
       !
!pss:-

       DO jst=1,nstm
          ! var_name= "cvs_over_veg_1" ... "cvs_over_veg_3"
          WRITE (var_name,"('cvs_over_veg_',i1)") jst
          CALL restput_p(rest_id, var_name, nbp_glo,  nvm, 1, kjit, cvs_over_veg(:,:,jst), 'scatter',  nbp_glo, index_g)
       END DO
       !
       IF ( check_waterbal ) THEN
          var_name= 'tot_water_beg'
          CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  tot_water_end, 'scatter', nbp_glo, index_g)
       ENDIF

       IF (ok_explicitsnow) THEN
          var_name= 'grndflux'
          CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  grndflux,'scatter',  nbp_glo, index_g)
          
          var_name = 'snowrho'
          CALL restput_p(rest_id, var_name, nbp_glo, nsnow, 1, kjit, snowrho, 'scatter', nbp_glo, index_g)
          
          var_name = 'snowheat'
          CALL restput_p(rest_id, var_name, nbp_glo, nsnow, 1, kjit, snowheat, 'scatter', nbp_glo, index_g)
          
          var_name = 'snowliq'
          CALL restput_p(rest_id, var_name, nbp_glo, nsnow, 1, kjit, snowliq, 'scatter', nbp_glo, index_g)
          
          var_name = 'snowtemp'
          CALL restput_p(rest_id, var_name, nbp_glo, nsnow, 1, kjit, snowtemp, 'scatter', nbp_glo, index_g)
          
          var_name = 'soiltemp'
          CALL restput_p(rest_id, var_name, nbp_glo, ngrnd, 1, kjit, soiltemp, 'scatter', nbp_glo, index_g)
          
          var_name = 'snowdz'
          CALL restput_p(rest_id, var_name, nbp_glo, nsnow, 1, kjit, snowdz, 'scatter', nbp_glo, index_g)

          var_name= 'snowflx'
          CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  snowflx, 'scatter',  nbp_glo, index_g)

          var_name= 'snowcap'
          CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  snowcap, 'scatter',  nbp_glo, index_g)

          var_name = 'cgrnd_snow'
          CALL restput_p(rest_id, var_name, nbp_glo, nsnow, 1, kjit, cgrnd_snow, 'scatter', nbp_glo, index_g)

          var_name = 'dgrnd_snow'
          CALL restput_p(rest_id, var_name, nbp_glo, nsnow, 1, kjit, dgrnd_snow, 'scatter', nbp_glo, index_g)

          var_name = 'zdz1_soil'
          CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, zdz1_soil, 'scatter', nbp_glo, index_g)

          var_name = 'zdz2_soil'
          CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, zdz2_soil, 'scatter', nbp_glo, index_g)

          var_name = 'cgrnd_soil'
          CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, cgrnd_soil, 'scatter', nbp_glo, index_g)

          var_name = 'dgrnd_soil'
          CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, dgrnd_soil, 'scatter', nbp_glo, index_g)

       END IF

       RETURN

    END IF

    !
    !! 3. Shared time step
    IF (long_print) WRITE (numout,*) 'hydrol pas de temps = ',dtradia


!pss:+ 
    !! 3.0 Calculate wetland fractions
        
    IF (TOPM_calcul) THEN
        CALL hydro_subgrid_main(kjpindex, ZTAB_FSAT, ZTAB_WTOP, humtot, profil_froz_hydro, fsat,&
           & ZTAB_FWET,ZTAB_WTOP_WET,fwet, dpu_max, &
           & 1000*(mcs_grid(:)-mcw_grid(:)), fwt1, fwt2, fwt3, fwt4, ZM, ZMIN, ZMAX, ZZPAS)

    ELSE
        fsat(:)=0.0
        fwet(:)=0.0
        fwt1(:)=0.0
        fwt2(:)=0.0
        fwt3(:)=0.0
        fwt4(:)=0.0
    ENDIF

    fwet_out(:)=fwet(:)
!pss:-

    ! 
    !! 3.1 Calculate snow processes with explicit method or bucket snow model
    IF (ok_explicitsnow) THEN
       ! Explicit snow model
       IF (long_print) WRITE (numout,*) ' ok_explicitsnow : use three-snow layer '
       
       CALL explicitsnow_main(kjpindex,dtradia,precip_rain,precip_snow,temp_air,pb,u,v,temp_sol_new,soilcap,&
             & pgflux,frac_nobio,totfrac_nobio,&
             & gtemp,gthick,gpkappa,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil,vevapsno, &
             & snow_age,snow_nobio_age,snow_nobio,snowrho,snowgrain,snowdz,snowtemp,snowheat,snowliq,&
             & snow,subsnownobio,grndflux,snowmelt,tot_melt,soilflxresid,subsinksoil,snowflx,snowcap,&
             & pkappa_snow,lambda_snow,cgrnd_snow,dgrnd_snow,temp_sol_add, veget_max)  !! Arsene 04-03-2015 Add for snowcompact
 
    ELSE
       ! Bucket snow model
       CALL hydrol_snow(kjpindex, dtradia, precip_rain, precip_snow, temp_sol_new, soilcap, &
            frac_nobio, totfrac_nobio, vevapnu, vevapsno, snow, snow_age, snow_nobio, snow_nobio_age, &
            tot_melt, snowdepth,snowmelt)
    END IF
        
    !
    !! 3.2 computes vegetations reservoirs  ==>hydrol_vegupd
    CALL hydrol_vegupd(kjpindex, veget, veget_max, soiltile, qsintveg,resdist)

    !
    !! 3.3 computes canopy  ==>hydrol_canop
    CALL hydrol_canop(kjpindex, precip_rain, vevapwet, veget_max, veget, qsintmax, qsintveg,precisol,tot_melt)

    !
    !! 3.4 computes surface reservoir  ==>hydrol_flood
    CALL hydrol_flood(kjpindex, dtradia, vevapflo, flood_frac, flood_res, floodout)

    !
    !! 3.5 computes soil hydrologie ==>hydrol_soil

    CALL hydrol_soil(kjpindex, dtradia, veget_max, soiltile, njsc, reinf_slope,  &
         transpir, vevapnu, evapot, evapot_penm, runoff, &
         drainage, returnflow, reinfiltration, irrigation, &
         tot_melt,evap_bare_lim, shumdiag, shumdiag_perma, &
         k_litt, litterhumdiag, humrel, vegstress, drysoil_frac,&
         stempdiag,snow,snowdz, & !pss:+
         fsat) !pss- !! Arsene 28-01-2016 - REMOVE drunoff_tot because never user and bug in sechiba_output.f90

    ! If we check the water balance we end with the comparison of total water change and fluxes
    !! X if check_waterbal ==> hydrol_waterbal
    IF (check_waterbal) THEN
       CALL hydrol_waterbal(kjpindex, index, .FALSE., dtradia, veget_max, totfrac_nobio, &
            & qsintveg, snow,snow_nobio, precip_rain, precip_snow, returnflow, reinfiltration, &
            & irrigation, tot_melt, vevapwet, transpir, vevapnu, vevapsno, vevapflo, floodout, runoff, drainage)
    ENDIF

    !! 4 write out file  ==> hydrol_alma/histwrite(*)
    !
    ! If we use the ALMA standards
    IF (almaoutput) THEN
       CALL hydrol_alma(kjpindex, index, .FALSE., qsintveg, snow, snow_nobio, soilwet)
    ENDIF


    DO jst = 1, nstm
       ! var_name= "mc_1" ... "mc_3"
       WRITE (var_name,"('moistc_',i1)") jst
       CALL xios_orchidee_send_field(TRIM(var_name),mc(:,:,jst))

       ! var_name= "kfactroot_1" ... "kfactroot_3"
       WRITE (var_name,"('kfactroot_',i1)") jst
       CALL xios_orchidee_send_field(TRIM(var_name),kfact_root(:,:,jst))

       ! var_name= "vegetsoil_1" ... "vegetsoil_3"
       WRITE (var_name,"('vegetsoil_',i1)") jst
       CALL xios_orchidee_send_field(TRIM(var_name),corr_veg_soil(:,:,jst))
    END DO

    CALL xios_orchidee_send_field("evapnu_soil",ae_ns*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("drainage_soil",dr_ns*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("transpir_soil",tr_ns*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("runoff_soil",ru_ns*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("humtot_soil",tmc)
    CALL xios_orchidee_send_field("humtot",humtot)
    njsc_tmp(:)=njsc(:)
    CALL xios_orchidee_send_field("soilindex",njsc_tmp)
    CALL xios_orchidee_send_field("humrel",humrel)     
    CALL xios_orchidee_send_field("drainage",drainage*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("runoff",runoff*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("precisol",precisol*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("rain",precip_rain*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("rain_alma",precip_rain*0.0005555556)
    CALL xios_orchidee_send_field("snowf",precip_snow)
    CALL xios_orchidee_send_field("snowf_alma",precip_snow*0.0005555556)
    CALL xios_orchidee_send_field("qsintmax",qsintmax)
    CALL xios_orchidee_send_field("qsintveg",qsintveg)
    CALL xios_orchidee_send_field("SWI",swi)

    IF ( control_in%do_floodplains ) THEN
       CALL xios_orchidee_send_field("floodout",floodout*one_day/dt_sechiba)
    END IF

    IF (check_waterbal) THEN
       CALL xios_orchidee_send_field("TotWater",tot_water_end)
       CALL xios_orchidee_send_field("TotWaterFlux",tot_flux*one_day/dt_sechiba)
    END IF

    CALL xios_orchidee_send_field("Qs",runoff/dt_sechiba)
    CALL xios_orchidee_send_field("Qsb",drainage/dt_sechiba)
    CALL xios_orchidee_send_field("Qsm",snowmelt/dt_sechiba)
    CALL xios_orchidee_send_field("SoilMoist",soilmoist)
    CALL xios_orchidee_send_field("SubSnow",vevapsno/dt_sechiba)
    CALL xios_orchidee_send_field("SnowDepth",snowdepth)

    IF (almaoutput) THEN
       CALL xios_orchidee_send_field("SoilWet",soilwet)
       CALL xios_orchidee_send_field("RootMoist",tot_watsoil_end)
       CALL xios_orchidee_send_field("DelSoilMoist",delsoilmoist)
       CALL xios_orchidee_send_field("DelSWE",delswe)
       CALL xios_orchidee_send_field("DelIntercept",delintercept)  
    END IF


    IF ( .NOT. almaoutput ) THEN
       DO jst=1,nstm
          ! var_name= "mc_1" ... "mc_3"
          WRITE (var_name,"('moistc_',i1)") jst
          CALL histwrite_p(hist_id, TRIM(var_name), kjit,mc(:,:,jst), kjpindex*nslm, indexlayer)

          ! var_name= "kfactroot_1" ... "kfactroot_3"
          WRITE (var_name,"('kfactroot_',i1)") jst
          CALL histwrite_p(hist_id, TRIM(var_name), kjit, kfact_root(:,:,jst), kjpindex*nslm, indexlayer)

          ! var_name= "vegetsoil_1" ... "vegetsoil_3"
          WRITE (var_name,"('vegetsoil_',i1)") jst
          CALL histwrite_p(hist_id, TRIM(var_name), kjit,corr_veg_soil(:,:,jst), kjpindex*nvm, indexveg)
       ENDDO
       CALL histwrite_p(hist_id, 'evapnu_soil', kjit, ae_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'drainage_soil', kjit, dr_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'transpir_soil', kjit, tr_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'runoff_soil', kjit, ru_ns, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'humtot_soil', kjit, tmc, kjpindex*nstm, indexsoil)
       CALL histwrite_p(hist_id, 'humtot', kjit, humtot, kjpindex, index)
       njsc_tmp(:)=njsc(:)
       CALL histwrite_p(hist_id, 'soilindex', kjit, njsc_tmp, kjpindex, index)
       CALL histwrite_p(hist_id, 'humrel',   kjit, humrel,   kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'drainage', kjit, drainage, kjpindex, index)
       CALL histwrite_p(hist_id, 'runoff', kjit, runoff, kjpindex, index)
       CALL histwrite_p(hist_id, 'precisol', kjit, precisol, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'rain', kjit, precip_rain, kjpindex, index)
       CALL histwrite_p(hist_id, 'snowf', kjit, precip_snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'qsintmax', kjit, qsintmax, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'qsintveg', kjit, qsintveg, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'SWI', kjit, swi, kjpindex, index)
       CALL histwrite_p(hist_id, 'snowmelt',kjit,snowmelt,kjpindex,index)
       CALL histwrite_p(hist_id, 'shumdiag_perma',kjit,shumdiag_perma,kjpindex*nbdl,indexnbdl)

!pss:+ ! write out wetland fraction and CTI parameters
       CALL histwrite_p(hist_id, 'fsat', kjit, fsat, kjpindex, index)
       CALL histwrite_p(hist_id, 'fwet', kjit, fwet, kjpindex, index)
       CALL histwrite_p(hist_id, 'fwt1', kjit, fwt1, kjpindex, index)
       CALL histwrite_p(hist_id, 'fwt2', kjit, fwt2, kjpindex, index)
       CALL histwrite_p(hist_id, 'fwt3', kjit, fwt3, kjpindex, index)
       CALL histwrite_p(hist_id, 'fwt4', kjit, fwt4, kjpindex, index)
       CALL histwrite_p(hist_id, 'ZMIN', kjit, ZMIN, kjpindex, index)
       CALL histwrite_p(hist_id, 'ZMAX', kjit, ZMAX, kjpindex, index)
       CALL histwrite_p(hist_id, 'ZMEAN', kjit, ZMEAN, kjpindex, index)
       !CALL histwrite_p(hist_id, 'NB_PIXE', kjit, NB_PIXE, kjpindex, index)
       CALL histwrite_p(hist_id, 'ZSTDT', kjit, ZSTDT, kjpindex, index)
       CALL histwrite_p(hist_id, 'ZSKEW', kjit, ZSKEW, kjpindex, index)
!       CALL histwrite_p(hist_id, 'dsg', kjit, dsg, kjpindex*nvm, indexveg)
!       CALL histwrite_p(hist_id, 'dsp', kjit, dsp, kjpindex*nvm, indexveg)
!       CALL histwrite_p(hist_id, 'ZWSAT', kjit, ZWSAT, kjpindex, index)
!       CALL histwrite_p(hist_id, 'ZWWILT', kjit, ZWWILT, kjpindex, index) 
!       CALL histwrite_p(hist_id, 'ZWFC', kjit, ZWFC, kjpindex, index) 
!       CALL histwrite_p(hist_id, 'RU', kjit, ruu_ch, kjpindex, index) 
!       CALL histwrite_p(hist_id, 'mx_eau_var', kjit, mx_eau_var, kjpindex, index)
!       CALL histwrite_p(hist_id, 'drunoff_tot', kjit, drunoff_tot, kjpindex, index) !! Arsene 28-01-2016 - REMOVE because never user and bug in sechiba_output.f90
!pss:-

       IF ( control_in%do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'floodout', kjit, floodout, kjpindex, index)
       ENDIF
       !
       IF ( hist2_id > 0 ) THEN
          DO jst=1,nstm
             ! var_name= "mc_1" ... "mc_3"
             WRITE (var_name,"('moistc_',i1)") jst
             CALL histwrite_p(hist2_id, TRIM(var_name), kjit,mc(:,:,jst), kjpindex*nslm, indexlayer)

             ! var_name= "kfactroot_1" ... "kfactroot_3"
             WRITE (var_name,"('kfactroot_',i1)") jst
             CALL histwrite_p(hist2_id, TRIM(var_name), kjit, kfact_root(:,:,jst), kjpindex*nslm, indexlayer)

             ! var_name= "vegetsoil_1" ... "vegetsoil_3"
             WRITE (var_name,"('vegetsoil_',i1)") jst
             CALL histwrite_p(hist2_id, TRIM(var_name), kjit,corr_veg_soil(:,:,jst), kjpindex*nvm, indexveg)
          ENDDO
          CALL histwrite_p(hist2_id, 'evapnu_soil', kjit, ae_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'drainage_soil', kjit, dr_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'transpir_soil', kjit, tr_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'runoff_soil', kjit, ru_ns, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'humtot_soil', kjit, tmc, kjpindex*nstm, indexsoil)
          CALL histwrite_p(hist2_id, 'humtot', kjit, humtot, kjpindex, index)
          njsc_tmp(:)=njsc(:)
          CALL histwrite_p(hist2_id, 'soilindex', kjit, njsc_tmp, kjpindex, index)
          CALL histwrite_p(hist2_id, 'humrel',   kjit, humrel,   kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'drainage', kjit, drainage, kjpindex, index)
          CALL histwrite_p(hist2_id, 'runoff', kjit, runoff, kjpindex, index)
          IF ( control_in%do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'floodout', kjit, floodout, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'precisol', kjit, precisol, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'rain', kjit, precip_rain, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowf', kjit, precip_snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowmelt',kjit,snowmelt,kjpindex,index)
          CALL histwrite_p(hist2_id, 'qsintmax', kjit, qsintmax, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'qsintveg', kjit, qsintveg, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'SWI', kjit, swi, kjpindex, index) 
          !
          IF (check_waterbal) THEN
             CALL histwrite_p(hist2_id, 'TotWater', kjit, tot_water_end, kjpindex, index)
             CALL histwrite_p(hist2_id, 'TotWaterFlux', kjit, tot_flux, kjpindex, index)
          ENDIF
       ENDIF
    ELSE
       CALL histwrite_p(hist_id, 'Snowf', kjit, precip_snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'Rainf', kjit, precip_rain, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qs', kjit, runoff, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qsb', kjit, drainage, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qsm', kjit, snowmelt, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelSoilMoist', kjit, delsoilmoist, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelSWE', kjit, delswe, kjpindex, index)
       CALL histwrite_p(hist_id, 'DelIntercept', kjit, delintercept, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'SoilMoist', kjit, soilmoist, kjpindex*nslm, indexlayer)
       CALL histwrite_p(hist_id, 'SoilWet', kjit, soilwet, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'RootMoist', kjit, tot_watsoil_end, kjpindex, index)
       CALL histwrite_p(hist_id, 'SubSnow', kjit, vevapsno, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'SnowDepth', kjit, snowdepth, kjpindex, index)
       !
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'Snowf', kjit, precip_snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Rainf', kjit, precip_rain, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qs', kjit, runoff, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qsb', kjit, drainage, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qsm', kjit, snowmelt, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSoilMoist', kjit, delsoilmoist, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSWE', kjit, delswe, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelIntercept', kjit, delintercept, kjpindex, index)
          !
          CALL histwrite_p(hist2_id, 'SoilMoist', kjit, soilmoist, kjpindex*nslm, indexlayer)
          CALL histwrite_p(hist2_id, 'SoilWet', kjit, soilwet, kjpindex, index)
          !
          CALL histwrite_p(hist2_id, 'RootMoist', kjit, tot_watsoil_end, kjpindex, index)
          CALL histwrite_p(hist2_id, 'SubSnow', kjit, vevapsno, kjpindex, index)
          !
          CALL histwrite_p(hist2_id, 'SnowDepth', kjit, snowdepth, kjpindex, index)
       ENDIF
    ENDIF

    IF (ok_freeze_cwrr) THEN
       CALL histwrite_p(hist_id, 'profil_froz_hydro', kjit,profil_froz_hydro , kjpindex*nslm, indexlayer)
       DO jst=1,nstm
          WRITE (var_name,"('profil_froz_hydro_',i1)") jst
          CALL histwrite_p(hist_id, TRIM(var_name), kjit, profil_froz_hydro_ns(:,:,jst), kjpindex*nslm, indexlayer)
       ENDDO
       CALL histwrite_p(hist_id, 'temp_hydro', kjit,temp_hydro , kjpindex*nslm, indexlayer)
       CALL histwrite_p(hist_id, 'kk_moy', kjit, kk_moy,kjpindex*nslm, indexlayer) ! averaged over soiltiles
    ENDIF

    IF (l_second_hydrol) THEN 
       l_second_hydrol=.FALSE.
    ENDIF
    IF (long_print) WRITE (numout,*) ' hydrol_main Done '

  END SUBROUTINE hydrol_main


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_init
!!
!>\BRIEF        Initializations and memory allocation   
!!
!! DESCRIPTION  :
!! - 1 Some initializations
!! - 2 make dynamic allocation with good dimension
!! - 2.1 array allocation for soil textur
!! - 2.2 Soil texture choice
!! - 3 Other array allocation
!! - 4 Open restart input file and read data for HYDROLOGIC process
!! - 5 get restart values if none were found in the restart file
!! - 6 Vegetation array      
!! - 7 set humrelv from us
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!!_ hydrol_init

  SUBROUTINE hydrol_init(kjit, ldrestart_read, kjpindex, index, rest_id, veget_max, soiltile, humrel,&
       vegstress, snow, snow_age, snow_nobio, snow_nobio_age, qsintveg, &
       grndflux,soiltemp,&
       snowdz,snowgrain,snowrho,snowtemp,snowliq,snowheat,&
       zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil,cgrnd_snow,dgrnd_snow,lambda_snow,snowflx,snowcap, & !pss+
       fwet_out) !pss-


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number 
    LOGICAL,INTENT (in)                                 :: ldrestart_read     !! Logical for _restart_ file to read
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! _Restart_ file identifier 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max          !! Carte de vegetation max
    REAL(r_std),DIMENSION (kjpindex,nstm), INTENT (in)  :: soiltile           !! Fraction of each soil tile (0-1, unitless)

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: humrel             !! Stress hydrique, relative humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: vegstress          !! Veg. moisture stress (only for vegetation growth)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: snow               !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: snow_age           !! Snow age
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: snow_nobio       !! Snow on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: snow_nobio_age   !! Snow age on ice, lakes, ...
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: qsintveg           !! Water on vegetation due to interception
!pss:+
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)  :: fwet_out            !! output wetland fraction to change energy or runoff ???!!!
!pss:-

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)       :: grndflux          !! Net flux into soil W/m2
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout) :: snowdz            !! Snow depth
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout) :: snowgrain         !! Snow grain size
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout) :: snowheat          !! Snow heat content
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout) :: snowtemp          !! Snow temperature
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout) :: snowliq           !! Liquid water content
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout) :: snowrho           !! Snow density
    REAL(r_std),DIMENSION (kjpindex,ngrnd),INTENT(inout) :: soiltemp          !! Soil temperature
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)       :: zdz1_soil
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)       :: zdz2_soil
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)       :: cgrnd_soil
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)       :: dgrnd_soil
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout) :: cgrnd_snow
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(inout) :: dgrnd_snow
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)       :: snowflx
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)       :: snowcap
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)       :: lambda_snow

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ier                   !! Error code
    INTEGER(i_std)                                     :: ji, jv, jst, jsl, jsc !! Indices
    INTEGER(i_std), PARAMETER                          :: error_level = 3       !! Error level for consistency check
                                                                                !! Switch to 2 tu turn fatal errors into warnings  

!_ ================================================================================================================================

    ! initialisation
    IF (l_first_hydrol) THEN 
       l_first_hydrol=.FALSE.
    ELSE 
       WRITE (numout,*) ' l_first_hydrol false . we stop '
       STOP 'hydrol_init'
    ENDIF

    !
    !! 1 Some initializations
    !
    !
    !Config Key   = CHECK_CWRR
    !Config Desc  = Should we check detailed CWRR water balance ?
    !Config Def   = n
    !Config If    = HYDROL_CWRR
    !Config Help  = This parameters allows the user to check
    !Config         the detailed water balance in each time step
    !Config         of CWRR.
    !Config Units = [FLAG]
    !
    check_cwrr = .FALSE.
    CALL getin_p('CHECK_CWRR', check_cwrr)
    !
    !Config Key   = DO_PONDS
    !Config Desc  = Should we include ponds 
    !Config Def   = n
    !Config If    = HYDROL_CWRR
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the ponds and return 
    !Config         the water into the soil moisture. If this is 
    !Config         activated, then there is no reinfiltration 
    !Config         computed inside the hydrol module.
    !Config Units = [FLAG]
    !
    doponds = .FALSE.
    CALL getin_p('DO_PONDS', doponds)


    !! 2 make dynamic allocation with good dimension

    !! 2.1 array allocation for soil texture

    ALLOCATE (nvan(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in nvan allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (avan(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in avan allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcr(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcr allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcs(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcs allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ks(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ks allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (pcent(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in pcent allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (free_drain_max(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in free_drain_max allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcf(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcf allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcw(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcw allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF
    
    ALLOCATE (mc_awet(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc_awet allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mc_adry(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc_adry allocation. We stop. We need nscm words = ',nscm
       STOP 'hydrol_init'
    END IF
       
    !!__2.2 Soil texture choose

    SELECTCASE (nscm)
    CASE (3)
              
       nvan(:) = nvan_fao(:)       
       avan(:) = avan_fao(:)
       mcr(:) = mcr_fao(:)
       mcs(:) = mcs_fao(:)
       ks(:) = ks_fao(:)
       pcent(:) = pcent_fao(:)
       free_drain_max(:) = free_drain_max_fao(:)
       mcf(:) = mcf_fao(:)
       mcw(:) = mcw_fao(:)
       mc_awet(:) = mc_awet_fao(:)
       mc_adry(:) = mc_adry_fao(:)
    CASE (12)
       
       nvan(:) = nvan_usda(:)
       avan(:) = avan_usda(:)
       mcr(:) = mcr_usda(:)
       mcs(:) = mcs_usda(:)
       ks(:) = ks_usda(:)
       pcent(:) = pcent_usda(:)
       free_drain_max(:) = free_drain_max_usda(:)
       mcf(:) = mcf_usda(:)
       mcw(:) = mcw_usda(:)
       mc_awet(:) = mc_awet_usda(:)
       mc_adry(:) = mc_adry_usda(:)
       
    CASE DEFAULT
       WRITE (numout,*) 'Unsupported soil type classification. Choose between zobler, fao and usda according to the map'
       STOP 'hydrol_init'
    ENDSELECT


    !! 2.3 Read in the run.def the parameters values defined by the user

    !Config Key   = CWRR_N_VANGENUCHTEN
    !Config Desc  = Van genuchten coefficient n
    !Config If    = HYDROL_CWRR
    !Config Def   = 1.89, 1.56, 1.31
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [-]
    CALL getin_p("CWRR_N_VANGENUCHTEN",nvan)

    !! Check parameter value (correct range)
    IF ( ANY(nvan(:) <= zero) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for CWRR_N_VANGENUCHTEN.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_A_VANGENUCHTEN
    !Config Desc  = Van genuchten coefficient a
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.0075, 0.0036, 0.0019
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [1/mm]  
    CALL getin_p("CWRR_A_VANGENUCHTEN",avan)

    !! Check parameter value (correct range)
    IF ( ANY(avan(:) <= zero) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for CWRR_A_VANGENUCHTEN.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_RESIDUAL
    !Config Desc  = Residual soil water content
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.065, 0.078, 0.095
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [mm]  
    CALL getin_p("VWC_RESIDUAL",mcr)

    !! Check parameter value (correct range)
    IF ( ANY(mcr(:) < zero) .OR. ANY(mcr(:) > 1.)  ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_RESIDUAL.", &
            &     "This parameter is ranged between 0 and 1. ", &
            &     "Please, check parameter value in run.def. ")
    END IF

    
    !Config Key   = VWC_SAT
    !Config Desc  = Saturated soil water content
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.41, 0.43, 0.41
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [-]  
    CALL getin_p("VWC_SAT",mcs)

    !! Check parameter value (correct range)
    IF ( ANY(mcs(:) < zero) .OR. ANY(mcs(:) > 1.) .OR. ANY(mcs(:) <= mcr(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_SAT.", &
            &     "This parameter should be greater than VWC_RESIDUAL and less than 1. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_KS 
    !Config Desc  = Hydraulic conductivity Saturation
    !Config If    = HYDROL_CWRR 
    !Config Def   = 1060.8, 249.6, 62.4
    !Config Help  = This parameter will be constant over the entire 
    !Config         simulated domain, thus independent from soil
    !Config         texture.   
    !Config Units = [mm/d]   
    CALL getin_p("CWRR_KS",ks)

    !! Check parameter value (correct range)
    IF ( ANY(ks(:) <= zero) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for CWRR_KS.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = WETNESS_TRANSPIR_MAX
    !Config Desc  = Soil moisture above which transpir is max
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.5, 0.5, 0.5
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [-]    
    CALL getin_p("WETNESS_TRANSPIR_MAX",pcent)

    !! Check parameter value (correct range)
    IF ( ANY(pcent(:) <= zero) .OR. ANY(pcent(:) > 1.) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for WETNESS_TRANSPIR_MAX.", &
            &     "This parameter should be positive and less or equals than 1. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = DRAINAGE_FACTOR_F
    !Config Desc  = Max value of the permeability coeff at the bottom of the soil 
    !Config If    = HYDROL_CWRR
    !Config Def   = 1.0, 1.0, 1.0 
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [-]   
    CALL getin_p("DRAINAGE_FACTOR_F",free_drain_max)

    !! Check parameter value (correct range)
    IF ( ANY(free_drain_max(:) < zero) .OR. ANY(free_drain_max(:) > 1.) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for DRAINAGE_FACTOR_F.", &
            &     "This parameter should be positive and not greater than 1. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_FC 
    !Config Desc  = Volumetric water content field capacity
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.32, 0.32, 0.32
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [-]   
    CALL getin_p("VWC_FC",mcf)

    !! Check parameter value (correct range)
    IF ( ANY(mcf(:) > mcs(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_FC.", &
            &     "This parameter should be less than VWC_SAT. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_WP
    !Config Desc  = Volumetric water content Wilting pt
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.10, 0.10, 0.10 
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [-]   
    CALL getin_p("VWC_WP",mcw)

    !! Check parameter value (correct range)
    IF ( ANY(mcw(:) > mcf(:)) .OR. ANY(mcw(:) < mcr(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_WP.", &
            &     "This parameter should be greater or equal than VWC_RESIDUAL and less or equal than VWC_SAT.", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_MIN_FOR_WET_ALB
    !Config Desc  = Vol. wat. cont. above which albedo is cst
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.25, 0.25, 0.25
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [-]  
    CALL getin_p("VWC_MIN_FOR_WET_ALB",mc_awet)

    !! Check parameter value (correct range)
    IF ( ANY(mc_awet(:) < 0) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_MIN_FOR_WET_ALB.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = VWC_MAX_FOR_DRY_ALB
    !Config Desc  = Vol. wat. cont. below which albedo is cst
    !Config If    = HYDROL_CWRR
    !Config Def   = 0.1, 0.1, 0.1
    !Config Help  = This parameter is independent from soil texture for
    !Config         the time being.
    !Config Units = [-]   
    CALL getin_p("VWC_MAX_FOR_DRY_ALB",mc_adry)

    !! Check parameter value (correct range)
    IF ( ANY(mc_adry(:) < 0) .OR. ANY(mc_adry(:) > mc_awet(:)) ) THEN
       CALL ipslerr_p(error_level, "hydrol_init.", &
            &     "Wrong parameter value for VWC_MAX_FOR_DRY_ALB.", &
            &     "This parameter should be positive and not greater than VWC_MIN_FOR_WET_ALB.", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !! 3 Other array allocation


    ALLOCATE (mask_veget(kjpindex,nvm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mask_veget allocation. We stop. We need kjpindex*nvm words = ',kjpindex*nvm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mask_soiltile(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mask_soiltile allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (humrelv(kjpindex,nvm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in humrelv allocation. We stop. We need kjpindex words = ',kjpindex*nvm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (vegstressv(kjpindex,nvm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vegstressv allocation. We stop. We need kjpindex words = ',kjpindex*nvm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (us(kjpindex,nvm,nstm,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in us allocation. We stop. We need kjpindex words = ',kjpindex*nvm*nstm*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (precisol(kjpindex,nvm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in precisol allocation. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (precisol_ns(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in precisol_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (free_drain_coef(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in free_drain_coef allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (frac_bare_ns(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in frac_bare_ns allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (water2infilt(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in water2infilt allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ae_ns(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ae_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (evap_bare_lim_ns(kjpindex,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in evap_bare_lim_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (rootsink(kjpindex,nslm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in rootsink allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (subsnowveg(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in subsnowveg allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (subsnownobio(kjpindex,nnobio),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in subsnownobio allocation. We stop. We need kjpindex words = ',kjpindex*nnobio
       STOP 'hydrol_init'
    END IF

    ALLOCATE (snowmelt(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in snowmelt allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrol_init'
    END IF

    ALLOCATE (icemelt(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in icemelt allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (subsinksoil(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in subsinksoil allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mx_eau_var(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mx_eau_var allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (vegtot(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vegtot allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (resdist(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in resdist allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (humtot(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in humtot allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (resolv(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in resolv allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (k(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in k allocation. We stop. We need kjpindex*nslm words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    IF (ok_freeze_cwrr) THEN
       ALLOCATE (kk_moy(kjpindex,nslm),stat=ier) 
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in kk allocation. We stop. We need kjpindex words = ',kjpindex*nslm
          STOP 'hydrol_init'
       END IF

       ALLOCATE (kk(kjpindex,nslm,nstm),stat=ier) 
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in kk allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
          STOP 'hydrol_init'
       END IF
    ENDIF

    ALLOCATE (a(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in a allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (b(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in b allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (d(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in d allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (e(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in e allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (f(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in f allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (g1(kjpindex,nslm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in g1 allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ep(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ep allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (fp(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in fp allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (gp(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gp allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (rhs(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in rhs allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (srhs(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in srhs allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmcs(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmcs allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmcr(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmcr allocation. We stop. We need kjpindex*nstm words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litt_mea(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litt_mea allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_res(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_res allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_wilt(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_wilt allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_field(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_field allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_sat(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_sat allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_awet(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_awet allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litter_adry(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litter_adry allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litt_wet_mea(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litt_wet_mea allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmc_litt_dry_mea(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmc_litt_dry_mea allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (v1(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in v1 allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (qflux00(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in qflux00 allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (ru_ns(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in ru_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF
    ru_ns(:,:) = zero

    ALLOCATE (dr_ns(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in dr_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF
    dr_ns(:,:) = zero

    ALLOCATE (tr_ns(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tr_ns allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (cvs_over_veg(kjpindex,nvm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in cvs_over_veg allocation. We stop. We need kjpindex*nvm*nstm words = ',kjpindex*nvm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (corr_veg_soil(kjpindex,nvm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in corr_veg_soil allocation. We stop. We need kjpindex*nvm*nstm words = ',kjpindex*nvm*nstm
       STOP 'hydrol_init'
    END IF


    ALLOCATE (mc(kjpindex,nslm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (soilmoist(kjpindex,nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soilmoist allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (soil_wet(kjpindex,nslm,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soil_wet allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (soil_wet_litter(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soil_wet allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (qflux(kjpindex,nslm,nstm),stat=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in qflux allocation. We stop. We need kjpindex words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (tmat(kjpindex,nslm,3),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tmat allocation. We stop. We need kjpindex words = ',kjpindex*nslm*trois
       STOP 'hydrol_init'
    END IF

    ALLOCATE (stmat(kjpindex,nslm,3),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in stmat allocation. We stop. We need kjpindex words = ',kjpindex*nslm*trois
       STOP 'hydrol_init'
    END IF

    ALLOCATE (nroot(nvm, nstm, nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in nroot allocation. We stop. We need nvm*nstm*nslm words = ',nvm * nstm * nslm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (kfact_root(kjpindex, nslm, nstm), stat=ier)
    IF (ier .NE. 0) THEN
       WRITE (numout,*) 'error in kfact_root allocation, We stop. We need kjpindex*nslm*nstm words = ',kjpindex*nslm*nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (kfact(nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in kfact allocation. We stop. We need nslm*nscm words = ',nslm * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (zz(nslm+1, nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in zz allocation. We stop. We need (nslm+1)*nstm words = ',(nslm+1) * nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (dz(nslm+1, nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in dz allocation. We stop. We need (nslm+1)*nstm words = ',(nslm+1) * nstm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mc_lin(imin:imax, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mc_lin allocation. We stop. We need (imax-imin)*nscm words = ',(imax-imin) * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (k_lin(imin:imax, nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in k_lin allocation. We stop. We need (imax-imin)*nslm*nscm words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (d_lin(imin:imax, nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in d_lin allocation. We stop. We need (imax-imin)*nslm*nscm words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (a_lin(imin:imax, nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in a_lin allocation. We stop. We need (imax-imin)*nslm*nscm words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF

    ALLOCATE (b_lin(imin:imax, nslm, nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in b_lin allocation. We stop. We need (imax-imin)*nslm*nscm words = ',(imax-imin) * nslm * nscm
       STOP 'hydrol_init'
    END IF

!pss+ ! WETALND variables allocation
!pss:+
    
    ALLOCATE (fsat(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in fsat allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (fwet(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in fwet allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (fwt1(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in fwt1 allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF
    
    ALLOCATE (fwt2(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in fwt2 allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF
   
    ALLOCATE (fwt3(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in fwt3 allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (fwt4(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in fwt4 allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (drunoff(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in drunoff allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF
       
    ALLOCATE (ZMEAN(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zmean allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF
!    ALLOCATE (NB_PIXE(kjpindex),stat=ier)
!    IF (ier.NE.0) THEN
!        WRITE (numout,*) ' error in mx_eau_var allocation. We stop. We need kjpindex words = ',kjpindex
!        STOP 'hydrolc_init'
!    END IF
    ALLOCATE (ZSTDT(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zstdt allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF
    
    ALLOCATE (ZSKEW(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zskew allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (ZMIN(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zmin allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF
    
    ALLOCATE (ZMAX(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zmax allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF
    
    ALLOCATE (ZM(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ZM allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (ZZPAS(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zzpas allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (ZTAB_FSAT(kjpindex,1000),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error inztab_fsat allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (ZTAB_WTOP(kjpindex,1000),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ztab_wtop allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (ZTAB_FWET(kjpindex,1000),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ztab_fwet allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

    ALLOCATE (ZTAB_WTOP_WET(kjpindex,1000),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ztab_wtop_wet allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'hydrolc_init'
    END IF

!pss+
    ALLOCATE (mcw_grid(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcw_grid allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF

    ALLOCATE (mcs_grid(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcs_grid allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    END IF
!pss-


    IF (ok_freeze_cwrr) THEN
       ALLOCATE (profil_froz_hydro(kjpindex, nslm),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in profil_froz_hydro allocation. We stop. ISA'
          STOP 'hydrol_init'
       END IF

       ALLOCATE (profil_froz_hydro_ns(kjpindex, nslm, nstm),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in profil_froz_hydro_ns allocation. We stop. ISA'
          STOP 'hydrol_init'
       END IF

       ALLOCATE (temp_hydro(kjpindex, nslm),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in temp_hydro allocation. We stop. ISA'
          STOP 'hydrol_init'
       END IF
    ENDIF

    ALLOCATE (mcl(kjpindex, nslm, nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in mcl allocation. We stop. ISA'
       STOP 'hydrol_init'
    END IF

    ALLOCATE (frac_hydro_diag(nslm, nbdl),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in frac_hydro_diag allocation. We stop'
       STOP 'hydrol_init'
    END IF

    !  If we check the water balance we need two more variables
    IF ( check_waterbal ) THEN

       ALLOCATE (tot_water_beg(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_water_beg allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_water_end(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_water_end allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_flux(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_flux allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

    ENDIF

    ! Soil Wetness Index
    ALLOCATE (swi(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in swi. We stop. We need kjpindex words = ',kjpindex
       STOP 'hydrol_init'
    ENDIF  
 
    !
    !  If we use the almaoutputs we need a few more variables
    !  tdo - they could be allocated only if alma_output, but then they should also be computed only if alma_output
    !
    IF ( almaoutput ) THEN

       ALLOCATE (tot_watveg_beg(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_watveg_beg allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_watveg_end(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_watveg_end allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_watsoil_beg(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_watsoil_beg allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (tot_watsoil_end(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in tot_watsoil_end allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (delsoilmoist(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in delsoilmoist allocation. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (delintercept(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in delintercept. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (delswe(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in delswe. We stop. We need kjpindex words = ',kjpindex
          STOP 'hydrol_init'
       ENDIF

       ALLOCATE (snow_beg(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in snow_beg allocation. We stop. We need kjpindex words =',kjpindex
          STOP 'hydrol_init'
       END IF

       ALLOCATE (snow_end(kjpindex),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in snow_end allocation. We stop. We need kjpindex words =',kjpindex
          STOP 'hydrol_init'
       END IF

    ENDIF

    !! 4 Open restart input file and read data for HYDROLOGIC process

    IF (ldrestart_read) THEN

       IF (long_print) WRITE (numout,*) ' we have to read a restart file for HYDROLOGIC variables'

       IF (is_root_prc) CALL ioconf_setatt_p('UNITS', '-')
       !
       DO jst=1,nstm
          ! var_name= "mc_1" ... "mc_3"
           WRITE (var_name,"('moistc_',I1)") jst
           IF (is_root_prc) CALL ioconf_setatt_p('LONG_NAME',var_name)
           CALL restget_p (rest_id, var_name, nbp_glo, nslm , 1, kjit, .TRUE., mc(:,:,jst), "gather", nbp_glo, index_g)
       END DO
       !
       IF (is_root_prc) CALL ioconf_setatt_p('UNITS', '-')
       DO jst=1,nstm
          DO jsl=1,nslm
             ! var_name= "us_1_01" ... "us_3_11"
             WRITE (var_name,"('us_',i1,'_',i2.2)") jst,jsl
             IF (is_root_prc) CALL ioconf_setatt_p('LONG_NAME',var_name)
             CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., us(:,:,jst,jsl), "gather", nbp_glo, index_g)
          END DO
       END DO
       !
       var_name= 'free_drain_coef'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Coefficient for free drainage at bottom of soil')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., free_drain_coef, "gather", nbp_glo, index_g)
       !
       var_name= 'water2infilt'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Remaining water to be infiltrated on top of the soil')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., water2infilt, "gather", nbp_glo, index_g)
       !
       var_name= 'ae_ns'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', 'kg/m^2')
          CALL ioconf_setatt_p('LONG_NAME','Bare soil evap on each soil type')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., ae_ns, "gather", nbp_glo, index_g)
       !
       var_name= 'snow'        
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', 'kg/m^2')
          CALL ioconf_setatt_p('LONG_NAME','Snow mass')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., snow, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_age'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', 'd')
          CALL ioconf_setatt_p('LONG_NAME','Snow age')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., snow_age, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_nobio'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', 'kg/m^2')
          CALL ioconf_setatt_p('LONG_NAME','Snow on other surface types')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nnobio  , 1, kjit, .TRUE., snow_nobio, "gather", nbp_glo, index_g)
       !
       var_name= 'snow_nobio_age'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', 'd')
          CALL ioconf_setatt_p('LONG_NAME','Snow age on other surface types')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nnobio  , 1, kjit, .TRUE., snow_nobio_age, "gather", nbp_glo, index_g)
       !
       var_name= 'vegstress'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Vegetation growth moisture stress')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., vegstress, "gather", nbp_glo, index_g)
       !
       var_name= 'qsintveg'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', 'kg/m^2')
          CALL ioconf_setatt_p('LONG_NAME','Intercepted moisture')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., qsintveg, "gather", nbp_glo, index_g)

!pss:+ !          
       var_name= 'fwet_out'      
       IF (is_root_prc) THEN
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','fwet pr autres routines')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE.,fwet_out , "gather", nbp_glo, index_g)
!pss:-

       IF (ok_explicitsnow) THEN
          var_name= 'grndflux'
          CALL ioconf_setatt('UNITS', 'W/m^2')
          CALL ioconf_setatt('LONG_NAME','ground heat flux')
          CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE.,grndflux,"gather", nbp_glo, index_g)
          CALL setvar_p (grndflux, val_exp, 'ground heat flux', 0.0)
          
          var_name = 'snowrho'
          CALL ioconf_setatt('UNITS', 'kg/m3')
          CALL ioconf_setatt('LONG_NAME','Snow Density profile')
          CALL restget_p (rest_id, var_name, nbp_glo, nsnow, 1, kjit,.TRUE.,snowrho, "gather", nbp_glo, index_g) !need to add veg dim
          CALL setvar_p (snowrho, val_exp, 'Snow Density profile', xrhosmin)
          
          var_name = 'snowtemp'
          CALL ioconf_setatt('UNITS', 'K')
          CALL ioconf_setatt('LONG_NAME','Snow Temperature profile')
          CALL restget_p (rest_id, var_name, nbp_glo,  nsnow, 1, kjit,.TRUE.,snowtemp, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (snowtemp, val_exp, 'Snow Temperature profile', tp_00)
          
          var_name = 'soiltemp'
          CALL ioconf_setatt('UNITS', 'K')
          CALL ioconf_setatt('LONG_NAME','Soil Temperature profile')
          CALL restget_p (rest_id, var_name, nbp_glo,  ngrnd, 1, kjit,.TRUE.,soiltemp, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (soiltemp, val_exp, 'Soil Temperature profile', tp_00)
          
          var_name = 'snowliq'
          CALL ioconf_setatt('UNITS', 'm')
          CALL ioconf_setatt('LONG_NAME','Snow liquid content profile')
          CALL restget_p (rest_id, var_name, nbp_glo,  nsnow, 1, kjit,.TRUE.,snowliq, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (snowliq, val_exp, 'Snow liquid content profile', 0.0)
          
          var_name = 'snowheat'
          CALL ioconf_setatt('UNITS', 'J/m2')
          CALL ioconf_setatt('LONG_NAME','Snow Heat profile')
          CALL restget_p (rest_id, var_name, nbp_glo,  nsnow, 1, kjit,.TRUE.,snowheat, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (snowheat, val_exp, 'Snow Heat profile', 0.0)
          
          var_name = 'snowdz'
          CALL ioconf_setatt('UNITS', 'm')
          CALL ioconf_setatt('LONG_NAME','Snow depth profile')
          CALL restget_p (rest_id, var_name, nbp_glo,  nsnow, 1, kjit,.TRUE.,snowdz, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (snowdz, val_exp, 'Snow depth profile', 0.0)
          
          var_name = 'snowgrain'
          CALL ioconf_setatt('UNITS', 'm')
          CALL ioconf_setatt('LONG_NAME','Snow grain profile')
          CALL restget_p (rest_id, var_name, nbp_glo,  nsnow, 1, kjit,.TRUE.,snowgrain, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (snowgrain, val_exp, 'Snow grain profile', 0.0)


          var_name = 'cgrnd_snow'
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','cgrnd snow coefficient')
          CALL restget_p (rest_id, var_name, nbp_glo,  nsnow, 1, kjit,.TRUE.,cgrnd_snow, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (cgrnd_snow, val_exp, 'cgrnd snow coefficient', 273.15)

          var_name = 'dgrnd_snow'
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','dgrnd snow coefficient')
          CALL restget_p (rest_id, var_name, nbp_glo,  nsnow, 1, kjit,.TRUE.,dgrnd_snow, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (dgrnd_snow, val_exp, 'dgrnd snow coefficient', 0.0)

          var_name = 'zdz1_soil'
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','zdz1 soil coefficient')
          CALL restget_p (rest_id, var_name, nbp_glo,  1,1, kjit,.TRUE.,zdz1_soil, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (zdz1_soil, val_exp, 'zdz1 soil coefficient', 0.0)

          var_name = 'zdz2_soil'
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','zdz2 soil coefficient')
          CALL restget_p (rest_id, var_name, nbp_glo,  1, 1, kjit,.TRUE.,zdz2_soil, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (zdz2_soil, val_exp, 'zdz2 soil coefficient', 0.0)

          var_name = 'cgrnd_soil'
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','cgrnd soil coefficient')
          CALL restget_p (rest_id, var_name, nbp_glo,  1,1, kjit,.TRUE.,cgrnd_soil, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (cgrnd_soil, val_exp, 'cgrnd soil coefficient', 0.0)

          var_name = 'dgrnd_soil'
          CALL ioconf_setatt('UNITS', '-')
          CALL ioconf_setatt('LONG_NAME','dgrnd soil coefficient')
          CALL restget_p (rest_id, var_name, nbp_glo,  1,1, kjit,.TRUE.,dgrnd_soil, "gather", nbp_glo,  index_g) !need to add veg dim
          CALL setvar_p (dgrnd_soil, val_exp, 'dgrnd soil coefficient', 0.0)

!          var_name = 'snowflx'
!          CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit,.TRUE.,snowflx,"gather", nbp_glo,  index_g) !need to add veg dim
!          CALL setvar_p (snowflx, val_exp, 'NO_KEYWORD', zero)
          
!          var_name = 'snowcap'
!          CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit,.TRUE.,snowcap,"gather", nbp_glo,  index_g) !need to add veg dim
!          CALL setvar_p (snowcap, val_exp, 'NO_KEYWORD', zero)
          
       END IF

       var_name= 'resdist'
       IF (is_root_prc) THEN
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Distribution of reservoirs')
       ENDIF
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., resdist, "gather", nbp_glo, index_g)
       
       IF (is_root_prc) CALL ioconf_setatt_p('UNITS', '-')
       DO jst=1,nstm
          ! var_name= "cvs_over_veg_1" ... "cvs_over_veg_3"
          WRITE (var_name,"('cvs_over_veg_',i1)") jst
          IF (is_root_prc) CALL ioconf_setatt_p('LONG_NAME',var_name)
          CALL restget_p (rest_id, var_name, nbp_glo,  nvm, 1, kjit, .TRUE., cvs_over_veg(:,:,jst), "gather",  nbp_glo, index_g)
       END DO
       
       IF ( check_waterbal ) THEN
          var_name= 'tot_water_beg'
          IF (is_root_prc) THEN
             CALL ioconf_setatt_p('UNITS', 'kg/m^2')
             CALL ioconf_setatt_p('LONG_NAME','Previous Total water')
          ENDIF
          CALL restget_p (rest_id, var_name, nbp_glo, 1  , 1, kjit, .TRUE., tot_water_beg, "gather", nbp_glo, index_g)
       ENDIF


    !! 5 get restart values if none were found in the restart file
       !
       !Config Key   = HYDROL_MOISTURE_CONTENT
       !Config Desc  = Soil moisture on each soil tile and levels
       !Config If    = HYDROL_CWRR       
       !Config Def   = 0.3
       !Config Help  = The initial value of mc if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (mc, val_exp, 'HYDROL_MOISTURE_CONTENT', 0.3_r_std)
       !
       !Config Key   = US_INIT
       !Config Desc  = US_NVM_NSTM_NSLM
       !Config If    = HYDROL_CWRR       
       !Config Def   = 0.0
       !Config Help  = The initial value of us (relative moisture) if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       DO jsl=1,nslm
          CALL setvar_p (us(:,:,:,jsl), val_exp, 'US_INIT', zero)
       ENDDO
       !
       !Config Key   = FREE_DRAIN_COEF
       !Config Desc  = Coefficient for free drainage at bottom
       !Config If    = HYDROL_CWRR       
       !Config Def   = 1.0, 1.0, 1.0
       !Config Help  = The initial value of free drainage if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (free_drain_coef, val_exp, 'FREE_DRAIN_COEF', free_drain_max)
       !
       !Config Key   = WATER_TO_INFILT
       !Config Desc  = Water to be infiltrated on top of the soil
       !Config If    = HYDROL_CWRR    
       !Config Def   = 0.0
       !Config Help  = The initial value of free drainage if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (water2infilt, val_exp, 'WATER_TO_INFILT', zero)
       !
       !Config Key   = EVAPNU_SOIL
       !Config Desc  = Bare soil evap on each soil if not found in restart
       !Config If    = HYDROL_CWRR  
       !Config Def   = 0.0
       !Config Help  = The initial value of bare soils evap if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (ae_ns, val_exp, 'EVAPNU_SOIL', zero)
       !
       !Config Key  = HYDROL_SNOW
       !Config Desc  = Initial snow mass if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow mass if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow, val_exp, 'HYDROL_SNOW', zero)
       !
       !Config Key   = HYDROL_SNOWAGE
       !Config Desc  = Initial snow age if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow age if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow_age, val_exp, 'HYDROL_SNOWAGE', zero)
       !
       !Config Key   = HYDROL_SNOW_NOBIO
       !Config Desc  = Initial snow amount on ice, lakes, etc. if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow_nobio, val_exp, 'HYDROL_SNOW_NOBIO', zero)
       !
       !Config Key   = HYDROL_SNOW_NOBIO_AGE
       !Config Desc  = Initial snow age on ice, lakes, etc. if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of snow age if its value is not found
       !Config         in the restart file. This should only be used if the model is 
       !Config         started without a restart file.
       !Config Units =
       !
       CALL setvar_p (snow_nobio_age, val_exp, 'HYDROL_SNOW_NOBIO_AGE', zero)
       !
       !Config Key   = HYDROL_QSV
       !Config Desc  = Initial water on canopy if not found in restart
       !Config If    = OK_SECHIBA
       !Config Def   = 0.0
       !Config Help  = The initial value of moisture on canopy if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units =
       !
       CALL setvar_p (qsintveg, val_exp, 'HYDROL_QSV', zero)

       IF (ok_freeze_cwrr) THEN  
          CALL setvar_p (profil_froz_hydro, val_exp, 'NO_KEYWORD', zero)
          CALL setvar_p (profil_froz_hydro_ns, val_exp, 'NO_KEYWORD', zero)
          CALL setvar_p (kk, val_exp, 'NO_KEYWORD', 276.48)
          CALL setvar_p (kk_moy, val_exp, 'NO_KEYWORD', 276.48)
          CALL setvar_p (temp_hydro, val_exp, 'NO_KEYWORD', 280.)
       ENDIF
       
!pss:+
       !
       !Config Key   = HYDROL_FWET
       !Config Desc  = Initial fwet_out if not found in restart
       !Config If    = TOPM_calcul
       !Config Def   = 0.0
       !Config Help  = The initial value of fwet_out if its value 
       !Config         is not found in the restart file. This should only be used if
       !Config         the model is started without a restart file. 
       !Config Units =
       CALL setvar_p (fwet_out, val_exp,'HYDROL_FWET', zero)

!pss:-

    !! 6 Vegetation array      
       !
       ! There is no need to configure the initialisation of resdist. If not available it is the vegetation map
       !
       IF ( MINVAL(resdist) .EQ.  MAXVAL(resdist) .AND. MINVAL(resdist) .EQ. val_exp) THEN
          resdist(:,:) = soiltile(:,:)
       ENDIF
       !
       !  Remember that it is only frac_nobio + SUM(veget_max(,:)) that is equal to 1. Thus we need vegtot
       !
       DO ji = 1, kjpindex
          vegtot(ji) = SUM(veget_max(ji,:))
       ENDDO
       !
       !
       ! compute the masks for veget


       mask_veget(:,:) = 0
       mask_soiltile(:,:) = 0

       DO jst=1,nstm
          DO ji = 1, kjpindex
             IF(soiltile(ji,jst) .GT. min_sechiba) THEN
                mask_soiltile(ji,jst) = 1
             ENDIF
          END DO
       ENDDO
          
       DO jv = 1, nvm
          DO ji = 1, kjpindex
             IF(veget_max(ji,jv) .GT. min_sechiba) THEN
                mask_veget(ji,jv) = 1
             ENDIF
          END DO
       END DO
          
    !! 7 set humrelv from us

       humrelv(:,:,:) = SUM(us,dim=4)
       vegstressv(:,:,:) = humrelv(:,:,:)
       ! set humrel from humrelv, assuming equi-repartition for the first time step

       humrel(:,:) = zero
       CALL setvar_p (cvs_over_veg, val_exp, 'NO_KEYWORD', un)

       DO jst=1,nstm
          DO jv=1,nvm
             DO ji=1,kjpindex

                vegstress(ji,jv)=vegstress(ji,jv) + vegstressv(ji,jv,jst) * &
                     & soiltile(ji,jst) * cvs_over_veg(ji,jv,jst) * vegtot(ji)

                humrel(ji,jv)=humrel(ji,jv) + humrelv(ji,jv,jst) * & 
                     & soiltile(ji,jst) * cvs_over_veg(ji,jv,jst) * vegtot(ji)
                humrel(ji,jv)=MAX(humrel(ji,jv), zero)* mask_veget(ji,jv)           
             END DO
          END DO
       END DO
    ENDIF
    !
    !
    IF (long_print) WRITE (numout,*) ' hydrol_init done '
    !
  END SUBROUTINE hydrol_init


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_clear
!!
!>\BRIEF        Deallocate arrays 
!!
!_ ================================================================================================================================
!_ hydrol_clear

  SUBROUTINE hydrol_clear()

    l_first_hydrol=.TRUE.

    ! Allocation for soiltile related parameters
    IF ( ALLOCATED (nvan)) DEALLOCATE (nvan)
    IF ( ALLOCATED (avan)) DEALLOCATE (avan)
    IF ( ALLOCATED (mcr)) DEALLOCATE (mcr)
    IF ( ALLOCATED (mcs)) DEALLOCATE (mcs)
    IF ( ALLOCATED (ks)) DEALLOCATE (ks)
    IF ( ALLOCATED (pcent)) DEALLOCATE (pcent)
    IF ( ALLOCATED (free_drain_max)) DEALLOCATE (free_drain_max)
    IF ( ALLOCATED (mcf)) DEALLOCATE (mcf)
    IF ( ALLOCATED (mcw)) DEALLOCATE (mcw)
    IF ( ALLOCATED (mc_awet)) DEALLOCATE (mc_awet)
    IF ( ALLOCATED (mc_adry)) DEALLOCATE (mc_adry)
!pss+
    IF ( ALLOCATED (mcs_grid)) DEALLOCATE (mcs_grid)
    IF ( ALLOCATED (mcw_grid)) DEALLOCATE (mcw_grid)
!pss-
    ! Other arrays
    IF (ALLOCATED (mask_veget)) DEALLOCATE (mask_veget)
    IF (ALLOCATED (mask_soiltile)) DEALLOCATE (mask_soiltile)
    IF (ALLOCATED (humrelv)) DEALLOCATE (humrelv)
    IF (ALLOCATED (vegstressv)) DEALLOCATE (vegstressv)
    IF (ALLOCATED (us)) DEALLOCATE (us)
    IF (ALLOCATED  (precisol)) DEALLOCATE (precisol)
    IF (ALLOCATED  (precisol_ns)) DEALLOCATE (precisol_ns)
    IF (ALLOCATED  (free_drain_coef)) DEALLOCATE (free_drain_coef)
    IF (ALLOCATED  (frac_bare_ns)) DEALLOCATE (frac_bare_ns)
    IF (ALLOCATED  (water2infilt)) DEALLOCATE (water2infilt)
    IF (ALLOCATED  (ae_ns)) DEALLOCATE (ae_ns)
    IF (ALLOCATED  (evap_bare_lim_ns)) DEALLOCATE (evap_bare_lim_ns)
    IF (ALLOCATED  (rootsink)) DEALLOCATE (rootsink)
    IF (ALLOCATED  (subsnowveg)) DEALLOCATE (subsnowveg)
    IF (ALLOCATED  (subsnownobio)) DEALLOCATE (subsnownobio)
    IF (ALLOCATED  (snowmelt)) DEALLOCATE (snowmelt)
    IF (ALLOCATED  (icemelt)) DEALLOCATE (icemelt)
    IF (ALLOCATED  (subsinksoil)) DEALLOCATE (subsinksoil)
    IF (ALLOCATED  (mx_eau_var)) DEALLOCATE (mx_eau_var)
    IF (ALLOCATED  (vegtot)) DEALLOCATE (vegtot)
    IF (ALLOCATED  (resdist)) DEALLOCATE (resdist)
    IF (ALLOCATED  (tot_water_beg)) DEALLOCATE (tot_water_beg)
    IF (ALLOCATED  (tot_water_end)) DEALLOCATE (tot_water_end)
    IF (ALLOCATED  (tot_flux)) DEALLOCATE (tot_flux)
    IF (ALLOCATED  (tot_watveg_beg)) DEALLOCATE (tot_watveg_beg)
    IF (ALLOCATED  (tot_watveg_end)) DEALLOCATE (tot_watveg_end)
    IF (ALLOCATED  (tot_watsoil_beg)) DEALLOCATE (tot_watsoil_beg)
    IF (ALLOCATED  (tot_watsoil_end)) DEALLOCATE (tot_watsoil_end)
    IF (ALLOCATED  (delsoilmoist)) DEALLOCATE (delsoilmoist)
    IF (ALLOCATED  (delintercept)) DEALLOCATE (delintercept)
    IF (ALLOCATED  (snow_beg)) DEALLOCATE (snow_beg)
    IF (ALLOCATED  (snow_end)) DEALLOCATE (snow_end)
    IF (ALLOCATED  (delswe)) DEALLOCATE (delswe)
    IF (ALLOCATED  (swi)) DEALLOCATE (swi)
    IF (ALLOCATED  (v1)) DEALLOCATE (v1)
    IF (ALLOCATED  (humtot)) DEALLOCATE (humtot)
    IF (ALLOCATED  (resolv)) DEALLOCATE (resolv)
    IF (ALLOCATED  (k)) DEALLOCATE (k)
    IF (ALLOCATED  (kk)) DEALLOCATE (kk)
    IF (ALLOCATED  (kk_moy)) DEALLOCATE (kk_moy)
    IF (ALLOCATED  (a)) DEALLOCATE (a)
    IF (ALLOCATED  (b)) DEALLOCATE (b)
    IF (ALLOCATED  (d)) DEALLOCATE (d)
    IF (ALLOCATED  (e)) DEALLOCATE (e)
    IF (ALLOCATED  (f)) DEALLOCATE (f)
    IF (ALLOCATED  (g1)) DEALLOCATE (g1)
    IF (ALLOCATED  (ep)) DEALLOCATE (ep)
    IF (ALLOCATED  (fp)) DEALLOCATE (fp)
    IF (ALLOCATED  (gp)) DEALLOCATE (gp)
    IF (ALLOCATED  (rhs)) DEALLOCATE (rhs)
    IF (ALLOCATED  (srhs)) DEALLOCATE (srhs)
    IF (ALLOCATED  (tmc)) DEALLOCATE (tmc)
    IF (ALLOCATED  (tmcs)) DEALLOCATE (tmcs)
    IF (ALLOCATED  (tmcr)) DEALLOCATE (tmcr)
    IF (ALLOCATED  (tmc_litter)) DEALLOCATE (tmc_litter)
    IF (ALLOCATED  (tmc_litt_mea)) DEALLOCATE (tmc_litt_mea)
    IF (ALLOCATED  (tmc_litter_res)) DEALLOCATE (tmc_litter_res)
    IF (ALLOCATED  (tmc_litter_wilt)) DEALLOCATE (tmc_litter_wilt)
    IF (ALLOCATED  (tmc_litter_field)) DEALLOCATE (tmc_litter_field)
    IF (ALLOCATED  (tmc_litter_sat)) DEALLOCATE (tmc_litter_sat)
    IF (ALLOCATED  (tmc_litter_awet)) DEALLOCATE (tmc_litter_awet)
    IF (ALLOCATED  (tmc_litter_adry)) DEALLOCATE (tmc_litter_adry)
    IF (ALLOCATED  (tmc_litt_wet_mea)) DEALLOCATE (tmc_litt_wet_mea)
    IF (ALLOCATED  (tmc_litt_dry_mea)) DEALLOCATE (tmc_litt_dry_mea)
    IF (ALLOCATED  (qflux00)) DEALLOCATE (qflux00)
    IF (ALLOCATED  (ru_ns)) DEALLOCATE (ru_ns)
    IF (ALLOCATED  (dr_ns)) DEALLOCATE (dr_ns)
    IF (ALLOCATED  (tr_ns)) DEALLOCATE (tr_ns)
    IF (ALLOCATED  (cvs_over_veg)) DEALLOCATE (cvs_over_veg)
    IF (ALLOCATED  (corr_veg_soil)) DEALLOCATE (corr_veg_soil)
    IF (ALLOCATED  (mc)) DEALLOCATE (mc)
    IF (ALLOCATED  (soilmoist)) DEALLOCATE (soilmoist)
    IF (ALLOCATED  (soil_wet)) DEALLOCATE (soil_wet)
    IF (ALLOCATED  (soil_wet_litter)) DEALLOCATE (soil_wet_litter)
    IF (ALLOCATED  (qflux)) DEALLOCATE (qflux)
    IF (ALLOCATED  (tmat)) DEALLOCATE (tmat)
    IF (ALLOCATED  (stmat)) DEALLOCATE (stmat)
    IF (ALLOCATED  (nroot)) DEALLOCATE (nroot)
    IF (ALLOCATED  (kfact_root)) DEALLOCATE (kfact_root)
    IF (ALLOCATED  (kfact)) DEALLOCATE (kfact)
    IF (ALLOCATED  (zz)) DEALLOCATE (zz)
    IF (ALLOCATED  (dz)) DEALLOCATE (dz)
    IF (ALLOCATED  (mc_lin)) DEALLOCATE (mc_lin)
    IF (ALLOCATED  (k_lin)) DEALLOCATE (k_lin)
    IF (ALLOCATED  (d_lin)) DEALLOCATE (d_lin)
    IF (ALLOCATED  (a_lin)) DEALLOCATE (a_lin)
    IF (ALLOCATED  (b_lin)) DEALLOCATE (b_lin)
    IF (ALLOCATED  (frac_hydro_diag)) DEALLOCATE (frac_hydro_diag)
    
!pss:+ !WETLAND variables
    IF (ALLOCATED  (fsat))  DEALLOCATE (fsat)
    IF (ALLOCATED  (fwet))  DEALLOCATE (fwet)
    IF (ALLOCATED  (fwt1))  DEALLOCATE (fwt1)
    IF (ALLOCATED  (fwt2))  DEALLOCATE (fwt2)
    IF (ALLOCATED  (fwt3))  DEALLOCATE (fwt3)
    IF (ALLOCATED  (fwt4))  DEALLOCATE (fwt4)
    IF (ALLOCATED  (drunoff))  DEALLOCATE (drunoff)
    IF (ALLOCATED  (ZMEAN)) DEALLOCATE (ZMEAN)
!    IF (ALLOCATED  (NB_PIXE)) DEALLOCATE (NB_PIXE)
    IF (ALLOCATED  (ZSTDT)) DEALLOCATE (ZSTDT)
    IF (ALLOCATED  (ZSKEW)) DEALLOCATE (ZSKEW)
    IF (ALLOCATED  (ZMIN)) DEALLOCATE (ZMIN)
    IF (ALLOCATED  (ZMAX)) DEALLOCATE (ZMAX)
    IF (ALLOCATED  (ZM)) DEALLOCATE (ZM)
    IF (ALLOCATED  (ZZPAS)) DEALLOCATE (ZZPAS)
    IF (ALLOCATED  (ZTAB_FSAT)) DEALLOCATE (ZTAB_FSAT)
    IF (ALLOCATED  (ZTAB_WTOP)) DEALLOCATE (ZTAB_WTOP)
    IF (ALLOCATED  (ZTAB_FWET)) DEALLOCATE (ZTAB_FWET)
    IF (ALLOCATED  (ZTAB_WTOP_WET)) DEALLOCATE (ZTAB_WTOP_WET)
!pss:-

  END SUBROUTINE hydrol_clear

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_tmc_update
!!
!>\BRIEF        This routine updates the soil moisture profiles when the vegetation fraction have changed. 
!!
!! DESCRIPTION  :
!! 
!!    This routine update tmc and mc with variation of veget_max (LAND_USE or DGVM activated)
!! 
!!
!!
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_tmc_update
  SUBROUTINE hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg, resdist )

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! domain size
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max     !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile      !! Fraction of each soil tile (0-1, unitless)

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)   :: qsintveg   !! Amount of water in the canopy interception 
    REAL(r_std), DIMENSION (kjpindex, nstm), INTENT(inout) :: resdist    !! Old vegetation map

    !! 0.4 Local variables
    INTEGER(i_std)                           :: ji, jv, jst,jsl
    LOGICAL                                  :: soil_upd        !! True if soiltile changed since last time step
    LOGICAL                                  :: error=.FALSE.   !! If true, exit in the end of subroutine
    REAL(r_std), DIMENSION(kjpindex,nstm)    :: vmr             !! Change in soiltile
    REAL(r_std), DIMENSION(kjpindex)         :: vmr_sum
    REAL(r_std), DIMENSION(kjpindex,nslm)    :: mc_dilu         !! Total loss of moisture content
    REAL(r_std), DIMENSION(kjpindex)         :: infil_dilu      !! Total loss for water2infilt
    REAL(r_std), DIMENSION(kjpindex,nstm)    :: tmc_old         !! tmc before calculations
    REAL(r_std), DIMENSION(kjpindex,nstm)    :: water2infilt_old!! water2infilt before calculations
    REAL(r_std), DIMENSION (kjpindex,nvm)    :: qsintveg_old    !! qsintveg before calculations
    REAL(r_std), DIMENSION(kjpindex)         :: test

    !! 0. Check if soiltiles changed since last time step
    soil_upd=SUM(ABS(soiltile(:,:)-resdist(:,:))) .GT. zero
!    WRITE(numout,*) 'zdcheck1 soil_upd',soil_upd,'soiltile',soiltile(1,:),'resdist',resdist(1,:),'soiltile-resdist',soiltile(1,:)-resdist(1,:)
    IF (long_print) WRITE (numout,*) 'soil_upd ', soil_upd

    IF (check_cwrr) THEN
       ! Save soil moisture for later use
       tmc_old(:,:) = tmc(:,:) 
       water2infilt_old(:,:) = water2infilt(:,:)
       qsintveg_old(:,:) = qsintveg(:,:)
    ENDIF

    !! 1. If a PFT has disapperead as result from a veget_max change, 
    !!    then add canopy water to surface water.

    DO ji=1,kjpindex
       DO jv=1,nvm
          IF ((veget_max(ji,jv).LT.min_sechiba).AND.(qsintveg(ji,jv).GT.0.)) THEN
             jst=pref_soil_veg(jv) ! soil tile index
             water2infilt(ji,jst) = water2infilt(ji,jst) + qsintveg(ji,jv)/resdist(ji,jst)
             qsintveg(ji,jv) = zero
          ENDIF
       ENDDO
    ENDDO
    
    !! 2. Compute new soil moisture if soiltile changed due to veget_max's change
    IF (soil_upd) THEN
       !! 2.1 Define the change in soiltile
       vmr(:,:) = soiltile(:,:) - resdist(:,:)  ! resdist is the previous values of soiltiles 

       ! Total area loss by the three soil tiles
       DO ji=1,kjpindex
          vmr_sum(ji)=SUM(vmr(ji,:),MASK=vmr(ji,:).LT.zero)
       ENDDO

       !! 2.2 Shrinking soil tiles
       !! 2.2.1 Total loss of moisture content from the shrinking soil tiles, expressed by soil layer
       mc_dilu(:,:)=zero
       DO jst=1,nstm
          DO jsl = 1, nslm
             DO ji=1,kjpindex
                IF ( vmr(ji,jst) < zero ) THEN
                   mc_dilu(ji,jsl) = mc_dilu(ji,jsl) + mc(ji,jsl,jst) * vmr(ji,jst) / vmr_sum(ji)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       !! 2.2.2 Total loss of water2inft from the shrinking soil tiles
       infil_dilu(:)=zero
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( vmr(ji,jst) < zero ) THEN
                infil_dilu(ji) = infil_dilu(ji) + water2infilt(ji,jst) * vmr(ji,jst) / vmr_sum(ji)
             ENDIF
          ENDDO
       ENDDO

       !! 2.3 Each gaining soil tile gets moisture proportionally to both the total loss and its areal increase 

       ! As the original mc from each soil tile are in [mcr,mcs] and we do weighted avrage, the new mc are in [mcr,mcs]
       ! The case where the soiltile is created (soiltile_old=0) works as the other cases

       ! 2.3.1 Update mc(kjpindex,nslm,nstm) !m3/m3
       DO jst=1,nstm
          DO jsl = 1, nslm
             DO ji=1,kjpindex
                IF ( vmr(ji,jst) > zero ) THEN
                   mc(ji,jsl,jst) = ( mc(ji,jsl,jst) * resdist(ji,jst) + mc_dilu(ji,jsl) * vmr(ji,jst) ) / soiltile(ji,jst)
                   ! NB : soiltile can not be zero for case vmr > zero, see slowproc_veget
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
!       DO jst=1,nstm
!          IF ( vmr(1,jst) > zero ) THEN
!             WRITE(numout,*) 'zdcheck2 jst=',jst,'soiltile,resdist,vmr',soiltile(1,jst),resdist(1,jst),vmr(1,jst)
!          ENDIF
!       ENDDO
       
       ! 2.3.2 Update water2inft
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( vmr(ji,jst) > zero ) THEN !donc soiltile>0     
                water2infilt(ji,jst) = ( water2infilt(ji,jst) * resdist(ji,jst) + infil_dilu(ji) * vmr(ji,jst) ) / soiltile(ji,jst)
             ENDIF !donc resdist>0
          ENDDO
       ENDDO

       ! 2.3.3 Case where soiltile < min_sechiba 
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF ( soiltile(ji,jst) .LT. min_sechiba ) THEN
                water2infilt(ji,jst) = zero
                mc(ji,:,jst) = zero
             ENDIF
          ENDDO
       ENDDO


    ENDIF ! soil_upd 


    !2.3.3 we compute tmc(kjpindex,nstm) and humtot!
    DO jst=1,nstm
       DO ji=1,kjpindex
             tmc(ji,jst) = dz(2,jst) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                     + dz(jsl+1,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
             ENDDO
             tmc(ji,jst) = tmc(ji,jst) + dz(nslm,jst) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
             tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
       ENDDO
    ENDDO

    humtot(:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          humtot(ji) = humtot(ji) + soiltile(ji,jst) * tmc(ji,jst)
       ENDDO
    ENDDO

    !! 4 check
    IF (check_cwrr) THEN
       DO ji=1,kjpindex
          test(ji) = ABS(SUM(tmc(ji,:)*soiltile(ji,:)) - SUM(tmc_old(ji,:)*resdist(ji,:)) + &
               SUM(qsintveg(ji,:)) - SUM(qsintveg_old(ji,:))) ! sum(soiltile)=1
          IF ( test(ji) .GT.  10.*allowed_err ) THEN
             WRITE(numout,*) 'tmc update WRONG: ji',ji
             WRITE(numout,*) 'tot water avant:',SUM(tmc_old(ji,:)*resdist(ji,:)) + SUM(qsintveg_old(ji,:))
             WRITE(numout,*) 'tot water apres:',SUM(tmc(ji,:)*soiltile(ji,:)) + SUM(qsintveg(ji,:))
             WRITE(numout,*) 'err:',test(ji)
             WRITE(numout,*) 'allowed_err:',allowed_err
             WRITE(numout,*) 'tmc:',tmc(ji,:)
             WRITE(numout,*) 'tmc_old:',tmc_old(ji,:)
             WRITE(numout,*) 'qsintveg:',qsintveg(ji,:)
             WRITE(numout,*) 'qsintveg_old:',qsintveg_old(ji,:)
             WRITE(numout,*) 'SUMqsintveg:',SUM(qsintveg(ji,:))
             WRITE(numout,*) 'SUMqsintveg_old:',SUM(qsintveg_old(ji,:))
             WRITE(numout,*) 'veget_max:',veget_max(ji,:)
             WRITE(numout,*) 'soiltile:',soiltile(ji,:)
             WRITE(numout,*) 'resdist:',resdist(ji,:)
             WRITE(numout,*) 'vmr:',vmr(ji,:)
             WRITE(numout,*) 'vmr_sum:',vmr_sum(ji)
             DO jst=1,nstm
                WRITE(numout,*) 'mc(',jst,'):',mc(ji,:,jst)
             ENDDO
             WRITE(numout,*) 'water2infilt:',water2infilt(ji,:)
             WRITE(numout,*) 'water2infilt_old:',water2infilt_old(ji,:)
             WRITE(numout,*) 'infil_dilu:',infil_dilu(ji)
             WRITE(numout,*) 'mc_dilu:',mc_dilu(ji,:)

             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_tmc_update', 'Error in water balance', 'We STOP in the end of this subroutine','')
          ENDIF
       ENDDO
    ENDIF

    !! Now that the work is done, update resdist
    resdist(:,:) = soiltile(:,:)

    !
    !!  Exit if error was found previously in this subroutine
    !
    IF ( error ) THEN
       WRITE(numout,*) 'One or more errors have been detected in hydrol_tmc_update. Model stops.'
       CALL ipslerr_p(3, 'hydrol_tmc_update', 'We will STOP now.',&
                  & 'One or several fatal errors were found previously.','')
    END IF

    IF (long_print) WRITE (numout,*) ' hydrol_tmc_update done '

  END SUBROUTINE hydrol_tmc_update

!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_var_init
!!
!>\BRIEF        This routine initializes HYDROLOGIC variables.  
!!
!! DESCRIPTION  :
!! - 1 compute the depths
!! - 2 compute the profile for roots
!! - 3 compute the profile for ksat, a and n Van Genuchten parameter
!! - 4 compute the linearized values of k, a, b and d for the resolution of Fokker Planck equation
!! - 5 water reservoirs initialisation
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_var_init

  SUBROUTINE hydrol_var_init (kjpindex, veget, veget_max, soiltile, njsc, &
       mx_eau_var, shumdiag, shumdiag_perma, k_litt, &
       litterhumdiag, drysoil_frac, evap_bare_lim,qsintveg) 

    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                          :: kjpindex      !! domain size
    ! input fields
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget_max     !! max fraction of vegetation type
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: veget         !! fraction of vegetation type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc          !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless) 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in) :: soiltile      !! Fraction of each soil tile (0-1, unitless)

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: mx_eau_var    !!
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out) :: shumdiag      !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out) :: shumdiag_perma!! Percent of porosity filled with water (mc/mcs) used for the thermal computations
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: k_litt        !! litter cond.
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: litterhumdiag !! litter humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: drysoil_frac  !! function of litter humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)       :: evap_bare_lim !! 

    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: qsintveg           !! Water on vegetation due to interception

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jv, jd, jst, jsc, jsl, i
    REAL(r_std)                                         :: m, frac
    REAL(r_std)                                         :: avan_mod, nvan_mod !! modified VG parameters with exponantial profile
    REAL(r_std), DIMENSION(nslm,nscm)                   :: afact, nfact  !! multiplicative factor for decay of a and n 
    ! parameters for "soil densification" with depth
    REAL(r_std)                                         :: dp_comp       !! depth at which the 'compacted value' (Van Genuchten) of 
                                                                         !! ksat is reached 
    REAL(r_std)                                         :: f_ks          !! exponential factor for decay of ksat with depth
    ! Fixed parameters from fitted relationships
    REAL(r_std)                                         :: n0            !! fitted value for relation log((n-n0)/(n_ref-n0)) = 
                                                                         !! nk_rel * log(k/k_ref) 
    REAL(r_std)                                         :: nk_rel        !! fitted value for relation log((n-n0)/(n_ref-n0)) = 
                                                                         !! nk_rel * log(k/k_ref) 
    REAL(r_std)                                         :: a0            !! fitted value for relation log((a-a0)/(a_ref-a0)) = 
                                                                         !! ak_rel * log(k/k_ref) 
    REAL(r_std)                                         :: ak_rel        !! fitted value for relation log((a-a0)/(a_ref-a0)) = 
                                                                         !! ak_rel * log(k/k_ref) 
    REAL(r_std)                                         :: kfact_max     !! Maximum factor for K due to depth
    REAL(r_std)                                         :: k_tmp, tmc_litter_ratio
    INTEGER(i_std), PARAMETER                           :: error_level = 3 !! Error level for consistency check
                                                                           !! Switch to 2 tu turn fatal errors into warnings

!_ ================================================================================================================================

!!??Aurelien: Les 3 parametres qui suient pourait peut-Ítre mis dans hydrol_init?
    !
    !
    !Config Key   = CWRR_NKS_N0 
    !Config Desc  = fitted value for relation log((n-n0)/(n_ref-n0)) = nk_rel * log(k/k_ref)
    !Config Def   = 0.95
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    n0 = 0.95
    CALL getin_p("CWRR_NKS_N0",n0)

    !! Check parameter value (correct range)
    IF ( n0 < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_NKS_N0.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_NKS_POWER
    !Config Desc  = fitted value for relation log((n-n0)/(n_ref-n0)) = nk_rel * log(k/k_ref)
    !Config Def   = 0.34
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    nk_rel = 0.34
    CALL getin_p("CWRR_NKS_POWER",nk_rel)

    !! Check parameter value (correct range)
    IF ( nk_rel < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_NKS_POWER.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_AKS_A0 
    !Config Desc  = fitted value for relation log((a-a0)/(a_ref-a0)) = ak_rel * log(k/k_ref)
    !Config Def   = 0.00012
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    a0 = 0.00012
    CALL getin_p("CWRR_AKS_A0",a0)

    !! Check parameter value (correct range)
    IF ( a0 < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_AKS_A0.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = CWRR_AKS_POWER
    !Config Desc  = fitted value for relation log((a-a0)/(a_ref-a0)) = ak_rel * log(k/k_ref)
    !Config Def   = 0.53
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    ak_rel = 0.53
    CALL getin_p("CWRR_AKS_POWER",ak_rel)

    !! Check parameter value (correct range)
    IF ( nk_rel < zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for CWRR_AKS_POWER.", &
            &     "This parameter should be non-negative. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_DECAY_RATE
    !Config Desc  = Factor for Ks decay with depth
    !Config Def   = 2.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =  
    !Config Units = [-]
    f_ks = 2.0
    CALL getin_p ("KFACT_DECAY_RATE", f_ks)

    !! Check parameter value (correct range)
    IF ( f_ks <= zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_DECAY_RATE.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_STARTING_DEPTH
    !Config Desc  = Depth for compacted value of Ks 
    !Config Def   = 0.3
    !Config If    = HYDROL_CWRR 
    !Config Help  =  
    !Config Units = [-]
    dp_comp = 0.3
    CALL getin_p ("KFACT_STARTING_DEPTH", dp_comp)

    !! Check parameter value (correct range)
    IF ( dp_comp <= zero ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_STARTING_DEPTH.", &
            &     "This parameter should be positive. ", &
            &     "Please, check parameter value in run.def. ")
    END IF


    !Config Key   = KFACT_MAX
    !Config Desc  = Maximum Factor for Ks increase due to vegetation
    !Config Def   = 10.0
    !Config If    = HYDROL_CWRR 
    !Config Help  =
    !Config Units = [-]
    kfact_max = 10.0
    CALL getin_p ("KFACT_MAX", kfact_max)

    !! Check parameter value (correct range)
    IF ( kfact_max < 10. ) THEN
       CALL ipslerr_p(error_level, "hydrol_var_init.", &
            &     "Wrong parameter value for KFACT_MAX.", &
            &     "This parameter should be greater than 10. ", &
            &     "Please, check parameter value in run.def. ")
    END IF

    !
    DO jst=1,nstm
       !-
    !! 1 compute the depths
       !-
       zz(1,jst) = zero
       dz(1,jst) = zero
       DO jsl=2,nslm
          zz(jsl,jst) = dpu_max* mille*((2**(jsl-1))-1)/ ((2**(nslm-1))-1)
          dz(jsl,jst) = zz(jsl,jst)-zz(jsl-1,jst)
       ENDDO
       zz(nslm+1,jst) = zz(nslm,jst)
       dz(nslm+1,jst) = zero

       !-
    !! 2 compute the profile for roots
       !-
       !! The three following equations concerning nroot computation are derived from the integrals 
       !! of equations C9 to C11 of De Rosnay's (1999) PhD thesis (page 158).
       !! The occasional absenceof minus sign before humcste parameter is correct.
       !! Arsene : sum(nroot)~=1  // For moss humscste = 50 >> 10 x grass. nroot max for nslm=4 !! Arsene 19-03-2014
       !! Arsene : Grass: 0-1-2-4-8-14-21-25-17-5-2                     !! Arsene 19-03-2014
       !! Arsene : Moss:  0-13-22-29-25-10-1... (humcste=50)            !! Arsene 19-03-2014
       !! Arsene : humcste=25: 0-7-13-20-26-23-9.5-1...                 !! Arsene 19-03-2014
       !! Arsene : zz = profondeur / dz = zz+1 - zz
       DO jv = 1,nvm
          DO jsl = 2, nslm-1
             nroot(jv,jst,jsl) = (EXP(-humcste(jv)*zz(jsl,jst)/mille)) * &
                     & (EXP(humcste(jv)*dz(jsl,jst)/mille/deux) - &
                     & EXP(-humcste(jv)*dz(jsl+1,jst)/mille/deux))/ &
                     & (EXP(-humcste(jv)*dz(2,jst)/mille/deux) &
                     & -EXP(-humcste(jv)*zz(nslm,jst)/mille))
          ENDDO
       ENDDO
       DO jv=1,nvm
          nroot(jv,jst,1) = zero
          nroot(jv,jst,nslm) = (EXP(humcste(jv)*dz(nslm,jst)/mille/deux) -un) * &
                  & EXP(-humcste(jv)*zz(nslm,jst)/mille) / &
                  & (EXP(-humcste(jv)*dz(2,jst)/mille/deux) &
                  & -EXP(-humcste(jv)*zz(nslm,jst)/mille))
       ENDDO
    ENDDO

       ! An additional exponential factor for ks depending on the amount of roots in the soil 
       ! through a geometric average over the vegets
       !!??Aurelien: Pkoi utiliser ks_usda?
!! Pas modif pour le moment... Arsene 18-02-2014
    kfact_root(:,:,:) = un
       DO jsl = 1, nslm
          DO jv = 2, nvm
             jst = pref_soil_veg(jv)
             DO ji = 1, kjpindex
                IF(soiltile(ji,jst) .GT. min_sechiba) THEN
                   kfact_root(ji,jsl,jst) = kfact_root(ji,jsl,jst) * &
                        & MAX((MAXVAL(ks_usda)/ks(njsc(ji)))**(- veget(ji,jv) * (humcste(jv)*zz(jsl,jst)/mille - un)/deux), &
                        un) 
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    !-
    !! 3 Compute the profile for ksat, a and n
    !-

       ! For every soil texture
       DO jsc = 1, nscm
          DO jsl=1,nslm
             kfact(jsl,jsc) = MIN(MAX(EXP(- f_ks * (zz(jsl,jsc)/mille - dp_comp)), un/kfact_max),un)
             !!??Aurelien: comment nk_rel et ak_rel ont-ils ete choisi? Et ce sont les memes pour chaque type de sol !?!
             nfact(jsl,jsc) = ( kfact(jsl,jsc) )**nk_rel
             afact(jsl,jsc) = ( kfact(jsl,jsc) )**ak_rel
          ENDDO
       ENDDO

    
    ! For every soil texture
    DO jsc = 1, nscm
       !-
    !! 4 compute the linearized values of k, a, b and d
       !-
    ! Calcul the matrix coef for dublin model:
    ! pice-wise linearised hydraulic conductivity k_lin=alin * mc_lin + b_lin
    ! and diffusivity d_lin in each interval of mc, called mc_lin,
    ! between imin, for residual mcr, 
    ! and imax for saturation mcs.

       mc_lin(imin,jsc)=mcr(jsc)
       mc_lin(imax,jsc)=mcs(jsc)
       DO ji= imin+1, imax-1 
          mc_lin(ji,jsc) = mcr(jsc) + (ji-imin)*(mcs(jsc)-mcr(jsc))/(imax-imin)
       ENDDO

       DO jsl = 1, nslm
          nvan_mod = n0 + (nvan(jsc)-n0) * nfact(jsl,jsc)
          avan_mod = a0 + (avan(jsc)-a0) * afact(jsl,jsc)
          !!??Aurelien: comment n0 et a0 ont-ils ete choisi? Et ce sont les memes pour chaque type de sol !?!
          m = un - un / nvan_mod
          ! Perhaps we may have problems with precision here for some machines (k being very small for ji=imin+1
          ! How can we handle it?
          DO ji = imax,imin+1,-1 
             frac=MIN(un,(mc_lin(ji,jsc)-mcr(jsc))/(mcs(jsc)-mcr(jsc)))
             k_lin(ji,jsl,jsc) = ks(jsc) * kfact(jsl,jsc) * (frac**0.5) * ( un - ( un - frac ** (un/m)) ** m )**2
          ENDDO
          ! We have to avoid k=0
          k_lin(imin,jsl,jsc) = k_lin(imin+1,jsl,jsc)/mille
          
          DO ji = imin,imax-1 
             a_lin(ji,jsl,jsc) = (k_lin(ji+1,jsl,jsc)-k_lin(ji,jsl,jsc)) / (mc_lin(ji+1,jsc)-mc_lin(ji,jsc))
             b_lin(ji,jsl,jsc)  = k_lin(ji,jsl,jsc) - a_lin(ji,jsl,jsc)*mc_lin(ji,jsc)

             IF(ji.NE.imin.AND.ji.NE.imax-1) THEN
                frac=MIN(un,(mc_lin(ji,jsc)-mcr(jsc))/(mcs(jsc)-mcr(jsc)))
                d_lin(ji,jsl,jsc) =(k_lin(ji,jsl,jsc) / (avan_mod*m*nvan_mod)) *  &
                     &  ( (frac**(-un/m))/(mc_lin(ji,jsc)-mcr(jsc)) ) * &
                  &  (  frac**(-un/m) -un ) ** (-m)
                frac=MIN(un,(mc_lin(ji+1,jsc)-mcr(jsc))/(mcs(jsc)-mcr(jsc)))
                d_lin(ji+1,jsl,jsc) =(k_lin(ji+1,jsl,jsc) / (avan_mod*m*nvan_mod))*&
                     &  ( (frac**(-un/m))/(mc_lin(ji+1,jsc)-mcr(jsc)) ) * &
                  &  (  frac**(-un/m) -un ) ** (-m)
                d_lin(ji,jsl,jsc) = undemi * (d_lin(ji,jsl,jsc)+d_lin(ji+1,jsl,jsc))
             ELSEIF(ji.EQ.imin) THEN
                d_lin(ji,jsl,jsc) = zero
             ELSEIF(ji.EQ.imax-1) THEN
                frac=MIN(un,(mc_lin(ji,jsc)-mcr(jsc))/(mcs(jsc)-mcr(jsc)))
                d_lin(ji,jsl,jsc) =(k_lin(ji,jsl,jsc) / (avan_mod*m*nvan_mod)) * &
                     & ( (frac**(-un/m))/(mc_lin(ji,jsc)-mcr(jsc)) ) *  &
                  & (  frac**(-un/m) -un ) ** (-m)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    !! 5 Water reservoir initialisation
    !
!!$    DO jst = 1,nstm
!!$       DO ji = 1, kjpindex
!!$          mx_eau_var(ji) = mx_eau_var(ji) + soiltile(ji,jst)*&
!!$               &   dpu_max*mille*mcs(njsc(ji))
!!$       END DO
!!$    END DO
!!$    IF (check_CWRR) THEN
!!$       IF ( ANY ( ABS( mx_eau_var(:) - dpu_max*mille*mcs(njsc(:)) ) > min_sechiba ) ) THEN
!!$          ji=MAXLOC ( ABS( mx_eau_var(:) - dpu_max*mille*mcs(njsc(:)) ) , 1)
!!$          WRITE(numout, *) "Erreur formule simplifi√©e mx_eau_var ! ", mx_eau_var(ji), dpu_max*mille*mcs(njsc(ji))
!!$          WRITE(numout, *) "err = ",ABS(mx_eau_var(ji) - dpu_max*mille*mcs(njsc(ji)))
!!$          STOP 1
!!$       ENDIF
!!$    ENDIF

    mx_eau_var(:) = zero
    mx_eau_var(:) = dpu_max*mille*mcs(njsc(:)) 

    DO ji = 1,kjpindex 
       IF (vegtot(ji) .LE. zero) THEN
          mx_eau_var(ji) = mx_eau_nobio*dpu_max
          !!Et Áa veut dire quoi vegtot=0? cad frac_nobio=1? Et si 0<frac_nobio<1 ???
       ENDIF

    END DO

    ! Compute the litter humidity, shumdiag and fry
    litterhumdiag(:) = zero
    k_litt(:) = zero
    tmc_litt_mea(:) = zero
    tmc_litt_wet_mea(:) = zero
    tmc_litt_dry_mea(:) = zero
    shumdiag(:,:) = zero
    shumdiag_perma(:,:) = zero
    humtot(:) = zero
    tmc(:,:) = zero
    swi(:) = zero

    ! Loop on soil types to compute the variables (ji,jst)
    DO jst=1,nstm 
       DO ji = 1, kjpindex
          tmcs(ji,jst)=dpu_max* mille*mcs(njsc(ji))
          tmcr(ji,jst)=dpu_max* mille*mcr(njsc(ji))
       ENDDO
    ENDDO
       
    ! The total soil moisture for each soil type:
    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc(ji,jst)= dz(2,jst) * ( trois*mc(ji,1,jst)+ mc(ji,2,jst))/huit
       END DO
    ENDDO

    DO jst=1,nstm 
       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) * ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO
    ENDDO

    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc(ji,jst) = tmc(ji,jst) +  dz(nslm,jst) * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
          tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
       ENDDO
    END DO

    ! If veget has been updated before restart (with LAND USE or DGVM),
    ! tmc and mc must be modified with respect to humtot conservation.
    CALL hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg, resdist )

    ! The litter variables:
    ! level 1
    DO jst=1,nstm 
       DO ji=1,kjpindex
          tmc_litter(ji,jst) = dz(2,jst) * (trois*mc(ji,1,jst)+mc(ji,2,jst))/huit
          tmc_litter_wilt(ji,jst) = dz(2,jst) * mcw(njsc(ji)) / deux
          tmc_litter_res(ji,jst) = dz(2,jst) * mcr(njsc(ji)) / deux
          tmc_litter_field(ji,jst) = dz(2,jst) * mcf(njsc(ji)) / deux
          tmc_litter_sat(ji,jst) = dz(2,jst) * mcs(njsc(ji)) / deux
          tmc_litter_awet(ji,jst) = dz(2,jst) * mc_awet(njsc(ji)) / deux
          tmc_litter_adry(ji,jst) = dz(2,jst) * mc_adry(njsc(ji)) / deux
       ENDDO
    END DO
    ! sum from level 2 to 4
    DO jst=1,nstm 
       DO jsl=2,4
          DO ji=1,kjpindex
             tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl,jst) * & 
                  & ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
             tmc_litter_wilt(ji,jst) = tmc_litter_wilt(ji,jst) + &
                  &(dz(jsl,jst)+ dz(jsl+1,jst))*& 
                  & mcw(njsc(ji))/deux
             tmc_litter_res(ji,jst) = tmc_litter_res(ji,jst) + &
                  &(dz(jsl,jst)+ dz(jsl+1,jst))*& 
                  & mcr(njsc(ji))/deux
             tmc_litter_sat(ji,jst) = tmc_litter_sat(ji,jst) + &
                  &(dz(jsl,jst)+ dz(jsl+1,jst))* & 
                  & mcs(njsc(ji))/deux
             tmc_litter_field(ji,jst) = tmc_litter_field(ji,jst) + &
                  & (dz(jsl,jst)+ dz(jsl+1,jst))* & 
                  & mcf(njsc(ji))/deux
             tmc_litter_awet(ji,jst) = tmc_litter_awet(ji,jst) + &
                  &(dz(jsl,jst)+ dz(jsl+1,jst))* & 
                  & mc_awet(njsc(ji))/deux
             tmc_litter_adry(ji,jst) = tmc_litter_adry(ji,jst) + &
                  & (dz(jsl,jst)+ dz(jsl+1,jst))* & 
                  & mc_adry(njsc(ji))/deux
          END DO
       END DO
    END DO

    ! subsequent calcul of soil_wet_litter (tmc-tmcw)/(tmcf-tmcw)
    DO jst=1,nstm 
       DO ji=1,kjpindex
          soil_wet_litter(ji,jst)=MIN(un, MAX(zero,&
               &(tmc_litter(ji,jst)-tmc_litter_wilt(ji,jst))/&
               & (tmc_litter_field(ji,jst)-tmc_litter_wilt(ji,jst)) ))
       END DO
    ENDDO

    ! Soil wetness profiles (mc-mcw)/(mcs-mcw)
    DO jst=1,nstm 
       DO ji=1,kjpindex
          soil_wet(ji,1,jst) = MIN(un, MAX(zero,&
               &(trois*mc(ji,1,jst) + mc(ji,2,jst) - quatre*mcw(njsc(ji)))&
               & /(quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
    !!??Aurelien: a quoi sert cette ligne?
          humrelv(ji,1,jst) = zero
       ENDDO
    END DO

    DO jst=1,nstm 
       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             soil_wet(ji,jsl,jst) = MIN(un, MAX(zero,&
                  & (trois*mc(ji,jsl,jst) + & 
                  & mc(ji,jsl-1,jst) *(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst))) &
                  & + mc(ji,jsl+1,jst)*(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst))) &
                  & - quatre*mcw(njsc(ji))) / (quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
          END DO
       END DO
    END DO

    DO jst=1,nstm 
       DO ji=1,kjpindex
          soil_wet(ji,nslm,jst) = MIN(un, MAX(zero,&
               & (trois*mc(ji,nslm,jst) &
               & + mc(ji,nslm-1,jst)-quatre*mcw(njsc(ji)))/(quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
       ENDDO
    END DO

    ! Compute the grid averaged values
    DO jst=1,nstm        
       DO ji=1,kjpindex
          !
          IF ( tmc_litter(ji,jst) < tmc_litter_res(ji,jst)) THEN
             i = imin
          ELSE
             tmc_litter_ratio = (tmc_litter(ji,jst)-tmc_litter_res(ji,jst)) / &
                  & (tmc_litter_sat(ji,jst)-tmc_litter_res(ji,jst))
             i= MAX(MIN(INT((imax-imin)*tmc_litter_ratio)+imin , imax-1), imin)
          ENDIF
          ! k_litt is an averaged conductivity for saturated infiltration in the 'litter' layer
          ! This is used for reinfiltration from surface water
          k_tmp = MAX(k_lin(i,1,njsc(ji))*ks(njsc(ji)), zero)
          k_litt(ji) = k_litt(ji) + soiltile(ji,jst) * SQRT(k_tmp)
       ENDDO
    ENDDO

    DO jst=1,nstm        
       DO ji=1,kjpindex
          litterhumdiag(ji) = litterhumdiag(ji) + &
               & soil_wet_litter(ji,jst) * soiltile(ji,jst)

          tmc_litt_wet_mea(ji) =  tmc_litt_wet_mea(ji) + & 
               & tmc_litter_awet(ji,jst)* soiltile(ji,jst)

          tmc_litt_dry_mea(ji) = tmc_litt_dry_mea(ji) + &
               & tmc_litter_adry(ji,jst) * soiltile(ji,jst) 

          tmc_litt_mea(ji) = tmc_litt_mea(ji) + &
               & tmc_litter(ji,jst) * soiltile(ji,jst) 
       ENDDO
    ENDDO

    ! Caluculate frac_hydro_diag for interpolation between hydrological and diagnostic axes
    CALL hydrol_calculate_frac_hydro_diag

    ! Calculate shumdiag and shumdiag_perma from soilwet
    ! Note : shumdiag and shumdiag_perma are at diagnostic levels
    DO jst=1,nstm 
       DO jd=1,nbdl
          DO ji=1,kjpindex
             DO jsl = 1, nslm
                shumdiag(ji,jd)= shumdiag(ji,jd) + soil_wet(ji,jsl,jst)  &
                     *frac_hydro_diag(jsl,jd)* &
                     ((mcs(njsc(ji))-mcw(njsc(ji))) &
                     /(mcf(njsc(ji))-mcw(njsc(ji)))) * &
                     soiltile(ji,jst)
             
                shumdiag_perma(ji,jd)= shumdiag_perma(ji,jd)  &
                     + mc(ji,jsl,jst) *frac_hydro_diag(jsl,jd) &
                     /mcs(jst)*soiltile(ji,jst)
             END DO
             shumdiag(ji,jd) = MAX(MIN(shumdiag(ji,jd), un), zero) 
             shumdiag_perma(ji,jd) = MAX(MIN(shumdiag_perma(ji,jd), un), zero) 
          END DO
       END DO
    END DO
    
    ! Calculate soilmoist
    soilmoist(:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
             soilmoist(ji,1) = soilmoist(ji,1) + soiltile(ji,jst) * &
                  dz(2,jst) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                soilmoist(ji,jsl) = soilmoist(ji,jsl) + soiltile(ji,jst) * &
                     ( dz(jsl,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                     + dz(jsl+1,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit )
             END DO
             soilmoist(ji,nslm) = soilmoist(ji,nslm) + soiltile(ji,jst) * &
                  dz(nslm,jst) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO
    END DO


    !
    !
    DO ji=1,kjpindex
       IF ( tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji) > zero ) THEN
          drysoil_frac(ji) = un + MAX( MIN( (tmc_litt_dry_mea(ji) - tmc_litt_mea(ji)) / &
               & (tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji)), zero), - un)
       ELSE
          drysoil_frac(ji) = zero
       ENDIF
    END DO

    evap_bare_lim = zero

    IF (long_print) WRITE (numout,*) ' hydrol_var_init done '

  END SUBROUTINE hydrol_var_init


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_snow
!!
!>\BRIEF        This routine computes snow processes. 
!!
!! DESCRIPTION  :
!! - 0 initialisation
!! - 1 On vegetation
!! - 1.1 Compute snow masse
!! - 1.2 Sublimation 
!! - 1.2.1 Check that sublimation on the vegetated fraction is possible.
!! - 1.3. snow melt only if temperature positive
!! - 1.3.1 enough snow for melting or not
!! - 1.3.2 not enough snow
!! - 1.3.3 negative snow - now snow melt
!! - 1.4 Snow melts only on weight glaciers
!! - 2 On Land ice
!! - 2.1 Compute snow
!! - 2.2 Sublimation 
!! - 2.3 Snow melt only for continental ice fraction
!! - 2.3.1 If there is snow on the ice-fraction it can melt
!! - 2.4 Snow melts only on weight glaciers 
!! - 3 On other surface types - not done yet
!! - 4 computes total melt (snow and ice)
!! - 5 computes snow age on veg and ice (for albedo)
!! - 5.1 Snow age on vegetation
!! - 5.2 Snow age on ice
!! - 6 Diagnose the depth of the snow layer
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_snow

  SUBROUTINE hydrol_snow (kjpindex, dtradia, precip_rain, precip_snow , temp_sol_new, soilcap,&
       & frac_nobio, totfrac_nobio, vevapnu, vevapsno, snow, snow_age, snow_nobio, snow_nobio_age, &
       & tot_melt, snowdepth,snowmelt)

    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex      !! Domain size
    REAL(r_std), INTENT (in)                                 :: dtradia       !! Time step in seconds
    ! input fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain   !! Rainfall
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_snow   !! Snow precipitation
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: temp_sol_new  !! New soil temperature
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: soilcap       !! Soil capacity
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(in)     :: frac_nobio    !! Fraction of continental ice, lakes, ...
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: totfrac_nobio !! Total fraction of continental ice+lakes+ ...

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: tot_melt      !! Total melt from snow and ice  
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: snowmelt      !! Snow melt
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: snowdepth     !! Snow depth

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu       !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapsno      !! Snow evaporation
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow          !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: snow_age      !! Snow age
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(inout)  :: snow_nobio    !! Ice water balance
    REAL(r_std), DIMENSION (kjpindex,nnobio), INTENT(inout)  :: snow_nobio_age!! Snow age on ice, lakes, ...

    !! 0.4 Local variables

    INTEGER(i_std)                               :: ji, jv
    REAL(r_std), DIMENSION (kjpindex)             :: d_age  !! Snow age change
    REAL(r_std), DIMENSION (kjpindex)             :: xx     !! temporary
    REAL(r_std)                                   :: snowmelt_tmp !! The name says it all !

!_ ================================================================================================================================

    !
    ! for continental points
    !

    !
    !!_0 initialisation
    !
    DO jv = 1, nnobio
       DO ji=1,kjpindex
          subsnownobio(ji,jv) = zero
       ENDDO
    ENDDO
    DO ji=1,kjpindex
       subsnowveg(ji) = zero
       snowmelt(ji) = zero
       icemelt(ji) = zero
       subsinksoil(ji) = zero
       tot_melt(ji) = zero
    ENDDO
    !
    !! 1 On vegetation
    !
    DO ji=1,kjpindex
       !
    !! 1.1 Compute snow masse
       !
       snow(ji) = snow(ji) + (un - totfrac_nobio(ji))*precip_snow(ji)
       !
       !
    !! 1.2 Sublimation 
       !      Separate between vegetated and no-veget fractions 
       !      Care has to be taken as we might have sublimation from the
       !      the frac_nobio while there is no snow on the rest of the grid.
       !
       IF ( snow(ji) > snowcri ) THEN
          subsnownobio(ji,iice) = frac_nobio(ji,iice)*vevapsno(ji)
          subsnowveg(ji) = vevapsno(ji) - subsnownobio(ji,iice)
       ELSE
          ! Correction Nathalie - Juillet 2006.
          ! On doit d'abord tester s'il existe un frac_nobio!
          ! Pour le moment je ne regarde que le iice
          IF ( frac_nobio(ji,iice) .GT. min_sechiba) THEN
             subsnownobio(ji,iice) = vevapsno(ji)
             subsnowveg(ji) = zero
          ELSE 
             subsnownobio(ji,iice) = zero
             subsnowveg(ji) = vevapsno(ji)
          ENDIF
       ENDIF
       !
       !
    !! 1.2.1 Check that sublimation on the vegetated fraction is possible.
       !
       IF (subsnowveg(ji) .GT. snow(ji)) THEN
          ! What could not be sublimated goes into soil evaporation
          !         vevapnu(ji) = vevapnu(ji) + (subsnowveg(ji) - snow(ji))
          IF( (un - totfrac_nobio(ji)).GT.min_sechiba) THEN
             subsinksoil (ji) = (subsnowveg(ji) - snow(ji))/ (un - totfrac_nobio(ji))
          END IF
          ! Sublimation is thus limited to what is available
          subsnowveg(ji) = snow(ji)
          snow(ji) = zero
          vevapsno(ji) = subsnowveg(ji) + subsnownobio(ji,iice)
       ELSE
          snow(ji) = snow(ji) - subsnowveg(ji)
       ENDIF
       !
    !! 1.3. snow melt only if temperature positive
       !
       IF (temp_sol_new(ji).GT.tp_00) THEN
          !
          IF (snow(ji).GT.sneige) THEN
             !
             snowmelt(ji) = (un - frac_nobio(ji,iice))*(temp_sol_new(ji) - tp_00) * soilcap(ji) / chalfu0
             !
    !! 1.3.1 enough snow for melting or not
             !
             IF (snowmelt(ji).LT.snow(ji)) THEN
                snow(ji) = snow(ji) - snowmelt(ji)
             ELSE
                snowmelt(ji) = snow(ji)
                snow(ji) = zero
             END IF
             !
          ELSEIF (snow(ji).GE.zero) THEN
             !
    !! 1.3.2 not enough snow
             !
             snowmelt(ji) = snow(ji)
             snow(ji) = zero
          ELSE
             !
    !! 1.3.3 negative snow - now snow melt
             !
             snow(ji) = zero
             snowmelt(ji) = zero
             WRITE(numout,*) 'hydrol_snow: WARNING! snow was negative and was reset to zero. '
             !
          END IF

       ENDIF
    !! 1.4 Snow melts only on weight glaciers
       ! Ice melt only if there is more than a given mass : maxmass_snow,
       ! Ajouts Edouard Davin / Nathalie de Noblet add extra to melting
       !
       IF ( snow(ji) .GT. maxmass_snow ) THEN
          snowmelt(ji) = snowmelt(ji) + (snow(ji) - maxmass_snow)
          snow(ji) = maxmass_snow
       ENDIF
       !
    END DO
    !
    !! 2 On Land ice
    !
    DO ji=1,kjpindex
       !
    !! 2.1 Compute snow
       !
       !!??Aurelien: pkoi mettre precip_rain en dessous? We considere liquid precipitations becomes instantly snow?  
       snow_nobio(ji,iice) = snow_nobio(ji,iice) + frac_nobio(ji,iice)*precip_snow(ji) + &
            & frac_nobio(ji,iice)*precip_rain(ji)
       !
    !! 2.2 Sublimation 
       !      Was calculated before it can give us negative snow_nobio but that is OK
       !      Once it goes below a certain values (-maxmass_snow for instance) we should kill
       !      the frac_nobio(ji,iice) !
       !
       snow_nobio(ji,iice) = snow_nobio(ji,iice) - subsnownobio(ji,iice)
       !
    !! 2.3 Snow melt only for continental ice fraction
       !
       snowmelt_tmp = zero
       IF (temp_sol_new(ji) .GT. tp_00) THEN
          !
    !! 2.3.1 If there is snow on the ice-fraction it can melt
          !
          snowmelt_tmp = frac_nobio(ji,iice)*(temp_sol_new(ji) - tp_00) * soilcap(ji) / chalfu0
          !
          IF ( snowmelt_tmp .GT. snow_nobio(ji,iice) ) THEN
             snowmelt_tmp = MAX( zero, snow_nobio(ji,iice))
          ENDIF
          snowmelt(ji) = snowmelt(ji) + snowmelt_tmp
          snow_nobio(ji,iice) = snow_nobio(ji,iice) - snowmelt_tmp
          !
       ENDIF
       !
    !! 2.4 Snow melts only on weight glaciers 
       !      Ice melt only if there is more than a given mass : maxmass_snow, 
       !
       IF ( snow_nobio(ji,iice) .GT. maxmass_snow ) THEN
          icemelt(ji) = snow_nobio(ji,iice) - maxmass_snow
          snow_nobio(ji,iice) = maxmass_snow
       ENDIF
       !
    END DO

    !
    !! 3 On other surface types - not done yet
    !
    IF ( nnobio .GT. 1 ) THEN
       WRITE(numout,*) 'WE HAVE',nnobio-1,' SURFACE TYPES I DO NOT KNOW'
       WRITE(numout,*) 'CANNOT TREAT SNOW ON THESE SURFACE TYPES'
       STOP 'in hydrol_snow' 
    ENDIF

    !
    !! 4 computes total melt (snow and ice)
    !
    DO ji = 1, kjpindex
       tot_melt(ji) = icemelt(ji) + snowmelt(ji)
    ENDDO

    !
    !! 5 computes snow age on veg and ice (for albedo)
    !
    DO ji = 1, kjpindex
       !
    !! 5.1 Snow age on vegetation
       !
       IF (snow(ji) .LE. zero) THEN
          snow_age(ji) = zero
       ELSE
          snow_age(ji) =(snow_age(ji) + (un - snow_age(ji)/max_snow_age) * dtradia/one_day) &
               & * EXP(-precip_snow(ji) / snow_trans)
       ENDIF
       !
    !! 5.2 Snow age on ice
       !
       ! age of snow on ice: a little bit different because in cold regions, we really
       ! cannot negect the effect of cold temperatures on snow metamorphism any more.
       !
       IF (snow_nobio(ji,iice) .LE. zero) THEN
          snow_nobio_age(ji,iice) = zero
       ELSE
          !
          d_age(ji) = ( snow_nobio_age(ji,iice) + &
               &  (un - snow_nobio_age(ji,iice)/max_snow_age) * dtradia/one_day ) * &
               &  EXP(-precip_snow(ji) / snow_trans) - snow_nobio_age(ji,iice)
          IF (d_age(ji) .GT. min_sechiba ) THEN
             xx(ji) = MAX( tp_00 - temp_sol_new(ji), zero )
             xx(ji) = ( xx(ji) / 7._r_std ) ** 4._r_std
             d_age(ji) = d_age(ji) / (un+xx(ji))
          ENDIF
          snow_nobio_age(ji,iice) = MAX( snow_nobio_age(ji,iice) + d_age(ji), zero )
          !
       ENDIF

    ENDDO

    !
    !! 6 Diagnose the depth of the snow layer
    !

    DO ji = 1, kjpindex
       snowdepth(ji) = snow(ji) /sn_dens
    ENDDO

    IF (long_print) WRITE (numout,*) ' hydrol_snow done '

  END SUBROUTINE hydrol_snow

   
!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_canop
!!
!>\BRIEF        This routine computes canopy processes.
!!
!! DESCRIPTION  :
!! - 1 evaporation off the continents
!! - 1.1 The interception loss is take off the canopy. 
!! - 1.2 precip_rain is shared for each vegetation type
!! - 1.3 Limits the effect and sum what receives soil
!! - 1.4 swap qsintveg to the new value
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_canop

  SUBROUTINE hydrol_canop (kjpindex, precip_rain, vevapwet, veget_max, veget, qsintmax, &
       & qsintveg,precisol,tot_melt)

    ! 
    ! interface description
    !

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex    !! Domain size
    ! input fields
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: precip_rain !! Rain precipitation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: vevapwet    !! Interception loss
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: veget_max   !! max fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: veget       !! Fraction of vegetation type 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: qsintmax    !! Maximum water on vegetation for interception
    REAL(r_std), DIMENSION  (kjpindex), INTENT (in)          :: tot_melt    !! Total melt

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)       :: precisol    !! Water fallen onto the ground (throughfall)

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(inout)     :: qsintveg    !! Water on vegetation due to interception

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv
    REAL(r_std), DIMENSION (kjpindex,nvm)                    :: zqsintvegnew, water_leack !! 20-01-2016 Arsene - ADD water_leack

!_ ================================================================================================================================

    ! boucle sur les points continentaux
    ! calcul de qsintveg au pas de temps suivant
    ! par ajout du flux interception loss
    ! calcule par enerbil en fonction
    ! des calculs faits dans diffuco
    ! calcul de ce qui tombe sur le sol
    ! avec accumulation dans precisol
    ! essayer d'harmoniser le traitement du sol nu
    ! avec celui des differents types de vegetation
    ! fait si on impose qsintmax ( ,1) = 0.0
    !
    ! loop for continental subdomain
    !
    !
    !! 1 evaporation off the continents
    !
    !! 1.1 The interception loss is take off the canopy. 
    DO jv=2,nvm                                        !! 20-01-2016 Arsene - BE CAREFUL: the result can be negative (qsintveg = 0 and vevapwet >0 )
       qsintveg(:,jv) = qsintveg(:,jv) - vevapwet(:,jv)
    END DO


    !! 20-01-2016 Arsene - Add - START
    !! 1.1.bis. For moss interception is greather, but some water is loss infilt at each time step
    !!          Before new rain ? 
    !!          The part of water loss (leak/fuite) have to be dependent of time step
    !!          int = (1 - gamma) x int , with gamma = min(dt / time, 1.) to lost all water. Time: 2 week ? 

    water_leack(:,:)=zero !! 20-01-2016 Arsene - ADD

    IF ( moss_water_leack_ok .AND. ANY(.NOT.vascular(:)) ) THEN
      DO jv=2,nvm
        IF (.NOT.vascular(jv)) THEN
         WHERE ( qsintveg(:,jv).GT.min_sechiba)    !! 22-01-2016 - Arsene - Yes..q the qsintveg can be negative...
           qsintveg(:,jv) = (1. - (dt_sechiba / (moss_water_leack * one_day))) * qsintveg(:,jv)     !! 20-01-2016 Arsene - dt_sechiba = dtradia... Why ? 
           water_leack(:,jv) = dt_sechiba / (moss_water_leack * one_day) * qsintveg(:,jv)           !! 20-01-2016 Arsene - ADD
         ENDWHERE
        ENDIF
      ENDDO
    ENDIF
    !! 20-01-2016 Arsene - Add - END

    !     It is raining :
    !! 1.2 precip_rain is shared for each vegetation type
    !     sum (veget (1,nvm)) must be egal to 1-totfrac_nobio.
    !     iniveget computes veget each day
    !
    qsintveg(:,1) = zero
    DO jv=2,nvm
       IF ( ok_throughfall_by_pft ) THEN
          ! Correction Nathalie - Juin 2006 - une partie de la pluie arrivera toujours sur le sol
          ! sorte de throughfall supplementaire
          qsintveg(:,jv) = qsintveg(:,jv) + veget(:,jv) * ((1-throughfall_by_pft(jv))*precip_rain(:)) !! 20-01-2016 Arsene - change throughfall_by_pft for moss: 30% to 15%
       ELSE
          qsintveg(:,jv) = qsintveg(:,jv) + veget(:,jv) * precip_rain(:)
       ENDIF
    END DO


    !
    !! 1.3 Limits the effect and sum what receives soil
    !
    precisol(:,1)=veget_max(:,1)*precip_rain(:)
    DO jv=2,nvm
       DO ji = 1, kjpindex
          zqsintvegnew(ji,jv) = MIN (qsintveg(ji,jv),qsintmax(ji,jv))                                 !! 20-01-2016 Arsene - NB: qsintmax is increase qsintmax 
          IF ( ok_throughfall_by_pft .AND. vascular(jv)) THEN
             ! correction throughfall Nathalie - Juin 2006
             precisol(ji,jv) = (veget(ji,jv)*throughfall_by_pft(jv)*precip_rain(ji)) + &
                  qsintveg(ji,jv) - zqsintvegnew (ji,jv) + &
                  (veget_max(ji,jv) - veget(ji,jv))*precip_rain(ji) + water_leack(ji,jv)              !! 20-01-2016 Arsene - ADD the part of water leak (fuite) have to be but in precisol (not lost) to be infilt
          ELSE
             precisol(ji,jv) = qsintveg(ji,jv) - zqsintvegnew (ji,jv) + &
                  (veget_max(ji,jv) - veget(ji,jv))*precip_rain(ji) + water_leack(ji,jv)              !! 20-01-2016 Arsene - ADD the part of water leak (fuite) have to be but in precisol (not lost) to be infilt
          ENDIF
       ENDDO
    END DO
    !    
    DO jv=1,nvm
       DO ji = 1, kjpindex
          IF (vegtot(ji).GT.min_sechiba) THEN
             precisol(ji,jv) = precisol(ji,jv)+tot_melt(ji)*veget_max(ji,jv)/vegtot(ji)
          ENDIF
       ENDDO
    END DO
    !  
    !
    !! 1.4 swap qsintveg to the new value
    !
    DO jv=2,nvm
       qsintveg(:,jv) = zqsintvegnew (:,jv)
    END DO


    IF (long_print) WRITE (numout,*) ' hydrol_canop done '

  END SUBROUTINE hydrol_canop


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_vegupd
!!
!>\BRIEF        Vegetation update   
!!
!! DESCRIPTION  :
!!   The vegetation cover has changed and we need to adapt the reservoir distribution 
!!   and the distribution of plants on different soil types.
!!   You may note that this occurs after evaporation and so on have been computed. It is
!!   not a problem as a new vegetation fraction will start with humrel=0 and thus will have no
!!   evaporation. If this is not the case it should have been caught above.
!!
!! - 1 Update of vegetation is it needed?
!! - 2 calculate water mass that we have to redistribute
!! - 3 put it into reservoir of plant whose surface area has grown
!! - 4 Soil tile gestion
!! - 5 update the corresponding masks
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_vegupd

  SUBROUTINE hydrol_vegupd(kjpindex, veget, veget_max, soiltile,qsintveg,resdist)


    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                            :: kjpindex 
    ! input fields
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(in)    :: veget            !! New vegetation map
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)     :: veget_max        !! Max. fraction of vegetation type
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)   :: soiltile         !! Fraction of each soil tile (0-1, unitless)

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)  :: qsintveg         !! Water on old vegetation 
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(inout) :: resdist          !! Old vegetation map

    !! 0.4 Local variables

    INTEGER(i_std)                                 :: ji,jv,jst
    REAL(r_std), DIMENSION(kjpindex)               :: tot_corr_veg_soil

!_ ================================================================================================================================

    !! 1 If veget has been updated at last time step (with LAND USE or DGVM),
    !! tmc and mc must be modified with respect to humtot conservation.
    CALL hydrol_tmc_update ( kjpindex, veget_max, soiltile, qsintveg, resdist )

    ! Remember that it is only frac_nobio + SUM(veget_max(,:)) that is equal to 1. Thus we need vegtot
    DO ji = 1, kjpindex
       vegtot(ji) = SUM(veget_max(ji,:))
    ENDDO

    ! Compute the masks for veget
    
    mask_veget(:,:) = 0
    mask_soiltile(:,:) = 0
    
    DO jst=1,nstm
       DO ji = 1, kjpindex
          IF(soiltile(ji,jst) .GT. min_sechiba) THEN
             mask_soiltile(ji,jst) = 1
          ENDIF
       END DO
    ENDDO
          
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          IF(veget_max(ji,jv) .GT. min_sechiba) THEN
             mask_veget(ji,jv) = 1
          ENDIF
       END DO
    END DO

    ! Compute corr_veg_soil 
    corr_veg_soil(:,:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji=1,kjpindex
          ! for veget distribution used in sechiba via humrel
          IF (mask_soiltile(ji,jst).GT.0 .AND. vegtot(ji) > min_sechiba) THEN
             corr_veg_soil(ji,jv,jst)=veget_max(ji,jv)/soiltile(ji,jst)
          ENDIF
       ENDDO
    ENDDO

    IF (check_cwrr .AND. l_second_hydrol) THEN
       ! somme(soiltile * corr_veg_soil ) = 1
       tot_corr_veg_soil(:)=zero
       DO jst = 1, nstm
          DO jv = 1,nvm
             DO ji=1,kjpindex
                tot_corr_veg_soil(ji)=tot_corr_veg_soil(ji)+soiltile(ji,jst)*corr_veg_soil(ji,jv,jst)
             ENDDO
          ENDDO
       ENDDO

       DO ji=1,kjpindex
          IF ( ABS( tot_corr_veg_soil(ji) - vegtot(ji) ) > 10*EPS1 ) THEN
             WRITE(numout,*) 'corr_veg_soil SPLIT FALSE:ji=',ji,&
                  tot_corr_veg_soil(ji)
             WRITE(numout,*) 'err',ABS( tot_corr_veg_soil(ji) - vegtot(ji) )
             WRITE(numout,*) 'vegtot',vegtot(ji)
             DO jv=1,nvm
                WRITE(numout,*) 'jv,veget_max,corr_veg_soil',jv,veget_max(ji,jv),corr_veg_soil(ji,jv,:)
             END DO
             CALL ipslerr_p(3, 'hydrol_vegupd', 'Error in tot_corr_veg_soil or vegtot','','')
          ENDIF
       ENDDO
    ENDIF

    ! Tout dans cette routine est maintenant certainement obsolete (veget_max etant constant) en dehors des lignes suivantes:
    frac_bare_ns(:,:) = zero
    DO jst = 1, nstm
       DO jv = 1, nvm
          DO ji =1, kjpindex
             IF(vegtot(ji) .GT. min_sechiba) THEN
                frac_bare_ns(ji,jst) = frac_bare_ns(ji,jst) + corr_veg_soil(ji,jv,jst) * frac_bare(ji,jv) / vegtot(ji)
             ENDIF
          END DO
       ENDDO
    END DO

    IF (long_print) WRITE (numout,*) ' hydrol_vegupd done '

    RETURN
    !
  END SUBROUTINE hydrol_vegupd


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_flood
!!
!>\BRIEF        This routine computes the evolution of the surface reservoir (floodplain).  
!!
!! DESCRIPTION  :
!! - 1 Take out vevapflo from the reservoir and transfer the remaining to subsinksoil
!! - 2 Compute the total flux from floodplain floodout (transfered to routing)
!! - 3 Discriminate between precip over land and over floodplain
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_flood

  SUBROUTINE hydrol_flood (kjpindex, dtradia, vevapflo, flood_frac, flood_res, floodout)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex         !!
    ! input fields
    REAL(r_std), INTENT (in)                                 :: dtradia          !! Time step in seconds
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: flood_frac       !! Fraction of floodplains in grid box

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: floodout         !! Flux to take out from floodplains

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: flood_res        !! Floodplains reservoir estimate
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapflo         !! Evaporation over floodplains

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv           !! Indices
    REAL(r_std), DIMENSION (kjpindex)                        :: temp             !! 

!_ ================================================================================================================================
    !- 
    !! 1 Take out vevapflo from the reservoir and transfer the remaining to subsinksoil 
    !-
    DO ji = 1,kjpindex
       temp(ji) = MIN(flood_res(ji), vevapflo(ji))
    ENDDO
    DO ji = 1,kjpindex
       flood_res(ji) = flood_res(ji) - temp(ji)
       subsinksoil(ji) = subsinksoil(ji) + vevapflo(ji) - temp(ji)
       vevapflo(ji) = temp(ji)
    ENDDO

    !- 
    !! 2 Compute the total flux from floodplain floodout (transfered to routing) 
    !-
    DO ji = 1,kjpindex
       floodout(ji) = vevapflo(ji) - flood_frac(ji) * SUM(precisol(ji,:))
    ENDDO

    !-
    !! 3 Discriminate between precip over land and over floodplain
    !-
    DO jv=1, nvm
       DO ji = 1,kjpindex
  precisol(ji,jv) = precisol(ji,jv) * (1 - flood_frac(ji))
ENDDO
ENDDO 

IF (long_print) WRITE (numout,*) ' hydrol_flood done'

END SUBROUTINE hydrol_flood


!! ================================================================================================================================
!! SUBROUTINE 	: hydrol_soil
!!
!>\BRIEF        This routine computes soil processes with CWRR scheme.
!!
!! DESCRIPTION  :
!! - 0 Arrays initialisation
!! -   for each soil type
!! - 1 We compare water2infilt and water2extract to keep only difference
!! - 1.1 add to the first layer
!! - 1.2 filling layers
!! - 2 Before diffusion scheme
!! - 2.1 Some initialisation necessary for the diffusion scheme to work
!! - 2.2 coefficients are computed for the profile of mc before infiltration:
!! - 3 The infiltration is computed 
!! - 4 Coefficient are recomputed for the profile of mc after infiltration:
!! - 5 Prepar the diffusion scheme
!! - 5.1 Set the values for diffusion scheme
!! - 5.2 verifications for a good soil humidity
!! - 5.3 compute matrix coefficients
!! - 6 Resolve diffusion scheme
!! - 6.1 solve equations assuming atmosphere limiting
!! - 6.2 check if really atmosphere limiting
!! - 6.3 Reset the coefficient for diffusion (only used if resolv(ji) = .TRUE.)
!! - 6.4 resolve the equations with new boundary conditions if necessary
!! - 7 close the water balance
!! - 7.1 compute dr_ns with the bottom boundary condition 
!! - 7.2 compute total soil moisture content
!! - 7.3 deduction of upper flux from soil moisture variation and bottom flux
!! - 7.4 deduction of ae_ns and ru_ns:
!! - 8 Special treatment for the unstable cases
!! - 9 Then compute the temporary surface water and correct the outgoing runoff
!! - 10 smooth again
!! - 11 Optional computation of the fluxes 
!! - 12 We make some useful output
!! - 13 before closing the soil water, we check the water balance of soil
!! - 14 sum 3d variables into 2d variables
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil

SUBROUTINE hydrol_soil (kjpindex, dtradia, veget_max, soiltile, njsc, reinf_slope, &
& transpir, vevapnu, evapot, evapot_penm, runoff, drainage, &
& returnflow, reinfiltration, irrigation, &
& tot_melt, evap_bare_lim, shumdiag, shumdiag_perma,&
& k_litt, litterhumdiag, humrel,vegstress, drysoil_frac, &
& stempdiag,snow, &
& snowdz, & !pss:+
& fsat) !pss:-  !! Arsene 28-01-2016 - REMOVE drunoff_tot because never user and bug in sechiba_output.f90
! 
! interface description

!! 0. Variable and parameter declaration

!! 0.1 Input variables

! input scalar 
INTEGER(i_std), INTENT(in)                               :: kjpindex 
! input fields
REAL(r_std), INTENT (in)                                 :: dtradia          !! Time step [s]
REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: veget_max        !! Map of max vegetation types [-]
INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile         !! Fraction of each soil tile (0-1, unitless)
REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(in)        :: transpir         !! Transpiration [mm]
REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: reinf_slope      !! Slope coef
REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: returnflow       !! Water returning to the deep reservoir [mm]
REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: reinfiltration   !! Water returning to the top of the soil [mm]
REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: irrigation       !! Irrigation [mm]
REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot           !! Potential evaporation [mm]
REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot_penm      !! Potential evaporation Penman [mm] 
REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_melt         !! Total melt from snow and ice [mm]
REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)       :: stempdiag        !! Diagnostic temp profile from thermosoil
REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: snow             !! Snow mass [Kg/m^2]
REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(in)       :: snowdz           !! Snow depth [m]
!pss:+
REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: fsat             !! fraction of saturation soil 
!pss:-

!! 0.2 Output variables

REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: runoff           !! complete runoff [mm]
REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: drainage         !! complete drainage [mm]
REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: evap_bare_lim    !! limitation of bare soil evaporation on each 
									 !! soil column [mm] 
REAL(r_std), DIMENSION (kjpindex,nbdl), INTENT (out)     :: shumdiag         !! Relative soil moisture
REAL(r_std), DIMENSION (kjpindex,nbdl), INTENT (out)     :: shumdiag_perma   !! Percent of porosity filled with water (mc/mcs) used for the thermal computations
REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: k_litt           !! Litter cond.
REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: litterhumdiag    !! Litter humidity
REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(out)      :: vegstress        !! Veg. moisture stress (only for vegetation 
									 !! growth) 
REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: drysoil_frac     !! Function of the litter humidity, that will be 
!pss:+
!REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: drunoff_tot      !! Dunne runoff !! Arsene 28-01-2016 - REMOVE because never user and bug in sechiba_output.f90
!pss:-

!! 0.3 Modified variables

REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu          !! 
REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (inout)    :: humrel           !! Relative humidity [0-1, dimensionless]

!! 0.4 Local variables

INTEGER(i_std)                                 :: ji, jv, jsl, jst           !! indices
REAL(r_std), PARAMETER                         :: frac_mcs = 0.66            !! temporary depth
!!??Aurelien: frac_ms inutile?
REAL(r_std)                                    :: us_tmp                     !! temporary stress
REAL(r_std), DIMENSION(kjpindex)               :: temp                       !! temporary value for fluxes
REAL(r_std), DIMENSION(kjpindex)               :: tmcold, tmcint             !!
REAL(r_std), DIMENSION(kjpindex,nslm,nstm)     :: moderwilt                  !!
REAL(r_std), DIMENSION(kjpindex,nslm)          :: mcint                      !! To save mc values for future use
LOGICAL, DIMENSION(kjpindex)                   :: is_under_mcr               !! Allows under residual soil moisture due to evap 
LOGICAL, DIMENSION(kjpindex)                   :: is_over_mcs                !! Allows over saturated soil moisture due to 
									 !! returnflow 
REAL(r_std), DIMENSION(kjpindex)               :: sum_rootsink               !! Sum of the root sink
REAL(r_std), DIMENSION(kjpindex)               :: deltahum,diff              !!
LOGICAL(r_std), DIMENSION(kjpindex)            :: test                       !!
REAL(r_std), DIMENSION(kjpindex)               :: tsink                      !!
REAL(r_std), DIMENSION(kjpindex)               :: water2extract              !! Temporary variable [mm]
REAL(r_std), DIMENSION(kjpindex)               :: returnflow_soil            !! Water from the routing back to the bottom of 
									 !! the soil [mm] 
REAL(r_std), DIMENSION(kjpindex)               :: reinfiltration_soil        !! Water from the routing back to the top of the 
									 !! soil [mm] 
REAL(r_std), DIMENSION(kjpindex)               :: irrigation_soil            !! Water from irrigation returning to soil 
									 !! moisture for each soil type [mm] 
REAL(r_std), DIMENSION(kjpindex)               :: flux_infilt                !!
REAL(r_std), DIMENSION(kjpindex)               :: flux_bottom                !! Flux at bottom of the soil column
REAL(r_std), DIMENSION(kjpindex)               :: flux_top                   !! Flux at top of the soil column
LOGICAL                                        :: error=.FALSE.              !! If true, exit in the end of subroutine

REAL(r_std)                                    :: x_moss, x_tot               ! Fraction of moss / total veget fract in a soil tile  !! 20-01-2016 Arsene - Add

!_ ================================================================================================================================

!! 0 Arrays initialisation

returnflow_soil(:) = zero
reinfiltration_soil(:) = zero
irrigation_soil(:) = zero
qflux(:,:,:) = zero
is_under_mcr(:) = .FALSE.
is_over_mcs(:) = .FALSE.
flux_infilt(:) = zero
IF (ok_freeze_cwrr) THEN
kk(:,:,:)=zero
kk_moy(:,:)=zero
ENDIF

IF (ok_freeze_cwrr) THEN
!
! 1.1 Calculate the temperature at hydrological levels
!
CALL hydrol_calculate_temp_hydro(kjpindex, stempdiag, snow,snowdz)
ENDIF
!
! split 2d variables to 3d variables, per soil tile
!
CALL hydrol_split_soil (kjpindex, veget_max, soiltile, vevapnu, transpir, humrel, evap_bare_lim)
!
! Common variables
!
DO ji=1,kjpindex
IF(vegtot(ji).GT.min_sechiba) THEN
  returnflow_soil(ji) = zero
  reinfiltration_soil(ji) = (returnflow(ji) + reinfiltration(ji))/vegtot(ji)
  irrigation_soil(ji) = irrigation(ji)/vegtot(ji)
ELSE
  returnflow_soil(ji) = zero
  reinfiltration_soil(ji) = zero
  irrigation_soil(ji) = zero
ENDIF
ENDDO

!
!!_  for each soil tile
!
DO jst = 1,nstm
!
!- Compute the sum of the sinks for future check-up
sum_rootsink(:)=SUM(rootsink(:,:,jst),dim=2)
DO ji=1,kjpindex
  tsink(ji) = sum_rootsink(ji)+MAX(ae_ns(ji,jst),zero)+subsinksoil(ji)
ENDDO
!
! The total moisture content (including water2infilt) is saved for balance checks at the end
tmcold(:) = tmc(:,jst)

!- The value of mc is kept in mcint, used in the flux computation after diffusion:
DO jsl = 1, nslm
  DO ji = 1, kjpindex
     mcint(ji,jsl) = mask_soiltile(ji,jst) * mc(ji,jsl,jst)
  ENDDO
ENDDO

DO ji = 1, kjpindex
  tmcint(ji) = dz(2,jst) * ( trois*mcint(ji,1) + mcint(ji,2) )/huit 
ENDDO

DO jsl = 2,nslm-1
  DO ji = 1, kjpindex
     tmcint(ji) = tmcint(ji) + dz(jsl,jst) &
	  & * (trois*mcint(ji,jsl)+mcint(ji,jsl-1))/huit &
	  & + dz(jsl+1,jst) * (trois*mcint(ji,jsl)+mcint(ji,jsl+1))/huit
  ENDDO
ENDDO

DO ji = 1, kjpindex
  tmcint(ji) = tmcint(ji) + dz(nslm,jst) &
       & * (trois * mcint(ji,nslm) + mcint(ji,nslm-1))/huit
ENDDO
!! 1 Compare water2infilt and water2extract to keep only difference

! The bare soil evaporation is substracted to the soil moisture profile and first to the water available:  
DO ji = 1, kjpindex
  water2extract(ji) = MIN(water2infilt(ji,jst) + irrigation_soil(ji) + reinfiltration_soil(ji), &
       & MAX(ae_ns(ji,jst),zero) + subsinksoil(ji))
ENDDO

! First we substract from the surface
DO ji = 1, kjpindex
  water2infilt(ji,jst) = water2infilt(ji,jst) + irrigation_soil(ji) + reinfiltration_soil(ji) - water2extract(ji)
ENDDO       
     
! Then we update the water to extract from the soil
DO ji = 1, kjpindex
  water2extract(ji) =  MAX(ae_ns(ji,jst),zero) + subsinksoil(ji) - water2extract(ji)
ENDDO

!! We add and substract components to the soil
!! by filling the layers from top to bottom (same for reinfiltration) - this is done by smooth below
!! 1.1 add to the first layer
DO ji = 1, kjpindex
  mc(ji,1,jst) = mc(ji,1,jst) - water2extract(ji) * deux / dz(2,jst)
ENDDO

!!??Aurelien: here, the first layer can be oversaturated, thats why we need hydrol_soil_smooth
!! 1.2 filling layers
CALL hydrol_soil_smooth(kjpindex,jst, njsc, is_under_mcr, is_over_mcs)


!! 2 Before diffusion scheme

!! 2.1 Some initialisation necessary for the diffusion scheme to work

DO ji = 1, kjpindex
  !- We correct rootsink for first two layers so that it is not too low in the first layer
  v1(ji,jst) = dz(2,jst)/huit * (trois * mc(ji,1,jst)+ mc(ji,2,jst))
  rootsink(ji,2,jst) = rootsink(ji,2,jst) + MAX(rootsink(ji,1,jst)-v1(ji,jst), zero) 
  rootsink(ji,1,jst) = MIN(rootsink(ji,1,jst),v1(ji,jst))
ENDDO
!- Flux at top of the soil column is zero
flux_top(:) = zero

DO ji = 1, kjpindex
  ! Initialise the flux to be infiltrated 
  flux_infilt(ji) = water2infilt(ji,jst) + precisol_ns(ji,jst) 
          !- The incoming flux is also first dedicated to fill the soil up to mcr (in case needed)
ENDDO


     !! 2.1.2 Treat problem with to low soil moisture content in soil column : is_under_mcr
       DO ji = 1, kjpindex
          IF (is_under_mcr(ji)) THEN
             ! Calculate minimum mc
             temp(ji)=mc(ji,1,jst)
             DO jsl=2,nslm
                IF (mc(ji,jsl,jst) .LT. temp(ji)) THEN
                   temp(ji)=mc(ji,jsl,jst)
                END IF
             END DO
      
             ! Calculate water deficit
             temp(ji)=mcr(njsc(ji))-temp(ji)

             ! Add to each soil layer
             DO jsl=1,nslm
                mc(ji,jsl,jst) = mc(ji,jsl,jst) + temp(ji)
             END DO

             ! Remove from flux_infilt
             flux_infilt(ji)=flux_infilt(ji)-temp(ji)*(dpu_max*mille)
             IF (flux_infilt(ji) .LT. 0.) THEN
                WRITE(numout,*) 'BIG PROBLEM at grid cell ',ji
                WRITE(numout,*) 'No available water in soil column!!!'
                flux_infilt(ji) = 0.
             END IF
          END IF
       END DO


    !! 2.2 Coefficients are computed for the profile of mc before infiltration
       CALL hydrol_soil_coef(kjpindex,jst,njsc)

    !! 3 The infiltration is computed 
       CALL hydrol_soil_infilt(kjpindex, jst, dtradia, njsc, flux_infilt)

    !! 4 Coefficient are recomputed for the profile of mc after infiltration
       CALL hydrol_soil_coef(kjpindex,jst,njsc)

    !! 5 Prepar the diffusion scheme
 
    !! 5.1 Set the values for diffusion scheme
       CALL hydrol_soil_setup(kjpindex,jst,dtradia)

    !! 5.2 verifications for a good soil humidity
    ! We only run the scheme in case we are not under mcr after precip and various reinfiltrations and not over mcs
       resolv(:) = (mask_soiltile(:,jst) .GT. 0) .AND. & 
               & (.NOT. is_under_mcr(:)) .AND. (.NOT. is_over_mcs(:))    

       ! In oversaturated case, we first take the evaporation from the outgoing water (the rest will be taken from the soil)
       sum_rootsink(:)=SUM(rootsink(:,:,jst),dim=2)
       WHERE (is_over_mcs(:))
          mc(:,1,jst) = mc(:,1,jst) - flux_top(:) * deux / dz(2,jst)
       ENDWHERE
       DO jsl=1, nslm
          WHERE (is_over_mcs(:))
             mc(:,jsl,jst) = mc(:,jsl,jst) - sum_rootsink(:) / (dpu_max*mille) 
          ENDWHERE
       ENDDO
       ! In under residual case, we equally spread the transpiration over the layers
       DO jsl = 1, nslm
          WHERE (is_under_mcr(:))
             mc(:,jsl,jst) = mc(:,jsl,jst) - sum_rootsink(:) / (dpu_max*mille) 
          ENDWHERE
       ENDDO

       IF (ok_freeze_cwrr) THEN
          DO ji =1, kjpindex
             DO jsl = 1, nslm
                mcl(ji,jsl,:)= MIN(mc(ji,jsl,:),mcr(njsc(ji))+(1-profil_froz_hydro_ns(ji,jsl, :))*(mc(ji,jsl,:)-mcr(njsc(ji))))
             ENDDO
          ENDDO
       ELSE
          mcl(:,:,:)=mc(:,:,:)
       ENDIF

    !! 5.3 compute matrix coefficients
       !- First layer
       DO ji = 1, kjpindex
          tmat(ji,1,1) = zero
          tmat(ji,1,2) = f(ji,1)
          tmat(ji,1,3) = g1(ji,1)
          rhs(ji,1)    = fp(ji,1) * mcl(ji,1,jst) + gp(ji,1)*mcl(ji,2,jst) &
               &  - flux_top(ji) - (b(ji,1)+b(ji,2))/deux *(dtradia/one_day) - rootsink(ji,1,jst)
       ENDDO

       !- soil body
       DO jsl=2, nslm-1
          DO ji = 1, kjpindex
             tmat(ji,jsl,1) = e(ji,jsl)
             tmat(ji,jsl,2) = f(ji,jsl)
             tmat(ji,jsl,3) = g1(ji,jsl)
             rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
                  & +  gp(ji,jsl) * mcl(ji,jsl+1,jst) & 
                  & + (b(ji,jsl-1) - b(ji,jsl+1)) * (dtradia/one_day) / deux & 
                  & - rootsink(ji,jsl,jst) 
          ENDDO
       ENDDO
       
       !- Last layer
       DO ji = 1, kjpindex
          jsl=nslm
          tmat(ji,jsl,1) = e(ji,jsl)
          tmat(ji,jsl,2) = f(ji,jsl)
          tmat(ji,jsl,3) = zero
          rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
               & + (b(ji,jsl-1) + b(ji,jsl)*(un-deux*free_drain_coef(ji,jst))) * (dtradia/one_day) / deux &
               & - rootsink(ji,jsl,jst)
       ENDDO

       !- store the equations in case needed again
       DO jsl=1,nslm
          DO ji = 1, kjpindex
             srhs(ji,jsl) = rhs(ji,jsl)
             stmat(ji,jsl,1) = tmat(ji,jsl,1)
             stmat(ji,jsl,2) = tmat(ji,jsl,2)
             stmat(ji,jsl,3) = tmat(ji,jsl,3) 
          ENDDO
       ENDDO

    !! 6 Resolve diffusion scheme

    !! 6.1 solve diffusion equations (update mcl)
       CALL hydrol_soil_tridiag(kjpindex,jst)

    !! 6.2 Correct bad moisture content due to numerical errors before water balance
       DO ji = 1, kjpindex
          DO jsl = 1,nslm
             IF (mcl(ji,jsl,jst) .LT. mcr(njsc(ji)) .AND. mcl(ji,jsl,jst) .NE. 0.) THEN
                mcl(ji,jsl,jst)=mcr(njsc(ji))
             END IF
             IF (mcl(ji,jsl,jst) .GT. mcs(njsc(ji)) .AND. mcl(ji,jsl,jst) .NE. 0.) THEN
                mcl(ji,jsl,jst)=mcs(njsc(ji))
             END IF
          ENDDO
       ENDDO

       ! Calculation of total soil moisture content (liquid + frozen)
       IF (ok_freeze_cwrr) THEN
          DO ji =1, kjpindex
             DO jsl = 1, nslm
                mc(ji,jsl,:)=MAX(mcl(ji,jsl,:), mcl(ji,jsl,:)+profil_froz_hydro_ns(ji,jsl,:)*(mc(ji,jsl,:)-mcr(njsc(ji))))
             ENDDO
          ENDDO
       ELSE
          mc(:,:,:)=mcl(:,:,:)
       ENDIF

    !! 7 close the water balance
    !! 7.1 compute dr_ns with the bottom boundary condition 
    !initialize qflux at bottom of diffusion and avoid over saturated or under residual soil moisture 

       DO ji = 1, kjpindex
          dr_ns(ji,jst)=zero
          jsl=nslm
          IF (.NOT. is_under_mcr(ji)) THEN
             dr_ns(ji,jst) = mask_soiltile(ji,jst)*k(ji,jsl)*free_drain_coef(ji,jst) * (dtradia/one_day)
          ENDIF
       ENDDO
          

    !! 7.2 compute total soil moisture content
       DO ji = 1, kjpindex
          tmc(ji,jst) = dz(2,jst) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit 
       ENDDO
          
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) &
                  & * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO

       DO ji = 1, kjpindex
          tmc(ji,jst) = tmc(ji,jst) + dz(nslm,jst) &
               & * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO

       DO ji = 1, kjpindex
    !! 7.3 deduction of upper flux from soil moisture variation and bottom flux
          qflux00(ji,jst) = mask_soiltile(ji,jst) * &
               & (MIN(tmcs(ji,jst),tmc(ji,jst))-tmcint(ji)+SUM(rootsink(ji,:,jst))+dr_ns(ji,jst)-returnflow_soil(ji))

    !! 7.4 deduction of ae_ns and ru_ns:
    ! ae_ns+ru_ns=precisol_ns+irrigation-q0  
          ru_ns(ji,jst) = (precisol_ns(ji,jst)  &
              & +water2infilt(ji,jst)-water2extract(ji)-qflux00(ji,jst)) * mask_soiltile(ji,jst)
          ae_ns(ji,jst) =ae_ns(ji,jst) + subsinksoil(ji)
!! Arsene 14-01-2016 - here the potentital runoff is comput. But one part can be conserve in water2infilt (with rein_slope)...

       ENDDO

    !! 8 Special treatment for the unstable cases
    !! 8.1 Negative runoff

       IF (long_print) THEN
          DO ji = 1, kjpindex
             IF (ru_ns(ji,jst) .LT. zero) THEN
                WRITE (numout,*) 'Negative runoff corrected', ji,jst,ru_ns(ji,jst), mc(ji,1,jst), tmc(ji,jst)
             ENDIF
          ENDDO
       ENDIF

      DO ji = 1, kjpindex
             IF (ru_ns(ji,jst) .LT. zero) THEN
                ! Available water in first layer
                temp(ji)=MAX(((mc(ji,1,jst)-mcr(njsc(ji))) * dz(2,jst) / deux), zero)

                ! Calculate and add maximum water to runoff
                temp(ji)=MIN(temp(ji),ABS(ru_ns(ji,jst)))
                ru_ns(ji,jst)=ru_ns(ji,jst)+temp(ji)

                ! Update water balance
                qflux00(ji,jst) = qflux00(ji,jst) - temp(ji)
                mc(ji,1,jst)=mc(ji,1,jst)-(temp(ji) * deux / dz(2,jst))
             END IF

             ! If still negative runoff, take water from deeper layers
             DO jsl = 2, nslm-1
                IF (ru_ns(ji,jst) .LT. zero) THEN
                   temp(ji)= MAX((mc(ji,jsl,jst)-mcr(njsc(ji)))*(dz(jsl,jst)+dz(jsl+1,jst)) / deux, zero)
                   temp(ji)=MIN(temp(ji),ABS(ru_ns(ji,jst)))
                   ru_ns(ji,jst)=ru_ns(ji,jst)+temp(ji)
                   qflux00(ji,jst) = qflux00(ji,jst) - temp(ji)
                   mc(ji,jsl,jst)=mc(ji,jsl,jst)-temp(ji) * deux / (dz(jsl,jst)+dz(jsl+1,jst))
                END IF
             ENDDO

             ! Last layer
             IF (ru_ns(ji,jst) .LT. zero) THEN
                temp(ji)=MAX((mc(ji,nslm,jst)-mcr(njsc(ji))) * dz(nslm,jst) / deux, zero)
                temp(ji)=MIN(temp(ji),ABS(ru_ns(ji,jst)))
                ru_ns(ji,jst)=ru_ns(ji,jst)+temp(ji)
                qflux00(ji,jst) = qflux00(ji,jst) - temp(ji)
                mc(ji,nslm,jst)=mc(ji,nslm,jst)-temp(ji) * deux / dz(nslm,jst)
             END IF

             ! If still negative runoff, take water from bottom drainage
             IF (ru_ns(ji,jst) .LT. zero) THEN
                IF (long_print)  WRITE (numout,*) 'runoff and drainage before correction',&
                     ru_ns(ji,jst),dr_ns(ji,jst)
                dr_ns(ji,jst)=dr_ns(ji,jst)+ru_ns(ji,jst)
                qflux00(ji,jst) = qflux00(ji,jst) + ru_ns(ji,jst)
                ru_ns(ji,jst)= 0.
             END IF
             
          ENDDO

    !! 8.2 Compute total soil moisture content
       DO ji = 1, kjpindex
          tmc(ji,jst) = dz(2,jst) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
       ENDDO

       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) &
                  & * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO

       DO ji = 1, kjpindex
          tmc(ji,jst) = tmc(ji,jst) + dz(nslm,jst) &
               & * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO

    !! 8.3 Avoid under-precision value for the 3 outward flux
       DO ji = 1, kjpindex  
          IF (ABS(ae_ns(ji,jst)).LT.min_sechiba) THEN
             ae_ns(ji,jst) = zero
          ENDIF

          IF(ABS(ru_ns(ji,jst)).LT.min_sechiba) THEN
             ru_ns(ji,jst) = zero
          ENDIF
       ENDDO

    !! 9 Compute temporary surface water and extract from runoff
       IF ( .NOT. doponds ) THEN 

         DO ji = 1, kjpindex

           IF ( reinf_slope_moss_ok .AND. ANY(.NOT.vascular(:)) ) THEN !! Arsene 23-02-2016 - change runoff is moss

             !! 9.1 Compute fraction of moss in each tile soil  !! Arsene 20-01-2016 - Add for moss: lower runoff
             x_moss = zero
             x_tot = zero
             DO jv = 1,nvm
               IF ( pref_soil_veg(jv).EQ.jst ) THEN
                 x_tot = veget_max(ji,jv) + x_tot
                 IF ( .NOT.vascular(jv) .AND. jv.NE.1 ) THEN 
                   x_moss = veget_max(ji,jv) + x_moss
                 ENDIF
               ENDIF
             ENDDO
             IF (x_tot .LT. min_sechiba ) THEN
               x_moss = zero
             ELSE
               x_moss = x_moss / x_tot
             ENDIF

             !! 9.2 Limit the runoff. reinf_slope topography dependent (if less land slope, less runof)
             !!     Take care: the water in modified for moss and grasses. Maybe a new tile soil is needed ?
             !!     More land slope -> less reinf_slope. For mosses: increase reinf_slope to decrease runoff

! Fisrt test:
!             water2infilt(ji,jst) = MIN(1.,(reinf_slope(ji) * (1. + x_moss*reinf_slope_moss_coef))) * ru_ns(ji,jst)

!Second test :
!             water2infilt(ji,jst) = MIN(50., ( 1. - (1.-MIN(0.99,(1-x_moss)*reinf_slope(ji)+x_moss*reinf_slope_moss)) &  ! With max limit for water2infilt
!                                         * dt_sechiba / (moss_water_runoff * one_day)) * ru_ns(ji,jst))                  ! ==> arround 0 runoff

             water2infilt(ji,jst) = ( 1. - (1.-((1-x_moss)*reinf_slope(ji)+x_moss*reinf_slope_moss)) & !! We need the (1 - (1 - x) * dt/dt ) to be adaptable to timestep changes
                                         * (dtradia / (one_day/48.))) * ru_ns(ji,jst)


           ELSE    !! Arsene 23-02-2016 - if reinf_slope_moss_ok --> original

!             water2infilt(ji,jst) = reinf_slope(ji) * ru_ns(ji,jst)  !! Arsene 23-02-2016 : need to change to be adaptable to timestep changes
             water2infilt(ji,jst) = ( 1. - (1. - reinf_slope(ji)) * (dtradia / (one_day/48.)) ) * ru_ns(ji,jst)

           ENDIF   !! Arsene 23-02-2016 - if reinf_slope_moss_ok  
         ENDDO

       ELSE
          DO ji = 1, kjpindex           
             water2infilt(ji,jst) = zero
          ENDDO
       ENDIF

      !
       DO ji = 1, kjpindex           
          ru_ns(ji,jst) = MAX(0., ru_ns(ji,jst) - water2infilt(ji,jst)) !! Arsene 09-02-2016 - Secure runoff to don't be negative (with new add)
       END DO

    !! 10 Smooth again
       ! Probably not necessary but harmless (Aurelien)
       CALL hydrol_soil_smooth(kjpindex, jst, njsc, is_under_mcr, is_over_mcs)

    !! 11 Optional computation of the fluxes 
       IF ( check_cwrr ) THEN
          CALL hydrol_soil_flux(kjpindex,jst,mcint,returnflow_soil)
       ENDIF

    !! 12 We make some useful output
    !- Total soil moisture, soil moisture at litter levels, soil wetness...

       !-total soil moisture:
       DO ji=1,kjpindex
          tmc(ji,jst)= dz(2,jst) * (trois*mc(ji,1,jst) + mc(ji,2,jst))/huit
       END DO

       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) * ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO

       DO ji=1,kjpindex
          tmc(ji,jst) = tmc(ji,jst) +  dz(nslm,jst) * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
          tmc(ji,jst) = tmc(ji,jst) + water2infilt(ji,jst)
       END DO

       ! The litter is the 4 top levels of the soil
       ! Compute various field of soil moisture for the litter (used for stomate and for albedo)

       DO ji=1,kjpindex
          tmc_litter(ji,jst) = dz(2,jst) * ( trois*mc(ji,1,jst)+ mc(ji,2,jst))/huit
          tmc_litter(ji,jst) = tmc_litter(ji,jst) 
       END DO

       ! sum from level 1 to 4

       DO jsl=2,4
          DO ji=1,kjpindex
             tmc_litter(ji,jst) = tmc_litter(ji,jst) + dz(jsl,jst) * & 
                  & ( trois*mc(ji,jsl,jst) + mc(ji,jsl-1,jst))/huit &
                  & + dz(jsl+1,jst)*(trois*mc(ji,jsl,jst) + mc(ji,jsl+1,jst))/huit
          END DO
       END DO

       ! subsequent calcul of soil_wet_litter (tmc-tmcw)/(tmcf-tmcw)

       DO ji=1,kjpindex
          soil_wet_litter(ji,jst) = MIN(un, MAX(zero,&
               & (tmc_litter(ji,jst)-tmc_litter_wilt(ji,jst)) / &
               & (tmc_litter_field(ji,jst)-tmc_litter_wilt(ji,jst)) ))
       END DO

       ! Soil wetness profiles (mc-mcw)/(mcs-mcw)
       ! soil_wet is the ratio of soil moisture to available soil moisture for plant
       ! (ie soil moisture at saturation minus soil moisture at wilting point).

       DO ji=1,kjpindex
          soil_wet(ji,1,jst) = MIN(un, MAX(zero,&
               & (trois*mc(ji,1,jst) + mc(ji,2,jst) - quatre*mcw(njsc(ji)))&
               & /(quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
          humrelv(ji,1,jst) = zero
       END DO

       DO jsl=2,nslm-1
          DO ji=1,kjpindex
             soil_wet(ji,jsl,jst) = MIN(un, MAX(zero,&
                  & (trois*mc(ji,jsl,jst) + & 
                  & mc(ji,jsl-1,jst) *(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst))) &
                  & + mc(ji,jsl+1,jst)*(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst))) &
                  & - quatre*mcw(njsc(ji))) / (quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
          END DO
       END DO

       DO ji=1,kjpindex
          soil_wet(ji,nslm,jst) = MIN(un, MAX(zero,&
               & (trois*mc(ji,nslm,jst) &
               & + mc(ji,nslm-1,jst)-quatre*mcw(njsc(ji)))/(quatre*(mcs(njsc(ji))-mcw(njsc(ji)))) ))
       END DO

       ! Compute the moderation of transpiration due to wilting point
       ! moderwilt is a factor which is zero if soil moisture is below the wilting point
       ! and is un if soil moisture is above the wilting point.

       DO jsl=1,nslm
          DO ji=1,kjpindex
             moderwilt(ji,jsl,jst) = INT( MAX(soil_wet(ji,jsl,jst), zero) + un - min_sechiba )
          END DO
       END DO

       ! Compute the new humrelv to use in sechiba:
       ! loop on each vegetation type
       humrelv(:,1,jst) = zero   

       ! calcul of us for each layer and vegetation type.
       DO jv = 2,nvm
          DO ji=1,kjpindex
             !- Here we make the assumption that roots do not take water from the 1st layer. 
             !- Comment the us=0 if you want to change this.
!             us(ji,jv,jst,1) = moderwilt(ji,1,jst)*MIN(un,((trois*mc(ji,1,jst) + mc(ji,2,jst)) &
!                  & /(quatre*mcs(jst)*pcent(jst))) )* (un-EXP(-humcste(jv)*dz(2,jst)/mille/deux)) &
!                  & /(un-EXP(-humcste(jv)*zz(nslm,jst)/mille))
             us(ji,jv,jst,1) = zero
             humrelv(ji,jv,jst) = MAX(us(ji,jv,jst,1),zero)
          END DO
       ENDDO


       DO jsl = 2,nslm-1
          DO jv = 2, nvm
             DO ji=1,kjpindex
                ! Influence of freezing in the soil on the soil moisture
                IF (ok_freeze_cwrr) THEN
                   us_tmp= (trois*(1-profil_froz_hydro_ns(ji,jsl, jst))*mc(ji,jsl,jst)+ &
                        (1-profil_froz_hydro_ns(ji,jsl-1, jst))*mc(ji,jsl-1,jst) &
                        *(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))+ &
                        (1-profil_froz_hydro_ns(ji,jsl+1, jst))*mc(ji,jsl+1,jst) &
                        *(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))) &
                        /(quatre*mcs(njsc(ji))*pcent(njsc(ji)))

                ELSE
                   us_tmp= (trois*mc(ji,jsl,jst) + &
                        mc(ji,jsl-1,jst)*(dz(jsl,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))+ &
                        mc(ji,jsl+1,jst)*(dz(jsl+1,jst)/(dz(jsl,jst)+dz(jsl+1,jst)))) &
                        /(quatre*mcs(njsc(ji))*pcent(njsc(ji)))
                ENDIF

                ! us_tmp should not be negative
                !!??Aurelien: Useless, here we are sure mc between mcr and mcs
                us_tmp=MAX(us_tmp, zero)
                ! us is computed with a SQRT in order for it to grow more rapidly with soil moisture.
                ! it is not essential
                us(ji,jv,jst,jsl) = moderwilt(ji,jsl,jst) * MIN( un, SQRT(us_tmp) ) * nroot(jv,jst,jsl)
                !
                us(ji,jv,jst,jsl) = MAX(us (ji,jv,jst,jsl), zero)
                
                humrelv(ji,jv,jst) = MAX((humrelv(ji,jv,jst) + us(ji,jv,jst,jsl)),zero)
             END DO
          END DO
       ENDDO


       DO jv = 2, nvm
          DO ji=1,kjpindex

             IF (ok_freeze_cwrr) THEN
                us_tmp = (trois*(1-profil_froz_hydro_ns(ji,nslm, jst))*mc(ji,nslm,jst)  &
                     + (1-profil_froz_hydro_ns(ji,nslm-1, jst))*mc(ji,nslm-1,jst))  &
                     / (quatre*mcs(njsc(ji))*pcent(njsc(ji)))
             ELSE 
                us_tmp = (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))  &
                     / (quatre*mcs(njsc(ji))*pcent(njsc(ji)))
             ENDIF

             ! us_tmp should not be negative
             us_tmp=MAX(us_tmp, zero)
             !
             us(ji,jv,jst,nslm) =moderwilt(ji,nslm,jst)*  &
                  & MIN( un, SQRT(us_tmp) ) * nroot(jv,jst,nslm)

             us(ji,jv,jst,nslm) = MAX(us(ji,jv,jst,nslm), zero)
             humrelv(ji,jv,jst) = MAX(zero,MIN(un, humrelv(ji,jv,jst) + us(ji,jv,jst,nslm)))
             vegstressv(ji,jv,jst) = humrelv(ji,jv,jst)
          END DO
       END DO


       DO jv = 2, nvm
          DO ji = 1, kjpindex
             IF (corr_veg_soil(ji,jv,jst) .LT. min_sechiba) THEN
                humrelv(ji,jv,jst) = zero  
             ENDIF
          END DO
       END DO


    !! 13 before closing the soil water, we check the water balance of soil

       IF(check_cwrr) THEN
          DO ji = 1,kjpindex
             
             deltahum(ji) = (tmc(ji,jst) - tmcold(ji))
             diff(ji)     = precisol_ns(ji,jst)-ru_ns(ji,jst)-dr_ns(ji,jst)-tsink(ji) &
                  & + irrigation_soil(ji) + returnflow_soil(ji) + reinfiltration_soil(ji)

             test(ji) = (ABS(deltahum(ji)-diff(ji))*mask_soiltile(ji,jst) .GT. allowed_err)
          ENDDO
             
          DO ji = 1,kjpindex
             IF(test(ji)) THEN
                
                WRITE (numout,*)'CWRR pat: bilan non nul',ji,jst,njsc(ji),deltahum(ji)-diff(ji)
                WRITE (numout,*)'tmc,tmcold,diff',tmc(ji,jst),tmcold(ji),deltahum(ji)
                WRITE(numout,*) 'evapot,evapot_penm,ae_ns',evapot(ji),evapot_penm(ji),ae_ns(ji,jst)
                WRITE (numout,*)'flux_top,ru_ns,qdrain,tsink,q0,precisol',flux_top(ji),ru_ns(ji,jst), &
                     &      dr_ns(ji,jst),tsink(ji),qflux00(ji,jst),precisol_ns(ji,jst)
                WRITE (numout,*)'water2infilt',water2infilt(ji,jst)
                WRITE (numout,*)'soiltile',soiltile(ji,jst)
                WRITE (numout,*)'irrigation, returnflow, reinfiltration', &
                     & irrigation_soil(ji),returnflow_soil(ji),reinfiltration_soil(ji)
                WRITE (numout,*)'mc',mc(ji,:,jst)
                WRITE (numout,*)'qflux',qflux(ji,:,jst)
                WRITE (numout,*)'veget_max', veget_max(ji,:)
                WRITE (numout,*)'k', k(ji,:)
                
                error=.TRUE.
                CALL ipslerr_p(2, 'hydrol_soil', 'We will STOP in the end of this subroutine.',&
                     & 'CWRR water balance check','')
             ENDIF
          ENDDO

          DO ji = 1,kjpindex
             
             IF(MINVAL(mc(ji,:,jst)).LT.-min_sechiba) THEN
                WRITE (numout,*)'CWRR MC NEGATIVE', &
                     & ji,lalo(ji,:),MINLOC(mc(ji,:,jst)),jst,mc(ji,:,jst)
                WRITE (numout,*)'evapot,evapot_penm,ae_ns',evapot(ji),evapot_penm(ji),ae_ns(ji,jst)
                WRITE (numout,*)'flux_top,ru_ns,qdrain,tsink,q0,precisol',flux_top(ji),ru_ns(ji,jst), &
                     &      dr_ns(ji,jst),tsink(ji),qflux00(ji,jst),precisol_ns(ji,jst)
                WRITE (numout,*)'water2infilt',water2infilt(ji,jst)
                WRITE (numout,*)'soiltile',soiltile(ji,jst)
                WRITE (numout,*)'irrigation, returnflow, reinfiltration', &
                     & irrigation_soil(ji),returnflow_soil(ji),reinfiltration_soil(ji)
                WRITE (numout,*)'mc',mc(ji,:,jst)
                WRITE (numout,*)'qflux',qflux(ji,:,jst)
                WRITE (numout,*)'veget_max', veget_max(ji,:)
                WRITE (numout,*)'k', k(ji,:)
                WRITE (numout,*)'soiltile',soiltile(ji,jst)

                error=.TRUE.
                CALL ipslerr_p(2, 'hydrol_soil', 'We will STOP in the end of this subroutine.',&
                     & 'CWRR MC NEGATIVE','')
             ENDIF
          END DO

          DO ji=1,kjpindex
             IF (ru_ns(ji,jst)*soiltile(ji,jst).LT.-min_sechiba) THEN
                WRITE (numout,*) 'Negative runoff', ji,jst, mask_soiltile(ji,jst) 
                WRITE (numout,*) 'mc1, mc2', mc(ji,1,jst), mc(ji,2,jst)
                WRITE (numout,*) 'mcint1, mcint2', mcint(ji,1), mcint(ji,2)
                WRITE (numout,*) 'qflux1, flux_top', qflux(ji,nslm,jst), flux_top(ji)
                WRITE (numout,*) 'is_over_mcs, is_under_mcr, test', &
                     & is_over_mcs(ji), is_under_mcr(ji), tmc(ji,jst)-tmcint(ji)+qflux(ji,nslm,jst)+SUM(rootsink(ji,:,jst))
                WRITE (numout,*)'mc', mc(ji,:,jst)
                WRITE (numout,*)'mcint', mcint(ji,:)
                WRITE (numout,*)'qflux', qflux(ji,:,jst)
                WRITE (numout,*)'rootsink1,evapot_penm,vegtot', rootsink(ji,1,jst), evapot_penm(ji), vegtot(ji)
                WRITE (numout,*) 'ae_ns, tsink, returnflow, reinfiltration, precisol_ns, irrigation, qflux0, ru_ns', &
                     & ae_ns(ji,jst), tsink(ji), returnflow_soil(ji), reinfiltration_soil(ji), &
                     & precisol_ns(ji,jst), irrigation_soil(ji), qflux00(ji,jst), ru_ns(ji,jst)

                error=.TRUE.
                CALL ipslerr_p(2, 'hydrol_soil', 'We will STOP in the end of this subroutine.',&
                     & 'Negative runoff, non-saturated soil','')
             ENDIF
          ENDDO
       ENDIF

       IF (long_print) WRITE (numout,*) ' hydrol_soil done for jst =', jst     

       IF (ok_freeze_cwrr) THEN
          DO ji = 1, kjpindex
             kk_moy(ji,:) =kk_moy(ji,:)+soiltile(ji,jst)*k(ji,:) 
             kk(ji,:,jst)=k(ji,:)
          ENDDO
       ENDIF

    END DO  ! end of loop on soiltile


    !
    !! 14 sum 3d variables into 2d variables
    !
    CALL hydrol_diag_soil (kjpindex, veget_max, soiltile, njsc, runoff, drainage, &
         & evapot, vevapnu, returnflow, reinfiltration, irrigation, &
         & shumdiag,shumdiag_perma, k_litt, litterhumdiag, humrel, vegstress, drysoil_frac,tot_melt) !! Arsene 28-01-2016 - REMOVE drunoff_tot because never user and bug in sechiba_output.f90
 !pss:+       !pss:-

    !
    !! 15 Calculation of evap_bare_limit : limitation factor for bare soil evaporation
    !
    evap_bare_lim(:) = zero
    evap_bare_lim_ns(:,:) = zero

    !!_  for each soil tile
    !   
    DO jst = 1,nstm

   !! 15.1 Save mc and tmc for use at the end of the time step
   !!      and calculate tmcint corresponding to mc without water2infilt
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mcint(ji,jsl) = mask_soiltile(ji,jst) * mc(ji,jsl,jst)
          ENDDO
       ENDDO

       DO ji = 1, kjpindex
          temp(ji) = tmc(ji,jst)
       ENDDO

       ! Calculate tmc corresponding to mc without water2infilt and save it for later use
       DO ji = 1, kjpindex
          tmcint(ji) = dz(2,jst) * ( trois*mcint(ji,1) + mcint(ji,2) )/huit
       ENDDO

       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmcint(ji) = tmcint(ji) + dz(jsl,jst) &
                  * (trois*mcint(ji,jsl)+mcint(ji,jsl-1))/huit &
                  + dz(jsl+1,jst) * (trois*mcint(ji,jsl)+mcint(ji,jsl+1))/huit
          ENDDO
       ENDDO

       DO ji = 1, kjpindex
          tmcint(ji) = tmcint(ji) + dz(nslm,jst) &
               * (trois * mcint(ji,nslm) + mcint(ji,nslm-1))/huit
       ENDDO



     !! 15.2  Coefficient are recomputed for the profile of mc
       CALL hydrol_soil_coef(kjpindex,jst,njsc)

     !! 15.3 Set the values for diffusion scheme
       CALL hydrol_soil_setup(kjpindex,jst,dtradia)
       resolv(:) = (mask_soiltile(:,jst) .GT. 0)

     !! 15.3 compute matrix coefficients 
     !- estimate maximum evaporation flux_top in mm/step, assuming the water is available
       DO ji = 1, kjpindex
          IF(vegtot(ji).GT.min_sechiba) THEN
             flux_top(ji) = evapot_penm(ji) * &
                  AINT(frac_bare_ns(ji,jst)+un-min_sechiba)
          ELSE
             flux_top(ji) = zero
          ENDIF
       ENDDO

       IF (ok_freeze_cwrr) THEN
          DO ji =1, kjpindex
             DO jsl = 1, nslm
                mcl(ji,jsl,:)= MIN(mc(ji,jsl,:),mcr(njsc(ji))+(1-profil_froz_hydro_ns(ji,jsl, :))*(mc(ji,jsl,:)-mcr(njsc(ji))))
             ENDDO
          ENDDO
       ELSE
          mcl(:,:,:)=mc(:,:,:)
       ENDIF

       !- First layer
       DO ji = 1, kjpindex
          tmat(ji,1,1) = zero
          tmat(ji,1,2) = f(ji,1)
          tmat(ji,1,3) = g1(ji,1)
          rhs(ji,1)    = fp(ji,1) * mcl(ji,1,jst) + gp(ji,1)*mcl(ji,2,jst) &
               - flux_top(ji) - (b(ji,1)+b(ji,2))/deux *(dtradia/one_day)
       ENDDO

       !- soil body
       DO jsl=2, nslm-1
          DO ji = 1, kjpindex
             tmat(ji,jsl,1) = e(ji,jsl)
             tmat(ji,jsl,2) = f(ji,jsl)
             tmat(ji,jsl,3) = g1(ji,jsl)
             rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
                  +  gp(ji,jsl) * mcl(ji,jsl+1,jst) &
                  + (b(ji,jsl-1) - b(ji,jsl+1)) * (dtradia/one_day) / deux
          ENDDO
       ENDDO

       !- Last layer
       DO ji = 1, kjpindex
          jsl=nslm
          tmat(ji,jsl,1) = e(ji,jsl)
          tmat(ji,jsl,2) = f(ji,jsl)
          tmat(ji,jsl,3) = zero
          rhs(ji,jsl) = ep(ji,jsl)*mcl(ji,jsl-1,jst) + fp(ji,jsl)*mcl(ji,jsl,jst) &
               + (b(ji,jsl-1) + b(ji,jsl)*(un-deux*free_drain_coef(ji,jst))) * (dtradia/one_day) / deux
       ENDDO

       !- Store the equations for later use
       DO jsl=1,nslm
          DO ji = 1, kjpindex
             srhs(ji,jsl) = rhs(ji,jsl)
             stmat(ji,jsl,1) = tmat(ji,jsl,1)
             stmat(ji,jsl,2) = tmat(ji,jsl,2)
             stmat(ji,jsl,3) = tmat(ji,jsl,3)
          ENDDO
       ENDDO

    !! 15.5.1 Resolution with flux=evapot_penm (update mcl)
       CALL hydrol_soil_tridiag(kjpindex,jst)

    !! 15.5.2 Resolution with mc(1)=mcr
       !! Prepare to rerun in case of under residual with evaporation 
       DO ji = 1, kjpindex
          resolv(ji) = (mcl(ji,1,jst).LT.(mcr(njsc(ji))).AND.flux_top(ji).GT.min_sechiba)
       ENDDO
       !! Reset the coefficient for diffusion (only used if resolv(ji) = .TRUE.)
       DO jsl=1,nslm
          !- The new condition is to put the upper layer at residual soil moisture
          DO ji = 1, kjpindex
             rhs(ji,jsl) = srhs(ji,jsl)
             tmat(ji,jsl,1) = stmat(ji,jsl,1)
             tmat(ji,jsl,2) = stmat(ji,jsl,2)
             tmat(ji,jsl,3) = stmat(ji,jsl,3)
          END DO
       END DO
       
       DO ji = 1, kjpindex
          tmat(ji,1,2) = un
          tmat(ji,1,3) = zero
          rhs(ji,1) = mcr(njsc(ji))
       ENDDO
       
       !! Resolve the equations with new boundary conditions if necessary (update mcl)
       CALL hydrol_soil_tridiag(kjpindex,jst)

       ! Calculation of total soil moisture content (liquid + frozen)
       IF (ok_freeze_cwrr) THEN
          DO ji =1, kjpindex
             DO jsl = 1, nslm
                mc(ji,jsl,:)=MAX(mcl(ji,jsl,:), mcl(ji,jsl,:)+profil_froz_hydro_ns(ji,jsl,:)*(mc(ji,jsl,:)-mcr(njsc(ji))))
             ENDDO
          ENDDO
       ELSE
          mc(:,:,:)=mcl(:,:,:)
       ENDIF


       !! Correct bad moisture content due to numerical errors before water balance
       DO ji = 1, kjpindex
          DO jsl = 1,nslm
             IF (mc(ji,jsl,jst) .LT. mcr(njsc(ji)) .AND. mc(ji,jsl,jst) .NE. 0.) THEN
                mc(ji,jsl,jst)=mcr(njsc(ji))
             END IF
             IF (mc(ji,jsl,jst) .GT. mcs(njsc(ji)) .AND. mc(ji,jsl,jst) .NE. 0.) THEN
                mc(ji,jsl,jst)=mcs(njsc(ji))
             END IF
          ENDDO
       ENDDO

    !! 15.6 Water balance

    !! 15.6.1 Compute dr_ns with the bottom boundary condition 

    !         Initialize qflux at bottom of diffusion
       DO ji = 1, kjpindex
          flux_bottom(ji) = mask_soiltile(ji,jst)*k(ji,nslm)*free_drain_coef(ji,jst) * (dtradia/one_day)
       ENDDO

    !! 15.6.2 Compute total soil moisture content
       DO ji = 1, kjpindex
          tmc(ji,jst) = dz(2,jst) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
       ENDDO
       
       DO jsl = 2,nslm-1
          DO ji = 1, kjpindex
             tmc(ji,jst) = tmc(ji,jst) + dz(jsl,jst) &
                  * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                  + dz(jsl+1,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit
          ENDDO
       ENDDO
       
       DO ji = 1, kjpindex
          tmc(ji,jst) = tmc(ji,jst) + dz(nslm,jst) &
               * (trois * mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO
    
    !! 15.6.3 Deducte upper flux from soil moisture variation and bottom flux
    !! TMCi-D-E=TMC   
       DO ji = 1, kjpindex
          evap_bare_lim_ns(ji,jst) = mask_soiltile(ji,jst) * &
               (tmcint(ji)-tmc(ji,jst)-flux_bottom(ji))
       END DO

    !! 15.7 Determination of evap_bar_lim_ns
       DO ji = 1, kjpindex
          ! Here we weight evap_bare_lim_ns by the fraction of bare evaporating soil. 
          ! This is given by frac_bare_ns, taking into account bare soil under vegetation
          IF(vegtot(ji) .GT. min_sechiba) THEN
             evap_bare_lim_ns(ji,jst) = evap_bare_lim_ns(ji,jst) * frac_bare_ns(ji,jst)
          ELSE
             evap_bare_lim_ns(ji,jst) = 0.
          ENDIF
       END DO

       DO ji=1,kjpindex
          IF ((evapot(ji).GT.min_sechiba) .AND. &
               (tmc_litter(ji,jst).GT.(tmc_litter_wilt(ji,jst)))) THEN
             evap_bare_lim_ns(ji,jst) = evap_bare_lim_ns(ji,jst) / evapot(ji)
          ELSEIF((evapot(ji).GT.min_sechiba).AND. &
               (tmc_litter(ji,jst).GT.(tmc_litter_res(ji,jst)))) THEN
             evap_bare_lim_ns(ji,jst) =  (un/deux) * evap_bare_lim_ns(ji,jst) / evapot(ji)
          END IF
          evap_bare_lim_ns(ji,jst)=MAX(MIN(evap_bare_lim_ns(ji,jst),1.),0.)
       END DO

    !! 15.8 Restore mc and tmc
       DO jsl = 1, nslm
          DO ji = 1, kjpindex
             mc(ji,jsl,jst) = mask_soiltile(ji,jst) * mcint(ji,jsl)
          ENDDO
       ENDDO

       DO ji = 1, kjpindex
          tmc(ji,jst) = temp(ji)
       ENDDO

    ENDDO !end loop on tiles

    !! 15.9 Deduction of evap_bar_lim
    DO ji = 1, kjpindex
       evap_bare_lim(ji) =  SUM(evap_bare_lim_ns(ji,:)*vegtot(ji)*soiltile(ji,:))
    ENDDO


    !
    !! 16 Exit if error was found previously in this subroutine
    !
    IF ( error ) THEN
       WRITE(numout,*) 'One or more errors have been detected in hydrol_soil. Model stops.'
       CALL ipslerr_p(3, 'hydrol_soil', 'We will STOP now.',&
                  & 'One or several fatal errors were found previously.','')
    END IF

  END SUBROUTINE hydrol_soil


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_infilt
!!
!>\BRIEF        Infiltration
!!
!! DESCRIPTION  :
!! - 1 First layer
!! - 2 Infiltration layer by layer 
!! - 2.1 Initialisation
!! - 2.2 Infiltrability of each layer if under a saturated one
!! - 2.3 Compute the mean rate at which water actually infiltrate:
!! - 2.4 From which we deduce the time it takes to fill up the layer or to end the time step...
!! - 2.5 The water enters in the layer
!! - 3 Verification
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_infilt

  SUBROUTINE hydrol_soil_infilt(kjpindex, ins, dtradia, njsc, flux_infilt)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! GLOBAL (in or inout)
    INTEGER(i_std), INTENT(in)                        :: kjpindex        !! Domain size
    REAL(r_std), INTENT (in)                          :: dtradia         !! Time step in seconds
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)  :: njsc            !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)    :: flux_infilt     !! Water to infiltrate

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                :: ji, jsl, ins        !! Indices
    REAL(r_std), DIMENSION (kjpindex)             :: wat_inf_pot         !! infiltrable water in the layer
    REAL(r_std), DIMENSION (kjpindex)             :: wat_inf             !! infiltrated water in the layer
    REAL(r_std), DIMENSION (kjpindex)             :: dt_tmp              !! time remaining before the end of the time step
    REAL(r_std), DIMENSION (kjpindex)             :: dt_inf              !! the time it takes to complete the infiltration in the 
                                                                         !! layer 
    REAL(r_std)                                   :: k_m                 !! the mean conductivity used for the saturated front
    REAL(r_std), DIMENSION (kjpindex)             :: infilt_tmp          !! infiltration rate for the considered layer 
    REAL(r_std), DIMENSION (kjpindex)             :: infilt_tot          !! total infiltration 
    REAL(r_std), DIMENSION (kjpindex)             :: flux_tmp            !! rate at which precip hits the ground

!_ ================================================================================================================================

    ! If data (or coupling with GCM) was available, a parameterization for subgrid rainfall could be performed

    DO ji = 1, kjpindex
       !-
    !_ 1 First layer
       !-
       ! First we fill up the first layer (about 1mm) without any resistance and quasi-immediately
       wat_inf_pot(ji) = MAX((mcs(njsc(ji))-mc(ji,1,ins)) * dz(2,ins) / deux, zero)
       wat_inf(ji) = MIN(wat_inf_pot(ji), flux_infilt(ji))
       mc(ji,1,ins) = mc(ji,1,ins) + wat_inf(ji) * deux / dz(2,ins)
       !
    ENDDO


    !-
    !! 2 Infiltration layer by layer 
    !! 2.1 Initialisation
    ! Initialize a countdown for infiltration during the time-step and the value of potential runoff
    dt_tmp(:) = dtradia / one_day
    infilt_tot(:) = wat_inf(:)
    ! Compute the rate at which water will try to infiltrate each layer

    flux_tmp(:) = (flux_infilt(:)-wat_inf(:)) / dt_tmp(:)
    !
    DO jsl = 2, nslm-1
       DO ji = 1, kjpindex
    !! 2.2 Infiltrability of each layer if under a saturated one
          ! This is computed by an simple arithmetic average because 
          ! the time step (30min) is not appropriate for a geometric average (advised by Haverkamp and Vauclin)
          k_m = (k(ji,jsl) + ks(njsc(ji))*kfact(jsl-1,njsc(ji))*kfact_root(ji,jsl,ins)) / deux 

          IF (ok_freeze_cwrr) THEN
             IF (temp_hydro(ji, jsl) .LT. ZeroCelsius) THEN
                k_m = k(ji,jsl)
             ENDIF
          ENDIF

    !! 2.3 We compute the mean rate at which water actually infiltrate:
          !- Subgrid: Exponential distribution of  k around i_m, but average p directly used 
          !!??Aurelien: A big black cloud for me
          infilt_tmp(ji) = k_m * (un - EXP(- flux_tmp(ji) / k_m)) 

    !! 2.4 From which we deduce the time it takes to fill up the layer or to end the time step...
          wat_inf_pot(ji) =  MAX((mcs(njsc(ji))-mc(ji,jsl,ins)) * (dz(jsl,ins) + dz(jsl+1,ins)) / deux, zero)
          IF ( infilt_tmp(ji) > min_sechiba) THEN
             dt_inf(ji) =  MIN(wat_inf_pot(ji)/infilt_tmp(ji), dt_tmp(ji))
             ! The water infiltration TIME has to limited by what is still available for infiltration.
             IF ( dt_inf(ji) * infilt_tmp(ji) > flux_infilt(ji)-infilt_tot(ji) ) THEN
                dt_inf(ji) = MAX(flux_infilt(ji)-infilt_tot(ji),zero)/infilt_tmp(ji)
             ENDIF
          ELSE
             dt_inf(ji) = dt_tmp(ji)
          ENDIF

    !! 2.5 The water enters in the layer
          wat_inf(ji) = dt_inf(ji) * infilt_tmp(ji)
          ! bviously the moisture content
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + &
               & wat_inf(ji) * deux / (dz(jsl,ins) + dz(jsl+1,ins))
          ! the time remaining before the next time step
          dt_tmp(ji) = dt_tmp(ji) - dt_inf(ji)
          ! and finally the infilt_tot (which is just used to check if there is a problem, below) 
          infilt_tot(ji) = infilt_tot(ji) + infilt_tmp(ji) * dt_inf(ji)
       ENDDO

    ENDDO


    !! 3 Verification
    DO ji = 1, kjpindex
       IF (infilt_tot(ji) .LT. -min_sechiba .OR. infilt_tot(ji) .GT. flux_infilt(ji) + min_sechiba) THEN
          WRITE (numout,*) 'Error in the calculation of infilt tot', infilt_tot(ji)
          WRITE (numout,*) 'k, ji, jst, mc', k(ji,1:2), ji, ins, mc(ji,1,ins)
          CALL ipslerr_p(3, 'hydrol_soil_infilt', 'We will STOP now.','Error in calculation of infilt tot','')
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE hydrol_soil_infilt


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_smooth
!!
!>\BRIEF        Smooth soil moisture values: avoid over-saturation or under-residual values.
!!
!! DESCRIPTION  :
!! - 1 Avoid over-saturation values
!! - 1.1 top to bottom
!! - 1.2 bottom to top
!! - 2 Avoid below residual values
!! - 2.1 top to bottom
!! - 2.2 bottom to top
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_smooth

  SUBROUTINE hydrol_soil_smooth(kjpindex, ins, njsc, is_under_mcr, is_over_mcs)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! number of soil type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: njsc            !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)

    !! 0.2 Output variables

    LOGICAL, DIMENSION(kjpindex), INTENT(out)          :: is_under_mcr    !! Allows under residual soil moisture due to evap 
    LOGICAL, DIMENSION(kjpindex), INTENT(out)          :: is_over_mcs     !! Allows over saturated soil moisture due to returnflow 

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                       :: ji,jsl
    REAL(r_std)                          :: excess
    REAL(r_std), DIMENSION(kjpindex)     :: excessji

!_ ================================================================================================================================       
    !-
    !! 1 Avoid over-saturation values
    !-
    
    ! in case of over-saturation we put the water where it is possible

    !! 1.1 top to bottom
    DO jsl = 1, nslm-2
       DO ji=1, kjpindex
          excess = MAX(mc(ji,jsl,ins)-mcs(njsc(ji)),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
          mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) + excess * &
               &  (dz(jsl,ins)+dz(jsl+1,ins))/(dz(jsl+1,ins)+dz(jsl+2,ins))
       ENDDO
    ENDDO

    jsl = nslm-1
    DO ji=1, kjpindex
       excess = MAX(mc(ji,jsl,ins)-mcs(njsc(ji)),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
       mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) + excess * &
            &  (dz(jsl,ins)+dz(jsl+1,ins))/dz(jsl+1,ins)
    ENDDO

    jsl = nslm
    DO ji=1, kjpindex
       excess = MAX(mc(ji,jsl,ins)-mcs(njsc(ji)),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
       mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) + excess * &
            &  dz(jsl,ins)/(dz(jsl-1,ins)+dz(jsl,ins))
    ENDDO
    !! 1.2 bottom to top
    DO jsl = nslm-1,2,-1
       DO ji=1, kjpindex
          excess = MAX(mc(ji,jsl,ins)-mcs(njsc(ji)),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excess
          mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) + excess * &
               &  (dz(jsl,ins)+dz(jsl+1,ins))/(dz(jsl-1,ins)+dz(jsl,ins))
       ENDDO
    ENDDO

    DO ji=1, kjpindex
       excessji(ji) = mask_soiltile(ji,ins) * MAX(mc(ji,1,ins)-mcs(njsc(ji)),zero)
    ENDDO

    DO ji=1, kjpindex
       mc(ji,1,ins) = mc(ji,1,ins) - excessji(ji)

       is_over_mcs(ji) = (excessji(ji) .GT. min_sechiba)
    ENDDO
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excessji(ji) * dz(2,ins) / (deux * dpu_max*mille)
       ENDDO
    ENDDO

    !-
    !! 2 Avoid below residual values
    !-
       
    ! Smooth the profile to avoid negative values of punctual soil moisture

    !! 2.1 top to bottom
    DO jsl = 1,nslm-2
       DO ji=1, kjpindex
          excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
          mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) - excess * &
               &  (dz(jsl,ins)+dz(jsl+1,ins))/(dz(jsl+1,ins)+dz(jsl+2,ins))
       ENDDO
    ENDDO

    jsl = nslm-1
    DO ji=1, kjpindex
       excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
       mc(ji,jsl+1,ins) = mc(ji,jsl+1,ins) - excess * &
            &  (dz(jsl,ins)+dz(jsl+1,ins))/dz(jsl+1,ins)
    ENDDO

    jsl = nslm
    DO ji=1, kjpindex
       excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
       mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
       mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) - excess * &
            &  dz(jsl,ins)/(dz(jsl-1,ins)+dz(jsl,ins))
    ENDDO

    !! 2.2 bottom to top
    DO jsl = nslm-1,2,-1
       DO ji=1, kjpindex
          excess = MAX(mcr(njsc(ji))-mc(ji,jsl,ins),zero)
          mc(ji,jsl,ins) = mc(ji,jsl,ins) + excess
          mc(ji,jsl-1,ins) = mc(ji,jsl-1,ins) - excess * &
               &  (dz(jsl,ins)+dz(jsl+1,ins))/(dz(jsl-1,ins)+dz(jsl,ins))
       ENDDO
    ENDDO

    DO ji=1, kjpindex
       excessji(ji) = mask_soiltile(ji,ins) * MAX(mcr(njsc(ji))-mc(ji,1,ins),zero)
    ENDDO
    DO ji=1, kjpindex
       mc(ji,1,ins) = mc(ji,1,ins) + excessji(ji)
       is_under_mcr(ji) = (excessji(ji) .GT. min_sechiba)
    ENDDO

    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mc(ji,jsl,ins) - excessji(ji) * dz(2,ins) / (deux * dpu_max*mille)
       ENDDO
    ENDDO

    ! We just get sure that mc remains at 0 where soiltile=0
    DO jsl = 1, nslm
       DO ji=1, kjpindex
          mc(ji,jsl,ins) = mask_soiltile(ji,ins) * mc(ji,jsl,ins)
       ENDDO
    ENDDO
    
    RETURN
  END SUBROUTINE hydrol_soil_smooth


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_flux
!!
!>\BRIEF        This subroutine computes the hydrological fluxes between the different soil layers.
!!
!! DESCRIPTION  :
!! - 1 Initialize qflux from the bottom, with dr_ns
!! - 2 between layer nslm and nslm-1   
!! - 3 we go to top and deduct qflux(1:nslm-2)   
!! - 4 Water balance verification
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_flux

  SUBROUTINE hydrol_soil_flux(kjpindex,ins,mcint,returnflow_soil)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! index of soil type
    REAL(r_std), DIMENSION (kjpindex,nslm), INTENT(in) :: mcint           !! mc values at the beginning of the time step
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)      :: returnflow_soil !! returnflow

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: jsl,ji
    REAL(r_std), DIMENSION(kjpindex)                   :: temp

!_ ================================================================================================================================    
    !- Compute the flux at every level from bottom to top (using mc and sink values)
    DO ji = 1, kjpindex
    !! 1 Initialize qflux from the bottom, with dr_ns
       jsl = nslm
       qflux(ji,jsl,ins) = dr_ns(ji,ins) - returnflow_soil(ji)
    !!_ between layer nslm and nslm-1   
       jsl = nslm-1
       qflux(ji,jsl,ins) = qflux(ji,jsl+1,ins) & 
            &  + (mc(ji,jsl,ins)-mcint(ji,jsl) &
            &  + trois*mc(ji,jsl+1,ins) - trois*mcint(ji,jsl+1)) &
            &  * (dz(jsl+1,ins)/huit) &
            &  + rootsink(ji,jsl+1,ins) 
    ENDDO
    !! 3 we go to top and deduct qflux(1:nslm-2)   
    DO jsl = nslm-2,1,-1
       DO ji = 1, kjpindex
          qflux(ji,jsl,ins) = qflux(ji,jsl+1,ins) & 
               &  + (mc(ji,jsl,ins)-mcint(ji,jsl) &
               &  + trois*mc(ji,jsl+1,ins) - trois*mcint(ji,jsl+1)) &
               &  * (dz(jsl+1,ins)/huit) &
               &  + rootsink(ji,jsl+1,ins) &
               &  + (dz(jsl+2,ins)/huit) &
               &  * (trois*mc(ji,jsl+1,ins) - trois*mcint(ji,jsl+1) &
               &  + mc(ji,jsl+2,ins)-mcint(ji,jsl+2)) 
       END DO
    ENDDO
    
    !! 4 Water balance verification  
    DO ji = 1, kjpindex
       temp(ji) =  qflux(ji,1,ins) + (dz(2,ins)/huit) &
            &  * (trois* (mc(ji,1,ins)-mcint(ji,1)) + (mc(ji,2,ins)-mcint(ji,2))) &
            &  + rootsink(ji,1,ins)
    ENDDO

    DO ji = 1, kjpindex
       IF (ABS(qflux00(ji,ins)-temp(ji)).GT. deux*min_sechiba) THEN
          WRITE(numout,*) 'Problem in the water balance, qflux computation', qflux00(ji,ins),temp(ji)
          WRITE (numout,*) 'returnflow_soil', returnflow_soil(ji)
          WRITE(numout,*) 'ji', ji, 'jsl',jsl,'ins',ins
          WRITE(numout,*) 'mcint', mcint(ji,:)
          WRITE(numout,*) 'mc', mc(ji,:,ins)
          WRITE (numout,*) 'rootsink', rootsink(ji,1,ins)
          CALL ipslerr_p(3, 'hydrol_soil_flux', 'We will STOP now.',&
               & 'Problem in the water balance, qflux computation','')
       ENDIF
    ENDDO

    RETURN
  END SUBROUTINE hydrol_soil_flux


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_tridiag
!!
!>\BRIEF        This subroutine solves a set of linear equations which has a tridiagonal coefficient matrix. 
!!
!! DESCRIPTION  : None
!! 
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : mcl (global module variable)
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_tridiag 

  SUBROUTINE hydrol_soil_tridiag(kjpindex,ins)

    !- arguments

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex        !! Domain size
    INTEGER(i_std), INTENT(in)                         :: ins             !! number of soil type

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ji,jsl
    REAL(r_std), DIMENSION(kjpindex)                   :: bet
    REAL(r_std), DIMENSION(kjpindex,nslm)              :: gam

!_ ================================================================================================================================
    DO ji = 1, kjpindex

       IF (resolv(ji)) THEN
          bet(ji) = tmat(ji,1,2)
          mcl(ji,1,ins) = rhs(ji,1)/bet(ji)
       ENDIF
    ENDDO

    DO jsl = 2,nslm
       DO ji = 1, kjpindex
          
          IF (resolv(ji)) THEN

             gam(ji,jsl) = tmat(ji,jsl-1,3)/bet(ji)
             bet(ji) = tmat(ji,jsl,2) - tmat(ji,jsl,1)*gam(ji,jsl)
             mcl(ji,jsl,ins) = (rhs(ji,jsl)-tmat(ji,jsl,1)*mcl(ji,jsl-1,ins))/bet(ji)
          ENDIF

       ENDDO
    ENDDO

    DO ji = 1, kjpindex
       IF (resolv(ji)) THEN
          DO jsl = nslm-1,1,-1
             mcl(ji,jsl,ins) = mcl(ji,jsl,ins) - gam(ji,jsl+1)*mcl(ji,jsl+1,ins)
          ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE hydrol_soil_tridiag


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_coef
!!
!>\BRIEF        Computes coef for the linearised hydraulic conductivity 
!! k_lin=a_lin mc_lin+b_lin and the linearised diffusivity d_lin. 
!!
!! DESCRIPTION  :
!! First, we identify the interval i in which the current value of mc is located.
!! Then, we give the values of the linearized parameters to compute 
!! conductivity and diffusivity as K=a*mc+b and d.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_soil_coef
 
  SUBROUTINE hydrol_soil_coef(kjpindex,ins,njsc)

    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                        :: kjpindex         !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins              !! Index of soil type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)  :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                    :: jsl,ji,i
    REAL(r_std)                                       :: mc_ratio
    REAL(r_std)                                       :: mc_used    !! Used liquid water content
    REAL(r_std)                                       :: x,m
    
!_ ================================================================================================================================

    IF (ok_freeze_cwrr) THEN
    
       ! Calculation of liquid and frozen saturation degrees with respect to residual
       ! x=liquid saturation degree/residual=(mcl-mcr)/(mcs-mcr)
       ! 1-x=frozen saturation degree/residual=(mcf-mcr)/(mcs-mcr) (=profil_froz_hydro)
       
       ! Van Genuchten parameter for thermodynamical calculation
       m = 1.-1./nvan(ins)
       
       DO jsl=1,nslm
          DO ji=1,kjpindex    
             IF ((.NOT. ok_thermodynamical_freezing).OR.(mc(ji,jsl, ins).LT.(mcr(njsc(ji))+min_sechiba))) THEN
                ! Linear soil freezing or soil moisture below residual
                IF (temp_hydro(ji, jsl).GE.273._r_std) THEN
                   x=1._r_std
                ELSE IF (273._r_std.GT.temp_hydro(ji, jsl).AND.temp_hydro(ji, jsl).GE.271._r_std) THEN 
                   x=(temp_hydro(ji, jsl)-271._r_std)/fr_dT
                ELSE 
                   x=0._r_std
                ENDIF
             ELSE IF (ok_thermodynamical_freezing) THEN
                ! Thermodynamical soil freezing
                IF (temp_hydro(ji, jsl).GE.273._r_std) THEN
                   x=1._r_std
                ELSE IF (273._r_std.GT.temp_hydro(ji, jsl).AND.temp_hydro(ji, jsl).GE.271._r_std) THEN 
                   x=MIN(((mcs(njsc(ji))-mcr(njsc(ji))) &
                        *((1000.*avan(njsc(ji))*(ZeroCelsius-temp_hydro(ji, jsl)) &
                        *lhf/ZeroCelsius/10.)**nvan(njsc(ji))+1.)**(-m))/(mc(ji,jsl, ins) &
                        -mcr(njsc(ji))),1._r_std)                
                ELSE
                   x=0._r_std
                ENDIF
             ENDIF
             profil_froz_hydro_ns(ji, jsl,ins)=1._r_std-x
             
             ! mc_used is used in the calculation of hydrological properties
             mc_used = mcr(njsc(ji))+x*(mc(ji,jsl, ins)-mcr(njsc(ji))) 
             !
             ! calcul de k based on mc_liq
             !
             i= MAX(imin, MIN(imax-1, INT(imin +(imax-imin)*(mc_used-mcr(njsc(ji)))/(mcs(njsc(ji))-mcr(njsc(ji))))))
             a(ji,jsl) = a_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
             b(ji,jsl) = b_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
             d(ji,jsl) = d_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
             k(ji,jsl) = MAX(k_lin(imin+1,jsl,njsc(ji)), &
                  a_lin(i,jsl,njsc(ji)) * mc_used + b_lin(i,jsl,njsc(ji)))
          ENDDO ! loop on grid
       ENDDO
             
    ELSE
       ! .NOT. ok_freeze_cwrr
       DO jsl=1,nslm
          DO ji=1,kjpindex 
             
             mc_ratio = MAX(mc(ji,jsl,ins)-mcr(njsc(ji)), zero)/(mcs(njsc(ji))-mcr(njsc(ji)))
             
             i= MAX(MIN(INT((imax-imin)*mc_ratio)+imin , imax-1), imin)
             a(ji,jsl) = a_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
             b(ji,jsl) = b_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
             d(ji,jsl) = d_lin(i,jsl,njsc(ji)) * kfact_root(ji,jsl,ins)
             k(ji,jsl) = MAX(k_lin(imin+1,jsl,njsc(ji)), &
                  a_lin(i,jsl,njsc(ji)) * mc(ji,jsl,ins) + b_lin(i,jsl,njsc(ji))) 
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE hydrol_soil_coef


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_soil_setup
!!
!>\BRIEF        This subroutine computes the matrix coef.  
!!
!! DESCRIPTION  : None 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : matrix coef
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_soil_setup(kjpindex,ins,dtradia)


    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    REAL(r_std), INTENT (in)                           :: dtradia          !! Time step in seconds
    ! parameters
    INTEGER(i_std), INTENT(in)                        :: kjpindex          !! Domain size
    INTEGER(i_std), INTENT(in)                        :: ins               !! index of soil type

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: jsl,ji
    REAL(r_std)                        :: temp3, temp4

!_ ================================================================================================================================
    !-we compute tridiag matrix coefficients (LEFT and RIGHT) 
    ! of the system to solve [LEFT]*mc_{t+1}=[RIGHT]*mc{t}+[add terms]: 
    ! e(nslm),f(nslm),g1(nslm) for the [left] vector
    ! and ep(nslm),fp(nslm),gp(nslm) for the [right] vector

    ! w_time=1 (in constantes_soil) indicates implicit computation for diffusion 
    temp3 = w_time*(dtradia/one_day)/deux
    temp4 = (un-w_time)*(dtradia/one_day)/deux

    ! Passage to arithmetic means for layer averages also in this subroutine : Aurelien 11/05/10

    !- coefficient for first layer
    DO ji = 1, kjpindex
       e(ji,1) = zero
       f(ji,1) = trois * dz(2,ins)/huit  + temp3 &
            & * ((d(ji,1)+d(ji,2))/(dz(2,ins))+a(ji,1))
       g1(ji,1) = dz(2,ins)/(huit)       - temp3 &
            & * ((d(ji,1)+d(ji,2))/(dz(2,ins))-a(ji,2))
       ep(ji,1) = zero
       fp(ji,1) = trois * dz(2,ins)/huit - temp4 &
            & * ((d(ji,1)+d(ji,2))/(dz(2,ins))+a(ji,1))
       gp(ji,1) = dz(2,ins)/(huit)       + temp4 &
            & * ((d(ji,1)+d(ji,2))/(dz(2,ins))-a(ji,2))
    ENDDO

    !- coefficient for medium layers

    DO jsl = 2, nslm-1
       DO ji = 1, kjpindex
          e(ji,jsl) = dz(jsl,ins)/(huit)                        - temp3 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl,ins))+a(ji,jsl-1))

          f(ji,jsl) = trois * (dz(jsl,ins)+dz(jsl+1,ins))/huit  + temp3 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl,ins)) + &
               & (d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1,ins)) )

          g1(ji,jsl) = dz(jsl+1,ins)/(huit)                     - temp3 &
               & * ((d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1,ins))-a(ji,jsl+1))

          ep(ji,jsl) = dz(jsl,ins)/(huit)                       + temp4 &
               & * ((d(ji,jsl)+d(ji,jsl-1))/(dz(jsl,ins))+a(ji,jsl-1))

          fp(ji,jsl) = trois * (dz(jsl,ins)+dz(jsl+1,ins))/huit - temp4 &
               & * ( (d(ji,jsl)+d(ji,jsl-1))/(dz(jsl,ins)) + &
               & (d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1,ins)) )

          gp(ji,jsl) = dz(jsl+1,ins)/(huit)                     + temp4 &
               & *((d(ji,jsl)+d(ji,jsl+1))/(dz(jsl+1,ins))-a(ji,jsl+1))
       ENDDO
    ENDDO

    !- coefficient for last layer
    DO ji = 1, kjpindex
       e(ji,nslm) = dz(nslm,ins)/(huit)        - temp3 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm,ins))+a(ji,nslm-1))
       f(ji,nslm) = trois * dz(nslm,ins)/huit  + temp3 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) / (dz(nslm,ins)) &
            & -a(ji,nslm)*(un-deux*free_drain_coef(ji,ins)))
       g1(ji,nslm) = zero
       ep(ji,nslm) = dz(nslm,ins)/(huit)       + temp4 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm,ins))+a(ji,nslm-1))
       fp(ji,nslm) = trois * dz(nslm,ins)/huit - temp4 &
            & * ((d(ji,nslm)+d(ji,nslm-1)) /(dz(nslm,ins)) &
            & -a(ji,nslm)*(un-deux*free_drain_coef(ji,ins)))
       gp(ji,nslm) = zero
    ENDDO

    RETURN
  END SUBROUTINE hydrol_soil_setup


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_split_soil
!!
!>\BRIEF        Splits 2d variables into 3d variables, per soil type. 
!!
!! DESCRIPTION  :
!! - 1 split 2d variables into 3d variables, per soil type
!! - 1.1 precipitation
!! - 1.2 evaporation
!! - 1.2.1 vevapnu_old
!! - 1.2.2 ae_ns new
!! - 1.3 transpiration
!! - 1.4 root sink
!! - 2 Verification
!! - 2.1 Check of mc
!! - 2.2 Check if the deconvolution is correct and conserves the fluxes
!! - 2.2.1 the precisol and evapnu
!! - 2.2.2 the transpiration and root sink
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_split_soil

  SUBROUTINE hydrol_split_soil (kjpindex, veget_max, soiltile, vevapnu, transpir, humrel,evap_bare_lim)
    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(in)       :: veget_max        !! max Vegetation map 
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile         !! Fraction of each soil tile (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: vevapnu          !! Bare soil evaporation
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: transpir         !! Transpiration
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)       :: humrel           !! Relative humidity
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evap_bare_lim    !!   

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                :: ji, jv, jsl, jst
    REAL(r_std), DIMENSION (kjpindex)             :: vevapnu_old
    REAL(r_std), DIMENSION (kjpindex)             :: tmp_check1
    REAL(r_std), DIMENSION (kjpindex)             :: tmp_check2
    REAL(r_std), DIMENSION (kjpindex,nstm)        :: tmp_check3
    LOGICAL                                       :: error=.FALSE. !! If true, exit in the end of subroutine

!_ ================================================================================================================================
    !
    !
    !! 1 split 2d variables into 3d variables, per soil type
    !
    !
    !! 1.1 precipitation
    precisol_ns(:,:)=zero
    DO jv=1,nvm
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF(veget_max(ji,jv).GT.min_sechiba) THEN
                precisol_ns(ji,jst)=precisol_ns(ji,jst)+precisol(ji,jv)* &
                     & corr_veg_soil(ji,jv,jst) /vegtot(ji) / veget_max(ji,jv)
             ENDIF
          END DO
       END DO
    END DO

    !
    !
    !! 1.2 evaporation
    !! 1.2.1 vevapnu_old
    vevapnu_old(:)=zero
    DO jst=1,nstm
       DO ji=1,kjpindex
          IF ( vegtot(ji) .GT. min_sechiba) THEN
             vevapnu_old(ji)=vevapnu_old(ji)+ &
                  & ae_ns(ji,jst)*soiltile(ji,jst)*vegtot(ji)
          ENDIF
       END DO
    END DO
    !
    !
    !! 1.2.2 ae_ns new
    DO jst=1,nstm
       DO ji=1,kjpindex
          IF (vevapnu_old(ji).GT.min_sechiba) THEN   
             IF(evap_bare_lim(ji).GT.min_sechiba) THEN       
                ae_ns(ji,jst) = vevapnu(ji) * evap_bare_lim_ns(ji,jst)/evap_bare_lim(ji)
             ELSE
                IF(vevapnu_old(ji).GT.min_sechiba) THEN  
                   ae_ns(ji,jst)=ae_ns(ji,jst) * vevapnu(ji)/vevapnu_old(ji)
                ELSE
                   ae_ns(ji,jst)=zero
                ENDIF
             ENDIF
          ELSEIF(frac_bare_ns(ji,jst).GT.min_sechiba) THEN
             IF(evap_bare_lim(ji).GT.min_sechiba) THEN  
                ae_ns(ji,jst) = vevapnu(ji) * evap_bare_lim_ns(ji,jst)/evap_bare_lim(ji)
             ELSE
                IF(tot_bare_soil(ji).GT.min_sechiba) THEN  
                   ae_ns(ji,jst) = vevapnu(ji) * frac_bare_ns(ji,jst)/tot_bare_soil(ji)
                ELSE
                   ae_ns(ji,jst) = zero
                ENDIF
             ENDIF
          ENDIF
          precisol_ns(ji,jst)=precisol_ns(ji,jst)+MAX(-ae_ns(ji,jst),zero)
       END DO
    END DO
    !
    !
    !! 1.3 transpiration
    tr_ns(:,:)=zero
    DO jv=1,nvm
       DO jst=1,nstm
          DO ji=1,kjpindex
             IF (humrel(ji,jv).GT.min_sechiba) THEN 
                tr_ns(ji,jst)=tr_ns(ji,jst)+ cvs_over_veg(ji,jv,jst)*humrelv(ji,jv,jst)* & 
                     & transpir(ji,jv)/humrel(ji,jv)
             ENDIF
          END DO
       END DO
    END DO

    !
    !
    !! 1.4 root sink
    rootsink(:,:,:)=zero
    DO jv=1,nvm
       DO jsl=1,nslm
          DO jst=1,nstm
             DO ji=1,kjpindex
                IF (humrel(ji,jv).GT.min_sechiba) THEN 
                   rootsink(ji,jsl,jst) = rootsink(ji,jsl,jst) &
                        & + cvs_over_veg(ji,jv,jst)* (transpir(ji,jv)*us(ji,jv,jst,jsl))/ &
                        & humrel(ji,jv)
                END IF
             END DO
          END DO
       END DO
    END DO

    !! 2 Verification
    !! 2.1 Check of mc
    IF(check_cwrr) THEN
       DO jsl=1,nslm
          DO jst=1,nstm
             DO ji=1,kjpindex
                IF(mc(ji,jsl,jst).LT.-0.05) THEN
                   WRITE(numout,*) 'CWRR split-----------------------------------------------'
                   WRITE(numout,*) 'ji,jst,jsl',ji,jst,jsl
                   WRITE(numout,*) 'mc',mc(ji,jsl,jst)
                   WRITE(numout,*) 'rootsink,us',rootsink(ji,:,jst),us(ji,:,jst,jsl)
                   WRITE(numout,*) 'corr_veg_soil',corr_veg_soil(ji,:,jst)
                   WRITE(numout,*) 'transpir',transpir(ji,:)
                   WRITE(numout,*) 'veget_max',veget_max(ji,:)
                   WRITE(numout,*) 'cvs_over_veg',cvs_over_veg(ji,:,jst)
                   WRITE(numout,*) 'humrel',humrel(ji,:)
                   WRITE(numout,*) 'humrelv (pour ce jst)',humrelv(ji,:,jst)
                   WRITE(numout,*) 'ae_ns',ae_ns(ji,jst)
                   WRITE(numout,*) 'tr_ns',tr_ns(ji,jst)
                   WRITE(numout,*) 'vevapnuold',vevapnu_old(ji)
                ENDIF
             END DO
          END DO
       END DO
    ENDIF


    !! 2.2 Check if the deconvolution is correct and conserves the fluxes

    IF (check_cwrr) THEN


       tmp_check1(:)=zero
       tmp_check2(:)=zero  

    !! 2.2.1 the precisol and evapnu

       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + &
                  & (precisol_ns(ji,jst)-MAX(-ae_ns(ji,jst),zero))* &
                  & soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO

       DO jv=1,nvm
          DO ji=1,kjpindex
             tmp_check2(ji)=tmp_check2(ji) + precisol(ji,jv)
          END DO
       END DO


       DO ji=1,kjpindex   

          IF(ABS(tmp_check1(ji)- tmp_check2(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'PRECISOL SPLIT FALSE:ji=',ji,tmp_check1(ji),tmp_check2(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- tmp_check2(ji))
             WRITE(numout,*) 'vegtot',vegtot(ji)

             DO jv=1,nvm
                WRITE(numout,'(a,i2.2,"|",F13.4,"|",F13.4,"|",3(F9.6))') 'jv,veget_max, precisol, corr_veg_soil ',&
                     jv,veget_max(ji,jv),precisol(ji,jv),corr_veg_soil(ji,jv,:)
             END DO

             DO jst=1,nstm
                WRITE(numout,*) 'jst,precisol_ns',jst,precisol_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO

             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','PRECISOL SPLIT FALSE')
          ENDIF

       END DO


       tmp_check1(:)=zero

       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + ae_ns(ji,jst)* &
                  & soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO

       DO ji=1,kjpindex   

          IF(ABS(tmp_check1(ji)- vevapnu(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'VEVAPNU SPLIT FALSE:ji, Sum(ae_ns), vevapnu =',ji,tmp_check1(ji),vevapnu(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- vevapnu(ji))
             WRITE(numout,*) 'ae_ns',ae_ns(ji,:)
             WRITE(numout,*) 'vegtot',vegtot(ji)
             WRITE(numout,*) 'evap_bare_lim, evap_bare_lim_ns',evap_bare_lim(ji), evap_bare_lim_ns(ji,:)
             WRITE(numout,*) 'tot_bare_soil,frac_bare_ns',tot_bare_soil(ji),frac_bare_ns(ji,:)
             WRITE(numout,*) 'vevapnu_old',vevapnu_old(ji)
             DO jst=1,nstm
                WRITE(numout,*) 'jst,ae_ns',jst,ae_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
                WRITE(numout,*) 'veget_max/vegtot/soiltile', veget_max(ji,:)/vegtot(ji)/soiltile(ji,jst)
                WRITE(numout,*) "corr_veg_soil",corr_veg_soil(ji,:,jst)
             END DO

             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','VEVAPNU SPLIT FALSE')
          ENDIF
       ENDDO

    !! 2.2.2 the transpiration and root sink

       tmp_check1(:)=zero
       tmp_check2(:)=zero  


       DO jst=1,nstm
          DO ji=1,kjpindex
             tmp_check1(ji)=tmp_check1(ji) + tr_ns(ji,jst)* &
                  & soiltile(ji,jst)*vegtot(ji)
          END DO
       END DO

       DO jv=1,nvm
          DO ji=1,kjpindex
             tmp_check2(ji)=tmp_check2(ji) + transpir(ji,jv)
          END DO
       END DO

       DO ji=1,kjpindex   

          IF(ABS(tmp_check1(ji)- tmp_check2(ji)).GT.allowed_err) THEN
             WRITE(numout,*) 'TRANSPIR SPLIT FALSE:ji=',ji,tmp_check1(ji),tmp_check2(ji)
             WRITE(numout,*) 'err',ABS(tmp_check1(ji)- tmp_check2(ji))
             WRITE(numout,*) 'vegtot',vegtot(ji)

             DO jv=1,nvm
                WRITE(numout,*) 'jv,veget_max, transpir',jv,veget_max(ji,jv),transpir(ji,jv)
                DO jst=1,nstm
                   WRITE(numout,*) 'corr_veg_soil:ji,jv,jst',ji,jv,jst,corr_veg_soil(ji,jv,jst)
                END DO
             END DO

             DO jst=1,nstm
                WRITE(numout,*) 'jst,tr_ns',jst,tr_ns(ji,jst)
                WRITE(numout,*) 'soiltile', soiltile(ji,jst)
             END DO

             error=.TRUE.
             CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','TRANSPIR SPLIT FALSE')
          ENDIF

       END DO


       tmp_check3(:,:)=zero

       DO jst=1,nstm
          DO jsl=1,nslm
             DO ji=1,kjpindex
                tmp_check3(ji,jst)=tmp_check3(ji,jst) + rootsink(ji,jsl,jst)
             END DO
          END DO
       ENDDO

       DO jst=1,nstm
          DO ji=1,kjpindex
             IF(ABS(tmp_check3(ji,jst)- tr_ns(ji,jst)).GT.allowed_err) THEN
                WRITE(numout,*) 'ROOTSINK SPLIT FALSE:ji,jst=', ji,jst,&
                     & tmp_check3(ji,jst),tr_ns(ji,jst)
                WRITE(numout,*) 'err',ABS(tmp_check3(ji,jst)- tr_ns(ji,jst))
                WRITE(numout,*) 'HUMREL(jv=1:13)',humrel(ji,:)
                WRITE(numout,*) 'TRANSPIR',transpir(ji,:)
                DO jv=1,nvm 
                   WRITE(numout,*) 'jv=',jv,'us=',us(ji,jv,jst,:)
                ENDDO

                error=.TRUE.
                CALL ipslerr_p(2, 'hydrol_split_soil', 'We will STOP in the end of this subroutine.',&
                  & 'check_CWRR','ROOTSINK SPLIT FALSE')
             ENDIF
          END DO
       END DO

    ENDIF

!! Exit if error was found previously in this subroutine
    IF ( error ) THEN
       WRITE(numout,*) 'One or more errors have been detected in hydrol_split_soil. Model stops.'
       CALL ipslerr_p(3, 'hydrol_split_soil', 'We will STOP now.',&
                  & 'One or several fatal errors were found previously.','')
    END IF

  END SUBROUTINE hydrol_split_soil


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_diag_soil
!!
!>\BRIEF        ??
!!
!! DESCRIPTION  :
!! - 1 Apply mask_soiltile
!! - 2 sum 3d variables in 2d variables with fraction of vegetation per soil type
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_diag_soil

  SUBROUTINE hydrol_diag_soil (kjpindex, veget_max, soiltile, njsc, runoff, drainage, &
       & evapot, vevapnu, returnflow, reinfiltration, irrigation, &
       & shumdiag,shumdiag_perma, k_litt, litterhumdiag, humrel, vegstress, drysoil_frac, tot_melt) !! Arsene 28-01-2016 - REMOVE drunoff_tot because never user and bug in sechiba_output.f90
 !pss:+ !pss:-

    ! 
    ! interface description

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)        :: veget_max       !! Max. vegetation type
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: njsc            !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT (in)      :: soiltile        !! Fraction of each soil tile (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: evapot          !! 
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: returnflow      !! Water returning to the deep reservoir
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: reinfiltration  !! Water returning to the top of the soil
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: irrigation      !! Water from irrigation
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)            :: tot_melt        !!

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT (out)          :: drysoil_frac    !! Function of litter wetness
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: runoff          !! complete runoff
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)           :: drainage        !! Drainage
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out)      :: shumdiag        !! relative soil moisture
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out)      :: shumdiag_perma  !! Percent of porosity filled with water (mc/mcs) used for the thermal computations
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: k_litt          !! litter cond.
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: litterhumdiag   !! litter humidity
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: humrel          !! Relative humidity
    REAL(r_std), DIMENSION (kjpindex, nvm), INTENT(out)      :: vegstress       !! Veg. moisture stress (only for vegetation growth)

!pss:+
!    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)           :: drunoff_tot          !! Dunne runoff !! Arsene 28-01-2016 - REMOVE because never user and bug in sechiba_output.f90
!pss:-

    !! 0.3 Modified variables 

    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)         :: vevapnu         !!

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jv, jsl, jst, i, jd
    REAL(r_std), DIMENSION (kjpindex)                        :: mask_vegtot
    REAL(r_std)                                              :: k_tmp, tmc_litter_ratio

!_ ================================================================================================================================
    !
    ! Put the prognostics variables of soil to zero if soiltype is zero

    !! 1 Apply mask_soiltile
    DO jst=1,nstm 

       DO ji=1,kjpindex

!!??Aurelien: who and why coment this line
!          IF(soiltile(ji,jst).EQ.zero) THEN

             ae_ns(ji,jst) = ae_ns(ji,jst) * mask_soiltile(ji,jst)
             dr_ns(ji,jst) = dr_ns(ji,jst) * mask_soiltile(ji,jst)
             ru_ns(ji,jst) = ru_ns(ji,jst) * mask_soiltile(ji,jst)
             tmc(ji,jst) =  tmc(ji,jst) * mask_soiltile(ji,jst)
             IF (ok_freeze_cwrr) THEN
                profil_froz_hydro_ns(ji,:,jst)=profil_froz_hydro_ns(ji,:,jst)*mask_soiltile(ji,jst)
             ENDIF !if (ok_freeze_cwrr) then

             DO jv=1,nvm
                humrelv(ji,jv,jst) = humrelv(ji,jv,jst) * mask_soiltile(ji,jst)
                DO jsl=1,nslm
                   us(ji,jv,jst,jsl) = us(ji,jv,jst,jsl)  * mask_soiltile(ji,jst)
                END DO
             END DO

             DO jsl=1,nslm          
                mc(ji,jsl,jst) = mc(ji,jsl,jst)  * mask_soiltile(ji,jst)
             END DO

!          ENDIF

       END DO
    END DO

    runoff(:) = zero
    drainage(:) = zero
    humtot(:) = zero
    shumdiag(:,:)= zero
    shumdiag_perma(:,:)=zero
    k_litt(:) = zero
    litterhumdiag(:) = zero
    tmc_litt_mea(:) = zero
    humrel(:,:) = zero
    vegstress(:,:) = zero
    swi(:) = zero
    IF (ok_freeze_cwrr) THEN
       profil_froz_hydro(:,:)=zero
    ENDIF
    
    !! 2 sum 3d variables in 2d variables with fraction of vegetation per soil type

    DO ji = 1, kjpindex
       mask_vegtot(ji) = 0
       IF(vegtot(ji) .GT. min_sechiba) THEN
          mask_vegtot(ji) = 1
       ENDIF
    END DO
    
    DO ji = 1, kjpindex 
       ! Here we weight ae_ns by the fraction of bare evaporating soil. 
       ! This is given by frac_bare_ns, taking into account bare soil under vegetation
       ae_ns(ji,:) = mask_vegtot(ji) * ae_ns(ji,:) * frac_bare_ns(ji,:)
    END DO

    DO jst = 1, nstm
       DO ji = 1, kjpindex 
          drainage(ji) = mask_vegtot(ji) * (drainage(ji) + vegtot(ji)*soiltile(ji,jst) * dr_ns(ji,jst))
          runoff(ji) = mask_vegtot(ji) *  (runoff(ji) +   vegtot(ji)*soiltile(ji,jst) * ru_ns(ji,jst)) &
               & + (1 - mask_vegtot(ji)) * (tot_melt(ji) + irrigation(ji) + returnflow(ji) + reinfiltration(ji))
          humtot(ji) = mask_vegtot(ji) * (humtot(ji) + soiltile(ji,jst) * tmc(ji,jst)) 
          IF (ok_freeze_cwrr) THEN 
             profil_froz_hydro(ji,:)=mask_vegtot(ji) * &
                  (profil_froz_hydro(ji,:) + soiltile(ji,jst) * profil_froz_hydro_ns(ji,:, jst))
          ENDIF
       END DO
    END DO

    ! we add the excess of snow sublimation to vevapnu

    DO ji = 1,kjpindex
       vevapnu(ji) = vevapnu (ji) + subsinksoil(ji)*vegtot(ji)
    END DO

    DO jst=1,nstm
       DO jv=1,nvm
          DO ji=1,kjpindex
             IF(veget_max(ji,jv).GT.min_sechiba) THEN
                vegstress(ji,jv)=vegstress(ji,jv)+vegstressv(ji,jv,jst)*soiltile(ji,jst) &
                     & * corr_veg_soil(ji,jv,jst) *vegtot(ji)/veget_max(ji,jv)
                vegstress(ji,jv)= MAX(vegstress(ji,jv),zero)
             ENDIF
          END DO
       END DO
    END DO

    cvs_over_veg(:,:,:) = zero
    DO jv=1,nvm
       DO ji=1,kjpindex
          IF(veget_max(ji,jv).GT.min_sechiba) THEN
             DO jst=1,nstm
                cvs_over_veg(ji,jv,jst) = corr_veg_soil(ji,jv,jst)/vegtot(ji) / veget_max(ji,jv)
             ENDDO
          ENDIF
       END DO
    END DO

    DO jst=1,nstm
       DO jv=1,nvm
          DO ji=1,kjpindex
             humrel(ji,jv)=humrel(ji,jv)+humrelv(ji,jv,jst)*soiltile(ji,jst) &
                  & * cvs_over_veg(ji,jv,jst)*vegtot(ji)
             humrel(ji,jv)=MAX(humrel(ji,jv),zero)
          END DO
       END DO
    END DO

    DO jst=1,nstm
       DO ji=1,kjpindex
          ! We compute here a mean k for the 'litter' used for reinfiltration from floodplains of ponds
          !
          IF ( tmc_litter(ji,jst) < tmc_litter_res(ji,jst)) THEN
             i = imin
          ELSE
             tmc_litter_ratio = (tmc_litter(ji,jst)-tmc_litter_res(ji,jst)) / &
                  & (tmc_litter_sat(ji,jst)-tmc_litter_res(ji,jst))
             i= MAX(MIN(INT((imax-imin)*tmc_litter_ratio)+imin, imax-1), imin)
          ENDIF
          !
          !
          k_tmp = MAX(k_lin(i,1,njsc(ji))*ks(njsc(ji)), zero)
          k_litt(ji) = k_litt(ji) + soiltile(ji,jst) * SQRT(k_tmp)
       ENDDO
    ENDDO

    DO jst=1,nstm        

       DO ji=1,kjpindex
          litterhumdiag(ji) = litterhumdiag(ji) + &
               & soil_wet_litter(ji,jst) * soiltile(ji,jst)

          tmc_litt_mea(ji) = tmc_litt_mea(ji) + &
               & tmc_litter(ji,jst) * soiltile(ji,jst) 

       END DO

       DO jd=1,nbdl
          DO ji=1,kjpindex
               DO jsl=1,nslm  
                  shumdiag(ji,jd)= shumdiag(ji,jd) + soil_wet(ji,jsl,jst)  &
                       *frac_hydro_diag(jsl,jd)* &
                       ((mcs(njsc(ji))-mcw(njsc(ji)))/(mcf(njsc(ji))-mcw(njsc(ji)))) * &
                       soiltile(ji,jst)
                  
                  shumdiag_perma(ji,jd)= shumdiag_perma(ji,jd)  &
                       + mc(ji,jsl,jst) *frac_hydro_diag(jsl,jd) &
                       /mcs(jst)*soiltile(ji,jst)
               ENDDO
               shumdiag(ji,jd) = MAX(MIN(shumdiag(ji,jd), un), zero) 
               shumdiag_perma(ji,jd) = MAX(MIN(shumdiag_perma(ji,jd), un), zero) 
          END DO
       END DO

    END DO

    ! Calculate soilmoist
    soilmoist(:,:) = zero
    DO jst=1,nstm
       DO ji=1,kjpindex
             soilmoist(ji,1) = soilmoist(ji,1) + soiltile(ji,jst) * &
                  dz(2,jst) * ( trois*mc(ji,1,jst) + mc(ji,2,jst) )/huit
             DO jsl = 2,nslm-1
                soilmoist(ji,jsl) = soilmoist(ji,jsl) + soiltile(ji,jst) * &
                     ( dz(jsl,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl-1,jst))/huit &
                     + dz(jsl+1,jst) * (trois*mc(ji,jsl,jst)+mc(ji,jsl+1,jst))/huit )
             END DO
             soilmoist(ji,nslm) = soilmoist(ji,nslm) + soiltile(ji,jst) * &
                  dz(nslm,jst) * (trois*mc(ji,nslm,jst) + mc(ji,nslm-1,jst))/huit
       END DO
    END DO


    ! First we compute swi (ALMIP requirement) - we assume here that dz is independant of jst (jst=1)
    jst=1
    DO ji=1,kjpindex
       swi(ji) = swi(ji) + shumdiag(ji,1) * (dz(2,jst))/(deux*dpu_max*mille)
       
       DO jsl=2,nbdl-1 
          swi(ji) = swi(ji) + shumdiag(ji,jsl) * (dz(jsl,jst)+dz(jsl+1,jst))/(deux*dpu_max*mille)
       ENDDO
       swi(ji) = swi(ji) + shumdiag(ji,nbdl) * (dz(nbdl,jst))/(deux*dpu_max*mille)
    END DO

    DO ji=1,kjpindex
       IF ( tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji) > zero ) THEN
          drysoil_frac(ji) = un + MAX( MIN( (tmc_litt_dry_mea(ji) - tmc_litt_mea(ji)) / &
               & (tmc_litt_wet_mea(ji) - tmc_litt_dry_mea(ji)), zero), - un)
       ELSE
          drysoil_frac(ji) = zero
       ENDIF
    END DO

  END SUBROUTINE hydrol_diag_soil


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_waterbal 
!!
!>\BRIEF        Checks the water balance.
!!
!! DESCRIPTION  :
!! This routine checks the water balance. First it gets the total
!! amount of water and then it compares the increments with the fluxes.
!! The computation is only done over the soil area as over glaciers (and lakes?)
!! we do not have water conservation.
!! This verification does not make much sense in REAL*4 as the precision is the same as some
!! of the fluxes
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_waterbal

  SUBROUTINE hydrol_waterbal (kjpindex, index, first_call, dtradia, veget_max, totfrac_nobio, &
       & qsintveg, snow,snow_nobio, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, tot_melt, &
       & vevapwet, transpir, vevapnu, vevapsno, vevapflo, floodout, runoff, drainage)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                        :: kjpindex     !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index        !! Indeces of the points on the map
    LOGICAL, INTENT (in)                               :: first_call   !! At which time is this routine called ?
    REAL(r_std), INTENT (in)                           :: dtradia      !! Time step in seconds
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: veget_max    !! Max Fraction of vegetation type 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: totfrac_nobio!! Total fraction of continental ice+lakes+...
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg     !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow         !! Snow mass [Kg/m^2]
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio !!Ice water balance
    !
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_rain  !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: precip_snow  !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: returnflow   !! Water to the bottom
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: reinfiltration !! Water to the top
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: irrigation   !! Water from irrigation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: tot_melt     !! Total melt
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: vevapwet     !! Interception loss
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: transpir     !! Transpiration
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: vevapnu      !! Bare soil evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: vevapsno     !! Snow evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: vevapflo     !! Floodplains evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: floodout     !! flow out of floodplains
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: runoff       !! complete runoff
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: drainage     !! Drainage

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: ji
    REAL(r_std) :: watveg, delta_water
    LOGICAL     :: error=.FALSE.  !! If true, exit in the end of subroutine

!_ ================================================================================================================================
    !
    !
    !
    IF ( ALL( tot_water_beg(:) == val_exp ) ) THEN

       DO ji = 1, kjpindex
          watveg = SUM(qsintveg(ji,:))
          tot_water_beg(ji) = humtot(ji)*vegtot(ji) + watveg + snow(ji)&
               & + SUM(snow_nobio(ji,:))
       ENDDO

       tot_water_end(:) = tot_water_beg(:)
       tot_flux(:) = zero

       RETURN
    ELSE IF ( first_call ) THEN
       tot_water_end(:) = tot_water_beg(:)
       tot_flux(:) = zero

       RETURN
    ENDIF

    tot_water_end(:) = zero
    tot_flux(:) = zero
    !
    DO ji = 1, kjpindex
       !
       ! If the fraction of ice, lakes, etc. does not complement the vegetation fraction then we do not
       ! need to go any further
       !
       IF ( ABS(un - (totfrac_nobio(ji) + vegtot(ji))) .GT. allowed_err ) THEN
          WRITE(numout,*) 'HYDROL problem in vegetation or frac_nobio on point ', ji
          WRITE(numout,*) 'totfrac_nobio : ', totfrac_nobio(ji)
          WRITE(numout,*) 'vegetation fraction : ', vegtot(ji)

          error=.TRUE.
          CALL ipslerr_p(2, 'hydrol_waterbal', 'We will STOP in the end of hydrol_waterbal.','','')
       ENDIF
    ENDDO

    DO ji = 1, kjpindex
       !
       watveg = SUM(qsintveg(ji,:))
       tot_water_end(ji) = humtot(ji)*vegtot(ji) + watveg + &
            & snow(ji) + SUM(snow_nobio(ji,:))
       !
       tot_flux(ji) =  precip_rain(ji) + precip_snow(ji) + irrigation (ji) - &
            & SUM(vevapwet(ji,:)) - SUM(transpir(ji,:)) - vevapnu(ji) - vevapsno(ji) - vevapflo(ji) + &
            & floodout(ji) - runoff(ji) - drainage(ji) + returnflow(ji) + reinfiltration(ji)
    ENDDO
    
    DO ji = 1, kjpindex
       !
       delta_water = tot_water_end(ji) - tot_water_beg(ji)
       !
       !
       !  Set some precision ! This is a wild guess and corresponds to what works on an IEEE machine
       !  under double precision (REAL*8).
       !
       !
       IF ( ABS(delta_water-tot_flux(ji)) .GT. deux*allowed_err ) THEN
          WRITE(numout,*) '------------------------------------------------------------------------- '
          WRITE(numout,*) 'HYDROL does not conserve water. The erroneous point is : ', ji
          WRITE(numout,*) 'Coord erroneous point', lalo(ji,:)
          WRITE(numout,*) 'The error in mm/s is :', (delta_water-tot_flux(ji))/dtradia, ' and in mm/dt : ', &
               & delta_water-tot_flux(ji)
          WRITE(numout,*) 'delta_water : ', delta_water, ' tot_flux : ', tot_flux(ji)
          WRITE(numout,*) 'Actual and allowed error : ', ABS(delta_water-tot_flux(ji)), allowed_err
          WRITE(numout,*) 'vegtot : ', vegtot(ji)
          WRITE(numout,*) 'precip_rain : ', precip_rain(ji)
          WRITE(numout,*) 'precip_snow : ',  precip_snow(ji)
          WRITE(numout,*) 'Water from routing. Reinfiltration/returnflow/irrigation : ', reinfiltration(ji), &
               & returnflow(ji),irrigation(ji)
          WRITE(numout,*) 'Total water in soil humtot:',  humtot(ji)
          WRITE(numout,*) 'mc:' , mc(ji,:,:)
          WRITE(numout,*) 'Water on vegetation watveg:', watveg
          WRITE(numout,*) 'Snow mass snow:', snow(ji)
          WRITE(numout,*) 'Snow mass on ice snow_nobio:', SUM(snow_nobio(ji,:))
          WRITE(numout,*) 'Melt water tot_melt:', tot_melt(ji)
          WRITE(numout,*) 'evapwet : ', vevapwet(ji,:)
          WRITE(numout,*) 'transpir : ', transpir(ji,:)
          WRITE(numout,*) 'evapnu, evapsno, evapflo: ', vevapnu(ji), vevapsno(ji), vevapflo(ji)
          WRITE(numout,*) 'drainage,runoff,floodout : ', drainage(ji),runoff(ji),floodout(ji)
         ! error=.TRUE.
         ! CALL ipslerr_p(2, 'hydrol_waterbal', 'We will STOP in the end of hydrol_waterbal.','','')
       ENDIF
       !
    ENDDO
    !
    ! Transfer the total water amount at the end of the current timestep top the begining of the next one.
    !
    tot_water_beg = tot_water_end
    !
    
    ! Exit if one or more errors were found
    IF ( error ) THEN
       WRITE(numout,*) 'One or more errors have been detected in hydrol_waterbal. Model stops.'
       CALL ipslerr_p(3, 'hydrol_waterbal', 'We will STOP now.',&
            'One or several fatal errors were found previously.','')
    END IF
    
  END SUBROUTINE hydrol_waterbal


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_alma 
!!
!>\BRIEF        This routine computes the changes in soil moisture and interception storage for the ALMA outputs.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!_ hydrol_alma

  SUBROUTINE hydrol_alma (kjpindex, index, first_call, qsintveg, snow, snow_nobio, soilwet)
    !
    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                        :: kjpindex     !! Domain size
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)   :: index        !! Indeces of the points on the map
    LOGICAL, INTENT (in)                              :: first_call   !! At which time is this routine called ?
    !
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)  :: qsintveg     !! Water on vegetation due to interception
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: snow         !! Snow water equivalent
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: snow_nobio     !! Water balance on ice, lakes, .. [Kg/m^2]

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: soilwet     !! Soil wetness

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std) :: ji
    REAL(r_std) :: watveg

!_ ================================================================================================================================
    !
    !
    IF ( first_call ) THEN

       tot_watveg_beg(:) = zero
       tot_watsoil_beg(:) = zero
       snow_beg(:)        = zero
       !
       DO ji = 1, kjpindex
          watveg = SUM(qsintveg(ji,:))
          tot_watveg_beg(ji) = watveg
          tot_watsoil_beg(ji) = humtot(ji)
          snow_beg(ji)        = snow(ji)+ SUM(snow_nobio(ji,:))
       ENDDO
       !
       tot_watveg_end(:) = tot_watveg_beg(:)
       tot_watsoil_end(:) = tot_watsoil_beg(:)
       snow_end(:)        = snow_beg(:)

       RETURN

    ENDIF
    !
    ! Calculate the values for the end of the time step
    !
    tot_watveg_end(:) = zero
    tot_watsoil_end(:) = zero
    snow_end(:) = zero
    delintercept(:) = zero
    delsoilmoist(:) = zero
    delswe(:) = zero
    !
    DO ji = 1, kjpindex
       watveg = SUM(qsintveg(ji,:))
       tot_watveg_end(ji) = watveg
       tot_watsoil_end(ji) = humtot(ji)
       snow_end(ji) = snow(ji)+ SUM(snow_nobio(ji,:))
       !
       delintercept(ji) = tot_watveg_end(ji) - tot_watveg_beg(ji)
       delsoilmoist(ji) = tot_watsoil_end(ji) - tot_watsoil_beg(ji)
       delswe(ji)       = snow_end(ji) - snow_beg(ji)
       !
       !
    ENDDO
    !
    !
    ! Transfer the total water amount at the end of the current timestep top the begining of the next one.
    !
    tot_watveg_beg = tot_watveg_end
    tot_watsoil_beg = tot_watsoil_end
    snow_beg(:) = snow_end(:)
    !
    DO ji = 1,kjpindex
       IF ( mx_eau_var(ji) > 0 ) THEN
          soilwet(ji) = tot_watsoil_end(ji) / mx_eau_var(ji)
       ELSE
          soilwet(ji) = zero
       ENDIF
    ENDDO
    !
  END SUBROUTINE hydrol_alma
  !


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_calculate_temp_hydro
!!
!>\BRIEF         Calculate the temperature at hydrological levels  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================


  SUBROUTINE hydrol_calculate_temp_hydro(kjpindex, stempdiag, snow,snowdz)

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                             :: kjpindex 
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)     :: stempdiag
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)          :: snow
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (in)    :: snowdz


    !! 0.2 Local variables
    
    INTEGER jh, jd, ji
    REAL(r_std) :: snow_h, profil_froz, m, mc_used
    REAL(r_std)  :: lev_diag, prev_diag, lev_prog, prev_prog
    REAL(r_std), DIMENSION(nslm,nbdl) :: intfactt
    
    
    DO ji=1,kjpindex
       IF (ok_explicitsnow) THEN 
          snow_h=SUM(snowdz(ji,:))
       ELSE  
          snow_h=snow(ji)/sn_dens
       ENDIF
       
       intfactt(:,:)=0.
       prev_diag = snow_h
       DO jh = 1, nslm
          IF (jh.EQ.1) THEN
             lev_diag = zz(2,nstm)/1000./2.+snow_h
          ELSEIF (jh.EQ.nslm) THEN
             lev_diag = zz(nslm,nstm)/1000.+snow_h
             
          ELSE
             lev_diag = zz(jh, nstm)/1000. &
                  & +(zz(jh+1,nstm)-zz(jh,nstm))/1000./2.+snow_h
             
          ENDIF
          prev_prog = 0.0
          DO jd = 1, nbdl
             lev_prog = diaglev(jd)
             IF ((lev_diag.GT.diaglev(nbdl).AND. &
                  & prev_diag.LT.diaglev(nbdl)-min_sechiba)) THEN
                lev_diag=diaglev(nbdl)          
             ENDIF
             intfactt(jh,jd) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog),&
                  & 0.0)/(lev_diag-prev_diag)
             prev_prog = lev_prog
          ENDDO
          IF (lev_diag.GT.diaglev(nbdl).AND. &
               & prev_diag.GE.diaglev(nbdl)-min_sechiba) intfactt(jh,nbdl)=1.
          prev_diag = lev_diag
       ENDDO
    ENDDO
    
    temp_hydro(:,:)=0.
    DO jd= 1, nbdl
       DO jh= 1, nslm
          DO ji = 1, kjpindex
             temp_hydro(ji,jh) = temp_hydro(ji,jh) + stempdiag(ji,jd)*intfactt(jh,jd)
          ENDDO
       ENDDO
    ENDDO
    
  END SUBROUTINE hydrol_calculate_temp_hydro


!! ================================================================================================================================
!! SUBROUTINE   : hydrol_calculate_frac_hydro_diag
!!
!>\BRIEF         Caluculate frac_hydro_diag for interpolation between hydrological and diagnostic axes
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE hydrol_calculate_frac_hydro_diag

    !! 0.1 Local variables

    INTEGER(i_std) :: jd, jh
    REAL(r_std)    :: prev_hydro, next_hydro, prev_diag, next_diag
    

    frac_hydro_diag(:,:)=0.
    prev_diag = 0.0
    
    DO jd = 1, nbdl 
       
       next_diag = diaglev(jd)
       prev_hydro = 0.0
       DO jh = 1, nslm
          IF (jh.EQ.1) THEN
             next_hydro = zz(2,nstm)/1000./2.
          ELSEIF (jh.EQ.nslm) THEN
             next_hydro = zz(nslm,nstm)/1000.
          ELSE
             next_hydro = zz(jh, nstm)/1000.+(zz(jh+1,nstm)-zz(jh,nstm))/1000./2.
          ENDIF
          frac_hydro_diag(jh,jd) = MAX(MIN(next_hydro, next_diag)-MAX(prev_hydro, prev_diag), 0.)/(next_diag - prev_diag)
          prev_hydro=next_hydro
       ENDDO
       
       prev_diag = next_diag
    ENDDO

  END SUBROUTINE hydrol_calculate_frac_hydro_diag
  
END MODULE hydrol
