! ================================================================================================================================
! MODULE       : stomate_lpj
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Main entry point for daily processes in STOMATE and LPJ (phenology, 
!! allocation, npp_calc, kill, turn, light, establish, crown, cover, lcchange)
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_lpj.f90 $
!! $Date: 2014-04-04 18:19:41 +0200 (Fri, 04 Apr 2014) $
!! $Revision: 2031 $
!! \n
!_ ================================================================================================================================

MODULE stomate_lpj

  ! modules used:

  USE ioipsl_para
  USE xios_orchidee
  USE grid
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE lpj_constraints
  USE lpj_pftinout
  USE lpj_kill
  USE lpj_crown
  USE lpj_fire
  USE lpj_gap
  USE lpj_light
  USE lpj_establish
  USE lpj_cover
  USE stomate_prescribe
  USE stomate_phenology
  USE stomate_alloc
  USE stomate_npp
  USE stomate_turnover
  USE stomate_litter
  USE stomate_soilcarbon
  USE stomate_vmax
  USE stomate_lcchange
!JCADD
  USE Grassland_Management
!ENDJCADD
!  USE Orch_Write_field_p

  !pss:+
!  USE stomate_cste_wetlands
  USE stomate_wet_ch4_pt_ter_0
  USE stomate_wet_ch4_pt_ter_wet1
  USE stomate_wet_ch4_pt_ter_wet2
  USE stomate_wet_ch4_pt_ter_wet3
  USE stomate_wet_ch4_pt_ter_wet4
  !pss:-



  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC StomateLpj,StomateLpj_clear

  LOGICAL, SAVE                         :: firstcall = .TRUE.             !! first call
!$OMP THREADPRIVATE(firstcall)
!JCADD
  ! flag that enable grazing
  LOGICAL, SAVE :: enable_grazing
  ! flag for cutting
  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE :: flag_cutting
  ! how many days ago was the beginning of the last cut
  REAL(r_std), DIMENSION(:,:), ALLOCATABLE, SAVE :: when_growthinit_cut
!ENDJCADD

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : StomateLpj_clear
!!
!>\BRIEF        Re-initialisation of variable
!!
!! DESCRIPTION  : This subroutine reinitializes variables. To be used if we want to relaunch 
!! ORCHIDEE but the routine is not used in current version.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE StomateLpj_clear

    CALL prescribe_clear
    CALL phenology_clear
    CALL npp_calc_clear
    CALL turn_clear
    CALL soilcarbon_clear
    CALL constraints_clear
    CALL establish_clear
    CALL fire_clear
    CALL gap_clear
    CALL light_clear
    CALL pftinout_clear
    CALL alloc_clear

    !pss:+
    CALL ch4_wet_flux_density_clear_0
    CALL ch4_wet_flux_density_clear_wet1
    CALL ch4_wet_flux_density_clear_wet2
    CALL ch4_wet_flux_density_clear_wet3
    CALL ch4_wet_flux_density_clear_wet4
    !pss:-

  END SUBROUTINE StomateLpj_clear


!! ================================================================================================================================
!! SUBROUTINE   : StomateLPJ
!!
!>\BRIEF        Main entry point for daily processes in STOMATE and LPJ, structures the call sequence 
!!              to the different processes such as dispersion, establishment, competition and mortality of PFT's.
!! 
!! DESCRIPTION  : This routine is the main entry point to all processes calculated on a 
!! daily time step. Is mainly devoted to call the different STOMATE and LPJ routines 
!! depending of the ok_dgvm (is dynamic veg used) and lpj_constant_mortality (is background mortality used).
!! It also prepares the cumulative 
!! fluxes or pools (e.g TOTAL_M TOTAL_BM_LITTER etc...)
!!
!! This routine makes frequent use of "weekly", "monthly" and "long term" variables. Quotion is used because
!! by default "weekly" denotes 7 days, by default "monthly" denotes 20 days and by default "Long term" denotes
!! 3 years. dtslow refers to 24 hours (1 day).
!!
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): All variables related to stomate and required for LPJ dynamic vegetation mode.
!!
!! REFERENCE(S) : 
!! - Krinner, G., N. Viovy, N. de Noblet-Ducoudré, J. Ogeé, J. Polcher, P. Friedlingstein, P. Ciais, S. Sitch, 
!! and I. C. Prentice. 2005. A dynamic global vegetation model for studies of the coupled atmosphere-biosphere 
!! system. Global Biogeochemical Cycles 19:GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, I. C. Prentice, A. Arneth, A. Bondeau, W. Cramer, J. O. Kaplan, S. Levis, W. Lucht, 
!! M. T. Sykes, K. Thonicke, and S. Venevsky. 2003. Evaluation of ecosystem dynamics, plant geography and 
!! terrestrial carbon cycling in the LPJ dynamic global vegetation model. Global Change Biology 9:161-185.
!!
!! FLOWCHART    : Update with existing flowchart from N Viovy (Jan 19, 2012)
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE StomateLpj (npts, dt_days, EndOfYear, EndOfMonth, &
       neighbours, resolution, &
       clay, herbivores, &
       tsurf_daily, tsoil_daily, t2m_daily, t2m_min_daily, &
       litterhum_daily, soilhum_daily, &
       maxmoiavail_lastyear, minmoiavail_lastyear, &
       gdd0_lastyear, precip_lastyear, &
       moiavail_month, moiavail_week, tlong_ref, t2m_month, t2m_week, &
       tsoil_month, soilhum_month, &
       gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       turnover_longterm, gpp_daily, &
       time_hum_min, maxfpc_lastyear, resp_maint_part, &
       PFTpresent, age, fireindex, firelitter, &
       leaf_age, leaf_frac, biomass, ind, adapted, regenerate, &
       senescence, when_growthinit, &
       litterpart, litter, litter_avail, litter_not_avail, litter_avail_frac, &
       dead_leaves, carbon,carbon_surf, black_carbon, lignin_struc, &
       veget_max, npp_longterm, lm_lastyearmax, veget_lastlight, &
       everywhere, need_adjacent, RIP_time, &
       lai, rprof,npp_daily, turnover_daily, turnover_time,&
       control_moist, control_temp, soilcarbon_input, &
       co2_to_bm, co2_fire, resp_hetero, resp_maint, resp_growth, &
       height, deadleaf_cover, vcmax, &
       bm_to_litter, &
       prod10,prod100,flux10, flux100, veget_max_new, &
       convflux,cflux_prod10,cflux_prod100, harvest_above, carb_mass_total, lcchange, &
       fpc_max, Tseason, Tseason_length, Tseason_tmp, &
       Tmin_spring, Tmin_spring_time, begin_leaves, onset_date, &
       MatrixA, npp0_cumul, snowtemp_min, snowdz_min, dia_cut, &  !! Arsene 25-06-2014 NPPcumul Add npp0_cumul !! Arsene 19-08-2014 Add snowtemp_min & snowdz_min) !! Arsene 27-08-2015 add dia_cut
       zz_coef_deep, deepC_a, deepC_s, deepC_p, & !pss:+
       ch4_flux_density_tot_0, ch4_flux_density_dif_0, ch4_flux_density_bub_0, ch4_flux_density_pla_0,&
       ch4_flux_density_tot_wet1,ch4_flux_density_dif_wet1,ch4_flux_density_bub_wet1,ch4_flux_density_pla_wet1,&
       ch4_flux_density_tot_wet2,ch4_flux_density_dif_wet2,ch4_flux_density_bub_wet2,ch4_flux_density_pla_wet2,&
       ch4_flux_density_tot_wet3,ch4_flux_density_dif_wet3,ch4_flux_density_bub_wet3,ch4_flux_density_pla_wet3,&
       ch4_flux_density_tot_wet4,ch4_flux_density_dif_wet4,ch4_flux_density_bub_wet4,ch4_flux_density_pla_wet4,&
       tsurf_year, &!) !pss:-
!JCADD
       wshtotsum, sr_ugb, compt_ugb, nb_ani, grazed_frac, &
       import_yield, sla_age1, t2m_14, sla_calc, snow, day_of_year,N_limfert)!, &
!       resp_hetero_litter_d,resp_hetero_soil_d)
!ENDJCADD
    
  !! 0. Variable and parameter declaration

    !! 0.1 input

    INTEGER(i_std), INTENT(in)                                 :: npts                 !! Domain size (unitless)
    REAL(r_std), INTENT(in)                                    :: dt_days              !! Time step of Stomate (days)
    INTEGER(i_std), DIMENSION(npts,8), INTENT(in)              :: neighbours           !! Indices of the 8 neighbours of each grid 
                                                                                       !! point [1=N, 2=NE, 3=E, 4=SE,
                                                                                       !!  5=S, 6=SW, 7=W, 8=NW] 
    REAL(r_std), DIMENSION(npts,2), INTENT(in)                 :: resolution           !! Resolution at each grid point (m)  
                                                                                       !! [1=E-W, 2=N-S] 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: clay                 !! Clay fraction (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: herbivores           !! Time constant of probability of a leaf to 
                                                                                       !! be eaten by a herbivore (days) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tsurf_daily          !! Daily surface temperatures (K)  !! Arsene 20-08-2014 N'est pas utilisé
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)              :: tsoil_daily          !! Daily soil temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_daily            !! Daily 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_min_daily        !! Daily minimum 2 meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: litterhum_daily      !! Daily litter humidity (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)              :: soilhum_daily        !! Daily soil humidity (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxmoiavail_lastyear !! Last year's maximum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: minmoiavail_lastyear !! Last year's minimum moisture availability 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: gdd0_lastyear        !! Last year's GDD0 (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: precip_lastyear      !! Lastyear's precipitation 
                                                                                       !! @tex $(mm year^{-1})$ @endtex
										       !! to determine if establishment possible
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_month       !! "Monthly" moisture availability (0 to 1, 
                                                                                       !! unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: moiavail_week        !! "Weekly" moisture availability 
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tlong_ref            !! "Long term" 2 meter reference 
                                                                                       !! temperatures (K) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_month            !! "Monthly" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: t2m_week             !! "Weekly" 2-meter temperatures (K)
    ! "seasonal" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: Tseason
    ! temporary variable to calculate Tseason
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: Tseason_length
    ! temporary variable to calculate Tseason
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: Tseason_tmp

    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: Tmin_spring
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: Tmin_spring_time
    REAL(r_std), DIMENSION(npts,nvm,2), INTENT(in)             :: onset_date

    !pss:+
    REAL(r_std), DIMENSION(npts), INTENT(in)                   :: tsurf_year           ! annual surface temperatures (K)
    !pss:-
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)              :: tsoil_month          !! "Monthly" soil temperatures (K)
    REAL(r_std), DIMENSION(npts,nbdl), INTENT(in)              :: soilhum_month        !! "Monthly" soil humidity
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: gdd_m5_dormance      !! Growing degree days (K), threshold -5 deg 
                                                                                       !! C (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: gdd_from_growthinit  !! growing degree days, since growthinit for crops
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: gdd_midwinter        !! Growing degree days (K), since midwinter 
                                                                                       !! (for phenology) - this is written to the history files 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: ncd_dormance         !! Number of chilling days (days), since 
                                                                                       !! leaves were lost (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: ngd_minus5           !! Number of growing days (days), threshold 
                                                                                       !! -5 deg C (for phenology) 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)        :: turnover_longterm    !! "Long term" turnover rate  
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: gpp_daily            !! Daily gross primary productivity  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: time_hum_min         !! Time elapsed since strongest moisture 
                                                                                       !! availability (days) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: maxfpc_lastyear      !! Last year's maximum foliage projected
                                                                                       !! coverage for each natural PFT,
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts), INTENT(in)        :: resp_maint_part      !! Maintenance respiration of different 
                                                                                       !! plant parts  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: fpc_max              !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    LOGICAL, INTENT(in)                                        :: lcchange             !! Land cover change flag
    LOGICAL, INTENT(in)                                        :: EndOfYear            !! Flag set on the last day of the year used 
                                                                                       !! to update "yearly variables". This 
                                                                                       !! variable must be .TRUE. once a year
    LOGICAL, INTENT(in)                                        :: EndOfMonth           !! Flag set at end of each month to update 
                                                                                       !! monthly variable 
    REAL(r_std), DIMENSION(ndeep),   INTENT (in)               :: zz_coef_deep         !! deep vertical profile

    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)               :: npp0_cumul           !! Arsene 25-06-2014 NPPcumul
!! Arsene 25-06-2014 NPPcumul - Variable to count number of days since npp =0 or <0. Could become DIMENSION(npts:nvm) if we add the same for others pft
    REAL(r_std), DIMENSION(npts,nsnow), INTENT(in)              :: snowtemp_min        !! Min daily snow layer temperature  !! Arsene 19-08-2014 Add
    REAL(r_std), DIMENSION(npts,nsnow), INTENT(in)              :: snowdz_min          !! Min daily snow layer thicknesse   !! Arsene 19-08-2014 Add



  !! 0.2 Output variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: npp_daily            !! Net primary productivity 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(out) :: turnover_daily       !! Turnover rates 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: co2_to_bm            !! CO2 taken up from atmosphere when 
                                                                                       !! introducing a new PFT (introduced for 
                                                                                       !! carbon balance closure) 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: co2_fire             !! Carbon emitted into the atmosphere by 
                                                                                       !! fire (living and dead biomass)  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: resp_hetero          !! Heterotrophic respiration
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: resp_maint           !! Maintenance respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: resp_growth          !! Growth respiration  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: deadleaf_cover       !! Fraction of soil covered by dead leaves 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)              :: vcmax                !! Maximum rate of carboxylation 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(out):: bm_to_litter      !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex
    LOGICAL, DIMENSION(npts,nvm), INTENT(out)                  :: begin_leaves     !! signal to start putting leaves on (true/false)

    ! Wetland CH4 methane density
    !pss:+
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_tot_0
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_dif_0
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_bub_0
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_pla_0

    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_tot_wet1
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_dif_wet1
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_bub_wet1
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_pla_wet1

    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_tot_wet2
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_dif_wet2
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_bub_wet2
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_pla_wet2

    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_tot_wet3
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_dif_wet3
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_bub_wet3
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_pla_wet3

    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_tot_wet4
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_dif_wet4
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_bub_wet4
    REAL(r_std), DIMENSION(npts), INTENT(in)             :: ch4_flux_density_pla_wet4
    !pss:-
 

    !! 0.3 Modified variables
    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: height               !! Height of vegetation (m) 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(inout)          :: control_moist        !! Moisture control of heterotrophic 
                                                                                       !! respiration (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nlevs), INTENT(inout)          :: control_temp         !! Temperature control of heterotrophic 
                                                                                       !! respiration, above and below 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: soilcarbon_input     !! Quantity of carbon going into carbon 
                                                                                       !! pools from litter decomposition  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lai                  !! Leaf area index OF AN INDIVIDUAL PLANT,
										       !! where a PFT contains n indentical plants
										       !! i.e., using the mean individual approach 
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: rprof                !! Prescribed root depth (m) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: PFTpresent           !! Tab indicating which PFTs are present in 
                                                                                       !! each pixel 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: age                  !! Age (years)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: dia_cut              !! Fix diameter of vegetation (for shrub) after loss biomass (above snow) !! Arsene 27-08-2015 add dia_cut
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: fireindex            !! Probability of fire (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: firelitter           !! Longer term litter above the ground that 
                                                                                       !! can be burned, @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_age             !! Leaf age (days)
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout)  :: leaf_frac            !! Fraction of leaves in leaf age class, 
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: ind                  !! Density of individuals 
                                                                                       !! @tex $(m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: adapted              !! Adaptation of PFT (killed if too cold) 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: regenerate           !! "Fitness": Winter sufficiently cold for 
                                                                                       !! PFT regeneration ? (0 to 1, unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: senescence           !! Flag for setting senescence stage (only 
                                                                                       !! for deciduous trees) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: when_growthinit      !! How many days ago was the beginning of 
                                                                                       !! the growing season (days) 
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: litterpart           !! Fraction of litter above the ground 
                                                                                       !! belonging to different PFTs
                                                                                       !! (0 to 1, unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter     !! Metabolic and structural litter, above 
                                                                                       !! and below ground 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
!JCADD for grazing litter
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(out):: litter_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm) , INTENT(out):: litter_not_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(in):: litter_avail_frac
!ENDJCADD
    REAL(r_std), DIMENSION(npts,nvm,nlitt), INTENT(inout)      :: dead_leaves          !! Dead leaves on ground, per PFT, metabolic 
                                                                                       !! and structural,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon               !! Carbon pool: active, slow, or passive, 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon_surf
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: black_carbon         !! Black carbon on the ground 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nlevs), INTENT(inout)      :: lignin_struc         !! Ratio of Lignin/Carbon in structural 
                                                                                       !! litter, above and below ground,  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_max            !! "Maximal" coverage fraction of a PFT (LAI 
                                                                                       !! -> infinity) on ground 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: npp_longterm         !! "Long term" mean yearly primary 
                                                                                       !! productivity 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: lm_lastyearmax       !! Last year's maximum leaf mass, for each 
                                                                                       !! PFT @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: veget_lastlight      !! Vegetation fractions (on ground) after 
                                                                                       !! last light competition  
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: everywhere           !! Is the PFT everywhere in the grid box or 
                                                                                       !! very localized (after its introduction) 
                                                                                       !! (unitless) 
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)                :: need_adjacent        !! In order for this PFT to be introduced, 
                                                                                       !! does it have to be present in an 
                                                                                       !! adjacent grid box? 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: RIP_time             !! How much time ago was the PFT eliminated 
                                                                                       !! for the last time (y) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)            :: turnover_time        !! Turnover_time of leaves for grasses 
                                                                                       !! (days)
    REAL(r_std), DIMENSION(npts,nvm),INTENT(inout)             :: veget_max_new        !! New "maximal" coverage fraction of a PFT 
                                                                                       !! (LAI -> infinity) (unitless) 
    REAL(r_std),DIMENSION(npts,0:10), INTENT(inout)            :: prod10               !! Products remaining in the 10
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (10
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,0:100), INTENT(inout)           :: prod100              !! Products remaining in the 100 
                                                                                       !! year-turnover pool after the annual 
                                                                                       !! release for each compartment (100 
                                                                                       !! + 1 : input from year of land cover 
                                                                                       !! change) @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,10), INTENT(inout)              :: flux10               !! Annual release from the 10
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts,100), INTENT(inout)             :: flux100              !! Annual release from the 100 
                                                                                       !! year-turnover pool compartments  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: convflux             !! Release during first year following land 
                                                                                       !! cover change @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: cflux_prod10         !! Total annual release from the 10 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts), INTENT(inout)                 :: cflux_prod100        !! Total annual release from the 100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: harvest_above        !! Harvest above ground biomass for 
                                                                                       !! agriculture @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)                :: carb_mass_total      !! Carbon Mass total (soil, litter, veg) 
                                                                                       !! @tex $(gC m^{-2})$ @endtex  
    REAL(r_std), DIMENSION(npts,nvm,nbpools,nbpools), INTENT(inout) :: MatrixA         !! Matrix containing the fluxes  
                                                                                       !! between the carbon pools
                                                                                       !! per sechiba time step 
                                                                                       !! @tex $(gC.m^2.day^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_a           !! permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_s           !! permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_p           !! permafrost soil carbon (g/m**3) passive
!JCADD
   ! snow mass (kg/m2)
    REAL(r_std), DIMENSION(npts), INTENT(in)         :: snow
    ! "14days" 2-meter temperatures (K)
    REAL(r_std), DIMENSION(npts), INTENT(in)         ::  t2m_14
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  sla_calc
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  wshtotsum
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  sr_ugb
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  compt_ugb
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  nb_ani
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  grazed_frac
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  import_yield
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  sla_age1
    INTEGER(i_std), INTENT(in)                       ::  day_of_year
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  N_limfert
!    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  resp_hetero_litter_d
!    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)  :: resp_hetero_soil_d
!ENDJCADD
    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_bm_to_litter    !! Total conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_live_biomass    !! Total living biomass  
                                                                                       !! @tex $(gC m{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements)           :: bm_alloc            !! Biomass increase, i.e. NPP per plant part 
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nelements)                  :: tot_turnover        !! Total turnover rate  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_litter_soil_carb!! Total soil and litter carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_litter_carb     !! Total litter carbon 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: tot_soil_carb       !! Total soil carbon  
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std), DIMENSION(npts)                                :: carb_mass_variation !! Carbon Mass variation  
                                                                                       !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: cn_ind              !! Crown area of individuals 
                                                                                       !! @tex $(m^{2})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm)                            :: woodmass_ind        !! Woodmass of individuals (gC) 
    REAL(r_std), DIMENSION(npts,nvm,nparts)                     :: f_alloc             !! Fraction that goes into plant part 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: avail_tree          !! Space availability for trees 
                                                                                       !! (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: avail_grass         !! Space availability for grasses 
                                                                                       !! (0 to 1, unitless) 
    INTEGER                                                     :: j,ji,i,m            !! Arsene 20-08-2014 Add ji & i & m
    REAL(r_std),DIMENSION(npts)                                 :: prod10_total        !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: prod100_total       !! Total products remaining in the pool 
                                                                                       !! after the annual release 
                                                                                       !! @tex $(gC m^{-2})$ @endtex 
    REAL(r_std),DIMENSION(npts)                                 :: cflux_prod_total    !! Total flux from conflux and the 10/100 
                                                                                       !! year-turnover pool 
                                                                                       !! @tex $(gC m^{-2} year^{-1})$ @endtex 
    REAL(r_std),DIMENSION(npts,nvm)                             :: veget_max_old       !! "Maximal" coverage fraction of a PFT  
                                                                                       !! (LAI-> infinity) on ground (unitless) 
    REAL(r_std), DIMENSION(npts,nvm)                            :: mortality           !! Fraction of individual dying this time 
                                                                                       !! step (0 to 1, unitless) 
    REAL(r_std), DIMENSION(npts)                                :: vartmp              !! Temporary variable used to add history
    REAL(r_std), DIMENSION(npts,nvm)                            :: histvar             !! History variables

    REAL(r_std)                  :: za, zb, ta, tb, slope, offset, biomass_loss        !! Arsene 20-08-2014 / 03-04-2015 local variables to obtain biomass_loss
    REAL(r_std), DIMENSION(npts,3)                              :: veget_layer         !! Arsene 20-08-2014  veget_max by layer : grasses + shrubs + trees

!JCADD lcchange of managed grassland
    ! "maximal" coverage fraction of a PFT (LAI -> infinity) on ground
    INTEGER(i_std)                       :: ier
    LOGICAL                               :: l_error =.FALSE.

!ENDJCADD
!_ ================================================================================================================================

    IF (bavard.GE.3) WRITE(numout,*) 'Entering stomate_lpj'

  
  !! 1. Initializations
!JCADD

    IF (firstcall) THEN

        firstcall = .FALSE.

        !Config  Key  = ENABLE_GRAZING
        !Config  Desc = grazing allowed
        !Config  Def  = n
        !Config  Help = flag for choose if you want animals or not.
        !
        enable_grazing = .FALSE.
        CALL getin_p('ENABLE_GRAZING',enable_grazing)
        WRITE (numout,*) 'enable_grazing',enable_grazing
        WRITE (numout,*) 'manage',is_grassland_manag
        WRITE (numout,*) 'cut',is_grassland_cut
        WRITE (numout,*) 'grazed',is_grassland_grazed
        ALLOCATE (flag_cutting        (npts,nvm), stat=ier) ; l_error=l_error .OR. (ier .NE. 0)
        ALLOCATE (when_growthinit_cut (npts,nvm), stat=ier) ; l_error=l_error .OR. (ier .NE. 0)
        IF (l_error) THEN
            STOP 'error allocation memory for flag_cutting in stomateLPJ'
        END IF

        flag_cutting(:,:) = 0
        when_growthinit_cut(:,:) = 20.0
    END IF
!ENDJCADD
    
    !! 1.1 Initialize variables to zero
    co2_to_bm(:,:) = zero
    co2_fire(:,:) = zero
    npp_daily(:,:) = zero
    resp_maint(:,:) = zero
    resp_growth(:,:) = zero
    harvest_above(:) = zero
    bm_to_litter(:,:,:,:) = zero
    cn_ind(:,:) = zero
    woodmass_ind(:,:) = zero
    turnover_daily(:,:,:,:) = zero
    
    !! 1.2  Initialize variables to veget_max
    veget_max_old(:,:) = veget_max(:,:)

    !! 1.3 Calculate some vegetation characteristics

    !! 1.3.1 Calculate some vegetation characteristics 
    !        Calculate cn_ind (individual crown mass) and individual height from
    !        state variables if running DGVM or dynamic mortality in static cover mode
    !??        Explain (maybe in the header once) why you mulitply with veget_max in the DGVM
    !??        and why you don't multiply with veget_max in stomate.

    IF ( control%ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN
       IF(control%ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_max(:,:))/ind(:,:)
          ENDWHERE
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon)+biomass(:,:,isapbelow,icarbon) &
                  +biomass(:,:,iheartabove,icarbon)+biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF

       !WRITE(numout,*) 'zd bio1 crown biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)
       CALL crown (npts,  PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_max, cn_ind, height, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zd bio2 crown biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)
    ENDIF

    !! 1.3.2 Prescribe characteristics if the vegetation is not dynamic
    !        IF the DGVM is not activated, the density of individuals and their crown
    !        areas don't matter, but they should be defined for the case we switch on
    !        the DGVM afterwards. At the first call, if the DGVM is not activated, 
    !        impose a minimum biomass for prescribed PFTs and declare them present.
    !WRITE(numout,*) 'zd leaffrac1 prescribe leaf_frac(1,10,:)', leaf_frac(1,10,:)
    CALL prescribe (npts, &
         veget_max, dt_days, PFTpresent, everywhere, when_growthinit, &
         biomass, leaf_frac, ind, cn_ind, co2_to_bm, height, dia_cut)   !! Arsene 04-09-2014 - Add height !! 27-08-2015 Arsene - Add dia_cut
    !WRITE(numout,*) 'zd leaffrac2 prescribe leaf_frac(1,10,:)', leaf_frac(1,10,:)
    !WRITE(numout,*) 'zd bio3 prescribe biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

  !! 2. Climatic constraints for PFT presence and regenerativeness

    !   Call this even when DGVM is not activated so that "adapted" and "regenerate"
    !   are kept up to date for the moment when the DGVM is activated.
    CALL constraints (npts, dt_days, &
         t2m_month, t2m_min_daily,when_growthinit, &
         adapted, regenerate, Tseason)  !!, snowdz_min)             !! Arsene 19-08-2014 Add snowdz_min

!! 2.1 Protection of shrub by snow. BY Arsene   08-2014
!! #############################     AVANT été placé ici la protection des SHRUBS     ###########################


  !! 3. Determine introduction and elimination of PTS based on climate criteria
 
    IF ( control%ok_dgvm ) THEN
      
       !! 3.1 Calculate introduction and elimination
       CALL pftinout (npts, dt_days, adapted, regenerate, &
            neighbours, veget_max, &
            biomass, ind, cn_ind, age, leaf_frac, npp_longterm, lm_lastyearmax, senescence, &
            PFTpresent, everywhere, when_growthinit, need_adjacent, RIP_time, &
            co2_to_bm, &
            avail_tree, avail_grass, &
!JCADD
            sla_calc)
!ENDJCADD 

       !WRITE(numout,*) 'zd leaffrac3 pftinout leaf_frac(1,10,:)', leaf_frac(1,10,:)
       !WRITE(numout,*) 'zd bio4 pftinout biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)
       !! 3.2 Reset attributes for eliminated PFTs.
       !     This also kills PFTs that had 0 leafmass during the last year. The message
       !     "... after pftinout" is misleading in this case.
       !WRITE(numout,*) 'zdcheck1 kill leaf_age(1,10,:)', leaf_age(1,10,:)
       CALL kill (npts, 'pftinout  ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zdcheck2 kill leaf_age(1,10,:)', leaf_age(1,10,:)
       !WRITE(numout,*) 'zd leaffrac4 kill leaf_frac(1,10,:)', leaf_frac(1,10,:)
       !WRITE(numout,*) 'zd bio5 kill biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

       
       !! 3.3 Calculate new crown area and diameter 
       !      Calculate new crown area, diameter and maximum vegetation cover**[No longer used in the subroutine]
       !      unsure whether this is really required
       !      - in theory this could ONLY be done at the END of stomate_lpj
       !      calculate woodmass of individual tree
       WHERE ((ind(:,:).GT.min_stomate))
          WHERE  ( veget_max(:,:) .GT. min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))*veget_max(:,:))/ind(:,:)
          ELSEWHERE
             woodmass_ind(:,:) =(biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE

       ENDWHERE
       
       ! Calculate crown area and diameter for all PFTs (including the newly established)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_max, cn_ind, height, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zd bio6 crown biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

    ENDIF

    
  !! 4. Phenology

    !! 4.1 Write values to history file
    !      Current values for ::when_growthinit 
    CALL xios_orchidee_send_field("WHEN_GROWTHINIT",when_growthinit)

    CALL histwrite_p (hist_id_stomate, 'WHEN_GROWTHINIT', itime, when_growthinit, npts*nvm, horipft_index)

    ! Set and write values for ::PFTpresent
    WHERE(PFTpresent)
       histvar=un
    ELSEWHERE
       histvar=zero
    ENDWHERE

    CALL xios_orchidee_send_field("PFTPRESENT",histvar)

    CALL histwrite_p (hist_id_stomate, 'PFTPRESENT', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_midwinter
    WHERE(gdd_midwinter.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_midwinter
    ENDWHERE

    CALL xios_orchidee_send_field("GDD_MIDWINTER",histvar)

    CALL histwrite_p (hist_id_stomate, 'GDD_MIDWINTER', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for gdd_m5_dormance
    WHERE(gdd_m5_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=gdd_m5_dormance
    ENDWHERE
    
    CALL xios_orchidee_send_field('GDD_M5_DORMANCE',histvar)
    CALL histwrite_p (hist_id_stomate, 'GDD_M5_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    ! Set and write values for ncd_dormance
    WHERE(ncd_dormance.EQ.undef)
       histvar=val_exp
    ELSEWHERE
       histvar=ncd_dormance
    ENDWHERE

    CALL xios_orchidee_send_field("NCD_DORMANCE",histvar)

    CALL histwrite_p (hist_id_stomate, 'NCD_DORMANCE', itime, histvar, npts*nvm, horipft_index)

    !! 4.2 Calculate phenology
    CALL phenology (npts, dt_days, PFTpresent, &
         veget_max, &
         tlong_ref, t2m_month, t2m_week, gpp_daily, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_month, moiavail_week, &
         gdd_m5_dormance, gdd_midwinter, ncd_dormance, ngd_minus5, &
         senescence, time_hum_min, &
         biomass, leaf_frac, leaf_age, &
         when_growthinit, co2_to_bm, begin_leaves, &!)
!JCADD
         sla_calc)
!ENDJCADD
    !WRITE(numout,*) 'zdcheck3 phenology leaf_age(1,10,:)', leaf_age(1,10,:)
    !WRITE(numout,*) 'zd leaffrac5 phenology leaf_frac(1,10,:)', leaf_frac(1,10,:)
    !WRITE(numout,*) 'zd bio7 phenology biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)
    
  !! 5. Allocate C to different plant parts
    
    CALL alloc (npts, dt_days, &
         lai, veget_max, senescence, when_growthinit, &
         moiavail_week, tsoil_month, soilhum_month, &
         biomass, age, leaf_age, leaf_frac, rprof, f_alloc, &!)
!JCADD
         sla_calc, when_growthinit_cut)
!ENDJCADD
    !WRITE(numout,*) 'zdcheck4 alloc leaf_age(1,10,:)', leaf_age(1,10,:)
    !WRITE(numout,*) 'zd leaffrac6 alloc leaf_frac(1,10,:)', leaf_frac(1,10,:)
    !WRITE(numout,*) 'zd bio8 alloc biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

  !! 6. NPP, maintenance and growth respiration

    !! 6.1 Calculate NPP and respiration terms
    CALL npp_calc (npts, dt_days, &
         PFTpresent, &
         tlong_ref, t2m_daily, tsoil_daily, lai, rprof, &
         gpp_daily, f_alloc, bm_alloc, resp_maint_part,&
         biomass, leaf_age, leaf_frac, age, &
         resp_maint, resp_growth, npp_daily, &!)
!JCADD
         sla_calc, sla_age1,N_limfert)
!ENDJCADD
    !WRITE(numout,*) 'zdcheck5 npp_calc leaf_age(1,10,:)', leaf_age(1,10,:)
    !WRITE(numout,*) 'zd leaffrac7 npp_calc leaf_frac(1,10,:)', leaf_frac(1,10,:)
    !WRITE(numout,*) 'zd bio9 npp_calc biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

    !! 6.2 Kill slow growing PFTs in DGVM or STOMATE with constant mortality
    IF ( control%ok_dgvm .OR. .NOT.lpj_gap_const_mort) THEN
       CALL kill (npts, 'npp       ', lm_lastyearmax,  &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zdcheck6 kill leaf_age(1,10,:)', leaf_age(1,10,:)
       !WRITE(numout,*) 'zd leaffrac8 kill leaf_frac(1,10,:)', leaf_frac(1,10,:)
       !WRITE(numout,*) 'zd bio10 kill biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

       !! 6.2.1 Update wood biomass      
       !        For the DGVM
       IF(control%ok_dgvm) THEN
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  ((biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon)) & 
                  *veget_max(:,:))/ind(:,:)
          ENDWHERE

       ! For all pixels with individuals
       ELSE
          WHERE (ind(:,:).GT.min_stomate)
             woodmass_ind(:,:) = &
                  (biomass(:,:,isapabove,icarbon) + biomass(:,:,isapbelow,icarbon) &
                  + biomass(:,:,iheartabove,icarbon) + biomass(:,:,iheartbelow,icarbon))/ind(:,:)
          ENDWHERE
       ENDIF ! control%ok_dgvm

       !! 6.2.2 New crown area and maximum vegetation cover after growth
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind,&
            veget_max, cn_ind, height, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zd bio11 crown biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

    ENDIF ! control%ok_dgvm
    
  !! 7. fire

    !! 7.1. Burn PFTs
    CALL fire (npts, dt_days, litterpart, &
         litterhum_daily, t2m_daily, lignin_struc, veget_max, &
         fireindex, firelitter, biomass, ind, &
         litter, dead_leaves, bm_to_litter, black_carbon, &
         co2_fire, MatrixA)
    !WRITE(numout,*) 'zd bio12 fire biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

!JCADD update available and not available litter for grazing litter
! after fire burning
  litter_avail(:,:,:) = litter(:,:,:,iabove,icarbon) * &
            litter_avail_frac(:,:,:)
  litter_not_avail(:,:,:) = litter(:,:,:,iabove,icarbon) * &
            (1.0 - litter_avail_frac(:,:,:))
!ENDJCADD

    !! 7.2 Kill PFTs in DGVM
    IF ( control%ok_dgvm ) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'fire      ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zdcheck7 kill leaf_age(1,10,:)', leaf_age(1,10,:)
       !WRITE(numout,*) 'zd leaffrac9 kill leaf_frac(1,10,:)', leaf_frac(1,10,:)
       !WRITE(numout,*) 'zd bio13 kill biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

    ENDIF ! control%ok_dgvm
 
  !! 8. Tree mortality
    ! Does not depend on age, therefore does not change crown area.
    CALL gap (npts, dt_days, &
         npp_longterm, turnover_longterm, lm_lastyearmax, &
         PFTpresent, biomass, ind, bm_to_litter, mortality, t2m_min_daily, Tmin_spring, Tmin_spring_time, &!)
!JCADD
         sla_calc, snowdz_min, snowtemp_min) !! Arsene 31-03-2015 Add snowdz_min & snowtemp_min
!ENDJCADD
    !WRITE(numout,*) 'zd bio14 gap biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)


  !! 8.1 Protection of shrub by snow - shrub mortality. BY Arsene   08-2014

!! Arsene 20-08-2014 Add protection of shrub by snow. START.
!!
!! #############################     protection of shrub by snow     ###########################
!! ## Ne prend pas en compte : -  Take care about "termal chock" : si T chute trop rapidement
!! ## Ne prend pas en compte : -  IF All(T) < T_crit THEN no loss of biomasse but mortality... 
    IF ( ANY(is_shrub(:)) ) THEN

        !! ## On calcul la fration de strate de végétation, afin de "simuler" des différence de dépôt de masse / pixel
        veget_layer(:,:)=zero
        DO ji = 1, npts

            DO j = 1, nvm
                IF ( is_tree(j) ) THEN
                    veget_layer(ji,3)=veget_layer(ji,3)+veget_max(ji,j)     !! Trees
                ELSEIF ( is_shrub(j) ) THEN
                    veget_layer(ji,2)=veget_layer(ji,2)+veget_max(ji,j)     !! Shrubs
                ELSE
                    veget_layer(ji,1)=veget_layer(ji,1)+veget_max(ji,j)     !! Grasse & bare soil
                ENDIF
            ENDDO
        ENDDO

        DO ji = 1, npts     !! Arsene 31-03-2015  On pourrait mettre à ce niveau le IF SUM(snowdz_min(ji,:)) .GT. min_stomate
           DO j = 2,nvm
              !! ## Si besoin : s'il y a de la neige, des buissons, de la vegetation,... On commence
              IF ( is_shrub(j) .AND. ( SUM(snowdz_min(ji,:)) .GT. min_stomate ).AND. &
                      ind(ji,j).GE.min_stomate .AND.  PFTpresent(ji,j) .AND. &
                      ((.NOT.control%ok_dgvm.AND.(height(ji,j).GT.(height_presc(j)/fact_min_height))) .OR. control%ok_dgvm)) THEN !! Arsene 03-04-2015 If no DGVM, need heigth >= min_height. Note: Pour le moment, height_max=height_presc & height_min = height_presc/10 (see stomate_lpj)
                      !! Arsene 27-08-2015 - Could be interesting to ad a crteria : if ANY(temp) < temp_crit. (but wrong in extrem case)
                 !! ## Calcul des courbes de température, de bas en haut. Afin de calculer tmin_crit - t
                 biomass_loss = zero
                 DO i = nsnow+1, 1, -1
                    !! ## Calcul de la temperature - Equation linéaires à 2 inconnues - on cherche les solutions de T(z)=tmin_crit
                    !! ## Afin de garder une logique "naturelle" (et autre), on résoud de bas en haut
                    !! ## Il y a 4 courbes avec :   - Comme points charnières sont dz= 0 ; puis au milieu de chaque couche ;
                    !! ##                                  et en surface (snowdz_min(1+2+3)) -> za et zb
                    !! ##                           - Les températures associées sont tsoil_daily(1) ; snowtemp_min(a) (avec a=3,2,1) ;
                    !! ##                                  T2m_min_daily -> ta et tb
                    !! For za & tb (compatible si nsow > 3)
                    IF ( i .EQ. nsnow+1 ) THEN
                       za = zero ; ta = tsoil_daily(ji,1)    !! Ce n'est la temperature minimum... Mais assez stable ==> à voir
                    ELSEIF ( i .EQ. nsnow ) THEN
                       za = snowdz_min(ji,i)/2 ; ta = snowtemp_min(ji,i)
                    ELSE
                       za =  ( SUM(snowdz_min(ji,i:nsnow)) + SUM(snowdz_min(ji,i+1:nsnow)) )/2 ; ta = snowtemp_min(ji,i)
                    ENDIF ! za & ta
                    !! For zb & tb (compatible si nsow > 3)
                    IF ( i .EQ. (nsnow+1) ) THEN
                       zb = snowdz_min(ji,i-1)/2 ; tb = snowtemp_min(ji,i-1)
                    ELSEIF ( i .EQ. 1) THEN
                       zb = SUM(snowdz_min(ji,i:nsnow)) ; tb = t2m_min_daily(ji)
                    ELSE
                       zb = ( SUM(snowdz_min(ji,i-1:nsnow)) + SUM(snowdz_min(ji,i:nsnow)) )/2 ;  tb = snowtemp_min(ji,i-1)
                    ENDIF !! zb & tb

                    !! ## On simule l'effet d'accumulation préférentiel de la neige sur les buissons via différence de dépôt de masse / pixel
                    !! ##    --> un facteur à z est ajouter afin de simuler le transfert de masse de neige au niveau des buissons (entre 0 et 0.5)
                    !! Arsene 28-01-2018 Change limit 0.5 by 0.2 and the equation in order to take into account reviewer comment - START
                    IF ( za.GT.min_stomate .AND. (veget_layer(ji,2)+veget_layer(ji,3)).GT.min_stomate & 
                                         & .AND. (veget_layer(ji,2)+veget_layer(ji,3)).LE.0.2 ) THEN
                        za = za * (1 + 4*(veget_layer(ji,2) + veget_layer(ji,3)))
                    ELSEIF ( za.GT.min_stomate .AND. (veget_layer(ji,2)+veget_layer(ji,3)).GT.0.2 & 
                                         & .AND. (veget_layer(ji,2)+veget_layer(ji,3)).LT.1 ) THEN
                        za = za * (2 - veget_layer(ji,2) - veget_layer(ji,3))
                    ENDIF

                    IF ( zb.GT.min_stomate .AND. (veget_layer(ji,2)+veget_layer(ji,3)).GT.min_stomate &
                                         & .AND. (veget_layer(ji,2)+veget_layer(ji,3)).LE.0.2 ) THEN
                        zb = zb * (1 + 4*(veget_layer(ji,2) + veget_layer(ji,3)))
                    ELSEIF ( zb.GT.min_stomate .AND. (veget_layer(ji,2)+veget_layer(ji,3)).GT.0.2 & 
                                         & .AND. (veget_layer(ji,2)+veget_layer(ji,3)).LT.1 ) THEN
                        zb = zb * (2 - veget_layer(ji,2) - veget_layer(ji,3))
                    ENDIF
                    !! Arsene 28-01-2018 Change limit 0.5 by 0.2 and the equation in order to take into account reviewer comment - END

                     !! ## Si On a des températures inférieur à tmin_crit ET que l'on se situe à z < height (hauteur du buisson)
                     IF ( ( (ta .LT. tmin_crit(j)) .OR. (tb .LT. tmin_crit(j) ) ) & 
                           & .AND. ((za .LT. height(ji,j)) .OR. (zb .LT. height(ji,j)) ) &  !! Arsene 01-04-2015 probably, don't need for zb vecause za<zb
                           & .AND. (control%ok_dgvm .OR. &
                           & (za.GT.height(ji,j)/fact_min_height .OR. zb.GT.height(ji,j)/fact_min_height))) THEN     !! Arsene 01-04-2015 probably, don't need for za vecause za<zb

                        !! ## Pour simplifier le calcul, on fait un changement de repère
                        !! ##      - On passe de z à z/height: on a donc des valeurs [0-1] et vb-va = % hauteur
                        !! ##      - On passe de t à "tmin_crit - t": On a l'écart de t à la température critique
                        !! ##             (Fonctionne si tmincrit < 0 (sinon prendre la valeur absolu de la diff)
                        
                        za = za/height(ji,j) ; zb = zb/height(ji,j)
                        ta = tmin_crit(j)-ta ; tb = tmin_crit(j)-tb
                         
                        !! ## On calcul l'équation de chacune des courbes
                        slope = ( tb - ta ) / (zb - za )
                        offset = ( zb * ta - tb * za ) / ( zb - za )

                        !! ## On se s'intéresse qu'aux parties où t > 0
                        !! ## ==> Si ta et tb <0, on regarde où t = 0
                        IF ( (ta .LT. zero ) .AND. (tb .GT. zero ) ) THEN
                            za = - offset / slope ; ta = zero
                        ELSEIF ( (ta .GT. zero ) .AND. (tb .LT. zero) ) THEN
                            zb = - offset / slope ; tb = zero
                        ENDIF
                        
                        !! ## Si le buisson est plus petit que SUM(snowdz) il faut vérifier que z<=height
                        IF ( zb .GT. height(ji,j) ) THEN 
                            zb = height(ji,j) ; tb = slope * zb + offset
                        ELSEIF ( za .GT. height(ji,j) ) THEN                !! Arsene 01-04-2015 probably, don't need for zb vecause za<zb
                            za = height(ji,j) ; ta = slope * za + offset    !! Arsene 01-04-2015 probably, don't need for zb vecause za<zb
                        ENDIF
                        !! ## Si hors DGVM ==> on ne touche pas au buisson si height < min_height
                        IF ( .NOT.control%ok_dgvm .AND. zb .LT. height(ji,j)/fact_min_height ) THEN
                            zb = height(ji,j)/fact_min_height ; tb = slope * zb + offset
                        ELSEIF ( .NOT.control%ok_dgvm .AND. za .LT. height(ji,j)/fact_min_height ) THEN                !! Arsene 01-04-2015 probably, don't need for za vecause za<zb
                            za = height(ji,j)/fact_min_height ; ta = slope * za + offset    !! Arsene 01-04-2015 probably, don't need for za vecause za<zb
                        ENDIF

                        !! ## Perte de biomasse fonction de l'aire sous la courbe
                        !! ## correspont à 4% par degrès de différence sur chaque hauteur
                        biomass_loss = 0.04 * ( (slope/2 * (zb**2-za**2)) + (offset*(zb-za)) ) + biomass_loss
                        !! Arsene 07-04-2014 Take care of snow temperature at first layer: with wind, it could be << T_air

                    ENDIF
                 ENDDO ! nsnow: for nsnow layers
                 biomass_loss = MIN(un,biomass_loss)

                 !! ## Maintenant, si biomasse loss > 0, alors on doit "enlever" de la bionass et de la hauteur (dia & nb ind inchangé)
                 !! ##  ==> on diminue la hauteur, à travers une diminution de biomasse.
                 IF ( biomass_loss.GT.min_stomate ) THEN

                     !! ## On calcul le diametre avant la coupe qui sera alors sauvegardé (afinde maintenir la "rupture" d'allometry) !! Arsene 27-08-2015
                     !! ## Vérifier si jamais que la hauteur est bonne (non modif prescedement...) - On part ici du principe que l'on travail sur la hauteur "du jour"...
                      IF (shrubs_like_trees) THEN
                         dia_cut(ji,j) = (height(ji,j)/pipe_tune2_for_shrub)**(1/pipe_tune3_for_shrub)
                      ELSE !! shrub and New Allometry
                         dia_cut(ji,j) = (height(ji,j)*height_presc(j) / (pipe_tune_shrub2*(height_presc(j)-height(ji,j))) ) &
                                 & **(1/pipe_tune_shrub3) /100
                      ENDIF
                      IF ( dia_cut(ji,j) .GE. maxdia(j)) dia_cut(ji,j)=zero !! Arsene 21-09-15 - If biomass accumulation... 


                     !! ## Calcul de la nouvelle heuteur du buissson (height)
                     !! ## Si DGVM non activé, heigth peut pas être inférieur à min heigh fixé ! Hmin = 50 cm ????????
                     !! ## Deux solution:  -soit on dis simplement que si ça réduit trop la taille, alors on fixe heigh = heigh_min
                     !! ##                      pb : cela indique que d'un seul coup, la partie inférieur de la plante deviens insensible au froid
                     !! ##                           Courbe non continue
                     !! ##                      ==> Choisi car même sensibilité
                     !! ##                 -soit on dis que la partie < height_min n'est pas sensible dès le début (za & zb > heigh_min)
                     !! ##                      pb : sensibilité différente avec ou sans DGVM

                     IF (.NOT.control%ok_dgvm .AND. &
                             & ( ((1-biomass_loss)*height(ji,j)) .LT. (height_presc(j)/fact_min_height)))   THEN !! Arsene 02-04-2015
                          biomass_loss = 1 - (height_presc(j)/fact_min_height) / height(ji,j)     !! Arsene 02-04-2015
                          height(ji,j) = height_presc(j)/fact_min_height                          !! Arsene 02-04-2015
                     ELSE
                          height(ji,j) = (1-biomass_loss)*height(ji,j)
                     ENDIF


                     !! Arsene 23-08-2016 ADD - start
                     !! ## If we delete a very important part of shrubs, kill the PFT via ind
                     IF ( control%ok_dgvm .AND. biomass_loss.gt.0.95 ) THEN
                         ind(ji,j) = zero
                         biomass_loss=1.
                     ENDIF
                     !! Arsene 23-08-2016 ADD - end

                     !! ## Simultanément, après avoir diminuer la hauteur on doit diminuer la biomasse (dans un cylindre, si uniquement height diminue ==> propotionnel à la biomasse
                     !! ## 2 méthodes : 1) on diminue seulement la biomasse en surface (ileaf, isapabove, iheartabove, ifruit et ?icarbres?): plus réaliste 
                     !! ##                  ==> pb: c'est woodmass_ind qui détermine le diamétre et la hauteur. Si en surface biomasse = 0, on a malgrès tout woodmass_ind non nul donc dia et heigth >0...
                     !! ##              2) on retire la avec une même proportion toute les biomasses ==> Sécurité mais moins réel.
                     !! ##                  ==> méthode choisi. Modifier l'équation des shrubs ?

                     !! ## On supprime les biomasses correspondantes, et on les rajoute dans la litière
                     DO m = 1,nelements
                           bm_to_litter(ji,j,ileaf,m) = bm_to_litter(ji,j,ileaf,m) + biomass_loss*biomass(ji,j,ileaf,m)
                           bm_to_litter(ji,j,isapabove,m) = bm_to_litter(ji,j,isapabove,m) + biomass_loss*biomass(ji,j,isapabove,m)
                           bm_to_litter(ji,j,isapbelow,m) = bm_to_litter(ji,j,isapbelow,m) + biomass_loss*biomass(ji,j,isapbelow,m)
                           bm_to_litter(ji,j,iheartabove,m) = bm_to_litter(ji,j,iheartabove,m) + &
                                                   biomass_loss*biomass(ji,j,iheartabove,m)
                           bm_to_litter(ji,j,iheartbelow,m) = bm_to_litter(ji,j,iheartbelow,m) + &
                                                   biomass_loss*biomass(ji,j,iheartbelow,m)
                           bm_to_litter(ji,j,iroot,m) = bm_to_litter(ji,j,iroot,m) + biomass_loss*biomass(ji,j,iroot,m)
                           bm_to_litter(ji,j,ifruit,m) = bm_to_litter(ji,j,ifruit,m) + biomass_loss*biomass(ji,j,ifruit,m)
                           bm_to_litter(ji,j,icarbres,m) = bm_to_litter(ji,j,icarbres,m) + biomass_loss*biomass(ji,j,icarbres,m)

                           biomass(ji,j,ileaf,m) = (1-biomass_loss) * biomass(ji,j,ileaf,m)
                           biomass(ji,j,isapabove,m) = (1-biomass_loss) * biomass(ji,j,isapabove,m)
                           biomass(ji,j,isapbelow,m) = (1-biomass_loss) * biomass(ji,j,isapbelow,m)
                                     biomass(ji,j,iheartabove,m) = (1-biomass_loss) * biomass(ji,j,iheartabove,m)
                           biomass(ji,j,iheartbelow,m) = (1-biomass_loss) * biomass(ji,j,iheartbelow,m)
                           biomass(ji,j,iroot,m) = (1-biomass_loss) * biomass(ji,j,iroot,m)
                           biomass(ji,j,ifruit,m) = (1-biomass_loss) * biomass(ji,j,ifruit,m)
                           biomass(ji,j,icarbres,m) = (1-biomass_loss) * biomass(ji,j,icarbres,m)
                     ENDDO

                     !! ## On met à jour la woodmass (au cas où ça serve rapidement... Mais possible que NON)
                     IF ( veget_max(ji,j) .GT. min_stomate) THEN
                           woodmass_ind(ji,j) = &
                             ((biomass(ji,j,isapabove,icarbon) + biomass(ji,j,isapbelow,icarbon) &
                             + biomass(ji,j,iheartabove,icarbon) + biomass(ji,j,iheartbelow,icarbon))*veget_max(ji,j))/ind(ji,j)
                     ELSE
                           woodmass_ind(ji,j) =(biomass(ji,j,isapabove,icarbon) + biomass(ji,j,isapbelow,icarbon) &
                             + biomass(ji,j,iheartabove,icarbon) + biomass(ji,j,iheartbelow,icarbon))/ind(ji,j)
                     ENDIF

                     !! ## si plus assez de biomasse, on "KILL"
                     IF ( SUM(biomass(ji,j,:,1)) .LT. min_stomate .AND. control%ok_dgvm) THEN
                          ind(ji,j) = zero
                     ENDIF

                ENDIF ! IF biomass_loss > 0 decrease biomass

              ENDIF ! IF shrub, if snow,... ==> when we use that !
         ENDDO ! on pft
      ENDDO ! on grill
   ENDIF ! if any shrub
!! Arsene 20-08-2014 Add protection of shrub by snow. END.


    IF ( control%ok_dgvm ) THEN

       ! reset attributes for eliminated PFTs
       CALL kill (npts, 'gap       ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zdcheck8 kill leaf_age(1,10,:)', leaf_age(1,10,:)
       !WRITE(numout,*) 'zd leaffrac10 kill leaf_frac(1,10,:)', leaf_frac(1,10,:)
       !WRITE(numout,*) 'zd bio15 kill biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

    ENDIF

  !! 9. Calculate vcmax 

    CALL vmax (npts, dt_days, &
         leaf_age, leaf_frac, &
         vcmax, moiavail_month, &      !! Arsene 24-06-2014 for dessication add soilhum_month 
!JCADD
         N_limfert)
!ENDJCADD
    !WRITE(numout,*) 'zdcheck9 vmax leaf_age(1,10,:)', leaf_age(1,10,:)
    !WRITE(numout,*) 'zd leaffrac11 vmax leaf_frac(1,10,:)', leaf_frac(1,10,:)

  !! 10. Leaf senescence, new lai and other turnover processes

    CALL turn (npts, dt_days, PFTpresent, &
         herbivores, &
         maxmoiavail_lastyear, minmoiavail_lastyear, &
         moiavail_week,  moiavail_month,tlong_ref, t2m_month, t2m_week, veget_max, &
         gdd_from_growthinit, leaf_age, leaf_frac, age, lai, biomass, &
         turnover_daily, senescence,turnover_time, npp0_cumul, &  !! Arsene 25-06-2014 NPPcumul : add of  npp0_cumul
!JCADD
         sla_calc)
!ENDJCADD
    !WRITE(numout,*) 'zdcheck10 turn leaf_age(1,10,:)', leaf_age(1,10,:)
    !WRITE(numout,*) 'zd leaffrac12 turn leaf_frac(1,10,:)', leaf_frac(1,10,:)
    !WRITE(numout,*) 'zd bio16 turn biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)


    !! 11. Light competition

    !! If not using constant mortality then kill with light competition
!    IF ( control%ok_dgvm .OR. .NOT.(lpj_gap_const_mort) ) THEN
    IF ( control%ok_dgvm ) THEN
 
       !! 11.1 Light competition
       CALL light (npts, dt_days, &
            veget_max, fpc_max, PFTpresent, cn_ind, lai, maxfpc_lastyear, &
            lm_lastyearmax, ind, biomass, veget_lastlight, bm_to_litter, mortality, &!)
!JCADD
         sla_calc)
!ENDJCADD
       !WRITE(numout,*) 'zd bio17 light biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

!! 23-08-2016 Arsene ADD - START
       !! 11.1.bis Give elimited attribute if nonsence PFT

       ! Update height and cn_ind
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_max, cn_ind, height, dia_cut) !! Arsene 27-08-2015 - Add dia_cut

       DO j = 2,nvm ! loop over PFTs

          ! Test the value of new theorical (lpj_cover)cf  veget_frac and chech if the difference is too important
!*!          IF ( natural(j) ) THEN    !! Cette partie est potentiellement a reconsidere, car cn_ind evolue encore apres... Mettre juste avant CALL gap ? (ralonge le temps de calcul)
!*!             WHERE ( veget_max_old(:,j)/(ind(:,j)*cn_ind(:,j)) .GT. 1./(min_stomate*10000.) ) !! il faut passer par un GT... comment le définir ?
!*!                ind(j,:)=zero
!*!             ENDWHERE
!*!          ENDIF ==> MARCHE PAS. Probablement à cause d'estabish ou kill... ou "ind" qui est trop faible Vérif que inf est pas trop faible avant 

          ! Test the value of height is compatible with survive  
          IF (  is_tree(j) .OR. is_shrub(j) ) THEN  ! ELSE KILL NVPS
             WHERE ( height(:,j) .LT. min_stomate*10000 )
                 ind(:,j)=zero
             ENDWHERE
          ENDIF
       ENDDO
!! 23-08-2016 Arsene ADD - END

       !! 11.2 Reset attributes for eliminated PFTs
       CALL kill (npts, 'light     ', lm_lastyearmax, &
            ind, PFTpresent, cn_ind, biomass, senescence, RIP_time, &
            lai, age, leaf_age, leaf_frac, npp_longterm, &
            when_growthinit, everywhere, veget_max, bm_to_litter, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zdcheck11 kill leaf_age(1,10,:)', leaf_age(1,10,:)
       !WRITE(numout,*) 'zd leaffrac13 kill leaf_frac(1,10,:)', leaf_frac(1,10,:)
       !WRITE(numout,*) 'zd bio18 kill biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

    ENDIF

    
  !! 12. Establishment of saplings
    
    IF ( control%ok_dgvm .OR. .NOT.lpj_gap_const_mort ) THEN

       !! 12.1 Establish new plants
       CALL establish (npts, dt_days, PFTpresent, regenerate, &
            neighbours, resolution, need_adjacent, herbivores, &
            precip_lastyear, gdd0_lastyear, lm_lastyearmax, &
            cn_ind, lai, avail_tree, avail_grass, npp_longterm, &
            leaf_age, leaf_frac, &
            ind, biomass, age, everywhere, co2_to_bm, veget_max, woodmass_ind, &!)
!JCADD
            sla_calc,height,dia_cut) !! Arsene 26-08-2015 - Add height (in)
!ENDJCADD

       !WRITE(numout,*) 'zdcheck12 establish leaf_age(1,10,:)', leaf_age(1,10,:)
       !WRITE(numout,*) 'zd leaffrac14 establish leaf_frac(1,10,:)', leaf_frac(1,10,:)
       !WRITE(numout,*) 'zd bio19 establish biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

       !! 12.2 Calculate new crown area (and maximum vegetation cover)
       CALL crown (npts, PFTpresent, &
            ind, biomass, woodmass_ind, &
            veget_max, cn_ind, height, dia_cut) !! Arsene 27-08-2015 - Add dia_cut
       !WRITE(numout,*) 'zd bio20 crown biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)
    ENDIF
!JCADD Grassland_management


    !
    ! 13 calculate grazing by animals or cutting for forage
    !
    IF (enable_grazing) THEN
      WRITE (numout, *) 'enter the grassland management process'
        CALL Main_Grassland_Management(&
           npts           , &
           dt_days        , &
           day_of_year    , &
           t2m_daily      , &
           t2m_min_daily  , &
           t2m_14         , &
           tsurf_daily    , &
           snow           , &
           biomass        , &
           bm_to_litter   , &
           litter         , &
           litter_avail   , &
           litter_not_avail , &
           .TRUE.         , &
           EndofYear      , &
!           ldrestart_read , &
!           ldrestart_write, &
!           index          , &
           flag_cutting   , &
           when_growthinit_cut , &
           lai,sla_calc,leaf_age,leaf_frac, &
           wshtotsum,sr_ugb,compt_ugb, &
           nb_ani,grazed_frac,import_yield,N_limfert)
    ENDIF
!ENDJCADD
  !! 13. Calculate final LAI and vegetation cover
    
    CALL cover (npts, cn_ind, ind, biomass, &
         veget_max, veget_max_old, lai, &
         litter, litter_avail, litter_not_avail, carbon, &
         turnover_daily, bm_to_litter,deepC_a, deepC_s,deepC_p)
    !WRITE(numout,*) 'zd bio21 cover biomass(1,10,ileaf,icarbon)',biomass(1,10,ileaf,icarbon)

  !! 14. Update litter pools to account for harvest
 
    ! the whole litter stuff:
    !    litter update, lignin content, PFT parts, litter decay, 
    !    litter heterotrophic respiration, dead leaf soil cover.
    !    No vertical discretisation in the soil for litter decay.\n
    ! added by shilong for harvest
    IF(harvest_agri) THEN
       CALL harvest(npts, dt_days, veget_max, &
            bm_to_litter, turnover_daily, &
            harvest_above)
    ENDIF

  !! 15. Land cover change

    !shilong adde turnover_daily
    IF(EndOfYear) THEN
       IF (lcchange) THEN
          CALL lcchange_main (npts, dt_days, veget_max, veget_max_new, &
               biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &
               co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind,flux10,flux100, &
!!$               prod10,prod100,prod10_total,prod100_total,&
!!$               convflux,cflux_prod_total,cflux_prod10,cflux_prod100,leaf_frac,&
               prod10,prod100,convflux,cflux_prod10,cflux_prod100,leaf_frac,&
               npp_longterm, lm_lastyearmax, litter, litter_avail, litter_not_avail, &
               carbon, &
               deepC_a, deepC_s, deepC_p)
       ENDIF
    ENDIF
    !MM déplacement pour initialisation correcte des grandeurs cumulées :
    cflux_prod_total(:) = convflux(:) + cflux_prod10(:) + cflux_prod100(:)
    prod10_total(:)=SUM(prod10,dim=2)
    prod100_total(:)=SUM(prod100,dim=2)
    
  !! 16. Total heterotrophic respiration

    tot_soil_carb(:,:) = zero
    tot_litter_carb(:,:) = zero
    DO j=2,nvm

       tot_litter_carb(:,j) = tot_litter_carb(:,j) + (litter(:,istructural,j,iabove,icarbon) + &
            &          litter(:,imetabolic,j,iabove,icarbon) + &
            &          litter(:,istructural,j,ibelow,icarbon) + litter(:,imetabolic,j,ibelow,icarbon))

       tot_soil_carb(:,j) = tot_soil_carb(:,j) + (carbon(:,iactive,j) + &
            &          carbon(:,islow,j)+  carbon(:,ipassive,j))

    ENDDO
    tot_litter_soil_carb(:,:) = tot_litter_carb(:,:) + tot_soil_carb(:,:)

!!$     DO k = 1, nelements ! Loop over # elements
!!$        tot_live_biomass(:,:,k) = biomass(:,:,ileaf,k) + biomass(:,:,isapabove,k) + biomass(:,:,isapbelow,k) +&
!!$             &                    biomass(:,:,iheartabove,k) + biomass(:,:,iheartbelow,k) + &
!!$             &                    biomass(:,:,iroot,k) + biomass(:,:,ifruit,k) + biomass(:,:,icarbres,k)
!!$    END DO ! Loop over # elements

    tot_live_biomass(:,:,:) = biomass(:,:,ileaf,:) + biomass(:,:,isapabove,:) + biomass(:,:,isapbelow,:) +&
             &                    biomass(:,:,iheartabove,:) + biomass(:,:,iheartbelow,:) + &
             &                    biomass(:,:,iroot,:) + biomass(:,:,ifruit,:) + biomass(:,:,icarbres,:)


    tot_turnover(:,:,:) = turnover_daily(:,:,ileaf,:) + turnover_daily(:,:,isapabove,:) + &
         &         turnover_daily(:,:,isapbelow,:) + turnover_daily(:,:,iheartabove,:) + &
         &         turnover_daily(:,:,iheartbelow,:) + turnover_daily(:,:,iroot,:) + &
         &         turnover_daily(:,:,ifruit,:) + turnover_daily(:,:,icarbres,:)

    tot_bm_to_litter(:,:,:) = bm_to_litter(:,:,ileaf,:) + bm_to_litter(:,:,isapabove,:) +&
         &             bm_to_litter(:,:,isapbelow,:) + bm_to_litter(:,:,iheartbelow,:) +&
         &             bm_to_litter(:,:,iheartabove,:) + bm_to_litter(:,:,iroot,:) + &
         &             bm_to_litter(:,:,ifruit,:) + bm_to_litter(:,:,icarbres,:)

    carb_mass_variation(:)=-carb_mass_total(:)
    carb_mass_total(:)=SUM((tot_live_biomass(:,:,icarbon)+tot_litter_carb+tot_soil_carb)*veget_max,dim=2) + &
         &                 (prod10_total + prod100_total)
    carb_mass_variation(:)=carb_mass_total(:)+carb_mass_variation(:)

    
  !! 17. Write history

    CALL xios_orchidee_send_field("RESOLUTION_X",resolution(:,1))
    CALL xios_orchidee_send_field("RESOLUTION_Y",resolution(:,2))
    CALL xios_orchidee_send_field("CONTFRAC_STOMATE",contfrac(:))
    CALL xios_orchidee_send_field("T2M_MONTH",t2m_month)
    CALL xios_orchidee_send_field("T2M_WEEK",t2m_week)
    CALL xios_orchidee_send_field("HET_RESP",resp_hetero(:,:))
!JCADD
!  WRITE(numout,*) 'write t2m_14'    ! Arsene 05-04-2015 Remove
    CALL xios_orchidee_send_field("T2M_14",t2m_14)
!    CALL xios_orchidee_send_field("LITTER_RESP",resp_hetero_litter_d(:,:))
!    CALL xios_orchidee_send_field("ACTIVE_RESP",resp_hetero_soil_d(:,iactive,:))
!    CALL xios_orchidee_send_field("SLOW_RESP",resp_hetero_soil_d(:,islow,:))
!    CALL xios_orchidee_send_field("PASSIVE_RESP",resp_hetero_soil_d(:,ipassive,:))
    CALL xios_orchidee_send_field("LITTER_STR_AVAIL",litter_avail(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_AVAIL",litter_avail(:,imetabolic,:))
    CALL xios_orchidee_send_field("LITTER_STR_NAVAIL",litter_not_avail(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_NAVAIL",litter_not_avail(:,imetabolic,:))
    CALL xios_orchidee_send_field("LITTER_STR_AVAILF",litter_avail_frac(:,istructural,:))
    CALL xios_orchidee_send_field("LITTER_MET_AVAILF",litter_avail_frac(:,imetabolic,:))
    CALL xios_orchidee_send_field("N_LIMFERT",N_limfert)
!ENDJCADD
    CALL xios_orchidee_send_field("CO2_FIRE",co2_fire)
    CALL xios_orchidee_send_field("CO2_TAKEN",co2_to_bm)
    CALL xios_orchidee_send_field("LAI",lai)
    CALL xios_orchidee_send_field("VEGET_MAX",veget_max)
    CALL xios_orchidee_send_field("NPP_STOMATE",npp_daily)
    CALL xios_orchidee_send_field("GPP",gpp_daily)
    CALL xios_orchidee_send_field("IND",ind)
    CALL xios_orchidee_send_field("CN_IND",cn_ind)
    CALL xios_orchidee_send_field("WOODMASS_IND",woodmass_ind)
    CALL xios_orchidee_send_field("TOTAL_M",tot_live_biomass)
    CALL xios_orchidee_send_field("MOISTRESS",moiavail_week)
    CALL xios_orchidee_send_field("LEAF_M",biomass(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("SAP_M_AB",biomass(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("SAP_M_BE",biomass(:,:,isapbelow,icarbon))
    CALL xios_orchidee_send_field("HEART_M_AB",biomass(:,:,iheartabove,icarbon))
    CALL xios_orchidee_send_field("HEART_M_BE",biomass(:,:,iheartbelow,icarbon))
    CALL xios_orchidee_send_field("ROOT_M",biomass(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_M",biomass(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("RESERVE_M",biomass(:,:,icarbres,icarbon))
    CALL xios_orchidee_send_field("TOTAL_TURN",tot_turnover)
    CALL xios_orchidee_send_field("LEAF_TURN",turnover_daily(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("MAINT_RESP",resp_maint)
    CALL xios_orchidee_send_field("GROWTH_RESP",resp_growth)
    CALL xios_orchidee_send_field("SAP_AB_TURN",turnover_daily(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("ROOT_TURN",turnover_daily(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_TURN",turnover_daily(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("TOTAL_BM_LITTER",tot_bm_to_litter(:,:,icarbon))
    CALL xios_orchidee_send_field("LEAF_BM_LITTER",bm_to_litter(:,:,ileaf,icarbon))
    CALL xios_orchidee_send_field("SAP_AB_BM_LITTER",bm_to_litter(:,:,isapabove,icarbon))
    CALL xios_orchidee_send_field("SAP_BE_BM_LITTER",bm_to_litter(:,:,isapbelow,icarbon))
    CALL xios_orchidee_send_field("HEART_AB_BM_LITTER",bm_to_litter(:,:,iheartabove,icarbon))
    CALL xios_orchidee_send_field("HEART_BE_BM_LITTER",bm_to_litter(:,:,iheartbelow,icarbon))
    CALL xios_orchidee_send_field("ROOT_BM_LITTER",bm_to_litter(:,:,iroot,icarbon))
    CALL xios_orchidee_send_field("FRUIT_BM_LITTER",bm_to_litter(:,:,ifruit,icarbon))
    CALL xios_orchidee_send_field("RESERVE_BM_LITTER",bm_to_litter(:,:,icarbres,icarbon))
    CALL xios_orchidee_send_field("LITTER_STR_AB",litter(:,istructural,:,iabove,icarbon))
    CALL xios_orchidee_send_field("LITTER_MET_AB",litter(:,imetabolic,:,iabove,icarbon))
    CALL xios_orchidee_send_field("LITTER_STR_BE",litter(:,istructural,:,ibelow,icarbon))
    CALL xios_orchidee_send_field("LITTER_MET_BE",litter(:,imetabolic,:,ibelow,icarbon))
    CALL xios_orchidee_send_field("DEADLEAF_COVER",deadleaf_cover)
    CALL xios_orchidee_send_field("TOTAL_SOIL_CARB",tot_litter_soil_carb)
    CALL xios_orchidee_send_field("CARBON_ACTIVE",carbon(:,iactive,:))
    CALL xios_orchidee_send_field("CARBON_SLOW",carbon(:,islow,:))
    CALL xios_orchidee_send_field("CARBON_PASSIVE",carbon(:,ipassive,:))
    CALL xios_orchidee_send_field("LITTERHUM",litterhum_daily)
    CALL xios_orchidee_send_field("TURNOVER_TIME",turnover_time)
    CALL xios_orchidee_send_field("PROD10",prod10)
    CALL xios_orchidee_send_field("FLUX10",flux10)
    CALL xios_orchidee_send_field("PROD100",prod100)
    CALL xios_orchidee_send_field("FLUX100",flux100)
    CALL xios_orchidee_send_field("CONVFLUX",convflux)
    CALL xios_orchidee_send_field("CFLUX_PROD10",cflux_prod10)
    CALL xios_orchidee_send_field("CFLUX_PROD100",cflux_prod100)
    CALL xios_orchidee_send_field("HARVEST_ABOVE",harvest_above)
    CALL xios_orchidee_send_field("VCMAX",vcmax)
    CALL xios_orchidee_send_field("AGE",age)
    CALL xios_orchidee_send_field("HEIGHT",height)
    CALL xios_orchidee_send_field("BLACK_CARBON",black_carbon)
    CALL xios_orchidee_send_field("FIREINDEX",fireindex(:,:))
! ipcc history
     CALL xios_orchidee_send_field("cVeg",SUM(tot_live_biomass(:,:,icarbon)*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cLitter",SUM(tot_litter_carb*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cSoil",SUM(tot_soil_carb*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cProduct",(prod10_total + prod100_total)/1e3)
     CALL xios_orchidee_send_field("cMassVariation",carb_mass_variation/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("lai_ipcc",SUM(lai*veget_max,dim=2)*contfrac)
     CALL xios_orchidee_send_field("gpp_ipcc",SUM(gpp_daily*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("ra",SUM((resp_maint+resp_growth)*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("npp_ipcc",SUM(npp_daily*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("rh",SUM(resp_hetero*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fFire",SUM(co2_fire*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fHarvest",harvest_above/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fLuc",cflux_prod_total/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("nbp",(SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire) &
            &        *veget_max,dim=2)-cflux_prod_total-harvest_above)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fVegLitter",SUM((tot_bm_to_litter(:,:,icarbon) + tot_turnover(:,:,icarbon))*&
          veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("fLitterSoil",SUM(SUM(soilcarbon_input,dim=2)*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("cLeaf",SUM(biomass(:,:,ileaf,icarbon)*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cWood",SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon))*&
          veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cRoot",SUM(( biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon) + &
          biomass(:,:,iheartbelow,icarbon) )*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cMisc",SUM(( biomass(:,:,icarbres,icarbon) + biomass(:,:,ifruit,icarbon))*&
          veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cLitterAbove",SUM((litter(:,istructural,:,iabove,icarbon)+&
          litter(:,imetabolic,:,iabove,icarbon))*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cLitterBelow",SUM((litter(:,istructural,:,ibelow,icarbon)+&
          litter(:,imetabolic,:,ibelow,icarbon))*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cSoilFast",SUM(carbon(:,iactive,:)*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cSoilMedium",SUM(carbon(:,islow,:)*veget_max,dim=2)/1e3*contfrac)
     CALL xios_orchidee_send_field("cSoilSlow",SUM(carbon(:,ipassive,:)*veget_max,dim=2)/1e3*contfrac)
       DO j=1,nvm
          histvar(:,j)=veget_max(:,j)*contfrac(:)*100
       ENDDO
     CALL xios_orchidee_send_field("landCoverFrac",histvar)
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_deciduous(j)) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
     CALL xios_orchidee_send_field("woodFracPrimDec",vartmp)         !! Arsene 31-07-2014 modifications name (old: treeFracPrimDec) ==> shrub+tree
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_evergreen(j)) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
     CALL xios_orchidee_send_field("woodFracPrimEver",vartmp)        !! Arsene 31-07-2014 modifications name (old: treeFracPrimEver) ==> shrub+tree
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( .NOT.(is_c4(j)) ) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
     CALL xios_orchidee_send_field("c3PftFrac",vartmp)
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( is_c4(j) ) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
     CALL xios_orchidee_send_field("c4PftFrac",vartmp)
     CALL xios_orchidee_send_field("rGrowth",SUM(resp_growth*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("rMaint",SUM(resp_maint*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("nppLeaf",SUM(bm_alloc(:,:,ileaf,icarbon)*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("nppWood",SUM(bm_alloc(:,:,isapabove,icarbon)*veget_max,dim=2)/1e3/one_day*contfrac)
     CALL xios_orchidee_send_field("nppRoot",SUM(( bm_alloc(:,:,isapbelow,icarbon) + bm_alloc(:,:,iroot,icarbon) )*&
          veget_max,dim=2)/1e3/one_day*contfrac)     


    CALL histwrite_p (hist_id_stomate, 'RESOLUTION_X', itime, &
         resolution(:,1), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'RESOLUTION_Y', itime, &
         resolution(:,2), npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CONTFRAC', itime, &
         contfrac(:), npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AB', itime, &
         litter(:,istructural,:,iabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AB', itime, &
         litter(:,imetabolic,:,iabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_BE', itime, &
         litter(:,istructural,:,ibelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_BE', itime, &
         litter(:,imetabolic,:,ibelow,icarbon), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'DEADLEAF_COVER', itime, &
         deadleaf_cover, npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'TOTAL_SOIL_CARB', itime, &
         tot_litter_soil_carb, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_ACTIVE', itime, &
         carbon(:,iactive,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_SLOW', itime, &
         carbon(:,islow,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_PASSIVE', itime, &
         carbon(:,ipassive,:), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'CARBON_ACTIVE_SURF', itime, &
         carbon_surf(:,iactive,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_SLOW_SURF', itime, &
         carbon_surf(:,islow,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CARBON_PASSIVE_SURF', itime, &
         carbon_surf(:,ipassive,:), npts*nvm, horipft_index)

!!!! Wetland CH4 methane
!pss:+
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_0', itime, &
                    ch4_flux_density_tot_0, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_0', itime, &
                    ch4_flux_density_dif_0, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_0', itime, &
                    ch4_flux_density_bub_0, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_0', itime, &
                    ch4_flux_density_pla_0, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_wet1', itime, &
                    ch4_flux_density_tot_wet1, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_wet1', itime, &
                    ch4_flux_density_dif_wet1, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_wet1', itime, &
                    ch4_flux_density_bub_wet1, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_wet1', itime, &
                    ch4_flux_density_pla_wet1, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_wet2', itime, &
                    ch4_flux_density_tot_wet2, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_wet2', itime, &
                    ch4_flux_density_dif_wet2, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_wet2', itime, &
                    ch4_flux_density_bub_wet2, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_wet2', itime, &
                    ch4_flux_density_pla_wet2, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_wet3', itime, &
                    ch4_flux_density_tot_wet3, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_wet3', itime, &
                    ch4_flux_density_dif_wet3, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_wet3', itime, &
                    ch4_flux_density_bub_wet3, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_wet3', itime, &
                    ch4_flux_density_pla_wet3, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_TOT_wet4', itime, &
                    ch4_flux_density_tot_wet4, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_DIF_wet4', itime, &
                    ch4_flux_density_dif_wet4, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_BUB_wet4', itime, &
                    ch4_flux_density_bub_wet4, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CH4_FLUX_PLA_wet4', itime, &
                    ch4_flux_density_pla_wet4, npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'TSURF_YEAR', itime, &
                    tsurf_year, npts, hori_index)
!pss:-


    CALL histwrite_p (hist_id_stomate, 'T2M_MONTH', itime, &
         t2m_month, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'T2M_WEEK', itime, &
         t2m_week, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'TSEASON', itime, &
         Tseason, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'TMIN_SPRING', itime, &
         Tmin_spring, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TMIN_SPRING_TIME', itime, &
         Tmin_spring_time, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ONSET_DATE', itime, &
         onset_date(:,:,2), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'HET_RESP', itime, &
         resp_hetero(:,:), npts*nvm, horipft_index)
! JCADD
    CALL histwrite_p(hist_id_stomate ,'T2M_14'   ,itime, &
         t2m_14, npts, hori_index)
!    CALL histwrite (hist_id_stomate, 'LITTER_RESP', itime, &
!         resp_hetero_litter_d(:,:), npts*nvm, horipft_index)
!    CALL histwrite (hist_id_stomate, 'ACTIVE_RESP', itime, &
!         resp_hetero_soil_d(:,iactive,:), npts*nvm, horipft_index)
!    CALL histwrite (hist_id_stomate, 'SLOW_RESP', itime, &
!         resp_hetero_soil_d(:,islow,:), npts*nvm, horipft_index)
!    CALL histwrite (hist_id_stomate, 'PASSIVE_RESP', itime, &
!         resp_hetero_soil_d(:,ipassive,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AVAIL', itime, &
         litter_avail(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AVAIL', itime, &
         litter_avail(:,imetabolic,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_NAVAIL', itime, &
         litter_not_avail(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_NAVAIL', itime, &
         litter_not_avail(:,imetabolic,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_STR_AVAILF', itime, &
         litter_avail_frac(:,istructural,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTER_MET_AVAILF', itime, &
         litter_avail_frac(:,imetabolic,:), npts*nvm, horipft_index)

    CALL histwrite_p (hist_id_stomate, 'N_LIMFERT', itime, &
         N_limfert, npts*nvm, horipft_index)
!WRITE (numout,*) 'endjcwrite'
! ENDJCADD
    CALL histwrite_p (hist_id_stomate, 'BLACK_CARBON', itime, &
         black_carbon, npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'FIREINDEX', itime, &
         fireindex(:,:), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LITTERHUM', itime, &
         litterhum_daily, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CO2_FIRE', itime, &
         co2_fire, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CO2_TAKEN', itime, &
         co2_to_bm, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite_p (hist_id_stomate, 'CONVFLUX', itime, &
         convflux, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD10', itime, &
         cflux_prod10, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'CFLUX_PROD100', itime, &
         cflux_prod100, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'HARVEST_ABOVE', itime, &
         harvest_above, npts, hori_index)

    CALL histwrite_p (hist_id_stomate, 'LAI', itime, &
         lai, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FPC_MAX', itime, &
         fpc_max, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MAXFPC_LASTYEAR', itime, &
         maxfpc_lastyear, npts*nvm, horipft_index) 
    CALL histwrite_p (hist_id_stomate, 'VEGET_MAX', itime, &
         veget_max, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'NPP', itime, &
         npp_daily, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'GPP', itime, &
         gpp_daily, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'IND', itime, &
         ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'CN_IND', itime, &
         cn_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'WOODMASS_IND', itime, &
         woodmass_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_M', itime, &
         tot_live_biomass(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_M', itime, &
         biomass(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_M_AB', itime, &
         biomass(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_M_BE', itime, &
         biomass(:,:,isapbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_M_AB', itime, &
         biomass(:,:,iheartabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_M_BE', itime, &
         biomass(:,:,iheartbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_M', itime, &
         biomass(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_M', itime, &
         biomass(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'RESERVE_M', itime, &
         biomass(:,:,icarbres,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_TURN', itime, &
         tot_turnover(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_TURN', itime, &
         turnover_daily(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_AB_TURN', itime, &
         turnover_daily(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_TURN', itime, &
         turnover_daily(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_TURN', itime, &
         turnover_daily(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TOTAL_BM_LITTER', itime, &
         tot_bm_to_litter(:,:,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'LEAF_BM_LITTER', itime, &
         bm_to_litter(:,:,ileaf,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_AB_BM_LITTER', itime, &
         bm_to_litter(:,:,isapabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'SAP_BE_BM_LITTER', itime, &
         bm_to_litter(:,:,isapbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_AB_BM_LITTER', itime, &
         bm_to_litter(:,:,iheartabove,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEART_BE_BM_LITTER', itime, &
         bm_to_litter(:,:,iheartbelow,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ROOT_BM_LITTER', itime, &
         bm_to_litter(:,:,iroot,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'FRUIT_BM_LITTER', itime, &
         bm_to_litter(:,:,ifruit,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'RESERVE_BM_LITTER', itime, &
         bm_to_litter(:,:,icarbres,icarbon), npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MAINT_RESP', itime, &
         resp_maint, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'GROWTH_RESP', itime, &
         resp_growth, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'AGE', itime, &
         age, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'HEIGHT', itime, &
         height, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MOISTRESS', itime, &
         moiavail_week, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'MOISTRESS_month', itime, &    !! Arsene 13-05-2014
         moiavail_month, npts*nvm, horipft_index)                     !! Arsene 13-05-2014
    CALL histwrite_p (hist_id_stomate, 'VCMAX', itime, &
         vcmax, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'TURNOVER_TIME', itime, &
         turnover_time, npts*nvm, horipft_index)
    ! land cover change
    CALL histwrite_p (hist_id_stomate, 'PROD10', itime, &
         prod10, npts*11, horip11_index)
    CALL histwrite_p (hist_id_stomate, 'PROD100', itime, &
         prod100, npts*101, horip101_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX10', itime, &
         flux10, npts*10, horip10_index)
    CALL histwrite_p (hist_id_stomate, 'FLUX100', itime, &
         flux100, npts*100, horip100_index)
!!DZADD
!    CALL histwrite_p (hist_id_stomate, 'LEAF_FRAC1', itime, leaf_frac(:,:,1), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'LEAF_FRAC2', itime, leaf_frac(:,:,2), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'LEAF_FRAC3', itime, leaf_frac(:,:,3), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'LEAF_FRAC4', itime, leaf_frac(:,:,4), npts*nvm, horipft_index)
!!ENDDZADD

    IF ( hist_id_stomate_IPCC > 0 ) THEN
       vartmp(:)=SUM(tot_live_biomass(:,:,icarbon)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cVeg", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_litter_carb*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(tot_soil_carb*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=(prod10_total + prod100_total)/1e3
       CALL histwrite_p (hist_id_stomate_IPCC, "cProduct", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=carb_mass_variation/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cMassVariation", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(lai*veget_max,dim=2)*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "lai", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(gpp_daily*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "gpp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((resp_maint+resp_growth)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "ra", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(npp_daily*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "npp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_hetero*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "rh", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(co2_fire*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "fFire", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=harvest_above/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "fHarvest", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=cflux_prod_total/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "fLuc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=(SUM((gpp_daily-(resp_maint+resp_growth+resp_hetero)-co2_fire) &
            &        *veget_max,dim=2)-cflux_prod_total-harvest_above)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "nbp", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((tot_bm_to_litter(:,:,icarbon) + tot_turnover(:,:,icarbon))*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "fVegLitter", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(SUM(soilcarbon_input,dim=2)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "fLitterSoil", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(biomass(:,:,ileaf,icarbon)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((biomass(:,:,isapabove,icarbon)+biomass(:,:,iheartabove,icarbon))*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cWood", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon) + biomass(:,:,iheartbelow,icarbon) ) &
            &        *veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cRoot", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( biomass(:,:,icarbres,icarbon) + biomass(:,:,ifruit,icarbon))*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cMisc", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,iabove,icarbon)+litter(:,imetabolic,:,iabove,icarbon))*&
            veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitterAbove", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM((litter(:,istructural,:,ibelow,icarbon)+litter(:,imetabolic,:,ibelow,icarbon))*&
            veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cLitterBelow", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,iactive,:)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilFast", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,islow,:)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilMedium", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(carbon(:,ipassive,:)*veget_max,dim=2)/1e3*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "cSoilSlow", itime, &
            vartmp, npts, hori_index)
       DO j=1,nvm
          histvar(:,j)=veget_max(:,j)*contfrac(:)*100
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "landCoverFrac", itime, &
            histvar, npts*nvm, horipft_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_deciduous(j)) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "woodFracPrimDec", itime, &    !! Arsene 31-07-2014 modifications name (old: treeFracPrimDec) ==> shrub+tree
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF (is_evergreen(j)) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "woodFracPrimEver", itime, &   !! Arsene 31-07-2014 modifications name (old: treeFracPrimEver) ==> shrub+tree
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( .NOT.(is_c4(j)) ) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "c3PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=zero
       DO j = 2,nvm
          IF ( is_c4(j) ) THEN
             vartmp(:) = vartmp(:) + veget_max(:,j)*contfrac*100
          ENDIF
       ENDDO
       CALL histwrite_p (hist_id_stomate_IPCC, "c4PftFrac", itime, &
            vartmp, npts, hori_index)
       !-
       vartmp(:)=SUM(resp_growth*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "rGrowth", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(resp_maint*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "rMaint", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(bm_alloc(:,:,ileaf,icarbon)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "nppLeaf", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(bm_alloc(:,:,isapabove,icarbon)*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "nppWood", itime, &
            vartmp, npts, hori_index)
       vartmp(:)=SUM(( bm_alloc(:,:,isapbelow,icarbon) + bm_alloc(:,:,iroot,icarbon) )*veget_max,dim=2)/1e3/one_day*contfrac
       CALL histwrite_p (hist_id_stomate_IPCC, "nppRoot", itime, &
            vartmp, npts, hori_index)

       CALL histwrite_p (hist_id_stomate_IPCC, 'RESOLUTION_X', itime, &
            resolution(:,1), npts, hori_index)
       CALL histwrite_p (hist_id_stomate_IPCC, 'RESOLUTION_Y', itime, &
            resolution(:,2), npts, hori_index)
       CALL histwrite_p (hist_id_stomate_IPCC, 'CONTFRAC', itime, &
            contfrac(:), npts, hori_index)

    ENDIF

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving stomate_lpj'

  END SUBROUTINE StomateLpj


!! ================================================================================================================================
!! SUBROUTINE   : harvest
!!
!>\BRIEF        Harvest of croplands
!!
!! DESCRIPTION  : To take into account biomass harvest from crop (mainly to take 
!! into account for the reduced litter input and then decreased soil carbon. it is a 
!! constant (40\%) fraction of above ground biomass.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): ::harvest_above the harvested biomass
!!
!! REFERENCE(S) :
!! - Piao, S., P. Ciais, P. Friedlingstein, N. de Noblet-Ducoudre, P. Cadule, N. Viovy, and T. Wang. 2009. 
!!   Spatiotemporal patterns of terrestrial carbon cycle during the 20th century. Global Biogeochemical 
!!   Cycles 23:doi:10.1029/2008GB003339.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE harvest(npts, dt_days, veget_max, &
       bm_to_litter, turnover_daily, &
       harvest_above)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER, INTENT(in)                                    :: npts            !! Domain size (unitless) 
    REAL(r_std), INTENT(in)                                :: dt_days         !! Time step (days)                               
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)           :: veget_max       !! new "maximal" coverage fraction of a PFT (LAI -> 
                                                                              !! infinity) on ground @tex $(m^2 m^{-2})$ @endtex 
    
   !! 0.2 Output variables
   
   !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! [DISPENSABLE] conversion of biomass to litter 
                                                                                     !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily   !! Turnover rates 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    REAL(r_std), DIMENSION(npts), INTENT(inout)            :: harvest_above    !! harvest above ground biomass for agriculture 
                                                                               !! @tex $(gC m^{-2} dtslow^{-1})$ @endtex 
    !! 0.4 Local variables

    INTEGER(i_std)                                         :: i, j, k, l, m    !! indices                       
    REAL(r_std)                                            :: above_old        !! biomass of previous time step 
                                                                               !! @tex $(gC m^{-2})$ @endtex 
!_ ================================================================================================================================

  !! 1. Yearly initialisation

    above_old             = zero
    harvest_above         = zero

    DO i = 1, npts
       DO j = 1,nvm
          IF (.NOT. natural(j)) THEN
             above_old = turnover_daily(i,j,ileaf,icarbon) + turnover_daily(i,j,isapabove,icarbon) + &
                  &       turnover_daily(i,j,iheartabove,icarbon) + turnover_daily(i,j,ifruit,icarbon) + &
                  &       turnover_daily(i,j,icarbres,icarbon) + turnover_daily(i,j,isapbelow,icarbon) + &
                  &       turnover_daily(i,j,iheartbelow,icarbon) + turnover_daily(i,j,iroot,icarbon)

             turnover_daily(i,j,ileaf,icarbon) = turnover_daily(i,j,ileaf,icarbon)*frac_turnover_daily
             turnover_daily(i,j,isapabove,icarbon) = turnover_daily(i,j,isapabove,icarbon)*frac_turnover_daily
             turnover_daily(i,j,isapbelow,icarbon) = turnover_daily(i,j,isapbelow,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iheartabove,icarbon) = turnover_daily(i,j,iheartabove,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iheartbelow,icarbon) = turnover_daily(i,j,iheartbelow,icarbon)*frac_turnover_daily
             turnover_daily(i,j,iroot,icarbon) = turnover_daily(i,j,iroot,icarbon)*frac_turnover_daily
             turnover_daily(i,j,ifruit,icarbon) = turnover_daily(i,j,ifruit,icarbon)*frac_turnover_daily
             turnover_daily(i,j,icarbres,icarbon) = turnover_daily(i,j,icarbres,icarbon)*frac_turnover_daily
             harvest_above(i)  = harvest_above(i) + veget_max(i,j) * above_old *(un - frac_turnover_daily)
          ENDIF
       ENDDO
    ENDDO

!!$    harvest_above = harvest_above
  END SUBROUTINE harvest
END MODULE stomate_lpj
