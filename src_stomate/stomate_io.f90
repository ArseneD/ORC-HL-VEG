!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_io.f90 $ 
!< $Date: 2013-07-26 18:45:06 +0200 (Fri, 26 Jul 2013) $
!< $Author: josefine.ghattas $
!< $Revision: 1392 $
! IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
MODULE stomate_io
  !---------------------------------------------------------------------
  !- Not all variables saved in the start files are absolutely necessary.
  !- However, Sechiba's and Stomate's PFTs are not necessarily identical,
  !- and for that case this information needs to be saved.
  !---------------------------------------------------------------------
  USE stomate_data
  USE constantes
  USE constantes_soil
  USE mod_orchidee_para
  USE ioipsl_para 

  !-
  IMPLICIT NONE
  !-
  PRIVATE
  PUBLIC readstart, writerestart, readbc,get_reftemp_clear
  !-
  ! first call?
  !-
  LOGICAL,SAVE :: firstcall = .TRUE.
!$OMP THREADPRIVATE(firstcall)
  !-
  ! reference temperature (K)
  !-
  REAL(r_std),ALLOCATABLE,DIMENSION(:),SAVE :: trefe
!$OMP THREADPRIVATE(trefe)
  !-
CONTAINS
  !-
  !===
  !-
  SUBROUTINE readstart &
       & (npts, index, lalo, resolution, day_counter, dt_days, date, &
       &  ind, adapted, regenerate, moiavail_daily, gdd_init_date, litterhum_daily, &
       &  t2m_daily, t2m_min_daily, tsurf_daily, tsoil_daily, &
       &  soilhum_daily, precip_daily, gpp_daily, npp_daily, &
       &  turnover_daily, moiavail_month, moiavail_week,&
       &  t2m_longterm, tlong_ref, t2m_month, t2m_week,&
       &  tsoil_month, soilhum_month, fireindex, firelitter,  &
       &  maxmoiavail_lastyear, maxmoiavail_thisyear, &
       &  minmoiavail_lastyear, minmoiavail_thisyear, &
       &  maxgppweek_lastyear, maxgppweek_thisyear, &
       &  gdd0_lastyear, gdd0_thisyear, precip_lastyear, precip_thisyear, &
       &  gdd_m5_dormance,  gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       &  PFTpresent, npp_longterm, lm_lastyearmax, lm_thisyearmax, &
       &  maxfpc_lastyear, maxfpc_thisyear, &
       &  turnover_longterm, gpp_week, biomass, resp_maint_part, &
       &  leaf_age, leaf_frac, senescence, when_growthinit, age, &
       &  resp_hetero, resp_maint, resp_growth, co2_fire, co2_to_bm_dgvm, &
       &  veget_lastlight, everywhere, need_adjacent, RIP_time, &
       &  time_hum_min, hum_min_dormance, &
       &  litterpart, litter, dead_leaves, &
       &  carbon, black_carbon, lignin_struc,turnover_time, &
       &  prod10,prod100,flux10, flux100, &
       &  convflux, cflux_prod10, cflux_prod100, bm_to_litter, carb_mass_total,&
       &  Tseason, Tseason_length, Tseason_tmp, &
       &  Tmin_spring, Tmin_spring_time, begin_leaves, onset_date, &
       &  global_years, ok_equilibrium, nbp_accu, nbp_flux, &
       &  MatrixV, VectorU, previous_stock, current_stock,&
       &  npp0_cumul,snowtemp_min,snowdz_min, dia_cut, &  !! Arsene 25-06-2014 NPPcumul ADD  !! Arsene 19-08-2014 Add snowtemp and snowdz !! Arsene 27-08-2015 add dia_cut
       &  deepC_a, deepC_s, deepC_p, O2_soil, CH4_soil, O2_snow, CH4_snow, &
       &  thawed_humidity, depth_organic_soil, altmax, fixed_cryoturbation_depth, & !pss+
       &  uo_0, uold2_0, uo_wet1, uold2_wet1, uo_wet2, uold2_wet2, uo_wet3, uold2_wet3, &
       &  uo_wet4, uold2_wet4, tsurf_year, &!) !pss-
!JCADD
       &  wshtotsum, sr_ugb, sla_calc, nb_ani, grazed_frac, &
       &  import_yield, t2m_14, litter_not_avail)
!ENDJCADD
    !---------------------------------------------------------------------
    !- read start file
    !---------------------------------------------------------------------
    !-
    ! 0 declarations
    !-
    ! 0.1 input
    !-
    ! Domain size
    INTEGER(i_std),INTENT(in) :: npts
    ! Indices of the points on the map
    INTEGER(i_std),DIMENSION(npts),INTENT(in) :: index
    ! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std),DIMENSION(npts,2),INTENT(in) :: lalo
    ! size in x an y of the grid (m)
    REAL(r_std),DIMENSION(npts,2),INTENT(in) :: resolution
    !-
    ! 0.2 output
    !-
    ! counts time until next STOMATE time step
    REAL(r_std),INTENT(out) :: day_counter
    ! time step of STOMATE in days
    REAL(r_std),INTENT(out) :: dt_days
    ! date (d)
    INTEGER(i_std),INTENT(out) :: date
    ! density of individuals (1/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: ind
    ! Winter too cold? between 0 and 1
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: adapted
    ! Winter sufficiently cold? between 0 and 1
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: regenerate
    ! daily moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: moiavail_daily
    ! date for beginning of gdd count
    REAL(r_std),DIMENSION(npts,2),INTENT(out) :: gdd_init_date
    ! daily litter humidity
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: litterhum_daily
    ! daily 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: t2m_daily
    ! daily minimum 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: t2m_min_daily
    ! daily surface temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: tsurf_daily
    ! daily soil temperatures (K)
    REAL(r_std),DIMENSION(npts,nbdl),INTENT(out) :: tsoil_daily
    ! daily soil humidity
    REAL(r_std),DIMENSION(npts,nbdl),INTENT(out) :: soilhum_daily
    ! daily precipitations (mm/day) (for phenology)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: precip_daily
    ! daily gross primary productivity (gC/m**2/day)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: gpp_daily
    ! daily net primary productivity (gC/m**2/day)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: npp_daily
    ! daily turnover rates (gC/m**2/day)
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(out) :: turnover_daily
    ! "monthly" moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: moiavail_month
    ! "weekly" moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: moiavail_week
    ! "long term" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: t2m_longterm
    ! "monthly" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: t2m_month
    ! "seasonal" 2 meter temperatures (K) 
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: Tseason
    ! temporary variable to calculate Tseason
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: Tseason_length
    ! temporary variable to calculate Tseason
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: Tseason_tmp

    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)  :: Tmin_spring
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)  :: Tmin_spring_time
    REAL(r_std),DIMENSION(npts,nvm,2),INTENT(out)  :: onset_date
    LOGICAL,DIMENSION(npts,nvm),INTENT(out)      :: begin_leaves

    ! "weekly" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: t2m_week
    ! "monthly" soil temperatures (K)
    REAL(r_std),DIMENSION(npts,nbdl),INTENT(out) :: tsoil_month
    ! "monthly" soil humidity
    REAL(r_std),DIMENSION(npts,nbdl),INTENT(out) :: soilhum_month
    ! Probability of fire
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: fireindex
    ! Longer term total litter above the ground, gC/m**2 of ground
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: firelitter
    ! last year's maximum moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: maxmoiavail_lastyear
    ! this year's maximum moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: maxmoiavail_thisyear
    ! last year's minimum moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: minmoiavail_lastyear
    ! this year's minimum moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: minmoiavail_thisyear
    ! last year's maximum weekly GPP
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: maxgppweek_lastyear
    ! this year's maximum weekly GPP
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: maxgppweek_thisyear
    ! last year's annual GDD0
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: gdd0_lastyear
    ! this year's annual GDD0
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: gdd0_thisyear
    ! last year's annual precipitation (mm/year)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: precip_lastyear
    ! this year's annual precipitation (mm/year)
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: precip_thisyear
    ! growing degree days, threshold -5 deg C (for phenology)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: gdd_m5_dormance
    ! growing degree days, from begin of season
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: gdd_from_growthinit
    ! growing degree days since midwinter (for phenology)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: gdd_midwinter
    ! number of chilling days since leaves were lost (for phenology)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: ncd_dormance
    ! number of growing days, threshold -5 deg C (for phenology)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: ngd_minus5
    ! PFT exists (equivalent to fpc_max > 0 for natural PFTs)
    LOGICAL,DIMENSION(npts,nvm),INTENT(out)    :: PFTpresent
    ! "long term" net primary productivity (gC/m**2/year)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: npp_longterm

    ! "counter" null or negative npp (days)                           !! Arsene 25-06-2014 NPPcumul
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)     :: npp0_cumul     !! Arsene 25-06-2014 NPPcumul
    ! Min daily snow layer temperature (K)                            !! Arsene 19-08-2014 Add
    REAL(r_std), DIMENSION(npts,nsnow), INTENT(out) :: snowtemp_min   !! Arsene 19-08-2014 Add
    ! Min daily snow layer thicknesse (cm?)                           !! Arsene 19-08-2014 Add
    REAL(r_std), DIMENSION(npts,nsnow), INTENT(out) :: snowdz_min     !! Arsene 19-08-2014 Add
    ! Fix diameter of vegetation (for shrub) after loss biomass (above snow) !! Arsene 27-08-2015 add dia_cut
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)   :: dia_cut               !! Arsene 27-08-2015 add dia_cut

    ! last year's maximum leaf mass, for each PFT (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: lm_lastyearmax
    ! this year's maximum leaf mass, for each PFT (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: lm_thisyearmax
    ! last year's maximum fpc for each natural PFT, on ground
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: maxfpc_lastyear
    ! this year's maximum fpc for each PFT,
    ! on *total* ground (see stomate_season)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: maxfpc_thisyear
    ! "long term" turnover rate (gC/m**2/year)
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(out) :: turnover_longterm
    ! "weekly" GPP (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: gpp_week
    ! biomass (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(out) :: biomass
    ! maintenance resp (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm,nparts),INTENT(out) :: resp_maint_part
    ! leaf age (days)
    REAL(r_std),DIMENSION(npts,nvm,nleafages),INTENT(out) :: leaf_age
    ! fraction of leaves in leaf age class
    REAL(r_std),DIMENSION(npts,nvm,nleafages),INTENT(out) :: leaf_frac
    ! is the plant senescent ? 
    !(only for deciduous trees - carbohydrate reserve)
    LOGICAL,DIMENSION(npts,nvm),INTENT(out) :: senescence
    ! how many days ago was the beginning of the growing season
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: when_growthinit
    ! mean age (years)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: age
    ! heterotrophic respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: resp_hetero
    ! maintenance respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: resp_maint
    ! growth respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: resp_growth
    ! carbon emitted into the atmosphere by fire (living and dead biomass)
    ! (in gC/m**2/time step)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: co2_fire
    ! biomass uptaken (gC/(m**2 of total ground)/day)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: co2_to_bm_dgvm
    ! vegetation fractions (on ground) after last light competition
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: veget_lastlight
    ! is the PFT everywhere in the grid box or very localized
    ! (after its introduction)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: everywhere
    ! in order for this PFT to be introduced,
    ! does it have to be present in an adjacent grid box?
    LOGICAL,DIMENSION(npts,nvm),INTENT(out) :: need_adjacent
    ! How much time ago was the PFT eliminated for the last time (y)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: RIP_time
    ! time elapsed since strongest moisture availability (d)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: time_hum_min
    ! minimum moisture during dormance
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: hum_min_dormance
    ! fraction of litter above the ground belonging to different PFTs
    ! separated for natural and agricultural PFTs.
    REAL(r_std),DIMENSION(npts,nvm,nlitt),INTENT(out) :: litterpart
    ! metabolic and structural litter, natural and agricultural,
    ! above and below ground (gC/m**2)
    REAL(r_std),DIMENSION(npts,nlitt,nvm,nlevs,nelements),INTENT(out):: litter
    ! dead leaves on ground, per PFT, metabolic and structural,
    ! in gC/(m**2 of ground)
    REAL(r_std),DIMENSION(npts,nvm,nlitt),INTENT(out) :: dead_leaves
    ! carbon pool: active, slow, or passive, (gC/m**2)
    REAL(r_std),DIMENSION(npts,ncarb,nvm),INTENT(out) :: carbon
    ! black carbon on the ground (gC/(m**2 of total ground))
    REAL(r_std),DIMENSION(npts),INTENT(out)                 :: black_carbon
    ! ratio Lignine/Carbon in structural litter, above and below ground,(gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm,nlevs),INTENT(out) :: lignin_struc
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out) :: turnover_time

    ! For Spinup matrix resolution
    INTEGER(i_std), INTENT(out) :: global_years   
    LOGICAL, DIMENSION(npts), INTENT(out) :: ok_equilibrium
    REAL(r_std), DIMENSION(npts), INTENT(out) :: nbp_accu  !! Accumulated Net Biospheric Production over the year
    REAL(r_std), DIMENSION(npts), INTENT(out) :: nbp_flux  !! Net Biospheric Production over the year
    !-
    REAL(r_std), DIMENSION(npts,nvm,nbpools,nbpools), INTENT(out) :: MatrixV
    REAL(r_std), DIMENSION(npts,nvm,nbpools), INTENT(out) :: VectorU
    REAL(r_std), DIMENSION(npts,nvm,nbpools), INTENT(out) :: previous_stock
    REAL(r_std), DIMENSION(npts,nvm,nbpools), INTENT(out) :: current_stock    

    !-
    REAL(r_std), DIMENSION(npts,ndeep,nvm),INTENT(inout) :: deepC_a
    REAL(r_std), DIMENSION(npts,ndeep,nvm),INTENT(inout) :: deepC_s
    REAL(r_std), DIMENSION(npts,ndeep,nvm),INTENT(inout) :: deepC_p
    REAL(r_std), DIMENSION(npts,ndeep,nvm),INTENT(inout) :: O2_soil
    REAL(r_std), DIMENSION(npts,ndeep,nvm),INTENT(inout) :: CH4_soil
    REAL(r_std), DIMENSION(npts,nsnow,nvm),INTENT(inout) :: O2_snow
    REAL(r_std), DIMENSION(npts,nsnow,nvm),INTENT(inout) :: CH4_snow
    REAL(r_std), DIMENSION(npts),INTENT(inout)           :: thawed_humidity
    REAL(r_std), DIMENSION(npts),INTENT(inout)           :: depth_organic_soil
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)      :: altmax     !! Active layer thickness (m) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)      :: fixed_cryoturbation_depth !! Depth to hold cryoturbation to for fixed runs 

    !pss:+
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out) :: uo_0
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out)  :: uold2_0
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out) :: uo_wet1
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out)  :: uold2_wet1
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out) :: uo_wet2
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out)  :: uold2_wet2
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out) :: uo_wet3
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out)  :: uold2_wet3
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out) :: uo_wet4
    REAL(r_std),DIMENSION(npts,nvert),INTENT(out)  :: uold2_wet4
    REAL(r_std),DIMENSION(npts),INTENT(out)      :: tsurf_year
    !pss:-



    !-
    ! 0.3 not necessarily output
    !-
    ! "long term" reference 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(inout) :: tlong_ref

!JCADD
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)    ::  sla_calc
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)    ::  wshtotsum
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)  ::  sr_ugb
!    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)  ::  compt_ugb
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)    ::  nb_ani
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)  ::  grazed_frac
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)  ::  import_yield
    REAL(r_std),DIMENSION(npts),INTENT(out)        ::  t2m_14
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(out)  ::  litter_not_avail
!ENDJCADD
    !-
    ! 0.4 local
    !-
    ! date, real
    REAL(r_std) :: date_real
    ! PFT exists (equivalent to fpc_max > 0 for natural PFTs), real
    REAL(r_std),DIMENSION(npts,nvm) :: PFTpresent_real
    ! is the plant senescent ?
    ! (only for deciduous trees - carbohydrate reserve), real
    REAL(r_std),DIMENSION(npts,nvm) :: senescence_real

    REAL(r_std),DIMENSION(npts,nvm) :: begin_leaves_real
    ! in order for this PFT to be introduced,
    ! does it have to be present in an adjacent grid box? - real
    REAL(r_std),DIMENSION(npts,nvm) :: need_adjacent_real
    ! To store variables names for I/O
    CHARACTER(LEN=80) :: var_name
    ! string suffix indicating an index
    CHARACTER(LEN=10) :: part_str
    ! string suffix indicating litter type
    CHARACTER(LEN=3),DIMENSION(nlitt) :: litter_str
    ! string suffix indicating level
    CHARACTER(LEN=2),DIMENSION(nlevs) :: level_str
    ! temporary storage
    REAL(r_std),DIMENSION(1) :: xtmp
    ! index
    INTEGER(i_std) :: j,k,l,m
    ! reference temperature (K)
    REAL(r_std),DIMENSION(npts) :: tref

    CHARACTER(LEN=1),DIMENSION(nelements) :: element_str   !! string suffix indicating element
    REAL(r_std), DIMENSION(1) :: temp_global_years
    CHARACTER(LEN=6), DIMENSION(nbpools) :: pools_str
    REAL(r_std), DIMENSION(npts) :: ok_equilibrium_real    

    ! Permafrost carbon processes
    LOGICAL :: read_input_deepC_a
    LOGICAL :: read_input_deepC_s
    LOGICAL :: read_input_deepC_p
    LOGICAL :: read_input_thawed_humidity
    LOGICAL :: read_input_depth_organic_soil
    real(r_std) :: deepC_a_init, deepC_s_init, deepC_p_init
    real(r_std) :: thawed_humidity_input = 0.5
    LOGICAL ::  reset_thawed_humidity = .FALSE.

    ! Wetland CH4
    INTEGER(i_std) :: nivo !pss:+-
    REAL(r_std) :: CH4atmo_CONC !pss:+-

    ! land cover change variables 
    ! products remaining in the 10/100 year-turnover pool after the annual release for each compartment
    ! (10 or 100 + 1 : input from year of land cover change)
    REAL(r_std),DIMENSION(npts,0:10),INTENT(out)                           :: prod10
    REAL(r_std),DIMENSION(npts,0:100),INTENT(out)                          :: prod100
    ! annual release from the 10/100 year-turnover pool compartments
    REAL(r_std),DIMENSION(npts,10),INTENT(out)                           :: flux10
    REAL(r_std),DIMENSION(npts,100),INTENT(out)                          :: flux100
    REAL(r_std), DIMENSION(npts), INTENT(out)                            :: convflux
    REAL(r_std), DIMENSION(npts), INTENT(out)                            :: cflux_prod10
    REAL(r_std), DIMENSION(npts), INTENT(out)                            :: cflux_prod100
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(out)         :: bm_to_litter
    REAL(r_std),DIMENSION(npts),INTENT(out)                              :: carb_mass_total
    !---------------------------------------------------------------------
    IF (bavard >= 3) WRITE(numout,*) 'Entering readstart'

    ! near surface CH4 concentration
    !pss:+
    CH4atmo_CONC=0.0017
    CALL GETIN('CH4atmo_CONC', CH4atmo_CONC)
    !pss:-

    !-
    ! 0 When the vegetation is dynamic,
    !   the long-term reference temperature is prognostic.
    !   In this case, it is read from the restart file.
    !   If the corresponding field does not exist in the restart file,
    !   read it from another file in order to initialize it correctly.
    !-
    CALL get_reftemp( npts, lalo, resolution, tref )
    !-
    ! 1 string definitions
    !-
    DO l=1,nlitt
       IF     (l == imetabolic) THEN
          litter_str(l) = 'met'
       ELSEIF (l == istructural) THEN
          litter_str(l) = 'str'
       ELSE
          STOP 'Define litter_str'
       ENDIF
    ENDDO
    !-
    DO l=1,nlevs
       IF     (l == iabove) THEN
          level_str(l) = 'ab'
       ELSEIF (l == ibelow) THEN
          level_str(l) = 'be'
       ELSE
          STOP 'Define level_str'
       ENDIF
    ENDDO

    pools_str(1:nbpools) =(/'str_ab','str_be','met_ab','met_be','actif ','slow  ','passif'/)

    !-
    DO l=1,nelements
       IF     (l == icarbon) THEN
          element_str(l) = ''
!!$       ELSEIF (l == initrogen) THEN
!!$          element_str(l) = '_n'
       ELSE
          STOP 'Define element_str'
       ENDIF
    ENDDO
    !-
    ! 2 run control
    !-
    ! 2.1 day counter
    !-
    IF (is_root_prc) THEN
       var_name = 'day_counter'
       CALL restget (rest_id_stomate, var_name, 1   , 1     , 1, itime, &
            &                 .TRUE., xtmp)
       day_counter = xtmp(1)
       IF (day_counter == val_exp) day_counter = un
    ENDIF
    CALL bcast(day_counter)
    !-
    ! 2.2 time step of STOMATE in days
    !-
    IF (is_root_prc) THEN
       var_name = 'dt_days'
       CALL restget (rest_id_stomate, var_name, 1   , 1     , 1, itime, &
            &                 .TRUE., xtmp)
       dt_days = xtmp(1)
       IF (dt_days == val_exp) dt_days = un
    ENDIF
    CALL bcast(dt_days)
    !-
    ! 2.3 date
    !-
    IF (is_root_prc) THEN
       var_name = 'date'
       CALL restget (rest_id_stomate, var_name, 1   , 1     , 1, itime, &
            &                 .TRUE., xtmp)
       date_real = xtmp(1)
       IF (date_real == val_exp) date_real = zero
       date = NINT(date_real)
    ENDIF
    CALL bcast(date)
    !-
    ! 3 daily meteorological variables
    !-
    moiavail_daily(:,:) = val_exp
    var_name = 'moiavail_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., moiavail_daily, 'gather', nbp_glo, index_g)
    IF (ALL(moiavail_daily(:,:) == val_exp)) moiavail_daily(:,:) = zero
    !-
    gdd_init_date(:,:) = val_exp
    var_name = 'gdd_init_date'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 2 , 1, itime, &
         &              .TRUE., gdd_init_date, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_init_date(:,1) == val_exp)) gdd_init_date(:,1) = 365.
    !-
    litterhum_daily(:) = val_exp
    var_name = 'litterhum_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., litterhum_daily, 'gather', nbp_glo, index_g)
    IF (ALL(litterhum_daily(:) == val_exp)) litterhum_daily(:) = zero
    !-
    t2m_daily(:) = val_exp
    var_name = 't2m_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., t2m_daily, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_daily(:) == val_exp)) t2m_daily(:) = zero
    !-
    t2m_min_daily(:) = val_exp
    var_name = 't2m_min_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., t2m_min_daily, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_min_daily(:) == val_exp)) t2m_min_daily(:) = large_value
    !-
    tsurf_daily(:) = val_exp
    var_name = 'tsurf_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., tsurf_daily, 'gather', nbp_glo, index_g)
    IF (ALL(tsurf_daily(:) == val_exp)) tsurf_daily(:) = tref(:)
    !-
    tsoil_daily(:,:) = val_exp
    var_name = 'tsoil_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nbdl, 1, itime, &
         &                .TRUE., tsoil_daily, 'gather', nbp_glo, index_g)
    IF (ALL(tsoil_daily(:,:) == val_exp)) tsoil_daily(:,:) = zero
    !-
    soilhum_daily(:,:) = val_exp
    var_name = 'soilhum_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nbdl, 1, itime, &
         &                .TRUE., soilhum_daily, 'gather', nbp_glo, index_g)
    IF (ALL(soilhum_daily(:,:) == val_exp)) soilhum_daily(:,:) = zero
    !-
    precip_daily(:) = val_exp
    var_name = 'precip_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., precip_daily, 'gather', nbp_glo, index_g)
    IF (ALL(precip_daily(:) == val_exp)) precip_daily(:) = zero
    !-
!! Arsene 19-08-2014 Add snowtemp_min and snowdz_min
    snowtemp_min(:,:) = val_exp                                             !! Arsene 19-08-2014 Add
    var_name = 'snowtemp_min'                                               !! Arsene 19-08-2014 Add
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nsnow , 1, itime, & !! Arsene 19-08-2014 Add
         &                .TRUE., snowtemp_min, 'gather', nbp_glo, index_g) !! Arsene 19-08-2014 Add
    IF (ALL(snowtemp_min(:,:) == val_exp)) snowtemp_min(:,:) = large_value  !! Arsene 19-08-2014 Add
    !-                                                                      !! Arsene 19-08-2014 Add
    snowdz_min(:,:) = val_exp                                               !! Arsene 19-08-2014 Add
    var_name = 'snowdz_min'                                                 !! Arsene 19-08-2014 Add
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nsnow , 1, itime, & !! Arsene 19-08-2014 Add
         &                .TRUE., snowtemp_min, 'gather', nbp_glo, index_g) !! Arsene 19-08-2014 Add
    IF (ALL(snowdz_min(:,:) == val_exp)) snowdz_min(:,:) = large_value      !! Arsene 19-08-2014 Add
!! Arsene 19-08-2014 Add snowtemp_min and snowdz_min
    !-
    ! 4 productivities
    !-
    gpp_daily(:,:) = val_exp
    var_name = 'gpp_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gpp_daily, 'gather', nbp_glo, index_g)
    IF (ALL(gpp_daily(:,:) == val_exp)) gpp_daily(:,:) = zero
    !-
    npp_daily(:,:) = val_exp
    var_name = 'npp_daily'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., npp_daily, 'gather', nbp_glo, index_g)
    IF (ALL(npp_daily(:,:) == val_exp)) npp_daily(:,:) = zero
    !-
    turnover_daily(:,:,:,:) = val_exp
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'turnover_daily_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                .TRUE., turnover_daily(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(turnover_daily(:,:,k,l) == val_exp)) &
               &       turnover_daily(:,:,k,l) = zero
       ENDDO
    END DO
    !-
    ! 5 monthly meteorological variables
    !-
    moiavail_month(:,:) = val_exp
    var_name = 'moiavail_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., moiavail_month, 'gather', nbp_glo, index_g)
    IF (ALL(moiavail_month(:,:) == val_exp)) moiavail_month(:,:) = zero
    !-
    moiavail_week(:,:) = val_exp
    var_name = 'moiavail_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., moiavail_week, 'gather', nbp_glo, index_g)
    IF (ALL(moiavail_week(:,:) == val_exp)) moiavail_week(:,:) = zero
    !-
    t2m_longterm(:) = val_exp
    var_name = 't2m_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_longterm(:) == val_exp)) t2m_longterm(:) = tref(:)
    !-
    ! the long-term reference temperature is a prognostic variable
    ! only in case the vegetation is dynamic
    !-
    IF (control%ok_dgvm) THEN
       tlong_ref(:) = val_exp
       var_name = 'tlong_ref'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
            &                .TRUE., tlong_ref, 'gather', nbp_glo, index_g)
       IF (ALL(tlong_ref(:) == val_exp)) tlong_ref(:) = tref(:)
    ENDIF
    !-
    t2m_month(:) = val_exp
    var_name = 't2m_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_month, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_month(:) == val_exp)) t2m_month(:) = tref(:)
    !-
    Tseason(:) = val_exp
    var_name = 'Tseason'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., Tseason, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason(:) == val_exp)) Tseason(:) = tref(:)
    !-
    Tseason_length(:) = val_exp
    var_name = 'Tseason_length'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., Tseason_length, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason_length(:) == val_exp)) Tseason_length(:) = tref(:)
    !-
    Tseason_tmp(:) = val_exp
    var_name = 'Tseason_tmp'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., Tseason_tmp, 'gather', nbp_glo, index_g)
    IF (ALL(Tseason_tmp(:) == val_exp)) Tseason_tmp(:) = tref(:)
    !-

    Tmin_spring(:,:) = 999
    var_name = 'Tmin_spring'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              .TRUE., Tmin_spring, 'gather', nbp_glo, index_g)
!    IF (ALL(Tmin_spring(:) == val_exp)) Tmin_spring(:) = tref(:)
    !-
    Tmin_spring_time(:,:) = val_exp
    var_name = 'Tmin_spring_time'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              .TRUE., Tmin_spring_time, 'gather', nbp_glo, index_g)
    IF (ALL(Tmin_spring_time(:,:) == val_exp)) Tmin_spring_time(:,:) = zero
    !-

    onset_date(:,:,:) = val_exp
    DO k=1,2
       WRITE(part_str,'(I2)') k
       IF ( k < 10 ) part_str(1:1) = '0'
       var_name = 'onset_date_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                .TRUE., onset_date(:,:,k), 'gather', nbp_glo, index_g)
       IF (ALL(onset_date(:,:,k) == val_exp)) onset_date(:,:,k) = zero
    ENDDO
    !-

    t2m_week(:) = val_exp
    var_name = 't2m_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., t2m_week, 'gather', nbp_glo, index_g)
    IF (ALL(t2m_week(:) == val_exp)) t2m_week(:) = tref(:)
    !-
    tsoil_month(:,:) = val_exp
    var_name = 'tsoil_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nbdl, 1, itime, &
         &              .TRUE., tsoil_month, 'gather', nbp_glo, index_g)

    IF (ALL(tsoil_month(:,:) == val_exp)) THEN
       DO l=1,nbdl
          tsoil_month(:,l) = tref(:)
       ENDDO
    ENDIF
    !-
    soilhum_month(:,:) = val_exp
    var_name = 'soilhum_month'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nbdl, 1, itime, &
         &              .TRUE., soilhum_month, 'gather', nbp_glo, index_g)
    IF (ALL(soilhum_month(:,:) == val_exp)) soilhum_month(:,:) = zero
    !-
    ! 6 fire probability
    !-
    fireindex(:,:) = val_exp
    var_name = 'fireindex'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              .TRUE., fireindex, 'gather', nbp_glo, index_g)
    IF (ALL(fireindex(:,:) == val_exp)) fireindex(:,:) = zero
    !-
    firelitter(:,:) = val_exp
    var_name = 'firelitter'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              .TRUE., firelitter, 'gather', nbp_glo, index_g)
    IF (ALL(firelitter(:,:) == val_exp)) firelitter(:,:) = zero
    !-
    ! 7 maximum and minimum moisture availabilities for tropic phenology
    !-
    maxmoiavail_lastyear(:,:) = val_exp
    var_name = 'maxmoistr_last'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxmoiavail_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxmoiavail_lastyear(:,:) == val_exp)) &
         &     maxmoiavail_lastyear(:,:) = zero
    !-
    maxmoiavail_thisyear(:,:) = val_exp
    var_name = 'maxmoistr_this'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxmoiavail_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxmoiavail_thisyear(:,:) == val_exp)) &
         &     maxmoiavail_thisyear(:,:) = zero
    !-
    minmoiavail_lastyear(:,:) = val_exp
    var_name = 'minmoistr_last'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., minmoiavail_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(minmoiavail_lastyear(:,:) == val_exp)) &
         &     minmoiavail_lastyear(:,:) = un
    !-
    minmoiavail_thisyear(:,:) = val_exp
    var_name = 'minmoistr_this'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., minmoiavail_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL( minmoiavail_thisyear(:,:) == val_exp)) &
         &     minmoiavail_thisyear(:,:) = un
    !-
    ! 8 maximum "weekly" GPP
    !-
    maxgppweek_lastyear(:,:) = val_exp
    var_name = 'maxgppweek_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxgppweek_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxgppweek_lastyear(:,:) == val_exp)) &
         &     maxgppweek_lastyear(:,:) = zero
    !-
    maxgppweek_thisyear(:,:) = val_exp
    var_name = 'maxgppweek_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxgppweek_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxgppweek_thisyear(:,:) == val_exp)) &
         &     maxgppweek_thisyear(:,:) = zero
    !-
    ! 9 annual GDD0
    !-
    gdd0_thisyear(:) = val_exp
    var_name = 'gdd0_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., gdd0_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(gdd0_thisyear(:) == val_exp)) gdd0_thisyear(:) = zero
    !-
    gdd0_lastyear(:) = val_exp
    var_name = 'gdd0_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., gdd0_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(gdd0_lastyear(:) == val_exp)) gdd0_lastyear(:) = gdd_crit_estab
    !-
    ! 10 annual precipitation
    !-
    precip_thisyear(:) = val_exp
    var_name = 'precip_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., precip_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(precip_thisyear(:) == val_exp)) precip_thisyear(:) = zero
    !-
    precip_lastyear(:) = val_exp
    var_name = 'precip_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., precip_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(precip_lastyear(:) == val_exp)) &
         &     precip_lastyear(:) = precip_crit
    !-
    ! 11 derived "biometeorological" variables
    !-
    gdd_m5_dormance(:,:) = val_exp
    var_name = 'gdd_m5_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_m5_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_m5_dormance(:,:) == val_exp)) &
         &     gdd_m5_dormance(:,:) = undef
    !-
    gdd_from_growthinit(:,:) = val_exp
    var_name = 'gdd_from_growthinit'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_from_growthinit, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_from_growthinit(:,:) == val_exp)) &
         &     gdd_from_growthinit(:,:) = zero
    !-
    gdd_midwinter(:,:) = val_exp
    var_name = 'gdd_midwinter'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gdd_midwinter, 'gather', nbp_glo, index_g)
    IF (ALL(gdd_midwinter(:,:) == val_exp)) gdd_midwinter(:,:) = undef
    !-
    ncd_dormance(:,:) = val_exp
    var_name = 'ncd_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., ncd_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(ncd_dormance(:,:) == val_exp)) ncd_dormance(:,:) = undef
    !-
    ngd_minus5(:,:) = val_exp
    var_name = 'ngd_minus5'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., ngd_minus5, 'gather', nbp_glo, index_g)
    IF (ALL(ngd_minus5(:,:) == val_exp)) ngd_minus5(:,:) = zero
    !-
    time_hum_min(:,:) = val_exp
    var_name = 'time_hum_min'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., time_hum_min, 'gather', nbp_glo, index_g)
    IF (ALL(time_hum_min(:,:) == val_exp)) time_hum_min(:,:) = undef
    !-
    hum_min_dormance(:,:) = val_exp
    var_name = 'hum_min_dormance'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., hum_min_dormance, 'gather', nbp_glo, index_g)
    IF (ALL(hum_min_dormance(:,:) == val_exp)) &
         &     hum_min_dormance(:,:) = undef
    !-
    ! 12 Plant status
    !-
    PFTpresent_real(:,:) = val_exp
    var_name = 'PFTpresent'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., PFTpresent_real, 'gather', nbp_glo, index_g)
    IF (ALL(PFTpresent_real(:,:) == val_exp)) PFTpresent_real(:,:) = zero
    WHERE (PFTpresent_real(:,:) >= .5)
       PFTpresent = .TRUE.
    ELSEWHERE
       PFTpresent = .FALSE.
    ENDWHERE
    !-
    ind(:,:) = val_exp
    var_name = 'ind'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., ind, 'gather', nbp_glo, index_g)
    IF (ALL(ind(:,:) == val_exp)) ind(:,:) = zero
    !-
    adapted(:,:) = val_exp
    var_name = 'adapted'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., adapted, 'gather', nbp_glo, index_g)
    IF (ALL(adapted(:,:) == val_exp)) adapted(:,:) = zero
    !-
    regenerate(:,:) = val_exp
    var_name = 'regenerate'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., regenerate, 'gather', nbp_glo, index_g)
    IF (ALL(regenerate(:,:) == val_exp)) regenerate(:,:) = zero
    !-
    npp_longterm(:,:) = val_exp
    var_name = 'npp_longterm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., npp_longterm, 'gather', nbp_glo, index_g)
    IF (ALL(npp_longterm(:,:) == val_exp)) npp_longterm(:,:) = zero
    !-
    !! Arsene 25-06-2014 NPPcumul
    npp0_cumul(:,:) = val_exp                                                  !! Arsene 25-06-2014 NPPcumul
    var_name = 'npp0_cumul'                                                    !! Arsene 25-06-2014 NPPcumul
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &     !! Arsene 25-06-2014 NPPcumul
         &              .TRUE., npp0_cumul, 'gather', nbp_glo, index_g)        !! Arsene 25-06-2014 NPPcumul
    IF (ALL(npp0_cumul(:,:) == val_exp)) npp0_cumul(:,:) = zero                !! Arsene 25-06-2014 NPPcumul
    !! Arsene 25-06-2014 NPPcumul
    !-
    lm_lastyearmax(:,:) = val_exp
    var_name = 'lm_lastyearmax'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., lm_lastyearmax, 'gather', nbp_glo, index_g)
    IF (ALL(lm_lastyearmax(:,:) == val_exp)) lm_lastyearmax(:,:) = zero
    !-
    lm_thisyearmax(:,:) = val_exp
    var_name = 'lm_thisyearmax'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., lm_thisyearmax, 'gather', nbp_glo, index_g)
    IF (ALL(lm_thisyearmax(:,:) == val_exp)) lm_thisyearmax(:,:) = zero
    !-
    maxfpc_lastyear(:,:) = val_exp
    var_name = 'maxfpc_lastyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxfpc_lastyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxfpc_lastyear(:,:) == val_exp)) maxfpc_lastyear(:,:) = zero
    !-
    maxfpc_thisyear(:,:) = val_exp
    var_name = 'maxfpc_thisyear'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., maxfpc_thisyear, 'gather', nbp_glo, index_g)
    IF (ALL(maxfpc_thisyear(:,:) == val_exp)) maxfpc_thisyear(:,:) = zero
    !-
    turnover_time(:,:) = val_exp
    var_name = 'turnover_time'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., turnover_time, 'gather', nbp_glo, index_g)
    IF ( ALL( turnover_time(:,:) == val_exp)) turnover_time(:,:) = 100.
    !-
    turnover_longterm(:,:,:,:) = val_exp
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'turnover_longterm_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &              .TRUE., turnover_longterm(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(turnover_longterm(:,:,k,l) == val_exp)) &
               &       turnover_longterm(:,:,k,l) = zero
       ENDDO
    END DO
    !-
    gpp_week(:,:) = val_exp
    var_name = 'gpp_week'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &              .TRUE., gpp_week, 'gather', nbp_glo, index_g)
    IF (ALL(gpp_week(:,:) == val_exp)) gpp_week(:,:) = zero
    !-
    biomass(:,:,:,:) = val_exp
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'biomass_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                   .TRUE., biomass(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(biomass(:,:,k,l) == val_exp)) biomass(:,:,k,l) = zero
       ENDDO
    END DO
    !-
    resp_maint_part(:,:,:) = val_exp
    DO k=1,nparts
       WRITE(part_str,'(I2)') k
       IF ( k < 10 ) part_str(1:1) = '0'
       var_name = 'maint_resp_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                   .TRUE., resp_maint_part(:,:,k), 'gather', nbp_glo, index_g)
       IF (ALL(resp_maint_part(:,:,k) == val_exp)) resp_maint_part(:,:,k) = zero
    ENDDO
    !-
    leaf_age(:,:,:) = val_exp
    DO m=1,nleafages
       WRITE (part_str,'(I2)') m
       IF ( m < 10 ) part_str(1:1) = '0'
       var_name = 'leaf_age_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                   .TRUE., leaf_age(:,:,m), 'gather', nbp_glo, index_g)
       IF (ALL(leaf_age(:,:,m) == val_exp)) leaf_age(:,:,m) = zero
    ENDDO
    !-
    leaf_frac(:,:,:) = val_exp
    DO m=1,nleafages
       WRITE(part_str,'(I2)') m
       IF ( m < 10 ) part_str(1:1) = '0'
       var_name = 'leaf_frac_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                  .TRUE., leaf_frac(:,:,m), 'gather', nbp_glo, index_g)
       IF (ALL(leaf_frac(:,:,m) == val_exp)) leaf_frac(:,:,m) = zero
    ENDDO
    !-
    senescence_real(:,:) = val_exp
    var_name = 'senescence'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., senescence_real, 'gather', nbp_glo, index_g)
    IF (ALL(senescence_real(:,:) == val_exp)) senescence_real(:,:) = zero
    WHERE ( senescence_real(:,:) >= .5 )
       senescence = .TRUE.
    ELSEWHERE
       senescence = .FALSE.
    ENDWHERE
    !-
    begin_leaves_real(:,:) = val_exp
    var_name = 'begin_leaves'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., begin_leaves_real, 'gather', nbp_glo, index_g)
    IF (ALL(begin_leaves_real(:,:) == val_exp)) begin_leaves_real(:,:) = zero
    WHERE ( begin_leaves_real(:,:) >= .5 )
       begin_leaves = .TRUE.
    ELSEWHERE
       begin_leaves = .FALSE.
    ENDWHERE
    !-
    when_growthinit(:,:) = val_exp
    var_name = 'when_growthinit'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., when_growthinit, 'gather', nbp_glo, index_g)
    IF (ALL(when_growthinit(:,:) == val_exp)) &
         &     when_growthinit(:,:) = zero
    !-
    age(:,:) = val_exp
    var_name = 'age'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., age, 'gather', nbp_glo, index_g)
    IF (ALL(age(:,:) == val_exp)) age(:,:) = zero
    !-
!! Arsene 27-08-2015 Add dia_cut
    dia_cut(:,:) = val_exp                                             !! Arsene 27-08-2015 Add
    var_name = 'dia_cut'                                               !! Arsene 27-08-2015 Add
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, & !! Arsene 27-08-2015 Add
         &                .TRUE., dia_cut, 'gather', nbp_glo, index_g) !! Arsene 27-08-2015 Add
    IF (ALL(dia_cut(:,:) == val_exp)) dia_cut(:,:) = zero           !! Arsene 27-08-2015 Add
!! Arsene 27-08-2015 Add dia_cut
    !-
    ! 13 CO2
    !-
    resp_hetero(:,:) = val_exp
    var_name = 'resp_hetero'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                .TRUE., resp_hetero, 'gather', nbp_glo, index_g)
    IF (ALL(resp_hetero(:,:) == val_exp)) resp_hetero(:,:) = zero
    !-
    resp_maint(:,:) = val_exp
    var_name = 'resp_maint'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., resp_maint, 'gather', nbp_glo, index_g)
    IF (ALL(resp_maint(:,:) == val_exp)) resp_maint(:,:) = zero
    !-
    resp_growth(:,:) = val_exp
    var_name = 'resp_growth'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., resp_growth, 'gather', nbp_glo, index_g)
    IF (ALL(resp_growth(:,:) == val_exp)) resp_growth(:,:) = zero
    !-
    co2_fire(:,:) = val_exp
    var_name = 'co2_fire'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &                .TRUE., co2_fire, 'gather', nbp_glo, index_g)
    IF (ALL(co2_fire(:,:) == val_exp)) co2_fire(:,:) = zero
    !-
    co2_to_bm_dgvm(:,:) = val_exp
    var_name = 'co2_to_bm_dgvm'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
         &                .TRUE., co2_to_bm_dgvm, 'gather', nbp_glo, index_g)
    IF (ALL(co2_to_bm_dgvm(:,:) == val_exp)) co2_to_bm_dgvm(:,:) = zero
    !-
    ! 14 vegetation distribution after last light competition
    !-
    veget_lastlight(:,:) = val_exp
    var_name = 'veget_lastlight'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., veget_lastlight, 'gather', nbp_glo, index_g)
    IF (ALL(veget_lastlight(:,:) == val_exp)) veget_lastlight(:,:) = zero
    !-
    ! 15 establishment criteria
    !-
    everywhere(:,:) = val_exp
    var_name = 'everywhere'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., everywhere, 'gather', nbp_glo, index_g)
    IF (ALL(everywhere(:,:) == val_exp)) everywhere(:,:) = zero
    !-
    need_adjacent_real(:,:) = val_exp
    var_name = 'need_adjacent'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., need_adjacent_real, 'gather', nbp_glo, index_g)
    IF (ALL(need_adjacent_real(:,:) == val_exp)) &
         &     need_adjacent_real(:,:) = zero
    WHERE ( need_adjacent_real(:,:) >= .5 )
       need_adjacent = .TRUE.
    ELSEWHERE
       need_adjacent = .FALSE.
    ENDWHERE
    !-
    RIP_time(:,:) = val_exp
    var_name = 'RIP_time'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
         &                .TRUE., RIP_time, 'gather', nbp_glo, index_g)
    IF (ALL(RIP_time(:,:) == val_exp)) RIP_time(:,:) = large_value
    !-
    ! 16 black carbon
    !-
    black_carbon(:) = val_exp
    var_name = 'black_carbon'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &                .TRUE., black_carbon, 'gather', nbp_glo, index_g)
    IF (ALL(black_carbon(:) == val_exp)) black_carbon(:) = zero
    !-
    ! 17 litter
    !-
    litterpart(:,:,:) = val_exp
    DO l=1,nlitt
       var_name = 'litterpart_'//litter_str(l)
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                   .TRUE., litterpart(:,:,l), 'gather', nbp_glo, index_g)
       IF (ALL(litterpart(:,:,l) == val_exp)) litterpart(:,:,l) = zero
    ENDDO
    !-
    litter(:,:,:,:,:) = val_exp
    DO k = 1,nelements
       DO l = 1,nlevs
          DO m = 1,nvm
             WRITE (part_str, '(I2)') m
             IF (m<10) part_str(1:1)='0'
             var_name = 'litter_'//part_str(1:LEN_TRIM(part_str))//'_'//level_str(l)//element_str(k)
             CALL restget_p (rest_id_stomate, var_name, nbp_glo, nlitt , 1, itime, &
                  &                     .TRUE., litter(:,:,m,l,k), 'gather', nbp_glo, index_g)
             IF (ALL(litter(:,:,m,l,k) == val_exp)) litter(:,:,m,l,k) = zero
          ENDDO
       ENDDO
    END DO
    !-
    dead_leaves(:,:,:) = val_exp
    DO l=1,nlitt
       var_name = 'dead_leaves_'//litter_str(l)
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                   .TRUE., dead_leaves(:,:,l), 'gather', nbp_glo, index_g)
       IF (ALL(dead_leaves(:,:,l) == val_exp)) dead_leaves(:,:,l) = zero
    ENDDO
    !-
    carbon(:,:,:) = val_exp
    DO m=1,nvm
       WRITE (part_str, '(I2)') m
       IF (m<10) part_str(1:1)='0'
       var_name = 'carbon_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, ncarb , 1, itime, &
            &                   .TRUE., carbon(:,:,m), 'gather', nbp_glo, index_g)
       IF (ALL(carbon(:,:,m) == val_exp)) carbon(:,:,m) = zero
    ENDDO
    !-
    lignin_struc(:,:,:) = val_exp
    DO l=1,nlevs
       var_name = 'lignin_struc_'//level_str(l)
       CALL restget_p &
            &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &     .TRUE., lignin_struc(:,:,l), 'gather', nbp_glo, index_g)
       IF (ALL(lignin_struc(:,:,l) == val_exp)) lignin_struc(:,:,l) = zero
    ENDDO
    !-
    ! 18 land cover change
    !-
    prod10(:,:) = val_exp
    var_name = 'prod10'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 11     , 1, itime, &
         &                .TRUE., prod10, 'gather', nbp_glo, index_g)
    IF (ALL(prod10(:,:) == val_exp)) prod10(:,:) = zero

    prod100(:,:) = val_exp
    var_name = 'prod100'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 101     , 1, itime, &
         &                .TRUE., prod100, 'gather', nbp_glo, index_g)
    IF (ALL(prod100(:,:) == val_exp)) prod100(:,:) = zero


    flux10(:,:) = val_exp
    var_name = 'flux10'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 10     , 1, itime, &
         &                .TRUE., flux10, 'gather', nbp_glo, index_g)
    IF (ALL(flux10(:,:) == val_exp)) flux10(:,:) = zero

    flux100(:,:) = val_exp
    var_name = 'flux100'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 100     , 1, itime, &
         &                .TRUE., flux100, 'gather', nbp_glo, index_g)
    IF (ALL(flux100(:,:) == val_exp)) flux100(:,:) = zero

    convflux(:) = val_exp
    var_name = 'convflux'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., convflux, 'gather', nbp_glo, index_g)
    IF (ALL(convflux(:) == val_exp)) convflux(:) = zero

    cflux_prod10(:) = val_exp
    var_name = 'cflux_prod10'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., cflux_prod10, 'gather', nbp_glo, index_g)
    IF (ALL(cflux_prod10(:) == val_exp)) cflux_prod10(:) = zero

    cflux_prod100(:) = val_exp
    var_name = 'cflux_prod100'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., cflux_prod100, 'gather', nbp_glo, index_g)
    IF (ALL(cflux_prod100(:) == val_exp)) cflux_prod100(:) = zero

    bm_to_litter(:,:,:,:) = val_exp
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF ( k < 10 ) part_str(1:1) = '0'
          var_name = 'bm_to_litter_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
               &                .TRUE., bm_to_litter(:,:,k,l), 'gather', nbp_glo, index_g)
          IF (ALL(bm_to_litter(:,:,k,l) == val_exp)) bm_to_litter(:,:,k,l) = zero
       ENDDO
    END DO

    carb_mass_total(:) = val_exp
    var_name = 'carb_mass_total'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
         &              .TRUE., carb_mass_total, 'gather', nbp_glo, index_g)
    IF (ALL(carb_mass_total(:) == val_exp)) carb_mass_total(:) = zero

    !Permafrost carbon related
    deepC_a(:,:,:) = val_exp
    DO m=1,nvm
       WRITE(part_str,'(I2)') m
       IF ( m < 10 ) part_str(1:1) = '0'
       var_name = 'deepC_a_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, ndeep, 1, itime, &
            &               .TRUE., deepC_a(:,:,m), 'gather', nbp_glo, index_g)
    END DO
    IF (ALL(deepC_a(:,:,:) == val_exp)) THEN
       deepC_a(:,:,:) = zero !deepC_a_init
    ENDIF

    deepC_s(:,:,:) = val_exp
    DO m=1,nvm
      WRITE(part_str,'(I2)') m
      IF ( m < 10 ) part_str(1:1) = '0'
      var_name = 'deepC_s_'//part_str(1:LEN_TRIM(part_str))
      CALL restget_p (rest_id_stomate,var_name, nbp_glo, ndeep, 1, itime, &
           &               .TRUE., deepC_s(:,:,m), 'gather', nbp_glo, index_g)
    END DO
    IF (ALL(deepC_s(:,:,:) == val_exp)) THEN
        deepC_s(:,:,:) = zero !deepC_s_init
    ENDIF

    deepC_p(:,:,:) = val_exp
    DO m=1,nvm
      WRITE(part_str,'(I2)') m
      IF ( m < 10 ) part_str(1:1) = '0'
      var_name = 'deepC_p_'//part_str(1:LEN_TRIM(part_str))
      CALL restget_p (rest_id_stomate,var_name, nbp_glo, ndeep, 1, itime, &
           &               .TRUE., deepC_p(:,:,m), 'gather', nbp_glo, index_g)
    END DO
    IF (ALL(deepC_p(:,:,:) == val_exp)) THEN
        deepC_p(:,:,:) = zero !deepC_p_init
    ENDIF
    O2_soil(:,:,:) = val_exp
    DO m=1,nvm
      WRITE(part_str,'(I2)') m
      IF ( m < 10 ) part_str(1:1) = '0'
      var_name = 'O2_soil_'//part_str(1:LEN_TRIM(part_str))
      CALL restget_p (rest_id_stomate,var_name, nbp_glo, ndeep, 1, itime, &
           &               .TRUE., O2_soil(:,:,m), 'gather', nbp_glo, index_g)
    END DO
    IF (ALL(O2_soil(:,:,:) == val_exp)) THEN
       O2_soil(:,:,:) = O2_init_conc
    ENDIF

    CH4_soil(:,:,:) = val_exp
    DO m=1,nvm
       WRITE(part_str,'(I2)') m
       IF ( m < 10 ) part_str(1:1) = '0'
       var_name = 'CH4_soil_'//part_str(1:LEN_TRIM(part_str))
       CALL restget_p (rest_id_stomate,var_name, nbp_glo, ndeep, 1, itime, &
            &               .TRUE., CH4_soil(:,:,m), 'gather', nbp_glo, index_g)
    END DO
    IF (ALL(CH4_soil(:,:,:) == val_exp)) THEN
       CH4_soil(:,:,:) = CH4_init_conc
    ENDIF

    O2_snow(:,:,:) = val_exp
    DO m=1,nvm
      WRITE(part_str,'(I2)') m
      IF ( m < 10 ) part_str(1:1) = '0'
      var_name = 'O2_snow_'//part_str(1:LEN_TRIM(part_str))
      CALL restget_p (rest_id_stomate,var_name, nbp_glo, nsnow, 1, itime, &
           &               .TRUE., O2_snow(:,:,m), 'gather', nbp_glo, index_g)
   END DO
   IF (ALL(O2_snow(:,:,:) == val_exp)) THEN
      O2_snow(:,:,:) = O2_init_conc
   ENDIF

   CH4_snow(:,:,:) = val_exp
   DO m=1,nvm
     WRITE(part_str,'(I2)') m
     IF ( m < 10 ) part_str(1:1) = '0'
     var_name = 'CH4_snow_'//part_str(1:LEN_TRIM(part_str))
     CALL restget_p (rest_id_stomate,var_name, nbp_glo, nsnow, 1, itime, &
          &               .TRUE., CH4_snow(:,:,m), 'gather', nbp_glo, index_g)
  END DO
  IF (ALL(CH4_snow(:,:,:) == val_exp)) THEN
     CH4_snow(:,:,:) = CH4_init_conc
  ENDIF
  thawed_humidity(:) = val_exp
  var_name = 'thawed_humidity'

  CALL getin('reset_thawed_humidity', reset_thawed_humidity)
  if ( reset_thawed_humidity ) then
     CALL getin('thawed_humidity_input', thawed_humidity_input)
     thawed_humidity(:) = thawed_humidity_input
  else
     CALL restget_p (rest_id_stomate,var_name, nbp_glo, 1, 1, itime, &
          &               .TRUE., thawed_humidity, 'gather', nbp_glo, index_g)
     IF (ALL(thawed_humidity(:) == val_exp)) THEN
        thawed_humidity(:) = thawed_humidity_input
        read_input_thawed_humidity = .TRUE.
     ENDIF
  endif


  depth_organic_soil(:) = val_exp
  var_name = 'depth_organic_soil'
  CALL restget_p (rest_id_stomate,var_name, nbp_glo, 1, 1, itime, &
       &               .TRUE., depth_organic_soil, 'gather', nbp_glo, index_g)
  IF (ALL(depth_organic_soil(:) == val_exp)) THEN
     depth_organic_soil(:) = 0.0
     read_input_depth_organic_soil = .TRUE.
  ENDIF

  altmax(:,:) = val_exp
  var_name = 'altmax'
  CALL restget_p (rest_id_stomate,var_name, nbp_glo, nvm, 1, itime, &
       &               .TRUE., altmax, 'gather', nbp_glo, index_g)
  IF (ALL(altmax(:,:) == val_exp)) THEN
     altmax(:,:) = 0.0
  ENDIF

  fixed_cryoturbation_depth(:,:) = val_exp
  var_name = 'fixed_cryoturb_depth'
  CALL restget_p (rest_id_stomate,var_name, nbp_glo, nvm, 1, itime, &
       &               .TRUE., fixed_cryoturbation_depth, 'gather', nbp_glo, index_g)
  IF (ALL(fixed_cryoturbation_depth(:,:) == val_exp)) THEN
     fixed_cryoturbation_depth(:,:) = 0.0
  ENDIF

  !-

!!!!! Wetland CH4 methane
    !pss:+
    uo_0(:,:) = val_exp
    var_name = 'uo_0'
    CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
         &              .TRUE.,uo_0 , 'gather', nbp_glo, index_g)
    IF (ALL(uo_0(:,:) == val_exp)) THEN
       DO nivo=1,nvert
          IF (nivo .LE. ns) THEN
            uo_0(:,nivo) = scmax
         ELSE
            uo_0(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uold2_0(:,:) = val_exp
   var_name = 'uold2_0'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_0 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_0(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns) THEN
            uold2_0(:,nivo) = scmax
         ELSE
            uold2_0(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
   
   uo_wet1(:,:) = val_exp
   var_name = 'uo_wet1'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
       &              .TRUE.,uo_wet1 , 'gather', nbp_glo, index_g)
   IF (ALL(uo_wet1(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uo_wet1(:,nivo) = scmax
         ELSE
            uo_wet1(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
   
   uold2_wet1(:,:) = val_exp
   var_name = 'uold2_wet1'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_wet1 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_wet1(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uold2_wet1(:,nivo) = scmax
         ELSE
            uold2_wet1(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
  
   uo_wet2(:,:) = val_exp
   var_name = 'uo_wet2'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uo_wet2 , 'gather', nbp_glo, index_g)
   IF (ALL(uo_wet2(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uo_wet2(:,nivo) = scmax
         ELSE
            uo_wet2(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uold2_wet2(:,:) = val_exp
   var_name = 'uold2_wet2'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_wet2 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_wet2(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uold2_wet2(:,nivo) = scmax
         ELSE
            uold2_wet2(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uo_wet3(:,:) = val_exp
   var_name = 'uo_wet3'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uo_wet3 , 'gather', nbp_glo, index_g)
   IF (ALL(uo_wet3(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uo_wet3(:,nivo) = scmax
         ELSE
            uo_wet3(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
   
   uold2_wet3(:,:) = val_exp
   var_name = 'uold2_wet3'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_wet3 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_wet3(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uold2_wet3(:,nivo) = scmax
         ELSE
            uold2_wet3(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uo_wet4(:,:) = val_exp
   var_name = 'uo_wet4'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uo_wet4 , 'gather', nbp_glo, index_g)
   IF (ALL(uo_wet4(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uo_wet4(:,nivo) = scmax
         ELSE
            uo_wet4(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF

   uold2_wet4(:,:) = val_exp
   var_name = 'uold2_wet4'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo,   nvert, 1, itime, &
        &              .TRUE.,uold2_wet4 , 'gather', nbp_glo, index_g)
   IF (ALL(uold2_wet4(:,:) == val_exp)) THEN
      DO nivo=1,nvert
         IF (nivo .LE. ns-10) THEN
            uold2_wet4(:,nivo) = scmax
         ELSE
            uold2_wet4(:,nivo) = CH4atmo_CONC
         ENDIF
      ENDDO
   ENDIF
   
   tsurf_year(:) = val_exp
   var_name = 'tsurf_year'
   CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
        &              .TRUE., tsurf_year, 'gather', nbp_glo, index_g)
   IF (ALL(tsurf_year(:) == val_exp)) tsurf_year(:) = tref(:)
   
!pss:-
!JCADD
!-
!  sla_calc(:,:) = val_exp
!  var_name = 'sla_calc'
!  CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm   , 1, itime, &
! &              .TRUE., sla_calc, 'gather', nbp_glo, index_g)
!  IF (ALL(sla_calc(:,:) == val_exp)) THEN
!     DO j=2,nvm
!        sla_calc(:,j) = sla(j)
!     END DO
!  END IF
  wshtotsum(:,:) = val_exp
  var_name = 'wshtotsum'
  CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
 &              .TRUE., wshtotsum, 'gather', nbp_glo, index_g)
    IF (ALL(wshtotsum(:,:) == val_exp)) wshtotsum(:,:) = zero
!-
  sr_ugb(:,:) = val_exp
  var_name = 'sr_ugb'
  CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
 &              .TRUE., sr_ugb, 'gather', nbp_glo, index_g)
    IF (ALL(sr_ugb(:,:) == val_exp)) sr_ugb(:,:) = zero
!-
  sla_calc(:,:) = val_exp
  var_name = 'sla_calc'
  CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
 &              .TRUE., sla_calc, 'gather', nbp_glo, index_g)
!    IF (ALL(sla_calc(:,:) == val_exp)) sla_calc(:,:) = zero
  IF (ALL(sla_calc(:,:) == val_exp)) THEN
     DO j=1,nvm
        sla_calc(:,j) = sla(j)
     END DO
  END IF
!-
  nb_ani(:,:) = val_exp
  var_name = 'nb_ani'
  CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
 &              .TRUE., nb_ani, 'gather', nbp_glo, index_g)
    IF (ALL(nb_ani(:,:) == val_exp)) nb_ani(:,:) = zero
!-
  grazed_frac(:,:) = val_exp
  var_name = 'grazed_frac'
  CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
 &              .TRUE., grazed_frac, 'gather', nbp_glo, index_g)
    IF (ALL(grazed_frac(:,:) == val_exp)) grazed_frac(:,:) = zero
!-
  import_yield(:,:) = val_exp
  var_name = 'import_yield'
  CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm     , 1, itime, &
 &              .TRUE., import_yield, 'gather', nbp_glo, index_g)
    IF (ALL(import_yield(:,:) == val_exp)) import_yield(:,:) = zero

!-
   t2m_14(:) = val_exp
  var_name = 't2m_14'
  CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
 &              .TRUE., t2m_14, 'gather', nbp_glo, index_g)
  IF (ALL(t2m_14(:) == val_exp)) t2m_14(:) = tref(:)
!

    litter_not_avail(:,:,:) = val_exp
    DO l=1,nlitt
       var_name = 'litter_not_avail_'//litter_str(l)
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm  , 1, itime, &
            &                   .TRUE., litter_not_avail(:,l,:), 'gather', nbp_glo, index_g)
       IF (ALL(litter_not_avail(:,l,:) == val_exp)) litter_not_avail(:,l,:) = zero
    ENDDO

!ENDJCADD
 
    !-
    ! 19. Spinup
    !-
!    IF (spinup_analytic) THEN

       IF (is_root_prc) THEN
          temp_global_years(1) = val_exp
          var_name = 'Global_years'
          CALL restget (rest_id_stomate, var_name, 1 ,1  , 1, itime, &
               &                .TRUE., temp_global_years)
          IF(temp_global_years(1) == val_exp) temp_global_years(1) = zero
          global_years = INT(temp_global_years(1))
       ENDIF
       CALL bcast(global_years)

       nbp_accu(:) = val_exp
       var_name = 'nbp_sum'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
            &              .TRUE., nbp_accu, 'gather', nbp_glo, index_g)
       IF (ALL(nbp_accu(:) == val_exp)) nbp_accu(:) = zero    

       nbp_flux(:) = val_exp
       var_name = 'nbp_flux'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo, 1     , 1, itime, &
            &              .TRUE., nbp_flux, 'gather', nbp_glo, index_g)
       IF (ALL(nbp_flux(:) == val_exp)) nbp_flux(:) = zero     

       !-
       ok_equilibrium_real(:) = val_exp
       var_name = 'ok_equilibrium'
       CALL restget_p (rest_id_stomate, var_name, nbp_glo , 1  , 1, itime, &
            &                .TRUE., ok_equilibrium_real,'gather', nbp_glo, index_g)
       IF (ALL(ok_equilibrium_real(:) == val_exp)) ok_equilibrium_real(:) = zero
       WHERE(ok_equilibrium_real(:) >= 0.5) 
          ok_equilibrium = .TRUE.
       ELSEWHERE
          ok_equilibrium = .FALSE.
       ENDWHERE

       MatrixV(:,:,:,:) = val_exp
       DO k = 1,nbpools
          DO j = 1,nbpools
             WRITE(part_str,'(I2)') k
             IF (k < 10) part_str(1:1) = '0'             
             var_name = 'MatrixV_'//part_str(1:LEN_TRIM(part_str))//'_'//TRIM(pools_str(j))
             CALL restget_p (rest_id_stomate, var_name, nbp_glo, nvm , 1, itime, &
                  &                     .TRUE., MatrixV(:,:,k,j), 'gather', nbp_glo, index_g)
          ENDDO
       ENDDO
       ! If nothing is found in the restart file, we initialize each submatrix by identity
       IF (ALL(MatrixV(:,:,:,:) == val_exp))  THEN 
          MatrixV(:,:,:,:) = zero
          DO l = 1,nbpools
             MatrixV(:,:,l,l) = un
          END DO
       END IF

       VectorU(:,:,:)  = val_exp
       DO k= 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'Vector_U_'//part_str(1:LEN_TRIM(part_str))
          CALL restget_p &
               &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     .TRUE., VectorU(:,:,k), 'gather', nbp_glo, index_g)
          IF (ALL(VectorU(:,:,k) == val_exp))  VectorU(:,:,k) = zero
       ENDDO
       
       previous_stock(:,:,:)  = val_exp
       DO k= 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'previous_stock_'//part_str(1:LEN_TRIM(part_str))
          CALL restget_p &
               &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     .TRUE., previous_stock(:,:,k), 'gather', nbp_glo, index_g)
          IF (ALL(previous_stock(:,:,k) == val_exp))  previous_stock(:,:,k) = undef_sechiba
       ENDDO
       
       current_stock(:,:,:)  = val_exp
       DO k= 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'current_stock_'//part_str(1:LEN_TRIM(part_str))
          CALL restget_p &
               &    (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &     .TRUE., current_stock(:,:,k), 'gather', nbp_glo, index_g)
          IF (ALL(current_stock(:,:,k) == val_exp))  current_stock(:,:,k) = zero
       ENDDO
          
!    ENDIF ! spinup_matrix_method
    IF (bavard >= 4) WRITE(numout,*) 'Leaving readstart'
    !-----------------------
  END SUBROUTINE readstart
  !-
  !===
  !-
  SUBROUTINE writerestart &
       & (npts, index, day_counter, dt_days, date, &
       &  ind, adapted, regenerate, moiavail_daily, gdd_init_date, litterhum_daily, &
       &  t2m_daily, t2m_min_daily, tsurf_daily, tsoil_daily, &
       &  soilhum_daily, precip_daily, gpp_daily, npp_daily, &
       &  turnover_daily, moiavail_month, moiavail_week, &
       &  t2m_longterm, tlong_ref, t2m_month, t2m_week, &
       &  tsoil_month, soilhum_month, fireindex, firelitter, &
       &  maxmoiavail_lastyear, maxmoiavail_thisyear, &
       &  minmoiavail_lastyear, minmoiavail_thisyear, &
       &  maxgppweek_lastyear, maxgppweek_thisyear, &
       &  gdd0_lastyear, gdd0_thisyear, precip_lastyear, precip_thisyear, &
       &  gdd_m5_dormance, gdd_from_growthinit, gdd_midwinter, ncd_dormance, ngd_minus5, &
       &  PFTpresent, npp_longterm, lm_lastyearmax, lm_thisyearmax, &
       &  maxfpc_lastyear, maxfpc_thisyear, &
       &  turnover_longterm, gpp_week, biomass, resp_maint_part, &
       &  leaf_age, leaf_frac, senescence, when_growthinit, age, &
       &  resp_hetero, resp_maint, resp_growth, co2_fire, co2_to_bm_dgvm, &
       &  veget_lastlight, everywhere, need_adjacent, RIP_time, &
       &  time_hum_min, hum_min_dormance, &
       &  litterpart, litter, dead_leaves, &
       &  carbon, black_carbon, lignin_struc, turnover_time, &
       &  prod10,prod100 ,flux10, flux100, &
       &  convflux, cflux_prod10, cflux_prod100, bm_to_litter, carb_mass_total, &
       &  Tseason, Tseason_length, Tseason_tmp, &
       &  Tmin_spring, Tmin_spring_time, begin_leaves, onset_date, &
       &  global_years, ok_equilibrium, nbp_accu, nbp_flux, &
       &  MatrixV, VectorU, previous_stock, current_stock,&
       &  npp0_cumul,snowtemp_min,snowdz_min, dia_cut, & !! Arsene 25-06-2014 NPPcumul ADD  !! Arsene 19-08-2014 Add snowtemp and snowdz !! Arsene 27-08-2015 add dia_cut
       &  deepC_a, deepC_s, deepC_p, O2_soil, CH4_soil, O2_snow, CH4_snow, &
       &  thawed_humidity, depth_organic_soil, altmax, fixed_cryoturbation_depth, & !pss+
       &  uo_0, uold2_0, uo_wet1, uold2_wet1, uo_wet2, uold2_wet2, uo_wet3, uold2_wet3, &
       &  uo_wet4, uold2_wet4, tsurf_year, &!) !pss:-
!JCADD
       &  wshtotsum, sr_ugb, sla_calc, nb_ani, grazed_frac, &
       &  import_yield, t2m_14, litter_not_avail)
!ENDJCADD
    !---------------------------------------------------------------------
    !- write restart file
    !---------------------------------------------------------------------
    !-
    ! 0 declarations
    !-
    ! 0.1 input
    !-
    ! Domain size
    INTEGER(i_std),INTENT(in) :: npts
    ! Indices of the points on the map
    INTEGER(i_std),DIMENSION(npts),INTENT(in) :: index
    ! counts time until next STOMATE time step
    REAL(r_std),INTENT(in) :: day_counter
    ! time step of STOMATE in days
    REAL(r_std),INTENT(in) :: dt_days
    ! date (d)
    INTEGER(i_std),INTENT(in) :: date
    ! density of individuals (1/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: ind
    ! Winter too cold? between 0 and 1
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: adapted
    ! Winter sufficiently cold? between 0 and 1
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: regenerate
    ! daily moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: moiavail_daily
    ! gdd init date
    REAL(r_std),DIMENSION(npts,2),INTENT(in) :: gdd_init_date
    ! daily litter humidity
    REAL(r_std),DIMENSION(npts),INTENT(in) :: litterhum_daily
    ! daily 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: t2m_daily
    ! daily minimum 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: t2m_min_daily
    ! daily surface temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: tsurf_daily
    ! daily soil temperatures (K)
    REAL(r_std),DIMENSION(npts,nbdl),INTENT(in) :: tsoil_daily
    ! daily soil humidity
    REAL(r_std),DIMENSION(npts,nbdl),INTENT(in) :: soilhum_daily
    ! daily precipitations (mm/day) (for phenology)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: precip_daily
    ! daily gross primary productivity (gC/m**2/day)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: gpp_daily
    ! daily net primary productivity (gC/m**2/day)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: npp_daily
    ! daily turnover rates (gC/m**2/day)
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: turnover_daily
    ! "monthly" moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: moiavail_month
    ! "weekly" moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: moiavail_week
    ! "long term" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: t2m_longterm
    ! "long term" reference 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: tlong_ref
    ! "monthly" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: t2m_month
    ! "seasonal" 2 meter temperatures (K) 
    REAL(r_std),DIMENSION(npts),INTENT(in)      :: Tseason
    ! temporary variable to calculate Tseason
    REAL(r_std),DIMENSION(npts),INTENT(in)      :: Tseason_length
    ! temporary variable to calculate Tseason
    REAL(r_std),DIMENSION(npts),INTENT(in)      :: Tseason_tmp

    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)  :: Tmin_spring
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)  :: Tmin_spring_time
    REAL(r_std),DIMENSION(npts,nvm,2),INTENT(in)  :: onset_date
    LOGICAL,DIMENSION(npts,nvm),INTENT(in)      :: begin_leaves

    ! "weekly" 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: t2m_week
    ! "monthly" soil temperatures (K)
    REAL(r_std),DIMENSION(npts,nbdl),INTENT(in) :: tsoil_month
    ! "monthly" soil humidity
    REAL(r_std),DIMENSION(npts,nbdl),INTENT(in) :: soilhum_month
    ! Probability of fire
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: fireindex
    ! Longer term total litter above the ground, gC/m**2 of ground
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: firelitter
    ! last year's maximum moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: maxmoiavail_lastyear
    ! this year's maximum moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: maxmoiavail_thisyear
    ! last year's minimum moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: minmoiavail_lastyear
    ! this year's minimum moisture availability
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: minmoiavail_thisyear
    ! last year's maximum weekly GPP
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: maxgppweek_lastyear
    ! this year's maximum weekly GPP
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: maxgppweek_thisyear
    ! last year's annual GDD0
    REAL(r_std),DIMENSION(npts),INTENT(in) :: gdd0_lastyear
    ! this year's annual GDD0
    REAL(r_std),DIMENSION(npts),INTENT(in) :: gdd0_thisyear
    ! last year's annual precipitation (mm/year)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: precip_lastyear
    ! this year's annual precipitation (mm/year)
    REAL(r_std),DIMENSION(npts),INTENT(in) :: precip_thisyear
    ! growing degree days, threshold -5 deg C (for phenology)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: gdd_m5_dormance
    ! growing degree days, from begin of season (crops)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: gdd_from_growthinit
    ! growing degree days since midwinter (for phenology)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: gdd_midwinter
    ! number of chilling days since leaves were lost (for phenology)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: ncd_dormance
    ! number of growing days, threshold -5 deg C (for phenology)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: ngd_minus5
    ! PFT exists (equivalent to fpc_max > 0 for natural PFTs)
    LOGICAL,DIMENSION(npts,nvm),INTENT(in) :: PFTpresent
    ! "long term" net primary productivity (gC/m**2/year)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: npp_longterm

    ! "counter" null or negative npp (days)                           !! Arsene 25-06-2014 NPPcumul
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)     :: npp0_cumul      !! Arsene 25-06-2014 NPPcumul
    ! Min daily snow layer temperature (K)                            !! Arsene 19-08-2014 Add
    REAL(r_std), DIMENSION(npts,nsnow), INTENT(in) :: snowtemp_min    !! Arsene 19-08-2014 Add
    ! Min daily snow layer thicknesse (cm?)                           !! Arsene 19-08-2014 Add
    REAL(r_std), DIMENSION(npts,nsnow), INTENT(in) :: snowdz_min      !! Arsene 19-08-2014 Add
    ! Fix diameter of vegetation (for shrub) after loss biomass (above snow) !! Arsene 27-08-2015 add dia_cut
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)   :: dia_cut               !! Arsene 27-08-2015 add dia_cut

    ! last year's maximum leaf mass, for each PFT (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: lm_lastyearmax
    ! this year's maximum leaf mass, for each PFT (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: lm_thisyearmax
    ! last year's maximum fpc for each natural PFT, on ground
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: maxfpc_lastyear
    ! this year's maximum fpc for each PFT,
    ! on *total* ground (see stomate_season)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: maxfpc_thisyear
    ! "long term" turnover rate (gC/m**2/year)
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: turnover_longterm
    ! "weekly" GPP (gC/day/(m**2 covered)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: gpp_week
    ! biomass (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: biomass
    ! maintenance respiration (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm,nparts),INTENT(in) :: resp_maint_part
    ! leaf age (days)
    REAL(r_std),DIMENSION(npts,nvm,nleafages),INTENT(in) :: leaf_age
    ! fraction of leaves in leaf age class
    REAL(r_std),DIMENSION(npts,nvm,nleafages),INTENT(in) :: leaf_frac
    ! is the plant senescent ?
    ! (only for deciduous trees - carbohydrate reserve)
    LOGICAL,DIMENSION(npts,nvm),INTENT(in) :: senescence
    ! how many days ago was the beginning of the growing season
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: when_growthinit
    ! mean age (years)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: age
    ! heterotrophic respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: resp_hetero
    ! maintenance respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: resp_maint
    ! growth respiration (gC/day/m**2)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: resp_growth
    ! carbon emitted into the atmosphere by fire (living and dead biomass)
    ! (in gC/m**2/time step)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: co2_fire
    ! biomass uptaken (gC/(m**2 of total ground)/day)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: co2_to_bm_dgvm
    ! vegetation fractions (on ground) after last light competition
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: veget_lastlight
    ! is the PFT everywhere in the grid box or very localized
    ! (after its introduction)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: everywhere
    ! in order for this PFT to be introduced,
    ! does it have to be present in an adjacent grid box?
    LOGICAL,DIMENSION(npts,nvm),INTENT(in) :: need_adjacent
    ! How much time ago was the PFT eliminated for the last time (y)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: RIP_time
    ! time elapsed since strongest moisture availability (d)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: time_hum_min
    ! minimum moisture during dormance
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: hum_min_dormance
    ! fraction of litter above the ground belonging to different PFTs
    REAL(r_std),DIMENSION(npts,nvm,nlitt),INTENT(in) :: litterpart
    ! metabolic and structural litter, above and below ground (gC/m**2)
    REAL(r_std),DIMENSION(npts,nlitt,nvm,nlevs,nelements),INTENT(in) :: litter
    ! dead leaves on ground, per PFT, metabolic and structural,
    ! in gC/(m**2 of ground)
    REAL(r_std),DIMENSION(npts,nvm,nlitt),INTENT(in) :: dead_leaves
    ! carbon pool: active, slow, or passive, (gC/m**2)
    REAL(r_std),DIMENSION(npts,ncarb,nvm),INTENT(in) :: carbon
    ! black carbon on the ground (gC/(m**2 of total ground))
    REAL(r_std),DIMENSION(npts),INTENT(in) :: black_carbon
    ! ratio Lignine/Carbon in structural litter, above and below ground, (gC/m**2)
    REAL(r_std),DIMENSION(npts,nvm,nlevs),INTENT(in) :: lignin_struc
    ! turnover_time of leaves
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in) :: turnover_time

    ! For Spinup matrix resolution
    INTEGER(i_std), INTENT(in) :: global_years   
    LOGICAL, DIMENSION(npts), INTENT(in) :: ok_equilibrium
    REAL(r_std), DIMENSION(npts), INTENT(in) :: nbp_accu  !! Accumulated Net Biospheric Production over the year 
    REAL(r_std), DIMENSION(npts), INTENT(in) :: nbp_flux  !! Net Biospheric Production over the year 
    !-
    REAL(r_std), DIMENSION(npts,nvm,nbpools,nbpools), INTENT(in) :: MatrixV
    REAL(r_std), DIMENSION(npts,nvm,nbpools), INTENT(in) :: VectorU
    REAL(r_std), DIMENSION(npts,nvm,nbpools), INTENT(in) :: previous_stock
    REAL(r_std), DIMENSION(npts,nvm,nbpools), INTENT(in) :: current_stock 
    ! Permafrost carbon related
    real(r_std), DIMENSION(npts,ndeep,nvm),INTENT(inout) :: deepC_a
    real(r_std), DIMENSION(npts,ndeep,nvm),intent(inout) :: deepC_s
    real(r_std), DIMENSION(npts,ndeep,nvm),intent(inout) :: deepC_p
    real(r_std), DIMENSION(npts,ndeep,nvm),intent(inout) :: O2_soil
    real(r_std), DIMENSION(npts,ndeep,nvm),intent(inout) :: CH4_soil
    real(r_std), DIMENSION(npts,nsnow,nvm),intent(inout) :: O2_snow
    real(r_std), DIMENSION(npts,nsnow,nvm),intent(inout) :: CH4_snow
    real(r_std), DIMENSION(npts),intent(inout)           :: thawed_humidity
    real(r_std), DIMENSION(npts),intent(inout)           :: depth_organic_soil
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)      :: altmax
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)      :: fixed_cryoturbation_depth

    ! Wetland CH4 methane
    !pss:+
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in) :: uo_0
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in)  :: uold2_0
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in) :: uo_wet1
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in)  :: uold2_wet1
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in) :: uo_wet2
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in)  :: uold2_wet2
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in) :: uo_wet3
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in)  :: uold2_wet3
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in) :: uo_wet4
    REAL(r_std),DIMENSION(npts,nvert),INTENT(in)  :: uold2_wet4

    REAL(r_std),DIMENSION(npts),INTENT(in) :: tsurf_year
    !pss:-
!JCADD
  REAL(r_std),DIMENSION(npts,nvm),INTENT(in)    :: sla_calc
  REAL(r_std),DIMENSION(npts,nvm),INTENT(in)    :: wshtotsum
  REAL(r_std), DIMENSION(npts,nvm), INTENT(in)  ::  sr_ugb
!  REAL(r_std), DIMENSION(npts,nvm), INTENT(in)  ::  compt_ugb
  REAL(r_std),DIMENSION(npts,nvm),INTENT(in)    :: nb_ani
  REAL(r_std), DIMENSION(npts,nvm), INTENT(in)  ::  grazed_frac
  REAL(r_std), DIMENSION(npts,nvm), INTENT(in)  ::  import_yield
  REAL(r_std),DIMENSION(npts),INTENT(in)        :: t2m_14
  REAL(r_std),DIMENSION(npts,nlitt,nvm),INTENT(in)    :: litter_not_avail
!ENDJCADD
    !-
    ! 0.2 local
    !-
    ! date, real
    REAL(r_std) :: date_real
    ! PFT exists (equivalent to fpc_max > 0 for natural PFTs), real
    REAL(r_std),DIMENSION(npts,nvm) :: PFTpresent_real
    ! is the plant senescent ?
    ! (only for deciduous trees - carbohydrate reserve), real
    REAL(r_std),DIMENSION(npts,nvm) :: senescence_real

    REAL(r_std),DIMENSION(npts,nvm) :: begin_leaves_real

    ! in order for this PFT to be introduced,
    ! does it have to be present in an adjacent grid box? - real
    REAL(r_std),DIMENSION(npts,nvm) :: need_adjacent_real
    ! To store variables names for I/O
    CHARACTER(LEN=80) :: var_name
    ! string suffix indicating an index
    CHARACTER(LEN=10) :: part_str
    ! string suffix indicating litter type
    CHARACTER(LEN=3),DIMENSION(nlitt) :: litter_str
    ! string suffix indicating level
    CHARACTER(LEN=2),DIMENSION(nlevs) :: level_str
    ! temporary storage
    REAL(r_std),DIMENSION(1) :: xtmp
    ! index
    INTEGER(i_std) :: j,k,l,m
    CHARACTER(LEN=1),DIMENSION(nelements) :: element_str  !! string suffix indicating element
    REAL(r_std), DIMENSION(1) :: temp_global_years
    CHARACTER(LEN=6),DIMENSION(nbpools) :: pools_str
    REAL(r_std), DIMENSION(npts) :: ok_equilibrium_real    

    ! land cover change variables 
    ! products remaining in the 10/100 year-turnover pool after the annual release for each compartment
    ! (10 or 100 + 1 : input from year of land cover change)
    REAL(r_std),DIMENSION(npts,0:10),INTENT(in)                           :: prod10
    REAL(r_std),DIMENSION(npts,0:100),INTENT(in)                          :: prod100
    ! annual release from the 10/100 year-turnover pool compartments
    REAL(r_std),DIMENSION(npts,10),INTENT(in)                           :: flux10
    REAL(r_std),DIMENSION(npts,100),INTENT(in)                          :: flux100
    REAL(r_std), DIMENSION(npts), INTENT(in)                            :: convflux
    REAL(r_std), DIMENSION(npts), INTENT(in)                            :: cflux_prod10
    REAL(r_std), DIMENSION(npts), INTENT(in)                            :: cflux_prod100
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(in)         :: bm_to_litter
    REAL(r_std),DIMENSION(npts),INTENT(in)                              :: carb_mass_total
    !---------------------------------------------------------------------
    IF (bavard >= 3) WRITE(numout,*) 'Entering writerestart'
    !-
    ! 1 string definitions
    !-
    DO l=1,nlitt
       IF     (l == imetabolic) THEN
          litter_str(l) = 'met'
       ELSEIF (l == istructural) THEN
          litter_str(l) = 'str'
       ELSE
          STOP 'Define litter_str'
       ENDIF
    ENDDO
    !-
    DO l=1,nlevs
       IF     (l == iabove) THEN
          level_str(l) = 'ab'
       ELSEIF (l == ibelow) THEN
          level_str(l) = 'be'
       ELSE
          STOP 'Define level_str'
       ENDIF
    ENDDO
    !-
    DO l=1,nelements
       IF     (l == icarbon) THEN
          element_str(l) = ''
!!$       ELSEIF (l == initrogen) THEN
!!$          element_str(l) = '_n'
       ELSE
          STOP 'Define element_str'
       ENDIF
    ENDDO
    !-
    pools_str(1:nbpools) =(/'str_ab','str_be','met_ab','met_be','actif ','slow  ','passif'/)
    !-
    IF (is_root_prc) THEN
       CALL ioconf_setatt_p ('UNITS','-')
       CALL ioconf_setatt_p ('LONG_NAME',' ')
    ENDIF
    !-
    ! 2 run control
    !-
    ! 2.1 day counter
    !-
    IF (is_root_prc) THEN
       var_name = 'day_counter'
       xtmp(1) = day_counter
       CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, xtmp)
    ENDIF
    !-
    ! 2.2 time step of STOMATE in days
    !-
    IF (is_root_prc) THEN
       var_name = 'dt_days'
       xtmp(1) = dt_days
       CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, xtmp)
    ENDIF
    !-
    ! 2.3 date
    !-
    IF (is_root_prc) THEN
       var_name = 'date'
       date_real = REAL(date,r_std)
       xtmp(1) = date_real
       CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, xtmp)
    ENDIF
    !-
    ! 3 daily meteorological variables
    !-
    var_name = 'moiavail_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                moiavail_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_init_date'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    2, 1, itime, &
         &              gdd_init_date, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'litterhum_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                litterhum_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_min_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_min_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsurf_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                tsurf_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsoil_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nbdl, 1, itime, &
         &                tsoil_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'soilhum_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nbdl, 1, itime, &
         &                soilhum_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'precip_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                precip_daily, 'scatter', nbp_glo, index_g)
    !-
!! Arsene 19-08-2014 Add snowtemp_min and snowdz_min
       var_name = 'snowtemp_min'                                           !! Arsene 19-08-2014 Add
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nsnow, 1, itime, & !! Arsene 19-08-2014 Add
         &                snowtemp_min, 'scatter', nbp_glo, index_g)       !! Arsene 19-08-2014 Add
    !-
       var_name = 'snowdz_min'                                             !! Arsene 19-08-2014 Add
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nsnow, 1, itime, & !! Arsene 19-08-2014 Add
         &                snowdz_min, 'scatter', nbp_glo, index_g)         !! Arsene 19-08-2014 Add
!! Arsene 19-08-2014 Add snowtemp_min and snowdz_min
    !-
    ! Wetland CH4 methane
    !pss:+
    var_name = 'uo_0'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uo_0, 'scatter', nbp_glo, index_g)
    
    var_name = 'uold2_0'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_0, 'scatter', nbp_glo, index_g)
    
    var_name = 'uo_wet1'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uo_wet1, 'scatter', nbp_glo, index_g)
    
    var_name = 'uold2_wet1'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_wet1, 'scatter', nbp_glo, index_g)
    
    var_name = 'uo_wet2'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uo_wet2, 'scatter', nbp_glo, index_g)
 
    var_name = 'uold2_wet2'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_wet2, 'scatter', nbp_glo, index_g)

    var_name = 'uo_wet3'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
       &              uo_wet3, 'scatter', nbp_glo, index_g)
    
    var_name = 'uold2_wet3'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_wet3, 'scatter', nbp_glo, index_g)
    
    var_name = 'uo_wet4'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uo_wet4, 'scatter', nbp_glo, index_g)
    
    var_name = 'uold2_wet4'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvert, 1, itime, &
         &              uold2_wet4, 'scatter', nbp_glo, index_g)
 
    var_name = 'tsurf_year'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &              tsurf_year, 'scatter', nbp_glo, index_g)
   
    !pss:-


    !-
    ! 4 productivities
    !-
    var_name = 'gpp_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_daily, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'npp_daily'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                npp_daily, 'scatter', nbp_glo, index_g)
    !-
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'turnover_daily_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   turnover_daily(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO
    !-
    ! 5 monthly meteorological variables
    !-
    var_name = 'moiavail_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                moiavail_month, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'moiavail_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                moiavail_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 't2m_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_longterm, 'scatter', nbp_glo, index_g)
    !-
    IF (control%ok_dgvm) THEN
       var_name = 'tlong_ref'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
            &                tlong_ref, 'scatter', nbp_glo, index_g)
    ENDIF
    !-
    var_name = 't2m_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_month, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tseason'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                Tseason, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tseason_length'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                Tseason_length, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tseason_tmp'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                Tseason_tmp, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tmin_spring'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                Tmin_spring, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'Tmin_spring_time'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                Tmin_spring_time, 'scatter', nbp_glo, index_g)
    !-
    DO k=1,2
       WRITE(part_str,'(I2)') k
       IF (k < 10) part_str(1:1) = '0'
       var_name = 'onset_date_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &                onset_date(:,:,k), 'scatter', nbp_glo, index_g)
    ENDDO
    !-

    var_name = 't2m_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                t2m_week, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'tsoil_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nbdl, 1, itime, &
         &                tsoil_month, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'soilhum_month'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nbdl, 1, itime, &
         &                soilhum_month, 'scatter', nbp_glo, index_g)
    !-
    ! 6 fire probability
    !-
    var_name = 'fireindex'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                fireindex, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'firelitter'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                firelitter, 'scatter', nbp_glo, index_g)
    !-
    ! 7 maximum and minimum moisture availabilities for tropic phenology
    !-
    var_name = 'maxmoistr_last'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxmoiavail_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxmoistr_this'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxmoiavail_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'minmoistr_last'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                minmoiavail_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'minmoistr_this'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                minmoiavail_thisyear, 'scatter', nbp_glo, index_g)
    !-
    ! 8 maximum "weekly" GPP
    !-
    var_name = 'maxgppweek_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxgppweek_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxgppweek_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxgppweek_thisyear, 'scatter', nbp_glo, index_g)
    !-
    ! 9 annual GDD0
    !-
    var_name = 'gdd0_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                gdd0_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd0_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                gdd0_lastyear, 'scatter', nbp_glo, index_g)
    !-
    ! 10 annual precipitation
    !-
    var_name = 'precip_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                precip_thisyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'precip_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &                precip_lastyear, 'scatter', nbp_glo, index_g)
    !-
    ! 11 derived "biometeorological" variables
    !-
    var_name = 'gdd_m5_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gdd_m5_dormance, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_from_growthinit'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &              gdd_from_growthinit, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'gdd_midwinter'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gdd_midwinter, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ncd_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                ncd_dormance, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ngd_minus5'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                ngd_minus5, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'time_hum_min'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                time_hum_min, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'hum_min_dormance'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                hum_min_dormance, 'scatter', nbp_glo, index_g)
    !-
    ! 12 Plant status
    !-
    var_name = 'PFTpresent'
    WHERE ( PFTpresent(:,:) )
       PFTpresent_real = un
    ELSEWHERE
       PFTpresent_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                PFTpresent_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'ind'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                ind, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'turnover_time'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                turnover_time, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'adapted'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                adapted, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'regenerate'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                regenerate, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'npp_longterm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                npp_longterm, 'scatter', nbp_glo, index_g)
    !-
   !! Arsene 25-06-2014 NPPcumul
    var_name = 'npp0_cumul'                                                !! Arsene 25-06-2014 NPPcumul
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &   !! Arsene 25-06-2014 NPPcumul
         &                npp0_cumul, 'scatter', nbp_glo, index_g)         !! Arsene 25-06-2014 NPPcumul
   !! Arsene 25-06-2014 NPPcumul 
    !-
    var_name = 'lm_lastyearmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                lm_lastyearmax, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'lm_thisyearmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                lm_thisyearmax, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxfpc_lastyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxfpc_lastyear, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'maxfpc_thisyear'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                maxfpc_thisyear, 'scatter', nbp_glo, index_g)
    !-
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'turnover_longterm_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   turnover_longterm(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO
    !-
    var_name = 'gpp_week'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                gpp_week, 'scatter', nbp_glo, index_g)
    !-
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'biomass_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                   biomass(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO
    !-
    DO k=1,nparts
       WRITE(part_str,'(I2)') k
       IF (k < 10) part_str(1:1) = '0'
       var_name = 'maint_resp_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &                   resp_maint_part(:,:,k), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    DO m=1,nleafages
       WRITE(part_str,'(I2)') m
       IF (m < 10) part_str(1:1) = '0'
       var_name = 'leaf_age_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &                  leaf_age(:,:,m), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    DO m=1,nleafages
       WRITE(part_str,'(I2)') m
       IF (m < 10) part_str(1:1) = '0'
       var_name = 'leaf_frac_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &                   leaf_frac(:,:,m), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    var_name = 'senescence'
    WHERE ( senescence(:,:) )
       senescence_real = un
    ELSEWHERE
       senescence_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                senescence_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'begin_leaves'
    WHERE ( begin_leaves(:,:) )
       begin_leaves_real = un
    ELSEWHERE
       begin_leaves_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                begin_leaves_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'when_growthinit'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                when_growthinit, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'age'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                age, 'scatter', nbp_glo, index_g)
!! Arsene 27-08-2015 Add dia_cut
    var_name = 'dia_cut'                                                 !! Arsene 27-08-2015 Add
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, & !! Arsene 27-08-2015 Add
         &                dia_cut, 'scatter', nbp_glo, index_g)          !! Arsene 27-08-2015 Add
!! Arsene 27-08-2015 Add dia_cut
    !-
    ! 13 CO2
    !-
    var_name = 'resp_hetero'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_hetero, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'resp_maint'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_maint, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'resp_growth'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                resp_growth, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'co2_fire'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, 1, itime, &
         &                co2_fire, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'co2_to_bm_dgvm'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                co2_to_bm_dgvm, 'scatter', nbp_glo, index_g)
    !-
    ! 14 vegetation distribution after last light competition
    !-
    var_name = 'veget_lastlight'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                veget_lastlight, 'scatter', nbp_glo, index_g)
    !-
    ! 15 establishment criteria
    !-
    var_name = 'everywhere'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                everywhere, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'need_adjacent'
    WHERE (need_adjacent(:,:))
       need_adjacent_real = un
    ELSEWHERE
       need_adjacent_real = zero
    ENDWHERE
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                need_adjacent_real, 'scatter', nbp_glo, index_g)
    !-
    var_name = 'RIP_time'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
         &                RIP_time, 'scatter', nbp_glo, index_g)
    !-
    ! 16 black carbon
    !-
    var_name = 'black_carbon'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
         &                black_carbon, 'scatter', nbp_glo, index_g)
    !-
    ! 17 litter
    !-
    DO l=1,nlitt
       var_name = 'litterpart_'//litter_str(l)
       CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, 1, itime, &
            &                   litterpart(:,:,l), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    DO k = 1,nelements
       DO l = 1,nlevs
          DO m = 1,nvm
             WRITE (part_str, '(I2)') m
             IF (m<10) part_str(1:1)='0'
             var_name = 'litter_'//part_str(1:LEN_TRIM(part_str))//'_'//level_str(l)//element_str(k)
             CALL restput_p (rest_id_stomate, var_name, nbp_glo, nlitt, 1, itime, &
                  &                     litter(:,:,m,l,k), 'scatter', nbp_glo, index_g)
          ENDDO
       ENDDO
    END DO
    !-
    DO l=1,nlitt
       var_name = 'dead_leaves_'//litter_str(l)
       CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, 1, itime, &
            &                   dead_leaves(:,:,l), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    DO m=1,nvm
       WRITE (part_str, '(I2)') m
       IF (m<10) part_str(1:1)='0'
       var_name = 'carbon_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ncarb, 1, itime, &
            &                   carbon(:,:,m), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    DO l=1,nlevs
       var_name = 'lignin_struc_'//level_str(l)
       CALL restput_p &
            &      (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
            &       lignin_struc(:,:,l), 'scatter', nbp_glo, index_g)
    ENDDO
    !-
    ! 18 land cover change
    !-
    var_name = 'prod10'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 11, 1, itime, &
         &                prod10, 'scatter', nbp_glo, index_g)
    var_name = 'prod100'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 101, 1, itime, &
         &                prod100, 'scatter', nbp_glo, index_g)
    var_name = 'flux10'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 10, 1, itime, &
         &                flux10, 'scatter', nbp_glo, index_g)
    var_name = 'flux100'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 100, 1, itime, &
         &                flux100, 'scatter', nbp_glo, index_g)

    var_name = 'convflux'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              convflux, 'scatter', nbp_glo, index_g)
    var_name = 'cflux_prod10'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              cflux_prod10, 'scatter', nbp_glo, index_g)
    var_name = 'cflux_prod100'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              cflux_prod100, 'scatter', nbp_glo, index_g)
    DO l = 1,nelements
       DO k = 1,nparts
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0'
          var_name = 'bm_to_litter_'//part_str(1:LEN_TRIM(part_str))//element_str(l)
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                bm_to_litter(:,:,k,l), 'scatter', nbp_glo, index_g)
       ENDDO
    END DO

    var_name = 'carb_mass_total'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
         &              carb_mass_total, 'scatter', nbp_glo, index_g)

    ! 19 Permafrost carbon related
    DO m=1,nvm
       WRITE(part_str,'(I2)') m
       IF (m < 10) part_str(1:1) = '0'
       var_name = 'deepC_a_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ndeep, 1, itime, &
            &               deepC_a(:,:,M), 'scatter', nbp_glo, index_g)
    END DO
    DO m=1,nvm
      WRITE(part_str,'(I2)') m
      IF (m < 10) part_str(1:1) = '0'
      var_name = 'deepC_s_'//part_str(1:LEN_TRIM(part_str))
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, ndeep, 1, itime, &
   &               deepC_s(:,:,m), 'scatter', nbp_glo, index_g)
    END DO
    DO m=1,nvm
       WRITE(part_str,'(I2)') m
       IF (m < 10) part_str(1:1) = '0'
       var_name = 'deepC_p_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ndeep, 1, itime, &
            &               deepC_p(:,:,m), 'scatter', nbp_glo, index_g)
    END DO
    DO m=1,nvm
       WRITE(part_str,'(I2)') m
       IF (m < 10) part_str(1:1) = '0'
       var_name = 'O2_soil_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ndeep, 1, itime, &
            &               O2_soil(:,:,m), 'scatter', nbp_glo, index_g)
    END DO
    DO m=1,nvm
       WRITE(part_str,'(I2)') m
       IF (m < 10) part_str(1:1) = '0'
       var_name = 'CH4_soil_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, ndeep, 1, itime, &
            &               CH4_soil(:,:,m), 'scatter', nbp_glo, index_g)
    END DO
    DO m=1,nvm
       WRITE(part_str,'(I2)') m
       IF (m < 10) part_str(1:1) = '0'
       var_name = 'O2_snow_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nsnow, 1, itime, &
            &               O2_snow(:,:,m), 'scatter', nbp_glo, index_g)
    END DO
    DO m=1,nvm
       WRITE(part_str,'(I2)') m
       IF (m < 10) part_str(1:1) = '0'
       var_name = 'CH4_snow_'//part_str(1:LEN_TRIM(part_str))
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, nsnow, 1, itime, &
            &               CH4_snow(:,:,m), 'scatter', nbp_glo, index_g)
    END DO

    var_name = 'thawed_humidity'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
   &               thawed_humidity, 'scatter', nbp_glo, index_g)

    var_name = 'depth_organic_soil'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
   &               depth_organic_soil, 'scatter', nbp_glo, index_g)

    var_name = 'altmax'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
   &               altmax, 'scatter', nbp_glo, index_g)
   !Isa dbg : fixed_cryoturbation_depth -> fixed_cryoturb_depth
    var_name = 'fixed_cryoturb_depth'
    CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
   &               fixed_cryoturbation_depth, 'scatter', nbp_glo, index_g)

!JCADD
!-
!  var_name = 'sla_calc'
!  CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
! &              sla_calc, 'scatter', nbp_glo, index_g)
!-
  var_name = 'wshtotsum'
  CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
 &              wshtotsum, 'scatter', nbp_glo, index_g)
!-
  var_name = 'sr_ugb'
  CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
 &              sr_ugb, 'scatter', nbp_glo, index_g)
!-
  var_name = 'sla_calc'
  CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
 &              sla_calc, 'scatter', nbp_glo, index_g)
!-
  var_name = 'nb_ani'
  CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
 &              nb_ani, 'scatter', nbp_glo, index_g)
!-
  var_name = 'grazed_frac'
  CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
 &              grazed_frac, 'scatter', nbp_glo, index_g)
!-
  var_name = 'import_yield'
  CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
 &              import_yield, 'scatter', nbp_glo, index_g)
!-
  var_name = 't2m_14'
  CALL restput_p (rest_id_stomate, var_name, nbp_glo,    1, 1, itime, &
 &              t2m_14, 'scatter', nbp_glo, index_g)

    DO l=1,nlitt
       var_name = 'litter_not_avail_'//litter_str(l)
       CALL restput_p (rest_id_stomate, var_name, nbp_glo,  nvm, 1, itime, &
            &                   litter_not_avail(:,l,:), 'scatter', nbp_glo, index_g)
    ENDDO
!ENDJCADD
 
    !-
    ! 19. Spinup
    !-
!    IF (spinup_analytic) THEN

       IF (is_root_prc) THEN
          temp_global_years(1) = REAL(global_years)
          var_name='Global_years'
          CALL restput (rest_id_stomate, var_name, 1, 1, 1, itime, temp_global_years)
       ENDIF
       
       var_name = 'nbp_sum'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &              nbp_accu, 'scatter', nbp_glo, index_g)

       var_name = 'nbp_flux'
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &              nbp_flux, 'scatter', nbp_glo, index_g)

       var_name = 'ok_equilibrium'
       WHERE(ok_equilibrium(:))
          ok_equilibrium_real = un
       ELSEWHERE
          ok_equilibrium_real = zero
       ENDWHERE
       CALL restput_p (rest_id_stomate, var_name, nbp_glo, 1, 1, itime, &
            &               ok_equilibrium_real, 'scatter', nbp_glo, index_g)
       
       DO k = 1,nbpools
          DO j = 1,nbpools
             WRITE(part_str,'(I2)') k
             IF (k < 10) part_str(1:1) = '0'             
             var_name = 'MatrixV_'//part_str(1:LEN_TRIM(part_str))//'_'//TRIM(pools_str(j))
             CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
                  &                MatrixV(:,:,k,j), 'scatter', nbp_glo, index_g)
          ENDDO
       ENDDO
          
       DO k = 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'Vector_U_'//part_str(1:LEN_TRIM(part_str))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                VectorU(:,:,k), 'scatter', nbp_glo, index_g)
       ENDDO
          
       DO k = 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'previous_stock_'//part_str(1:LEN_TRIM(part_str))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                previous_stock(:,:,k), 'scatter', nbp_glo, index_g)
       ENDDO
          
       DO k = 1,nbpools
          WRITE(part_str,'(I2)') k
          IF (k < 10) part_str(1:1) = '0' 
          var_name = 'current_stock_'//part_str(1:LEN_TRIM(part_str))
          CALL restput_p (rest_id_stomate, var_name, nbp_glo, nvm, 1, itime, &
               &                current_stock(:,:,k), 'scatter', nbp_glo, index_g)
       ENDDO
       
!    ENDIF !(spinup_analytic)
    !-

    IF (bavard >= 4) WRITE(numout,*) 'Leaving writerestart'
    !--------------------------
  END SUBROUTINE writerestart
  !-
  !===
  !-
  SUBROUTINE readbc (npts, lalo, resolution, tref)
    !---------------------------------------------------------------------
    !-
    ! 0.1 input
    !-
    ! Domain size
    INTEGER(i_std),INTENT(in) :: npts
    ! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std),DIMENSION (npts,2),INTENT(in) :: lalo
    ! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (npts,2),INTENT(in) :: resolution
    !-
    ! 0.2 not necessarily output
    !-
    ! "long term" reference 2 meter temperatures (K)
    REAL(r_std),DIMENSION(npts),INTENT(inout) :: tref
    !---------------------------------------------------------------------
    !-
    ! If the vegetation is static, then the long-term reference
    ! temperature is a boundary condition.
    !-
    IF ( .NOT. control%ok_dgvm ) THEN
       CALL get_reftemp (npts, lalo, resolution, tref)
    ENDIF
    !--------------------
  END SUBROUTINE readbc
  !-
  !===
  !-
  SUBROUTINE get_reftemp_clear
    !---------------------------------------------------------------------
    firstcall=.TRUE.
    IF (ALLOCATED (trefe)) DEALLOCATE( trefe )
    !-------------------------------
  END SUBROUTINE get_reftemp_clear
  !-
  !===
  !-
  SUBROUTINE get_reftemp (npts, lalo, resolution, tref_out)
    !---------------------------------------------------------------------
    !- read the long-term reference temperature from a boundary condition
    !- file. If the vegetation is dynamic, this field is used to
    !- initialize correctly the (prognostic) long-term reference
    !- temperature (in the case it is not found in the restart file).
    !- If the vegetation is static, the field read here is a real boundary
    !- condition that is not modified by the model.
    !---------------------------------------------------------------------
    !-
    ! 0 declarations
    !-
    ! 0.1 input
    !-
    ! Domain size
    INTEGER(i_std),INTENT(in) :: npts
    ! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std),DIMENSION (npts,2),INTENT(in) :: lalo
    ! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (npts,2),INTENT(in) :: resolution
    !-
    ! 0.2 output
    !-
    ! reference temperature (K)
    REAL(r_std), DIMENSION(npts),INTENT(out) :: tref_out
    !-
    ! 0.3 local
    !-
    INTEGER(i_std),PARAMETER :: nbvmax=200
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: &
         &  iml, jml, lml, tml, fid, ib, ip, jp, fopt, ilf, lastjp
    REAL(r_std) :: lev(1), date, dt, coslat
    INTEGER(i_std)                                 :: itau(1)
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: &
         &  lat_rel, lon_rel, lat_ful, lon_ful, tref_file
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: &
         &  loup_rel, lolow_rel, laup_rel, lalow_rel
    REAL(r_std) :: lon_up, lon_low, lat_up, lat_low
    REAL(r_std) :: ax, ay, sgn
    REAL(r_std),DIMENSION(nbvmax) :: area
    REAL(r_std),DIMENSION(nbvmax) :: tt
    REAL(r_std) :: resx, resy
    LOGICAL :: do_again
    !---------------------------------------------------------------------
    !-
    ! 1 If this is the first call, calculate the reference temperature
    !   and keep it in memory
    !-
    IF (firstcall) THEN
       !---
       !-- 1.1 only do this once
       !---
       firstcall = .FALSE.
       !---
       !-- 1.2 allocate the field
       !---
       ALLOCATE( trefe(npts) )
       !---
       !-- 1.3 read and interpolate the temperature file
       !---
       !-- Needs to be a configurable variable
       !---
       !Config Key   = REFTEMP_FILE
       !Config Desc  = Name of file from which the reference temperature is read
       !Config If    = OK_STOMATE
       !Config Def   = reftemp.nc
       !Config Help  = The name of the file to be opened to read
       !Config         the reference surface temperature.
       !Config         The data from this file is then interpolated
       !Config         to the grid of of the model.
       !Config         The aim is to get a reference temperature either
       !Config         to initialize the corresponding prognostic model
       !Config         variable correctly (ok_dgvm=TRUE) or to impose it
       !Config         as boundary condition (ok_dgvm=FALSE)
       !Config Units = [FILE]
       !---
       filename = 'reftemp.nc'
       CALL getin_p('REFTEMP_FILE',filename)
       !---
       IF (is_root_prc) CALL flininfo(filename,iml, jml, lml, tml, fid)
       CALL bcast(iml)
       CALL bcast(jml)
       CALL bcast(lml)
       CALL bcast(tml)
       !---
       ALLOCATE (lat_rel(iml,jml))
       ALLOCATE (lon_rel(iml,jml))
       ALLOCATE (laup_rel(iml,jml))
       ALLOCATE (loup_rel(iml,jml))
       ALLOCATE (lalow_rel(iml,jml))
       ALLOCATE (lolow_rel(iml,jml))
       ALLOCATE (lat_ful(iml+2,jml+2))
       ALLOCATE (lon_ful(iml+2,jml+2))
       ALLOCATE (tref_file(iml,jml))
       !---
       IF (is_root_prc) CALL flinopen (filename, .FALSE., iml, jml, lml, &
            &                                   lon_rel, lat_rel, lev, tml, itau, date, dt, fid)
       CALL bcast(lon_rel)
       CALL bcast(lat_rel)
       CALL bcast(itau)
       CALL bcast(date)
       CALL bcast(dt)

       !---
       IF (is_root_prc) CALL flinget (fid, 'temperature', iml, jml, lml, tml, &
            &                                  1, 1, tref_file)
       CALL bcast(tref_file)
       !---
       IF (is_root_prc) CALL flinclo (fid)
       !---
       !-- Duplicate the border assuming we have a global grid
       !-- going from west to east
       !---
       lon_ful(2:iml+1,2:jml+1) = lon_rel(1:iml,1:jml)
       lat_ful(2:iml+1,2:jml+1) = lat_rel(1:iml,1:jml)
       !---
       IF ( lon_rel(iml,1) < lon_ful(2,2)) THEN
          lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)
          lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
       ELSE
          lon_ful(1,2:jml+1) = lon_rel(iml,1:jml)-360
          lat_ful(1,2:jml+1) = lat_rel(iml,1:jml)
       ENDIF
       !---
       IF ( lon_rel(1,1) > lon_ful(iml+1,2)) THEN
          lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)
          lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
       ELSE
          lon_ful(iml+2,2:jml+1) = lon_rel(1,1:jml)+360
          lat_ful(iml+2,2:jml+1) = lat_rel(1,1:jml)
       ENDIF
       !---
       sgn = lat_rel(1,1)/ABS(lat_rel(1,1))
       lat_ful(2:iml+1,1) = sgn*180 - lat_rel(1:iml,1)
       sgn = lat_rel(1,jml)/ABS(lat_rel(1,jml))
       lat_ful(2:iml+1,jml+2) = sgn*180 - lat_rel(1:iml,jml)
       lat_ful(1,1) = lat_ful(iml+1,1)
       lat_ful(iml+2,1) = lat_ful(2,1)
       lat_ful(1,jml+2) = lat_ful(iml+1,jml+2)
       lat_ful(iml+2,jml+2) = lat_ful(2,jml+2)
       !---
       !-- Add the longitude lines to the top and bottom
       !---
       lon_ful(:,1) = lon_ful(:,2)
       lon_ful(:,jml+2) = lon_ful(:,jml+1)
       !---
       !-- Get the upper and lower limits of each grid box
       !---
       DO ip=1,iml
          DO jp=1,jml
             loup_rel(ip,jp) = &
                  &        MAX(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), &
                  &            0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
             lolow_rel(ip,jp) = &
                  &        MIN(0.5*(lon_ful(ip,jp+1)+lon_ful(ip+1,jp+1)), &
                  &            0.5*(lon_ful(ip+1,jp+1)+lon_ful(ip+2,jp+1)))
             laup_rel(ip,jp) = &
                  &        MAX(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), &
                  &            0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
             lalow_rel(ip,jp) = &
                  &        MIN(0.5*(lat_ful(ip+1,jp)+lat_ful(ip+1,jp+1)), &
                  &            0.5*(lat_ful(ip+1,jp+1)+lat_ful(ip+1,jp+2)))
          ENDDO
       ENDDO
       !---
       !-- Now we take each grid point and find out which values
       !-- from the forcing we need to average
       !---
       DO ib=1,npts
          !-----
          resx = resolution(ib,1)
          resy = resolution(ib,2)
          !-----
          do_again = .TRUE.
          !-----
          DO WHILE (do_again)
             !-----
             do_again = .FALSE.
             !-------
             !------ We find the 4 limits of the grid-box.
             !------ As we transform the resolution of the model into longitudes
             !------ and latitudes we do not have the problem of periodicity.
             !------ coslat is a help variable here !
             !-------
             coslat = MAX(COS(lalo(ib,1)*pi/180.),mincos)*pi/180.*R_Earth
             !-------
             lon_up  = lalo(ib,2)+resx/(2.0*coslat)
             lon_low = lalo(ib,2)-resx/(2.0*coslat)
             !-------
             coslat  = pi/180.*R_Earth
             !-------
             lat_up  = lalo(ib,1)+resy/(2.0*coslat)
             lat_low = lalo(ib,1)-resy/(2.0*coslat)
             !-------
             !------ Find the grid boxes from the data that go into
             !------ the model's boxes.
             !------ We still work as if we had a regular grid !
             !------ Well it needs to be localy regular so that
             !------ the longitude at the latitude of the last found point
             !------ is close to the one of the next point.
             !-------
             fopt = 0
             lastjp = 1
             DO ip=1,iml
                !---------
                !-------- Either the center of the data grid point is in the interval
                !-------- of the model grid or the East and West limits of the data
                !-------- grid point are on either sides of the border of the data grid
                !---------
                IF (      lon_rel(ip,lastjp) > lon_low &
                     &            .AND. lon_rel(ip,lastjp) < lon_up &
                     &             .OR. lolow_rel(ip,lastjp) < lon_low &
                     &            .AND. loup_rel(ip,lastjp) > lon_low &
                     &             .OR. lolow_rel(ip,lastjp) < lon_up &
                     &            .AND. loup_rel(ip,lastjp) > lon_up ) THEN
                   DO jp=1,jml
                      !-------------
                      !------------ Now that we have the longitude let us find the latitude
                      !-------------
                      IF (      lat_rel(ip,jp) > lat_low &
                           &                 .AND. lat_rel(ip,jp) < lat_up &
                           &                  .OR. lalow_rel(ip,jp) < lat_low &
                           &                 .AND. laup_rel(ip,jp) > lat_low &
                           &                  .OR. lalow_rel(ip,jp) < lat_up &
                           &                 .AND. laup_rel(ip,jp) > lat_up) THEN
                         lastjp = jp
                         !---------------
                         fopt = fopt + 1
                         IF ( fopt > nbvmax) THEN
                            WRITE(numout,*) &
                                 &                       'Please increase nbvmax in subroutine get_reftemp',ib
                            STOP
                         ELSE
                            !-----------------
                            !---------------- Get the area of the fine grid in the model grid
                            !-----------------
                            coslat = MAX(COS(lat_rel(ip,jp)*pi/180.),mincos)
                            ax =  ( MIN(lon_up,loup_rel(ip,jp)) &
                                 &                       -MAX(lon_low,lolow_rel(ip,jp))) &
                                 &                     *pi/180.*R_Earth*coslat
                            ay =  ( MIN(lat_up,laup_rel(ip,jp)) &
                                 &                       -MAX(lat_low,lalow_rel(ip,jp))) &
                                 &                     *pi/180.*R_Earth
                            area(fopt) = ax*ay
                            tt(fopt) = tref_file(ip,jp)
                         ENDIF
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
             !-------
             !------ Check that we found some points
             !-------
             trefe(ib) = zero
             !-------
             IF (fopt == 0) THEN
                do_again = .TRUE.
                !-------
                !------ increase search radius
                !-------
                resx = resx*2.
                resy = resy*2.
                IF ( resx > 2.*pi*R_Earth .OR. resy > pi*R_Earth ) THEN
                   STOP 'get_reftemp: found no point'
                ENDIF
             ELSE
                sgn = zero
                !-------
                !------ Compute the average surface air temperature
                !-------
                DO ilf=1,fopt
                   trefe(ib) = trefe(ib) + tt(ilf) * area(ilf)
                   sgn = sgn + area(ilf)
                ENDDO
                !-------
                !------ Normalize the surface
                !-------
                IF (sgn < min_sechiba) THEN
                   do_again = .TRUE.
                   !---------
                   !-------- increase search radius
                   !---------
                   resx = resx * 2.
                   resy = resy * 2.
                   IF ( resx > 2.*pi*R_Earth .OR. resy > pi*R_Earth ) THEN
                      STOP 'get_reftemp: found no point'
                   ENDIF
                ELSE
                   trefe(ib) = trefe(ib) / sgn
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       !-
       ! transform into Kelvin
       !-
       trefe(:) = trefe(:) + ZeroCelsius
       !-
       ! deallocate
       !-
       DEALLOCATE (lat_rel)
       DEALLOCATE (lon_rel)
       DEALLOCATE (laup_rel)
       DEALLOCATE (loup_rel)
       DEALLOCATE (lalow_rel)
       DEALLOCATE (lolow_rel)
       DEALLOCATE (lat_ful)
       DEALLOCATE (lon_ful)
       DEALLOCATE (tref_file)
    ENDIF
    !-
    ! 2 output the reference temperature
    !-
    tref_out(:) = trefe(:)
    !-------------------------
  END SUBROUTINE get_reftemp
  !-
  !===
  !-
END MODULE stomate_io
