! =================================================================================================================================
! MODULE       : slowproc
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Groups the subroutines that: (1) initialize all variables used in 
!! slowproc_main, (2) prepare the restart file for the next simulation, (3) Update the 
!! vegetation cover if needed, and (4) handle all slow processes if the carbon
!! cycle is activated (call STOMATE) or update the vegetation properties (LAI and 
!! fractional cover) in the case of a run with only SECHIBA.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_sechiba/slowproc.f90 $
!! $Date: 2014-05-23 12:34:33 +0200 (Fri, 23 May 2014) $
!! $Revision: 2135 $
!! \n
!_ ================================================================================================================================

MODULE slowproc

  ! modules used:

  USE defprec
  USE constantes 
  USE constantes_soil
  USE pft_parameters
  USE ioipsl
  USE xios_orchidee
  USE ioipsl_para
  USE sechiba_io
  USE interpol_help
  USE stomate
  USE stomate_data
  USE grid
  USE mod_orchidee_para

  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC slowproc_main,slowproc_clear

!!??VB Are these old comments? Questions you have yourself. I'm not clear what it's about anyway.
  ! Specific To use OLD or NEW interpolation schemes depending on the given arguments
  ! For interpol the OLD routines correspond to the tag 1.6 and should be dropped ???
  ! For interpol ONLY the subroutines *_g are used ?? (should delete the others ?)
  INTERFACE slowproc_interlai
     MODULE PROCEDURE slowproc_interlai_OLD, slowproc_interlai_NEW
  END INTERFACE
  INTERFACE slowproc_interpol
     MODULE PROCEDURE slowproc_interpol_OLD, slowproc_interpol_NEW
  END INTERFACE
  INTERFACE slowproc_interpol_g
     MODULE PROCEDURE slowproc_interpol_OLD_g, slowproc_interpol_NEW_g
  END INTERFACE

  !
  ! variables used inside slowproc module : declaration and initialisation
  !
!!?? Many units are lacking (even flags should be stated as "unitless")
!!?? You should use a sign (eg. !£!) which can be automatically searched and points to obsolete/dispensable variables
  LOGICAL, SAVE                                   :: l_first_slowproc = .TRUE.  !! Flag to do the initialisation only one time
!$OMP THREADPRIVATE(l_first_slowproc)
  REAL(r_std), SAVE                               :: dt_slow                    !! time step of slow processes and STOMATE
!$OMP THREADPRIVATE(dt_slow)
  !
  REAL(r_std), SAVE                                :: slope_default = 0.1
!$OMP THREADPRIVATE(slope_default)
  !
  INTEGER(i_std) , SAVE                           :: veget_update=0             !! update frequency in years for landuse (nb of years)
!$OMP THREADPRIVATE(veget_update)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: clayfraction            !! Clayfraction (0-1, unitless)
!$OMP THREADPRIVATE(clayfraction)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: laimap                  !! LAI map when the LAI is prescribed and not calculated by STOMATE
!$OMP THREADPRIVATE(laimap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: soilclass_default
!$OMP THREADPRIVATE(soilclass_default)
  !
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: veget_nextyear          !! next year fraction of vegetation type (0-1, unitless)
!$OMP THREADPRIVATE(veget_nextyear)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_nobio_nextyear     !! next year fraction of ice+lakes+cities+... (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_nextyear)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: totfrac_nobio_nextyear  !! next year total fraction of ice+lakes+cities+... (0-1, unitless)
!$OMP THREADPRIVATE(totfrac_nobio_nextyear)

  ! Update vegetation and bare soil fraction
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: tot_bare_soil    !! total evaporating bare soil fraction 
!$OMP THREADPRIVATE(tot_bare_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_bare        !! Fraction (of veget_max) of bare soil in each vegetation type           
!$OMP THREADPRIVATE(frac_bare)
  PUBLIC tot_bare_soil,frac_bare

CONTAINS

!!?? The "brief" section should be more precise (check and drop the "possibly").
!!?? Rather than "main routine that managed variable utilisation", be more precise:
!!?? Calls slowproc_init (which intialize variables) ...
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_main
!!
!>\BRIEF         Main routine that manage variable initialisation (slowproc_init), 
!! prepare the restart file with the slowproc variables, update the time variables 
!! for slow processes, and possibly update the vegetation cover, before calling 
!! STOMATE in the case of the carbon cycle activated or just update LAI (and possibly
!! the vegetation cover) for simulation with only SECHIBA   
!!
!!
!! DESCRIPTION  : (definitions, functional, design, flags): The subroutine manages 
!! diverses tasks:
!! (1) Initializing all variables of slowproc (first call)
!! (2) Preparation of the restart file for the next simulation with all prognostic variables
!! (3) Compute and update time variable for slow processes
!! (4) Update the vegetation cover if there is some land use change (only every years)
!! (5) Call STOMATE for the runs with the carbone cycle activated (ok_stomate) and compute the respiration
!!     and the net primary production
!! (6) Compute the LAI and possibly update the vegetation cover for run without STOMATE 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):  ::co2_flux, ::fco2_lu, ::lai, ::height, ::veget, ::frac_nobio,  
!! ::veget_max, ::totfrac_nobio, ::soiltype, ::assim_param, ::deadleaf_cover, ::qsintmax,
!! and resp_maint, resp_hetero, resp_growth, npp that are calculated and stored
!! in stomate is activated.  
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : 
! \latexonly 
! \includegraphics(scale=0.5){SlowprocMainFlow.eps} !PP to be finalize!!)
! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_main (kjit, kjpij, kjpindex, dtradia, date0, &
       ldrestart_read, ldrestart_write, ldforcing_write, ldcarbon_write, &
       IndexLand, indexveg, lalo, neighbours, resolution, contfrac, soiltile, reinf_slope, &
       t2m, t2m_min, temp_sol, stempdiag, &
       humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
       deadleaf_cover, &
       assim_param, &
       lai, frac_age, height, veget, frac_nobio, njsc, veget_max, totfrac_nobio, qsintmax, &
       rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
       co2_flux, fco2_lu, temp_growth,&
       tdeep, hsdeep_long, snow, heat_Zimov, pb, &
       sfluxCH4_deep, sfluxCO2_deep, &
       thawed_humidity, depth_organic_soil, zz_deep, zz_coef_deep, &
       soilc_total,snowdz,snowrho, snowtemp)    !! Arsene 19-08-2014 Add snowtemp (and snowdz) to stomate
  
!! INTERFACE DESCRIPTION

!! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjit                !! Time step number
    INTEGER(i_std), INTENT(in)                          :: kjpij               !! Total size of the un-compressed grid
    INTEGER(i_std),INTENT(in)                           :: kjpindex            !! Domain size - terrestrial pixels only
    REAL(r_std),INTENT(in)                              :: dtradia             !! Time step of SECHIBA
    REAL(r_std),INTENT (in)                             :: date0               !! Initial date of what ???
    LOGICAL, INTENT(in)                                 :: ldrestart_read      !! Logical for _restart_ file to read
    LOGICAL, INTENT(in)                                 :: ldrestart_write     !! Logical for _restart_ file to write
    LOGICAL, INTENT(in)                                 :: ldforcing_write     !! Logical for _forcing_ file to write
    LOGICAL, INTENT(in)                                 :: ldcarbon_write      !! Logical for _carbon_forcing_ file to write
    INTEGER(i_std),INTENT (in)                          :: rest_id,hist_id     !! _Restart_ file and _history_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist2_id            !! _history_ file 2 identifier
    INTEGER(i_std),INTENT (in)                          :: rest_id_stom        !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist_id_stom        !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                           :: hist_id_stom_IPCC   !! STOMATE's IPCC _history_ file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: IndexLand           !! Indices of the points on the land map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg            !! Indices of the points on the vegetation (3D map ???) 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo                !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,8), INTENT(in)  :: neighbours          !! neighoring grid points if land. In what ??? indices or geographical coordinates (lat/lon) ??? 
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution          !! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: contfrac            !! Fraction of continent in the grid (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)  :: humrel              !! Relative humidity ("moisture stress") (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: t2m                 !! 2 m air temperature (K)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: t2m_min             !! min. 2 m air temp. during forcing time step (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_sol            !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)  :: stempdiag           !! Soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)  :: shumdiag            !! Relative soil moisture (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: litterhumdiag       !! Litter humidity  (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: precip_rain         !! Rain precipitation (mm dt_slow^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: precip_snow         !! Snow precipitation (mm dt_slow^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: gpp                 !! GPP of total ground area (gC m^{-2} time step^{-1}). 
                                                                               !! Calculated in sechiba, account for vegetation cover and 
                                                                               !! effective time step to obtain gpp_d   
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (in)    :: tdeep      !! deep temperature profile
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (in)    :: hsdeep_long!! deep long term soil humidity profile
    REAL(r_std), DIMENSION(kjpindex),         INTENT (in)        :: snow       !! Snow mass [Kg/m^2]
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)               :: pb         !! Lowest level pressure
    REAL(r_std), DIMENSION(ndeep),   INTENT (in)                 :: zz_deep    !! deep vertical profile
    REAL(r_std), DIMENSION(ndeep),   INTENT (in)                 :: zz_coef_deep!! deep vertical profile   
    REAL(r_std), DIMENSION(kjpindex,nsnow),INTENT(in)            :: snowdz     !! snow depth for each layer
    REAL(r_std), DIMENSION(kjpindex,nsnow),INTENT(in)            :: snowrho    !! snow density for each layer    !! Arsene 19-08-2014 Need
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)          :: snowtemp   !! snow temperature      !! Arsene 19-08-2014 Add

 
!! 0.2 Output variables 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)  :: co2_flux            !! CO2 flux per average ground area (gC m^{-2} dt_slow^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: fco2_lu             !! CO2 flux from land-use (without forest management) (gC m^{-2} dt_slow^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: temp_growth         !! Growth temperature (°C) - Is equal to t2m_month 
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (out)   :: heat_Zimov !! heating associated with decomposition
    REAL(r_std), DIMENSION(kjpindex),     INTENT (out)           :: sfluxCH4_deep      !! surface flux of CH4 to atmosphere from permafrost
    REAL(r_std), DIMENSION(kjpindex),     INTENT (out)           :: sfluxCO2_deep      !! surface flux of CO2 to atmosphere from permafrost
 
!! 0.3 Modified variables
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(inout)       :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: lai            !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: height         !! height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(inout):: frac_age   !! Age efficacity from STOMATE for isoprene
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: veget          !! Fraction of vegetation type including none biological fraction (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout)  :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: veget_max      !! Maximum fraction of vegetation type including none biological fraction (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(inout)    :: soiltile       !! Fraction of each soil tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)          :: reinf_slope    !! slope coef for reinfiltration
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (inout):: assim_param    !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: deadleaf_cover !! Fraction of soil covered by dead leaves (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: qsintmax       !! Maximum water storage on vegetation from interception (mm)
    REAL(r_std), DIMENSION(kjpindex),   INTENT (inout)           :: thawed_humidity    !! specified humidity of thawed soil
    REAL(r_std), DIMENSION(kjpindex),   INTENT (inout)           :: depth_organic_soil !! how deep is the organic soil?
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (inout) :: soilc_total        !! total soil carbon for use in thermal

!! 0.4 Local variables
    INTEGER(i_std)                                     :: j, jv                !! indices 
    INTEGER(i_std), SAVE                               :: lcanop               !! canopy levels used for LAI
!$OMP THREADPRIVATE(lcanop)
    REAL(r_std)                                        :: tmp_day(1)           !! temporary variable for I/O
    CHARACTER(LEN=80)                                  :: var_name             !! To store variables names for I/O
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_maint           !! Maitanance component of autotrophic respiration in (gC m^{-2} dt_slow^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_hetero          !! heterotrophic resp. (gC/(m**2 of total ground)/time step)
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_growth          !! Growth component of autotrophic respiration in gC m^{-2} dt_slow^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: npp                  !! Net Ecosystem Exchange (gC/(m**2 of total ground)/time step)
    INTEGER(i_std) , SAVE                              :: veget_year           !! year for landuse
!$OMP THREADPRIVATE(veget_year)
    LOGICAL                                            :: land_use_updated     !! logical for landuse
    LOGICAL, PARAMETER                                 :: check = .FALSE.      !! logical 
    REAL(r_std), SAVE                                  :: sec_old = zero       !! Time counter for to save the old  
!$OMP THREADPRIVATE(sec_old)

!_ ================================================================================================================================

! 1 Initialize all variables at the first call 
!
    IF (l_first_slowproc) THEN

       ! 1.1 Perform the allocation of all variables, define some files and some flags. 
       !     Restart file read for Sechiba.
       IF (long_print) WRITE (numout,*) ' l_first_slowproc : call slowproc_init '
       CALL slowproc_init (kjit, ldrestart_read, dtradia, date0, kjpindex, IndexLand, lalo, neighbours, resolution, contfrac, &
            & rest_id, read_lai, lai, frac_age, veget, frac_nobio, totfrac_nobio, soiltile, reinf_slope, veget_max, njsc, &
            & height, lcanop, veget_update, veget_year)

       !
       ! 1.2 Define Time step in days for stomate
       dt_days = dt_slow / one_day

       ! 1.3 Set to zero all respiration fluxes
       resp_maint(:,:)  = zero
       resp_hetero(:,:) = zero
       resp_growth(:,:) = zero

       ! 1.4 check time step coherence between slow processes and fast processes
       IF ( dt_slow .LT. dtradia ) THEN
          WRITE(numout,*) 'slow_processes: time step smaller than forcing time step.'
          STOP 'slowproc_main'
       ENDIF

       ! 1.5 Call stomate to initialize all variables manadged in stomate,
       !     when Stomate is used or when the "watchout"?? diagnostics are requested 
       IF ( control%stomate_watchout .OR. control%ok_stomate ) THEN
          CALL stomate_main (kjit, kjpij, kjpindex, dtradia, dt_slow, &
               ldrestart_read, ldrestart_write, ldforcing_write, ldcarbon_write, &
               IndexLand, lalo, neighbours, resolution, contfrac, totfrac_nobio, clayfraction, &
               t2m, t2m_min, temp_sol, stempdiag, &
               humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
               deadleaf_cover, &
               assim_param, &
               lai, frac_age, height, veget, veget_max, &
               veget_nextyear, totfrac_nobio_nextyear, &
               hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
               co2_flux, fco2_lu, resp_maint,resp_hetero,resp_growth,temp_growth,&
               tdeep, hsdeep_long, snow, heat_Zimov, pb, &
               sfluxCH4_deep, sfluxCO2_deep, &
               thawed_humidity, depth_organic_soil, zz_deep, &
               zz_coef_deep, soilc_total,snowdz,snowrho, snowtemp)    !! Arsene 19-08-2014 Add snowtemp (and snowdz) to stomate
       ENDIF

       ! 1.6 Specific run without the carbon cycle (STOMATE not called): 
       !     Need to initialize some variables that will be used in SECHIBA:
       !     height, deadleaf_cover, assim_param, lai, qsintmax.
       IF ( .NOT. control%ok_stomate ) THEN

          CALL slowproc_derivvar (kjpindex, veget, lai, &
               qsintmax, deadleaf_cover, assim_param, height, temp_growth, humrel)    !! Arsene 28-05-2014 Add humrel
       ELSE
          qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
          qsintmax(:,1) = zero
       ENDIF

       RETURN

    ENDIF

    ! 2 preparation of a restart file for the next simulation 
    !  that will contain all pronostic variables to be used
    !  use restput routine for scalar and restput_p routine for arrays
    IF (ldrestart_write) THEN

       IF (long_print) WRITE (numout,*) ' we have to complete restart file with SLOWPROC variables '

       ! 2.1 Write a series of variables controled by slowproc: day
       ! counter, vegetation fraction, max vegetation fraction, LAI
       ! variable from stomate, fraction of bare soil, soiltype
       ! fraction, clay fraction, height of vegetation, map of LAI

       ! write day counter
       tmp_day(1) = day_counter
       var_name= 'day_counter'
       IF (is_root_prc) CALL restput (rest_id, var_name, 1 , 1  , 1, kjit, tmp_day)

!!$       var_name= 'veget'
!!$       CALL restput_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, veget, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'veget_max'
       CALL restput_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, veget_max, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'lai'
       CALL restput_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, lai, 'scatter',  nbp_glo, index_g)
       !
       var_name= 'frac_nobio'
       CALL restput_p (rest_id, var_name, nbp_glo, nnobio, 1, kjit, frac_nobio, 'scatter',  nbp_glo, index_g)

       IF ( control%ok_sechiba ) THEN

          var_name= 'soiltile_frac'
          CALL restput_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, soiltile, 'scatter',  nbp_glo, index_g)
          !
          var_name= 'njsc'
          CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, REAL(njsc, r_std), 'scatter',  nbp_glo, index_g)

          IF ( control%hydrol_cwrr ) THEN
             var_name= 'reinf_slope'
             CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, reinf_slope, 'scatter',  nbp_glo, index_g)
          END IF

       END IF !( control%ok_sechiba )

       var_name= 'clay_frac'
       CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, clayfraction, 'scatter',  nbp_glo, index_g)
       !
       ! The height of the vegetation could in principle be recalculated at the beginning of the run.
       ! However, this is very tedious, as many special cases have to be taken into account. This variable
       ! is therefore saved in the restart file.
       var_name= 'height'
       CALL restput_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, height, 'scatter',  nbp_glo, index_g)
       !
       ! Specific case where the LAI is read and not calculated by STOMATE: need to be saved
       IF (read_lai) THEN     
          var_name= 'laimap'
          CALL restput_p (rest_id, var_name, nbp_glo, nvm, 12, kjit, laimap)
       ENDIF
       !
       ! If there is some land use change, write the year for the land use ??? 
       IF (land_use) THEN
          tmp_day(1) = REAL(veget_year,r_std)
          var_name='veget_year'
          IF (is_root_prc) CALL restput (rest_id, var_name, 1 , 1  , 1, kjit, tmp_day)
       ENDIF
       
       ! 2.2 call STOMATE to write specific variables managed by STOMATE in the RESTART files 
       ! when Stomate is used or when the "watchout"?? diagnostics are requested 
       IF ( control%stomate_watchout .OR. control%ok_stomate ) THEN
          CALL stomate_main (kjit, kjpij, kjpindex, dtradia, dt_slow, &
               ldrestart_read, ldrestart_write, ldforcing_write, ldcarbon_write, &
               IndexLand, lalo, neighbours, resolution, contfrac, totfrac_nobio, clayfraction, &
               t2m, t2m_min, temp_sol, stempdiag, &
               humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
               deadleaf_cover, &
               assim_param, &
               lai, frac_age, height, veget, veget_max, &
               veget_nextyear, totfrac_nobio_nextyear, &
               hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
               co2_flux, fco2_lu, resp_maint,resp_hetero,resp_growth,temp_growth,&
               tdeep, hsdeep_long, snow, heat_Zimov, pb, &
               sfluxCH4_deep, sfluxCO2_deep, &
               thawed_humidity, depth_organic_soil, zz_deep, &
               zz_coef_deep, soilc_total,snowdz,snowrho, snowtemp)    !! Arsene 19-08-2014 Add snowtemp (and snowdz) to stomate
       ENDIF

       RETURN

    ENDIF


    ! 3 Compute and update all variables linked to the date and time
    !   counter for the call of STOMATE
    !
    ! 3.1 update day counter
    day_counter = day_counter + dtradia
    IF (check) WRITE(numout,*) "slowproc: day_counter 3",day_counter

    ! 3.2 Check if we are turning from one day to the next day 
    !     Assert all slow processes (days and years)
    IF ( sec_old >= one_day - dtradia .AND.  sec >= zero ) THEN

       ! 3.2.1 reset daily counter
       day_counter = zero
       IF (check) WRITE(numout,*) "slowproc: day_counter 2",day_counter

       ! 3.2.2 Activate slow processes
       do_slow = .TRUE.

       ! 3.2.3 Count the number of days 
       date = date + nint(dt_days)
       IF (check) WRITE(numout,*) "New date : ",date, 'year_length ',year_length,kjit

       ! 3.2.4 Calculate if we have just turned into a new year
       !       Note that: EndOfYear must be true once per year
       !       during a call of stomate_season ??
       IF ( month == 1 .AND. day == 1 .AND. sec .LT. dtradia ) THEN
          EndOfYear = .TRUE.
       ELSE
          EndOfYear = .FALSE.
       ENDIF
    
    ! 3.3 We are not a new day so that slow processes in STOMATE will
    ! not be done (but why end of year set to false ???)
    ELSE
       do_slow = .FALSE.
       EndOfYear = .FALSE.
    ENDIF

    ! 3.4 Reset old time counter to the new value
    sec_old = sec

    IF ( EndOfYear ) THEN
       WRITE(numout,*)  'slowproc: EndOfYear is activated.'
    ENDIF

    ! 4 Define the Land USE and update the vegetation for the next year
    !   This is done only every "veget_update" years (veget_update correspond
    !   to a number of years between each vegetation updates)
    land_use_updated=.FALSE.
    IF ( (land_use) .AND. (veget_update .GT. 0) ) THEN

       ! Check if we are effectively at the end of the year 
       IF ( EndOfYear ) THEN
          !
          veget_year = veget_year + 1

          ! Update of the vegetation cover with Land Use only if 
          ! the current year match the requested condition (a multiple of
          ! "veget_update")
          IF ( MOD(veget_year - veget_year_orig, veget_update) == 0 ) THEN
          
             WRITE(numout,*)  'We are updating land use veget for year =' , veget_year
             
             ! Call the routine to update the vegetation (output is veget_nextyear)
             CALL slowproc_update(kjpindex, lalo, neighbours, resolution, contfrac, &
               &               veget_max, frac_nobio, veget_nextyear, frac_nobio_nextyear, veget_year)
           
             !! this is added by cdk when he incorporating this permafrost
             !carbon into ORCHIDEE. This is because some pft disappeared in the
             !new map
             WHERE(veget_nextyear(:,:).LT.veget(:,:))
                veget(:,:) = veget_nextyear(:,:)
             ENDWHERE
 
             DO j = 1, nnobio
                totfrac_nobio_nextyear(:) = totfrac_nobio_nextyear(:) + frac_nobio_nextyear(:,j)
             ENDDO
             land_use_updated=.TRUE.
             !
          ENDIF
          !
       ENDIF
    ENDIF 
    IF (EndOfYear .AND. .NOT. land_use_updated) THEN
       lcchange=.FALSE.
    ENDIF

    ! 5 call STOMATE, either because we want to keep track of
    !   long-term variables (WATCHOUT case) or just because STOMATE is
    !   activated
    IF ( control%stomate_watchout .OR. control%ok_stomate ) THEN

       ! 5.1 call stomate main routine that will call all c-cycle routines       !
       CALL stomate_main (kjit, kjpij, kjpindex, dtradia, dt_slow, &
            ldrestart_read, ldrestart_write, ldforcing_write, ldcarbon_write, &
            IndexLand, lalo, neighbours, resolution, contfrac, totfrac_nobio, clayfraction, &
            t2m, t2m_min, temp_sol, stempdiag, &
            humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
            deadleaf_cover, &
            assim_param, &
            lai, frac_age, height, veget, veget_max, &
            veget_nextyear, totfrac_nobio_nextyear, &
            hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
            co2_flux, fco2_lu, resp_maint,resp_hetero,resp_growth,temp_growth,&
            tdeep, hsdeep_long, snow, heat_Zimov, pb, &
            sfluxCH4_deep, sfluxCO2_deep, &
            thawed_humidity, depth_organic_soil, zz_deep, &
            zz_coef_deep, soilc_total,snowdz,snowrho, snowtemp)    !! Arsene 19-08-2014 Add snowtemp (and snowdz) to stomate

       ! 5.2 output the respiration terms and the net primary
       !     production (NPP) that are calculated in STOMATE


       IF ( control%ok_stomate .AND. control%ok_sechiba ) THEN
          ! 5.2.1 Output the 3 respiration terms
          CALL xios_orchidee_send_field("maint_resp",resp_maint/dt_sechiba)
          CALL xios_orchidee_send_field("hetero_resp",resp_hetero/dt_sechiba)
          CALL xios_orchidee_send_field("growth_resp",resp_growth/dt_sechiba)

          CALL histwrite_p(hist_id, 'maint_resp', kjit, resp_maint, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist_id, 'hetero_resp', kjit, resp_hetero, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist_id, 'growth_resp', kjit, resp_growth, kjpindex*nvm, indexveg)

          ! 5.2.2 Compute the net primary production as the diff from
          ! Gross primary productin and the growth and maintenance
          ! respirations (WHY making separate case for index 1 (bare
          ! soil) ?)
          npp(:,1)=zero
          DO j = 2,nvm
             npp(:,j) = gpp(:,j) - resp_growth(:,j) - resp_maint(:,j)
          ENDDO

          CALL xios_orchidee_send_field("npp",npp/dt_sechiba)

          CALL histwrite_p(hist_id, 'npp', kjit, npp, kjpindex*nvm, indexveg)

       ENDIF

       IF ( hist2_id > 0 ) THEN
          IF ( control%ok_stomate ) THEN
             CALL histwrite_p(hist2_id, 'maint_resp', kjit, resp_maint, kjpindex*nvm, indexveg)
             CALL histwrite_p(hist2_id, 'hetero_resp', kjit, resp_hetero, kjpindex*nvm, indexveg)
             CALL histwrite_p(hist2_id, 'growth_resp', kjit, resp_growth, kjpindex*nvm, indexveg)
             CALL histwrite_p(hist2_id, 'npp', kjit, npp, kjpindex*nvm, indexveg)
          ENDIF
       ENDIF

    ENDIF

    ! 6 STOMATE is not activated: we have to compute the LAI  
    !   and possibly update the vegetation fraction if there was
    !   updated land use 
    IF ( .NOT. control%ok_stomate ) THEN

       IF (land_use_updated) THEN
          veget_max(:,:)=veget_nextyear(:,:)
          frac_nobio(:,:)=frac_nobio_nextyear(:,:)
       ENDIF

       IF ( do_slow ) THEN
          !
          ! 2.2.1 do daily processes if necessary
          !
          IF (long_print) WRITE (numout,*) 'slowproc_main : We update the daily variables'

          ! 6.1.1 Updates the LAI each dt_slow (day)
          CALL slowproc_lai (kjpindex, lcanop,stempdiag, &
               lalo,resolution,lai,month,day,read_lai,laimap)

          frac_age(:,:,1) = un
          frac_age(:,:,2) = zero
          frac_age(:,:,3) = zero
          frac_age(:,:,4) = zero
       ENDIF

       ! 6.3 Define the CO2 flux from the grid point to zero (no carbone cycle)
       co2_flux(:,:) = zero

    ELSE IF (land_use_updated) THEN
       !
       frac_nobio(:,:)=frac_nobio_nextyear(:,:)
       !
    ENDIF

    !
    ! 2.4 do daily processes if necessary
    !
    IF ( do_slow .OR. land_use_updated ) THEN
       !
       ! 2.4.1 updates veget
       CALL slowproc_veget (kjpindex, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)
       
       ! 2.4.2 updates qsintmax and other derived variables
       IF ( .NOT. control%ok_stomate ) THEN
          CALL slowproc_derivvar (kjpindex, veget, lai, &
               qsintmax, deadleaf_cover, assim_param, height, temp_growth, humrel)    !! Arsene 28-05-2014 Add humrel
       ELSE
          qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
          qsintmax(:,1) = zero
       ENDIF
       !
    END IF

    !
    ! Do some basic tests on the surface fractions updated above
    !
    CALL slowproc_checkveget(kjpindex, frac_nobio, veget_max, veget, soiltile)

    IF ( control%ok_sechiba ) THEN
       CALL xios_orchidee_send_field("tot_bare_soil",tot_bare_soil)
       CALL xios_orchidee_send_field("frac_bare",frac_bare)

       IF ( .NOT. almaoutput) THEN
          CALL histwrite_p(hist_id, 'tot_bare_soil', kjit, tot_bare_soil, kjpindex, IndexLand)
          CALL histwrite_p(hist_id, 'frac_bare', kjit, frac_bare, kjpindex*nvm, indexveg)
       ENDIF
    END IF

    IF (long_print) WRITE (numout,*) ' slowproc_main done '

  END SUBROUTINE slowproc_main


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_init
!!
!>\BRIEF         Initialisation of all variables linked to SLOWPROC
!!
!! DESCRIPTION  : (definitions, functional, design, flags): The subroutine manages 
!! diverses tasks:
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::lcanop, ::veget_update, ::veget_year,
!! ::lai, ::veget, ::frac_nobio, ::totfrac_nobio, ::veget_max, ::height, ::soiltype
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

SUBROUTINE slowproc_init (kjit, ldrestart_read, dtradia, date0, kjpindex, IndexLand, lalo, neighbours, resolution, contfrac, &
       & rest_id, read_lai, lai, frac_age, veget, frac_nobio, totfrac_nobio, soiltile, reinf_slope, veget_max, njsc, &
       & height, lcanop, veget_update, veget_year)

  !! INTERFACE DESCRIPTION

  !! 0. Parameters values
  LOGICAL, PARAMETER                                    :: check = .FALSE.!! Flag to check what???

  !! 0.1 input variables and scalar
  INTEGER(i_std), INTENT (in)                           :: kjit           !! Time step number
  LOGICAL, INTENT (in)                                  :: ldrestart_read !! Logical for _restart_ file to read
  REAL(r_std),INTENT (in)                               :: dtradia        !! Time step for SECHIBA (fast processes)
  REAL(r_std), INTENT (in)                              :: date0          !! intial date of the simulation ???
  INTEGER(i_std), INTENT (in)                           :: kjpindex       !! Domain size - Terrestrial pixels only 
  INTEGER(i_std), INTENT (in)                           :: rest_id        !! Restart file identifier

  INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: IndexLand      !! Indices of the land points on the map
  REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)       :: lalo           !! Geogr. coordinates (latitude,longitude) (degrees)
  INTEGER(i_std), DIMENSION (kjpindex,8), INTENT(in)    :: neighbours     !! Vector of neighbours for each grid point 1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
  REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)       :: resolution     !! size in x and y of the grid (m)
  REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: contfrac       !! Fraction of continent in the grid (unitless)
  LOGICAL, INTENT(in)                                   :: read_lai       !! Flag related to the read of laimap 

  !! 0.2 output variables and scalar
  INTEGER(i_std), INTENT(out)                           :: lcanop         !! Number of Canopy level used to compute LAI
  INTEGER(i_std), INTENT(out)                           :: veget_update   !! update frequency in timesteps (years or days ??) for landuse
  INTEGER(i_std), INTENT(out)                           :: veget_year     !! first year for landuse   (year or index ???)

  REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: lai            !! Leaf Area index (m^2 / m^2)
  REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: veget          !! Fraction of vegetation type (unitless)
  REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: frac_nobio     !! Fraction of ice,lakes,cities, ... (unitless)
  REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: totfrac_nobio  !! Total fraction of ice+lakes+cities+... (unitless)
  REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: veget_max      !! Max fraction of vegetation type (unitless)
  REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: height         !! Height of vegetation or surface in genral ??? (m)
  REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT (out):: frac_age !! Age efficacity from STOMATE for isoprene
  REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)   :: soiltile       !! Fraction of each soil tile (0-1, unitless)
  REAL(r_std), DIMENSION (kjpindex), INTENT(out)        :: reinf_slope    !! slope coef for reinfiltration
  INTEGER(i_std), DIMENSION(kjpindex), INTENT(out)      :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    
  !! 0.4 local scalar and variables 
  REAL(r_std)                                           :: tmp_day(1)        !! temporary variable 
  REAL(r_std)                                           :: tmp_veget_year(1) !! temporary variable
  REAL(r_std)                                           :: zcanop            !! ???? soil depth taken for canopy
  INTEGER(i_std)                                        :: vtmp(1)           !! temporary variable
  REAL(r_std), DIMENSION(nbdl)                          :: zsoil             !! soil depths at diagnostic levels
  INTEGER(i_std)                                        :: j,l               !! Indices
  CHARACTER(LEN=80)                                     :: var_name          !! To store variables names for I/O
  INTEGER(i_std)                                        :: ji, jv, ier       !! Indices 
  LOGICAL                                               :: get_slope
  REAL(r_std)                                           :: frac_nobio1       !! temporary variable for frac_nobio(see above)
  REAL(r_std), DIMENSION(kjpindex)                      :: tmp_real
  REAL(r_std), DIMENSION(kjpindex,nbdl)                 :: stempdiag2_bid    !! matrix to store stempdiag_bid
  REAL(r_std), DIMENSION (kjpindex,nscm)                :: soilclass         !! Fractions of each soil textural class in the grid cell (0-1, unitless)
  CHARACTER(LEN=4)                                      :: vegsoil_dist      !! Flag to choose the soil/vegetation distribution
  !
  CHARACTER(LEN=30), SAVE                               :: veget_str         !! update frequency for landuse
!$OMP THREADPRIVATE(veget_str)
  REAL(r_std),DIMENSION (nbp_glo,nnobio)                :: frac_nobio_g      !! Fraction of ice, lakes, cities etc. in the mesh (global)
  REAL(r_std),DIMENSION (nbp_glo,nvm)                   :: veget_max_g       !! Fraction of vegetation type at global scale
  REAL(r_std), DIMENSION(kjpindex)                      :: frac_crop_tot     !! Total fraction occupied by crops (0-1, unitless)
  LOGICAL, SAVE                                         :: read_veg_map_fr_restfile = .FALSE. 
!_ ================================================================================================================================

    ! 1 Allocation and initialisation of the different arrays: clayfraction, veget_nextyear,
    !   frac_nobio_nextyear
    !
    IF (long_print) WRITE (numout,*) "In slowproc_init"

    ALLOCATE (clayfraction(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in clayfraction allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'slowproc_init'
    END IF
    !? why there is a specific undef for sechiba ?
    clayfraction(:)=undef_sechiba

    ! Initialisation of the fraction of the different vegetation: Start with 100% of bare soil
    !? (RISKY to have bare soil index fixed to 1)
    ALLOCATE (soilclass_default(nscm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexveg allocation. We stop. We need nstm words = ',nscm
       STOP 'hydrol_init'
    END IF
    soilclass_default(:)=undef_sechiba
    !
    veget_max(:,1) = un
    veget_max(:,2:nvm) = zero
    frac_nobio(:,:) = zero
    totfrac_nobio(:) = zero

    ! Allocation of next year vegetation fraction in case of land use change
    ier=-1
    ALLOCATE(veget_nextyear(kjpindex, nvm), STAT=ier)
    IF (ier/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of veget_nextyear : ",ier
      STOP 
    ENDIF
    veget_nextyear(:,1) = un
    veget_nextyear(:,2:nvm) = zero

    ! Allocation of the fraction of non biospheric areas 
    ier=-1
    ALLOCATE(frac_nobio_nextyear(kjpindex, nnobio), STAT=ier)
    IF (ier/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of frac_nobio_nextyear : ",ier
      STOP 
    ENDIF
    frac_nobio_nextyear(:,:) = zero
    
    ier=-1
    ALLOCATE(totfrac_nobio_nextyear(kjpindex), STAT=ier)
    IF (ier/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of totfrac_nobio_nextyear : ",ier
      STOP 
    ENDIF
    !
    ALLOCATE (tot_bare_soil(kjpindex),STAT=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tot_bare_soil allocation. We stop. We need kjpindex words = ',kjpindex
       STOP
    END IF
    !
    ALLOCATE (frac_bare(kjpindex,nvm),STAT=ier) 
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in frac_bare allocation. We stop. We need kjpindex*nvm words = ',kjpindex*nvm
       STOP
    END IF

    IF (long_print) WRITE (numout,*) 'slowproc_init : End of allocations'

    !? Danger: it should be rather divided by nnobio if nnobio >1 
!MM must be corrected when nnobio > 1 !!
    totfrac_nobio_nextyear(:) = nnobio*un

    ! 2 read the day counter from the restart file: 
    ! val_exp is a default value to test if the variable exist in the restart file
    ! tmp_day is the current value of day_counter from the restart file: only readed
    ! by the main processor (parallelisation)
    ! Day_counter is set to zero if we are at the first time step of the simulation    
    !?? BUT DAY_counter is not used anymore: we should remove this variable
    var_name= 'day_counter'
    CALL ioconf_setatt_p('UNITS', 'd')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of computed day')
    IF (is_root_prc) THEN
       CALL restget (rest_id, var_name, 1 , 1  , 1, kjit, .TRUE., tmp_day)
       IF (tmp_day(1) == val_exp) THEN
          day_counter = zero
       ELSE
          day_counter = tmp_day(1)
       ENDIF
    ENDIF
    ! specific for parallelisation (synchronize day_counter)
    CALL bcast(day_counter)

    IF (long_print) WRITE (numout,*) 'slowproc_init : End of allocationsGot restart file set'

    ! get restart value if none were found in the restart file
    !

    !Config Key   = SECHIBA_DAY
    !Config Desc  = Time within the day simulated
    !Config if    = OK_SECHIBA
    !Config Def   = 0.0
    !Config Help  = This is the time spent simulating the current day. This variable is
    !Config         prognostic as it will trigger all the computations which are
    !Config         only done once a day.
    !Config Units = [days]
    !
    CALL setvar_p (day_counter, val_exp, 'SECHIBA_DAY', zero)

!!$    var_name= 'veget'
!!$    CALL ioconf_setatt_p('UNITS', '-')
!!$    CALL ioconf_setatt_p('LONG_NAME','Vegetation fraction')
!!$    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., veget, "gather", nbp_glo, index_g)
    !
    var_name= 'veget_max'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Maximum vegetation fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., veget_max, "gather", nbp_glo, index_g)
   
    ! Get frac_nobio from the restart file 
    frac_nobio(:,:) = val_exp
    var_name= 'frac_nobio'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Special soil type fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nnobio, 1, kjit, .TRUE., frac_nobio, "gather", nbp_glo, index_g)

    veget_update=0

    !
    ! control%do_land_use is not activate if land use map won't be updated
    control%do_land_use = .FALSE.
    !
    IF (land_use) THEN
       !
       var_name= 'veget_year'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Last year get in Land Use file.')
       IF (is_root_prc) THEN
          CALL restget (rest_id, var_name, 1       , 1  , 1, kjit, .TRUE., tmp_veget_year)
          !
          IF (tmp_veget_year(1) == val_exp) THEN
             veget_year=veget_year_orig
          ELSE
             IF (veget_reinit) THEN
                veget_year=veget_year_orig
             ELSE
                veget_year=INT(tmp_veget_year(1))
             ENDIF
          ENDIF
       ENDIF
       CALL bcast(veget_year)

       !
       !Config Key   = VEGET_UPDATE
       !Config Desc  = Update vegetation frequency
       !Config If    = LAND_USE
       !Config Def   = 0Y
       !Config Help  = The veget datas will be update each this time step.
       !Config Units = [years]
       !
       veget_update=0
       WRITE(veget_str,'(a)') '0Y'
     !?? danger : VEGET_UPDATE now called VEGET_LENGTH in Thomas L. run.def ??!!  
       CALL getin_p('VEGET_UPDATE', veget_str)
       !
       !
       l=INDEX(TRIM(veget_str),'Y')
       READ(veget_str(1:(l-1)),"(I2.2)") veget_update
       WRITE(numout,*) "Update frequency for land use in years :",veget_update
       !
       ! control%do_land_use is activated if update frequency is not zero.
       IF ( veget_update > 0 ) THEN
          control%do_land_use = .TRUE.
       ENDIF

       IF ( veget_update == 0 .AND. lcchange ) THEN
          CALL ipslerr_p (2,'slowproc_init', &
             &     'You have asked for LAND_COVER_CHANGE activated with VEGET_UPDATE = 0Y.',&
             &     'We can''t use this land cover change model if veget is not updated.', &
             &     'We have disabled it.')
          lcchange=.FALSE.
       ENDIF

    ENDIF
    !
    IF (long_print) WRITE (numout,*) 'slowproc_init : End of Land_Use configuration'

    IF ( control%ok_sechiba ) THEN
       
       var_name= 'soiltile_frac'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Fraction of each soil type')
       CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., soiltile, "gather", nbp_glo, index_g)

       ! Test if the soiltile_frac was in restartfile. If not, try the old name 'soiltype_frac'.
       IF ( ALL(soiltile(:,:) == val_exp ) ) THEN
          WRITE(numout,*) 'The variable soiltile_frac was not in restart file. Try to read old variable name soiltype_frac instead.'
          var_name= 'soiltype_frac'
          CALL restget_p (rest_id, var_name, nbp_glo, nstm, 1, kjit, .TRUE., soiltile, "gather", nbp_glo, index_g)
       END IF

       IF ( control%hydrol_cwrr ) THEN
          var_name= 'reinf_slope'
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Slope coef for reinfiltration')
          CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., reinf_slope, "gather", nbp_glo, index_g)
       END IF
       var_name= 'njsc'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Index of soil type')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., tmp_real, "gather", nbp_glo, index_g)
       WHERE ( tmp_real .LT. val_exp )
          njsc = NINT(tmp_real)
       ENDWHERE
       
    END IF !( control%ok_sechiba ) 
    !
    var_name= 'clay_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of clay in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., clayfraction, "gather", nbp_glo, index_g)
    !
    IF (long_print) WRITE (numout,*) 'slowproc_init : End CWRR configuration'
    !
    var_name= 'lai'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Leaf area index')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., lai, "gather", nbp_glo, index_g)
    !
    ! The height of the vegetation could in principle be recalculated at the beginning of the run.
    ! However, this is very tedious, as many special cases have to be taken into account. This variable
    ! is therefore saved in the restart file.
    var_name= 'height'
    CALL ioconf_setatt_p('UNITS', 'm')
    CALL ioconf_setatt_p('LONG_NAME','Height of vegetation')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., height, "gather", nbp_glo, index_g)
    
    !
    IF (read_lai)THEN
       !
       ALLOCATE (laimap(kjpindex,nvm,12),stat=ier)
       IF (ier.NE.0) THEN
          WRITE (numout,*) ' error in laimap allocation. We stop. We need kjpindex*nvm*12 words = ',kjpindex*nvm*12
          STOP 'slowproc_init'
       END IF
       laimap(:,:,:) = val_exp
       !
       var_name= 'laimap'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Leaf area index read')
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 12, kjit, .TRUE., laimap)
       !
    ELSE
       !
       ALLOCATE (laimap(1,1,1))
    ENDIF
    !
    !Config Key   = SECHIBA_ZCANOP
    !Config Desc  = Soil level used for canopy development (if STOMATE disactivated)
    !Config If    = OK_SECHIBA and .NOT. OK_STOMATE  
    !Config Def   = 0.5
    !Config Help  = The temperature at this soil depth is used to determine the LAI when
    !Config         STOMATE is not activated.
    !Config Units = [m]
    !
    zcanop = 0.5_r_std
    CALL setvar_p (zcanop, val_exp, 'SECHIBA_ZCANOP', 0.5_r_std)
    !
    ! depth at center of the levels
    zsoil(1) = diaglev(1) / 2.
    DO l = 2, nbdl
       zsoil(l) = ( diaglev(l) + diaglev(l-1) ) / 2.
    ENDDO

    ! index of this level
    vtmp = MINLOC ( ABS ( zcanop - zsoil(:) ) )
    lcanop = vtmp(1)

    !
    !  Interception reservoir coefficient
    !
    !Config Key   = SECHIBA_QSINT 
    !Config Desc  = Interception reservoir coefficient
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.1
    !Config Help  = Transforms leaf area index into size of interception reservoir
    !Config         for slowproc_derivvar or stomate
    !Config Units = [m]
    CALL getin_p('SECHIBA_QSINT', qsintcst)
    WRITE(numout, *)' SECHIBA_QSINT, qsintcst = ', qsintcst

    !
    ! Time step of STOMATE and LAI update
    !
    !Config Key   = DT_SLOW
    !Config Desc  = Time step of STOMATE and other slow processes
    !Config If    = OK_STOMATE
    !Config Def   = one_day
    !Config Help  = Time step (s) of regular update of vegetation
    !Config         cover, LAI etc. This is also the time step
    !Config         of STOMATE.
    !Config Units = [seconds]

    dt_slow = one_day
    CALL getin_p('DT_SLOW', dt_slow)

    !
    ! get restart value if none were found in the restart file
    !
    !
    IF ( impveg ) THEN
       !
       !  We are on a point and thus we can read the information from the run.def
       !
       !
       !Config Key   = SECHIBA_VEGMAX
       !Config Desc  = Maximum vegetation distribution within the mesh (0-dim mode)
       !Config If    = IMPOSE_VEG
       !Config Def   = 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0
       !Config Help  = The fraction of vegetation is read from the restart file. If
       !Config         it is not found there we will use the values provided here.
       !Config Units = [-]
       !
       CALL setvar_p (veget_max, val_exp, 'SECHIBA_VEGMAX', veget_ori_fixed_test_1)
       IF (SUM(veget_ori_fixed_test_1) > un) THEN
          CALL ipslerr_p (2,'slowproc_init', &
             &     'The sum of the fractions of the array SECHIBA_VEGMAX is greater than 1.', &
             &     'The sum should be equal to 1.', &
             &     'Check your configuration file.')
       ENDIF

       !
       !Config Key   = SECHIBA_FRAC_NOBIO
       !Config Desc  = Fraction of other surface types within the mesh (0-dim mode)
       !Config If    = IMPOSE_VEG
       !Config Def   = 0.0
       !Config Help  = The fraction of ice, lakes, etc. is read from the restart file. If
       !Config         it is not found there we will use the values provided here.
       !Config         For the moment, there is only ice.
       !Config Units = [-]
       !
     !?? problem with frac_nobio : everything intialized as if every frac_nobio point is 
        ! covered with ice sheets (continental ice) 
     !?? direct calculation possible from vegetmax not possible ?
       ! laisser ca tant qu'il n'y a que la glace. Pb avec setvar?
       frac_nobio1 = frac_nobio(1,1)
       CALL setvar_p (frac_nobio1, val_exp, 'SECHIBA_FRAC_NOBIO', frac_nobio_fixed_test_1)
       frac_nobio(:,:) = frac_nobio1
       ! CALL setvar (frac_nobio, val_exp, 'SECHIBA_FRAC_NOBIO', frac_nobio_fixed_test_1)

       !
       !Config Key   = SECHIBA_LAI
       !Config Desc  = LAI for all vegetation types (0-dim mode)
       !Config Def   = 0., 8., 8., 4., 4.5, 4.5, 4., 4.5, 4., 2., 2., 2., 2.
       !Config If    = IMPOSE_VEG
       !Config Help  = The maximum LAI used in the 0dim mode. The values should be found
       !Config         in the restart file. The new values of LAI will be computed anyway
       !Config         at the end of the current day. The need for this variable is caused
       !Config         by the fact that the model may stop during a day and thus we have not
       !Config         yet been through the routines which compute the new surface conditions.
       !Config Units = [-]
       !
       CALL setvar_p (lai, val_exp, 'SECHIBA_LAI', llaimax)

       IF (impsoilt) THEN

          ! If njsc is not in restart file, then initialize soilclass from values 
          ! from run.def file and recalculate njsc
          IF ( ALL(njsc(:) .EQ. undef_int )) THEN
             !Config Key   = SOIL_FRACTIONS
             !Config Desc  = Fraction of the 3 soil types (0-dim mode)
             !Config Def   = undef_sechiba
             !Config If    = IMPOSE_VEG and IMPOSE_SOILT
             !Config Help  = Determines the fraction for the 3 soil types
             !Config         in the mesh in the following order : sand loam and clay.
             !Config Units = [-]
          
             soilclass(1,:) = soilclass_default(:)
             CALL getin_p('SOIL_FRACTIONS',soilclass(1,:))
             DO ji=2,kjpindex
                soilclass(ji,:) = soilclass(1,:)
             ENDDO

             njsc(:) = 0
             DO ji = 1, kjpindex
                njsc(ji) = MAXLOC(soilclass(ji,:),1)
             ENDDO
          END IF

          !Config Key   = CLAY_FRACTION
          !Config Desc  = Fraction of the clay fraction (0-dim mode)
          !Config Def   = 0.2
          !Config If    = IMPOSE_VEG and IMPOSE_SOIL
          !Config Help  = Determines the fraction of clay in the grid box.
          !Config Units = [-] 
          !
          CALL setvar_p (clayfraction, val_exp, 'CLAY_FRACTION', clayfraction_default)
       ELSE
          IF ( MINVAL(clayfraction) .EQ. MAXVAL(clayfraction) .AND. MAXVAL(clayfraction) .EQ. val_exp .OR. &
               MINVAL(njsc) .EQ. MAXVAL(njsc) .AND. MAXVAL(njsc) .EQ. undef_int ) THEN
             
             CALL slowproc_soilt(kjpindex, lalo, neighbours, resolution, contfrac, soilclass, clayfraction)
             njsc(:) = 0
             DO ji = 1, kjpindex
                njsc(ji) = MAXLOC(soilclass(ji,:),1)
             ENDDO
          ENDIF
       ENDIF

       !Config Key   = REINF_SLOPE
       !Config Desc  = Slope coef for reinfiltration 
       !Config Def   = 0.1
       !Config If    = IMPOSE_VEG
       !Config Help  = Determines the reinfiltration ratio in the grid box due to flat areas
       !Config Units = [-]
       !
       slope_default=0.1
       CALL setvar_p (reinf_slope, val_exp, 'SLOPE', slope_default)

       !Config Key   = SLOWPROC_HEIGHT
       !Config Desc  = Height for all vegetation types 
       !Config Def   = 0., 30., 30., 20., 20., 20., 15., 15., 15., .5, .6, 1.0, 1.0
       !Config If    = OK_SECHIBA
       !Config Help  = The height used in the 0dim mode. The values should be found
       !Config         in the restart file. The new values of height will be computed anyway
       !Config         at the end of the current day. The need for this variable is caused
       !Config         by the fact that the model may stop during a day and thus we have not
       !Config         yet been through the routines which compute the new surface conditions.
       !Config Units = [m]
       !
       CALL setvar_p (height, val_exp, 'SLOWPROC_HEIGHT', height_presc)

    ELSE
       !
       !  We are in the full 2-D case thus we need to do the interpolation if the potential vegetation
       !  is not available
       !
       IF ( ALL( veget_max(:,:) .EQ. val_exp ) .OR. ALL( frac_nobio(:,:) .EQ. val_exp ) ) THEN

          IF ( .NOT. land_use ) THEN

             !Config Key  = SLOWPROC_VEGET_OLD_INTERPOL
             !Config Desc = Flag to use old "interpolation" of vegetation map.
             !Config If   = NOT IMPOSE_VEG and NOT LAND_USE
             !Config Def  = FALSE
             !Config Help = If you want to recover the old (ie orchidee_1_2
             !branch) 
             !Config        "interpolation" of vegetation map.
             !
             old_veget = .FALSE.
             CALL getin_p('SLOWPROC_VEGET_OLD_INTERPOL',old_veget)

             ! The interpolation of vegetation has changed.
             IF (is_root_prc) THEN
               CALL getin('read_veg_map_fr_restfile', read_veg_map_fr_restfile)
               IF ( .NOT. read_veg_map_fr_restfile ) THEN
                IF ( .NOT. old_veget ) THEN
                   ! NEW slowproc interpol :
                   CALL slowproc_interpol_g(nbp_glo, lalo_g, neighbours_g, resolution_g, contfrac_g, veget_max_g, frac_nobio_g)
                ELSE
                   ! OLD slowproc interpol :
                   CALL slowproc_interpol_g(nbp_glo, lalo_g, neighbours_g, resolution_g, veget_max_g, frac_nobio_g)
                ENDIF
               ELSE
                   CALL slowproc_read_veg_restfile(nbp_glo, lalo_g, neighbours_g, resolution_g, contfrac_g, &
                        veget_max_g, frac_nobio_g)
               ENDIF

             ENDIF

             CALL scatter(veget_max_g,veget_max)
             CALL scatter(frac_nobio_g, frac_nobio)
             !
             IF ( control%ok_dgvm ) THEN
                !
                ! If we are dealing with dynamic vegetation then all
                ! natural PFTs should be set to veget_max = 0
                !  And in case no agriculture is desired, agriculture PFTS should be
                ! set to 0 as well
             
                IF (agriculture) THEN
                   DO jv=2,nvm
                      IF (natural(jv)) THEN 
                         veget_max(:,jv)=zero
                      ENDIF
                   ENDDO
                   !
                   ! 1. Here we calculate the fraction of crop for each point
                   frac_crop_tot(:) = zero
                   ! 2. we sum only on the indexes corresponding to the non_natural pfts
                   DO jv = 2,nvm
                      IF (.NOT. natural(jv)) THEN
                         DO ji = 1, kjpindex               
                            frac_crop_tot(ji) = frac_crop_tot(ji) + veget_max(ji,jv)
                         ENDDO
                      ENDIF
                   ENDDO
                   ! 3. we calculate the fraction of bare soil
                   DO ji = 1, kjpindex
                      veget_max(ji,1) = un - frac_crop_tot(ji) - SUM(frac_nobio(ji,:))
                   ENDDO
                   !
                ELSE
                   veget_max(:,:) = zero
                   DO ji = 1, kjpindex
                      veget_max(ji,1) = un  - SUM(frac_nobio(ji,:))
                   ENDDO
                ENDIF
                !
             ENDIF
          ELSE
             WRITE(numout,*)  'We initialize land use veget for year =' , veget_year
             ! If restart doesn''t contain veget, then it is the first computation
             CALL slowproc_update(kjpindex, lalo, neighbours, resolution, contfrac, &
               &               veget_nextyear, frac_nobio_nextyear, veget_max, frac_nobio, veget_year, init=.TRUE.)
             !
             IF ( control%ok_dgvm  ) THEN
                !
                ! If we are dealing with dynamic vegetation then all
                ! natural PFTs should be set to veget_max = 0
                !  And in case no agriculture is desired, agriculture PFTS should be
                ! set to 0 as well
                
                IF (agriculture) THEN
                   !
                   DO jv = 2, nvm
                      IF (natural(jv)) THEN 
                         veget_max(:,jv)=zero
                      ENDIF
                   ENDDO
                   !
                   ! 1. Here we calculate the fraction of crop for each point
                   frac_crop_tot(:) = zero
                   ! 2. we sum only on the indexes corresponding to the non_natural pfts
                   DO jv = 2, nvm
                      IF(.NOT. natural(jv)) THEN
                         DO ji = 1, kjpindex
                            frac_crop_tot(ji) = frac_crop_tot(ji) + veget_max(ji,jv)
                         ENDDO
                      ENDIF
                   ENDDO
                   ! 3. we calculate the fraction of bare soil
                   DO ji = 1, kjpindex
                      veget_max(ji,1) = un - frac_crop_tot(ji) - SUM(frac_nobio(ji,:))   
                   ENDDO
                   !
                ELSE
                   veget_max(:,:) = zero
                   DO ji = 1, kjpindex
                      veget_max(ji,1) = un  - SUM(frac_nobio(ji,:))
                   ENDDO
                ENDIF
                !
             ENDIF
             !
          ENDIF
             !
       ELSE
          ! WITH restarts for vegetation and DGVM and NO AGRICULTURE
          IF ( control%ok_dgvm  ) THEN
             !
             ! If we are dealing with dynamic vegetation then all
             ! natural PFTs should be set to veget_max = 0
             !  And in case no agriculture is desired, agriculture PFTS should be
             ! set to 0 as well
             !
             IF (.NOT. agriculture) THEN
                ! 1 .We calculate the total fraction of crops for each point
                frac_crop_tot(:) = zero
                DO jv = 2, nvm  
                   IF ( .NOT. natural (jv))  THEN
                      DO ji = 1, kjpindex
                         frac_crop_tot(ji) = frac_crop_tot(ji) + veget_max(ji,jv)
                      ENDDO
                   ENDIF
                ENDDO
                ! 2. We add to the fraction of bare soil the fraction of crops
                DO ji = 1, kjpindex
                   veget_max(ji,1) = veget_max(ji,1) + frac_crop_tot(ji)
                ENDDO
                ! 3. We set the fractions of crops to zero
                DO jv = 2, nvm                  
                   IF ( .NOT. natural (jv))  THEN
                      veget_max(:,jv) = zero
                   ENDIF
                ENDDO
                !-
             ENDIF
             !
          ENDIF
          !
       ENDIF
       !
       !
       IF (read_lai) THEN

          !
          !  In case the LAI map was not found in the restart then we need to recompute it
          !
          IF ( ALL( laimap(:,:,:) .EQ. val_exp) ) THEN
             ! The interpolation of LAI has changed.
             IF ( .NOT. old_lai ) THEN
                ! NEW slowproc interlai :
                CALL slowproc_interlai (kjpindex, lalo, resolution,  neighbours, contfrac, laimap)
             ELSE
                ! OLD slowproc interlai :
                CALL slowproc_interlai(kjpindex, lalo,  resolution, laimap)
             ENDIF
             !
             ! Compute the current LAI just to be sure 
             !
             stempdiag2_bid(1:kjpindex,1:nbdl) = stempdiag_bid
             CALL slowproc_lai (kjpindex, lcanop, stempdiag2_bid, &
                  lalo,resolution,lai,month,day,read_lai,laimap)
             !
             !  Make sure that we redo the computation when we will be back in slowproc_main
             day_counter = dt_slow - dtradia
             !
          ENDIF
          !
       ENDIF
       !
       IF ( MINVAL(lai) .EQ. MAXVAL(lai) .AND. MAXVAL(lai) .EQ. val_exp) THEN
          !
          !  Get a first guess at the LAI
          !
          IF ( read_lai ) THEN
             IF ( ALL( laimap(:,:,:) .EQ. val_exp) ) THEN
                ! The interpolation of LAI has changed.
                IF ( .NOT. old_lai ) THEN
                   ! NEW slowproc interlai :
                   CALL slowproc_interlai (kjpindex, lalo, resolution,  neighbours, contfrac, laimap)
                ELSE
                   ! OLD slowproc interlai :
                   CALL slowproc_interlai(kjpindex, lalo,  resolution, laimap)
                ENDIF
             ENDIF
          ENDIF
          !
          stempdiag2_bid(1:kjpindex,1:nbdl) = stempdiag_bid
          CALL slowproc_lai (kjpindex, lcanop, stempdiag2_bid, &
               lalo,resolution,lai,month,day,read_lai,laimap)
          !
          ! If we start from scratch, we set lai to zero for consistency with stomate
          IF ( .NOT. read_lai ) THEN
             lai(:,:) = zero
          ENDIF
          !   
          frac_age(:,:,1) = un
          frac_age(:,:,2) = zero
          frac_age(:,:,3) = zero
          frac_age(:,:,4) = zero

          ! Make sure that we redo the computation when we will be back in slowproc_main
          day_counter = dt_slow - dtradia

       ENDIF

       IF ( MINVAL(height) .EQ. MAXVAL(height) .AND. MAXVAL(height) .EQ. val_exp) THEN
       
          ! Impose height
          DO jv = 1, nvm
             height(:,jv) = height_presc(jv)
          ENDDO

       ENDIF

       IF ( MINVAL(clayfraction) .EQ. MAXVAL(clayfraction) .AND. MAXVAL(clayfraction) .EQ. val_exp .OR. &
            MINVAL(njsc) .EQ. MAXVAL(njsc) .AND. MAXVAL(njsc) .EQ. undef_int ) THEN

          CALL slowproc_soilt(kjpindex, lalo, neighbours, resolution, contfrac, soilclass, clayfraction)

          njsc(:) = 0
          DO ji = 1, kjpindex
             njsc(ji) = MAXLOC(soilclass(ji,:),1)
          ENDDO
       ENDIF


       !Config Key   = GET_SLOPE
       !Config Desc  = Read slopes from file and do the interpolation
       !Config Def   = n
       !Config If    =
       !Config Help  = Needed for reading the slopesfile and doing the interpolation. This will be
       !               used by the re-infiltration parametrization
       !Config Units = [FLAG]
       get_slope = .FALSE.
       CALL getin_p('GET_SLOPE',get_slope)
       
       IF ( control%hydrol_cwrr ) THEN
          IF ( MINVAL(reinf_slope) .EQ. MAXVAL(reinf_slope) .AND. MAXVAL(reinf_slope) .EQ. val_exp .OR. get_slope) THEN
             
             CALL slowproc_slope(kjpindex, lalo, neighbours, resolution, contfrac, reinf_slope)
             
          ENDIF
       END IF
    ENDIF

    ! Have a first guess at the vegetation fraction
    CALL slowproc_veget (kjpindex, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)

    l_first_slowproc = .FALSE.

    IF (long_print) WRITE (numout,*) ' slowproc_init done '

  END SUBROUTINE slowproc_init

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_clear
!!
!>\BRIEF          Clear all variables related to slowproc and stomate modules  
!!
!_ ================================================================================================================================

  SUBROUTINE slowproc_clear 

  ! 1 Reset l_first_slowproc and clear all the variables defined as common for the routines in slowproc 

    l_first_slowproc = .TRUE.
    
    IF (ALLOCATED (clayfraction)) DEALLOCATE (clayfraction)
    IF (ALLOCATED (laimap)) DEALLOCATE (laimap)
    IF (ALLOCATED (veget_nextyear)) DEALLOCATE (veget_nextyear)
    IF (ALLOCATED (frac_nobio_nextyear)) DEALLOCATE (frac_nobio_nextyear)
    IF (ALLOCATED (totfrac_nobio_nextyear)) DEALLOCATE (totfrac_nobio_nextyear)
    !
    IF ( ALLOCATED (soilclass_default)) DEALLOCATE (soilclass_default)
    !
    IF (ALLOCATED  (tot_bare_soil)) DEALLOCATE (tot_bare_soil)
    IF (ALLOCATED  (frac_bare)) DEALLOCATE (frac_bare)
    !

 ! 2. Clear all the variables in stomate 

    CALL stomate_clear 
    !
  END SUBROUTINE slowproc_clear

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_derivvar
!!
!>\BRIEF         Initializes variables related to the
!! parameters to be assimilated, the maximum water on vegetation, the vegetation height, 
!! and the fraction of soil covered by dead leaves and the vegetation height 
!!
!! DESCRIPTION  : (definitions, functional, design, flags):
!! (1) Initialization of the variables relevant for the assimilation parameters  
!! (2) Intialization of the fraction of soil covered by dead leaves
!! (3) Initialization of the Vegetation height per PFT
!! (3) Initialization the maximum water on vegetation for interception with a particular treatement of the PFT no.1
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::qsintmax, ::deadleaf_cover, ::assim_param, ::height  
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_derivvar (kjpindex, veget, lai, &
       qsintmax, deadleaf_cover, assim_param, height, temp_growth, humrel)     !! Arsene 28-05-2014 Add humrel

    !! INTERFACE DESCRIPTION

    !! 0.1 Input scalar and fields 
    INTEGER(i_std), INTENT (in)                                :: kjpindex       !! Domain size - terrestrial pixels only
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)          :: veget          !! Fraction of pixel covered by PFT. Fraction accounts for none-biological land covers (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)          :: lai            !! PFT leaf area index (m^{2} m^{-2})

    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)          :: humrel         !! Arsene 28-05-2014 ADD Relative humidity ("moisture stress") (0-1, unitless)

    !! 0.2. Output scalar and fields 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)          :: qsintmax       !! Maximum water on vegetation for interception(mm)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)              :: deadleaf_cover !! fraction of soil covered by dead leaves (unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm,npco2), INTENT (out)   :: assim_param    !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)          :: height         !! height of the vegetation or surface in general ??? (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)              :: temp_growth    !! growth temperature (°C)  
    !
    !! 0.3 Local declaration
    INTEGER(i_std)                                              :: ji, jv         !! Local indices
!_ ================================================================================================================================

    !
    ! 1. Initialize (why here ??) the variables revelant for the assimilation parameters
    !
    DO jv = 1, nvm
       !assim_param(:,jv,ivcmax) = vcmax_fix(jv)
       assim_param(:,jv,ivcmax) = vcmax_fix(jv) * (0.5 + 0.5 *  humrel(:,jv))     !! Arsene 23-06-2014 For dessication. CARE. It's a little modification for the next time step
    ENDDO

    !
    ! 2. Intialize the fraction of soil covered by dead leaves 
    !
    deadleaf_cover(:) = zero

    !
    ! 3. Initialize the Vegetation height per PFT
    !
    DO jv = 1, nvm
       height(:,jv) = height_presc(jv)
    ENDDO
    !
    ! 4. Initialize the maximum water on vegetation for interception
    !
    qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)

    ! Added by Nathalie - July 2006
    !  Initialize the case of the PFT no.1 to zero 
    qsintmax(:,1) = zero

    temp_growth(:)=25.

  END SUBROUTINE slowproc_derivvar


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_mean
!!
!>\BRIEF          Accumulates field_in over a period of dt_tot.
!! Has to be called at every time step (dt). 
!! Mean value is calculated if ldmean=.TRUE.
!! field_mean must be initialized outside of this routine! 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) AcumAcuumlm 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::field_main
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_mean (npts, n_dim2, dt_tot, dt, ldmean, field_in, field_mean)

    !
    !! 0 declarations

    !! 0.1 input scalar and variables 
    INTEGER(i_std), INTENT(in)                           :: npts     !! Domain size- terrestrial pixels only 
    INTEGER(i_std), INTENT(in)                           :: n_dim2   !! Number of PFTs 
    REAL(r_std), INTENT(in)                              :: dt_tot   !! Time step of stomate (in days). The period over which the accumulation or the mean is computed 
    REAL(r_std), INTENT(in)                              :: dt       !! Time step in days 
    LOGICAL, INTENT(in)                                  :: ldmean   !! Flag to calculate the mean after the accumulation ???
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(in)      :: field_in !! Daily field 

    !! 0.3 Modified field; The computed sum or mean field over dt_tot time period depending on the flag ldmean 
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(inout)   :: field_mean !! Accumulated field at dt_tot time period or mean field over dt_tot 
 

!_ ================================================================================================================================

    !
    ! 1. Accumulation the field over dt_tot period 
    !
    field_mean(:,:) = field_mean(:,:) + field_in(:,:) * dt

    !
    ! 2. If the flag ldmean set, the mean field is computed over dt_tot period  
    !
    IF (ldmean) THEN
       field_mean(:,:) = field_mean(:,:) / dt_tot
    ENDIF

  END SUBROUTINE slowproc_mean


  
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_long
!!
!>\BRIEF        Calculates a temporally smoothed field (field_long) from
!! instantaneous input fields.Time constant tau determines the strength of the smoothing.
!! For tau -> infinity??, field_long becomes the true mean value of field_inst
!! (but  the spinup becomes infinietly long, too).
!! field_long must be initialized outside of this routine! 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) Testing the time coherence betwen the time step dt and the time tau over which
!! the rescaled of the mean is performed   
!!  (2) Computing the rescaled mean over tau period 
!! MAIN OUTPUT VARIABLE(S): field_long  
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::field_long
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_long (npts, n_dim2, dt, tau, field_inst, field_long)

    !
    ! 0 declarations
    !

    ! 0.1 input scalar and fields 

    INTEGER(i_std), INTENT(in)                                 :: npts        !! Domain size- terrestrial pixels only
    INTEGER(i_std), INTENT(in)                                 :: n_dim2      !! Second dimension of the fields, which represents the number of PFTs
    REAL(r_std), INTENT(in)                                    :: dt          !! Time step in days   
    REAL(r_std), INTENT(in)                                    :: tau         !! Integration time constant (has to have same unit as dt!)  
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(in)            :: field_inst  !! Instantaneous field 


    ! 0.2 modified field

    ! Long-term field
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(inout)         :: field_long  !! Mean value of the instantaneous field rescaled at tau time period 

!_ ================================================================================================================================

    !
    ! 1 test coherence of the time 

    IF ( ( tau .LT. dt ) .OR. ( dt .LE. zero ) .OR. ( tau .LE. zero ) ) THEN
       WRITE(numout,*) 'slowproc_long: Problem with time steps'
       WRITE(numout,*) 'dt=',dt
       WRITE(numout,*) 'tau=',tau
    ENDIF

    !
    ! 2 integration of the field over tau 

    field_long(:,:) = ( field_inst(:,:)*dt + field_long(:,:)*(tau-dt) ) / tau

  END SUBROUTINE slowproc_long


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_veget
!!
!>\BRIEF        Sum up veget_max and frac_nobio and test if sum is equal to 1.
!! In case there is a problem, a correction is performed when possible, otherwise
!! it stops.    
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) Summing up frac_nobio and veget_max. Test if the sum is equal to unity.
!!  If not, when possible a solution is proposed, otherwise the routine stops 
!! (2) Setting veget to veget_max 
!! (3) Testing the  lai of a vegetation type. If small the soil part of the pixel is increased
!! (4) Summing the surface fractions and test if it is equal to 1. If not, when
!! possible, a solution is given, otherwise the routine stops.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::veget
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_veget (kjpindex, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile)
    !
    ! 0. Declarations
    !

    ! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                             :: kjpindex    !! Domain size - terrestrial pixels only
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)       :: lai         !! PFT leaf area index (m^{2} m^{-2})

    ! 0.2 Modified variables 
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(inout) :: frac_nobio  !! Fraction of the mesh which is covered by ice, lakes, ...
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)    :: veget_max   !! Maximum fraction of vegetation type including none biological fraction (unitless)

    ! 0.3 Output variables 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)      :: veget       !! Fraction of pixel covered by PFT. Fraction accounts for none-biological land covers (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)    :: totfrac_nobio
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)   :: soiltile     !! Fraction of each soil tile (0-1, unitless)

    ! 0.4 Local scalar and varaiables 
    REAL(r_std), DIMENSION(kjpindex)                       :: fracsum     !! Sum of both fracnobio and veget_max 
    INTEGER(i_std)                                         :: nbad        !!
    INTEGER(i_std)                                         :: ji, jv, jst !! indices 

    ! Test Nathalie	
    REAL(r_std)                                            :: SUMveg      !!
!!!$    REAL(r_std)                                            :: trans_veg
!_ ================================================================================================================================

    ! 0. Set to zero fractions of frac_nobio and veget_max smaller than min_vegfrac
!    DO ji = 1, kjpindex
!       IF ( SUM(frac_nobio(ji,:)) .LT. min_vegfrac ) THEN
!          frac_nobio(ji,:) = zero
!       ENDIF
!       DO jv = 1, nvm
!          IF ( veget_max(ji,jv) .LT. min_vegfrac ) THEN
!             veget_max(ji,jv) = zero
!          ENDIF
!       ENDDO
!       SUMveg = SUM(frac_nobio(ji,:))+SUM(veget_max(ji,:))
!       frac_nobio(ji,:) = frac_nobio(ji,:)/SUMveg
!       veget_max(ji,:) = veget_max(ji,:)/SUMveg
!    ENDDO

    !
    ! 1. Sum up veget_max and frac_nobio and test if sum is equal to 1
    !
    !
    ! 1.1 Sum up
    !
    fracsum(:) = zero
    DO jv = 1, nnobio
       DO ji = 1, kjpindex
          fracsum(ji) = fracsum(ji) + frac_nobio(ji,jv)
       ENDDO
    ENDDO
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          fracsum(ji) = fracsum(ji) + veget_max(ji,jv)
       ENDDO
    ENDDO
    !
    ! 1.2 stop if there is a severe problem, that is an error above 0.01%
    !     The rest will be scaled
    !
    nbad = COUNT( ABS(fracsum(:)-un) .GT. 0.0001 )
    IF ( nbad .GT. 0 ) THEN
      WRITE(numout,*) 'Problem with total surface areas.'
      DO ji = 1, kjpindex
        IF ( ABS(fracsum(ji)-un) .GT. 0.0001 ) THEN
          WRITE(numout,*) 'Point :', ji
          WRITE(numout,*) '  frac_nobio :', frac_nobio(ji,:)
          WRITE(numout,*) '  Veget_max :', veget_max(ji,:)
          WRITE(numout,*) '  Fracsum :', fracsum(ji), EPSILON(un)
          WRITE(numout,*) '  The error is :', un - ( SUM(frac_nobio(ji,:)) + SUM(veget_max(ji,:)) )
          STOP 'slowproc_veget'
        ENDIF
      ENDDO
    ENDIF
    !
    ! 1.3 correction at places where the problem is surely precision-related
    !
    nbad = COUNT( ABS(fracsum(:)-un) .GT. EPSILON(un) )
    !
    IF ( nbad .GT. 0 ) THEN
       !
       IF ( long_print ) THEN
          WRITE(numout,*) 'WARNING! scaling frac_nobio and veget_max at', nbad, ' points'
       ENDIF
       !
       DO jv = 1, nnobio
          DO ji = 1, kjpindex
             IF ( ABS(fracsum(ji)-un) .GT. EPSILON(un) ) THEN
                frac_nobio(ji,jv) = frac_nobio(ji,jv)/fracsum(ji)
             ENDIF
          ENDDO
       ENDDO
       !
       DO jv = 1, nvm
          DO ji = 1, kjpindex
             IF ( ABS(fracsum(ji)-un) .GT. EPSILON(un) ) THEN
                veget_max(ji,jv) = veget_max(ji,jv)/fracsum(ji)
             ENDIF
          ENDDO
       ENDDO
       !
    ENDIF

    !
    ! 2. Set veget equal to veget_max
    !
    veget(:,1)=veget_max(:,1)
    !
    ! 3. if lai of a vegetation type (jv > 1) is small, increase soil part
    !
    ! Ajout Nouveau calcul (stomate-like) 
    DO ji = 1, kjpindex
       SUMveg = zero
       tot_bare_soil(ji) = veget_max(ji,1)
       DO jv = 2, nvm
          veget(ji,jv) = veget_max(ji,jv) * ( un - exp( - lai(ji,jv) * ext_coeff(jv) ) )
          tot_bare_soil(ji) = tot_bare_soil(ji) + (veget_max(ji,jv) - veget(ji,jv))
          SUMveg = SUMveg + veget(ji,jv)
       ENDDO
       SUMveg = SUMveg + tot_bare_soil(ji) + SUM(frac_nobio(ji,:)) 
       IF (SUMveg .LT. 0.99999) THEN
          WRITE(numout,*)' ATTENTION, en ji, SUMveg LT 1: ', ji, SUMveg
          WRITE(numout,*)'   frac_nobio = ',SUM(frac_nobio(ji,:))
          WRITE(numout,*)  '              ',veget(ji,:)
          WRITE(numout,*)  '              ',tot_bare_soil(ji)
       ENDIF
    ENDDO


    !
    ! 4. Sum up surface fractions and test if the sum is equal to 1
    !

    !
    ! 4.1 Sum up
    !
    fracsum(:) = zero
    DO jv = 1, nnobio
       DO ji = 1, kjpindex
          fracsum(ji) = fracsum(ji) + frac_nobio(ji,jv)
       ENDDO
    ENDDO
    DO jv = 1, nvm
       DO ji = 1, kjpindex
          fracsum(ji) = fracsum(ji) + veget_max(ji,jv)
       ENDDO
    ENDDO

    !
    ! 4.2 stop if there is a severe problem
    !
    nbad = COUNT( ABS(fracsum(:)-un) .GT. (REAL(nvm+nnobio,r_std)*EPSILON(un)) )
    !
    IF ( nbad .GT. 0 ) THEN
       WRITE(numout,*) 'Problem with veget or frac_nobio.'
       DO ji = 1, kjpindex
          IF ( ABS(fracsum(ji)-un) .GT. (10.*EPSILON(un)) ) THEN
             WRITE(numout,*) 'Point :', ji
             WRITE(numout,*) '  frac_nobio :', frac_nobio(ji,:)
             WRITE(numout,*) '  Veget :', veget(ji,:)
             WRITE(numout,*) '  The error is :', un - (SUM(frac_nobio(ji,:)) + SUM(veget(ji,:)))
             STOP 'slowproc_veget'
          ENDIF
       ENDDO
    ENDIF

    ! 
    ! 4.3 correction at places where the problem is surely precision-related
    !
    nbad = COUNT( ABS(fracsum(:)-un) .GT. EPSILON(un) )
    !
    IF ( nbad .GT. 0 ) THEN
       !
       IF ( long_print ) THEN
          WRITE(numout,*) 'slowproc_veget: scaling veget at', nbad, ' points'
       ENDIF
       !
       DO ji = 1, kjpindex
          IF ( ABS(fracsum(ji)-un) .GT. EPSILON(un) ) THEN
             veget_max(ji,1) = veget_max(ji,1) / fracsum(ji)
          ENDIF
       ENDDO
       veget(:,1)=veget_max(:,1)
       !
       DO jv = 2, nvm
          DO ji = 1, kjpindex
             IF ( ABS(fracsum(ji)-un) .GT. EPSILON(un) ) THEN
                veget_max(ji,jv) = veget_max(ji,jv) / fracsum(ji)
                veget(ji,jv) = veget(ji,jv) / fracsum(ji)
             ENDIF
          ENDDO
       ENDDO
       !
       DO jv = 1, nnobio
          DO ji = 1, kjpindex
             IF ( ABS(fracsum(ji)-un) .GT. EPSILON(un) ) THEN
                frac_nobio(ji,jv) = frac_nobio(ji,jv) / fracsum(ji)
             ENDIF
          ENDDO
       ENDDO

    ENDIF

    !
    tot_bare_soil(:) = veget_max(:,1)
    DO jv = 2, nvm
       DO ji =1, kjpindex
          tot_bare_soil(ji) = tot_bare_soil(ji) + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
    END DO
    
    DO ji =1, kjpindex
       IF( veget_max(ji,1) .GT. min_sechiba ) THEN
          frac_bare(ji,1) = un
       ELSE
          frac_bare(ji,1) = zero
       ENDIF
    ENDDO
    DO jv = 2, nvm
       DO ji =1, kjpindex
          IF( veget_max(ji,jv) .GT. min_sechiba ) THEN
             frac_bare(ji,jv) = un - veget(ji,jv)/veget_max(ji,jv)
          ELSE
             frac_bare(ji,jv) = zero
          ENDIF
       ENDDO
    ENDDO

    totfrac_nobio(:) = zero
    DO jv = 1, nnobio
       totfrac_nobio(:) = totfrac_nobio(:) + frac_nobio(:,jv)
    ENDDO
    
    ! Soiltiles are only used in hydrol, but we fix them in here because some time it might depend
    ! on a changing vegetation (but then some adaptation should be made to hydrol) and be also used
    ! in the other modules to perform separated energy balances
    soiltile(:,:) = zero
    soiltile(:,1) = totfrac_nobio(:)
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji = 1, kjpindex
          soiltile(ji,jst) = soiltile(ji,jst) + veget_max(ji,jv)
!!$          IF (jv==1) THEN
!!$             soiltile(ji,jst) = soiltile(ji,jst) + tot_bare_soil(ji)
!!$          ENDIF
       ENDDO
    ENDDO
!!$    ! Avoid soiltile < 0.01
!!$    DO jst = 1, nstm
!!$       DO ji = 1, kjpindex
!!$          IF (soiltile(ji,jst) .LT. 0.01) THEN
!!$             soiltile(ji,MAXLOC(soiltile(ji,:),1)) = soiltile(ji,MAXLOC(soiltile(ji,:),1)) + soiltile(ji,jst)
!!$             soiltile(ji,jst) = zero
!!$          ENDIF
!!$       ENDDO
!!$    ENDDO
    
  END SUBROUTINE slowproc_veget
 
 
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_lai
!!
!>\BRIEF        Do the interpolation of lai for the PFTs in case the laimap is not read   
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) Interplation by using the mean value of laimin and laimax for the PFTs    
!! (2) Interpolation between laimax and laimin values by using the temporal
!!  variations 
!! (3) If problem occurs during the interpolation, the routine stops 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_lai (kjpindex,lcanop,stempdiag,lalo,resolution,lai,mm,dd,read_lai,laimap)
    !
    ! 0. Declarations
    !
    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                          :: kjpindex   !! Domain size - terrestrial pixels only
    INTEGER(i_std), INTENT(in)                          :: lcanop     !! soil level used for LAI
    INTEGER(i_std), INTENT(in)                          :: mm, dd     !! Number of the month in the year and number of day of the month 
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)  :: stempdiag  !! Soil temperature (K) ???
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo       !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution !! Size in x an y of the grid (m) - surface area of the gridbox
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: laimap     !! map of lai read 
    LOGICAL, INTENT(in)                                 :: read_lai   !! Flag for the reading of laimap 

    !! 0.2 Output
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)   :: lai        !! PFT leaf area index (m^{2} m^{-2})LAI

    !! 0.4 Local
    INTEGER(i_std)                                      :: ji,jv      !! Local indices 
!_ ================================================================================================================================

    !
    IF  ( .NOT. read_lai ) THEN
    
       lai(: ,1) = zero
       ! On boucle sur 2,nvm au lieu de 1,nvm
       DO jv = 2,nvm
          SELECT CASE (type_of_lai(jv))
             
          CASE ("mean ")
             !
             ! 1. do the interpolation between laimax and laimin
             !
             lai(:,jv) = undemi * (llaimax(jv) + llaimin(jv))
             !
          CASE ("inter")
             !
             ! 2. do the interpolation between laimax and laimin
             !
             DO ji = 1,kjpindex
                lai(ji,jv) = llaimin(jv) + tempfunc(stempdiag(ji,lcanop)) * (llaimax(jv) - llaimin(jv))
             ENDDO
             !
          CASE default
             !
             ! 3. Problem
             !
             WRITE (numout,*) 'This kind of lai choice is not possible. '// &
                  ' We stop with type_of_lai ',jv,' = ', type_of_lai(jv) 
             STOP 'slowproc_lai'
             
          END SELECT
          
       ENDDO
       !
    ELSE
       lai(: ,1) = zero
       ! On boucle sur 2,nvm au lieu de 1,nvm
       DO jv = 2,nvm

          SELECT CASE (type_of_lai(jv))
             
          CASE ("mean ")
             !
             ! 1. force MAXVAL of laimap on lai on this PFT
             !
             DO ji = 1,kjpindex
                lai(ji,jv) = MAXVAL(laimap(ji,jv,:))
             ENDDO
             !
          CASE ("inter")
             !
             ! 2. do the interpolation between laimax and laimin
             !
             !
             ! If January
             !
             IF (mm .EQ. 1 ) THEN
                IF (dd .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,12)*(1-(dd+15)/30.) + laimap(:,jv,1)*((dd+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,1)*(1-(dd-15)/30.) + laimap(:,jv,2)*((dd-15)/30.)
                ENDIF
                !
                ! If December
                !
             ELSE IF (mm .EQ. 12) THEN
                IF (dd .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,11)*(1-(dd+15)/30.) + laimap(:,jv,12)*((dd+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,12)*(1-(dd-15)/30.) + laimap(:,jv,1)*((dd-15)/30.)
                ENDIF
          !
          ! ELSE
          !
             ELSE
                IF (dd .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,mm-1)*(1-(dd+15)/30.) + laimap(:,jv,mm)*((dd+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,mm)*(1-(dd-15)/30.) + laimap(:,jv,mm+1)*((dd-15)/30.)
                ENDIF
             ENDIF
             !
          CASE default
             !
             ! 3. Problem
             !
             WRITE (numout,*) 'This kind of lai choice is not possible. '// &
                  ' We stop with type_of_lai ',jv,' = ', type_of_lai(jv) 
             STOP 'slowproc_lai'
             
          END SELECT
          
       ENDDO
    ENDIF

  END SUBROUTINE slowproc_lai

! SHOULD BE DROPED (not commented) ???

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interlai_OLD
!!
!>\BRIEF         Interpolate the LAI map to the grid of the model 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::laimap
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_interlai_OLD(nbpt, lalo,  resolution, laimap)
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          !! Vector of latitude and longitudes (beware of the order !)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  laimap(nbpt,nvm,12)          !! lai read variable and re-dimensioned
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, ijml, i, j, ik, lml, tml, fid, ib, jb,ip, jp, vid, ai,iki,jkj
    REAL(r_std) :: lev(1), date, dt, coslat
    INTEGER(i_std) :: itau(1)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) ::  mask_lu
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_lu, lon_lu, mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_ful, lon_ful
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: laimaporig
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: laimap_lu
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lon_up, lon_low, lat_up, lat_low
    INTEGER, DIMENSION(nbpt) :: n_origlai
    INTEGER, DIMENSION(nbpt) :: n_found
    REAL(r_std), DIMENSION(nbpt,nvm) :: frac_origlai

    CHARACTER(LEN=80) :: meter
    REAL(r_std) :: prog, sumf
    LOGICAL :: found
    INTEGER :: idi,jdi, ilast, jlast, jj, ii, jv, inear, iprog
    REAL(r_std) :: domaine_lon_min, domaine_lon_max, domaine_lat_min, domaine_lat_max

!_ ================================================================================================================================
    !
    !
    !Config Key   = LAI_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = NOT(LAI_MAP)
    !Config Def   = lai2D.nc
    !Config Help  = The name of the file to be opened to read the LAI
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from a Nicolas VIOVY one. 
    !Config Units = [FILE]
    !
    filename = 'lai2D.nc'
    CALL getin_p('LAI_FILE',filename)
    !
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    !
    !
    ALLOCATE(lon_lu(iml))
    ALLOCATE(lat_lu(jml))
    ALLOCATE(laimap_lu(iml,jml,nvm,tml))
    ALLOCATE(mask_lu(iml,jml))
    !
    WRITE(numout,*) 'slowproc_interlai : Reading the LAI file'
    !
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
       CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
       CALL flinget(fid, 'LAI', iml, jml, nvm, tml, 1, 12, laimap_lu)
       CALL flinget(fid, 'mask', iml, jml, 0, 0, 0, 1, mask_lu)
       !
       CALL flinclo(fid)
    ENDIF
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    CALL bcast(laimap_lu)
    CALL bcast(mask_lu)
    !
    WRITE(numout,*) 'slowproc_interlai : ', lon_lu(1), lon_lu(iml),lat_lu(1), lat_lu(jml)
    !
    !
    ijml=iml*jml
    ALLOCATE(lon_ful(ijml))
    ALLOCATE(lat_ful(ijml))
    ALLOCATE(laimaporig(ijml,nvm,tml))
    ALLOCATE(mask(ijml))

    DO i=1,iml
       DO j=1,jml
          iki=(j-1)*iml+i
          lon_ful(iki)=lon_lu(i)
          lat_ful(iki)=lat_lu(j)
          laimaporig(iki,:,:)=laimap_lu(i,j,:,:)
          mask(iki)=mask_lu(i,j)
       ENDDO
    ENDDO
    !
    WHERE  ( laimaporig(:,:,:) .LT. 0 )
       laimaporig(:,:,:) = zero
    ENDWHERE
    !
    !
    ALLOCATE(lon_up(nbpt)) 
    ALLOCATE(lon_low(nbpt))
    ALLOCATE(lat_up(nbpt))
    ALLOCATE(lat_low(nbpt))
    !
    DO ib =1, nbpt
       !
       !  We find the 4 limits of the grid-box. As we transform the resolution of the model
       !  into longitudes and latitudes we do not have the problem of periodicity.
       !  coslat is a help variable here !
       !
       coslat = MAX(COS(lalo(ib,1) * pi/180. ), mincos )*pi/180. * R_Earth
       !
       lon_up(ib) = lalo(ib,2) + resolution(ib,1)/(2.0*coslat) 
       lon_low(ib) = lalo(ib,2) - resolution(ib,1)/(2.0*coslat) 
       !
       coslat = pi/180. * R_Earth
       !
       lat_up(ib) = lalo(ib,1) + resolution(ib,2)/(2.0*coslat) 
       lat_low(ib) = lalo(ib,1) - resolution(ib,2)/(2.0*coslat) 
       !
       !
       !
    ENDDO
    lon_up = NINT( lon_up * 1E6 ) * 1E-6
    lon_low = NINT( lon_low * 1E6 ) * 1E-6
    lat_up = NINT( lat_up * 1E6 ) * 1E-6
    lat_low = NINT( lat_low * 1E6 ) * 1E-6
    !
    !  Get the limits of the integration domaine so that we can speed up the calculations
    !
    domaine_lon_min = MINVAL(lon_low)
    domaine_lon_max = MAXVAL(lon_up)
    domaine_lat_min = MINVAL(lat_low)
    domaine_lat_max = MAXVAL(lat_up)
    !
!!$    WRITE(*,*) 'DOMAINE lon :', domaine_lon_min, domaine_lon_max
!!$    WRITE(*,*) 'DOMAINE lat :', domaine_lat_min, domaine_lat_max
    !
    ! Ensure that the fine grid covers the whole domain
    WHERE ( lon_ful(:) .LT. domaine_lon_min )
      lon_ful(:) = lon_ful(:) + 360.
    ENDWHERE
    !
    WHERE ( lon_ful(:) .GT. domaine_lon_max )
      lon_ful(:) = lon_ful(:) - 360.
    ENDWHERE
    lon_ful = NINT( lon_ful * 1E6 ) * 1E-6
    lat_ful = NINT( lat_ful * 1E6 ) * 1E-6
    !
    WRITE(numout,*) 'Interpolating the lai map :'
    WRITE(numout,'(2a40)')'0%--------------------------------------', &
                   & '------------------------------------100%'
    !
    ilast = 1
    n_origlai(:) = 0
    laimap(:,:,:) = zero   
    !
    DO ip=1,ijml
       !
       !   Give a progress meter
       !
       ! prog = ip/float(ijml)*79.
       !  IF ( ABS(prog - NINT(prog)) .LT. 1/float(ijml)*79. ) THEN
       !   meter(NINT(prog)+1:NINT(prog)+1) = 'x'
       !   WRITE(numout, advance="no", FMT='(a)') ACHAR(13)
       !   WRITE(numout, advance="no", FMT='(a80)') meter
       ! ENDIF
       iprog = NINT(float(ip)/float(ijml)*79.) - NINT(float(ip-1)/float(ijml)*79.)
       IF ( iprog .NE. 0 ) THEN
          WRITE(numout,'(a1,$)') 'y'
       ENDIF
       !
       !  Only start looking for its place in the smaler grid if we are within the domaine
       !  That should speed up things !
       !
       IF ( ( lon_ful(ip) .GE. domaine_lon_min ) .AND. &
            ( lon_ful(ip) .LE. domaine_lon_max ) .AND. &
            ( lat_ful(ip) .GE. domaine_lat_min ) .AND. &
            ( lat_ful(ip) .LE. domaine_lat_max )        ) THEN
          !
          ! look for point on GCM grid which this point on fine grid belongs to.
          ! First look at the point on the model grid where we arrived just before. There is 
          ! a good chace that neighbouring points on the fine grid fall into the same model
          ! grid box.
          !
          IF ( ( lon_ful(ip) .GE. lon_low(ilast) ) .AND. &
               ( lon_ful(ip) .LT. lon_up(ilast) ) .AND. &
               ( lat_ful(ip) .GE. lat_low(ilast) ) .AND. &
               ( lat_ful(ip) .LT. lat_up(ilast) )         ) THEN
             !
             ! We were lucky
             !
             IF (mask(ip) .GT. 0) THEN
                n_origlai(ilast) =  n_origlai(ilast) + 1  
                DO i=1,nvm
                   DO j=1,12   
                      laimap(ilast,i,j) = laimap(ilast,i,j) + laimaporig(ip,i,j)
                   ENDDO
                ENDDO
             ENDIF
             !
          ELSE
             !
             ! Otherwise, look everywhere.
             ! Begin close to last grid point.
             !
             found = .FALSE. 
             idi = 1
             !
             DO WHILE ( (idi .LT. nbpt) .AND. ( .NOT. found ) )

                !
                ! forward and backward
                !
                DO ii = 1,2
                   !
                   IF ( ii .EQ. 1 ) THEN
                      ib = ilast - idi
                   ELSE
                      ib = ilast + idi
                   ENDIF
                   !
                   IF ( ( ib .GE. 1 ) .AND. ( ib .LE. nbpt ) ) THEN 
                      IF ( ( lon_ful(ip) .GE. lon_low(ib) ) .AND. &
                           ( lon_ful(ip) .LT. lon_up(ib) ) .AND. &
                           ( lat_ful(ip) .GE. lat_low(ib) ) .AND. &
                           ( lat_ful(ip) .LT. lat_up(ib) )         ) THEN
                         !
                         IF (mask(ip) .gt. 0) THEN
                            DO i=1,nvm
                               DO j=1,12                
                                  laimap(ib,i,j) = laimap(ib,i,j) + laimaporig(ip,i,j) 
                               ENDDO
                            ENDDO
                            n_origlai(ib) =  n_origlai(ib) + 1
                         ENDIF
                         ilast = ib
                         found = .TRUE.
                         !
                      ENDIF
                   ENDIF
                   !
                ENDDO
                !
                idi = idi + 1
                !
             ENDDO
             !
          ENDIF ! lucky/not lucky
          !
       ENDIF     ! in the domain
    ENDDO


    ! determine fraction of LAI points in each box of the coarse grid
    DO ip=1,nbpt
       IF ( n_origlai(ip) .GT. 0 ) THEN
          DO jv =1,nvm
             laimap(ip,jv,:) = laimap(ip,jv,:)/REAL(n_origlai(ip),r_std)
          ENDDO
       ELSE
          !
          ! Now we need to handle some exceptions
          !
          IF ( lalo(ip,1) .LT. -56.0) THEN
             ! Antartica
             DO jv =1,nvm
                laimap(ip,jv,:) = zero
             ENDDO
             !
          ELSE IF ( lalo(ip,1) .GT. 70.0) THEN
             ! Artica
             DO jv =1,nvm
                laimap(ip,jv,:) = zero
             ENDDO
             !
          ELSE IF ( lalo(ip,1) .GT. 55.0 .AND. lalo(ip,2) .GT. -65.0 .AND. lalo(ip,2) .LT. -20.0) THEN
             ! Greenland
             DO jv =1,nvm
                laimap(ip,jv,:) = zero
             ENDDO
             !
          ELSE
             !
             WRITE(numout,*) 'PROBLEM, no point in the lai map found for this grid box'
             WRITE(numout,*) 'Longitude range : ', lon_low(ip), lon_up(ip)
             WRITE(numout,*) 'Latitude range : ', lat_low(ip), lat_up(ip)
             !
             WRITE(numout,*) 'Looking for nearest point on the lai map file'
             CALL slowproc_nearest (ijml, lon_ful, lat_ful, &
                  lalo(ip,2), lalo(ip,1), inear)
             WRITE(numout,*) 'Coordinates of the nearest point, ',inear,' :', &
                  lon_ful(inear),lat_ful(inear)
             !
             DO jv =1,nvm
                laimap(ip,jv,:) = laimaporig(inear,jv,:)
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    ! 
    WRITE(numout,*) 'slowproc_interlai : Interpolation Done'
    !
    !  
    !
    DEALLOCATE(lon_up)
    DEALLOCATE(lon_low)
    DEALLOCATE(lat_up)
    DEALLOCATE(lat_low)
    DEALLOCATE(lat_ful)
    DEALLOCATE(lon_ful)
    DEALLOCATE(lat_lu)
    DEALLOCATE(lon_lu)
    DEALLOCATE(laimap_lu)
    DEALLOCATE(laimaporig)
    DEALLOCATE(mask_lu)
    DEALLOCATE(mask)
    !
    RETURN
    !
  END SUBROUTINE slowproc_interlai_OLD

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interlai_NEW
!!
!>\BRIEF         Interpolate the LAI map to the grid of the model 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::laimap
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_interlai_NEW(nbpt, lalo,  resolution, neighbours, contfrac, laimap)
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          !! Vector of latitude and longitudes 
                                                                 !! (beware of the order = 1 : latitude, 2 : longitude)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    !
    INTEGER(i_std), INTENT(in)         :: neighbours(nbpt,8)     !! Vector of neighbours for each grid point 1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW) 
                            
    REAL(r_std), INTENT(in)             :: contfrac(nbpt)        !! Fraction of land in each grid box.
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  laimap(nbpt,nvm,12)          !! lai read variable and re-dimensioned
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename                               !! name of the LAI map read
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, jp, it, jj, jv
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_lu, lon_lu    !! latitude and
                                                                !! longitude, extract from LAI map
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat, lon        !! en 2D ???
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area     !! the area of the fine grid in the model grid ???
                                                                !! cf src_global/interpol_help.f90, line 377, called "areaoverlap"
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index !! the indexes from the grid boxes from the data that go into the 
                                                                !! model's boxes  
                                                                !! cf src_global/interpol_help.f90,line 300, called "ip"

    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: laimap_lu   !! value in LAIMAP
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: resol_lu
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    !
    REAL(r_std) :: coslat, lmax, lmin, ldelta
    INTEGER(i_std) :: nix, njx
    REAL(r_std) :: totarea
    INTEGER(i_std) :: idi, nbvmax                               !! nbvmax : number of maximum vegetation map
                                                                !! points in the GCM grid ; idi : its counter
    CHARACTER(LEN=30) :: callsign                               !! Allows to specify which variable is beeing treated
    !
    LOGICAL ::           renormelize_lai  ! flag to force LAI renormelization
    LOGICAL ::           ok_interpol                            !! optionnal return of aggregate_2d
    !
    INTEGER                  :: ALLOC_ERR 
!_ ================================================================================================================================
    !
    !Config Key   = LAI_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = LAI_MAP
    !Config Def   = lai2D.nc
    !Config Help  = The name of the file to be opened to read the LAI
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from a Nicolas VIOVY one. 
    !Config Units = [FILE]
    !
    filename = 'lai2D.nc'
    CALL getin_p('LAI_FILE',filename)
    !
    !
    !Config Key   = RENORM_LAI
    !Config Desc  = flag to force LAI renormelization
    !Config If    = LAI_MAP
    !Config Def   = n
    !Config Help  = If true, the laimap will be renormalize between llaimin and llaimax parameters.
    !Config Units = [FLAG]
    !
    renormelize_lai = .FALSE.
    CALL getin_p('RENORM_LAI',renormelize_lai)

    !
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    !
    !
    ALLOC_ERR=-1
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(laimap_lu(iml,jml,nvm,tml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of laimap_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(resol_lu(iml,jml,2), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of resol_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    !
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
       CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
       CALL flinget(fid, 'LAI', iml, jml, nvm, tml, 1, 12, laimap_lu)
       !
       WHERE (laimap_lu(:,:,:,:) < zero )
          laimap_lu(:,:,:,:) = zero
       ENDWHERE
       !
       CALL flinclo(fid)
    ENDIF
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    CALL bcast(laimap_lu)
    !
    ALLOC_ERR=-1
    ALLOCATE(lon(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lat(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    DO ip=1,iml
       lat(ip,:) = lat_lu(:)
    ENDDO
    DO jp=1,jml
       lon(:,jp) = lon_lu(:)
    ENDDO
    !
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    ! Consider all points a priori
    !
    mask(:,:) = 0
    !
    DO ip=1,iml
       DO jp=1,jml
          !
          ! Exclude the points where there is never a LAI value. It is probably 
          ! an ocean point.
          !
          IF ( ANY(laimap_lu(ip,jp,:,:) < 20.) ) THEN
             mask(ip,jp) = 1
          ENDIF
          !
          ! Resolution in longitude
          !
          coslat = MAX( COS( lat(ip,jp) * pi/180. ), mincos )     
          IF ( ip .EQ. 1 ) THEN
             resol_lu(ip,jp,1) = ABS( lon(ip+1,jp) - lon(ip,jp) ) * pi/180. * R_Earth * coslat
          ELSEIF ( ip .EQ. iml ) THEN
             resol_lu(ip,jp,1) = ABS( lon(ip,jp) - lon(ip-1,jp) ) * pi/180. * R_Earth * coslat
          ELSE
             resol_lu(ip,jp,1) = ABS( lon(ip+1,jp) - lon(ip-1,jp) )/2. * pi/180. * R_Earth * coslat
          ENDIF
          !
          ! Resolution in latitude
          !
          IF ( jp .EQ. 1 ) THEN
             resol_lu(ip,jp,2) = ABS( lat(ip,jp) - lat(ip,jp+1) ) * pi/180. * R_Earth
          ELSEIF ( jp .EQ. jml ) THEN
             resol_lu(ip,jp,2) = ABS( lat(ip,jp-1) - lat(ip,jp) ) * pi/180. * R_Earth
          ELSE
             resol_lu(ip,jp,2) =  ABS( lat(ip,jp-1) - lat(ip,jp+1) )/2. * pi/180. * R_Earth
          ENDIF
          !
       ENDDO
    ENDDO
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution_g(:,1))/MAXVAL(resol_lu(:,:,1)))+2
       njx=INT(MAXVAL(resolution_g(:,2))/MAXVAL(resol_lu(:,:,2)))+2
       nbvmax = nix*njx
    ENDIF
    CALL bcast(nbvmax)
    !
    callsign = 'LAI map'
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(nbpt, nbvmax, 2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_index(:,:,:)=0
       ALLOC_ERR=-1
       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon, lat, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
       ENDIF
       !
       nbvmax = nbvmax * 2
    ENDDO
    !
    laimap(:,:,:) = zero
    !
    DO ib=1,nbpt
       idi = COUNT(sub_area(ib,:) > zero)
       IF ( idi > 0 ) THEN
          totarea = zero
          DO jj=1,idi
             ip = sub_index(ib,jj,1)
             jp = sub_index(ib,jj,2)
             DO jv=1,nvm
                DO it=1,12
                   laimap(ib,jv,it) = laimap(ib,jv,it) + laimap_lu(ip,jp,jv,it)*sub_area(ib,jj)
                ENDDO
             ENDDO
             totarea = totarea + sub_area(ib,jj)
          ENDDO
          !
          ! Normalize
          !
          laimap(ib,:,:) = laimap(ib,:,:)/totarea
          !
       ELSE
          WRITE(numout,*) 'On point ', ib, ' no points where found for interpolating the LAI map.'
          WRITE(numout,*) 'Location : ', lalo(ib,2), lalo(ib,1)
          DO jv=1,nvm
             laimap(ib,jv,:) = (llaimax(jv)+llaimin(jv))/deux
          ENDDO
          WRITE(numout,*) 'Solved by putting the average LAI for the PFT all year long'
       ENDIF
    ENDDO
    !
    ! Normelize the read LAI by the values SECHIBA is used to
    !
    IF ( renormelize_lai ) THEN
       DO ib=1,nbpt
          DO jv=1,nvm
             lmax = MAXVAL(laimap(ib,jv,:))
             lmin = MINVAL(laimap(ib,jv,:))
             ldelta = lmax-lmin
             IF ( ldelta < min_sechiba) THEN
                ! LAI constante ... keep it constant
                laimap(ib,jv,:) = (laimap(ib,jv,:)-lmin)+(llaimax(jv)+llaimin(jv))/deux
             ELSE
                laimap(ib,jv,:) = (laimap(ib,jv,:)-lmin)/(lmax-lmin)*(llaimax(jv)-llaimin(jv))+llaimin(jv)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
    WRITE(numout,*) 'slowproc_interlai : Interpolation Done'
    !
    !  
    !
    DEALLOCATE(lat_lu)
    DEALLOCATE(lon_lu)
    DEALLOCATE(lon)
    DEALLOCATE(lat)
    DEALLOCATE(laimap_lu)
    DEALLOCATE(mask)
    DEALLOCATE(sub_area)
    DEALLOCATE(sub_index)
    DEALLOCATE (resol_lu)
    !
    RETURN
    !
  END SUBROUTINE slowproc_interlai_NEW

! NOT COMMENTED YET!

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_update
!!
!>\BRIEF          Interpolate a vegetation map (by pft) NEED TO BE DONE!!!!!!!!
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_update(nbpt, lalo, neighbours,  resolution, contfrac, & 
       &       veget_last, frac_nobio_last, & 
       &       veget_next, frac_nobio_next, veget_year, init)
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs 
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,8), INTENT(in)         :: neighbours       !! Vector of neighbours for each grid point 
                                                                              !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    REAL(r_std), DIMENSION(nbpt,nvm), INTENT(in)           :: veget_last      !! old max vegetfrac
    REAL(r_std), DIMENSION(nbpt,nnobio), INTENT(in)        :: frac_nobio_last !! old fraction of the mesh which is 
                                                                              !! covered by ice, lakes, ...
    !
    INTEGER(i_std), INTENT(in)         :: veget_year            !! first year for landuse (0 == NO TIME AXIS)
    LOGICAL, OPTIONAL, INTENT(in)      :: init                  !! initialisation : in case of dgvm, it forces update of all PFTs
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt,nvm), INTENT(out)          :: veget_next       !! new max vegetfrac
    REAL(r_std), DIMENSION(nbpt,nnobio), INTENT(out)       :: frac_nobio_next  !! new fraction of the mesh which is 
    !! covered by ice, lakes, ...
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, jp, inobio, jv
    INTEGER(i_std) :: nb_coord,nb_var, nb_gat,nb_dim
    REAL(r_std) :: date, dt
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:)  :: itau
    REAL(r_std), DIMENSION(1)  :: time_counter
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_lu, lon_lu
    INTEGER,DIMENSION(flio_max_var_dims) :: l_d_w, i_d_w
    LOGICAL :: exv, l_ex
    !
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: vegmap            ! last coord is time with only one value = 1
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: vegmap_1            ! last coord is time with only one value = 1  (IF VEGET_YEAR == 0 , NO TIME AXIS)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_ful, lon_ful
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    !
    REAL(r_std) :: sumf, err, norm
    REAL(r_std) :: totarea
    INTEGER(i_std) :: idi, nbvmax
    CHARACTER(LEN=30) :: callsign
    !
    LOGICAL ::           ok_interpol ! optionnal return of aggregate_2d
    !
    ! for DGVM case :
    REAL(r_std)                 :: sum_veg                     ! sum of vegets
    REAL(r_std)                 :: sum_nobio                   ! sum of nobios
    REAL(r_std)                 :: sumvAnthro_old, sumvAnthro  ! last an new sum of antrhopic vegets
    REAL(r_std)                 :: rapport                     ! (S-B) / (S-A)
    LOGICAL                    :: partial_update              ! if TRUE, partialy update PFT (only anthropic ones) 
                                                              ! e.g. in case of DGVM and not init (optional parameter)
    !
    LOGICAL, PARAMETER         :: debug = .FALSE.
    !
    INTEGER                  :: ALLOC_ERR

!_ ================================================================================================================================
    !
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = LAND_USE
    !Config Def   = PFTmap.nc
    !Config Help  = The name of the file to be opened to read a vegetation
    !Config         map (in pft) is to be given here. 
    !Config Units = [FILE]
    !
    filename = 'PFTmap.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    !
    IF (is_root_prc) THEN
       IF (debug) THEN
          WRITE(numout,*) "Entering slowproc_update. Debug mode."
          WRITE (*,'(/," --> fliodmpf")')
          CALL fliodmpf (TRIM(filename))
          WRITE (*,'(/," --> flioopfd")')
       ENDIF
       CALL flioopfd (TRIM(filename),fid,nb_dim=nb_coord,nb_var=nb_var,nb_gat=nb_gat)
       IF (debug) THEN
          WRITE (*,'(" Number of coordinate        in the file : ",I2)') nb_coord
          WRITE (*,'(" Number of variables         in the file : ",I2)') nb_var
          WRITE (*,'(" Number of global attributes in the file : ",I2)') nb_gat
       ENDIF
    ENDIF
    CALL bcast(nb_coord)
    CALL bcast(nb_var)
    CALL bcast(nb_gat)
    IF ( veget_year > 0 ) THEN
       IF (is_root_prc) &
         CALL flioinqv (fid,v_n="time_counter",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
       CALL bcast(l_d_w)
       tml=l_d_w(1)
    ELSE
       tml=0
    END IF
    WRITE(numout,*) "tml =",tml

    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lon",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    iml=l_d_w(1)
    WRITE(numout,*) "iml =",iml
    
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="lat",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    jml=l_d_w(1)
    WRITE(numout,*) "jml =",jml
    
    IF (is_root_prc) &
         CALL flioinqv (fid,v_n="veget",l_ex=l_ex,nb_dims=nb_dim,len_dims=l_d_w) 
    CALL bcast(l_d_w)
    lml=l_d_w(1)

    IF (lml /= nvm) &
         CALL ipslerr_p (3,'slowproc_update', &
               &          'Problem with vegetation file for Land Use','lml /= nvm', &
               &          '(number of pft must be equal)')
    !
    ALLOC_ERR=-1
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_lu : ",ALLOC_ERR
      STOP 
    ENDIF

    IF ( veget_year > 0 ) THEN
       IF (tml > 0) THEN
          ALLOC_ERR=-1
          ALLOCATE(itau(tml), STAT=ALLOC_ERR)
          IF (ALLOC_ERR/=0) THEN
             WRITE(numout,*) "ERROR IN ALLOCATION of itau : ",ALLOC_ERR
             STOP 
          ENDIF
       ENDIF
    END IF
       
    IF (is_root_prc) THEN
       IF (tml > 0) THEN
          CALL fliogstc (fid, t_axis=itau,x_axis=lon_lu,y_axis=lat_lu)
          IF (veget_year <= tml) THEN
             CALL fliogetv (fid,"time_counter",time_counter, start=(/ veget_year /), count=(/ 1 /))
             WRITE(numout,*) "slowproc_update LAND_USE : the date read for vegetmax is ",time_counter
          ELSE
             CALL fliogetv (fid,"time_counter",time_counter, start=(/ tml /), count=(/ 1 /))
             WRITE(numout,*) "slowproc_update LAND_USE : You try to update vegetmax with a the date greater than in the file."
             WRITE(numout,*) "                           We will keep the last one :",time_counter
          ENDIF
       ELSE
          CALL fliogstc (fid, x_axis=lon_lu,y_axis=lat_lu)
       ENDIF
    ENDIF
    
    IF (tml > 0) THEN
       CALL bcast(itau)
    ENDIF
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    !
    ALLOC_ERR=-1
    ALLOCATE(lat_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_ful : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_ful(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_ful : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    DO ip=1,iml
       lon_ful(ip,:)=lon_lu(ip)
    ENDDO
    DO jp=1,jml
       lat_ful(:,jp)=lat_lu(jp)
    ENDDO
    !
    !
    WRITE(numout,*) 'Reading the LAND USE vegetation file'
    !
    ALLOC_ERR=-1
    ALLOCATE(vegmap(iml,jml,nvm,1), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of vegmap : ",ALLOC_ERR
      STOP 
    ENDIF
    IF ( veget_year == 0 ) THEN
       IF (is_root_prc) THEN
          ALLOC_ERR=-1
          ALLOCATE(vegmap_1(iml,jml,nvm), STAT=ALLOC_ERR)
          IF (ALLOC_ERR/=0) THEN
             WRITE(numout,*) "ERROR IN ALLOCATION of vegmap_1 : ",ALLOC_ERR
             STOP 
          ENDIF
       ENDIF
    ENDIF
    !
!!!$    CALL flinopen &
!!!$       &  (filename, .FALSE., iml, jml, lml, lon_ful, lat_ful, &
!!!$       &   lev_ful, tml, itau, date, dt, fid)
!!!=> FATAL ERROR FROM ROUTINE flinopen
!!! --> No time axis found
 !!
!!!MM modif TAG 1.4 : 
!!!    CALL flinget(fid, 'lon', iml, 0, 0, 0, 1, 1, lon_lu)
!!!    CALL flinget(fid, 'lat', jml, 0, 0, 0, 1, 1, lat_lu)
!!!    CALL flinget(fid, 'maxvegetfrac', iml, jml, nvm, tml, 1, 293, vegmap_lu)
!!!FATAL ERROR FROM ROUTINE flinopen
!!! --> No variable lon
!!!   We get only the right year
!!!$    CALL flinget(fid, 'maxvegetfrac', iml, jml, nvm, tml, veget_year, veget_year, vegmap)
!!!$    !
!!!$    CALL flinclo(fid)

    IF (is_root_prc) &
         CALL flioinqv (fid,"maxvegetfrac",exv,nb_dims=nb_dim,len_dims=l_d_w,id_dims=i_d_w)
    CALL bcast(exv)
    CALL bcast(nb_dim)
    CALL bcast(l_d_w)
    CALL bcast(i_d_w)

    IF (exv) THEN
       IF (debug) THEN
          WRITE (numout,'(" Number of dimensions : ",I2)') nb_dim
          WRITE (numout,'(" Dimensions :",/,5(1X,I7,:))') l_d_w(1:nb_dim)
          WRITE (numout,'(" Identifiers :",/,5(1X,I7,:))') i_d_w(1:nb_dim)
       ENDIF
       !
       IF ( veget_year > 0 ) THEN
          IF (is_root_prc) THEN
             IF (veget_year <= tml) THEN
                CALL fliogetv (fid,"maxvegetfrac", vegmap, start=(/ 1, 1, 1, veget_year /), count=(/ iml, jml, nvm, 1 /))
             ELSE
                CALL fliogetv (fid,"maxvegetfrac", vegmap, start=(/ 1, 1, 1, tml /), count=(/ iml, jml, nvm, 1 /))
             ENDIF
          ENDIF
       ELSE
          IF (is_root_prc) THEN
             CALL fliogetv (fid,"maxvegetfrac", vegmap_1, start=(/ 1, 1, 1 /), count=(/ iml, jml, nvm /))
             vegmap(:,:,:,1)=vegmap_1(:,:,:)
             DEALLOCATE(vegmap_1)
          ENDIF
       ENDIF
       CALL bcast(vegmap)
       IF (is_root_prc) CALL flioclo (fid)
    ELSE
       CALL ipslerr_p (3,'slowproc_update', &
            &          'Problem with vegetation file for Land Use.', &
            &          "Variable maxvegetfrac doesn't exist.", &
            &          '(verify your land use file.)')
    ENDIF
    !
    ! Mask of permitted variables.
    !
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    mask(:,:) = 0
    DO ip=1,iml
       DO jp=1,jml
          sum_veg=SUM(vegmap(ip,jp,:,1))
          IF ( sum_veg .GE. min_sechiba .AND. sum_veg .LE. 1.-1.e-7) THEN
             mask(ip,jp) = 1
             IF (debug) THEN
                WRITE(numout,*) "update : SUM(vegmap(",ip,jp,")) = ",sum_veg
             ENDIF
          ELSEIF ( sum_veg .GT. 1.-1.e-7 .AND. sum_veg .LE. 2.) THEN
             ! normalization
             vegmap(ip,jp,:,1) = vegmap(ip,jp,:,1) / sum_veg
             mask(ip,jp) = 1
             IF (debug) THEN
                WRITE(numout,*) "update : SUM(vegmap(",ip,jp,"))_c = ",SUM(vegmap(ip,jp,:,1))
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    !
    ! The number of maximum vegetation map points in the GCM grid should
    ! also be computed and not imposed here.
    !
    nbvmax = 200
    !
    callsign="Land Use Vegetation map"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(nbpt, nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_index(:,:,:)=0

       ALLOC_ERR=-1
       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon_ful, lat_ful, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
       ENDIF
       !
       nbvmax = nbvmax * 2
    ENDDO
    !
    ! Compute the logical for partial (only anthropic) PTFs update
    IF (PRESENT(init)) THEN
       partial_update = control%ok_dgvm  .AND. .NOT. init
    ELSE
       partial_update = control%ok_dgvm  
    ENDIF
    !
    IF ( .NOT. partial_update ) THEN
       !
       veget_next(:,:)=zero
       !
       DO ib = 1, nbpt
          sumf=zero
          DO idi=1, nbvmax
             ! Leave the do loop if all sub areas are treated, sub_area <= 0
             IF ( sub_area(ib,idi) <= zero ) EXIT
             ip = sub_index(ib,idi,1)
             jp = sub_index(ib,idi,2)
             veget_next(ib,:) = veget_next(ib,:) + vegmap(ip,jp,:,1)*sub_area(ib,idi)
             sumf=sumf + sub_area(ib,idi)
          ENDDO
!!!$          !
!!!$          !  Limit the smalest vegetation fraction to 0.5%
!!!$          !
!!!$          DO jv = 1, nvm
!!!$             IF ( veget_next(ib,jv) .LT. min_vegfrac ) THEN
!!!$                veget_next(ib,jv) = zero
!!!$             ENDIF
!!!$          ENDDO
          !
          ! Normalize
          !
          IF (sumf > min_sechiba) THEN
             veget_next(ib,:) = veget_next(ib,:) / sumf
          ELSE
             WRITE(numout,*) "slowproc_update : No land point in the map for point ",&
                  ib,",(",lalo(ib,1),",",lalo(ib,2),")" 
             CALL ipslerr_p (2,'slowproc_update', &
                  &          'Problem with vegetation file for Land Use.', &
                  &          "No land point in the map for point", &
                  &          'Keep old values. (verify your land use file.)')
!!$             CALL slowproc_nearest (iml, lon_ful, lat_ful, &
!!$                  lalo(ib,2), lalo(ib,1), inear)
             IF (PRESENT(init)) THEN
                IF (init) THEN
                    veget_next(ib,1) = un
                    veget_next(ib,2:nvm) = zero
                ELSE
                   veget_next(ib,:) = veget_last(ib,:)
                ENDIF
             ELSE
                veget_next(ib,:) = veget_last(ib,:)
             ENDIF
          ENDIF
          !
          IF (debug) THEN
             WRITE(numout,*) "SUM(veget_next(",ib,")) = ",SUM(veget_next(ib,:))
          ENDIF
       ENDDO
    ELSE
       DO ib = 1, nbpt
          ! last veget for this point
          sum_veg=SUM(veget_last(ib,:))
          IF (debug) THEN
             WRITE(numout,*) "SUM(veget_last(",ib,")) = ",sum_veg
          ENDIF
          !
          ! If the DGVM is activated, only anthropiques PFT are utpdated,
          ! other are copied 
          veget_next(ib,:) = veget_last(ib,:)
          !
          ! natural ones are initialized to zero.
          DO jv = 2, nvm
             ! If the DGVM is activated, only anthropiques PFT are utpdated 
             IF ( .NOT. natural(jv) ) THEN
                veget_next(ib,jv) = zero
             ENDIF
          ENDDO
          !
          sumf=zero
          DO idi=1, nbvmax
             ! Leave the do loop if all sub areas are treated, sub_area <= 0
             IF ( sub_area(ib,idi) <= zero ) EXIT
             ip = sub_index(ib,idi,1)
             jp = sub_index(ib,idi,2)
             ! If the DGVM is activated, only anthropic PFTs are utpdated 
             DO jv = 2, nvm
                IF ( .NOT. natural(jv) ) THEN       
                   veget_next(ib,jv) = veget_next(ib,jv) + vegmap(ip,jp,jv,1)*sub_area(ib,idi)
                ENDIF
             ENDDO
             sumf=sumf + sub_area(ib,idi)
          ENDDO
!!!$          !
!!!$          !  Limit the smalest vegetation fraction to 0.5%
!!!$          !
!!!$          DO jv = 2, nvm
!!!$             ! On anthropic and natural PFTs ? 
!!!$             IF ( veget_next(ib,jv) .LT. min_vegfrac ) THEN
!!!$                veget_next(ib,jv) = zero
!!!$             ENDIF
!!!$          ENDDO
          !
          ! Normalize
          !
          ! Proposition de Pierre :
          ! apres modification de la surface des PFTs anthropiques,
          ! on doit conserver la proportion des PFTs naturels.
          ! ie la somme des vegets est conservee
          !    et PFT naturel / (somme des vegets - somme des vegets anthropiques)
          !       est conservee.
          ! Modification de Nathalie : 
          ! Si les PFTs anthropique diminue, on les remplace plutôt par du sol nu.
          ! Le DGVM est chargé de ré-introduire les PFTs naturels.
          IF (sumf > min_sechiba) THEN
             sumvAnthro_old = zero
             sumvAnthro     = zero
             DO jv = 2, nvm
                IF ( .NOT. natural(jv) ) THEN
                   veget_next(ib,jv) = veget_next(ib,jv) / sumf
                   sumvAnthro = sumvAnthro + veget_next(ib,jv)
                   sumvAnthro_old = sumvAnthro_old + veget_last(ib,jv)
                ENDIF
             ENDDO

             IF ( sumvAnthro_old < sumvAnthro ) THEN
                ! deforestation
                ! conservation :
                rapport = ( sum_veg - sumvAnthro ) / ( sum_veg - sumvAnthro_old )
                DO jv = 1, nvm
                   IF ( natural(jv) ) THEN
                      veget_next(ib,jv) = veget_last(ib,jv) * rapport
                   ENDIF
                ENDDO
             ELSE
                ! reforestation
                DO jv = 1, nvm
                   IF ( natural(jv) ) THEN
                      veget_next(ib,jv) = veget_last(ib,jv)
                   ENDIF
                ENDDO
                veget_next(ib,1) = veget_next(ib,1) + sumvAnthro_old - sumvAnthro
             ENDIF

             ! test
             IF ( ABS( SUM(veget_next(ib,:)) - sum_veg ) > 10*EPSILON(un) ) THEN
                WRITE(numout,*) "No conservation of sum of veget for point ",ib,",(",lalo(ib,1),",",lalo(ib,2),")" 
                WRITE(numout,*) "last sum of veget ",sum_veg," new sum of veget ",SUM(veget_next(ib,:))," error : ",&
                     &                         SUM(veget_next(ib,:)) - sum_veg
                WRITE(numout,*) "Anthropic modifications : last ",sumvAnthro_old," new ",sumvAnthro     
                CALL ipslerr_p (3,'slowproc_update', &
                     &          'No conservation of sum of veget_next', &
                     &          "The sum of veget_next is different after reading Land Use map.", &
                     &          '(verify the dgvm case model.)')
             ENDIF
          ELSE
             WRITE(numout,*) "No land point in the map for point ",ib,",(",lalo(ib,1),",",lalo(ib,2),")" 
!             CALL ipslerr_p (3,'slowproc_update', &
             CALL ipslerr_p (2,'slowproc_update', &
                  &          'Problem with vegetation file for Land Use.', &
                  &          "No land point in the map for point", &
                  &          '(verify your land use file.)')
             veget_next(ib,:) = veget_last(ib,:)
          ENDIF
          
       ENDDO       
    ENDIF
    !
    frac_nobio_next (:,:) = un
    !
!MM
    ! Work only for one nnobio !! (ie ice)
    DO inobio=1,nnobio
       DO jv=1,nvm
          !
          DO ib = 1, nbpt
             frac_nobio_next(ib,inobio) = frac_nobio_next(ib,inobio) - veget_next(ib,jv)
          ENDDO
       ENDDO
       !
    ENDDO
    !
    DO ib = 1, nbpt
       sum_veg = SUM(veget_next(ib,:))
       sum_nobio = SUM(frac_nobio_next(ib,:))
       IF (sum_nobio < 0.) THEN
          frac_nobio_next(ib,:) = zero
          veget_next(ib,1) = veget_next(ib,1) - sum_nobio
          sum_veg = SUM(veget_next(ib,:))
       ENDIF
       sumf = sum_veg + sum_nobio
       IF (sumf > min_sechiba) THEN
          veget_next(ib,:) = veget_next(ib,:) / sumf
          frac_nobio_next(ib,:) = frac_nobio_next(ib,:) / sumf
          norm=SUM(veget_next(ib,:))+SUM(frac_nobio_next(ib,:))
          err=norm-un
          IF (debug) &
             WRITE(numout,*) "ib ",ib," SUM(veget_next(ib,:)+frac_nobio_next(ib,:))-un, sumf",err,sumf
          IF (abs(err) > -EPSILON(un)) THEN
!MM 1.9.3
!          IF (abs(err) > 0.) THEN
             IF ( SUM(frac_nobio_next(ib,:)) > min_sechiba ) THEN
                frac_nobio_next(ib,1) = frac_nobio_next(ib,1) - err
             ELSE
                veget_next(ib,1) = veget_next(ib,1) - err
             ENDIF
             norm=SUM(veget_next(ib,:))+SUM(frac_nobio_next(ib,:))
             err=norm-un
             IF (debug) &
                  WRITE(numout,*) "ib ",ib," SUM(veget_next(ib,:)+frac_nobio_next(ib,:))-un",err
             IF (abs(err) > EPSILON(un)) THEN
!MM 1.9.3
!             IF (abs(err) > 0.) THEN
                WRITE(numout,*) "update : Problem with point ",ib,",(",lalo(ib,1),",",lalo(ib,2),")" 
                WRITE(numout,*) "         err(sum-1.) = ",abs(err)
                CALL ipslerr_p (2,'slowproc_update', &
                     &          'Problem with sum vegetation + sum fracnobio for Land Use.', &
                     &          "sum not equal to 1.", &
                     &          '(verify your land use file.)')
             ENDIF
          ENDIF
       ELSE
          WRITE(numout,*) "No vegetation nor frac_nobio for point ",ib,",(",lalo(ib,1),",",lalo(ib,2),")" 
          WRITE(numout,*)"Replaced by bare_soil !! "
          veget_next(ib,1) = un
          veget_next(ib,2:nvm) = zero
          frac_nobio_next(ib,:) = zero
!!!$          CALL ipslerr_p (3,'slowproc_update', &
!!!$               &          'Problem with vegetation file for Land Use.', &
!!!$               &          "No vegetation nor frac_nobio for point ", &
!!!$               &          '(verify your land use file.)')
       ENDIF
    ENDDO
    !
    WRITE(numout,*) 'slowproc_update : Interpolation Done'
    !
    DEALLOCATE(vegmap)
    DEALLOCATE(lat_lu,lon_lu)
    DEALLOCATE(lat_ful,lon_ful)
    DEALLOCATE(mask)
    DEALLOCATE(sub_index,sub_area)
    !
    RETURN
    !
  END SUBROUTINE slowproc_update


! NEED TO GET RID OF SLOWPROC_INTERPOL_OLD because it is never called!
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interpol_OLD
!!
!>\BRIEF         Interpolate the IGBP vegetation map to the grid of the model
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  !MM TAG 1.6 model !
  SUBROUTINE slowproc_interpol_OLD(nbpt, lalo, neighbours, resolution, veget, frac_nobio )
  !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,8)    !! Vector of neighbours for each grid point 
                                                                 !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  veget(nbpt,nvm)         !! Vegetation fractions
    REAL(r_std), INTENT(out)    ::  frac_nobio(nbpt,nnobio) !! Fraction of the mesh which is covered by ice, lakes, ...
    !
    !  0.3 LOCAL
    !
    INTEGER(i_std), PARAMETER                       :: nolson = 94  !! Number of Olson classes
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, jp, vid
    REAL(r_std) :: lev(1), date, dt, coslat
    INTEGER(i_std) :: itau(1)
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_ful, lon_ful, vegmap
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lon_up, lon_low, lat_up, lat_low
    INTEGER, DIMENSION(nbpt,nolson) :: n_origveg
    INTEGER, DIMENSION(nbpt) :: n_found
    REAL(r_std), DIMENSION(nbpt,nolson) :: frac_origveg
    REAL(r_std) :: vegcorr(nolson,nvm)
    REAL(r_std) :: nobiocorr(nolson,nnobio)
    CHARACTER(LEN=80) :: meter
    REAL(r_std) :: prog, sumf
    LOGICAL :: found
    INTEGER :: idi, ilast, ii, jv, inear, iprog
    REAL(r_std) :: domaine_lon_min, domaine_lon_max, domaine_lat_min, domaine_lat_max
!_ ================================================================================================================================
    !
    CALL get_vegcorr (nolson,vegcorr,nobiocorr)
    !
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = NOT(IMPOSE_VEG)
    !Config Def   = carteveg5km.nc
    !Config Help  = The name of the file to be opened to read the vegetation
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from the IGBP one. We assume that we have
    !Config         a classification in 87 types. This is Olson modified by Viovy.
    !Config Units = [FILE]
    !
    filename = 'carteveg5km.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    !
    if (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    !
    !
    ALLOCATE(lat_ful(iml))
    ALLOCATE(lon_ful(iml))
    ALLOCATE(vegmap(iml))
    !
    WRITE(numout,*) 'Reading the vegetation file'
    !
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, jml, lml, tml, 1, 1, lon_ful)
       CALL flinget(fid, 'latitude', iml, jml, lml, tml, 1, 1, lat_ful)
       CALL flinget(fid, 'vegetation_map', iml, jml, lml, tml, 1, 1, vegmap)
       !
       CALL flinclo(fid)
    ENDIF
    
    CALL bcast(lon_ful)
    CALL bcast(lat_ful)
    CALL bcast(vegmap)
    
    !
    IF (MAXVAL(vegmap) .LT. nolson) THEN
       WRITE(numout,*) 'WARNING -- WARNING'
       WRITE(numout,*) 'The vegetation map has to few vegetation types.'
       WRITE(numout,*) 'If you are lucky it will work but please check'
    ELSE IF ( MAXVAL(vegmap) .GT. nolson) THEN
       WRITE(numout,*) 'More vegetation types in file than the code can'
       WRITE(numout,*) 'deal with.: ',  MAXVAL(vegmap),  nolson
       STOP 'slowproc_interpol'
    ENDIF
    !
    ALLOCATE(lon_up(nbpt)) 
    ALLOCATE(lon_low(nbpt))
    ALLOCATE(lat_up(nbpt))
    ALLOCATE(lat_low(nbpt))
    !
    DO ib =1, nbpt
       !
       !  We find the 4 limits of the grid-box. As we transform the resolution of the model
       !  into longitudes and latitudes we do not have the problem of periodicity.
       !  coslat is a help variable here !
       !
       coslat = MAX(COS(lalo(ib,1) * pi/180. ), mincos )*pi/180. * R_Earth
       !
       lon_up(ib) = lalo(ib,2) + resolution(ib,1)/(2.0*coslat) 
       lon_low(ib) = lalo(ib,2) - resolution(ib,1)/(2.0*coslat) 
       !
       coslat = pi/180. * R_Earth
       !
       lat_up(ib) = lalo(ib,1) + resolution(ib,2)/(2.0*coslat) 
       lat_low(ib) = lalo(ib,1) - resolution(ib,2)/(2.0*coslat) 
       !
       !
       veget(ib,:) = zero
       frac_nobio (ib,:) = zero
       !
    ENDDO
    !
    !  Get the limits of the integration domaine so that we can speed up the calculations
    !
    domaine_lon_min = MINVAL(lon_low)
    domaine_lon_max = MAXVAL(lon_up)
    domaine_lat_min = MINVAL(lat_low)
    domaine_lat_max = MAXVAL(lat_up)
    !
!!$    WRITE(*,*) 'DOMAINE lon :', domaine_lon_min, domaine_lon_max
!!$    WRITE(*,*) 'DOMAINE lat :', domaine_lat_min, domaine_lat_max
    !
    ! Ensure that the fine grid covers the whole domain
    WHERE ( lon_ful(:) .LT. domaine_lon_min )
      lon_ful(:) = lon_ful(:) + 360.
    ENDWHERE
    !
    WHERE ( lon_ful(:) .GT. domaine_lon_max )
      lon_ful(:) = lon_ful(:) - 360.
    ENDWHERE
    !
    WRITE(numout,*) 'Interpolating the vegetation map :'
    WRITE(numout,'(2a40)')'0%--------------------------------------', &
                   & '------------------------------------100%'
    !
    ilast = 1
    n_origveg(:,:) = 0
    !
    DO ip=1,iml
      !
      !   Give a progress meter
      !
      ! prog = ip/float(iml)*79.
      !  IF ( ABS(prog - NINT(prog)) .LT. 1/float(iml)*79. ) THEN
      !   meter(NINT(prog)+1:NINT(prog)+1) = 'x'
      !   WRITE(numout, advance="no", FMT='(a)') ACHAR(13)
      !   WRITE(numout, advance="no", FMT='(a80)') meter
      ! ENDIF
      iprog = NINT(float(ip)/float(iml)*79.) - NINT(float(ip-1)/float(iml)*79.)
      IF ( iprog .NE. 0 ) THEN
        WRITE(numout,'(a1,$)') 'x'
      ENDIF
      !
      !  Only start looking for its place in the smaler grid if we are within the domaine
      !  That should speed up things !
      !
      IF ( ( lon_ful(ip) .GE. domaine_lon_min ) .AND. &
           ( lon_ful(ip) .LE. domaine_lon_max ) .AND. &
           ( lat_ful(ip) .GE. domaine_lat_min ) .AND. &
           ( lat_ful(ip) .LE. domaine_lat_max )        ) THEN
        !
        ! look for point on GCM grid which this point on fine grid belongs to.
        ! First look at the point on the model grid where we arrived just before. There is 
        ! a good chace that neighbouring points on the fine grid fall into the same model
        ! grid box.
        !
        !
        ! THERE IS A BUG HERE !!! IF THE GCM GRID SITS ON THE DATE LINE WE WILL HAVE FOR INSTANCE
        ! LON_LOW = -182 AND LON_UP = -178. THUS WE WILL ONLY PICK UP HALF THE POINTS NEEDED.
        !
        IF ( ( lon_ful(ip) .GT. lon_low(ilast) ) .AND. &
             ( lon_ful(ip) .LT. lon_up(ilast) ) .AND. &
             ( lat_ful(ip) .GT. lat_low(ilast) ) .AND. &
             ( lat_ful(ip) .LT. lat_up(ilast) )         ) THEN
          !
          ! We were lucky
          !
          n_origveg(ilast,NINT(vegmap(ip))) = n_origveg(ilast,NINT(vegmap(ip))) + 1
          !
        ELSE
          !
          ! Otherwise, look everywhere.
          ! Begin close to last grid point.
          !
          found = .FALSE. 
          idi = 1
          !
          DO WHILE ( (idi .LT. nbpt) .AND. ( .NOT. found ) )
            !
            ! forward and backward
            !
            DO ii = 1,2
              !
              IF ( ii .EQ. 1 ) THEN
                ib = ilast - idi
              ELSE
                ib = ilast + idi
              ENDIF
              !
              IF ( ( ib .GE. 1 ) .AND. ( ib .LE. nbpt ) ) THEN
                IF ( ( lon_ful(ip) .GT. lon_low(ib) ) .AND. &
                     ( lon_ful(ip) .LT. lon_up(ib) ) .AND. &
                     ( lat_ful(ip) .GT. lat_low(ib) ) .AND. &
                     ( lat_ful(ip) .LT. lat_up(ib) )         ) THEN
                  !
                  n_origveg(ib,NINT(vegmap(ip))) = n_origveg(ib,NINT(vegmap(ip))) + 1
                  ilast = ib
                  found = .TRUE.
                  !
                ENDIF
              ENDIF
              !
            ENDDO
            !
            idi = idi + 1
            !
          ENDDO
          !
        ENDIF ! lucky/not lucky
        !
      ENDIF     ! in the domain
    ENDDO

    !
    ! Now we know how many points of which Olson type from the fine grid fall
    ! into each box of the (coarse) model grid: n_origveg(nbpt,nolson)
    !

    !
    ! determine number of points of the fine grid which fall into each box of the 
    ! coarse grid
    !
    DO ib = 1, nbpt
      n_found(ib) = SUM( n_origveg(ib,:) )
    ENDDO

    !
    ! determine fraction of Olson vegetation type in each box of the coarse grid
    !
    DO vid = 1, nolson
      WHERE ( n_found(:) .GT. 0 ) 
        frac_origveg(:,vid) =  REAL(n_origveg(:,vid),r_std) /  REAL(n_found(:),r_std)
      ELSEWHERE
         frac_origveg(:,vid) = zero
      ENDWHERE
    ENDDO

    !
    ! now finally calculate coarse vegetation map
    ! Find which model vegetation corresponds to each Olson type 
    !
    DO vid = 1, nolson
      !
      DO jv = 1, nvm
        veget(:,jv) = veget(:,jv) + frac_origveg(:,vid) * vegcorr(vid,jv)
      ENDDO
      !
      DO jv = 1, nnobio
        frac_nobio(:,jv) = frac_nobio(:,jv) + frac_origveg(:,vid) * nobiocorr(vid,jv)
      ENDDO
      !
    ENDDO
    !
    !
    WRITE(numout,*)
    WRITE(numout,*) 'Interpolation Done'
    !
    !   Clean up the point of the map
    !
    DO ib = 1, nbpt
       !
       !  Let us see if all points found something in the 5km map !
       !
       IF ( n_found(ib) .EQ. 0 ) THEN
          !
          ! Now we need to handle some exceptions
          !
          IF ( lalo(ib,1) .LT. -56.0) THEN
             ! Antartica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 70.0) THEN
             ! Artica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 55.0 .AND. lalo(ib,2) .GT. -65.0 .AND. lalo(ib,2) .LT. -20.0) THEN
             ! Greenland
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE
             !
             WRITE(numout,*) 'PROBLEM, no point in the 5km map found for this grid box'
             WRITE(numout,*) 'Longitude range : ', lon_low(ib), lon_up(ib)
             WRITE(numout,*) 'Latitude range : ', lat_low(ib), lat_up(ib)
             !
             WRITE(numout,*) 'Looking for nearest point on the 5 km map'
             CALL slowproc_nearest (iml, lon_ful, lat_ful, &
                                    lalo(ib,2), lalo(ib,1), inear)
             WRITE(numout,*) 'Coordinates of the nearest point:', &
                              lon_ful(inear),lat_ful(inear)
             !
             DO jv = 1, nvm
               veget(ib,jv) = vegcorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
             DO jv = 1, nnobio
               frac_nobio(ib,jv) = nobiocorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
          ENDIF
          !
       ENDIF
       !
       !
       !  Limit the smalest vegetation fraction to 0.5%
       !
       DO vid = 1, nvm
          IF ( veget(ib,vid) .LT. min_vegfrac ) THEN
             veget(ib,vid) = zero
          ENDIF
       ENDDO
       !
       sumf = SUM(frac_nobio(ib,:))+SUM(veget(ib,:))
       frac_nobio(ib,:) = frac_nobio(ib,:)/sumf
       veget(ib,:) = veget(ib,:)/sumf
       !
       !       
    ENDDO
    !
    DEALLOCATE(lon_up)
    DEALLOCATE(lon_low)
    DEALLOCATE(lat_up)
    DEALLOCATE(lat_low)
    DEALLOCATE(lat_ful)
    DEALLOCATE(lon_ful)
    DEALLOCATE(vegmap)
    !
    RETURN
    !
 END SUBROUTINE slowproc_interpol_OLD


! NEED TO GET RID OF SLOWPROC_INTERPOL_NEW because it is never called!
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interpol_NEW
!!
!>\BRIEF         Interpolate the IGBP vegetation map to the grid of the model
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_interpol_NEW(nbpt, lalo, neighbours, resolution, contfrac, veget, frac_nobio )
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)         !! Vector of latitude and longitudes (beware of the order!)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,8)    !! Vector of neighbours for each grid point 
                                                                 !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)   !! The size in km of each grid-box in X and Y
    REAL(r_std),DIMENSION (nbpt), INTENT (in) :: contfrac        !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  veget(nbpt,nvm)              !! Vegetation fractions
    REAL(r_std), INTENT(out)    ::  frac_nobio(nbpt,nnobio)      !! Fraction of the mesh which is covered by ice, lakes, ...
    !
    LOGICAL ::           ok_interpol                             !! optionnal return of aggregate_vec
    !
    !  0.3 LOCAL
    !
    INTEGER(i_std), PARAMETER  :: nolson = 94                    !! Number of Olson classes
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, vid
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_ful, lon_ful, vegmap
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: sub_area
    INTEGER(i_std),ALLOCATABLE, DIMENSION(:,:) :: sub_index
    REAL(r_std), DIMENSION(nbpt,nolson) :: n_origveg
    REAL(r_std), DIMENSION(nbpt) :: n_found
    REAL(r_std), DIMENSION(nbpt,nolson) :: frac_origveg
    REAL(r_std) :: vegcorr(nolson,nvm)
    REAL(r_std) :: nobiocorr(nolson,nnobio)
    CHARACTER(LEN=40) :: callsign
    REAL(r_std) :: sumf, resol_lon, resol_lat
    INTEGER(i_std) :: idi, jv, inear, nbvmax, nix, njx
    !
    INTEGER                  :: ALLOC_ERR

!_ ================================================================================================================================
    !
    n_origveg(:,:) = zero
    n_found(:) = zero
    !
    CALL get_vegcorr (nolson,vegcorr,nobiocorr)
    !
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = NOT(IMPOSE_VEG) and NOT(LAND_USE)
    !Config Def   = carteveg5km.nc
    !Config Help  = The name of the file to be opened to read the vegetation
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from the IGBP one. We assume that we have
    !Config         a classification in 87 types. This is Olson modified by Viovy.
    !Config Units = [FILE]
    !
    filename = 'carteveg5km.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    !
    if (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    !
    !
    ALLOC_ERR=-1
    ALLOCATE(lat_ful(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_ful : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_ful(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_ful : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(vegmap(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of vegmap : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    WRITE(numout,*) 'Reading the OLSON type vegetation file'
    !
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, jml, lml, tml, 1, 1, lon_ful)
       CALL flinget(fid, 'latitude', iml, jml, lml, tml, 1, 1, lat_ful)
       CALL flinget(fid, 'vegetation_map', iml, jml, lml, tml, 1, 1, vegmap)
       !
       CALL flinclo(fid)
    ENDIF
    
    CALL bcast(lon_ful)
    CALL bcast(lat_ful)
    CALL bcast(vegmap)
    
    !
    IF (MAXVAL(vegmap) .LT. nolson) THEN
       WRITE(numout,*) 'WARNING -- WARNING'
       WRITE(numout,*) 'The vegetation map has to few vegetation types.'
       WRITE(numout,*) 'If you are lucky it will work but please check'
    ELSE IF ( MAXVAL(vegmap) .GT. nolson) THEN
       WRITE(numout,*) 'More vegetation types in file than the code can'
       WRITE(numout,*) 'deal with.: ',  MAXVAL(vegmap),  nolson
       STOP 'slowproc_interpol'
    ENDIF
    !
    ! Some assumptions on the vegetation file. This information should be
    ! be computed or read from the file. 
    ! It is the reolution in meters of the grid of the vegetation file.
    !
    resol_lon = 5000.
    resol_lat = 5000.
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution_g(:,1))*2/resol_lon)+2
       njx=INT(MAXVAL(resolution_g(:,2))*2/resol_lon)+2
       nbvmax = nix*njx
    ENDIF
    CALL bcast(nbvmax)
    !
    callsign="Vegetation map"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
    WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_index(:,:)=0
       ALLOC_ERR=-1
       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_area(:,:)=zero
       !
    WRITE(numout,*) 'Carteveg range LON:', MINVAL(lon_ful), MAXVAL(lon_ful)
    WRITE(numout,*) 'Carteveg range LAT:', MINVAL(lat_ful), MAXVAL(lat_ful)
    !
       CALL aggregate_p (nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, lon_ful, lat_ful, resol_lon, resol_lat, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          !
          nbvmax = nbvmax * 2
       ELSE
          !
          DO ib = 1, nbpt
             DO idi=1, nbvmax
                ! Leave the do loop if all sub areas are treated, sub_area <= 0
                IF ( sub_area(ib,idi) <= zero ) EXIT 
                ip = sub_index(ib,idi)
                n_origveg(ib,NINT(vegmap(ip))) = n_origveg(ib,NINT(vegmap(ip))) + sub_area(ib,idi)
                n_found(ib) =  n_found(ib) + sub_area(ib,idi)
             ENDDO
          ENDDO
          !
       ENDIF
    ENDDO
    !
    ! Now we know how many points of which Olson type from the fine grid fall
    ! into each box of the (coarse) model grid: n_origveg(nbpt,nolson)
    !
    !
    ! determine fraction of Olson vegetation type in each box of the coarse grid
    !
    DO vid = 1, nolson
       WHERE ( n_found(:) .GT. 0 ) 
          frac_origveg(:,vid) = n_origveg(:,vid) / n_found(:)
       ELSEWHERE
          frac_origveg(:,vid) = zero
       ENDWHERE
    ENDDO
    !
    ! now finally calculate coarse vegetation map
    ! Find which model vegetation corresponds to each Olson type 
    !
    veget(:,:) = zero
    frac_nobio(:,:) = zero
    !
    DO vid = 1, nolson
       !
       DO jv = 1, nvm
          veget(:,jv) = veget(:,jv) + frac_origveg(:,vid) * vegcorr(vid,jv)
       ENDDO
       !
       DO jv = 1, nnobio
          frac_nobio(:,jv) = frac_nobio(:,jv) + frac_origveg(:,vid) * nobiocorr(vid,jv)
       ENDDO
       !
    ENDDO
    !
    WRITE(numout,*) 'slowproc_interpol : Interpolation Done'
    !
    !   Clean up the point of the map
    !
    DO ib = 1, nbpt
       !
       !  Let us see if all points found something in the 5km map !
       !
       IF ( n_found(ib) .EQ. 0 ) THEN
          !
          ! Now we need to handle some exceptions
          !
          IF ( lalo(ib,1) .LT. -56.0) THEN
             ! Antartica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 70.0) THEN
             ! Artica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 55.0 .AND. lalo(ib,2) .GT. -65.0 .AND. lalo(ib,2) .LT. -20.0) THEN
             ! Greenland
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE
             !
             WRITE(numout,*) 'PROBLEM, no point in the 5km map found for this grid box',ib
             WRITE(numout,*) 'Longitude range : ', lalo(ib,2)
             WRITE(numout,*) 'Latitude range : ', lalo(ib,1)
             !
             WRITE(numout,*) 'Looking for nearest point on the 5 km map'
             CALL slowproc_nearest (iml, lon_ful, lat_ful, &
                  lalo(ib,2), lalo(ib,1), inear)
             WRITE(numout,*) 'Coordinates of the nearest point:', &
                  lon_ful(inear),lat_ful(inear)
             !
             DO jv = 1, nvm
                veget(ib,jv) = vegcorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
             DO jv = 1, nnobio
                frac_nobio(ib,jv) = nobiocorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
          ENDIF
          !
       ENDIF
       !
       !
       !  Limit the smalest vegetation fraction to 0.5%
       !
       DO vid = 1, nvm
          IF ( veget(ib,vid) .LT. min_vegfrac ) THEN
             veget(ib,vid) = zero
          ENDIF
       ENDDO
       !
       sumf = SUM(frac_nobio(ib,:))+SUM(veget(ib,:))
       frac_nobio(ib,:) = frac_nobio(ib,:)/sumf
       veget(ib,:) = veget(ib,:)/sumf
       !
       !       
    ENDDO
    !
    DEALLOCATE(vegmap)
    DEALLOCATE(lat_ful, lon_ful)
    DEALLOCATE(sub_index)
    DEALLOCATE(sub_area)

    !
    RETURN
    !
  END SUBROUTINE slowproc_interpol_NEW

! NEED TO GET RID OF SLOWPROC_INTERPOL_OLD_g because it correspond to version 1.6 of the code
! for CMIP3 simulations (too old)
!
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interpol_OLD_g
!!
!>\BRIEF         Interpolate the IGBP vegetation map to the grid of the model
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  !MM TAG 1.6 model !
  SUBROUTINE slowproc_interpol_OLD_g(nbpt, lalo, neighbours, resolution, veget, frac_nobio )
  !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,8)    !! Vector of neighbours for each grid point 
                                                                 !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  veget(nbpt,nvm)         !! Vegetation fractions
    REAL(r_std), INTENT(out)    ::  frac_nobio(nbpt,nnobio) !! Fraction of the mesh which is covered by ice, lakes, ...
    !
    !  0.3 LOCAL
    !
    INTEGER(i_std), PARAMETER                       :: nolson = 94      !! Number of Olson classes
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, jp, vid
    REAL(r_std) :: lev(1), date, dt, coslat
    INTEGER(i_std) :: itau(1)
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_ful, lon_ful, vegmap
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lon_up, lon_low, lat_up, lat_low
    INTEGER, DIMENSION(nbpt,nolson) :: n_origveg
    INTEGER, DIMENSION(nbpt) :: n_found
    REAL(r_std), DIMENSION(nbpt,nolson) :: frac_origveg
    REAL(r_std) :: vegcorr(nolson,nvm)
    REAL(r_std) :: nobiocorr(nolson,nnobio)
    CHARACTER(LEN=80) :: meter
    REAL(r_std) :: prog, sumf
    LOGICAL :: found
    INTEGER :: idi, ilast, ii, jv, inear, iprog
    REAL(r_std) :: domaine_lon_min, domaine_lon_max, domaine_lat_min, domaine_lat_max

!_ ================================================================================================================================
    !
    !
    CALL get_vegcorr (nolson,vegcorr,nobiocorr)
    !
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = NOT(IMPOSE_VEG)
    !Config Def   = carteveg5km.nc
    !Config Help  = The name of the file to be opened to read the vegetation
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from the IGBP one. We assume that we have
    !Config         a classification in 87 types. This is Olson modified by Viovy.
    !Config Units = [FILE]
    !
    filename = 'carteveg5km.nc'
    CALL getin('VEGETATION_FILE',filename)
    !
    CALL flininfo(filename, iml, jml, lml, tml, fid)
    !
    !
    ALLOCATE(lat_ful(iml))
    ALLOCATE(lon_ful(iml))
    ALLOCATE(vegmap(iml))
    !
    WRITE(numout,*) 'Reading the vegetation file'
    !
    CALL flinget(fid, 'longitude', iml, jml, lml, tml, 1, 1, lon_ful)
    CALL flinget(fid, 'latitude', iml, jml, lml, tml, 1, 1, lat_ful)
    CALL flinget(fid, 'vegetation_map', iml, jml, lml, tml, 1, 1, vegmap)
    !
    CALL flinclo(fid)
    
    !
    IF (MAXVAL(vegmap) .LT. nolson) THEN
      WRITE(*,*) 'WARNING -- WARNING'
      WRITE(*,*) 'The vegetation map has to few vegetation types.'
      WRITE(*,*) 'If you are lucky it will work but please check'
    ELSE IF ( MAXVAL(vegmap) .GT. nolson) THEN
      WRITE(*,*) 'More vegetation types in file than the code can'
      WRITE(*,*) 'deal with.: ',  MAXVAL(vegmap),  nolson
      STOP 'slowproc_interpol'
    ENDIF
    !
    ALLOCATE(lon_up(nbpt)) 
    ALLOCATE(lon_low(nbpt))
    ALLOCATE(lat_up(nbpt))
    ALLOCATE(lat_low(nbpt))
    !
    DO ib =1, nbpt
       !
       !  We find the 4 limits of the grid-box. As we transform the resolution of the model
       !  into longitudes and latitudes we do not have the problem of periodicity.
       !  coslat is a help variable here !
       !
       coslat = MAX(COS(lalo(ib,1) * pi/180. ), mincos )*pi/180. * R_Earth
       !
       lon_up(ib) = lalo(ib,2) + resolution(ib,1)/(2.0*coslat) 
       lon_low(ib) = lalo(ib,2) - resolution(ib,1)/(2.0*coslat) 
       !
       coslat = pi/180. * R_Earth
       !
       lat_up(ib) = lalo(ib,1) + resolution(ib,2)/(2.0*coslat) 
       lat_low(ib) = lalo(ib,1) - resolution(ib,2)/(2.0*coslat) 
       !
       !
       veget(ib,:) = zero
       frac_nobio (ib,:) = zero
       !
    ENDDO
    !
    !  Get the limits of the integration domaine so that we can speed up the calculations
    !
    domaine_lon_min = MINVAL(lon_low)
    domaine_lon_max = MAXVAL(lon_up)
    domaine_lat_min = MINVAL(lat_low)
    domaine_lat_max = MAXVAL(lat_up)
    !
!!$    WRITE(*,*) 'DOMAINE lon :', domaine_lon_min, domaine_lon_max
!!$    WRITE(*,*) 'DOMAINE lat :', domaine_lat_min, domaine_lat_max
    !
    ! Ensure that the fine grid covers the whole domain
    WHERE ( lon_ful(:) .LT. domaine_lon_min )
      lon_ful(:) = lon_ful(:) + 360.
    ENDWHERE
    !
    WHERE ( lon_ful(:) .GT. domaine_lon_max )
      lon_ful(:) = lon_ful(:) - 360.
    ENDWHERE
    !
    WRITE(numout,*) 'Interpolating the vegetation map :'
    WRITE(numout,'(2a40)')'0%--------------------------------------', &
                   & '------------------------------------100%'
    !
    ilast = 1
    n_origveg(:,:) = 0
    !
    DO ip=1,iml
      !
      !   Give a progress meter
      !
      ! prog = ip/float(iml)*79.
      !  IF ( ABS(prog - NINT(prog)) .LT. 1/float(iml)*79. ) THEN
      !   meter(NINT(prog)+1:NINT(prog)+1) = 'x'
      !   WRITE(numout, advance="no", FMT='(a)') ACHAR(13)
      !   WRITE(numout, advance="no", FMT='(a80)') meter
      ! ENDIF
      iprog = NINT(float(ip)/float(iml)*79.) - NINT(float(ip-1)/float(iml)*79.)
      IF ( iprog .NE. 0 ) THEN
        WRITE(numout,'(a1,$)') 'x'
      ENDIF
      !
      !  Only start looking for its place in the smaler grid if we are within the domaine
      !  That should speed up things !
      !
      IF ( ( lon_ful(ip) .GE. domaine_lon_min ) .AND. &
           ( lon_ful(ip) .LE. domaine_lon_max ) .AND. &
           ( lat_ful(ip) .GE. domaine_lat_min ) .AND. &
           ( lat_ful(ip) .LE. domaine_lat_max )        ) THEN
        !
        ! look for point on GCM grid which this point on fine grid belongs to.
        ! First look at the point on the model grid where we arrived just before. There is 
        ! a good chace that neighbouring points on the fine grid fall into the same model
        ! grid box.
        !
        !
        ! THERE IS A BUG HERE !!! IF THE GCM GRID SITS ON THE DATE LINE WE WILL HAVE FOR INSTANCE
        ! LON_LOW = -182 AND LON_UP = -178. THUS WE WILL ONLY PICK UP HALF THE POINTS NEEDED.
        !
        IF ( ( lon_ful(ip) .GT. lon_low(ilast) ) .AND. &
             ( lon_ful(ip) .LT. lon_up(ilast) ) .AND. &
             ( lat_ful(ip) .GT. lat_low(ilast) ) .AND. &
             ( lat_ful(ip) .LT. lat_up(ilast) )         ) THEN
          !
          ! We were lucky
          !
          n_origveg(ilast,NINT(vegmap(ip))) = n_origveg(ilast,NINT(vegmap(ip))) + 1
          !
        ELSE
          !
          ! Otherwise, look everywhere.
          ! Begin close to last grid point.
          !
          found = .FALSE. 
          idi = 1
          !
          DO WHILE ( (idi .LT. nbpt) .AND. ( .NOT. found ) )
            !
            ! forward and backward
            !
            DO ii = 1,2
              !
              IF ( ii .EQ. 1 ) THEN
                ib = ilast - idi
              ELSE
                ib = ilast + idi
              ENDIF
              !
              IF ( ( ib .GE. 1 ) .AND. ( ib .LE. nbpt ) ) THEN
                IF ( ( lon_ful(ip) .GT. lon_low(ib) ) .AND. &
                     ( lon_ful(ip) .LT. lon_up(ib) ) .AND. &
                     ( lat_ful(ip) .GT. lat_low(ib) ) .AND. &
                     ( lat_ful(ip) .LT. lat_up(ib) )         ) THEN
                  !
                  n_origveg(ib,NINT(vegmap(ip))) = n_origveg(ib,NINT(vegmap(ip))) + 1
                  ilast = ib
                  found = .TRUE.
                  !
                ENDIF
              ENDIF
              !
            ENDDO
            !
            idi = idi + 1
            !
          ENDDO
          !
        ENDIF ! lucky/not lucky
        !
      ENDIF     ! in the domain
    ENDDO

    !
    ! Now we know how many points of which Olson type from the fine grid fall
    ! into each box of the (coarse) model grid: n_origveg(nbpt,nolson)
    !

    !
    ! determine number of points of the fine grid which fall into each box of the 
    ! coarse grid
    !
    DO ib = 1, nbpt
      n_found(ib) = SUM( n_origveg(ib,:) )
    ENDDO

    !
    ! determine fraction of Olson vegetation type in each box of the coarse grid
    !
    DO vid = 1, nolson
      WHERE ( n_found(:) .GT. 0 ) 
        frac_origveg(:,vid) =  REAL(n_origveg(:,vid),r_std) /  REAL(n_found(:),r_std)
      ELSEWHERE
         frac_origveg(:,vid) = zero
      ENDWHERE
    ENDDO

    !
    ! now finally calculate coarse vegetation map
    ! Find which model vegetation corresponds to each Olson type 
    !
    DO vid = 1, nolson
      !
      DO jv = 1, nvm
        veget(:,jv) = veget(:,jv) + frac_origveg(:,vid) * vegcorr(vid,jv)
      ENDDO
      !
      DO jv = 1, nnobio
        frac_nobio(:,jv) = frac_nobio(:,jv) + frac_origveg(:,vid) * nobiocorr(vid,jv)
      ENDDO
      !
    ENDDO
    !
    !
    WRITE(numout,*)
    WRITE(numout,*) 'Interpolation Done'
    !
    !   Clean up the point of the map
    !
    DO ib = 1, nbpt
       !
       !  Let us see if all points found something in the 5km map !
       !
       IF ( n_found(ib) .EQ. 0 ) THEN
          !
          ! Now we need to handle some exceptions
          !
          IF ( lalo(ib,1) .LT. -56.0) THEN
             ! Antartica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 70.0) THEN
             ! Artica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 55.0 .AND. lalo(ib,2) .GT. -65.0 .AND. lalo(ib,2) .LT. -20.0) THEN
             ! Greenland
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE
             !
             WRITE(numout,*) 'PROBLEM, no point in the 5km map found for this grid box'
             WRITE(numout,*) 'Longitude range : ', lon_low(ib), lon_up(ib)
             WRITE(numout,*) 'Latitude range : ', lat_low(ib), lat_up(ib)
             !
             WRITE(numout,*) 'Looking for nearest point on the 5 km map'
             CALL slowproc_nearest (iml, lon_ful, lat_ful, &
                                    lalo(ib,2), lalo(ib,1), inear)
             WRITE(numout,*) 'Coordinates of the nearest point:', &
                              lon_ful(inear),lat_ful(inear)
             !
             DO jv = 1, nvm
               veget(ib,jv) = vegcorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
             DO jv = 1, nnobio
               frac_nobio(ib,jv) = nobiocorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
          ENDIF
          !
       ENDIF
       !
       !
       !  Limit the smalest vegetation fraction to 0.5%
       !
       DO vid = 1, nvm
          IF ( veget(ib,vid) .LT. min_vegfrac ) THEN
             veget(ib,vid) = zero
          ENDIF
       ENDDO
       !
       sumf = SUM(frac_nobio(ib,:))+SUM(veget(ib,:))
       frac_nobio(ib,:) = frac_nobio(ib,:)/sumf
       veget(ib,:) = veget(ib,:)/sumf
       !
       !       
    ENDDO
    !
    DEALLOCATE(lon_up)
    DEALLOCATE(lon_low)
    DEALLOCATE(lat_up)
    DEALLOCATE(lat_low)
    DEALLOCATE(lat_ful)
    DEALLOCATE(lon_ful)
    DEALLOCATE(vegmap)
    !
    RETURN
    !
  END SUBROUTINE slowproc_interpol_OLD_g


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interpol_NEW_g
!!
!>\BRIEF         Interpolate the IGBP vegetation map to the grid of the model
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::veget, ::frac_nobio
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_interpol_NEW_g(nbpt, lalo, neighbours, resolution, contfrac, veget, frac_nobio )
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)           :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)          !! Vector of latitude and longitudes 
                                                                  !! (beware of the order : 1=latitude ; 2=longitude)
    INTEGER(i_std), INTENT(in)           :: neighbours(nbpt,8)    !! Vector of neighbours for each grid point 
                                                                  !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    REAL(r_std),DIMENSION (nbpt), INTENT (in) :: contfrac         !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  veget(nbpt,nvm)               !! Vegetation fractions
    REAL(r_std), INTENT(out)    ::  frac_nobio(nbpt,nnobio)       !! Fraction of the mesh which is covered by ice, lakes, ...
    !
    LOGICAL ::           ok_interpol                              !! optionnal return of aggregate_vec
    !
    !  0.3 LOCAL
    !
    INTEGER(i_std), PARAMETER                       :: nolson = 94      !! Number of Olson classes
    !
    !
    CHARACTER(LEN=80) :: filename                                       !!vegetation map filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, vid              
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lat_ful, lon_ful, vegmap  !! for 5km vegetation map 
                                                                        !! latitude vector, longitude vector, and 
                                                                        !! value of Olson's classes for each location
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: sub_area                !! the area of the fine grid in the model grid ??? 
                                                                        !! cf src_global/interpol_help.f90, line 377, called "areaoverlap"
    INTEGER(i_std),ALLOCATABLE, DIMENSION(:,:) :: sub_index             !! the indexes from the grid boxes from the data that go 
                                                                        !! into the model's boxes 
                                                                        !! cf src_global/interpol_help.f90,line 300, called "ip"
    REAL(r_std), DIMENSION(nbpt,nolson) :: n_origveg                    !! number of points of each Olson type from the fine grid 
                                                                        !! in each box of the (coarse) model grid 
    REAL(r_std), DIMENSION(nbpt) :: n_found                             !! total number of different Olson types found in each 
                                                                        !! box of the (coarse) model grid
    REAL(r_std), DIMENSION(nbpt,nolson) :: frac_origveg                 !! fraction of each Olson type in each box of the (coarse) model grid
    REAL(r_std) :: vegcorr(nolson,nvm)                                  !! correspondance table between Olson and the following SECHIBA Classes.
                                                                        !!   vegcorr(i,:)+nobiocorr(i,:) = 1.  for all i 
                                                                        !! see each class in src_parameters/constantes_veg.f90

    REAL(r_std) :: nobiocorr(nolson,nnobio)                             !! non-biospheric surface typesi
    CHARACTER(LEN=40) :: callsign                                       !! Allows to specify which variable is beeing treated
    REAL(r_std) :: sumf, resol_lon, resol_lat                           !! sumf = sum veget + sum nobio
                                                                        !! resol_lon, resol_lat  reolution in meters of the grid of the vegetation file
    INTEGER(i_std) :: idi, jv, inear, nbvmax                            !! idi : counter for nbvmax, see below   
                                                                        !! jv : counter for nvm, number of PFT
                                                                        !! inear : location of the point of vegmap, which is the closest from the modelled point
                                                                        !! nbvmax : number of maximum vegetation map points in the GCM grid 
    INTEGER(i_std) :: nix, njx
    !
    INTEGER                  :: ALLOC_ERR                               !! location of the eventual missing value in vegmap

!_ ================================================================================================================================
    !
    n_origveg(:,:) = zero
    n_found(:) = zero
    !
    CALL get_vegcorr (nolson,vegcorr,nobiocorr)
    !
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = NOT(IMPOSE_VEG) and NOT(LAND_USE)
    !Config Def   = carteveg5km.nc
    !Config Help  = The name of the file to be opened to read the vegetation
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from the IGBP one. We assume that we have
    !Config         a classification in 87 types. This is Olson modified by Viovy.
    !Config Units = [FILE]
    !
    filename = 'carteveg5km.nc'
    CALL getin('VEGETATION_FILE',filename)  ! GETIN_P !!
    !
    CALL flininfo(filename, iml, jml, lml, tml, fid)   
    !
    ! see IOIPSL/src/flincom.f90, line 665
    ! fid      : File ID
    !- iml      | These 4 variables give the size of the variables
    !- jml      | to be read. It will be verified that the variables
    !- lml      | fits in there.
    !- tml     
    ! iml, jml : horizontal size of the grid, lml = vertical size
    ! tml : size of time axis

    ! TL : pourquoi 2 variables pour la taille horizontale ? cf
    ! IOIPSL/src/flincom.f90 , line 160 

    ALLOC_ERR=-1
    ALLOCATE(lat_ful(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_ful : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_ful(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_ful : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(vegmap(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of vegmap : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    WRITE(numout,*) 'Reading the OLSON type vegetation file'
    !
    CALL flinget(fid, 'longitude', iml, jml, lml, tml, 1, 1, lon_ful)
    CALL flinget(fid, 'latitude', iml, jml, lml, tml, 1, 1, lat_ful)
    CALL flinget(fid, 'vegetation_map', iml, jml, lml, tml, 1, 1, vegmap)
    !
    WRITE(numout,*) 'File name : ', filename
    WRITE(numout,*) 'Min and max vegetation numbers : ', MINVAL(vegmap), MAXVAL(vegmap)
    !
    CALL flinclo(fid)
    !
    IF (MAXVAL(vegmap) .LT. nolson) THEN
       WRITE(numout,*) 'WARNING -- WARNING'
       WRITE(numout,*) 'The vegetation map has too few vegetation types.'
       WRITE(numout,*) 'If you are lucky it will work but please check'
    ELSE IF ( MAXVAL(vegmap) .GT. nolson) THEN
       WRITE(numout,*) 'More vegetation types in file than the code can'
       WRITE(numout,*) 'deal with.: ',  MAXVAL(vegmap),  nolson
       STOP 'slowproc_interpol'
    ENDIF
    !
    ! Some assumptions on the vegetation file. This information should be
    ! be computed or read from the file. 
    ! It is the reolution in meters of the grid of the vegetation file.
    !
    
    !TL : CODE EN DUR ????? 
    resol_lon = 5000.
    resol_lat = 5000.
    !
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some margin is taken.
    !
    nix=INT(MAXVAL(resolution_g(:,1)*2)/resol_lon)+1
    njx=INT(MAXVAL(resolution_g(:,2)*2)/resol_lon)+1
    nbvmax = nix*njx
    !
    ! No need to broadcast as this routine is only called on root_proc
    !
    callsign="Vegetation map"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
    WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_index(:,:)=0
       ALLOC_ERR=-1
       ALLOCATE(sub_area(nbpt, nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate (nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, lon_ful, lat_ful, resol_lon, resol_lat, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       ! Defined as aggregate_2d or aggregate_vec in src_global/interpol_help.f90, depending
       ! on the dimensions (2D region or vector)i. 
       ! This routing will get for each point of the coarse grid the
       ! indexes of the finer grid and the area of overlap.
       ! This routine is designed for a fine grid which is regular in lat/lon.

       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          !
          nbvmax = nbvmax * 2
       ELSE
          !
          DO ib = 1, nbpt
             DO idi=1, nbvmax
                ! Leave the do loop if all sub areas are treated, sub_area <= 0
                IF ( sub_area(ib,idi) <= zero ) EXIT

                ip = sub_index(ib,idi)
                n_origveg(ib,NINT(vegmap(ip))) = n_origveg(ib,NINT(vegmap(ip))) + sub_area(ib,idi)
                n_found(ib) =  n_found(ib) + sub_area(ib,idi)
             ENDDO
          ENDDO
          !
       ENDIF
    ENDDO
    !
    ! Now we know how many points of which Olson type from the fine grid fall
    ! into each box of the (coarse) model grid: n_origveg(nbpt,nolson)
    !
    !
    ! determine fraction of Olson vegetation type in each box of the coarse grid
    !
    DO vid = 1, nolson
       WHERE ( n_found(:) .GT. 0 ) 
          frac_origveg(:,vid) = n_origveg(:,vid) / n_found(:)
       ELSEWHERE
          frac_origveg(:,vid) = zero
       ENDWHERE
    ENDDO
    !
    ! now finally calculate coarse vegetation map
    ! Find which model vegetation corresponds to each Olson type 
    !
    veget(:,:) = zero
    frac_nobio(:,:) = zero
    !
    DO vid = 1, nolson
       !
       DO jv = 1, nvm
          veget(:,jv) = veget(:,jv) + frac_origveg(:,vid) * vegcorr(vid,jv)
       ENDDO
       !
       DO jv = 1, nnobio
          frac_nobio(:,jv) = frac_nobio(:,jv) + frac_origveg(:,vid) * nobiocorr(vid,jv)
       ENDDO
       !
    ENDDO
    !
    WRITE(numout,*) 'slowproc_interpol : Interpolation Done'
    !
    !   Clean up the point of the map
    !
    DO ib = 1, nbpt
       !
       !  Let us see if all points found something in the 5km map !
       !
       IF ( n_found(ib) .EQ. 0 ) THEN
          !
          ! Now we need to handle some exceptions
          !
          IF ( lalo(ib,1) .LT. -56.0) THEN
             ! Antartica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 70.0) THEN
             ! Artica
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE IF ( lalo(ib,1) .GT. 55.0 .AND. lalo(ib,2) .GT. -65.0 .AND. lalo(ib,2) .LT. -20.0) THEN
             ! Greenland
             frac_nobio(ib,:) = zero
             frac_nobio(ib,iice) = un
             veget(ib,:) = zero
             !
          ELSE
             !
             WRITE(numout,*) 'PROBLEM, no point in the 5km map found for this grid box',ib
             WRITE(numout,*) 'Longitude range : ', lalo(ib,2)
             WRITE(numout,*) 'Latitude range : ', lalo(ib,1)
             !
             WRITE(numout,*) 'Looking for nearest point on the 5 km map'
             CALL slowproc_nearest (iml, lon_ful, lat_ful, &
                  lalo(ib,2), lalo(ib,1), inear)
             WRITE(numout,*) 'Coordinates of the nearest point:', &
                  lon_ful(inear),lat_ful(inear)
             !
             DO jv = 1, nvm
                veget(ib,jv) = vegcorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
             DO jv = 1, nnobio
                frac_nobio(ib,jv) = nobiocorr(NINT(vegmap(inear)),jv)
             ENDDO
             !
          ENDIF
          !
       ENDIF
       !
       !
       !  Limit the smallest vegetation fraction to 0.5%
       !
       DO vid = 1, nvm
          IF ( veget(ib,vid) .LT. min_vegfrac ) THEN  ! min_vegfrac=0.001 in constantes_veg.f90
             veget(ib,vid) = zero
          ENDIF
       ENDDO
       !
       sumf = SUM(frac_nobio(ib,:))+SUM(veget(ib,:))
       frac_nobio(ib,:) = frac_nobio(ib,:)/sumf
       veget(ib,:) = veget(ib,:)/sumf
       !
       !       
    ENDDO
    !
    DEALLOCATE(vegmap)
    DEALLOCATE(lat_ful, lon_ful)
    DEALLOCATE(sub_index)
    DEALLOCATE(sub_area)

    !
    RETURN
    !
  END SUBROUTINE slowproc_interpol_NEW_g


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_nearest
!!
!>\BRIEF         looks for nearest grid point on the fine map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::inear
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_nearest(iml, lon5, lat5, lonmod, latmod, inear)

    !! INTERFACE DESCRIPTION
    
    !! 0.1 input variables

    INTEGER(i_std), INTENT(in)                   :: iml             !! size of the vector
    REAL(r_std), DIMENSION(iml), INTENT(in)      :: lon5, lat5      !! longitude and latitude vector, for the 5km vegmap
    REAL(r_std), INTENT(in)                      :: lonmod, latmod  !! longitude  and latitude modelled

    !! 0.2 output variables
    
    INTEGER(i_std), INTENT(out)                  :: inear           !! location of the grid point from the 5km vegmap grid
                                                                    !! closest from the modelled grid point

    !! 0.4 Local variables

    REAL(r_std)                                  :: pa, p
    REAL(r_std)                                  :: coscolat, sincolat
    REAL(r_std)                                  :: cospa, sinpa
    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: cosang
    INTEGER(i_std)                               :: i
    INTEGER(i_std), DIMENSION(1)                 :: ineartab
    INTEGER                                      :: ALLOC_ERR

!_ ================================================================================================================================

    ALLOC_ERR=-1
    ALLOCATE(cosang(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of cosang : ",ALLOC_ERR
      STOP 
    ENDIF

    pa = pi/2.0 - latmod*pi/180.0 ! dist. between north pole and the point a 
                                                      !! COLATITUDE, in radian
    cospa = COS(pa)
    sinpa = SIN(pa)

    DO i = 1, iml

       sincolat = SIN( pi/2.0 - lat5(i)*pi/180.0 ) !! sinus of the colatitude
       coscolat = COS( pi/2.0 - lat5(i)*pi/180.0 ) !! cosinus of the colatitude

       p = (lonmod-lon5(i))*pi/180.0 !! angle between a & b (between their meridian)in radians

       !! dist(i) = ACOS( cospa*coscolat + sinpa*sincolat*COS(p))
       cosang(i) = cospa*coscolat + sinpa*sincolat*COS(p) !! TL : cosang is maximum when angle is at minimal value  
!! orthodromic distance between 2 points : cosang = cosinus (arc(AB)/R), with
!R = Earth radius, then max(cosang) = max(cos(arc(AB)/R)), reached when arc(AB)/R is minimal, when
! arc(AB) is minimal, thus when point B (corresponding grid point from LAI MAP) is the nearest from
! modelled A point
    ENDDO

    ineartab = MAXLOC( cosang(:) )
    inear = ineartab(1)

    DEALLOCATE(cosang)
  END SUBROUTINE slowproc_nearest

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_soilt
!!
!>\BRIEF         Interpolate the Zobler soil type map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::soiltype, ::clayfraction
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_soilt(nbpt, lalo, neighbours, resolution, contfrac, soilclass, clayfraction)
    !
    !
    !   This subroutine should read the Zobler map and interpolate to the model grid. The method
    !   is to get fraction of the three main soiltypes for each grid box.
    !   The soil fraction are going to be put into the array soiltype in the following order :
    !   coarse, medium and fine.
    !
    !
    !!  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)    :: nbpt                   !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)       :: lalo(nbpt,2)           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)    :: neighbours(nbpt,8)     !! Vector of neighbours for each grid point 
                                                            !! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)       :: resolution(nbpt,2)     !! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT(in)       :: contfrac(nbpt)         !! Fraction of land in each grid box.
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)      :: soilclass(nbpt, nstm)   !! Soil type map to be created from the Zobler map
    REAL(r_std), INTENT(out)      :: clayfraction(nbpt)     !! The fraction of clay as used by STOMATE
    !
    !
    !  0.3 LOCAL
    !
    INTEGER(i_std)               :: nbvmax
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, ip, jp, fopt, ilf, nbexp
    REAL(r_std) :: lev(1), date, dt
    INTEGER(i_std) :: itau(1)
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: lat_rel, lon_rel
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: soiltext, soiltext2
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)  :: sub_area
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:)  :: sub_index
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: resol_lu
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: solt, solt2
    REAL(r_std) ::  sgn, coslat
    CHARACTER(LEN=30) :: callsign
    INTEGER(i_std) :: nix, njx
    !
    ! Number of texture classes in Zobler
    !
    INTEGER(i_std), PARAMETER :: nzobler = 7
    REAL(r_std),ALLOCATABLE   :: textfrac_table(:,:)
    !   
    LOGICAL                  :: ok_interpol  ! optionnal return of aggregate_2d
    !   
    INTEGER                  :: ALLOC_ERR

!_ ================================================================================================================================
    !
    !
    !  Needs to be a configurable variable
    !
    !
    !Config Key   = SOILCLASS_FILE
    !Config Desc  = Name of file from which soil types are read
    !Config Def   = soils_param.nc
    !Config If    = NOT(IMPOSE_VEG)
    !Config Help  = The name of the file to be opened to read the soil types. 
    !Config         The data from this file is then interpolated to the grid of
    !Config         of the model. The aim is to get fractions for sand loam and
    !Config         clay in each grid box. This information is used for soil hydrology
    !Config         and respiration.
    !Config Units = [FILE]
    !
    filename = 'soils_param.nc'
    CALL getin_p('SOILCLASS_FILE',filename)
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
    ! soils_param.nc file is 1° soit texture file.
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
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(soiltext(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of soiltext : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(soiltext2(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of soiltext2 : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(resol_lu(iml,jml,2), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of resol_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    IF (is_root_prc) CALL flinopen(filename, .FALSE., iml, jml, lml, lon_rel, lat_rel, lev, tml, itau, date, dt, fid)
    CALL bcast(lon_rel)
    CALL bcast(lat_rel)
    CALL bcast(itau)
    CALL bcast(date)
    CALL bcast(dt)
    
    !
    IF (is_root_prc) CALL flinget(fid, 'soiltext', iml, jml, lml, tml, 1, 1, soiltext)
    CALL bcast(soiltext)
    !
    !
    IF (is_root_prc) CALL flinclo(fid)
    !
    nbexp = 0
    !
    !
    ! Mask of permitted variables.
    !
    mask(:,:) = zero
    DO ip=1,iml
       DO jp=1,jml
          !
          IF (soiltext(ip,jp) .GT. min_sechiba) THEN
             mask(ip,jp) = un
          ENDIF
          !
          ! Resolution in longitude
          !
          coslat = MAX( COS( lat_rel(ip,jp) * pi/180. ), mincos )     
          IF ( ip .EQ. 1 ) THEN
             resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip,jp) ) * pi/180. * R_Earth * coslat
          ELSEIF ( ip .EQ. iml ) THEN
             resol_lu(ip,jp,1) = ABS( lon_rel(ip,jp) - lon_rel(ip-1,jp) ) * pi/180. * R_Earth * coslat
          ELSE
             resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip-1,jp) )/2. * pi/180. * R_Earth * coslat
          ENDIF
          !
          ! Resolution in latitude
          !
          IF ( jp .EQ. 1 ) THEN
             resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp) - lat_rel(ip,jp+1) ) * pi/180. * R_Earth
          ELSEIF ( jp .EQ. jml ) THEN
             resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp) ) * pi/180. * R_Earth
          ELSE
             resol_lu(ip,jp,2) =  ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp+1) )/2. * pi/180. * R_Earth
          ENDIF
          !
       ENDDO
    ENDDO
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution_g(:,1))/MAXVAL(resol_lu(:,:,1)))+2
       njx=INT(MAXVAL(resolution_g(:,2))/MAXVAL(resol_lu(:,:,2)))+2
       nbvmax = nix*njx
    ENDIF
    CALL bcast(nbvmax)
    !
    callsign = "Soil types"
    !
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       WRITE(numout,*) "Projection arrays for ",callsign," : "
    WRITE(numout,*) "nbvmax = ",nbvmax, nix, njx
       !
       ALLOC_ERR=-1
       ALLOCATE(solt(nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of solt : ",ALLOC_ERR
          STOP 
       ENDIF
       ALLOC_ERR=-1
    ALLOCATE(solt2(nbvmax), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
       WRITE(numout,*) "ERROR IN ALLOCATION of solt2 : ",ALLOC_ERR
       STOP 
    ENDIF
    ALLOC_ERR=-1
       ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_index(:,:,:)=0
       ALLOC_ERR=-1
       ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon_rel, lat_rel, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       !
       IF ( .NOT. ok_interpol ) THEN
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          DEALLOCATE(solt)
          DEALLOCATE(solt2)
          !
          nbvmax = nbvmax * 2
       ENDIF
    ENDDO
    !
    !
    SELECTCASE(soil_classif)
    CASE('none')
       ALLOCATE(textfrac_table(nscm,ntext))
       DO ib=1, nbpt
          soilclass(ib,:) = soilclass_default_fao
          clayfraction(ib) = clayfraction_default
       ENDDO
    CASE('zobler')
       !
       soilclass_default=soilclass_default_fao
       !
       WRITE(numout,*) "Using a soilclass map with Zobler classification"
       !
       ALLOCATE(textfrac_table(nzobler,ntext))
       !
       CALL get_soilcorr (nzobler, textfrac_table)
       !
       !
       DO ib =1, nbpt
          !
          ! GO through the point we have found
          !
          !
          fopt = COUNT(sub_area(ib,:) > zero)
          !
          !    Check that we found some points
          !
          soilclass(ib,:) = zero
          clayfraction(ib) = zero
          !
          IF ( fopt .EQ. 0) THEN
             nbexp = nbexp + 1
             soilclass(ib,:) = soilclass_default(:)
             clayfraction(ib) = clayfraction_default
          ELSE
             !
             DO ilf = 1,fopt
                solt(ilf) = soiltext(sub_index(ib,ilf,1),sub_index(ib,ilf,2))
             ENDDO
             !
             sgn = zero
             !
             !   Compute the average bare soil albedo parameters
             !
             DO ilf = 1,fopt
                !
                ! We have to take care of two exceptions here : type 6 = glacier and type 0 = ocean
                !
                IF ( (solt(ilf) .LE. nzobler) .AND. (solt(ilf) .GT. 0) .AND.&
                     & (solt(ilf) .NE. 6)) THEN
                   SELECTCASE(solt(ilf))
                   CASE(1)
                      soilclass(ib,1) = soilclass(ib,1) + sub_area(ib,ilf)
                   CASE(2)
                      soilclass(ib,2) = soilclass(ib,2) + sub_area(ib,ilf)
                   CASE(3)
                      soilclass(ib,2) = soilclass(ib,2) + sub_area(ib,ilf)
                   CASE(4)
                      soilclass(ib,2) = soilclass(ib,2) + sub_area(ib,ilf)
                   CASE(5)
                      soilclass(ib,3) = soilclass(ib,3) + sub_area(ib,ilf)
                   CASE(7)
                      soilclass(ib,2) = soilclass(ib,2) + sub_area(ib,ilf)
                   CASE DEFAULT
                      WRITE(numout,*) 'We should not be here, an impossible case appeared'
                      STOP 'slowproc_soilt'
                   END SELECT
                   clayfraction(ib) = clayfraction(ib) + &
                        & textfrac_table(solt(ilf),3) * sub_area(ib,ilf)
                   sgn = sgn + sub_area(ib,ilf)
                ELSE
                   IF (solt(ilf) .GT. nzobler) THEN
                      WRITE(numout,*) 'The file contains a soil color class which is incompatible with this program'
                      STOP 'slowproc_soilt'
                   ENDIF
                ENDIF
                !
             ENDDO
             !
             ! Normalize the surface
             !
             IF ( sgn .LT. min_sechiba) THEN
                nbexp = nbexp + 1
                soilclass(ib,:) = soilclass_default(:)
                clayfraction(ib) = clayfraction_default
             ELSE
                soilclass(ib,:) = soilclass(ib,:)/sgn
                clayfraction(ib) = clayfraction(ib)/sgn
             ENDIF
             !
          ENDIF
          !
       ENDDO
   
    !
    CASE("fao")
       !
       soilclass_default=soilclass_default_fao
       !
       WRITE(numout,*) "Using a soilclass map with fao classification"
       !
       ALLOCATE(textfrac_table(nscm,ntext))
       !
       CALL get_soilcorr (nscm, textfrac_table)
       !
       DO ib =1, nbpt
          !
          ! GO through the point we have found
          !
          !
          fopt = COUNT(sub_area(ib,:) > zero)
          !
          !    Check that we found some points
          !
          soilclass(ib,:) = 0.0
          clayfraction(ib) = 0.0
          !
          IF ( fopt .EQ. 0) THEN
             nbexp = nbexp + 1
             soilclass(ib,:) = soilclass_default(:)
             clayfraction(ib) = clayfraction_default
          ELSE
             !
             DO ilf = 1,fopt
                solt(ilf) = soiltext(sub_index(ib,ilf,1),sub_index(ib,ilf,2))
             ENDDO
             !
             !
             !   Compute the average bare soil albedo parameters
             !
             sgn = zero
             !
             DO ilf = 1,fopt
                !
                ! 
                !
                IF ( (solt(ilf) .LE. nscm) .AND. (solt(ilf) .GT. 0) ) THEN
                   soilclass(ib,solt(ilf)) = soilclass(ib,solt(ilf)) + sub_area(ib,ilf)
                   clayfraction(ib) = clayfraction(ib) + textfrac_table(solt(ilf),3) * sub_area(ib,ilf)
                   sgn = sgn + sub_area(ib,ilf)
                ELSE
                   IF (solt(ilf) .GT. nscm) THEN
                      WRITE(*,*) 'The file contains a soil color class which is incompatible with this program'
                      STOP 'slowproc_soilt'
                   ENDIF
                ENDIF
                !
             ENDDO
             !
             ! Normalize the surface
             !
             IF ( sgn .LT. min_sechiba) THEN
                nbexp = nbexp + 1
                soilclass(ib,:) = soilclass_default(:)
                clayfraction(ib) = clayfraction_default
             ELSE
                soilclass(ib,:) = soilclass(ib,:)/sgn
                clayfraction(ib) = clayfraction(ib)/sgn
             ENDIF
             !
          ENDIF
          !
       ENDDO
       !
       !
    CASE("usda")
       !
       soilclass_default=soilclass_default_usda
       !
       WRITE(numout,*) "Using a soilclass map with usda classification"
       !
       ALLOCATE(textfrac_table(nscm,ntext))
       !
       CALL get_soilcorr (nscm, textfrac_table)
       !
       DO ib =1, nbpt
          !
          ! GO through the point we have found
          !
          !
          fopt = COUNT(sub_area(ib,:) > zero)
          !
          !    Check that we found some points
          !
          soilclass(ib,:) = 0.0
          clayfraction(ib) = 0.0
          !
          IF ( fopt .EQ. 0) THEN
             nbexp = nbexp + 1
             soilclass(ib,:) = soilclass_default
             clayfraction(ib) = clayfraction_default
          ELSE
             !
             DO ilf = 1,fopt
                solt(ilf) = soiltext(sub_index(ib,ilf,1),sub_index(ib,ilf,2))
             ENDDO
             !
             !
             !   Compute the average bare soil albedo parameters
             !
             sgn = zero
             !
             DO ilf = 1,fopt
                !
                ! 
                !
                IF ( (solt(ilf) .LE. nscm) .AND. (solt(ilf) .GT. 0) ) THEN
                   soilclass(ib,solt(ilf)) = soilclass(ib,solt(ilf)) + sub_area(ib,ilf)
                   clayfraction(ib) = clayfraction(ib) + textfrac_table(solt(ilf),3) * sub_area(ib,ilf)
                   sgn = sgn + sub_area(ib,ilf)
                ELSE
                   IF (solt(ilf) .GT. nscm) THEN
                      WRITE(*,*) 'The file contains a soil color class which is incompatible with this program'
                      STOP 'slowproc_soilt'
                   ENDIF
                ENDIF
                !
             ENDDO
             !
             ! Normalize the surface
             !
             IF ( sgn .LT. min_sechiba) THEN
                nbexp = nbexp + 1
                soilclass(ib,:) = soilclass_default(:)
                clayfraction(ib) = clayfraction_default
             ELSE
                soilclass(ib,:) = soilclass(ib,:)/sgn
                clayfraction(ib) = clayfraction(ib)/sgn
             ENDIF
             !
          ENDIF
          !
       ENDDO
       !
    CASE DEFAULT
       WRITE(*,*) 'A non supported soil type classification has been chosen'
       STOP 'slowproc_soilt'
    ENDSELECT
    !
    WRITE(numout,*) 'Interpolation Done'
       !
    IF ( nbexp .GT. 0 ) THEN
       WRITE(numout,*) 'slowproc_soilt : The interpolation of the bare soil albedo had ', nbexp
       WRITE(numout,*) 'slowproc_soilt : points without data. This are either coastal points or'
       WRITE(numout,*) 'slowproc_soilt : ice covered land.'
       WRITE(numout,*) 'slowproc_soilt : The problem was solved by using the default soil types.'
    ENDIF
    !
    DEALLOCATE (lat_rel)
    DEALLOCATE (lon_rel)
    DEALLOCATE (mask)
    DEALLOCATE (sub_area)
    DEALLOCATE (sub_index)
    DEALLOCATE (soiltext)
    DEALLOCATE (solt)
    !
    DEALLOCATE (soiltext2)
    DEALLOCATE (solt2)
    DEALLOCATE (textfrac_table)
    !
    DEALLOCATE (resol_lu)
    !
    !
    RETURN
    !
  END SUBROUTINE slowproc_soilt
 
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_slope
!!
!>\BRIEF         Calculate mean slope coef in each  model grid box from the slope map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::reinf_slope
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_slope(nbpt, lalo, neighbours, resolution, contfrac, reinf_slope)
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)          ! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,8)    ! Vector of neighbours for each grid point 
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)    ! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT (in)             :: contfrac(nbpt)         !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  reinf_slope(nbpt)                   ! slope coef 
    !
    !  0.3 LOCAL
    !
    !
    REAL(r_std)  :: slope_noreinf                 ! Slope above which runoff is maximum
    CHARACTER(LEN=80) :: filename
    CHARACTER(LEN=30) :: callsign
    INTEGER(i_std)    :: iml, jml, lml, tml, fid, ib, ip, jp, vid
    INTEGER(i_std)    :: idi, idi_last, nbvmax
    REAL(r_std)        :: slopecoef, coslat
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:) :: mask
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:,:) :: sub_index
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: lat_rel, lon_rel, slopemap
    REAL(r_std), ALLOCATABLE, DIMENSION(:)      :: lat_lu, lon_lu
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)    :: sub_area
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: resol_lu
    INTEGER(i_std) :: nix, njx
    !
    LOGICAL ::           ok_interpol = .FALSE. ! optionnal return of aggregate_2d
    !
    INTEGER                  :: ALLOC_ERR
!_ ================================================================================================================================

    !
    !Config Key   = SLOPE_NOREINF
    !Config Desc  = See slope_noreinf above
    !Config If    = 
    !Config Def   = 0.5
    !Config Help  = The slope above which there is no reinfiltration
    !Config Units = [-]
    !
    slope_noreinf = 0.5
    !
    CALL getin_p('SLOPE_NOREINF',slope_noreinf)
    !
    !Config Key   = TOPOGRAPHY_FILE
    !Config Desc  = Name of file from which the topography map is to be read
    !Config If    = 
    !Config Def   = cartepente2d_15min.nc
    !Config Help  = The name of the file to be opened to read the orography
    !Config         map is to be given here. Usualy SECHIBA runs with a 2'
    !Config         map which is derived from the NGDC one. 
    !Config Units = [FILE]
    !
    filename = 'cartepente2d_15min.nc'
    CALL getin_p('TOPOGRAPHY_FILE',filename)
    !
    IF (is_root_prc) CALL flininfo(filename, iml, jml, lml, tml, fid)
    CALL bcast(iml)
    CALL bcast(jml)
    CALL bcast(lml)
    CALL bcast(tml)
    !
    ALLOC_ERR=-1
    ALLOCATE(lat_lu(jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lon_lu(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(slopemap(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of slopemap : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(resol_lu(iml,jml,2), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of resol_lu : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    WRITE(numout,*) 'Reading the topography file'
    !
    IF (is_root_prc) THEN
       CALL flinget(fid, 'longitude', iml, 0, 0, 0, 1, 1, lon_lu)
       CALL flinget(fid, 'latitude', jml, 0, 0, 0, 1, 1, lat_lu)
       CALL flinget(fid, 'pente', iml, jml, 0, 0, 1, 1, slopemap)
       !
       CALL flinclo(fid)
    ENDIF
    CALL bcast(lon_lu)
    CALL bcast(lat_lu)
    CALL bcast(slopemap)
    !
    ALLOC_ERR=-1
    ALLOCATE(lon_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lon_rel : ",ALLOC_ERR
      STOP 
    ENDIF
    ALLOC_ERR=-1
    ALLOCATE(lat_rel(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of lat_rel : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    DO ip=1,iml
       lat_rel(ip,:) = lat_lu(:)
    ENDDO
    DO jp=1,jml
       lon_rel(:,jp) = lon_lu(:)
    ENDDO
    !
    !
    ! Mask of permitted variables.
    !
    ALLOC_ERR=-1
    ALLOCATE(mask(iml,jml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) THEN
      WRITE(numout,*) "ERROR IN ALLOCATION of mask : ",ALLOC_ERR
      STOP 
    ENDIF
    !
    mask(:,:) = zero
    DO ip=1,iml
       DO jp=1,jml
          IF (slopemap(ip,jp) .GT. min_sechiba) THEN
             mask(ip,jp) = un
          ENDIF
          !
          ! Resolution in longitude
          !
          coslat = MAX( COS( lat_rel(ip,jp) * pi/180. ), mincos )     
          IF ( ip .EQ. 1 ) THEN
             resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip,jp) ) * pi/180. * R_Earth * coslat
          ELSEIF ( ip .EQ. iml ) THEN
             resol_lu(ip,jp,1) = ABS( lon_rel(ip,jp) - lon_rel(ip-1,jp) ) * pi/180. * R_Earth * coslat
          ELSE
             resol_lu(ip,jp,1) = ABS( lon_rel(ip+1,jp) - lon_rel(ip-1,jp) )/2. * pi/180. * R_Earth * coslat
          ENDIF
          !
          ! Resolution in latitude
          !
          IF ( jp .EQ. 1 ) THEN
             resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp) - lat_rel(ip,jp+1) ) * pi/180. * R_Earth
          ELSEIF ( jp .EQ. jml ) THEN
             resol_lu(ip,jp,2) = ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp) ) * pi/180. * R_Earth
          ELSE
             resol_lu(ip,jp,2) =  ABS( lat_rel(ip,jp-1) - lat_rel(ip,jp+1) )/2. * pi/180. * R_Earth
          ENDIF
          !
       ENDDO
    ENDDO
    !
    !
    ! The number of maximum vegetation map points in the GCM grid is estimated.
    ! Some lmargin is taken.
    !
    IF (is_root_prc) THEN
       nix=INT(MAXVAL(resolution_g(:,1))/MAXVAL(resol_lu(:,:,1)))+2
       njx=INT(MAXVAL(resolution_g(:,2))/MAXVAL(resol_lu(:,:,2)))+2
       nbvmax = nix*njx
    ENDIF
    CALL bcast(nbvmax)
    !
    callsign="Slope map"
    ok_interpol = .FALSE.
    DO WHILE ( .NOT. ok_interpol )
       !
       WRITE(numout,*) "Projection arrays for ",callsign," : "
       WRITE(numout,*) "nbvmax = ",nbvmax
       !
       ALLOC_ERR=-1
       ALLOCATE(sub_index(nbpt,nbvmax,2), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_index : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_index(:,:,:)=0
       ALLOC_ERR=-1
       ALLOCATE(sub_area(nbpt,nbvmax), STAT=ALLOC_ERR)
       IF (ALLOC_ERR/=0) THEN
          WRITE(numout,*) "ERROR IN ALLOCATION of sub_area : ",ALLOC_ERR
          STOP 
       ENDIF
       sub_area(:,:)=zero
       !
       CALL aggregate_p(nbpt, lalo, neighbours, resolution, contfrac, &
            &                iml, jml, lon_rel, lat_rel, mask, callsign, &
            &                nbvmax, sub_index, sub_area, ok_interpol)
       
       IF (.NOT. ok_interpol ) THEN
          IF (long_print) WRITE(numout,*) 'nbvmax will be increased from ',nbvmax,' to ', nbvmax*2
          DEALLOCATE(sub_area)
          DEALLOCATE(sub_index)
          nbvmax = nbvmax * 2
       END IF
    END DO
    !
    !
    DO ib = 1, nbpt
      !-
      !- Reinfiltration coefficient due to the slope: Calculation with parameteres maxlope_ro 
      !-
      slopecoef = zero

      ! Initialize last index to the highest possible 
      idi_last=nbvmax
      DO idi=1, nbvmax
         ! Leave the do loop if all sub areas are treated, sub_area <= 0
         IF ( sub_area(ib,idi) <= zero ) THEN
            ! Set last index to the last one used
            idi_last=idi-1
            ! Exit do loop
            EXIT
         END IF

         ip = sub_index(ib,idi,1)
         jp = sub_index(ib,idi,2)

         slopecoef = slopecoef + MIN(slopemap(ip,jp)/slope_noreinf, un) * sub_area(ib,idi)
      ENDDO

      IF ( idi_last >= 1 ) THEN
         reinf_slope(ib) = un - slopecoef / SUM(sub_area(ib,1:idi_last)) 
      ELSE
         reinf_slope(ib) = slope_default
      ENDIF
    ENDDO
    !
    !
    WRITE(numout,*) 'Interpolation Done'
    !
    !
    DEALLOCATE(slopemap)
    DEALLOCATE(sub_index)
    DEALLOCATE(sub_area)
    DEALLOCATE(mask)
    DEALLOCATE(lon_lu)
    DEALLOCATE(lat_lu)
    DEALLOCATE(lon_rel)
    DEALLOCATE(lat_rel)
    !
    !
    RETURN
    !
  END SUBROUTINE slowproc_slope

!! ================================================================================================================================
!! SUBROUTINE 	: get_vegcorr
!!
!>\BRIEF         The "get_vegcorr" routine defines the table of correspondence
!!               between the 94 Olson vegetation types and the 13 Plant Functional Types known 
!!               by SECHIBA and STOMATE. Used by slowproc for the old interpolation.
!!
!!\DESCRIPTION : get_vegcorr is needed if you use the old_map carteveg5km.nc. \n
!!               Usually SECHIBA can run with a 5kmx5km map which is derived from the IGBP one. \n
!!               We assume that we have a classification in 94 types. This is Olson one modified by Nicolas Viovy.\n
!!               ORCHIDEE has to convert the Olson vegetation types into PFTs for the run (interpolation step).\n
!!               Each Olson matches to a combination of fractions of one or several PFTs.\n
!!               This routine uses the same process for the non-biospheric map (not used).\n
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::vegcorr, ::nobiocorr.
!!
!! REFERENCE(S)	: 
!! - Olson, J.S., J.A. Watts, and L.J. Allison., 1983. 
!! "Carbon in Live Vegetation of Major World Ecosystems." 
!! Report ORNL-5862. Oak Ridge National Laboratory, Oak Ridge, Tennessee. 
!! - Olson, J.S., J.A. Watts, and L.J. Allison., 1985. 
!! "Major World Ecosystem Complexes Ranked by Carbon in Live Vegetation: A Database." 
!! NDP-017. Carbon Dioxide Information Center, Oak Ridge National Laboratory, Oak Ridge, Tennessee.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_vegcorr (nolson,vegcorr,nobiocorr)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    INTEGER(i_std),PARAMETER :: nolson94 = 94                       !! Number of Olson vegetation types (unitless)
    INTEGER(i_std),PARAMETER :: nvm13 = 15     !!!Arsene 13-02-2014 !! Number of PFTS of ORCHIDEE (unitless) 

    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in) :: nolson                             !! Number of Olson vegetation types (unitless)
    
    !! 0.2 Output variables 

    REAL(r_std),DIMENSION(nolson,nvm),INTENT(out) :: vegcorr        !! Correspondence array between Olson types and PFTS 
                                                                    !! (0-1, unitless)
    REAL(r_std),DIMENSION(nolson,nnobio),INTENT(out) :: nobiocorr   !! Correspondence array between non-vegetation types and nobio 
                                                                    !! types (lake,etc..) (0-1, unitless)

    !! 0.4 Local variable
    
    INTEGER(i_std) :: ib                                            !! Indice (unitless)
    
 !_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-
    IF (nolson /= nolson94) THEN
       WRITE(numout,*) nolson,nolson94
       CALL ipslerr_p(3,'get_vegcorr', '', '',&
            &                 'wrong number of OLSON vegetation types.') ! Fatal error
    ENDIF !(nolson /= nolson94)
    
    IF (nvm /= nvm13) THEN
       WRITE(numout,*) nvm,nvm13
       CALL ipslerr_p(3,'get_vegcorr', '', '',&
            &                 'wrong number of SECHIBA vegetation types.') ! Fatal error
    ENDIF !(nvm /= nvm13)

!! Arsene 04-01-13 REMOVE   ??????
    ! The carteveg5km cannot be used if the PFTs are not in the standard order
    DO ib = 1,nvm
       IF (pft_to_mtc(ib) /= ib ) THEN
          CALL ipslerr_p(3,'get_vegcorr','You have redefined the order of the 13 PFTS', & 
               &          'You can not use carteveg5km', 'Use the standard configuration of PFTS' )
       ENDIF
    ENDDO

    !-
    ! 1 set the indices of non-biospheric surface types to 0.
    !-
    nobiocorr(:,:) = zero
    !-
    ! 2 Here we construct the correspondance table
    !   between Olson and the following SECHIBA Classes.
    !   vegcorr(i,:)+nobiocorr(i,:) = 1.  for all i.
    !-
    ! The modified OLSON types found in file carteveg5km.nc
    ! created by Nicolas Viovy :
    !! Arsene 27-11-13 Add 2 PFT: Moss, Schrubs

    !  1 Urban
    vegcorr( 1,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  2 Cool low sparse grassland
    vegcorr( 2,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  3 Cold conifer forest
    vegcorr( 3,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  4 Cold deciduous conifer forest
    vegcorr( 4,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  5 Cool Deciduous broadleaf forest
    vegcorr( 5,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  6 Cool evergreen broadleaf forests
    vegcorr( 6,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !  7 Cool tall grasses and shrubs
    vegcorr( 7,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.3/)!
    !  8 Warm C3 tall grasses and shrubs
    vegcorr( 8,:) = &
         & (/0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.0, 0.3/)!
    !  9 Warm C4 tall grases and shrubs
    vegcorr( 9,:) = &
         & (/0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.0, 0.3/)!
    ! 10 Bare desert
    vegcorr(10,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 11 Cold upland tundra
    vegcorr(11,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.1, 0.2/)!!
    ! 12 Cool irrigated grassland
    vegcorr(12,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 13 Semi desert
    vegcorr(13,:) = &
         & (/0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 14 Glacier ice
    vegcorr(14,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    nobiocorr(14,iice) = 1.
    ! 15 Warm wooded wet swamp
    vegcorr(15,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.2, 0.0/)!!
    ! 16 Inland water
    vegcorr(16,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 17 sea water
    vegcorr(17,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 18 cool shrub evergreen
    vegcorr(18,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/)!
    ! 19 cold shrub deciduous
    vegcorr(19,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/)!
    ! 20 Cold evergreen forest and fields
    vegcorr(20,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.3/)!!
    ! 21 cool rain forest
    vegcorr(21,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1, 0.0/)!!
    ! 22 cold conifer boreal forest
    vegcorr(22,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1, 0.0/)!!
    ! 23 cool conifer forest
    vegcorr(23,:) = &
         & (/0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1, 0.0/)!!
    ! 24 warm mixed forest
    vegcorr(24,:) = &
         & (/0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0/)
    ! 25 cool mixed forest
    vegcorr(25,:) = &
         & (/0.0, 0.0, 0.0, 0.4, 0.0, 0.4, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1, 0.0/)!
    ! 26 cool broadleaf forest
    vegcorr(26,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 27 cool deciduous broadleaf forest
    vegcorr(27,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.3, 0.5, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.1, 0.0/)!
    ! 28 warm montane tropical forest
    vegcorr(28,:) = &
         & (/0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0/)
    ! 29 warm seasonal tropical forest
    vegcorr(29,:) = &
         & (/0.0, 0.5, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0/)
    ! 30 cool crops and towns
    vegcorr(30,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    ! 31 warm crops and towns
    vegcorr(31,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0/)
    ! 32 cool crops and towns
    vegcorr(32,:) = &
         & (/0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0/)
    ! 33 warm dry tropical woods
    vegcorr(33,:) = &
         & (/0.2, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 34 warm tropical rain forest
    vegcorr(34,:) = &
         & (/0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 35 warm tropical degraded forest
    vegcorr(35,:) = &
         & (/0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.2, 0.0/)!!
    ! 36 warm corn and beans cropland
    vegcorr(36,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0/)
    ! 37 cool corn and bean cropland
    vegcorr(37,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0/)
    ! 38 warm rice paddy and field
    vegcorr(38,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0/)
    ! 39 hot irrigated cropland
    vegcorr(39,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0/)
    ! 40 cool irrigated cropland
    vegcorr(40,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0/)
    ! 41 cold irrigated cropland
    vegcorr(41,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0/)
    ! 42 cool grasses and shrubs
    vegcorr(42,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5/)!!
    ! 43 hot and mild grasses and shrubs
    vegcorr(43,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0, 0.5/)
    ! 44 cold grassland
    vegcorr(44,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 45 Savanna (woods) C3
    vegcorr(45,:) = &
         & (/0.1, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.4/)!
    ! 46 Savanna woods C4
    vegcorr(46,:) = &
         & (/0.1, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.3/)!
    ! 47 Mire, bog, fen
    vegcorr(47,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.4, 0.3/)!!
    ! 48 Warm marsh wetland
    vegcorr(48,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.5, 0.2/)!!
    ! 49 cold marsh wetland
    vegcorr(49,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.5, 0.2/)!!
    ! 50 mediteraean scrub
    vegcorr(50,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.8/)!
    ! 51 Cool dry woody scrub
    vegcorr(51,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.7/)!
    ! 52 Warm dry evergreen woods
    vegcorr(52,:) = &
         & (/0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1/)!
    ! 53 Volcanic rocks
    vegcorr(53,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 54 sand desert
    vegcorr(54,:) = &
         & (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 55 warm semi desert shrubs
    vegcorr(55,:) = &
         & (/0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4/)!
    ! 56 cool semi desert shrubs
    vegcorr(56,:) = &
         & (/0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.4/)!
    ! 57 semi desert sage
    vegcorr(57,:) = &
         & (/0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.5/)!
    ! 58 Barren tundra
    vegcorr(58,:) = &
         & (/0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.3, 0.1/)!
    ! 59 cool southern hemisphere mixed forest
    vegcorr(59,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.1/)!
    ! 60 cool fields and woods
    vegcorr(60,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.2/)!
    ! 61 warm forest and filed
    vegcorr(61,:) = &
         & (/0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.1/)!
    ! 62 cool forest and field
    vegcorr(62,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.2/)!
    ! 63 warm C3 fields and woody savanna
    vegcorr(63,:) = &
         & (/0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.3/)!!
    ! 64 warm C4 fields and woody savanna
    vegcorr(64,:) = &
         & (/0.1, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.2/)!
    ! 65 cool fields and woody savanna
    vegcorr(65,:) = &
         & (/0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.3/)!
    ! 66 warm succulent and thorn scrub
    vegcorr(66,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.4, 0.2/)!
    ! 67 cold small leaf mixed woods
    vegcorr(67,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.3, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.3/)!
    ! 68 cold deciduous and mixed boreal fores
    vegcorr(68,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.2/)!
    ! 69 cold narrow conifers
    vegcorr(69,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0/)!
    ! 70 cold wooded tundra
    vegcorr(70,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.4, 0.0, 0.0, 0.0, 0.1, 0.3/)!
    ! 71 cold heath scrub
    vegcorr(71,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.8/)!!
    ! 72 Polar and alpine desert
    vegcorr(72,:) = &
         & (/0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 73 warm Mangrove
    vegcorr(73,:) = &
         & (/0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 74 cool crop and water mixtures
    vegcorr(74,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.0, 0.0, 0.0/)
    ! 75 cool southern hemisphere mixed forest
    vegcorr(75,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.1/)!
    ! 76 cool moist eucalyptus
    vegcorr(76,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 77 warm rain green tropical forest
    vegcorr(77,:) = &
         & (/0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    ! 78 warm C3 woody savanna
    vegcorr(78,:) = &
         & (/0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.3/)!
    ! 79 warm C4 woody savanna
    vegcorr(79,:) = &
         & (/0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.2/)!
    ! 80 cool woody savanna
    vegcorr(80,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.4/)!!
    ! 81 cold woody savanna
    vegcorr(81,:) = &
         & (/0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.3/)!
    ! 82 warm broadleaf crops
    vegcorr(82,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0/)
    ! 83 warm C3 grass crops
    vegcorr(83,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0/)
    ! 84 warm C4 grass crops
    vegcorr(84,:) = &
         & (/0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0/)
    ! 85 cool grass crops
    vegcorr(85,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0/)
    ! 86 warm C3 crops grass,shrubs
    vegcorr(86,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.2/)!
    ! 87 cool crops,grass,shrubs
    vegcorr(87,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.5, 0.0, 0.0, 0.2/)!
    ! 88 warm evergreen tree crop
    vegcorr(88,:) = &
         & (/0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0/)
    ! 89 cool evergreen tree crop
    vegcorr(89,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 90 cold evergreen tree crop
    vegcorr(90,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 91 warm deciduous tree crop
    vegcorr(91,:) = &
         & (/0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0/)
    ! 92 cool deciduous tree crop
    vegcorr(92,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 93 cold deciduous tree crop
    vegcorr(93,:) = &
         & (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0/)
    ! 94 wet sclerophylic forest
    vegcorr(94,:) = &
         & (/0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    !-
    ! 3 Check the mapping for the Olson types which are going into the
    !   the veget and nobio array.
    !-
    DO ib=1,nolson
       !
       IF ( ABS(SUM(vegcorr(ib,:))+SUM(nobiocorr(ib,:))-1.0) &
            &       > EPSILON(1.0)) THEN
          WRITE(numout,*) 'Wrong correspondance for Olson type :', ib
          CALL ipslerr_p(3,'get_vegcorr', '', '',&
               &                 'Wrong correspondance for Olson type.') ! Fatal error
       ENDIF
       !
    ENDDO ! Loop over the # Olson type


  END SUBROUTINE get_vegcorr

!! ================================================================================================================================
!! SUBROUTINE 	: get_soilcorr
!!
!>\BRIEF         The "get_soilcorr" routine defines the table of correspondence
!!               between the Zobler types and the three texture types known by SECHIBA and STOMATE :
!!               silt, sand and clay. 
!!
!! DESCRIPTION : get_soilcorr is needed if you use soils_param.nc .\n
!!               The data from this file is then interpolated to the grid of the model. \n
!!               The aim is to get fractions for sand loam and clay in each grid box.\n
!!               This information is used for soil hydrology and respiration.
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S) : ::texfrac_table
!!
!! REFERENCE(S)	: 
!! - Zobler L., 1986, A World Soil File for global climate modelling. NASA Technical memorandum 87802. NASA 
!!   Goddard Institute for Space Studies, New York, U.S.A.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_soilcorr (nzobler,textfrac_table)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    INTEGER(i_std),PARAMETER :: nbtypes_zobler = 7                    !! Number of Zobler types (unitless)

    !! 0.1  Input variables
    
    INTEGER(i_std),INTENT(in) :: nzobler                              !! Size of the array (unitless)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION(nzobler,ntext),INTENT(out) :: textfrac_table !! Table of correspondence (0-1, unitless)
    
    !! 0.4 Local variables
    
    INTEGER(i_std) :: ib                                              !! Indice (unitless)
    
!_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-  
    IF (nzobler /= nbtypes_zobler) THEN 
       CALL ipslerr_p(3,'get_soilcorr', 'nzobler /= nbtypes_zobler',&
          &   'We do not have the correct number of classes', &
          &                 ' in the code for the file.')  ! Fatal error
    ENDIF

    !-
    ! 1. Textural fraction for : silt        sand         clay
    !-
    textfrac_table(1,:) = (/ 0.12, 0.82, 0.06 /)
    textfrac_table(2,:) = (/ 0.32, 0.58, 0.10 /)
    textfrac_table(3,:) = (/ 0.39, 0.43, 0.18 /)
    textfrac_table(4,:) = (/ 0.15, 0.58, 0.27 /)
    textfrac_table(5,:) = (/ 0.34, 0.32, 0.34 /)
    textfrac_table(6,:) = (/ 0.00, 1.00, 0.00 /)
    textfrac_table(7,:) = (/ 0.39, 0.43, 0.18 /)


    !-
    ! 2. Check the mapping for the Zobler types which are going into the ORCHIDEE textures classes 
    !-
    DO ib=1,nzobler ! Loop over # classes soil
       
       IF (ABS(SUM(textfrac_table(ib,:))-1.0) > EPSILON(1.0)) THEN ! The sum of the textural fractions should not exceed 1 !
          WRITE(numout,*) &
               &     'Error in the correspondence table', &
               &     ' sum is not equal to 1 in', ib
          WRITE(numout,*) textfrac_table(ib,:)
          CALL ipslerr_p(3,'get_soilcorr', 'SUM(textfrac_table(ib,:)) /= 1.0',&
               &                 '', 'Error in the correspondence table') ! Fatal error
       ENDIF
       
    ENDDO ! Loop over # classes soil

    
  END SUBROUTINE get_soilcorr

!! ================================================================================================================================
!! FUNCTION 	: tempfunc
!!
!>\BRIEF        ! This function interpolates value between ztempmin and ztempmax
!! used for lai detection. 
!!
!! DESCRIPTION   : This subroutine calculates a scalar between 0 and 1 with the following equation :\n
!!                 \latexonly
!!                 \input{constantes_veg_tempfunc.tex}
!!                 \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : tempfunc_result
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  FUNCTION tempfunc (temp_in) RESULT (tempfunc_result)


    !! 0. Variables and parameters declaration

    REAL(r_std),PARAMETER    :: ztempmin=273._r_std   !! Temperature for laimin (K)
    REAL(r_std),PARAMETER    :: ztempmax=293._r_std   !! Temperature for laimax (K)
    REAL(r_std)              :: zfacteur              !! Interpolation factor   (K^{-2})

    !! 0.1 Input variables

    REAL(r_std),INTENT(in)   :: temp_in               !! Temperature (K)

    !! 0.2 Result

    REAL(r_std)              :: tempfunc_result       !! (unitless)
    
!_ ================================================================================================================================

    !! 1. Define a coefficient
    zfacteur = un/(ztempmax-ztempmin)**2
    
    !! 2. Computes tempfunc
    IF     (temp_in > ztempmax) THEN
       tempfunc_result = un
    ELSEIF (temp_in < ztempmin) THEN
       tempfunc_result = zero
    ELSE
       tempfunc_result = un-zfacteur*(ztempmax-temp_in)**2
    ENDIF !(temp_in > ztempmax)


  END FUNCTION tempfunc

!!
!================================================================================================================================
!! SUBROUTINE   : slowproc_read_veg_restfile
!!
!>\BRIEF         read a vegetation map from a restart file, and interpolate if
!necessary to current model grid
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: none
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE slowproc_read_veg_restfile(nbpt, lalo, neighbours, resolution, contfrac, veget, frac_nobio )
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)          ! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,8)    ! Vector of neighbours for each grid point 
    ! (1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)    ! The size in km of each grid-box in X and Y
    REAL(r_std),DIMENSION (nbpt), INTENT (in) :: contfrac         !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  veget(nbpt,nvm)         ! Vegetation fractions
    REAL(r_std), INTENT(out)    ::  frac_nobio(nbpt,nnobio) ! Fraction of the mesh which is covered by ice, lakes, ...
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid, ib, jp, it, jj, jv, im, jm
    INTEGER(i_std) :: il, ils, ip, ix, iy, imin, jmin, ier
    REAL(r_std) :: dlon, dlonmin, dlat, dlatmin
    REAL(r_std), PARAMETER :: maxmargin = 5.
    !
    !  0.4 allocatable
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:) :: veget_max_in
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) ::  frac_nobio_in
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:) :: xx, yy
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: x, y

    ! get filename
    CALL getin_p('vegetation_map_restartfile_name', filename)

    ! load the sechiba restart file, and read in values of veget_max and
    ! frac_nobio
    CALL flininfo(filename, iml, jml, lml, tml, fid)
    !! for now, assume hard-coded  96 x 72 gridpoints because iml and jml
    !readings are wrong
    im = 96
    jm = 72
    ALLOCATE(xx(im, jm))
    ALLOCATE(yy(im,jm))
    ALLOCATE (x(im),y(jm))
    ALLOCATE(veget_max_in(im,jm,nvm))
    ALLOCATE(frac_nobio_in(im,jm))
    WRITE(*,*) 'cdk debug: '
    WRITE(*,*) 'filename: ', filename
    WRITE(*,*) 'im, jm, lml, tml, fid: ', im, jm, lml, tml, fid

    CALL flinget(fid, 'nav_lon', im, jm, 1, 1, 1, 1, xx)
    WRITE(*,*) 'cdk debug read nav_lon '
    CALL flinget(fid, 'nav_lat', im, jm, 1, 1, 1, 1, yy)
    WRITE(*,*) 'cdk debug read nav_lat '
    CALL flinget(fid, 'veget_max', im, jm, nvm, 1, 1, 1, veget_max_in)  !! these need to be double-precision
    WRITE(*,*) 'cdk debug read veget_max '
    CALL flinget(fid, 'frac_nobio', im, jm, 1, 1, 1, 1, frac_nobio_in)  !! these need to be double-precision
    WRITE(*,*) 'cdk debug read frac_nobio '

    x(:) = xx(:,1)
    y(:) = yy(1,:)

    ! prendre la valeur la plus proche
    DO ip = 1, nbpt
       dlonmin = HUGE(1.)
       DO ix = 1,im
          dlon = MIN( ABS(lalo(ip,2)-x(ix)), ABS(lalo(ip,2)+360.-x(ix)), ABS(lalo(ip,2)-360.-x(ix)) )
          IF ( dlon .LT. dlonmin ) THEN
             imin = ix
             dlonmin = dlon
          ENDIF
       ENDDO
       dlatmin = HUGE(1.)
       DO iy = 1,jm
          dlat = ABS(lalo(ip,1)-y(iy))
          IF ( dlat .LT. dlatmin ) THEN
             jmin = iy
             dlatmin = dlat
          ENDIF
       ENDDO
       WRITE(*,*) 'cdk debug found: ip,  dlonmin, dlatmin: ', ip, dlonmin,dlatmin
       ! imin, jmin, x(imin), y(jmin):,
       ! if nothing is close, then set as all nobio
       IF ( ( dlonmin .LE. maxmargin ) .AND. ( dlatmin .LE. maxmargin ) ) THEN
          veget(ip,:) = veget_max_in(imin,jmin,:)
          frac_nobio(ip,:) = frac_nobio_in(imin,jmin)
       ELSE
          veget(ip,:) = zero
          frac_nobio(ip,:) = un
       ENDIF
       WRITE(*,*) 'cdk debug: ip,  veget, frac_nobio: ', veget(ip,:),frac_nobio(ip,:)

    ENDDO

    DEALLOCATE(xx)
    DEALLOCATE(yy)
    DEALLOCATE(x)
    DEALLOCATE(y)
    DEALLOCATE(veget_max_in)
    DEALLOCATE(frac_nobio_in)
    !
    RETURN
    !
  END SUBROUTINE slowproc_read_veg_restfile


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_checkveget
!!
!>\BRIEF         To verify the consistency of the various fractions defined within the grid box after having been
!!               been updated by STOMATE or the standard procedures.
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: none
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!
  SUBROUTINE slowproc_checkveget(nbpt, frac_nobio, veget_max, veget, soiltile)
    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                      :: nbpt       ! Number of points for which the data needs to be interpolated
    REAL(r_std),DIMENSION (nbpt,nnobio), INTENT(in) :: frac_nobio ! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (nbpt,nvm), INTENT(in)    :: veget_max  ! Maximum fraction of vegetation type including none biological fraction (unitless) 
    REAL(r_std),DIMENSION (nbpt,nvm), INTENT(in)    :: veget      ! Vegetation fractions
    REAL(r_std),DIMENSION (nbpt,nstm), INTENT(in)   :: soiltile   ! Fraction of soil tiles in the gridbox (unitless)
    !
    !  0.2 OUTPUT
    !

    !
    !  0.3 LOCAL
    !
    INTEGER(i_std) :: ji, jn, jv, js
    REAL(r_std)  :: epsilocal
    REAL(r_std)  :: totfrac
    CHARACTER(len=80) :: str1, str2
    INTEGER(i_std), PARAMETER :: ibaresoil=1
    !
!_ ================================================================================================================================
    !
    !
    ! There is some margin added as the computing errors might bring us above EPSILON(un)
    !
    epsilocal = EPSILON(un)*1000.
    !
    ! 1.0 Verify that none of the fractions are smaller than min_vegfrac, without beeing zero.
    !
!    DO ji=1,nbpt
!       DO jn=1,nnobio
!          IF ( frac_nobio(ji,jn) > epsilocal .AND. frac_nobio(ji,jn) < min_vegfrac ) THEN
!             WRITE(str1,'("Occurs on grid box", I8," and nobio type ",I3 )') ji, jn
!             WRITE(str2,'("The small value obtained is ", E14.4)') frac_nobio(ji,jn)
!             CALL ipslerr_p (3,'slowproc_checkveget', &
!                   "frac_nobio is larger than zero but smaller than min_vegfrac.", str1, str2)
!          ENDIF
!       ENDDO
!       DO jv=1,nvm
!          IF ( veget_max(ji,jv) > epsilocal .AND. veget_max(ji,jv) < min_vegfrac ) THEN
!             WRITE(str1,'("Occurs on grid box", I8," and nobio type ",I3 )') ji, jn
!             WRITE(str2,'("The small value obtained is ", E14.4)') veget_max(ji,jv)
!             CALL ipslerr_p (3,'slowproc_checkveget', &
!                   "veget_max is larger than zero but smaller than min_vegfrac.", str1, str2)
!          ENDIF
!       ENDDO
!    ENDDO
    !
    ! 2.0 verify that with all the fractions we cover the entire grid box  
    !
    DO ji=1,nbpt
       totfrac = zero
       DO jn=1,nnobio
          totfrac = totfrac + frac_nobio(ji,jn)
       ENDDO
       DO jv=1,nvm
          totfrac = totfrac + veget_max(ji,jv)
       ENDDO
       IF ( ABS(totfrac - un) > epsilocal) THEN
             WRITE(str1,'("This occurs on grid box", I8)') ji
             WRITE(str2,'("The sum over all fraction and error are ", E14.4, E14.4)') totfrac, ABS(totfrac - un)
             CALL ipslerr_p (3,'slowproc_checkveget', &
                   "veget_max + frac_nobio is not equal to 1.", str1, str2)
             WRITE(*,*) "EPSILON =", epsilocal 
       ENDIF
    ENDDO
    !
    ! 3.0 Verify that veget is smaller or equal to veget_max
    !
    DO ji=1,nbpt
       DO jv=1,nvm
          IF ( jv == ibaresoil ) THEN
             IF ( ABS(veget(ji,jv) - veget_max(ji,jv)) > epsilocal ) THEN
                WRITE(str1,'("This occurs on grid box", I8)') ji
                WRITE(str2,'("The difference is ", E14.4)') veget(ji,jv) - veget_max(ji,jv)
                CALL ipslerr_p (3,'slowproc_checkveget', &
                     "veget is not equal to veget_max on bare soil.", str1, str2)
             ENDIF
          ELSE
             IF ( veget(ji,jv) > veget_max(ji,jv) ) THEN
                WRITE(str1,'("This occurs on grid box", I8)') ji
                WRITE(str2,'("The values for veget and veget_max :", F8.4, F8.4)') veget(ji,jv), veget_max(ji,jv)
                CALL ipslerr_p (3,'slowproc_checkveget', &
                     "veget is greater than veget_max.", str1, str2)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    ! 4.0 Test tot_bare_soil in relation to the other variables
    !
    DO ji=1,nbpt
       totfrac = zero
       DO jv=1,nvm
          totfrac = totfrac + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
       ! add the bare soil fraction to totfrac
       totfrac = totfrac + veget(ji,ibaresoil)
       ! do the test
       IF ( ABS(totfrac - tot_bare_soil(ji)) > epsilocal ) THEN
          WRITE(str1,'("This occurs on grid box", I8)') ji
          WRITE(str2,'("The values for tot_bare_soil, tot frac and error :", F8.4, F8.4, E14.4)') &
               &  tot_bare_soil(ji), totfrac, ABS(totfrac - tot_bare_soil(ji))
          CALL ipslerr_p (3,'slowproc_checkveget', &
               "tot_bare_soil does not correspond to the total bare soil fraction.", str1, str2)
       ENDIF
    ENDDO
    !
    ! 5.0 Test that soiltile has the right sum
    !
    DO ji=1,nbpt
       totfrac = zero
       DO js=1,nstm
          totfrac = totfrac + soiltile(ji,js)
       ENDDO
       IF ( ABS(totfrac - un) > epsilocal ) THEN
          WRITE(str1,'("This occurs on grid box", I8)') ji
          WRITE(str2,'("The sum of soiltile and error are :", F8.4, E14.4)') soiltile, ABS(totfrac - un)
          CALL ipslerr_p (3,'slowproc_checkveget', &
               "soiltile does not sum-up to one.", str1, str2)
       ENDIF
    ENDDO
    !
  END SUBROUTINE slowproc_checkveget
    

END MODULE slowproc
