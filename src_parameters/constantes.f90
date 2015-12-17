! =================================================================================================================================
! MODULE       : constantes
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        "constantes" module contains subroutines to initialize most of the exernalized parameters. This module
!!              also make a use to the module constantes_var where the parameters are declared.
!!
!!\n DESCRIPTION: This module contains subroutines to initialize most of the exernalized parameters. This module
!!                also make a use to the module constantes_var where the parameters are declared.\n
!!                This module can be used to acces the subroutines and the constantes. The constantes declarations
!!                can also be used seperatly with "USE constantes_var".
!!
!! RECENT CHANGE(S): Didier Solyga : This module contains now all the externalized parameters of ORCHIDEE 
!!                   listed by modules which are not pft-dependent  
!!                   Josefine Ghattas 2013 : The declaration part has been extracted and moved to module constates_var
!!
!! REFERENCE(S)	: 
!! - Louis, Jean-Francois (1979), A parametric model of vertical eddy fluxes in the atmosphere. 
!! Boundary Layer Meteorology, 187-202.
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2014-09-04 14:46:14 +0200 (Thu, 04 Sep 2014) $
!! $Revision: 2282 $
!! \n
!_ ================================================================================================================================

MODULE constantes

  USE constantes_var
  USE defprec
  USE ioipsl_para, ONLY : getin_p, ipslerr_p
  USE mod_orchidee_para, ONLY : numout

  IMPLICIT NONE

CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : activate_sub_models
!!
!>\BRIEF         This subroutine reads the flags in the configuration file to
!! activate some sub-models like routing, irrigation, fire, herbivory, ...  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

   SUBROUTINE activate_sub_models(active_flags)

     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.1 Input variables

     TYPE(control_type), INTENT(in) :: active_flags     !! What parts of the code are activated ?

     !! 0.4 Local variables

     LOGICAL, SAVE ::  first_call = .TRUE.             !! To keep first call trace (true/false)
!$OMP THREADPRIVATE(first_call)

!_ ================================================================================================================================

     IF (first_call) THEN 
           
        IF (active_flags%ok_stomate) THEN

           !Config Key   = HERBIVORES
           !Config Desc  = herbivores allowed?
           !Config If    = OK_STOMATE 
           !Config Def   = n
           !Config Help  = With this variable, you can determine
           !Config         if herbivores are activated
           !Config Units = [FLAG]
           CALL getin_p('HERBIVORES', ok_herbivores)
           !
           !Config Key   = TREAT_EXPANSION
           !Config Desc  = treat expansion of PFTs across a grid cell?
           !Config If    = OK_STOMATE 
           !Config Def   = n
           !Config Help  = With this variable, you can determine
           !Config         whether we treat expansion of PFTs across a
           !Config         grid cell.
           !Config Units = [FLAG]
           CALL getin_p('TREAT_EXPANSION', treat_expansion)
           !
           !Config Key   = LPJ_GAP_CONST_MORT
           !Config Desc  = prescribe mortality if not using DGVM?
           !Config If    = OK_STOMATE
           !Config Def   = y
           !Config Help  = set to TRUE if constant mortality is to be activated
           !Config         ignored if DGVM=true!
           !Config Units = [FLAG]
           CALL getin_p('LPJ_GAP_CONST_MORT', lpj_gap_const_mort)
           !
           !Config Key   = HARVEST_AGRI
           !Config Desc  = Harvest model for agricultural PFTs.
           !Config If    = OK_STOMATE 
           !Config Def   = y
           !Config Help  = Compute harvest above ground biomass for agriculture.
           !Config         Change daily turnover.
           !Config Units = [FLAG]
           CALL getin_p('HARVEST_AGRI', harvest_agri)
           !
           !Config Key   = FIRE_DISABLE
           !Config Desc  = no fire allowed
           !Config If    = OK_STOMATE 
           !Config Def   = n
           !Config Help  = With this variable, you can allow or not
           !Config         the estimation of CO2 lost by fire
           !Config Units = [FLAG]
           CALL getin_p('FIRE_DISABLE', disable_fire)
           !
           !Config Key   = SPINUP_ANALYTIC
           !Config Desc  = Activation of the analytic resolution of the spinup.
           !Config If    = OK_STOMATE
           !Config Def   = n
           !Config Help  = Activate this option if you want to solve the spinup by the Gauss-Jordan method.
           !Config Units = BOOLEAN    
           CALL getin_p('SPINUP_ANALYTIC',spinup_analytic)

        ENDIF

        !
        ! Check consistency (see later)
        !
!!$        IF(.NOT.(ok_routing) .AND. (doirrigation .OR. dofloodplains)) THEN
!!$           CALL ipslerr_p(2,'activate_sub_models', &
!!$               &     'Problem :you tried to activate the irrigation and floodplains without activating the routing',&
!!$               &     'Are you sure ?', &
!!$               &     '(check your parameters).')
!!$        ENDIF
       
!!$        IF(.NOT.(ok_stomate) .AND. (ok_herbivores .OR. treat_expansion .OR. lpj_gap_const_mort &
!!$            & .OR. harvest_agri .OR. disable_fire)) THEN
!!$          CALL ipslerr_p(2,'activate_sub_models', &
!!$               &     'Problem : try to activate the following options : herbivory, treat_expansion, fire,',&
!!$               &     'harvest_agri and constant mortality without stomate activated.',&
!!$               &     '(check your parameters).')
!!$        ENDIF
            
        first_call =.FALSE.

     ENDIF

   END SUBROUTINE activate_sub_models
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : veget_config
!!
!>\BRIEF         This subroutine reads the flags controlling the configuration for
!! the vegetation : impose_veg, veget_mpa, lai_map, etc...       
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

   SUBROUTINE veget_config

     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.4 Local variables  

     LOGICAL, SAVE ::  first_call = .TRUE.        !! To keep first call trace (true/false)  
!$OMP THREADPRIVATE(first_call)

!_ ================================================================================================================================
     
     IF (first_call) THEN 

        !Config Key   = AGRICULTURE
        !Config Desc  = agriculture allowed?
        !Config If    = OK_SECHIBA or OK_STOMATE
        !Config Def   = y
        !Config Help  = With this variable, you can determine
        !Config         whether agriculture is allowed
        !Config Units = [FLAG]
        CALL getin_p('AGRICULTURE', agriculture)
        !
        !Config Key   = IMPOSE_VEG
        !Config Desc  = Should the vegetation be prescribed ?
        !Config If    = OK_SECHIBA or OK_STOMATE
        !Config Def   = n
        !Config Help  = This flag allows the user to impose a vegetation distribution
        !Config         and its characteristics. It is espacially interesting for 0D
        !Config         simulations. On the globe it does not make too much sense as
        !Config         it imposes the same vegetation everywhere
        !Config Units = [FLAG]
        CALL getin_p('IMPOSE_VEG', impveg)

        IF (impveg) THEN
           !Config Key   = IMPOSE_SOILT
           !Config Desc  = Should the soil type be prescribed ?
           !Config Def   = n
           !Config If    = IMPOSE_VEG
           !Config Help  = This flag allows the user to impose a soil type distribution.
           !Config         It is espacially interesting for 0D
           !Config         simulations. On the globe it does not make too much sense as
           !Config         it imposes the same soil everywhere
           !Config Units = [FLAG]
           CALL getin_p('IMPOSE_SOILT', impsoilt)     
        ENDIF

        !Config Key   = LAI_MAP
        !Config Desc  = Read the LAI map
        !Config If    = OK_SECHIBA or OK_STOMATE
        !Config Def   = n
        !Config Help  = It is possible to read a 12 month LAI map which will
        !Config         then be interpolated to daily values as needed.
        !Config Units = [FLAG]
        CALL getin_p('LAI_MAP',read_lai)

        IF(read_lai) THEN
           !Config Key   = SLOWPROC_LAI_OLD_INTERPOL
           !Config Desc  = Flag to use old "interpolation" of LAI
           !Config If    = LAI_MAP
           !Config Def   = n
           !Config Help  = If you want to recover the old (ie orchidee_1_2 branch) 
           !Config         "interpolation" of LAI map.
           !Config Units = [FLAG]
           CALL getin_p('SLOWPROC_LAI_OLD_INTERPOL',old_lai)
        ENDIF
 
        !Config Key   = LAND_USE
        !Config Desc  = Read a land_use vegetation map
        !Config If    = OK_SECHIBA or OK_STOMATE
        !Config Def   = y
        !Config Help  = pft values are needed, max time axis is 293
        !Config Units = [FLAG]
        CALL getin_p('LAND_USE',land_use)

        IF(land_use) THEN
           !Config Key   = VEGET_REINIT
           !Config Desc  = booleen to indicate that a new LAND USE file will be used.
           !Config If    = LAND_USE
           !Config Def   = y
           !Config Help  = The parameter is used to bypass veget_year count 
           !Config         and reinitialize it with VEGET_YEAR parameter.
           !Config         Then it is possible to change LAND USE file.
           !Config Units = [FLAG] 
           CALL getin_p('VEGET_REINIT', veget_reinit)
           !
           !Config Key   = LAND_COVER_CHANGE
           !Config Desc  = treat land use modifications
           !Config If    = LAND_USE
           !Config Def   = n
           !Config Help  = With this variable, you can use a Land Use map
           !Config         to simulate anthropic modifications such as 
           !Config         deforestation.
           !Config Units = [FLAG] 
           CALL getin_p('LAND_COVER_CHANGE', lcchange)
           !
           !Config Key   = VEGET_YEAR
           !Config Desc  = Year of the land_use vegetation map to be read
           !Config If    = LAND_USE
           !Config Def   = 1
           !Config Help  = First year for landuse vegetation (2D map by pft).
           !Config         If VEGET_YEAR is set to 0, this means there is no time axis.
           !Config Units = [FLAG] 
           CALL getin_p('VEGET_YEAR', veget_year_orig)
        ENDIF

        IF(.NOT. impveg .AND. .NOT. land_use) THEN
           !Config Key   = SLOWPROC_VEGET_OLD_INTERPOL
           !Config Desc  = Flag to use old "interpolation" of vegetation map.
           !Config If    = NOT(IMPOSE_VEG) and NOT(LAND_USE)
           !Config Def   = n
           !Config Help  = If you want to recover the old (ie orchidee_1_2 branch) 
           !Config         "interpolation" of vegetation map.
           !Config Units = [FLAG] 
           CALL getin_p('SLOWPROC_VEGET_OLD_INTERPOL',old_veget)
         ENDIF  

         !
         ! Check consistency
         !
         ! 1. You have to activate agriculture and land_use
         IF ( .NOT. agriculture .AND. land_use ) THEN 
            CALL ipslerr_p(2,'veget_config', &
                 &     'Problem with agriculture desactivated and Land Use activated.',&
                 &     'Are you sure ?', &
                 &     '(check your parameters).')
         ENDIF


        first_call = .FALSE.

     ENDIF

!!$        ! DS : Add warning in case of a wrong configuration (need to be discussed)
!!$        ! 2. 
!!$        IF (.NOT.(read_lai) .AND. old_lai) THEN
!!$           CALL ipslerr_p(2,'veget_config', &
!!$               &     'Problem with lai_map desactivated and old_lai activated.',&
!!$               &     'Are you sure ?', &
!!$               &     '(check your parameters).')
!!$        ENDIF
!!$    
!!$        ! 3.
!!$        IF ((impveg .OR. land_use) .AND. old_veget) THEN
!!$           CALL ipslerr_p(2,'veget_config', &
!!$                &     'Problem : try to use the old interpolation with a land use map or in impose_veg.',&
!!$                &     'Are you sure ?', &
!!$                &     '(check your parameters).')
!!$        ENDIF
!!$
!!$        ! 4.
!!$        IF ( .NOT.(impveg) .AND. impsoilt) THEN
!!$           CALL ipslerr_p(2,'veget_config', &
!!$               &     'Problem : try to activate impose_soilt without activating impose_veg.',&
!!$               &     'Are you sure ?', &
!!$               &     '(check your parameters).')
!!$        ENDIF
!!$
!!$        ! 5.
!!$        IF (.NOT.(land_use) .AND. (veget_reinit)) THEN
!!$           CALL ipslerr_p(2,'veget_config', &
!!$                &     'Problem : try to use a land_use map without activating land_use.',&
!!$                &     'Are you sure ?', &
!!$                &     '(check your parameters).')        
!!$        ENDIF
!!$
!!$        ! 6.
!!$        IF (.NOT.(land_use) .AND. lcchange) THEN
!!$           CALL ipslerr_p(2,'veget_config', &
!!$                &     'Problem : lcchange is activated without activating land_use.',&
!!$                &     'Are you sure ?', &
!!$                &     '(check your parameters).')        
!!$        ENDIF
           
   END SUBROUTINE veget_config
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : veget_config
!!
!>\BRIEF         This subroutine reads in the configuration file the imposed values of the parameters for all SECHIBA modules.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

   SUBROUTINE config_sechiba_parameters

     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.4 Local variables
     
     LOGICAL, SAVE ::  first_call = .TRUE.    !! To keep first call trace (true/false)
!$OMP THREADPRIVATE(first_call)

!_ ================================================================================================================================
     
     IF(first_call) THEN 
        
        ! Global : parameters used by many modules
        !
        !Config Key   = MAXMASS_SNOW
        !Config Desc  = The maximum mass of a snow
        !Config If    = OK_SECHIBA or HYDROL_CWRR
        !Config Def   = 3000.
        !Config Help  = 
        !Config Units = [kg/m^2]  
        CALL getin_p('MAXMASS_SNOW',maxmass_snow)
        !
        !Config Key   = SNOWCRI
        !Config Desc  = Sets the amount above which only sublimation occures 
        !Config If    = OK_SECHIBA or HYDROL_CWRR
        !Config Def   = 1.5
        !Config Help  = 
        !Config Units = [kg/m^2]  
        CALL getin_p('SNOWCRI',snowcri)
        !
        !! Initialization of sneige
        sneige = snowcri/mille
        !
        !Config Key   = MIN_WIND
        !Config Desc  = Minimum wind speed
        !Config If    = OK_SECHIBA
        !Config Def   = 0.1
        !Config Help  = 
        !Config Units = [m/s]
        CALL getin_p('MIN_WIND',min_wind)
        !
        !Config Key   = MAX_SNOW_AGE
        !Config Desc  = Maximum period of snow aging 
        !Config If    = OK_SECHIBA
        !Config Def   = 50.
        !Config Help  = 
        !Config Units = [days?]
        CALL getin_p('MAX_SNOW_AGE',max_snow_age)
        !
        !Config Key   = SNOW_TRANS
        !Config Desc  = Transformation time constant for snow
        !Config If    = OK_SECHIBA
        !Config Def   = 0.3
        !Config Help  = 
        !Config Units = [m]   
        CALL getin_p('SNOW_TRANS',snow_trans)
 
        !-
        ! condveg
        !-
        !
        !Config Key   = Z0_OVER_HEIGHT
        !Config Desc  = to get z0 from height 
        !Config If    = OK_SECHIBA 
        !Config Def   = 1/16.
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('Z0_OVER_HEIGHT',z0_over_height)
        !
        !Config Key   = HEIGHT_DISPLACEMENT
        !Config Desc  = Magic number which relates the height to the displacement height.
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.75
        !Config Help  = 
        !Config Units = [m]  
        CALL getin_p('HEIGHT_DISPLACEMENT',height_displacement)
        !
        !Config Key   = Z0_BARE
        !Config Desc  = bare soil roughness length
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.01 
        !Config Help  = 
        !Config Units = [m]   
        CALL getin_p('Z0_BARE',z0_bare)
        !
        !Config Key   = Z0_ICE
        !Config Desc  = ice roughness length
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.001
        !Config Help  = 
        !Config Units = [m]   
        CALL getin_p('Z0_ICE',z0_ice)
        !
        !Config Key   = TCST_SNOWA
        !Config Desc  = Time constant of the albedo decay of snow
        !Config If    = OK_SECHIBA 
        !Config Def   = 5.0 
        !Config Help  = 
        !Config Units = [days]
        CALL getin_p('TCST_SNOWA',tcst_snowa)
        !
        !Config Key   = SNOWCRI_ALB
        !Config Desc  = Critical value for computation of snow albedo
        !Config If    = OK_SECHIBA
        !Config Def   = 10. 
        !Config Help  = 
        !Config Units = [cm]  
        CALL getin_p('SNOWCRI_ALB',snowcri_alb)
        !
        !
        !Config Key   = VIS_DRY
        !Config Desc  = The correspondance table for the soil color numbers and their albedo 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.27
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('VIS_DRY',vis_dry)
        !
        !Config Key   = NIR_DRY
        !Config Desc  = The correspondance table for the soil color numbers and their albedo 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.48, 0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20, 0.55
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('NIR_DRY',nir_dry)
        !
        !Config Key   = VIS_WET 
        !Config Desc  = The correspondance table for the soil color numbers and their albedo
        !Config If    = OK_SECHIBA  
        !Config Def   = 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('VIS_WET',vis_wet)
        !
        !Config Key   = NIR_WET
        !Config Desc  = The correspondance table for the soil color numbers and their albedo 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.31
        !Config Help  = 
        !Config Units = [-]    
        CALL getin_p('NIR_WET',nir_wet)
        !
        !Config Key   = ALBSOIL_VIS
        !Config Desc  = 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.18, 0.16, 0.16, 0.15, 0.12, 0.105, 0.09, 0.075, 0.25
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('ALBSOIL_VIS',albsoil_vis)
        !
        !Config Key   = ALBSOIL_NIR 
        !Config Desc  = 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.36, 0.34, 0.34, 0.33, 0.30, 0.25, 0.20, 0.15, 0.45
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('ALBSOIL_NIR',albsoil_nir)
        !-
        !
        !Config Key   = ALB_DEADLEAF 
        !Config Desc  = albedo of dead leaves, VIS+NIR 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.12, 0.35
        !Config Help  = 
        !Config Units = [-]     
        CALL getin_p('ALB_DEADLEAF',alb_deadleaf)
        !
        !Config Key   = ALB_ICE
        !Config Desc  = albedo of ice, VIS+NIR
        !Config If    = OK_SECHIBA
        !Config Def   = 0.60, 0.20
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('ALB_ICE',alb_ice)
        !
        ! Get the fixed snow albedo if needed
        !
        !Config Key   = CONDVEG_SNOWA
        !Config Desc  = The snow albedo used by SECHIBA
        !Config Def   = 1.E+20
        !Config if    = OK_SECHIBA
        !Config Help  = This option allows the user to impose a snow albedo.
        !Config         Default behaviour is to use the model of snow albedo
        !Config         developed by Chalita (1993).
        !Config Units = [-]
        CALL getin_p('CONDVEG_SNOWA',fixed_snow_albedo)
        !
        !Config Key   = ALB_BARE_MODEL
        !Config Desc  = Switch bare soil albedo dependent (if TRUE) on soil wetness
        !Config Def   = n
        !Config if    = OK_SECHIBA
        !Config Help  = If TRUE, the model for bare soil albedo is the old formulation.
        !Config         Then it depend on the soil dry or wetness. If FALSE, it is the 
        !Config         new computation that is taken, it is the mean of soil albedo.
        !Config Units = [FLAG]
        CALL getin_p('ALB_BARE_MODEL',alb_bare_model)
        !
        !Config Key   = Z0CDRAG_AVE
        !Config Desc  = Average method for z0
        !Config Def   = y
        !Config if    = OK_SECHIBA
        !Config Help  = If this flag is set to true (y) then the neutral Cdrag
        !Config         is averaged instead of the log(z0). This should be
        !Config         the prefered option. We still wish to keep the other
        !Config         option so we can come back if needed. If this is
        !Config         desired then one should set Z0CDRAG_AVE=n
        !Config Units = [FLAG]
        CALL getin_p('Z0CDRAG_AVE',z0cdrag_ave)
        !
        !Config Key   = IMPOSE_AZE
        !Config Desc  = Should the surface parameters be prescribed
        !Config Def   = n
        !Config if    = OK_SECHIBA
        !Config Help  = This flag allows the user to impose the surface parameters
        !Config         (Albedo Roughness and Emissivity). It is espacially interesting for 0D
        !Config         simulations. On the globe it does not make too much sense as
        !Config         it imposes the same vegetation everywhere
        !Config Units = [FLAG]
        CALL getin_p('IMPOSE_AZE',impaze)
        !
        IF(impaze) THEN
           !
           !Config Key   = CONDVEG_Z0
           !Config Desc  = Surface roughness
           !Config Def   = 0.15
           !Config If    = IMPOSE_AZE
           !Config Help  = Surface rougness to be used on the point if a 0-dim version
           !Config         of SECHIBA is used. Look at the description of the forcing  
           !Config         data for the correct value.
           !Config Units = [m]
           CALL getin_p('CONDVEG_Z0', z0_scal)  
           !
           !Config Key   = ROUGHHEIGHT
           !Config Desc  = Height to be added to the height of the first level
           !Config Def   = 0.0
           !Config If    = IMPOSE_AZE
           !Config Help  = ORCHIDEE assumes that the atmospheric level height is counted
           !Config         from the zero wind level. Thus to take into account the roughness
           !Config         of tall vegetation we need to correct this by a certain fraction
           !Config         of the vegetation height. This is called the roughness height in
           !Config         ORCHIDEE talk.
           !Config Units = [m] 
           CALL getin_p('ROUGHHEIGHT', roughheight_scal)
           ! 
           !Config Key   = CONDVEG_ALBVIS
           !Config Desc  = SW visible albedo for the surface
           !Config Def   = 0.25
           !Config If    = IMPOSE_AZE
           !Config Help  = Surface albedo in visible wavelengths to be used 
           !Config         on the point if a 0-dim version of SECHIBA is used. 
           !Config         Look at the description of the forcing data for 
           !Config         the correct value.
           !Config Units = [-]
           CALL getin_p('CONDVEG_ALBVIS', albedo_scal(ivis))
           !
           !Config Key   = CONDVEG_ALBNIR
           !Config Desc  = SW near infrared albedo for the surface
           !Config Def   = 0.25
           !Config If    = IMPOSE_AZE
           !Config Help  = Surface albedo in near infrared wavelengths to be used 
           !Config         on the point if a 0-dim version of SECHIBA is used. 
           !Config         Look at the description of the forcing data for 
           !Config         the correct value.
           !Config Units = [-]  
           CALL getin_p('CONDVEG_ALBNIR', albedo_scal(inir))
           !
           !Config Key   = CONDVEG_EMIS
           !Config Desc  = Emissivity of the surface for LW radiation
           !Config Def   = 1.0
           !Config If    = IMPOSE_AZE
           !Config Help  = The surface emissivity used for compution the LE emission
           !Config         of the surface in a 0-dim version. Values range between 
           !Config         0.97 and 1.. The GCM uses 0.98.
           !Config Units = [-] 
           CALL getin_p('CONDVEG_EMIS', emis_scal)
        ENDIF
        !- 
        ! Variables related to the explicitsnow module
        !-
        !Config Key = xansmax 
        !Config Desc = maximum snow albedo
        !Config If = OK_SECHIBA
        !Config Def = 0.85
        !Config Help = 
        !Config Units = [-] 
        CALL getin_p('XANSMAX',xansmax)
        !
        !Config Key = xansmin 
        !Config Desc = minimum snow albedo
        !Config If = OK_SECHIBA
        !Config Def = 0.50
        !Config Help = 
        !Config Units = [-] 
        CALL getin_p('XANSMIN',xansmin)
        !
        !Config Key = xans_todry 
        !Config Desc = albedo decay rate for the dry snow 
        !Config If = OK_SECHIBA
        !Config Def = 0.008
        !Config Help = 
        !Config Units = [S-1] 
        CALL getin_p('XANSDRY',xans_todry)
        !
        !Config Key = xans_t 
        !Config Desc = albedo decay rate for the wet snow 
        !Config If = OK_SECHIBA
        !Config Def = 0.24
        !Config Help = 
        !Config Units = [S-1] 
        CALL getin_p('XANS_T',xans_t)

        !Config Key = xrhosmax 
        !Config Desc = maximum snow density 
        !Config If = OK_SECHIBA
        !Config Def = 750
        !Config Help = 
        !Config Units = [-] 
        CALL getin_p('XRHOSMAX',xrhosmax)
        !
        !Config Key = xwsnowholdmax1
        !Config Desc = snow holding capacity 1
        !Config If = OK_SECHIBA
        !Config Def = 0.03
        !Config Help = 
        !Config Units = [-] 
        CALL getin_p('XWSNOWHOLDMAX1',xwsnowholdmax1)
        !
        !Config Key = xwsnowholdmax2
        !Config Desc = snow holding capacity 2
        !Config If = OK_SECHIBA
        !Config Def = 0.10
        !Config Help = 
        !Config Units = [-] 
        CALL getin_p('XWSNOWHOLDMAX2',xwsnowholdmax2)
        !
        !Config Key = xsnowrhohold 
        !Config Desc = snow density 
        !Config If = OK_SECHIBA
        !Config Def = 200.0
        !Config Help = 
        !Config Units = [kg/m3] 
        CALL getin_p('XSNOWRHOHOLD',xsnowrhohold)
        ! 
        !Config Key = ZSNOWTHRMCOND1
        !Config Desc = Thermal conductivity Coef 1
        !Config If = OK_SECHIBA
        !Config Def = 0.02 
        !Config Help = 
        !Config Units = [W/m/K] 
        CALL getin_p('ZSNOWTHRMCOND1',ZSNOWTHRMCOND1)
        !
        !Config Key = ZSNOWTHRMCOND2
        !Config Desc = Thermal conductivity Coef 2
        !Config If = OK_SECHIBA
        !Config Def = 2.5E-6
        !Config Help = 
        !Config Units = [W m5/(kg2 K)] 
        CALL getin_p('ZSNOWTHRMCOND2',ZSNOWTHRMCOND2)
        !
        !Config Key = ZSNOWTHRMCOND_AVAP
        !Config Desc = Thermal conductivity Coef 1 water vapor
        !Config If = OK_SECHIBA
        !Config Def = -0.06023
        !Config Help = 
        !Config Units = [W/m/K] 
        CALL getin_p('ZSNOWTHRMCOND_AVAP',ZSNOWTHRMCOND_AVAP)
        !
        !Config Key = ZSNOWTHRMCOND_BVAP
        !Config Desc = Thermal conductivity Coef 2 water vapor
        !Config If = OK_SECHIBA
        !Config Def = -2.5425
        !Config Help = 
        !Config Units = [W/m] 
        CALL getin_p('ZSNOWTHRMCOND_BVAP',ZSNOWTHRMCOND_BVAP)
        !
        !Config Key = ZSNOWTHRMCOND_CVAP
        !Config Desc = Thermal conductivity Coef 3 water vapor
        !Config If = OK_SECHIBA
        !Config Def = -289.99
        !Config Help = 
        !Config Units = [K] 
        CALL getin_p('ZSNOWTHRMCOND_CVAP',ZSNOWTHRMCOND_CVAP)

        !Snow compaction factors
        !Config Key = ZSNOWCMPCT_RHOD
        !Config Desc = Snow compaction coefficent
        !Config If = OK_SECHIBA
        !Config Def = 150.0
        !Config Help = 
        !Config Units = [kg/m3]
        CALL getin_p('ZSNOWCMPCT_RHOD',ZSNOWCMPCT_RHOD)

        !Config Key = ZSNOWCMPCT_ACM
        !Config Desc = Coefficent for the thermal conductivity
        !Config If = OK_SECHIBA
        !Config Def = 2.8e-6
        !Config Help = 
        !Config Units = [1/s]
        CALL getin_p('ZSNOWCMPCT_ACM',ZSNOWCMPCT_ACM)

        !Config Key = ZSNOWCMPCT_BCM
        !Config Desc = Coefficent for the thermal conductivity
        !Config If = OK_SECHIBA
        !Config Def = 0.04
        !Config Help = 
        !Config Units = [1/K]
        CALL getin_p('ZSNOWCMPCT_BCM',ZSNOWCMPCT_BCM)

        !Config Key = ZSNOWCMPCT_CCM
        !Config Desc = Coefficent for the thermal conductivity
        !Config If = OK_SECHIBA
        !Config Def = 460.
        !Config Help = 
        !Config Units = [m3/kg] 
        CALL getin_p('ZSNOWCMPCT_CCM',ZSNOWCMPCT_CCM)

        !Config Key = ZSNOWCMPCT_V0
        !Config Desc = Vapor coefficent for the thermal conductivity
        !Config If = OK_SECHIBA
        !Config Def = 3.7e7
        !Config Help = 
        !Config Units = [Pa/s]
        CALL getin_p('ZSNOWCMPCT_V0',ZSNOWCMPCT_V0)

        !Config Key = ZSNOWCMPCT_VT
        !Config Desc = Vapor coefficent for the thermal conductivity
        !Config If = OK_SECHIBA
        !Config Def = 0.081
        !Config Help = 
        !Config Units = [1/K]
        CALL getin_p('ZSNOWCMPCT_VT',ZSNOWCMPCT_VT)

        !Config Key = ZSNOWCMPCT_VR
        !Config Desc = Vapor coefficent for the thermal conductivity
        !Config If = OK_SECHIBA
        !Config Def = 0.018
        !Config Help = 
        !Config Units = [m3/kg]
        CALL getin_p('ZSNOWCMPCT_VR',ZSNOWCMPCT_VR)


        !Surface resistance
        !
        !Config Key = CB
        !Config Desc = Constant of the Louis scheme 
        !Config If = OK_SECHIBA
        !Config Def = 5.0
        !Config Help = 
        !Config Units = [-] 
        CALL getin_p('CB',cb)
        !
        !Config Key = CC
        !Config Desc = Constant of the Louis scheme 
        !Config If = OK_SECHIBA
        !Config Def = 5.0
        !Config Help = 
        !Config Units = [-] 
        CALL getin_p('CC',cc)
        !
        !Config Key = CD
        !Config Desc = Constant of the Louis scheme 
        !Config If = OK_SECHIBA
        !Config Def = 5.0
        !Config Help = 
        !Config Units = [-] 
        CALL getin_p('CD',cd)
        !
        !Config Key = RAYT_CSTE
        !Config Desc = Constant in the computation of surface resistance  
        !Config If = OK_SECHIBA
        !Config Def = 125
        !Config Help = 
        !Config Units = [W.m^{-2}] 
        CALL getin_p('RAYT_CSTE',rayt_cste)
        !
        !Config Key = DEFC_PLUS
        !Config Desc = Constant in the computation of surface resistance  
        !Config If = OK_SECHIBA
        !Config Def = 23.E-3
        !Config Help = 
        !Config Units = [K.W^{-1}] 
        CALL getin_p('DEFC_PLUS',defc_plus)
        !
        !Config Key = DEFC_MULT
        !Config Desc = Constant in the computation of surface resistance  
        !Config If = OK_SECHIBA
        !Config Def = 1.5
        !Config Help = 
        !Config Units = [K.W^{-1}] 
        CALL getin_p('DEFC_MULT',defc_mult)
        !

        !
        !-
        ! diffuco 
        !-
        !
        !Config Key   = NLAI
        !Config Desc  = Number of LAI levels
        !Config If    = OK_SECHIBA
        !Config Def   = 20
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('NLAI',nlai)
        !
        !Config Key   = LAIMAX
        !Config Desc  = Maximum LAI
        !Config If    = OK_SECHIBA
        !Config Def   = 
        !Config Help  = 
        !Config Units = [m^2/m^2]   
        CALL getin_p('LAIMAX',laimax)
        !
        !Config Key   = DEW_VEG_POLY_COEFF
        !Config Desc  = coefficients of the polynome of degree 5 for the dew
        !Config If    = OK_SECHIBA
        !Config Def   = 0.887773, 0.205673, 0.110112, 0.014843, 0.000824, 0.000017 
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('DEW_VEG_POLY_COEFF',dew_veg_poly_coeff)
        !
        !Config Key   = DOWNREGULATION_CO2
        !Config Desc  = Activation of CO2 downregulation
        !Config If    = OK_SECHIBA
        !Config Def   = n 
        !Config Help  = 
        !Config Units = [FLAG]   
        CALL getin_p('DOWNREGULATION_CO2',downregulation_co2)
        !
        !Config Key   = DOWNREGULATION_CO2_BASELEVEL
        !Config Desc  = CO2 base level
        !Config If    = OK_SECHIBA 
        !Config Def   = 280.
        !Config Help  = 
        !Config Units = [ppm]   
        CALL getin_p('DOWNREGULATION_CO2_BASELEVEL',downregulation_co2_baselevel)
        !-
        ! slowproc
        !-
        !
        !Config Key   = CLAYFRACTION_DEFAULT
        !Config Desc  = default fraction of clay
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.2 
        !Config Help  = 
        !Config Units = [-]   
        CALL getin_p('CLAYFRACTION_DEFAULT',clayfraction_default)
        !
        !Config Key   = MIN_VEGFRAC 
        !Config Desc  = Minimal fraction of mesh a vegetation type can occupy 
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.001 
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('MIN_VEGFRAC',min_vegfrac)
        !
        !Config Key   = STEMPDIAG_BID 
        !Config Desc  = only needed for an initial LAI if there is no restart file
        !Config If    = OK_SECHIBA 
        !Config Def   = 280.
        !Config Help  = 
        !Config Units = [K]
        CALL getin_p('STEMPDIAG_BID',stempdiag_bid)

!pss+
        !
        !Config Key   = SHIFT_FSAT_FWET
        !Config Desc  = shift saturation fraction to wetland fraction
        !Config If    = TOPM_calcul
        !Config Def   = 5.
        !Config Help  = 
        !Config Units = [-
        CALL getin_p('SHIFT_FSAT_FWET',SHIFT_fsat_fwet)
        !Config Key   = WTD1_borne 
        !Config Desc  = depth of subsurface saturation for wetland1 fraction
        !Config If    = TOPM_calcul 
        !Config Def   = 0.06
        !Config Help  = 
        !Config Units = [m]
        CALL getin_p('WTD1_BORNE',WTD1_borne)
        !Config Key   = WTD2_borne 
        !Config Desc  = depth of subsurface saturation for wetland1 fraction
        !Config If    = TOPM_calcul 
        !Config Def   = 0.12
        !Config Help  = 
        !Config Units = [m]
        CALL getin_p('WTD2_BORNE',WTD2_borne)
        !Config Key   = WTD3_borne 
        !Config Desc  = depth of subsurface saturation for wetland1 fraction
        !Config If    = TOPM_calcul 
        !Config Def   = 0.18
        !Config Help  = 
        !Config Units = [m]
        CALL getin_p('WTD3_BORNE',WTD3_borne)
        !Config Key   = WTD4_borne 
        !Config Desc  = depth of subsurface saturation for wetland1 fraction
        !Config If    = TOPM_calcul 
        !Config Def   = 0.24
        !Config Help  = 
        !Config Units = [m]
        CALL getin_p('WTD4_BORNE',WTD4_borne)
!pss-

        !
        first_call =.FALSE.
        
     ENDIF
     
   END SUBROUTINE config_sechiba_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_co2_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!! needed when OK_CO2 is set to true. (ie : when the photosynthesis is activated) 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

   SUBROUTINE config_co2_parameters
     
     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.4 Local variables
     
     LOGICAL, SAVE ::  first_call = .TRUE.      !! To keep first call trace (true/false)
!$OMP THREADPRIVATE(first_call)

!_ ================================================================================================================================
     
     IF(first_call) THEN
        
        !
        !Config Key   = LAI_LEVEL_DEPTH
        !Config Desc  = 
        !Config If    = OK_CO2
        !Config Def   = 0.15
        !Config Help  = 
        !Config Units = [-]  
        CALL getin_p('LAI_LEVEL_DEPTH',lai_level_depth)
        !
        !Config Key   = Oi
        !Config Desc  = Intercellular oxygen partial pressure
        !Config If    = OK_CO2
        !Config Def   = 210000.
        !Config Help  = See Legend of Figure 6 of Yin et al. (2009)
        !Config Units = [ubar]  
        CALL getin_p('Oi',Oi)
        
        first_call = .FALSE.
        
     ENDIF
     
   END SUBROUTINE config_co2_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_stomate_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!! needed when stomate is activated (ie : when OK_STOMATE is set to true).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

   SUBROUTINE config_stomate_parameters
     
    IMPLICIT NONE
    
    !! 0. Variables and parameters declaration

    !! 0.4 Local variables   

    LOGICAL, SAVE ::  first_call = .TRUE.  !! To keep first call trace (true/false)
!$OMP THREADPRIVATE(first_call)

!_ ================================================================================================================================
    
    IF(first_call) THEN
       !-
       ! constraints_parameters
       !-
       !
       !Config Key   = TOO_LONG 
       !Config Desc  = longest sustainable time without regeneration (vernalization)
       !Config If    = OK_STOMATE
       !Config Def   = 5.
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TOO_LONG',too_long)

       !-
       ! fire parameters
       !-
       !
       !Config Key   = TAU_FIRE 
       !Config Desc  = Time scale for memory of the fire index (days). Validated for one year in the DGVM. 
       !Config If    = OK_STOMATE 
       !Config Def   = 30.
       !Config Help  = 
       !Config Units = [days]    
       CALL getin_p('TAU_FIRE',tau_fire)
       !
       !Config Key   = LITTER_CRIT
       !Config Desc  = Critical litter quantity for fire
       !Config If    = OK_STOMATE 
       !Config Def   = 200.
       !Config Help  = 
       !Config Units = [gC/m^2]  
       CALL getin_p('LITTER_CRIT',litter_crit)
       !
       !Config Key   = FIRE_RESIST_STRUCT
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('FIRE_RESIST_STRUCT',fire_resist_struct)
       !
       !
       !Config Key   = CO2FRAC
       !Config Desc  = What fraction of a burned plant compartment goes into the atmosphere
       !Config If    = OK_STOMATE 
       !Config Def   = 0.95, 0.95, 0., 0.3, 0., 0., 0.95, 0.95
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('CO2FRAC',co2frac)
       !
       !Config Key   = BCFRAC_COEFF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3, 1.3, 88.2 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('BCFRAC_COEFF',bcfrac_coeff)
       !
       !Config Key   = FIREFRAC_COEFF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.45, 0.8, 0.6, 0.13
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FIREFRAC_COEFF',firefrac_coeff)

       !-
       ! gap parameters (+ lpj_const_mort)
       !-
       !
       !Config Key   = AVAILABILITY_FACT 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1
       !Config Help  = 
       !Config Units = [-]   
!       CALL getin_p('AVAILABILITY_FACT', availability_fact)  
       !
       !Config Key   = REF_GREFF
       !Config Desc  = Asymptotic maximum mortality rate
       !Config If    = OK_STOMATE 
       !Config Def   = 0.035
       !Config Help  = Set asymptotic maximum mortality rate from Sitch 2003
       !Config         (they use 0.01) (year^{-1})
       !Config Units = [1/year]  
       CALL getin_p('REF_GREFF',ref_greff)
       !-
       ! allocation parameters
       !-
       !
       !Config Key   = OK_MINRES
       !Config Desc  = Do we try to reach a minimum reservoir even if we are severely stressed?
       !Config If    = OK_STOMATE 
       !Config Def   = y
       !Config Help  = 
       !Config Units = [FLAG]
       CALL getin_p('OK_MINRES',ok_minres)
       !
       !Config Key   = RESERVE_TIME_TREE 
       !Config Desc  = maximum time during which reserve is used (trees) 
       !Config If    = OK_STOMATE 
       !Config Def   = 30.
       !Config Help  = 
       !Config Units = [days]    
       CALL getin_p('RESERVE_TIME_TREE',reserve_time_tree)
       !
!! Arsene 31-07-2014 modifications
       !Config Key   = RESERVE_TIME_SHRUB 
       !Config Desc  = maximum time during which reserve is used (shrub) 
       !Config If    = OK_STOMATE 
       !Config Def   = 25.
       !Config Help  = 
       !Config Units = [days]    
       CALL getin_p('RESERVE_TIME_SHRUB',reserve_time_shrub)
!! Arsene 31-07-2014 modifications
       !
       !Config Key   = RESERVE_TIME_GRASS 
       !Config Desc  = maximum time during which reserve is used (grasses) 
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('RESERVE_TIME_GRASS',reserve_time_grass)
       !
       !Config Key   = F_FRUIT
       !Config Desc  = Standard fruit allocation
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('F_FRUIT',f_fruit)
       !
       !Config Key   = ALLOC_SAP_ABOVE_GRASS 
       !Config Desc  = fraction of sapwood allocation above ground 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALLOC_SAP_ABOVE_GRASS',alloc_sap_above_grass)
       !
       !Config Key   = MIN_LTOLSR 
       !Config Desc  = extrema of leaf allocation fraction 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.2
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MIN_LTOLSR',min_LtoLSR)
       !
       !Config Key   = MAX_LTOLSR
       !Config Desc  = extrema of leaf allocation fraction
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MAX_LTOLSR',max_LtoLSR)
       !
       !Config Key   = Z_NITROGEN
       !Config Desc  = scaling depth for nitrogen limitation 
       !Config If    = OK_STOMATE
       !Config Def   = 0.2 
       !Config Help  =
       !Config Units = [m]  
       CALL getin_p('Z_NITROGEN',z_nitrogen)
       !
       !Config Key   = NLIM_TREF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 25. 
       !Config Help  = 
       !Config Units = [C]  
       CALL getin_p('NLIM_TREF',Nlim_tref)

!! Arsene 03-08-2015 - Add shrubs_like_trees
       !Config Key   = SHRUBS_LIKE_TREES
       !Config Desc  = Use tree sheme for shrubland (like geometry...)
       !Config If    = OK_STOMATE
       !Config Def   = FALSE
       !Config Help  = When to way are possible to modelise shrub, use the closest to tree
       !Config Units = [FLAG]
       CALL getin_p('SHRUBS_LIKE_TREES', shrubs_like_trees)
!! Arsene 03-08-2015 - Add shrubs_like_trees 
  
       !-
       ! data parameters
       !-
       !
       !Config Key   = PIPE_TUNE1
       !Config Desc  = crown area = pipe_tune1. stem diameter**(1.6) (Reinicke's theory)
       !Config If    = OK_STOMATE 
       !Config Def   = 100.0
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('PIPE_TUNE1',pipe_tune1)
       !
       !Config Key   = PIPE_TUNE2 
       !Config Desc  = height=pipe_tune2 * diameter**pipe_tune3
       !Config If    = OK_STOMATE 
       !Config Def   = 40.0 
       !Config Help  = 
       !Config Units = [-]      
       CALL getin_p('PIPE_TUNE2',pipe_tune2) 
        !
       !Config Key   = PIPE_TUNE3
       !Config Desc  = height=pipe_tune2 * diameter**pipe_tune3
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('PIPE_TUNE3',pipe_tune3)
       !
!! Arsene 11-08-2015 - New maxdia ==> pipe_tune4 not use
!!       !Config Key   = PIPE_TUNE4
!!       !Config Desc  = needed for stem diameter
!!       !Config If    = OK_STOMATE 
!!       !Config Def   = 0.3 
!!       !Config Help  = 
!!       !Config Units = [-]  
!!       CALL getin_p('PIPE_TUNE4',pipe_tune4)
!! Arsene 03-08-2015 - Note pipe_tune4 use only for max_dia... non use (for shrub ?) 
!! Arsene 11-08-2015 - New maxdia ==> pipe_tune4 not use
       !
       !Config Key   = PIPE_DENSITY 
       !Config Desc  = Density
       !Config If    = OK_STOMATE 
       !Config Def   = 2.e5 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('PIPE_DENSITY',pipe_density)
       !
       !Config Key   = PIPE_K1 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 8.e3 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('PIPE_K1',pipe_k1)
!! Arsene 03-08-2015 - Note: pipe_k1 use only to calcul first biomass... not need special for shrub ?
       !
       !Config Key   = PIPE_TUNE_EXP_COEFF 
       !Config Desc  = pipe tune exponential coeff 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.6 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('PIPE_TUNE_EXP_COEFF',pipe_tune_exp_coeff)
       !
       !
!! Arsene 03-08-2015 - New parametrisation for shrub... if like tree -START
       !Config Key   = PIPE_TUNE1_FOR_SHRUB
       !Config Desc  = crown area = pipe_tune1_for_shrub. stem diameter**(1.6) (Reinicke's theory)
       !Config If    = OK_STOMATE .AND. SHRUBS_LIKE_TREES
       !Config Def   = 217.
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('PIPE_TUNE1_FOR_SHRUB',pipe_tune1_for_shrub)
       !
       !Config Key   = PIPE_TUNE2_FOR_SHRUB
       !Config Desc  = height=pipe_tune2_for_shrub * diameter**pipe_tune3_for_shrub
       !Config If    = OK_STOMATE .AND. SHRUBS_LIKE_TREES
       !Config Def   = 8. 
       !Config Help  = 
       !Config Units = [-]      
       CALL getin_p('PIPE_TUNE2_FOR_SHRUB',pipe_tune2_for_shrub) 
       !
       !Config Key   = PIPE_TUNE3_FOR_SHRUB
       !Config Desc  = height=pipe_tune2_for_shrub * diameter**pipe_tune3_for_shrub
       !Config If    = OK_STOMATE .AND. SHRUBS_LIKE_TREES
       !Config Def   = 0.55 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('PIPE_TUNE3_FOR_SHRUB',pipe_tune3_for_shrub)
       !
       !Config Key   = PIPE_TUNE_EXP_COEFF_FOR_SHRUB 
       !Config Desc  = pipe tune exponential coeff for shrubs
       !Config If    = OK_STOMATE .AND. SHRUBS_LIKE_TREES
       !Config Def   = 2.8 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('PIPE_TUNE_EXP_COEFF_FOR_SHRUB',pipe_tune_exp_coeff_for_shrub)
!! Arsene 03-08-2015 - New parametrisation for shrub... if like tree -END
       !
       !
!! Arsene 03-08-2015 - New parametrisation for shrub...
       !Config Key   = PIPE_DENSITY_SHRUB 
       !Config Desc  = Density for shrubs... like tree for the moment 
       !Config If    = OK_STOMATE
       !Config Def   = 2.e5 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('PIPE_DENSITY_SHRUB',pipe_density_shrub)
       !
       !Config Key   = PIPE_K1_SHRUB
       !Config Desc  = for shrubs...
       !Config If    = OK_STOMATE 
       !Config Def   = 8.e3 ???
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('PIPE_K1_SHRUB',pipe_k1_shrub)
!! Arsene 03-08-2015 - New parametrisation for shrub...
       !
       !
!! Arsene 03-08-2015 - New parametrisation for shrub... if NOT like tree -START
       !Config Key   = PIPE_TUNE_SHRUB1
       !Config Desc  = total crown area = pipe_tune_shrub1. tatal basal area**(pipe_tune_shrub_exp_coeff) (Aiba & Kohyama (1996))
       !Config If    = OK_STOMATE .AND. .NOT.SHRUBS_LIKE_TREES
       !Config Def   = 10**2.42
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('PIPE_TUNE_SHRUB1',pipe_tune_shrub1)
       !
       !Config Key   = PIPE_TUNE_SHRUB2
       !Config Desc  = 1/height=1/(pipe_tune_shrub2 * (100*diameter)**pipe_tune_shrub3) + 1/height_presc
       !Config If    = OK_STOMATE .AND. .NOT.SHRUBS_LIKE_TREES
       !Config Def   = 0.75
       !Config Help  = 
       !Config Units = [-]      
       CALL getin_p('PIPE_TUNE_SHRUB2',pipe_tune_shrub2)
       !
       !Config Key   = PIPE_TUNE_SHRUB3
       !Config Desc  = 1/height=1/(pipe_tune_shrub2 * (100*diameter)**pipe_tune_shrub3) + 1/height_presc 
       !Config If    = OK_STOMATE .AND. .NOT.SHRUBS_LIKE_TREES
       !Config Def   = 1.15 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('PIPE_TUNE_SHRUB3',pipe_tune_shrub3)
       !
       !Config Key   = PIPE_TUNE_SHRUB_EXP_COEFF 
       !Config Desc  = pipe tune exponential coeff for shrubs (new equation)
       !Config If    = OK_STOMATE .AND. .NOT.SHRUBS_LIKE_TREES
       !Config Def   = 0.8 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('PIPE_TUNE_SHRUB_EXP_COEFF',pipe_tune_shrub_exp_coeff)
!! Arsene 03-08-2015 - New parametrisation for shrub... if like tree -END
       !
       !
!! Arsene 03-09-2015 - For iteration (new allomerty for shrubs)
       !Config Key   = SHRUB_IT_OK
       !Config Desc  = Active shrub iteration for allometry (with DGVM) - else Look-Up Table
       !Config If    = .NOT.SHRUBS_LIKE_TREES & DGVM_OK
       !Config Def   = .FALSE.
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('SHRUB_IT_OK',shrub_it_ok)
       !
       !Config Key   = ACCEPT_SIGMA_IT
       !Config Desc  = Define the precision wanted for the iteration result 
       !Config If    = .NOT.SHRUBS_LIKE_TREES
       !Config Def   = 0.01
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('ACCEPT_SIGMA_IT',accept_sigma_it)
       !
       !Config Key   = FACTOR_DIV_IT
       !Config Desc  = Define the factor division bewtween iteration 
       !Config If    = .NOT.SHRUBS_LIKE_TREES
       !Config Def   = 5.
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('FACTOR_DIV_IT',factor_div_it)
!! Arsene 03-09-2015 - For iteration (new allomerty for shrubs)
       !
!! Arsene 16-10-2015 - For shrub allometry array
       !Config Key   = SHRUB_ALLOM_LIG
       !Config Desc  = Define the number of line for shrub_allom_array
       !Config If    = .NOT.SHRUBS_LIKE_TREES
       !Config Def   = 1000
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('SHRUB_ALLOM_LIG', shrub_allom_lig)
       !
       !Config Key   = SHRUB_LIM_MAXDIA
       !Config Desc  = Define the proportion of real maxdia (= maxdia * shrub_lim_maxdia ) for shrub
       !Config If    = .NOT.SHRUBS_LIKE_TREES
       !Config Def   = 0.97
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('SHRUB_LIM_MAXDIA', shrub_lim_maxdia)
!! Arsene 16-10-2015 - For shrub allometry array
       !
       !
       !Config Key   = PRECIP_CRIT 
       !Config Desc  = minimum precip
       !Config If    = OK_STOMATE 
       !Config Def   = 100.
       !Config Help  = 
       !Config Units = [mm/year]  
       CALL getin_p('PRECIP_CRIT',precip_crit)
       !
       !Config Key   = GDD_CRIT_ESTAB
       !Config Desc  = minimum gdd for establishment of saplings
       !Config If    = OK_STOMATE 
       !Config Def   = 150. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('GDD_CRIT_ESTAB',gdd_crit_estab)
        !
       !Config Key   = FPC_CRIT
       !Config Desc  = critical fpc, needed for light competition and establishment
       !Config If    = OK_STOMATE 
       !Config Def   = 0.95
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('FPC_CRIT',fpc_crit)
       !
       !Config Key   = ALPHA_GRASS
       !Config Desc  = sapling characteristics : alpha's
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALPHA_GRASS',alpha_grass)
       !
       !Config Key   = ALPHA_TREE
       !Config Desc  = sapling characteristics : alpha's 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALPHA_TREE',alpha_tree)
       !
!! Arsene 31-07-2014 modifications
       !Config Key   = ALPHA_SHRUB
       !Config Desc  = sapling characteristics : alpha's 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.8
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALPHA_SHRUB',alpha_shrub)
!! Arsene 31-07-2014 modifications
       !-
       !
       !Config Key   = MASS_RATIO_HEART_SAP
       !Config Desc  = mass ratio (heartwood+sapwood)/sapwood
       !Config If    = OK_STOMATE 
       !Config Def   = 3.
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MASS_RATIO_HEART_SAP',mass_ratio_heart_sap)
       !
!!! Arsene 22-05-2015 modifications
!       !Config Key   = SHRUB_IND_FRAC
!       !Config Desc  = Ratio between number's of branches (ind) and real individual number 
!       !Config If    = OK_STOMATE 
!       !Config Def   = 2000.
!       !Config Help  = 
!       !Config Units = [-]   
!       CALL getin_p('SHRUB_IND_FRAC',shrub_ind_frac)
!!! Arsene 22-05-2015 modifications
       !-
       !Config Key   = TAU_HUM_MONTH
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('TAU_HUM_MONTH',tau_hum_month)
       !
       !Config Key   = TAU_HUM_WEEK
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 7.
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_HUM_WEEK',tau_hum_week)
       !
       !Config Key   = TAU_T2M_MONTH
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 20.
       !Config Help  = 
       !Config Units = [days]     
       CALL getin_p('TAU_T2M_MONTH',tau_t2m_month)
       !
       !Config Key   = TAU_T2M_WEEK
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 7.
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_T2M_WEEK',tau_t2m_week)
       !
       !Config Key   = TAU_TSOIL_MONTH 
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]     
       CALL getin_p('TAU_TSOIL_MONTH',tau_tsoil_month)
       !
       !Config Key   = TAU_SOILHUM_MONTH
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_SOILHUM_MONTH',tau_soilhum_month)
       !
       !Config Key   = TAU_GPP_WEEK 
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 7. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_GPP_WEEK',tau_gpp_week)
       !
       !Config Key   = TAU_GDD
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 40. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_GDD',tau_gdd)
       !
       !Config Key   = TAU_NGD
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 50.
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('TAU_NGD',tau_ngd)
       !
       !Config Key   = COEFF_TAU_LONGTERM
       !Config Desc  = time scales for phenology and other processes
       !Config If    = OK_STOMATE 
       !Config Def   = 3. 
       !Config Help  = 
       !Config Units = [days]   
       CALL getin_p('COEFF_TAU_LONGTERM',coeff_tau_longterm)
       !-
       !
       !Config Key   = BM_SAPL_CARBRES 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 5. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('BM_SAPL_CARBRES',bm_sapl_carbres)
       !
       !Config Key   = BM_SAPL_SAPABOVE
       !Config Desc  = 
       !Config If    = OK_STOMATE
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('BM_SAPL_SAPABOVE',bm_sapl_sapabove)
       !
       !Config Key   = BM_SAPL_HEARTABOVE 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 2.
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('BM_SAPL_HEARTABOVE',bm_sapl_heartabove)
       !
       !Config Key   = BM_SAPL_HEARTBELOW 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 2. 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('BM_SAPL_HEARTBELOW',bm_sapl_heartbelow)
       !
       !Config Key   = INIT_SAPL_MASS_LEAF_NAT
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_LEAF_NAT',init_sapl_mass_leaf_nat)
       !
       !Config Key   = INIT_SAPL_MASS_LEAF_AGRI
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 1. 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_LEAF_AGRI',init_sapl_mass_leaf_agri)
       !
!! Arsene 18-08-2015 - Add vor non vascular plant (important or bug with DGVM)
       !Config Key   = INIT_SAPL_MASS_LEAF_NOVASC
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.01 ?
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_LEAF_NOVASC', init_sapl_mass_leaf_novasc)
!! Arsene 18-08-2015 - Add for non vascular plant
       !
       !Config Key   = INIT_SAPL_MASS_CARBRES
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 5. 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_CARBRES',init_sapl_mass_carbres)
       !
       !Config Key   = INIT_SAPL_MASS_ROOT
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('INIT_SAPL_MASS_ROOT',init_sapl_mass_root)
       !
       !Config Key   = INIT_SAPL_MASS_FRUIT
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3 
       !Config Help  = 
       !Config Units = [-]    
       CALL getin_p('INIT_SAPL_MASS_FRUIT',init_sapl_mass_fruit)
       !
       !Config Key   = CN_SAPL_INIT 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('CN_SAPL_INIT',cn_sapl_init)
       !
       !Config Key   = MIGRATE_TREE 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 10000.
       !Config Help  = 
       !Config Units = [m/year]   
       CALL getin_p('MIGRATE_TREE',migrate_tree)
       !
!! Arsene 31-07-2014 modifications
       !Config Key   = MIGRATE_SHRUB
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 10000.
       !Config Help  = 
       !Config Units = [m/year]   
       CALL getin_p('MIGRATE_SHRUB',migrate_shrub)
!! Arsene 31-07-2014 modifications
       !
       !Config Key   = MIGRATE_GRASS
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 10000.
       !Config Help  = 
       !Config Units = [m/year]   
       CALL getin_p('MIGRATE_GRASS',migrate_grass)
       !
       !Config Key   = LAI_INITMIN_TREE
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3
       !Config Help  = 
       !Config Units = [m^2/m^2]  
       CALL getin_p('LAI_INITMIN_TREE',lai_initmin_tree)
       !
!! Arsene 31-07-2014 modifications
       !Config Key   = LAI_INITMIN_SHRUB
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.2
       !Config Help  = 
       !Config Units = [m^2/m^2]  
       CALL getin_p('LAI_INITMIN_SHRUB',lai_initmin_shrub)
!! Arsene 31-07-2014 modifications
       !
       !Config Key   = LAI_INITMIN_GRASS 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.1
       !Config Help  = 
       !Config Units = [m^2/m^2]    
       CALL getin_p('LAI_INITMIN_GRASS',lai_initmin_grass)
       !
       !Config Key   = DIA_COEFF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 4., 0.5
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('DIA_COEFF',dia_coeff)
       !
!! Arsene 11-08-2015 - New maxdia ==> maxdia_coeff not use
!!       !Config Key   = MAXDIA_COEFF
!!       !Config Desc  = 
!!       !Config If    = OK_STOMATE 
!!       !Config Def   = 100., 0.01 
!!       !Config Help  = 
!!       !Config Units = [-]   
!!       CALL getin_p('MAXDIA_COEFF',maxdia_coeff)
!! Arsene 11-08-2015 - New maxdia ==> maxdia_coeff not use
       !
!! Arsene 11-08-2015 - New calcul of bm_sapl ==> not use
!!       !Config Key   = BM_SAPL_LEAF
!!       !Config Desc  = 
!!       !Config If    = OK_STOMATE 
!!       !Config Def   = 4., 4., 0.8, 5. 
!!       !Config Help  = 
!!       !Config Units = [-]  
!!       CALL getin_p('BM_SAPL_LEAF',bm_sapl_leaf)
!! Arsene 11-08-2015 - New calcul of bm_sapl ==> not use

       !-
       ! litter parameters
       !-
       !
       !Config Key   = METABOLIC_REF_FRAC
       !Config Desc  =
       !Config If    = OK_STOMATE 
       !Config Def   = 0.85  
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('METABOLIC_REF_FRAC',metabolic_ref_frac)
       !
       !Config Key   = Z_DECOMP
       !Config Desc  = scaling depth for soil activity
       !Config If    = OK_STOMATE 
       !Config Def   = 0.2
       !Config Help  = 
       !Config Units = [m]   
       CALL getin_p('Z_DECOMP',z_decomp)
       !
       !Config Key   = CN
       !Config Desc  = C/N ratio
       !Config If    = OK_STOMATE 
       !Config Def   = 40., 40., 40., 40., 40., 40., 40., 40.
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('CN',CN)
       !
       !Config Key   = LC 
       !Config Desc  = Lignine/C ratio of the different plant parts
       !Config If    = OK_STOMATE 
       !Config Def   = 0.22, 0.35, 0.35, 0.35, 0.35, 0.22, 0.22, 0.22
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('LC',LC)
       !
       !Config Key   = FRAC_SOIL_STRUCT_AA
       !Config Desc  = frac_soil(istructural,iactive,iabove)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.55
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('FRAC_SOIL_STRUCT_AA',frac_soil_struct_aa)
       !
       !Config Key   = FRAC_SOIL_STRUCT_A 
       !Config Desc  = frac_soil(istructural,iactive,ibelow)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.45
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('FRAC_SOIL_STRUCT_AB',frac_soil_struct_ab)
       !
       !Config Key   = FRAC_SOIL_STRUCT_SA
       !Config Desc  = frac_soil(istructural,islow,iabove)
       !Config If    = OK_STOMATE
       !Config Def   = 0.7  
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FRAC_SOIL_STRUCT_SA',frac_soil_struct_sa)
       !
       !Config Key   = FRAC_SOIL_STRUCT_SB
       !Config Desc  = frac_soil(istructural,islow,ibelow) 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.7  
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FRAC_SOIL_STRUCT_SB',frac_soil_struct_sb)
       !
       !Config Key   = FRAC_SOIL_METAB_AA 
       !Config Desc  = frac_soil(imetabolic,iactive,iabove) 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.45 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FRAC_SOIL_METAB_AA',frac_soil_metab_aa)
       !
       !Config Key   = FRAC_SOIL_METAB_AB 
       !Config Desc  = frac_soil(imetabolic,iactive,ibelow)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.45  
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('FRAC_SOIL_METAB_AB',frac_soil_metab_ab)
       !
       !
       !Config Key   = METABOLIC_LN_RATIO
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.018  
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('METABOLIC_LN_RATIO',metabolic_LN_ratio) 
       !
       !Config Key   = TAU_METABOLIC
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.066
       !Config Help  = 
       !Config Units = [days] 
       CALL getin_p('TAU_METABOLIC',tau_metabolic)
       !
       !Config Key   = TAU_STRUCT 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.245 
       !Config Help  = 
       !Config Units = [days]
       CALL getin_p('TAU_STRUCT',tau_struct)
       !
       !Config Key   = SOIL_Q10
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.69 (=ln2)
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('SOIL_Q10',soil_Q10)
       !
       !Config Key   = TSOIL_REF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 30. 
       !Config Help  = 
       !Config Units = [C]   
       CALL getin_p('TSOIL_REF',tsoil_ref)
       !
       !Config Key   = LITTER_STRUCT_COEF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 3. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('LITTER_STRUCT_COEF',litter_struct_coef)
       !
       !Config Key   = MOIST_COEFF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.1, 2.4, 0.29
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MOIST_COEFF',moist_coeff)
       !
       !Config Key   = MOISTCONT_MIN
       !Config Desc  = minimum soil wetness to limit the heterotrophic respiration
       !Config If    = OK_STOMATE 
       !Config Def   = 0.25
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('MOISTCONT_MIN',moistcont_min)

       !-
       ! lpj parameters
       !-
       !
       !Config Key   = FRAC_TURNOVER_DAILY 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.55
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('FRAC_TURNOVER_DAILY',frac_turnover_daily)   
       !
!! Arsene 20-05-2015 fact_min_height
       !Config Key   = FACT_MIN_HEIGHT
       !Config Desc  = Coefficient: min_height = height_presc * fact_min_height
       !Config If    = OK_STOMATE
       !Config Def   = 10.
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('FACT_MIN_HEIGHT',fact_min_height)
!! Arsene 20-05-2015 fact_min_height


       !-
       ! npp parameters
       !-
       !
       !Config Key   = TAX_MAX
       !Config Desc  = maximum fraction of allocatable biomass used for maintenance respiration
       !Config If    = OK_STOMATE 
       !Config Def   = 0.8
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('TAX_MAX',tax_max) 

       !-
       ! phenology parameters
       !-
       !
       !Config Key   = ALWAYS_INIT
       !Config Desc  = take carbon from atmosphere if carbohydrate reserve too small? 
       !Config If    = OK_STOMATE 
       !Config Def   = n 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALWAYS_INIT',always_init)
       !
       !Config Key   = MIN_GROWTHINIT_TIME 
       !Config Desc  = minimum time since last beginning of a growing season
       !Config If    = OK_STOMATE 
       !Config Def   = 300. 
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('MIN_GROWTHINIT_TIME',min_growthinit_time)
       !
       !Config Key   = MOIAVAIL_ALWAYS_TREE
       !Config Desc  = moisture availability above which moisture tendency doesn't matter 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MOIAVAIL_ALWAYS_TREE',moiavail_always_tree)
       !
!! Arsene 31-07-2014 modifications
       !Config Key   = MOIAVAIL_ALWAYS_SHRUB
       !Config Desc  = moisture availability above which moisture tendency doesn't matter 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.9 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MOIAVAIL_ALWAYS_SHRUB',moiavail_always_shrub)
!! Arsene 31-07-2014 modifications
       !
       !Config Key   = MOIAVAIL_ALWAYS_GRASS 
       !Config Desc  = moisture availability above which moisture tendency doesn't matter
       !Config If    = OK_STOMATE 
       !Config Def   = 0.6 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('MOIAVAIL_ALWAYS_GRASS',moiavail_always_grass)
       !
       !Config Key   = T_ALWAYS_ADD
       !Config Desc  = monthly temp. above which temp. tendency doesn't matter 
       !Config If    = OK_STOMATE 
       !Config Def   = 10.
       !Config Help  = 
       !Config Units = [C]    
       CALL getin_p('T_ALWAYS_ADD',t_always_add)
       !
       !
       !Config Key   = GDDNCD_REF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 603. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('GDDNCD_REF',gddncd_ref)
       !
       !Config Key   = GDDNCD_CURVE
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.0091 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('GDDNCD_CURVE',gddncd_curve)
       !
       !Config Key   = GDDNCD_OFFSET
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 64. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('GDDNCD_OFFSET',gddncd_offset)
       !-
       ! prescribe parameters
       !-
       !
       !Config Key   = BM_SAPL_RESCALE 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 40. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('BM_SAPL_RESCALE',bm_sapl_rescale)

       !-
       ! respiration parameters
       !-
       !
       !Config Key   = MAINT_RESP_MIN_VMAX
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('MAINT_RESP_MIN_VMAX',maint_resp_min_vmax)  
       !
       !Config Key   = MAINT_RESP_COEFF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.4 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('MAINT_RESP_COEFF',maint_resp_coeff)

       !-
       ! soilcarbon parameters 
       !-
       !
       !Config Key   = FRAC_CARB_AP
       !Config Desc  = frac carb coefficients from active pool: depends on clay content
       !Config if    = OK_STOMATE 
       !Config Def   = 0.004
       !Config Help  = fraction of the active pool going into the passive pool
       !Config Units = [-]
       CALL getin_p('FRAC_CARB_AP',frac_carb_ap)  
       !
       !Config Key   = FRAC_CARB_SA
       !Config Desc  = frac_carb_coefficients from slow pool
       !Config if    = OK_STOMATE 
       !Config Def   = 0.42
       !Config Help  = fraction of the slow pool going into the active pool
       !Config Units = [-]
       CALL getin_p('FRAC_CARB_SA',frac_carb_sa)
       !
       !Config Key   = FRAC_CARB_SP
       !Config Desc  = frac_carb_coefficients from slow pool
       !Config if    = OK_STOMATE 
       !Config Def   = 0.03
       !Config Help  = fraction of the slow pool going into the passive pool
       !Config Units = [-] 
       CALL getin_p('FRAC_CARB_SP',frac_carb_sp)
       !
       !Config Key   = FRAC_CARB_PA
       !Config Desc  = frac_carb_coefficients from passive pool
       !Config if    = OK_STOMATE 
       !Config Def   = 0.45
       !Config Help  = fraction of the passive pool going into the active pool
       !Config Units = [-]
       CALL getin_p('FRAC_CARB_PA',frac_carb_pa)
       !
       !Config Key   = FRAC_CARB_PS
       !Config Desc  = frac_carb_coefficients from passive pool
       !Config if    = OK_STOMATE 
       !Config Def   = 0.0
       !Config Help  = fraction of the passive pool going into the slow pool
       !Config Units = [-]
       CALL getin_p('FRAC_CARB_PS',frac_carb_ps)
       !
       !Config Key   = ACTIVE_TO_PASS_CLAY_FRAC
       !Config Desc  = 
       !Config if    = OK_STOMATE 
       !Config Def   = 0.68  
       !Config Help  =
       !Config Units = [-] 
       CALL getin_p('ACTIVE_TO_PASS_CLAY_FRAC',active_to_pass_clay_frac)
       !
       !Config Key   = CARBON_TAU_IACTIVE
       !Config Desc  = residence times in carbon pools
       !Config if    = OK_STOMATE 
       !Config Def   = 0.149
       !Config Help  =
       !Config Units =  [days] 
       CALL getin_p('CARBON_TAU_IACTIVE',carbon_tau_iactive)
       !
       !Config Key   = CARBON_TAU_ISLOW
       !Config Desc  = residence times in carbon pools
       !Config if    = OK_STOMATE 
       !Config Def   = 5.48
       !Config Help  =
       !Config Units = [days]
       CALL getin_p('CARBON_TAU_ISLOW',carbon_tau_islow)
       !
       !Config Key   = CARBON_TAU_IPASSIVE
       !Config Desc  = residence times in carbon pools
       !Config if    = OK_STOMATE 
       !Config Def   = 241.
       !Config Help  = residence time in the passive pool
       !Config Units = [days] 
       CALL getin_p('CARBON_TAU_IPASSIVE',carbon_tau_ipassive)
       !
       !Config Key   = FLUX_TOT_COEFF
       !Config Desc  =
       !Config if    = OK_STOMATE 
       !Config Def   = 1.2, 1.4,.75
       !Config Help  =
       !Config Units = [days] 
       CALL getin_p('FLUX_TOT_COEFF',flux_tot_coeff)

       !-
       ! turnover parameters
       !-
       !
       !Config Key   = NEW_TURNOVER_TIME_REF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('NEW_TURNOVER_TIME_REF',new_turnover_time_ref)
       !
       !Config Key   = DT_TURNOVER_TIME 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 10.
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('DT_TURNOVER_TIME',dt_turnover_time)
       !
       !Config Key   = LEAF_AGE_CRIT_TREF
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 20. 
       !Config Help  = 
       !Config Units = [days]  
       CALL getin_p('LEAF_AGE_CRIT_TREF',leaf_age_crit_tref)
       !
       !Config Key   = LEAF_AGE_CRIT_COEFF 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.5, 0.75, 10. 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('LEAF_AGE_CRIT_COEFF',leaf_age_crit_coeff)
       !
!! Arsene 25-06-2014 NPPcumul parameters
       !Config Key   = NPP0_CUMUL_PARAM 
       !Config Desc  = 
       !Config If    = OK_STOMATE
       !Config Def   = 20., 70., 130., 0.01 
       !Config Help  = 
       !Config Units = [days,days,days,-] 
       CALL getin_p('NPP0_CUMUL_PARAM',npp0_c)
!! Arsene 25-06-2014 NPPcumul

!! Arsene 29-04-2015 LLAI_COEF for NVP turnover
       !Config Key   = LLAI_COEF 
       !Config Desc  = 
       !Config If    = OK_STOMATE
       !Config Def   = 0.06 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('LLAI_COEF',llai_coef)
!! Arsene 29-04-2015 LLAI_COEF for NVP turnover

       !-
       ! vmax parameters
       !-
       !
       !Config Key   = VMAX_OFFSET 
       !Config Desc  = offset (minimum relative vcmax)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.3
       !Config Help  = offset (minimum vcmax/vmax_opt)
       !Config Units = [-]  
       CALL getin_p('VMAX_OFFSET',vmax_offset)
       !
       !Config Key   = LEAFAGE_FIRSTMAX
       !Config Desc  = leaf age at which vmax attains vcmax_opt (in fraction of critical leaf age)
       !Config If    = OK_STOMATE 
       !Config Def   = 0.03 
       !Config Help  = relative leaf age at which vmax attains vcmax_opt
       !Config Units = [-] 
       CALL getin_p('LEAFAGE_FIRSTMAX',leafage_firstmax)
       !
       !Config Key   = LEAFAGE_LASTMAX 
       !Config Desc  = leaf age at which vmax falls below vcmax_opt (in fraction of critical leaf age) 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = relative leaf age at which vmax falls below vcmax_opt
       !Config Units = [-]  
       CALL getin_p('LEAFAGE_LASTMAX',leafage_lastmax)
       !
       !Config Key   = LEAFAGE_OLD 
       !Config Desc  = leaf age at which vmax attains its minimum (in fraction of critical leaf age)
       !Config If    = OK_STOMATE 
       !Config Def   = 1.
       !Config Help  = relative leaf age at which vmax attains its minimum
       !Config Units = [-]  
       CALL getin_p('LEAFAGE_OLD',leafage_old)
       !
!! Arsene 30-04-2015 Add vcmax_offset & humrel_mmin
       !Config Key   = VCMAX_OFFSET 
       !Config Desc  = offset of vcmax reduce by humrel_month (unitless) [0-1]
       !Config If    = OK_STOMATE
       !Config Def   = 0.3
       !Config Help  = impact of humrel_mont can't be below
       !Config Units = [-]  
       CALL getin_p('VCMAX_OFFSET',vcmax_offset)
       !
       !Config Key   = HUMREL_MMIN 
       !Config Desc  = Below this value, impact of humrel_month on NVP (unitless) [0-1]
       !Config If    = OK_STOMATE
       !Config Def   = 0.8
       !Config Help  = impact of humrel_month on vcmax below this value (for NVP)
       !Config Units = [-]  
       CALL getin_p('HUMREL_MMIN',humrel_mmin)
!! Arsene 30-04-2015 Add vcmax_offset & humrel_mmin
       !
!! Arsene 04-11-2015 Add vcmax_offset2 & humrel_mmax
       !Config Key   = VCMAX_OFFSET2 
       !Config Desc  = offset of vcmax reduce by humrel_month on grasses (unitless) [0-1]
       !Config If    = OK_STOMATE
       !Config Def   = 0.2
       !Config Help  = impact of humrel_month can't be below
       !Config Units = [-]  
       CALL getin_p('VCMAX_OFFSET2',vcmax_offset2)
       !
       !Config Key   = HUMREL_MMAX 
       !Config Desc  =  Above this value, impact of humrel_month on grasses (unitless) [0-1]
       !Config If    = OK_STOMATE
       !Config Def   = 0.97
       !Config Help  = impact of humrel_month on vcmax above this value (for grasses)
       !Config Units = [-]  
       CALL getin_p('HUMREL_MMAX',humrel_mmax)
!! Arsene 04-11-2015 Add vcmax_offset2 & humrel_mmax


       !-
       ! season parameters
       !-
       !
       !Config Key   = GPPFRAC_DORMANCE 
       !Config Desc  = rapport maximal GPP/GGP_max pour dormance
       !Config If    = OK_STOMATE 
       !Config Def   = 0.2 
       !Config Help  = 
       !Config Units = [-]
       CALL getin_p('GPPFRAC_DORMANCE',gppfrac_dormance)
       !
       !Config Key   = TAU_CLIMATOLOGY
       !Config Desc  = tau for "climatologic variables 
       !Config If    = OK_STOMATE 
       !Config Def   = 20
       !Config Help  = 
       !Config Units = [days]
       CALL getin_p('TAU_CLIMATOLOGY',tau_climatology)
       !
       !Config Key   = HVC1 
       !Config Desc  = parameters for herbivore activity
       !Config If    = OK_STOMATE 
       !Config Def   = 0.019
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('HVC1',hvc1)
       !
       !Config Key   = HVC2 
       !Config Desc  = parameters for herbivore activity 
       !Config If    = OK_STOMATE 
       !Config Def   = 1.38
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('HVC2',hvc2)
       !
       !Config Key   = LEAF_FRAC_HVC
       !Config Desc  = parameters for herbivore activity 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.33
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('LEAF_FRAC_HVC',leaf_frac_hvc)
       !
       !Config Key   = TLONG_REF_MAX
       !Config Desc  = maximum reference long term temperature 
       !Config If    = OK_STOMATE 
       !Config Def   = 303.1
       !Config Help  = 
       !Config Units = [K]  
       CALL getin_p('TLONG_REF_MAX',tlong_ref_max)
       !
       !Config Key   = TLONG_REF_MIN 
       !Config Desc  = minimum reference long term temperature 
       !Config If    = OK_STOMATE 
       !Config Def   = 253.1
       !Config Help  = 
       !Config Units = [K]  
       CALL getin_p('TLONG_REF_MIN',tlong_ref_min)
       !
       !Config Key   = NCD_MAX_YEAR
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 3. 
       !Config Help  = NCD : Number of Chilling Days
       !Config Units = [days]
       CALL getin_p('NCD_MAX_YEAR',ncd_max_year)
       !
       !Config Key   = GDD_THRESHOLD 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 5. 
       !Config Help  = GDD : Growing-Degree-Day
       !Config Units = [days] 
       CALL getin_p('GDD_THRESHOLD',gdd_threshold)
       !
       !Config Key   = GREEN_AGE_EVER 
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 2. 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('GREEN_AGE_EVER',green_age_ever)
       !
       !Config Key   = GREEN_AGE_DEC
       !Config Desc  = 
       !Config If    = OK_STOMATE 
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('GREEN_AGE_DEC',green_age_dec)

       !-
       ! establish WETLAND CH4 methane parameters
       !-
       !pss+

       !Config Key   = nvert
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 171 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NVERT',nvert)

       !Config Key   = ns
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 151 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NS',ns)

       !Config Key   = nday
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 24
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('NDAY',nday)

       !Config Key   = h
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.1
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('H',h)

       !Config Key   = rk
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 1
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RK',rk)

       !Config Key   = diffair
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 7.2 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('DIFFAIR',diffair)

       !Config Key   = pox
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('POX',pox)

       !Config Key   = dveg
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.001 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('DVEG',dveg)

       !Config Key   = rkm
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 5.0
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RKM',rkm)

       !Config Key   = xvmax
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 20.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('XVMAX',xvmax)

       !Config Key   = oxq10
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 2.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('OXQ10',oxq10)

       !Config Key   = scmax
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 500. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('SCMAX',scmax)

       !Config Key   = sr0pl
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 600. 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('SR0PL',sr0pl)

       !Config Key   = pwater_wet1
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -3 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_WET1',pwater_wet1)

       !Config Key   = pwater_wet2
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -9 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_WET2',pwater_wet2)

       !Config Key   = pwater_wet3
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -15 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_WET3',pwater_wet3)

       !Config Key   = pwater_wet4
       !Config Desc  = depth where saturation: definition for wetland 1  
       !Config If    = CH4_CALCUL
       !Config Def   = -21 
       !Config Help  = 
       !Config Units = [cm]   
       CALL getin_p('PWATER_WET4',pwater_wet4)

       !Config Key   = rpv
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 0.5 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RPV',rpv)

       !Config Key   = iother
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = -1.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('IOTHER',iother)

       !Config Key   = rq10
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = 3.0 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('RQ10',rq10)

       !Config Key   = alpha_CH4
       !Config Desc  = nb of vertical layers for CH4 diffusion 
       !Config If    = CH4_CALCUL
       !Config Def   = /0.009,0.004,0.021/
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ALPHA_CH4',alpha_CH4)
       !pss-

       
       first_call = .FALSE.
       
    ENDIF
    
  END SUBROUTINE config_stomate_parameters
!
!=
!

!! ================================================================================================================================
!! SUBROUTINE   : config_dgvm_parameters 
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters 
!! needed when the DGVM model is activated (ie : when ok_dgvm is set to true).
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_dgvm_parameters   
    
    IMPLICIT NONE
    
    !! 0. Variables and parameters declaration

    !! 0.4 Local variables

    LOGICAL, SAVE ::  first_call = .TRUE.         !! To keep first call trace (true/false)
!$OMP THREADPRIVATE(first_call)

!_ ================================================================================================================================    

    IF(first_call) THEN
  
       !-
       ! establish parameters
       !-
       !
       !Config Key   = ESTAB_MAX_TREE
       !Config Desc  = Maximum tree establishment rate 
       !Config If    = OK_DGVM
       !Config Def   = 0.12 
       !Config Help  = 
       !Config Units = [-]   
       CALL getin_p('ESTAB_MAX_TREE',estab_max_tree)
       !
       !Config Key   = ESTAB_MAX_GRASS
       !Config Desc  = Maximum grass establishment rate
       !Config If    = OK_DGVM
       !Config Def   = 0.12 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('ESTAB_MAX_GRASS',estab_max_grass)
       !
       !Config Key   = ESTABLISH_SCAL_FACT
       !Config Desc  = 
       !Config If    = OK_DGVM 
       !Config Def   = 5.
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('ESTABLISH_SCAL_FACT',establish_scal_fact)
       !
       !Config Key   = MAX_TREE_COVERAGE 
       !Config Desc  = 
       !Config If    = OK_DGVM 
       !Config Def   = 0.98
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('MAX_TREE_COVERAGE',max_tree_coverage)
       !
       !Config Key   = IND_0_ESTAB
       !Config Desc  = 
       !Config If    = OK_DGVM 
       !Config Def   = 0.2
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('IND_0_ESTAB',ind_0_estab)

       !-
       ! light parameters
       !-
       !
       !Config Key   = ANNUAL_INCREASE
       !Config Desc  = for diagnosis of fpc increase, compare today's fpc to last year's maximum (T) or to fpc of last time step (F)?
       !Config If    = OK_DGVM
       !Config Def   = y
       !Config Help  = 
       !Config Units = [FLAG]
       CALL getin_p('ANNUAL_INCREASE',annual_increase)
       !
       !Config Key   = MIN_COVER 
       !Config Desc  = For trees, minimum fraction of crown area occupied 
       !Config If    = OK_DGVM
       !Config Def   = 0.05 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('MIN_COVER',min_cover)

       !-
       ! pftinout parameters
       !
       !Config Key   = IND_0 
       !Config Desc  = initial density of individuals
       !Config If    = OK_DGVM
       !Config Def   = 0.02 
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('IND_0',ind_0)
       !
       !Config Key   = MIN_AVAIL
       !Config Desc  = minimum availability
       !Config If    = OK_DGVM
       !Config Def   = 0.01
       !Config Help  = 
       !Config Units = [-]  
       CALL getin_p('MIN_AVAIL',min_avail)
       !
       !Config Key   = RIP_TIME_MIN
       !Config Desc  = 
       !Config If    = OK_DGVM
       !Config Def   = 1.25 
       !Config Help  = 
       !Config Units = [year]  
       CALL getin_p('RIP_TIME_MIN',RIP_time_min)
       !
       !Config Key   = NPP_LONGTERM_INIT
       !Config Desc  = 
       !Config If    = OK_DGVM
       !Config Def   = 10.
       !Config Help  = 
       !Config Units = [gC/m^2/year]
       CALL getin_p('NPP_LONGTERM_INIT',npp_longterm_init)
       !
       !Config Key   = EVERYWHERE_INIT
       !Config Desc  = 
       !Config If    = OK_DGVM
       !Config Def   = 0.05 
       !Config Help  = 
       !Config Units = [-] 
       CALL getin_p('EVERYWHERE_INIT',everywhere_init)
       
       first_call = .FALSE.
       
    ENDIF
    
    
  END SUBROUTINE config_dgvm_parameters


END MODULE constantes
