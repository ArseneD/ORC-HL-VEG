! =================================================================================================================================
! MODULE 	: constantes_soil
!
! CONTACT       : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "constantes_soil" module contains subroutine to initialize the parameters related to soil and hydrology.
!!
!!\n DESCRIPTION : "constantes_soil" module contains subroutine to initialize the parameters related to soil and hydrology.
!!                 This module alos USE constates_soil and can therfor be used to acces the subroutines and the constantes.
!!                 The constantes declarations can also be used seperatly with "USE constantes_soil_var".
!!
!! RECENT CHANGE(S): 
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================

MODULE constantes_soil

  USE constantes_soil_var
  USE constantes
  USE ioipsl_para 

  IMPLICIT NONE

CONTAINS

  
!! ================================================================================================================================
!! SUBROUTINE   : config_soil_parameters
!!
!>\BRIEF        This subroutine reads in the configuration file all the parameters related to soil and hydrology. 
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

   SUBROUTINE config_soil_parameters(active_flags, impose_param)
     
     USE ioipsl

     IMPLICIT NONE

     !! 0. Variables and parameters declaration

     !! 0.1 Input variables

     TYPE(control_type), INTENT(in) :: active_flags            !! What parts of the code are activated ? (true/false)  
     LOGICAL, INTENT(IN)            :: impose_param            !! Flag for imposing parameters in run.def
    !! 0.4 Local variables 

     INTEGER(i_std), PARAMETER      :: error_level = 3         !! Switch to 2 to turn fatal errors into warnings.(1-3, unitless)
     LOGICAL, SAVE                  :: first_call = .TRUE.     !! To keep first call trace (true/false)
!$OMP THREADPRIVATE(first_call)
     LOGICAL                        :: ok_freeze               !! Local variable used to set default values for all flags 
                                                               !! controling the soil freezing scheme
!_ ================================================================================================================================
     
     IF ( first_call ) THEN

        !Config Key   = THERMOSOIL_NBLEV
        !Config Desc  = Number of soil level
        !Config If    = 
        !Config Def   = 7
        !Config Help  = Use at least 11 for long term simulation where soil thermal inertia matters
        !Config Units = (-)
        ngrnd=7
        CALL getin_p("THERMOSOIL_NBLEV",ngrnd)
        

        ! Following initializations are only done for option impose_param
        IF ( active_flags%ok_sechiba .AND. impose_param ) THEN

        !Config Key   = DRY_SOIL_HEAT_CAPACITY
        !Config Desc  = Dry soil Heat capacity of soils
        !Config If    = OK_SECHIBA 
        !Config Def   = 1.80e+6
        !Config Help  = Values taken from : PIELKE,'MESOSCALE METEOROLOGICAL MODELING',P.384.
        !Config Units = [J.m^{-3}.K^{-1}] 
        CALL getin_p("DRY_SOIL_HEAT_CAPACITY",so_capa_dry)
        
        !! Check parameter value (correct range)
        IF ( so_capa_dry <= zero ) THEN
           CALL ipslerr_p(error_level, "config_soil_parameters.", &
                &     "Wrong parameter value for DRY_SOIL_HEAT_CAPACITY.", &
                &     "This parameter should be positive. ", &
                &     "Please, check parameter value in run.def. ")
        END IF
        

        !Config Key   = DRY_SOIL_HEAT_COND
        !Config Desc  = Dry soil Thermal Conductivity of soils
        !Config If    = OK_SECHIBA
        !Config Def   = 0.40 
        !Config Help  = Values taken from : PIELKE,'MESOSCALE METEOROLOGICAL MODELING',P.384.
        !Config Units = [W.m^{-2}.K^{-1}] 
        CALL getin_p("DRY_SOIL_HEAT_COND",so_cond_dry)

        !! Check parameter value (correct range)
        IF ( so_cond_dry <= zero ) THEN
           CALL ipslerr_p(error_level, "config_soil_parameters.", &
                &     "Wrong parameter value for DRY_SOIL_HEAT_COND.", &
                &     "This parameter should be positive. ", &
                &     "Please, check parameter value in run.def. ")
        END IF


        !Config Key   = WET_SOIL_HEAT_CAPACITY
        !Config Desc  = Wet soil Heat capacity of soils 
        !Config If    = OK_SECHIBA
        !Config Def   = 3.03e+6
        !Config Help  = 
        !Config Units = [J.m^{-3}.K^{-1}]
        CALL getin_p("WET_SOIL_HEAT_CAPACITY",so_capa_wet)

        !! Check parameter value (correct range)
        IF ( so_capa_wet <= zero ) THEN
           CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for WET_SOIL_HEAT_CAPACITY.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
        END IF


        !Config Key   = WET_SOIL_HEAT_COND
        !Config Desc  = Wet soil Thermal Conductivity of soils
        !Config If    = OK_SECHIBA 
        !Config Def   = 1.89 
        !Config Help  = 
        !Config Units = [W.m^{-2}.K^{-1}]
        CALL getin_p("WET_SOIL_HEAT_COND",so_cond_wet)

        !! Check parameter value (correct range)
        IF ( so_cond_wet <= zero ) THEN
           CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for WET_SOIL_HEAT_COND.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
        END IF


        !Config Key   = SNOW_HEAT_COND
        !Config Desc  = Thermal Conductivity of snow
        !Config If    = OK_SECHIBA  
        !Config Def   = 0.3
        !Config Help  = 
        !Config Units = [W.m^{-2}.K^{-1}]
        CALL getin_p("SNOW_HEAT_COND",sn_cond)

        !! Check
        IF ( sn_cond <= zero ) THEN
           CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for SNOW_HEAT_COND.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
        END IF


        !Config Key   = SNOW_DENSITY
        !Config Desc  = Snow density for the soil thermodynamics 
        !Config If    = OK_SECHIBA 
        !Config Def   = 330.0
        !Config Help  = 
        !Config Units = [-] 
        CALL getin_p("SNOW_DENSITY",sn_dens)
        
        !! Check parameter value (correct range)
        IF ( sn_dens <= zero ) THEN
          CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for SNOW_DENSITY.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
        END IF


        !! Calculation of snow capacity
        !! If sn_dens is redefined by the user, sn_capa needs to be reset
        sn_capa = 2100.0_r_std*sn_dens


        !Config Key   = NOBIO_WATER_CAPAC_VOLUMETRI
        !Config Desc  = 
        !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
        !Config Def   = 150.
        !Config Help  = 
        !Config Units = [s/m^2]
        CALL getin_p('NOBIO_WATER_CAPAC_VOLUMETRI',mx_eau_nobio)

       !! Check parameter value (correct range)
        IF ( mx_eau_nobio <= zero ) THEN
           CALL ipslerr_p(error_level, "config_soil_parameters.", &
               &     "Wrong parameter value for NOBIO_WATER_CAPAC_VOLUMETRI.", &
               &     "This parameter should be positive. ", &
               &     "Please, check parameter value in run.def. ")
        END IF


        !Config Key   = SECHIBA_QSINT 
        !Config Desc  = Interception reservoir coefficient
        !Config If    = OK_SECHIBA 
        !Config Def   = 0.1
        !Config Help  = Transforms leaf area index into size of interception reservoir
        !Config         for slowproc_derivvar or stomate
        !Config Units = [m]
        CALL getin_p('SECHIBA_QSINT',qsintcst)

        !! Check parameter value (correct range)
        IF ( qsintcst <= zero ) THEN
           CALL ipslerr_p(error_level, "config_soil_parameters.", &
                &     "Wrong parameter value for SECHIBA_QSINT.", &
                &     "This parameter should be positive. ", &
                &     "Please, check parameter value in run.def. ")
        END IF


        IF ( .NOT.(active_flags%hydrol_cwrr) ) THEN
           
           !Config Key   = CHOISNEL_DIFF_MIN
           !Config Desc  = Diffusion constant for the slow regime
           !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
           !Config Def   = 0.001
           !Config Help  = 
           !Config Units = [kg/m^2/dt]
           CALL getin_p('CHOISNEL_DIFF_MIN',min_drain)

           !! Check parameter value (correct range)
           IF ( min_drain <= zero ) THEN
              CALL ipslerr_p(error_level, "config_soil_parameters.", &
                   &     "Wrong parameter value for CHOISNEL_DIFF_MIN.", &
                   &     "This parameter should be positive. ", &
                   &     "Please, check parameter value in run.def. ")
            END IF


           !Config Key   = CHOISNEL_DIFF_MAX
           !Config Desc  = Diffusion constant for the fast regime
           !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
           !Config Def   = 0.1
           !Config Help  = 
           !Config Units = [kg/m^2/dt]
           CALL getin_p('CHOISNEL_DIFF_MAX',max_drain)

           !! Check parameter value (correct range)
           IF (  ( max_drain <= zero ) .OR. ( max_drain <= min_drain ) ) THEN
              CALL ipslerr_p(error_level, "config_soil_parameters.", &
                   &     "Wrong parameter value for CHOISNEL_DIFF_MAX.", &
                   &     "This parameter should be positive or greater than CHOISNEL_DIFF_MIN.", &
                   &     "Please, check parameter value in run.def. ")
           END IF


           !Config Key   = CHOISNEL_DIFF_EXP
           !Config Desc  = The exponential in the diffusion law
           !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
           !Config Def   = 1.5
           !Config Help  = 
           !Config Units = [-]
           CALL getin_p('CHOISNEL_DIFF_EXP',exp_drain)
           
           !! Check parameter value (correct range)
           IF ( exp_drain <= zero ) THEN
              CALL ipslerr_p(error_level, "config_soil_parameters.", &
                   &     "Wrong parameter value for CHOISNEL_DIFF_EXP.", &
                   &     "This parameter should be positive. ", &
                   &     "Please, check parameter value in run.def. ")
           END IF


           !Config Key   = CHOISNEL_RSOL_CSTE
           !Config Desc  = Constant in the computation of resistance for bare  soil evaporation 
           !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR)
           !Config Def   = 33.E3
           !Config Help  = 
           !Config Units = [s/m^2]
           CALL getin_p('CHOISNEL_RSOL_CSTE',rsol_cste)

           !! Check parameter value (correct range)
           IF ( rsol_cste <= zero ) THEN
              CALL ipslerr_p(error_level, "config_soil_parameters.", &
                   &     "Wrong parameter value for CHOISNEL_RSOL_CSTE.", &
                   &     "This parameter should be positive. ", &
                   &     "Please, check parameter value in run.def. ")
           END IF


           !Config Key   = HCRIT_LITTER
           !Config Desc  = Scaling depth for litter humidity
           !Config If    = OK_SECHIBA and .NOT.(HYDROL_CWRR) 
           !Config Def   = 0.08 
           !Config Help  = 
           !Config Units = [m]
           CALL getin_p('HCRIT_LITTER',hcrit_litter)

           !! Check parameter value (correct range)
           IF ( hcrit_litter <= zero ) THEN
              CALL ipslerr_p(error_level, "config_soil_parameters.", &
                   &     "Wrong parameter value for HCRIT_LITTER.", &
                   &     "This parameter should be positive. ", &
                   &     "Please, check parameter value in run.def. ")
           END IF

        END IF
     
        END IF ! IF ( active_flags%ok_sechiba .AND. impose_param ) THEN



        !! Variables related to permafrost carbon

        !! Variables related to soil freezing in thermosoil module
        !
        !Config Key  = OK_FREEZE
        !Config Desc = Activate the complet soil freezing scheme
        !Config If   = OK_SECHIBA 
        !Config Def  = FALSE
        !Config Help = Activate soil freezing thermal effects. Activates soil freezing hydrological effects in CWRR scheme.
        !Config Units= [FLAG]

        ! ok_freeze is a flag that controls the default values for several flags controling 
        ! the different soil freezing processes
        ! Set ok_freeze=true for the complete soil freezing scheme
        ! ok_freeze is a local variable only used in this subroutine
        ok_freeze = .FALSE.
        CALL getin_p('OK_FREEZE',ok_freeze)


        !Config Key  = OK_ECORR
        !Config Desc = Energy correction for freezing
        !Config If   = 
        !Config Def  = True if OK_FREEZE else false
        !Config Help = Energy conservation : Correction to make sure that the same latent heat is 
        !Config        released and consumed during freezing and thawing
        !Config Units= [FLAG]

        IF (ok_freeze) THEN
           ok_Ecorr = .TRUE.
        ELSE
           ok_Ecorr = .FALSE.
        END IF
        CALL getin_p ('OK_ECORR',ok_Ecorr)

        !Config Key  = READ_PERMAFROST_MAP
        !Config Desc = Read information about ice content, overburden and permafrost type from IPA map
        !Config If   = 
        !Config Def  = FALSE
        !Config Help = 
        !Config Units= [FLAG]

        read_permafrost_map = .FALSE.
        CALL getin_p ('READ_PERMAFROST_MAP',read_permafrost_map)       

        !Config Key  = READ_REFTEMP
        !Config Desc = Initialize soil temperature using climatological temperature
        !Config If   = 
        !Config Def  = True if OK_FREEZE else false
        !Config Help = 
        !Config Units= [FLAG]

        IF (ok_freeze) THEN
           read_reftemp = .FALSE.
        ELSE
           read_reftemp = .FALSE.
        END IF
        CALL getin_p ('READ_REFTEMP',read_reftemp)

        !Config Key  = OK_FREEZE_THERMIX
        !Config Desc = Activate thermal part of the soil freezing scheme
        !Config If   = 
        !Config Def  = True if OK_FREEZE else false
        !Config Help = 
        !Config Units= [FLAG]

        IF (ok_freeze) THEN
           ok_freeze_thermix = .TRUE.
        ELSE
           ok_freeze_thermix = .FALSE.
        END IF
        CALL getin_p ('OK_FREEZE_THERMIX',ok_freeze_thermix)

        !! Coherence check for number of thermosoil levels for long term simulation where soil thermal inertia matters
        !! It is highly recommnaded to use at least ngrnd=11 when soil freezing is activated
        IF (ok_freeze_thermix .AND. ngrnd < 11) THEN
           WRITE(numout,*) 'ERROR : Incoherence between ok_freeze_thermix activated and ngrnd to small. Here used ngrnd=',ngrnd
           WRITE(numout,*) 'Set THERMOSOIL_NBLEV=11 or higher in run.def parameter file or deactivate soil freezing'
           CALL ipslerr_p(3,'thermosoil_init','Not enough thermodynamic soil levels for soil freezing', &
                'Adapt run.def with at least THERMOSOIL_NBLEV=11','')
        END IF


        !Config Key = POROS
        !Config Desc = Soil porosity 
        !Config If = OK_SECHIBA
        !Config Def = 0.41
        !Config Help = From USDA classification, mean value
        !Config Units = [-] 
        poros=0.41
        CALL getin_p('POROS',poros)


        !Config Key = fr_dT
        !Config Desc = Freezing window    
        !Config If = OK_SECHIBA
        !Config Def = 2.0
        !Config Help = 
        !Config Units = [K] 
        fr_dT=2.0
        CALL getin_p('FR_DT',fr_dT)
              

        !! Variables related to soil Freezing in diffuco module

        !Config Key  = OK_SNOWFACT
        !Config Desc = Activates the smoothering of landscapes by snow,
        !       e.g. reduces of the surface roughness length when snow is present.
        !Config If   = 
        !Config Def  = True if OK_FREEZE else false
        !Config Help = 
        !Config Units= [FLAG]

        IF (ok_freeze) THEN
           ok_snowfact = .TRUE.
        ELSE
           ok_snowfact = .FALSE.
        END IF
        CALL getin_p('OK_SNOWFACT', ok_snowfact)


        !! Variables related to soil Freezing in hydrol module

        !Config Key  = OK_FREEZE_CWRR
        !Config Desc = CWRR freezing scheme by I. Gouttevin
        !Config If   = 
        !Config Def  = True if OK_FREEZE else false
        !Config Help =
        !Config Units= [FLAG]

        IF (ok_freeze) THEN
           ok_freeze_cwrr = .TRUE.
        ELSE
           ok_freeze_cwrr = .FALSE.
        END IF
        CALL getin_p('OK_FREEZE_CWRR',ok_freeze_cwrr)


        IF (ok_freeze_cwrr) THEN
           !Config Key  = OK_THERMODYNAMICAL_FREEZING
           !Config Desc = Calculate frozen fraction thermodynamically 
           !Config If   = HYDROL_CWRR .AND. OK_FREEZE_CWRR
           !Config Def  = True
           !Config Help = Calculate frozen fraction thermodynamically if true,
           !Config      = else calculate frozen fraction linearly 
           !Config Units= [FLAG]
           ok_thermodynamical_freezing = .TRUE.
           CALL getin_p('OK_THERMODYNAMICAL_FREEZING',ok_thermodynamical_freezing)
        END IF
     
        first_call =.FALSE.
        
     ENDIF
     
   END SUBROUTINE config_soil_parameters
   

END MODULE constantes_soil
