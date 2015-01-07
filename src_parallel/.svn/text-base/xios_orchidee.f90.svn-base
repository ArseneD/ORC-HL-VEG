! ================================================================================================================================
!  MODULE       : xios_orchidee
!
!  CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   This module contains the initialization and interface to the XIOS code.
!!
!!\n DESCRIPTION: This module contains the interface for the use of the XIOS code. All call to XIOS is done in this module.
!!
!!                Summury of subroutines
!!                      xios_orchidee_comm_init       : First call to XIOS to get the MPI communicator 
!!                      xios_orchidee_init            : Initialize variables needed for use of XIOS 
!!                                                      Deactivation of fields not calculated due specific run options
!!                      xios_orchidee_update_calendar : Update the calandar in XIOS
!!                      xios_orchidee_finalize        : Last call to XIOS for finalization
!!                      xios_orchidee_send_field      : Interface to send fields with 1, 2 or 3 dimensions to XIOS
!!                      xios_orchidee_send_field_r1d  : Internal subroutine for 1D(array) fields
!!                      xios_orchidee_send_field_r2d  : Internal subroutine for 2D fields
!!                      xios_orchidee_send_field_r3d  : Internal subroutine for 3D fields
!!
!!                Note that compilation must be done with the preporcessing key XIOS and CPP_PARA. Compiling without these 
!!                keys makes it impossible to activate XIOS. To activate run using XIOS, the flag XIOS_ORCHIDEE_OK=y must 
!!                be set in run.def and the file iodef.xml must exist.  
!!
!! RECENT CHANGE(S): Created by Arnaud Caubel(LSCE), Josefine Ghattas (IPSL) 2013
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================

MODULE xios_orchidee

#ifdef XIOS
  USE xios
#endif
  USE defprec
  USE pft_parameters_var, ONLY : nvm
  USE constantes_var, ONLY : nlai, nnobio, nsnow, control_type, long_print, almaoutput, dt_sechiba
  USE constantes_soil_var, ONLY : ngrnd, nstm, nslm, nbdl, check_waterbal, diaglev
  USE IOIPSL, ONLY : ioget_calendar, ju2ymds, getin
  USE mod_orchidee_para_var
  USE mod_orchidee_transfert_para
  USE ioipsl_para

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xios_orchidee_comm_init, xios_orchidee_init, xios_orchidee_change_context, &
            xios_orchidee_update_calendar, xios_orchidee_context_finalize, xios_orchidee_finalize, &
            xios_orchidee_send_field

  !
  !! Declaration of public variables
  !
  LOGICAL, PUBLIC, SAVE           :: xios_orchidee_ok=.FALSE.     !! Use XIOS for diagnostic files
  !$OMP THREADPRIVATE(xios_orchidee_ok)

  !
  !! Declaration of internal variables
  !
#ifdef XIOS
  TYPE(xios_context)              :: ctx_hdl_orchidee      !! Handel for ORCHIDEE
  !$OMP THREADPRIVATE(ctx_hdl_orchidee)
#endif
  CHARACTER(len=*),PARAMETER      :: id="client"           !! Id for initialization of ORCHIDEE in XIOS



  !! ==============================================================================================================================
  !! INTERFACE   : xios_orchidee_send_field
  !!
  !>\BRIEF         Send a field to XIOS.
  !!
  !! DESCRIPTION  :\n Send a field to XIOS. The field can have 1, 2 or 3 dimensions.
  !!                  This interface should be called at each time-step for each output varaiables.
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE xios_orchidee_send_field
     MODULE PROCEDURE xios_orchidee_send_field_r1d, xios_orchidee_send_field_r2d, xios_orchidee_send_field_r3d
  END INTERFACE


CONTAINS
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_comm_init 
  !!
  !>\BRIEF         Get the MPI communicator.
  !!
  !! DESCRIPTION  :\n First call to XIOS to get the MPI communicator. 
  !!                  Note that it is XIOS that initialize the MPI communicator.
  !!                  This subroutine is only called in ORCHIDEE offline mode. When running in coupled mode, the 
  !!                  atmospheric model must initlialize XIOS at the same time as initializing MPI. 
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_comm_init(comm_local)
    !
    !! 0. Variable and parameter declaration
    !
    !!    Output variables
    INTEGER, INTENT(OUT) :: comm_local

    !_ ================================================================================================================================

    IF (is_omp_root) THEN
#ifdef XIOS
       CALL xios_initialize(id,return_comm=comm_local)
#else
       ! Write error messages and stop the model
       WRITE(numout,*) 'Preprocessing key XIOS is missing to run ORCHIDEE with XIOS'
       WRITE(numout,*) 'Recompile with preprocessing flag XIOS or set XIOS_ORCHIDEE_OK=n in run.def'
       WRITE(numout,*) 'Fatal error from ORCHIDEE. STOP in xios_orchidee_comm_init'
#ifdef CPP_PARA
       CALL MPI_ABORT(3)
#endif     
       
#endif
    
    END IF
  END SUBROUTINE xios_orchidee_comm_init


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_init 
  !!
  !>\BRIEF         Initialize variables needed for use of XIOS.
  !!
  !! DESCRIPTION  :\n Initialization of specific varaiables needed to use XIOS such as model domain and time step. 
  !!
  !!                  In this subroutine also a section containg deactivation of some fields is found. The variables are 
  !!                  deactivated of not according to the corresponding control flag. For exemple the variables cacluated by the 
  !!                  routing scheme will be deactivated if the routing is deactivated. This is done to be able to keep the same 
  !!                  iodef.xml input file for several options without geting empty fields in the output file. Note that a field that
  !!                  is activated in the code can always be deactivated from the iodef.xml external file. 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_init(MPI_COMM_ORCH,               &
       date0,    year,      month,             day,          &
       lon_mpi,  lat_mpi,   soilth_lev,        control_flags )

    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                            :: MPI_COMM_ORCH    !! Orchidee MPI communicator (from module mod_orchidee_mpi_data)
    REAL(r_std), INTENT(in)                               :: date0            !! Julian day at first time step
    INTEGER(i_std), INTENT(in)                            :: year, month, day !! Current date information
    REAL(r_std),DIMENSION (iim_g,jj_nb), INTENT(in)       :: lon_mpi, lat_mpi !! Longitudes and latitudes on MPI local domain 2D domain
    REAL(r_std),DIMENSION (ngrnd), INTENT(in)             :: soilth_lev       !! Vertical soil levels for thermal scheme (m)
    TYPE(control_type)                                    :: control_flags    !! Flags that activate parts of the model

    !
    !! 0.2 Local variables
    !
#ifdef XIOS
    TYPE(xios_time)                :: dtime_xios
    TYPE(xios_fieldgroup)          :: fieldgroup_handle
    TYPE(xios_field)               :: field_handle
    TYPE(xios_file)                :: file_handle
#endif
    INTEGER(i_std)                 :: i
    INTEGER(i_std)                 :: year0, month0, day0 !! Time origin date information
    REAL(r_std)                    :: sec0                !! Time origin date information
    CHARACTER(LEN=20)              :: calendar_str        !! Name of current calendar
    CHARACTER(LEN=30)              :: start_str           !! Current date as character string
    CHARACTER(LEN=30)              :: startorig_str       !! Time origin date as character string
    !_ ================================================================================================================================
    
    
    IF (long_print) WRITE(numout,*) 'Entering xios_orchidee_init'

    !Config Key   = XIOS_ORCHIDEE_OK
    !Config Desc  = Use XIOS for writing diagnostics file
    !Config If    = 
    !Config Def   = n
    !Config Help  = 
    !Config Units = [FLAG]
    CALL getin_p('XIOS_ORCHIDEE_OK',xios_orchidee_ok)
    WRITE(numout,*)'In xios_orchidee_init, xios_orchidee_ok=',xios_orchidee_ok

    ! Coherence test between flag and preprocessing key
#ifndef XIOS
    IF (xios_orchidee_ok) THEN
       ! Write error messages and stop the model
       WRITE(numout,*) 'Preprocessing key XIOS is missing to run ORCHIDEE with XIOS'
       WRITE(numout,*) 'Recompile with preprocessing flag XIOS or set XIOS_ORCHIDEE_OK=n in run.def'
       WRITE(numout,*) 'Fatal error from ORCHIDEE. STOP in xios_orchidee_init'
#ifdef CPP_PARA
       CALL MPI_ABORT(3)
#endif     
    END IF
#endif


    !
    !! 1. Set date and calendar information on the format needed by XIOS
    !

    ! Get the calendar from IOIPSL and captialize the string to be readable for XIOS
    CALL ioget_calendar(calendar_str)
    IF (calendar_str == 'gregorian') THEN
       calendar_str='Gregorian'
    ELSE IF (calendar_str == 'noleap') THEN
       calendar_str='NoLeap'
    ELSE IF (calendar_str == '360d') THEN
       calendar_str='D360'
    END IF

    ! Transform the time origin from julian days into year, month, day and seconds
    CALL ju2ymds(date0, year0, month0, day0, sec0)

    ! Transform the current date into character string as XIOS need
    WRITE(start_str,"(I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") year,month,day,0,0,0

    ! Transform the time origin date into character string as XIOS need
    WRITE(startorig_str,"(I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") year0,month0,day0,0,0,0

    IF (long_print) THEN
       WRITE(numout,*) 'In xios_orchidee_init : calendar_str=', calendar_str
       WRITE(numout,*) 'In xios_orchidee_init : start_date=',   start_str
       WRITE(numout,*) 'In xios_orchidee_init : time_origin=',  startorig_str
    END IF

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       !
       !! 2. Context initialization
       !
       CALL xios_context_initialize("orchidee",MPI_COMM_ORCH)
       CALL xios_get_handle("orchidee",ctx_hdl_orchidee)
       CALL xios_set_current_context(ctx_hdl_orchidee)

       !
       !! 2. Calendar and date definition
       !
       CALL xios_set_context_attr("orchidee",calendar_type=calendar_str)
       CALL xios_set_context_attr("orchidee",start_date=start_str)
       CALL xios_set_context_attr("orchidee",time_origin=startorig_str)

       !
       !! 3. Domain definition
       !
       ! Global domain
       CALL xios_set_domain_attr("domain_landpoints",ni_glo=iim_g, nj_glo=jjm_g)

       ! Local MPI domain
       CALL xios_set_domain_attr("domain_landpoints",ibegin=1, iend=iim_g, jbegin=jj_begin, jend=jj_end)

       ! Define how data is stored on memory : 1D array for only continental points
       CALL xios_set_domain_attr("domain_landpoints",data_dim=1, data_ibegin=0, data_ni=nbp_mpi)
       CALL xios_set_domain_attr("domain_landpoints",data_n_index=nbp_mpi, data_i_index=kindex_mpi)     

       ! Define longitudes and latitudes on local MPI domain
       CALL xios_set_domain_attr("domain_landpoints",lonvalue=lon_mpi(:,1),latvalue=lat_mpi(1,:))

       !
       !! 4. Axis definition
       !
       CALL xios_set_axis_attr("veget",size=nvm ,VALUE=(/(REAL(i,r_std),i=1,nvm)/))
       CALL xios_set_axis_attr("laiax", size=nlai+1,VALUE=(/(REAL(i,r_std),i=1,nlai+1)/))
       CALL xios_set_axis_attr("solth",size=ngrnd ,VALUE=soilth_lev(:))
       CALL xios_set_axis_attr("soiltyp", size=nstm,VALUE=(/(REAL(i,r_std),i=1,nstm)/))
       CALL xios_set_axis_attr("nobio", size=nnobio,VALUE=(/(REAL(i,r_std),i=1,nnobio)/))
       CALL xios_set_axis_attr("albtyp", size=2,VALUE=(/(REAL(i,r_std),i=1,2)/))
       CALL xios_set_axis_attr("solay", size=nslm,VALUE=(/(REAL(i,r_std),i=1,nslm)/))
       CALL xios_set_axis_attr("soildiag", size=nbdl,VALUE=diaglev(:))
       CALL xios_set_axis_attr("snowlev", size=nsnow,VALUE=(/(REAL(i,r_std),i=1,nsnow)/))
       CALL xios_set_axis_attr("PFT", size=nvm,VALUE=(/(REAL(i,r_std), i=1,nvm)/))
       CALL xios_set_axis_attr("P10", size=10,VALUE=(/(REAL(i,r_std), i=1,10)/))
       CALL xios_set_axis_attr("P100", size=100,VALUE=(/(REAL(i,r_std), i=1,100)/))
       CALL xios_set_axis_attr("P11", size=11,VALUE=(/(REAL(i,r_std), i=1,11)/))
       CALL xios_set_axis_attr("P101", size=101,VALUE=(/(REAL(i,r_std), i=1,101)/))

       !
       !! 5. Send the time-step length to XIOS
       !
       dtime_xios%second=dt_sechiba
       CALL xios_set_timestep(dtime_xios) 

       !
       !! 6. Deactivation of some fields if they are not calculated
       !
       IF ( .NOT. control_flags%river_routing ) THEN
          CALL xios_set_field_attr("basinmap",enabled=.FALSE.)
          CALL xios_set_field_attr("nbrivers",enabled=.FALSE.)
          CALL xios_set_field_attr("riversret",enabled=.FALSE.)
          CALL xios_set_field_attr("hydrographs",enabled=.FALSE.)
          CALL xios_set_field_attr("fastr",enabled=.FALSE.)
          CALL xios_set_field_attr("slowr",enabled=.FALSE.)
          CALL xios_set_field_attr("streamr",enabled=.FALSE.)
          CALL xios_set_field_attr("lakevol",enabled=.FALSE.)
          CALL xios_set_field_attr("pondr",enabled=.FALSE.)
       END IF


       IF (control_flags%hydrol_cwrr ) THEN
          CALL xios_set_field_attr("dss",enabled=.FALSE.)
          CALL xios_set_field_attr("gqsb",enabled=.FALSE.)
          CALL xios_set_field_attr("bqsb",enabled=.FALSE.)
          CALL xios_set_field_attr("rsol",enabled=.FALSE.)
          CALL xios_set_field_attr("mrsos",enabled=.FALSE.) 
          CALL xios_set_field_attr("mrso",enabled=.FALSE.) 
          CALL xios_set_field_attr("mrros",enabled=.FALSE.) 
          CALL xios_set_field_attr("mrro",enabled=.FALSE.) 
          CALL xios_set_field_attr("prveg",enabled=.FALSE.) 
       ELSE
          CALL xios_set_field_attr("reinf_slope",enabled=.FALSE.)
          CALL xios_set_field_attr("soilindex",enabled=.FALSE.)
          CALL xios_set_field_attr("evapnu_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("drainage_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("transpir_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("runoff_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("humtot",enabled=.FALSE.)
          CALL xios_set_field_attr("humtot_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("SWI",enabled=.FALSE.)
          CALL xios_set_field_attr("njsc",enabled=.FALSE.)
          CALL xios_set_field_attr("k_litt",enabled=.FALSE.)
          CALL xios_set_field_attr("SoilMoist",enabled=.FALSE.)
          CALL xios_set_field_attr("moistc_1",enabled=.FALSE.)
          CALL xios_set_field_attr("moistc_2",enabled=.FALSE.)
          CALL xios_set_field_attr("moistc_3",enabled=.FALSE.)
          CALL xios_set_field_attr("kfactroot_1",enabled=.FALSE.)
          CALL xios_set_field_attr("kfactroot_2",enabled=.FALSE.)
          CALL xios_set_field_attr("kfactroot_3",enabled=.FALSE.)
          CALL xios_set_field_attr("vegetsoil_1",enabled=.FALSE.)
          CALL xios_set_field_attr("vegetsoil_2",enabled=.FALSE.)
          CALL xios_set_field_attr("vegetsoil_3",enabled=.FALSE.)
       END IF

       IF ( .NOT. control_flags%do_floodplains ) THEN
          CALL xios_set_field_attr("flood_frac",enabled=.FALSE.)
          CALL xios_set_field_attr("reinfiltration",enabled=.FALSE.)
          CALL xios_set_field_attr("floodmap",enabled=.FALSE.)
          CALL xios_set_field_attr("floodh",enabled=.FALSE.)       
          CALL xios_set_field_attr("floodr",enabled=.FALSE.)       
          CALL xios_set_field_attr("floodout",enabled=.FALSE.)       
          CALL xios_set_field_attr("evapflo",enabled=.FALSE.) 
          CALL xios_set_field_attr("evapflo_alma",enabled=.FALSE.) 
       END IF

       ! Deactivate some stomate fields. 
       ! These fields were traditionally added in sechiba_history.nc output file.
       IF ( .NOT. control_flags%ok_stomate ) THEN
          CALL xios_set_field_attr("nee",enabled=.FALSE.)
          CALL xios_set_field_attr("maint_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("hetero_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("growth_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("npp",enabled=.FALSE.)
       END IF

       IF ( .NOT. control_flags%do_irrigation ) THEN
          CALL xios_set_field_attr("irrigation",enabled=.FALSE.)
          CALL xios_set_field_attr("netirrig",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigation_alma",enabled=.FALSE.)
          CALL xios_set_field_attr("netirrig_alma",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigmap",enabled=.FALSE.)
       END IF

       IF ( .NOT. control_flags%ok_co2)THEN
          CALL xios_set_field_attr("vbetaco2",enabled=.FALSE.)
          CALL xios_set_field_attr("cimean",enabled=.FALSE.)
          CALL xios_set_field_attr("gpp",enabled=.FALSE.)
       END IF

       IF ( .NOT. control_flags%ok_inca)THEN
          CALL xios_set_field_attr("PAR",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_fertil_no",enabled=.FALSE.)
          CALL xios_set_field_attr("CRF",enabled=.FALSE.)
          CALL xios_set_field_attr("ptnlev1",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_iso",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_mono",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_ORVOC",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_MBO",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_methanol",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetone",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetal",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_formal",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetic",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_formic",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_no_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_no",enabled=.FALSE.)
       END IF

       IF ( .NOT. control_flags%ok_radcanopy)THEN
          CALL xios_set_field_attr("PARsun",enabled=.FALSE.)
          CALL xios_set_field_attr("PARsh",enabled=.FALSE.)
          CALL xios_set_field_attr("laisun",enabled=.FALSE.)
          CALL xios_set_field_attr("laish",enabled=.FALSE.)
          CALL xios_set_field_attr("Fdf",enabled=.FALSE.)
          CALL xios_set_field_attr("PARsuntab",enabled=.FALSE.)
          CALL xios_set_field_attr("PARshtab",enabled=.FALSE.)
          CALL xios_set_field_attr("Sinang",enabled=.FALSE.)
          CALL xios_set_field_attr("PARdf",enabled=.FALSE.)
          CALL xios_set_field_attr("PARdr",enabled=.FALSE.)
          CALL xios_set_field_attr("Trans",enabled=.FALSE.)
          CALL xios_set_field_attr("Day",enabled=.FALSE.)
          CALL xios_set_field_attr("Year_length",enabled=.FALSE.)
       END IF

       IF ( .NOT. control_flags%ok_bbgfertil_Nox) THEN
          CALL xios_set_field_attr("flx_co2_bbg_year",enabled=.FALSE.)
       END IF

       IF ( .NOT. control_flags%ok_cropsfertil_Nox) THEN
          CALL xios_set_field_attr("N_qt_WRICE_year",enabled=.FALSE.)
          CALL xios_set_field_attr("N_qt_OTHER_year",enabled=.FALSE.)
       END IF

       IF (.NOT. check_waterbal) THEN
          CALL xios_set_field_attr("TotWater",enabled=.FALSE.)
          CALL xios_set_field_attr("TotWaterFlux",enabled=.FALSE.)
       END IF


       !
       !! 7. Close context
       !
       CALL xios_close_context_definition()      


       !
       !! 8. Activate almaoutput if needed 
       !! Some extra calculations have to be done for the variables tot_watsoil_end("RootMoist"), 
       !! delsoilmoist("DelSoilMoist"), delintercept, delswe("DelSWE") and soilwet("SoilWet").
       !! Set almaoutput=true if at least one of these variables are defined in an output file. 
       !! If not, keep the initial value of almaoutput. 
       !

       IF (xios_field_is_active("RootMoist") .OR. xios_field_is_active("DelSoilMoist") .OR. &
            xios_field_is_active("DelIntercept") .OR. xios_field_is_active("DelSWE") .OR. &
            xios_field_is_active("SoilWet")) THEN

          almaoutput=.TRUE.
          WRITE(numout,*) 'almaoutput has been set to true in xios_orchidee_init'
       END IF
#endif
    END IF

    IF (xios_orchidee_ok) THEN
       ! Send variable almaoutput to all processes
       CALL bcast(almaoutput)
    END IF

    IF (long_print) WRITE(numout,*) 'Exit xios_orchidee_init'
  END SUBROUTINE xios_orchidee_init


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_change_context
  !!
  !>\BRIEF         Use this subroutine to switch between different context.
  !!               This subroutine must be called when running in coupled mode at each time ORCHIDEE is called, in the
  !!               begining and end of intersurf_gathered. First call is done after xios_orchidee_init is done. 
  !!
  !! DESCRIPTION  :\n 
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_change_context(new_context)
    !
    !! 0. Variable and parameter declaration
    !
    !!    Input variable
    CHARACTER(LEN=*),INTENT(IN)              :: new_context

    !! Local variables
#ifdef XIOS
    TYPE(xios_context) :: ctx_hdl
#endif
    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       CALL xios_get_handle(new_context,ctx_hdl)
       CALL xios_set_current_context(ctx_hdl)
#endif
    END IF
    
  END SUBROUTINE xios_orchidee_change_context

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_update_calendar
  !!
  !>\BRIEF          Update the calandar in XIOS.
  !!
  !! DESCRIPTION  :\n Update the calendar in XIOS : let XIOS know that ORCHIDEE avanced one time-step.
  !!                  This subroutine should be called in the beginning of each time-step. The first 
  !!                  time-step in a new execution should always start at 1. Therefore, first calculate
  !!                  an offset that is substracted to the current time step in sechiba. 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_update_calendar(itau_sechiba)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(IN) :: itau_sechiba    !! Current time step of the model
    !
    !! 0.2 Local variables
    !
    LOGICAL, SAVE         :: first=.TRUE.         !! Flag for first entering in subroutine
    INTEGER(i_std), SAVE  :: offset               !! Offset to substract from itau_sechiba
    INTEGER(i_std)        :: itau_xios            !! Current time step for XIOS

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       ! Calculate the offset
       IF (first) THEN
          offset=itau_sechiba-1
          first=.FALSE.
       END IF

       ! Substract the offset to the current time step in sechiba
       itau_xios=itau_sechiba-offset

       ! Send the new time step to XIOS
       IF (long_print) WRITE(numout,*) 'xios_orchidee_update_calendar: itau_sechiba, itau_xios=',itau_sechiba,itau_xios
       CALL xios_update_calendar(itau_xios)
#endif
    END IF
  END SUBROUTINE xios_orchidee_update_calendar
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_context_finalize
  !!
  !>\BRIEF         Finalize orchidee context.
  !!
  !! DESCRIPTION  :\n This subroutine finalizes the orchidee context without finalizing XIOS. In coupled mode, the atmospheric
  !!                  modele must finalize XIOS. This subroutine is called in the end of the execution of ORCHIDEE only in 
  !!                  coupeld mode.
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_context_finalize

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
       IF (long_print) WRITE(numout,*) 'Entering xios_orchidee_context_finalize'
#ifdef XIOS
       CALL xios_context_finalize()
#endif
    END IF
  END SUBROUTINE xios_orchidee_context_finalize


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_finalize
  !!
  !>\BRIEF         Last call to XIOS for finalization.
  !!
  !! DESCRIPTION  :\n Last call to XIOS for finalization of the orchidee context and XIOS.
  !!                  This subroutine is called only when ORCHIDEE is run in offline mode. In coupled mode it is the atmospheric
  !!                  model that finalizes XIOS. In that case, the context orchidee must be finalized using the 
  !!                  subroutine xios_orchidee_context_finalize
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_finalize

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
       IF (long_print) WRITE(numout,*) 'Entering xios_orchidee_finalize'
#ifdef XIOS
       CALL xios_context_finalize()
       CALL xios_finalize()
#endif
    END IF
  END SUBROUTINE xios_orchidee_finalize


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r1d
  !!
  !>\BRIEF          Subroutine for sending 1D (array) fields to XIOS.
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 1D fields (array).
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r1d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)          :: field_id
    REAL(r_std), DIMENSION(:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (long_print) WRITE(numout,*) 'Entering xios_orchidee_send_field_r1d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r1d


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r2d
  !!
  !>\BRIEF          Subroutine for sending 2D fields to XIOS.
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 2D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r2d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)            :: field_id
    REAL(r_std), DIMENSION(:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (long_print) WRITE(numout,*) 'Entering xios_orchidee_send_field_r2d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r2d


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r3d
  !!
  !>\BRIEF          Subroutine for sending 3D fields to XIOS. 
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 3D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r3d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (long_print) WRITE(numout,*) 'Entering xios_orchidee_send_field_r3d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r3d
 
END MODULE xios_orchidee

