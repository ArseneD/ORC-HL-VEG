!! This subroutine is the interface between the main program 
!! (LMDZ or dim2_driver) and SECHIBA.
!! - Input fields are gathered to keep just continental points
!! - call sechiba_main That's SECHIBA process.
!! - Output fields are scattered to complete global fields
!!
!! @call sechiba_main
!! @Version : $Revision: 2291 $, $Date: 2014-09-08 15:25:21 +0200 (Mon, 08 Sep 2014) $
!!
!! @author Marie-Alice Foujols and Jan Polcher
!! 
!< $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_sechiba/intersurf.f90 $
!< $Date: 2014-09-08 15:25:21 +0200 (Mon, 08 Sep 2014) $
!< $Author: josefine.ghattas $
!< $Revision: 2291 $
!! IPSL (2006)
!!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!!
!f90doc MODULEintersurf
MODULE intersurf

  USE IOIPSL
  USE xios_orchidee
  USE ioipsl_para 
  USE defprec
  USE sechiba
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE mod_orchidee_para
  USE watchout
  USE solar
  USE grid
  USE Orch_Write_field_p
  USE thermosoil, ONLY : thermosoil_levels


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Init_intersurf, intersurf_main,  intersurf_gathered, stom_define_history, stom_IPCC_define_history, intsurf_time


  INTERFACE intersurf_main
    MODULE PROCEDURE intersurf_main_2d, intersurf_main_1d, intersurf_gathered
  END INTERFACE



  !
  !  Global variables
  !
  INTEGER(i_std),PARAMETER                           :: max_hist_level = 11
  !
  LOGICAL, SAVE                                     :: l_first_intersurf=.TRUE. !! Initialisation has to be done one time
!$OMP THREADPRIVATE(l_first_intersurf)
  !
  INTEGER(i_std), SAVE                               :: hist_id, rest_id        !! IDs for history and restart files
!$OMP THREADPRIVATE(hist_id, rest_id)
  INTEGER(i_std), SAVE                               :: hist2_id                !! ID for the second history files (Hi-frequency ?)
!$OMP THREADPRIVATE(hist2_id)
  INTEGER(i_std), SAVE                               :: hist_id_stom, hist_id_stom_IPCC, rest_id_stom !! Dito for STOMATE
!$OMP THREADPRIVATE(hist_id_stom, hist_id_stom_IPCC, rest_id_stom)
  REAL(r_std), SAVE                                  :: dw                      !! frequency of history write (sec.)
!$OMP THREADPRIVATE(dw)
  !
  INTEGER(i_std), SAVE                               :: itau_offset  !! This offset is used to phase the 
  !                                                                 !! calendar of the GCM or the driver.
!$OMP THREADPRIVATE(itau_offset)
  REAL(r_std), SAVE                                  :: date0_shifted
!$OMP THREADPRIVATE(date0_shifted)
  !
  TYPE(control_type), SAVE                          :: control_flags !! Flags that (de)activate parts of the model
!$OMP THREADPRIVATE(control_flags)
  !
  !
  !! first day of this year
  REAL(r_std), SAVE :: julian0
!$OMP THREADPRIVATE(julian0)
  !
  LOGICAL, PARAMETER :: check_INPUTS = .false.         !! (very) long print of INPUTs in intersurf 
  LOGICAL, SAVE :: check_time = .FALSE.
!$OMP THREADPRIVATE(check_time)
  LOGICAL, SAVE :: impose_param = .TRUE.  !! Flag impos_param : should we read all the parameters in the run.def file ?
!$OMP THREADPRIVATE(impose_param)
  LOGICAL, SAVE :: ok_histsync            !! Flag activate syncronization of IOIPSL output
!$OMP THREADPRIVATE(ok_histsync)
  

  PUBLIC check_time, l_first_intersurf, impose_param
  !
CONTAINS
  !
  !f90doc CONTAINS
  ! 

  SUBROUTINE intersurf_main_2d (kjit, iim, jjm, kjpindex, kindex, xrdt, &
     & lrestart_read, lrestart_write, lon, lat, zcontfrac, zneighbours, zresolution, date0, &
! First level conditions
     & zlev, u, v, qair, temp_air, epot_air, ccanopy, &
! Variables for the implicit coupling
     & cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
! Rain, snow, radiation and surface pressure
     & precip_rain, precip_snow, lwdown, swnet, swdown, pb, &
! Output : Fluxes
     & vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
! Surface temperatures and surface properties
     & tsol_rad, temp_sol_new, qsurf, albedo, emis, z0, &
! For VOC  radiation  
     & sinang)
    ! routines called : sechiba_main
    !
    IMPLICIT NONE
    !    
    ! interface description for dummy arguments
    ! input scalar 
    INTEGER(i_std),INTENT (in)                            :: kjit            !! Time step number
    INTEGER(i_std),INTENT (in)                            :: iim, jjm        !! Dimension of input fields
    INTEGER(i_std),INTENT (in)                            :: kjpindex        !! Number of continental points
    REAL(r_std),INTENT (in)                               :: xrdt            !! Time step in seconds
    LOGICAL, INTENT (in)                                 :: lrestart_read    !! Logical for _restart_ file to read
    LOGICAL, INTENT (in)                                 :: lrestart_write   !! Logical for _restart_ file to write'
    REAL(r_std), INTENT (in)                              :: date0           !! Date at which kjit = 0
    ! input fields
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: kindex          !! Index for continental points
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: u             !! Lowest level wind speed
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: v             !! Lowest level wind speed 
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: zlev          !! Height of first layer
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: qair          !! Lowest level specific humidity
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: precip_rain   !! Rain precipitation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: precip_snow   !! Snow precipitation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: lwdown        !! Down-welling long-wave flux 
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: swnet         !! Net surface short-wave flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: swdown        !! Downwelling surface short-wave flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: sinang        !!
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: temp_air      !! Air temperature in Kelvin
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: epot_air      !! Air potential energy
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: ccanopy       !! CO2 concentration in the canopy
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: petAcoef      !! Coeficients A from the PBL resolution
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: peqAcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: petBcoef      !! Coeficients B from the PBL resolution
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: peqBcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (iim,jjm), INTENT(inout)          :: cdrag         !! Cdrag
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: pb            !! Surface pressure
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: lon, lat      !! Geographical coordinates
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in)             :: zcontfrac     !! Fraction of continent in the grid
    INTEGER, DIMENSION (iim,jjm,8), INTENT(in)              :: zneighbours   !! land neighbours
    REAL(r_std),DIMENSION (iim,jjm,2), INTENT(in)           :: zresolution   !! resolution in x and y dimensions
    ! output fields
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: z0            !! Surface roughness
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: coastalflow   !! Diffuse flow of water into the ocean (m^3/dt)
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: riverflow     !! Largest rivers flowing into the ocean (m^3/dt)
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: tsol_rad      !! Radiative surface temperature
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: vevapp        !! Total of evaporation
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: temp_sol_new  !! New soil temperature
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: qsurf         !! Surface specific humidity
    REAL(r_std),DIMENSION (iim,jjm,2), INTENT(out)          :: albedo        !! Albedo
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: fluxsens      !! Sensible chaleur flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: fluxlat       !! Latent chaleur flux
    REAL(r_std),DIMENSION (iim,jjm), INTENT(out)            :: emis          !! Emissivity
    ! LOCAL declaration
    ! work arrays to scatter and/or gather information just before/after sechiba_main call's
    ! and to keep output value for next call
    REAL(r_std),DIMENSION (kjpindex)                      :: zu            !! Work array to keep u
    REAL(r_std),DIMENSION (kjpindex)                      :: zv            !! Work array to keep v
    REAL(r_std),DIMENSION (kjpindex)                      :: zzlev         !! Work array to keep zlev
    REAL(r_std),DIMENSION (kjpindex)                      :: zqair         !! Work array to keep qair
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_rain  !! Work array to keep precip_rain
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_snow  !! Work array to keep precip_snow
    REAL(r_std),DIMENSION (kjpindex)                      :: zlwdown       !! Work array to keep lwdown
    REAL(r_std),DIMENSION (kjpindex)                      :: zswnet        !! Work array to keep swnet
    REAL(r_std),DIMENSION (kjpindex)                      :: zswdown       !! Work array to keep swdown
    REAL(r_std),DIMENSION (kjpindex)                      :: zsinang       !! Work array to keep sinang
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_air     !! Work array to keep temp_air
    REAL(r_std),DIMENSION (kjpindex)                      :: zepot_air     !! Work array to keep epot_air
    REAL(r_std),DIMENSION (kjpindex)                      :: zccanopy      !! Work array to keep ccanopy
    REAL(r_std),DIMENSION (kjpindex)                      :: zpetAcoef     !! Work array to keep petAcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpeqAcoef     !! Work array to keep peqAcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpetBcoef     !! Work array to keep petBcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpeqBcoef     !! Work array to keep peqVcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zcdrag        !! Work array to keep cdrag
    REAL(r_std),DIMENSION (kjpindex)                      :: zpb           !! Work array to keep pb
    REAL(r_std),DIMENSION (kjpindex)                      :: zz0           !! Work array to keep z0
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoastal      !! Work array to keep coastalflow
    REAL(r_std),DIMENSION (kjpindex)                      :: zriver        !! Work array to keep riverflow
    REAL(r_std),DIMENSION (kjpindex)                      :: dcoastal      !! Work array to keep coastalflow
    REAL(r_std),DIMENSION (kjpindex)                      :: driver        !! Work array to keep riverflow
    REAL(r_std),DIMENSION (kjpindex)                      :: znetco2       !! Work array to keep netco2flux
    REAL(r_std),DIMENSION (kjpindex)                      :: zcarblu       !! Work array to keep fco2_land_use
    REAL(r_std),DIMENSION (kjpindex)                      :: ztsol_rad     !! Work array to keep tsol_rad
    REAL(r_std),DIMENSION (kjpindex)                      :: zvevapp       !! Work array to keep vevapp
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_sol_new !! Work array to keep temp_sol_new
    REAL(r_std),DIMENSION (kjpindex)                      :: zqsurf        !! Work array to keep qsurf
    REAL(r_std),DIMENSION (kjpindex,2)                    :: zalbedo       !! Work array to keep albedo
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxsens     !! Work array to keep fluxsens
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxlat      !! Work array to keep fluxlat
    REAL(r_std),DIMENSION (kjpindex)                      :: zemis         !! Work array to keep emis
    REAL(r_std),ALLOCATABLE, DIMENSION (:)                :: soilth_lev    !! Vertical soil axis for thermal scheme (m)
    !
    INTEGER(i_std)                                       :: i, j, ik
    INTEGER(i_std)                                       :: ier
    INTEGER(i_std)                                       :: itau_sechiba
    REAL(r_std)                                          :: zlev_mean
    LOGICAL                                              :: do_watch      !! if it's time, write watchout
    INTEGER                                              :: old_fileout   !! old Logical Int for std IO output
    LOGICAL :: check = .FALSE.
    !
    CALL ipslnlf(new_number=numout,old_number=old_fileout)
    !
    IF (l_first_intersurf) THEN
       !
       CALL intsurf_time( kjit, date0, xrdt )
       !
       IF ( check ) WRITE(numout,*) 'Initialisation of intersurf_main_2d'
       !
       OFF_LINE_MODE = .TRUE. 
       !
       DO ik=1,kjpindex
          
          j = ((kindex(ik)-1)/iim) + 1
          i = (kindex(ik) - (j-1)*iim)

          !- Create the internal coordinate table
          !-
          lalo(ik,1) = lat(i,j)
          lalo(ik,2) = lon(i,j)
          !
          !- Store the fraction of the continents only once so that the user
          !- does not change them afterwards.
          !-
          contfrac(ik) = zcontfrac(i,j)
       ENDDO
       CALL gather(contfrac,contfrac_g)
       CALL gather(lalo,lalo_g)
       CALL gather2D_mpi(lon,lon_g)
       CALL gather2D_mpi(lat,lat_g)
       CALL gather2D_mpi(zlev,zlev_g)
       !
       !  Configuration of SSL specific parameters
       !
       CALL intsurf_config(control_flags, xrdt)
       !
       CALL intsurf_restart(kjit, date0, xrdt, control_flags, rest_id, rest_id_stom, itau_offset)
       itau_sechiba = kjit + itau_offset
       !


       !!- Initialize module for output with XIOS
       !
       ! Get the vertical soil levels for the thermal scheme, to be used in xios_orchidee_init
       ! Get the vertical soil levels for the thermal scheme, to be used in xios_orchidee_init
       ALLOCATE(soilth_lev(ngrnd), stat=ier)
       IF (ier /= 0) THEN
          CALL ipslerr_p(3,'intersurf_main_2d', 'Error in allocation of soilth_lev','','')
       END IF
       soilth_lev(1:ngrnd) = thermosoil_levels()
       
       CALL xios_orchidee_init( MPI_COMM_ORCH,                &
            date0,    year,    month,           day,          &
            lon,      lat,     soilth_lev,      control_flags )

       !- Initialize IOIPSL output files
       CALL intsurf_history(iim, jjm, lon, lat,  kindex, kjpindex, itau_sechiba, date0_shifted, xrdt, control_flags, hist_id, &
            & hist2_id, hist_id_stom, hist_id_stom_IPCC)

       IF ( ok_watchout ) THEN
          IF (is_root_prc) THEN
             zlev_mean = zero
             DO ik=1, nbp_glo
                j = ((index_g(ik)-1)/iim_g) + 1
                i = (index_g(ik) - (j-1)*iim_g)
       
                zlev_mean = zlev_mean + zlev_g(i,j)
             ENDDO
             zlev_mean = zlev_mean / REAL(nbp_glo,r_std)
          ENDIF

          last_action_watch = itau_sechiba
          last_check_watch  = last_action_watch

          ! Only root proc write watchout file
          CALL watchout_init (iim_g, jjm_g, kjpindex, nbp_glo, &
               & date0_shifted, last_action_watch, dt_watch, index_g, lon_g, lat_g, zlev_mean)
       ENDIF
       !
       IF ( check ) WRITE(numout,*) 'End of Initialisation of intersurf'
       !
    ENDIF
    !
    !  Shift the time step to phase the two models
    !
    itau_sechiba = kjit + itau_offset
    !

    ! Update the calendar in xios by sending the new time step
    IF (l_first_intersurf) THEN
       ! Special case : the model is only in initialization phase and the current itau_sechiba is not a real time step. 
       ! Therefor give itau_sechiba+1 to xios to have a correct time axis in output files. 
       CALL xios_orchidee_update_calendar(itau_sechiba+1)
    ELSE
       CALL xios_orchidee_update_calendar(itau_sechiba)
    ENDIF

    CALL intsurf_time( itau_sechiba, date0_shifted, xrdt )
    !
    ! 1. gather input fields from kindex array
    !    Warning : I'm not sure this interface with one dimension array is the good one
    !
    DO ik=1, kjpindex
      
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)
       
       zu(ik)           = u(i,j)
       zv(ik)           = v(i,j)
       zzlev(ik)        = zlev(i,j)
       zqair(ik)        = qair(i,j)
       zprecip_rain(ik) = precip_rain(i,j)*xrdt
       zprecip_snow(ik) = precip_snow(i,j)*xrdt
       zlwdown(ik)      = lwdown(i,j)
       zswnet(ik)       = swnet(i,j)
       zswdown(ik)      = swdown(i,j)
       zsinang(ik)      = sinang(i,j)
       ztemp_air(ik)    = temp_air(i,j)
       zepot_air(ik)    = epot_air(i,j)
       zccanopy(ik)     = ccanopy(i,j)
       zpetAcoef(ik)    = petAcoef(i,j)
       zpeqAcoef(ik)    = peqAcoef(i,j)
       zpetBcoef(ik)    = petBcoef(i,j)
       zpeqBcoef(ik)    = peqBcoef(i,j)
       zcdrag(ik)       = cdrag(i,j)
       zpb(ik)          = pb(i,j)
       
    ENDDO
    !
    IF (check_INPUTS) THEN
       WRITE(numout,*) "Intersurf_main_2D :"
       WRITE(numout,*) "Time step number = ",kjit
       WRITE(numout,*) "Dimension of input fields = ",iim, jjm
       WRITE(numout,*) "Number of continental points = ",kjpindex
       WRITE(numout,*) "Time step in seconds = ",xrdt
       WRITE(numout,*) "Logical for _restart_ file to read, write = ",lrestart_read,lrestart_write
       WRITE(numout,*) "Date at which kjit = 0  =  ",date0
       WRITE(numout,*) "Index for continental points = ",kindex
       WRITE(numout,*) "Lowest level wind speed North = ",zu
       WRITE(numout,*) "Lowest level wind speed East = ",zv
       WRITE(numout,*) "Height of first layer = ",zzlev
       WRITE(numout,*) "Lowest level specific humidity = ",zqair
       WRITE(numout,*) "Rain precipitation = ",zprecip_rain
       WRITE(numout,*) "Snow precipitation = ",zprecip_snow
       WRITE(numout,*) "Down-welling long-wave flux = ",zlwdown
       WRITE(numout,*) "Net surface short-wave flux = ",zswnet
       WRITE(numout,*) "Downwelling surface short-wave flux = ",zswdown
       WRITE(numout,*) "Air temperature in Kelvin = ",ztemp_air
       WRITE(numout,*) "Air potential energy = ",zepot_air
       WRITE(numout,*) "CO2 concentration in the canopy = ",zccanopy
       WRITE(numout,*) "Coeficients A from the PBL resolution = ",zpetAcoef
       WRITE(numout,*) "One for T and another for q = ",zpeqAcoef
       WRITE(numout,*) "Coeficients B from the PBL resolution = ",zpetBcoef
       WRITE(numout,*) "One for T and another for q = ",zpeqBcoef
       WRITE(numout,*) "Cdrag = ",zcdrag
       WRITE(numout,*) "Surface pressure = ",zpb
       WRITE(numout,*) "Geographical coordinates lon = ", (/ ( lon(ilandindex(ik), jlandindex(ik)), ik=1,kjpindex ) /)
       WRITE(numout,*) "Geographical coordinates lat = ", (/ ( lat(ilandindex(ik), jlandindex(ik)), ik=1,kjpindex ) /) 
       WRITE(numout,*) "Fraction of continent in the grid = ",contfrac
    ENDIF
    !
    ! 2. save the grid
    !
    IF ( check ) WRITE(numout,*) 'Save the grid'
    !
    IF (l_first_intersurf) THEN
       CALL histwrite_p(hist_id, 'LandPoints',  itau_sechiba+1, (/ ( REAL(ik), ik=1,kjpindex ) /), kjpindex, kindex)
       CALL histwrite_p(hist_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       IF ( control_flags%ok_stomate ) THEN
          CALL histwrite_p(hist_id_stom, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
          CALL histwrite_p(hist_id_stom_IPCC, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       ENDIF

       CALL histwrite_p(hist_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)
       IF ( is_omp_root .AND. hist_id > 0 ) THEN
          ! Always syncronize output after initialization 
          CALL histsync(hist_id)
       END IF

       CALL histwrite_p(hist2_id, 'LandPoints',  itau_sechiba+1, (/ ( REAL(ik), ik=1,kjpindex ) /), kjpindex, kindex)
       CALL histwrite_p(hist2_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       CALL histwrite_p(hist2_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)
       IF ( is_omp_root .AND. hist2_id > 0 ) THEN
          ! Always syncronize output after initialization 
          CALL histsync(hist2_id)
       ENDIF
       !
    ENDIF
    !
    ! 3. call sechiba for continental points only
    !
    IF ( check ) WRITE(numout,*) 'Calling sechiba'
    !
    CALL sechiba_main (itau_sechiba, iim*jjm, kjpindex, kindex, xrdt, date0_shifted, &
       & lrestart_read, lrestart_write, control_flags, &
       & lalo, contfrac, neighbours, resolution, &
! First level conditions
! Ajout Nathalie - Juin 2006 - passage q2m/t2m pour calcul rveget
!       & zzlev, zu, zv, zqair, ztemp_air, zepot_air, zccanopy, &
       & zzlev, zu, zv, zqair, zqair, ztemp_air, ztemp_air, zepot_air, zccanopy, &
! Variables for the implicit coupling
       & zcdrag, zpetAcoef, zpeqAcoef, zpetBcoef, zpeqBcoef, &
! Rain, snow, radiation and surface pressure
       & zprecip_rain ,zprecip_snow,  zlwdown, zswnet, zswdown, zsinang, zpb, &
! Output : Fluxes
       & zvevapp, zfluxsens, zfluxlat, zcoastal, zriver, znetco2, zcarblu, &
! Surface temperatures and surface properties
       & ztsol_rad, ztemp_sol_new, zqsurf, zalbedo, zemis, zz0, &
! File ids
       & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC ) 
    
    !
    IF ( check ) WRITE(numout,*) 'out of SECHIBA'
    !
    ! 4. save watchout
    !
    IF ( ok_watchout .AND. .NOT. l_first_intersurf ) THEN
       ! Accumulate last time step
       sum_zlev(:) = sum_zlev(:) + zzlev(:)
       sum_u(:) = sum_u(:) + zu(:)
       sum_v(:) = sum_v(:) + zv(:)
       sum_qair(:) = sum_qair(:) + zqair(:) 
       sum_temp_air(:) = sum_temp_air(:) + ztemp_air(:)
       sum_epot_air(:) = sum_epot_air(:) + zepot_air(:)
       sum_ccanopy(:) = sum_ccanopy(:) + zccanopy(:)
       sum_cdrag(:) = sum_cdrag(:) + zcdrag(:)
       sum_petAcoef(:) = sum_petAcoef(:) + zpetAcoef(:)
       sum_peqAcoef(:) = sum_peqAcoef(:) + zpeqAcoef(:)
       sum_petBcoef(:) = sum_petBcoef(:) + zpetBcoef(:)
       sum_peqBcoef(:) = sum_peqBcoef(:) + zpeqBcoef(:)
       sum_rain(:) = sum_rain(:) + zprecip_rain(:)
       sum_snow(:) = sum_snow(:) + zprecip_snow(:)
       sum_lwdown(:) = sum_lwdown(:) + zlwdown(:)
       sum_pb(:) = sum_pb(:) + zpb(:)

!!$       IF ( dt_watch > 3600 ) THEN
!!$          julian_watch = date0_shifted+((itau_sechiba-0.5)/dt_split_watch)*dt_watch/one_day
!!$          WRITE(numout, *) "WATCH register : julian_watch ",julian_watch, " julian0",julian0,"date0_shifted ",date0_shifted, &
!!$               "itau_sechiba ",itau_sechiba, &
!!$               dt_split_watch,dt_watch,one_day
!!$          CALL solarang (julian_watch, julian0, iim, jjm, lon, lat, sinang)
!!$          WHERE ( sinang(:,:) .LT. EPSILON(un) ) 
!!$             isinang(:,:) = isinang(:,:) - 1
!!$          ENDWHERE
!!$          mean_sinang(:,:) = mean_sinang(:,:)+sinang(:,:)
!!$          WRITE(numout, *) "WATCH sinang : ",sinang, mean_sinang
!!$          WRITE(numout,*) "sum_swdown",sum_swdown
!!$          !
!!$          DO ik=1,kjpindex          
!!$             j = ((kindex(ik)-1)/iim) + 1
!!$             i = (kindex(ik) - (j-1)*iim)
!!$             
!!$             sum_swnet(ik) = sum_swnet(ik) + sinang(i,j)*zswnet(ik)
!!$             sum_swdown(ik) = sum_swdown(ik) + sinang(i,j)*zswdown(ik)
!!$          ENDDO
!!$       ELSE
          sum_swnet(:) = sum_swnet(:) + zswnet(:)
          sum_swdown(:) = sum_swdown(:) + zswdown(:)
!!$       ENDIF

       do_watch = .FALSE.
       call isittime &
            &  (itau_sechiba,date0_shifted,xrdt,dt_watch,&
            &   last_action_watch,last_check_watch,do_watch)
       last_check_watch = itau_sechiba
       IF (do_watch) THEN
          !
          IF ( check ) WRITE(numout,*) 'save watchout'
          !
          IF (long_print) THEN
             WRITE(numout,*) "intersurf : write watchout for date ",date0,date0_shifted,itau_sechiba,&
                  & last_action_watch, last_check_watch
          ENDIF
          last_action_watch = itau_sechiba

          sum_zlev(:) = sum_zlev(:) / dt_split_watch
          sum_u(:) = sum_u(:) / dt_split_watch
          sum_v(:) = sum_v(:) / dt_split_watch
          sum_qair(:) = sum_qair(:) / dt_split_watch
          sum_temp_air(:) = sum_temp_air(:) / dt_split_watch
          sum_epot_air(:) = sum_epot_air(:) / dt_split_watch
          sum_ccanopy(:) = sum_ccanopy(:) / dt_split_watch
          sum_cdrag(:) = sum_cdrag(:) / dt_split_watch
          sum_petAcoef(:) = sum_petAcoef(:) / dt_split_watch
          sum_peqAcoef(:) = sum_peqAcoef(:) / dt_split_watch
          sum_petBcoef(:) = sum_petBcoef(:) / dt_split_watch
          sum_peqBcoef(:) = sum_peqBcoef(:) / dt_split_watch
          sum_rain(:) = sum_rain(:) / dt_split_watch
          sum_snow(:) = sum_snow(:) / dt_split_watch
          sum_lwdown(:) = sum_lwdown(:) / dt_split_watch
          sum_pb(:) = sum_pb(:) / dt_split_watch

!!$          IF ( dt_watch > 3600 ) THEN 
!!$             WRITE(numout, *) "WATCH mean_sinang before norm : ",mean_sinang,isinang
!!$             WHERE ( isinang(:,:) .LT. dt_split_watch )
!!$                mean_sinang(:,:) = mean_sinang(:,:) / isinang(:,:)
!!$             ENDWHERE
!!$             WRITE(numout, *) "WATCH mean_sinang norm : ",mean_sinang
!!$             WRITE(numout,*) "SWDOWN 0 : ",sum_swdown(:)
!!$             !
!!$             DO ik=1,kjpindex          
!!$                j = ((kindex(ik)-1)/iim) + 1
!!$                i = (kindex(ik) - (j-1)*iim)
!!$                IF (mean_sinang(i,j) > zero) THEN
!!$                   sum_swdown(ik) = sum_swdown(ik)/mean_sinang(i,j)
!!$                   sum_swnet(ik) =  sum_swnet(ik)/mean_sinang(i,j)
!!$                ELSE
!!$                   sum_swdown(ik) = zero
!!$                   sum_swnet(ik) =  zero
!!$                ENDIF
!!$             ENDDO
!!$          ELSE
             sum_swnet(:) = sum_swnet(:) / dt_split_watch
             sum_swdown(:) = sum_swdown(:) / dt_split_watch
!!$          ENDIF

          CALL watchout_write_p(kjpindex, itau_sechiba, xrdt, sum_zlev, sum_swdown, sum_rain, &
               &   sum_snow, sum_lwdown, sum_pb, sum_temp_air, sum_epot_air, sum_qair, &
               &   sum_u, sum_v, sum_swnet, sum_petAcoef, sum_peqAcoef, sum_petBcoef, sum_peqBcoef, &
               &   sum_cdrag, sum_ccanopy )
       ENDIF
    ENDIF
    !
    ! 5. scatter output fields
    !
    z0(:,:)           = undef_sechiba
    coastalflow(:,:)  = undef_sechiba
    riverflow(:,:)    = undef_sechiba
    tsol_rad(:,:)     = undef_sechiba
    vevapp(:,:)       = undef_sechiba
    temp_sol_new(:,:) = undef_sechiba 
    qsurf(:,:)        = undef_sechiba 
    albedo(:,:,:)     = undef_sechiba
    fluxsens(:,:)     = undef_sechiba
    fluxlat(:,:)      = undef_sechiba
    emis(:,:)         = undef_sechiba 
    cdrag(:,:)        = undef_sechiba 
    !
    DO ik=1, kjpindex
      
    
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)

       z0(i,j)           = zz0(ik)
       coastalflow(i,j)  = zcoastal(ik)
       riverflow(i,j)    = zriver(ik)
       tsol_rad(i,j)     = ztsol_rad(ik)
       vevapp(i,j)       = zvevapp(ik)
       temp_sol_new(i,j) = ztemp_sol_new(ik)
       qsurf(i,j)        = zqsurf(ik)
       albedo(i,j,1)     = zalbedo(ik,1)
       albedo(i,j,2)     = zalbedo(ik,2)
       fluxsens(i,j)     = zfluxsens(ik)
       fluxlat(i,j)      = zfluxlat(ik)
       emis(i,j)         = zemis(ik)
       cdrag(i,j)        = zcdrag(ik)

    ENDDO
    !
    ! Modified fields for variables scattered during the writing
    !
    dcoastal(:) = zcoastal(:)
    driver(:)   = zriver(:)

    IF ( .NOT. l_first_intersurf) THEN
       
       CALL xios_orchidee_send_field("LandPoints" ,(/ ( REAL(ik), ik=1,kjpindex ) /))
       CALL xios_orchidee_send_field("Areas", area)
       CALL xios_orchidee_send_field("Contfrac",contfrac)
       CALL xios_orchidee_send_field("evap",zvevapp*one_day/dt_sechiba)
       CALL xios_orchidee_send_field("evap_alma",zvevapp/dt_sechiba)
       CALL xios_orchidee_send_field("coastalflow",dcoastal/dt_sechiba)
       CALL xios_orchidee_send_field("riverflow",driver/dt_sechiba)
       CALL xios_orchidee_send_field("temp_sol_C",ztemp_sol_new-ZeroCelsius)
       CALL xios_orchidee_send_field("temp_sol_K",ztemp_sol_new)
       CALL xios_orchidee_send_field("fluxsens",zfluxsens)
       CALL xios_orchidee_send_field("fluxlat",zfluxlat)
       CALL xios_orchidee_send_field("alb_vis",zalbedo(:,1))
       CALL xios_orchidee_send_field("alb_nir",zalbedo(:,2))
       CALL xios_orchidee_send_field("tair",ztemp_air)
       CALL xios_orchidee_send_field("qair",zqair)
       CALL xios_orchidee_send_field("q2m",zqair)
       CALL xios_orchidee_send_field("t2m",ztemp_air)
       CALL xios_orchidee_send_field("swnet",zswnet)
       CALL xios_orchidee_send_field("swdown",zswdown)
       

       IF ( .NOT. almaoutput ) THEN
       !
       !  scattered during the writing
       ! 
          CALL histwrite_p (hist_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'coastalflow',itau_sechiba, dcoastal, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'riverflow',itau_sechiba, driver, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'temp_sol', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'tsol_max', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'tsol_min', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'fluxlat',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'swnet',    itau_sechiba, zswnet, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'swdown',   itau_sechiba, zswdown, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
          CALL histwrite_p (hist_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
          CALL histwrite_p (hist_id, 'tair',     itau_sechiba, ztemp_air, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'qair',     itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'q2m',     itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist_id, 't2m',     itau_sechiba, ztemp_air, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'coastalflow',itau_sechiba, dcoastal, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'riverflow',itau_sechiba, driver, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'temp_sol', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tsol_max', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tsol_min', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'fluxlat',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'swnet',    itau_sechiba, zswnet, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'swdown',   itau_sechiba, zswdown, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tair',     itau_sechiba, ztemp_air, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'qair',     itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'q2m',     itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 't2m',     itau_sechiba, ztemp_air, kjpindex, kindex)
       ELSE
          CALL histwrite_p (hist_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'SWnet',    itau_sechiba, zswnet, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'AvgSurfT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'RadT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Tair', itau_sechiba, ztemp_air, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Qair', itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'SWnet',    itau_sechiba, zswnet, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'AvgSurfT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'RadT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       ENDIF
       !
       IF ( is_omp_root ) THEN
          IF ( (dw .EQ. xrdt) .AND. hist_id > 0 ) THEN
             ! Syncronize output but only if flag ok_histsync is set to true
             IF (ok_histsync) CALL histsync(hist_id)
          ENDIF
       END IF
   ENDIF
    !
    ! 6. Transform the water fluxes into Kg/m^2s and m^3/s
    !
    DO ik=1, kjpindex
    
       j = ((kindex(ik)-1)/iim) + 1
       i = (kindex(ik) - (j-1)*iim)

       vevapp(i,j) = vevapp(i,j)/xrdt
       coastalflow(i,j) = coastalflow(i,j)/xrdt
       riverflow(i,j) = riverflow(i,j)/xrdt

    ENDDO
    !
    IF ( lrestart_write .AND. ok_watchout .AND. is_root_prc ) THEN
       CALL watchout_close()
    ENDIF
    !
    ! Uncomment those lines to test all parallel getin in developments
!    IF(l_first_intersurf) THEN
!       CALL ipsldbg (new_status=.TRUE.)
!       CALL getin_dump_para("test",mpi_rank) 
!       CALL ipsldbg (new_status=.FALSE.)
!    ENDIF
    IF(l_first_intersurf .AND. is_root_prc) CALL getin_dump
    ! 
    l_first_intersurf = .FALSE.
    !
    IF (long_print) WRITE (numout,*) ' intersurf_main done '
    !
    CALL ipslnlf(new_number=old_fileout)
    !
  END SUBROUTINE intersurf_main_2d
!
  SUBROUTINE intersurf_main_1d (kjit, iim, jjm, kjpindex, kindex, xrdt, &
     & lrestart_read, lrestart_write, lon, lat, zcontfrac, zneighbours, zresolution, date0, &
! First level conditions
     & zlev,  u, v, qair, temp_air, epot_air, ccanopy, &
! Variables for the implicit coupling
     & cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
! Rain, snow, radiation and surface pressure
     & precip_rain, precip_snow, lwdown, swnet, swdown, pb, &
! Output : Fluxes
     & vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
! Surface temperatures and surface properties
     & tsol_rad, temp_sol_new, qsurf, albedo, emis, z0, &
! For VOC  radiation  
    & sinang)    
    ! routines called : sechiba_main
    !
    IMPLICIT NONE
    !    
    ! interface description for dummy arguments
    ! input scalar 
    INTEGER(i_std),INTENT (in)                            :: kjit          !! Time step number
    INTEGER(i_std),INTENT (in)                            :: iim, jjm      !! Dimension of input fields
    INTEGER(i_std),INTENT (in)                            :: kjpindex      !! Number of continental points
    REAL(r_std),INTENT (in)                               :: xrdt          !! Time step in seconds
    LOGICAL, INTENT (in)                                 :: lrestart_read !! Logical for _restart_ file to read
    LOGICAL, INTENT (in)                                 :: lrestart_write!! Logical for _restart_ file to write'
    REAL(r_std), INTENT (in)                              :: date0         !! Date at which kjit = 0
    ! input fields
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: kindex        !! Index for continental points
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: u             !! Lowest level wind speed
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: v             !! Lowest level wind speed 
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: zlev          !! Height of first layer
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: qair          !! Lowest level specific humidity
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: precip_rain   !! Rain precipitation
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: precip_snow   !! Snow precipitation
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: lwdown        !! Down-welling long-wave flux 
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: swnet         !! Net surface short-wave flux
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: swdown        !! Downwelling surface short-wave flux
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: sinang        !! Downwelling surface short-wave flux
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: temp_air      !! Air temperature in Kelvin
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: epot_air      !! Air potential energy
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: ccanopy       !! CO2 concentration in the canopy
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: petAcoef      !! Coeficients A from the PBL resolution
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: peqAcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: petBcoef      !! Coeficients B from the PBL resolution
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: peqBcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (iim*jjm), INTENT(inout)          :: cdrag         !! Cdrag
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: pb            !! Surface pressure
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: lon, lat      !! Geographical coordinates
    REAL(r_std),DIMENSION (iim*jjm), INTENT(in)             :: zcontfrac     !! Fraction of continent
    INTEGER, DIMENSION (iim*jjm,8), INTENT(in)             :: zneighbours   !! land neighbours
    REAL(r_std),DIMENSION (iim*jjm,2), INTENT(in)           :: zresolution   !! resolution in x and y dimensions
    ! output fields
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: z0            !! Surface roughness
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: coastalflow   !! Diffuse flow of water into the ocean (m^3/dt)
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: riverflow     !! Largest rivers flowing into the ocean (m^3/dt)
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: tsol_rad      !! Radiative surface temperature
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: vevapp        !! Total of evaporation
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: temp_sol_new  !! New soil temperature
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: qsurf         !! Surface specific humidity
    REAL(r_std),DIMENSION (iim*jjm,2), INTENT(out)          :: albedo        !! Albedo
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: fluxsens      !! Sensible chaleur flux
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: fluxlat       !! Latent chaleur flux
    REAL(r_std),DIMENSION (iim*jjm), INTENT(out)            :: emis          !! Emissivity
    ! LOCAL declaration
    ! work arrays to scatter and/or gather information just before/after sechiba_main call's
    ! and to keep output value for next call
    REAL(r_std),DIMENSION (kjpindex)                      :: zu            !! Work array to keep u
    REAL(r_std),DIMENSION (kjpindex)                      :: zv            !! Work array to keep v
    REAL(r_std),DIMENSION (kjpindex)                      :: zzlev         !! Work array to keep zlev
    REAL(r_std),DIMENSION (kjpindex)                      :: zqair         !! Work array to keep qair
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_rain  !! Work array to keep precip_rain
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_snow  !! Work array to keep precip_snow
    REAL(r_std),DIMENSION (kjpindex)                      :: zlwdown       !! Work array to keep lwdown
    REAL(r_std),DIMENSION (kjpindex)                      :: zswnet        !! Work array to keep swnet
    REAL(r_std),DIMENSION (kjpindex)                      :: zswdown       !! Work array to keep swdown
    REAL(r_std),DIMENSION (kjpindex)                      :: zsinang       !! Work array to keep sinang
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_air     !! Work array to keep temp_air
    REAL(r_std),DIMENSION (kjpindex)                      :: zepot_air     !! Work array to keep epot_air
    REAL(r_std),DIMENSION (kjpindex)                      :: zccanopy      !! Work array to keep ccanopy
    REAL(r_std),DIMENSION (kjpindex)                      :: zpetAcoef     !! Work array to keep petAcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpeqAcoef     !! Work array to keep peqAcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpetBcoef     !! Work array to keep petBcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zpeqBcoef     !! Work array to keep peqVcoef
    REAL(r_std),DIMENSION (kjpindex)                      :: zcdrag        !! Work array to keep cdrag
    REAL(r_std),DIMENSION (kjpindex)                      :: zpb           !! Work array to keep pb
    REAL(r_std),DIMENSION (kjpindex)                      :: zz0           !! Work array to keep z0
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoastal      !! Work array to keep coastal flow
    REAL(r_std),DIMENSION (kjpindex)                      :: zriver        !! Work array to keep river out flow
    REAL(r_std),DIMENSION (kjpindex)                      :: dcoastal      !! Work array to keep coastal flow
    REAL(r_std),DIMENSION (kjpindex)                      :: driver        !! Work array to keep river out flow
    REAL(r_std),DIMENSION (kjpindex)                      :: znetco2       !! Work array to keep netco2flux
    REAL(r_std),DIMENSION (kjpindex)                      :: zcarblu       !! Work array to keep fco2_land_use
    REAL(r_std),DIMENSION (kjpindex)                      :: ztsol_rad     !! Work array to keep tsol_rad
    REAL(r_std),DIMENSION (kjpindex)                      :: zvevapp       !! Work array to keep vevapp
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_sol_new !! Work array to keep temp_sol_new
    REAL(r_std),DIMENSION (kjpindex)                      :: zqsurf        !! Work array to keep qsurf
    REAL(r_std),DIMENSION (kjpindex,2)                    :: zalbedo       !! Work array to keep albedo
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxsens     !! Work array to keep fluxsens
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxlat      !! Work array to keep fluxlat
    REAL(r_std),DIMENSION (kjpindex)                      :: zemis         !! Work array to keep emis
    REAL(r_std),ALLOCATABLE, DIMENSION (:)                :: soilth_lev    !! Vertical soil axis for thermal scheme (m)
    !
    INTEGER(i_std)                                        :: i, j, ik
    INTEGER(i_std)                                        :: ier
    INTEGER(i_std)                                        :: itau_sechiba
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)              :: tmp_lon, tmp_lat, tmp_lev
    REAL(r_std)                                           :: zlev_mean
    LOGICAL                                               :: do_watch      !! if it's time, write watchout
    INTEGER                                               :: old_fileout   !! old Logical Int for std IO output
    LOGICAL :: check = .FALSE.
    !
    CALL ipslnlf(new_number=numout,old_number=old_fileout)
    !
    IF (l_first_intersurf) THEN
       !
       CALL intsurf_time( kjit, date0, xrdt )
       !
       IF ( check ) WRITE(numout,*) 'Initialisation of intersurf_main_1d'
       WRITE(numout,*) 'Initialisation of intersurf_main_1d'   
    !
       OFF_LINE_MODE = .TRUE. 
       !
       !  Create the internal coordinate table
       !
       IF ( (.NOT.ALLOCATED(tmp_lon))) THEN
          ALLOCATE(tmp_lon(iim,jjm))
       ENDIF
       IF ( (.NOT.ALLOCATED(tmp_lat))) THEN
          ALLOCATE(tmp_lat(iim,jjm))
       ENDIF
       IF ( (.NOT.ALLOCATED(tmp_lev))) THEN
          ALLOCATE(tmp_lev(iim,jjm))
       ENDIF
       !
       DO i=1,iim
          DO j=1,jjm
             ik = (j-1)*iim + i
             tmp_lon(i,j) = lon(ik)
             tmp_lat(i,j) = lat(ik)
             tmp_lev(i,j) = zlev(kindex(ik)) 
          ENDDO
       ENDDO
       !
       lalo(:,1) = lat(:)
       lalo(:,2) = lon(:)
       !
       !- Store the fraction of the continents only once so that the user
       !- does not change them afterwards.
       !
       DO ik=1,kjpindex

          contfrac(ik) = zcontfrac(kindex(ik))

       ENDDO
       contfrac_g(:) = contfrac(:)
       lalo_g(:,:) = lalo(:,:)
       lon_g(:,:) = tmp_lon(:,:)
       lat_g(:,:) = tmp_lat(:,:)
       zlev_g(:,:) = tmp_lev(:,:)
       !
       !  Configuration of SSL specific parameters
       !
       CALL intsurf_config(control_flags, xrdt)
       !
       CALL intsurf_restart(kjit, date0, xrdt, control_flags, rest_id, rest_id_stom, itau_offset)
       itau_sechiba = kjit + itau_offset
       !


       !!- Initialize module for output with XIOS
       !
       ! Get the vertical soil levels for the thermal scheme, to be used in xios_orchidee_init
       ALLOCATE(soilth_lev(ngrnd), stat=ier)
       IF (ier /= 0) THEN
          CALL ipslerr_p(3,'intersurf_main_1d', 'Error in allocation of soilth_lev','','')
       END IF
       soilth_lev(1:ngrnd) = thermosoil_levels()
       
       CALL xios_orchidee_init( MPI_COMM_ORCH,               &
            date0,    year,    month,          day,          &
            tmp_lon,  tmp_lat, soilth_lev,     control_flags )

       !- Initialize IOIPSL output files
       CALL intsurf_history(iim, jjm, tmp_lon, tmp_lat,  kindex, kjpindex, itau_sechiba, &
            date0_shifted, xrdt, control_flags, hist_id, &
            & hist2_id, hist_id_stom, hist_id_stom_IPCC)
       !
       IF ( ok_watchout ) THEN
          zlev_mean = zero
          DO ik=1, kjpindex

             zlev_mean = zlev_mean + zlev(ik)
          ENDDO
          ! Divide by one
          zlev_mean = zlev_mean / REAL(kjpindex,r_std)

          last_action_watch = itau_sechiba
          last_check_watch  = last_action_watch

          CALL watchout_init(iim, jjm, kjpindex, kjpindex, &
               & date0_shifted, last_action_watch, dt_watch, index_g, lon_g, lat_g, zlev_mean)
       ENDIF
       !
       IF ( check ) WRITE(numout,*) 'End of Initialisation of intersurf'
       !
    ENDIF
    !
    !  Shift the time step to phase the two models
    !
    itau_sechiba = kjit + itau_offset
    
    ! Update the calendar in xios by sending the new time step
    IF (l_first_intersurf) THEN
       ! Special case : the model is only in initialization phase and the current itau_sechiba is not a real time step. 
       ! Therefor give itau_sechiba+1 to xios to have a correct time axis in output files. 
       CALL xios_orchidee_update_calendar(itau_sechiba+1)
    ELSE
       CALL xios_orchidee_update_calendar(itau_sechiba)
    ENDIF

    CALL intsurf_time( itau_sechiba, date0_shifted, xrdt )
    !
    ! 1. gather input fields from kindex array
    !
    DO ik=1, kjpindex
       
       zu(ik)           = u(kindex(ik))
       zv(ik)           = v(kindex(ik))
       zzlev(ik)        = zlev(kindex(ik))
       zqair(ik)        = qair(kindex(ik))
       zprecip_rain(ik) = precip_rain(kindex(ik))*xrdt
       zprecip_snow(ik) = precip_snow(kindex(ik))*xrdt
       zlwdown(ik)      = lwdown(kindex(ik))
       zswnet(ik)       = swnet(kindex(ik))
       zswdown(ik)      = swdown(kindex(ik))
       zsinang(ik)      = sinang(kindex(ik))
       ztemp_air(ik)    = temp_air(kindex(ik))
       zepot_air(ik)    = epot_air(kindex(ik))
       zccanopy(ik)     = ccanopy(kindex(ik))
       zpetAcoef(ik)    = petAcoef(kindex(ik))
       zpeqAcoef(ik)    = peqAcoef(kindex(ik))
       zpetBcoef(ik)    = petBcoef(kindex(ik))
       zpeqBcoef(ik)    = peqBcoef(kindex(ik))
       zcdrag(ik)       = cdrag(kindex(ik))
       zpb(ik)          = pb(kindex(ik))
       
    ENDDO
    !
    ! 2. save the grid
    !
    IF ( check ) WRITE(numout,*) 'Save the grid'
    !
    IF (l_first_intersurf) THEN
       CALL histwrite_p(hist_id, 'LandPoints',  itau_sechiba+1, (/ REAL(kindex) /), kjpindex, kindex)
       CALL histwrite_p(hist_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       
       IF ( control_flags%ok_stomate ) THEN
          CALL histwrite_p(hist_id_stom, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
          CALL histwrite_p(hist_id_stom_IPCC, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       ENDIF
       
       CALL histwrite_p(hist_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)

       ! Syncronize output but only if flag ok_histsync is set to true
       IF (ok_histsync) THEN
          IF ( is_omp_root .AND. hist_id > 0 ) THEN
             CALL histsync(hist_id)
          END IF
       END IF
       

       CALL histwrite_p(hist2_id, 'LandPoints',  itau_sechiba+1, (/ REAL(kindex) /), kjpindex, kindex)
       CALL histwrite_p(hist2_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       CALL histwrite_p(hist2_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)
       ! Syncronize output but only if flag ok_histsync is set to true
       IF (ok_histsync) THEN
          IF ( is_omp_root .AND. hist2_id > 0 ) THEN
             CALL histsync(hist2_id)
          ENDIF
       END IF
    ENDIF
    !
    ! 3. call sechiba
    !
    IF ( check ) WRITE(numout,*) 'Calling sechiba'
    !
    CALL sechiba_main (itau_sechiba, iim*jjm, kjpindex, kindex, xrdt, date0_shifted, &
       & lrestart_read, lrestart_write, control_flags, &
       & lalo, contfrac, neighbours, resolution, &
! First level conditions
! Ajout Nathalie - Juin 2006 - passage q2m/t2m pour calcul rveget
!       & zzlev, zu, zv, zqair, ztemp_air, zepot_air, zccanopy, &
       & zzlev, zu, zv, zqair, zqair, ztemp_air, ztemp_air, zepot_air, zccanopy, &
! Variables for the implicit coupling
       & zcdrag, zpetAcoef, zpeqAcoef, zpetBcoef, zpeqBcoef, &
! Rain, snow, radiation and surface pressure
       & zprecip_rain ,zprecip_snow,  zlwdown, zswnet, zswdown, zsinang, zpb, &
! Output : Fluxes
       & zvevapp, zfluxsens, zfluxlat, zcoastal, zriver, znetco2, zcarblu, &
! Surface temperatures and surface properties
       & ztsol_rad, ztemp_sol_new, zqsurf, zalbedo, zemis, zz0, &
! File ids
       & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC ) 
    
    !
    IF ( check ) WRITE(numout,*) 'out of SECHIBA'
    !
    ! 4. save watchout
    !
    IF ( ok_watchout .AND. .NOT. l_first_intersurf ) THEN
       ! Accumulate last time step
       sum_zlev(:) = sum_zlev(:) + zzlev(:)
       sum_u(:) = sum_u(:) + zu(:)
       sum_v(:) = sum_v(:) + zv(:)
       sum_qair(:) = sum_qair(:) + zqair(:) 
       sum_temp_air(:) = sum_temp_air(:) + ztemp_air(:)
       sum_epot_air(:) = sum_epot_air(:) + zepot_air(:)
       sum_ccanopy(:) = sum_ccanopy(:) + zccanopy(:)
       sum_cdrag(:) = sum_cdrag(:) + zcdrag(:)
       sum_petAcoef(:) = sum_petAcoef(:) + zpetAcoef(:)
       sum_peqAcoef(:) = sum_peqAcoef(:) + zpeqAcoef(:)
       sum_petBcoef(:) = sum_petBcoef(:) + zpetBcoef(:)
       sum_peqBcoef(:) = sum_peqBcoef(:) + zpeqBcoef(:)
       sum_rain(:) = sum_rain(:) + zprecip_rain(:)
       sum_snow(:) = sum_snow(:) + zprecip_snow(:)
       sum_lwdown(:) = sum_lwdown(:) + zlwdown(:)
       sum_pb(:) = sum_pb(:) + zpb(:)
       
!!$       IF ( dt_watch > 3600 ) THEN
!!$          julian_watch = date0_shifted+((itau_sechiba-0.5)/dt_split_watch)*dt_watch/one_day
!!$          CALL solarang (julian_watch, julian0, iim, jjm, lon, lat, sinang)
!!$          WHERE ( sinang(:,:) .LT. EPSILON(un) ) 
!!$             isinang(:,:) = isinang(:,:) - 1
!!$          ENDWHERE
!!$          mean_sinang(:,:) = mean_sinang(:,:)+sinang(:,:)
!!$          !
!!$          DO ik=1,kjpindex          
!!$             j = ((kindex(ik)-1)/iim) + 1
!!$             i = (kindex(ik) - (j-1)*iim)
!!$             
!!$             sum_swnet(ik) = sum_swnet(ik) + sinang(i,j)*zswnet(ik)
!!$             sum_swdown(ik) = sum_swdown(ik) + sinang(i,j)*zswdown(ik)
!!$          ENDDO
!!$       ELSE
          sum_swnet(:) = sum_swnet(:) + zswnet(:)
          sum_swdown(:) = sum_swdown(:) + zswdown(:)
!!$       ENDIF
          
       do_watch = .FALSE.
       call isittime &
            &  (itau_sechiba,date0_shifted,xrdt,dt_watch,&
            &   last_action_watch,last_check_watch,do_watch)
       last_check_watch = itau_sechiba
       IF (do_watch) THEN
          !
          IF ( check ) WRITE(numout,*) 'save watchout'
          !
          IF (long_print) THEN
             WRITE(numout,*) "intersurf : write watchout for date ",date0,date0_shifted,itau_sechiba,&
                  & last_action_watch, last_check_watch
          ENDIF
          last_action_watch = itau_sechiba

          sum_zlev(:) = sum_zlev(:) / dt_split_watch
          sum_u(:) = sum_u(:) / dt_split_watch
          sum_v(:) = sum_v(:) / dt_split_watch
          sum_qair(:) = sum_qair(:) / dt_split_watch
          sum_temp_air(:) = sum_temp_air(:) / dt_split_watch
          sum_epot_air(:) = sum_epot_air(:) / dt_split_watch
          sum_ccanopy(:) = sum_ccanopy(:) / dt_split_watch
          sum_cdrag(:) = sum_cdrag(:) / dt_split_watch
          sum_petAcoef(:) = sum_petAcoef(:) / dt_split_watch
          sum_peqAcoef(:) = sum_peqAcoef(:) / dt_split_watch
          sum_petBcoef(:) = sum_petBcoef(:) / dt_split_watch
          sum_peqBcoef(:) = sum_peqBcoef(:) / dt_split_watch
          sum_rain(:) = sum_rain(:) / dt_split_watch
          sum_snow(:) = sum_snow(:) / dt_split_watch
          sum_lwdown(:) = sum_lwdown(:) / dt_split_watch
          sum_pb(:) = sum_pb(:) / dt_split_watch

!!$          IF ( dt_watch > 3600 ) THEN 
!!$             WHERE ( isinang(:,:) .GT. 0 )
!!$                mean_sinang(:,:) = mean_sinang(:,:) / isinang(:,:)
!!$             ENDWHERE
!!$             !
!!$             DO ik=1,kjpindex          
!!$                j = ((kindex(ik)-1)/iim) + 1
!!$                i = (kindex(ik) - (j-1)*iim)
!!$                IF (mean_sinang(i,j) > zero) THEN
!!$                   sum_swdown(ik) = sum_swdown(ik)/mean_sinang(i,j)
!!$                   sum_swnet(ik) =  sum_swnet(ik)/mean_sinang(i,j)
!!$                ELSE
!!$                   sum_swdown(ik) = zero
!!$                   sum_swnet(ik) =  zero
!!$                ENDIF
!!$             ENDDO
!!$          ELSE
             sum_swnet(:) = sum_swnet(:) / dt_split_watch
             sum_swdown(:) = sum_swdown(:) / dt_split_watch
!!$          ENDIF

          CALL watchout_write_p(kjpindex, itau_sechiba, xrdt, sum_zlev, sum_swdown, sum_rain, &
               &   sum_snow, sum_lwdown, sum_pb, sum_temp_air, sum_epot_air, sum_qair, &
               &   sum_u, sum_v, sum_swnet, sum_petAcoef, sum_peqAcoef, sum_petBcoef, sum_peqBcoef, &
               &   sum_cdrag, sum_ccanopy )
       ENDIF
    ENDIF
    !
    ! 5. scatter output fields
    !
    !
    z0(:)           = undef_sechiba
    coastalflow(:)  = undef_sechiba
    riverflow(:)    = undef_sechiba
    tsol_rad(:)     = undef_sechiba
    vevapp(:)       = undef_sechiba
    temp_sol_new(:) = undef_sechiba 
    qsurf(:)        = undef_sechiba 
    albedo(:,:)     = undef_sechiba
    fluxsens(:)     = undef_sechiba
    fluxlat(:)      = undef_sechiba
    emis(:)         = undef_sechiba 
    cdrag(:)        = undef_sechiba 
    !
    DO ik=1, kjpindex
       
       z0(kindex(ik))           = zz0(ik)
       coastalflow(kindex(ik))  = zcoastal(ik)
       riverflow(kindex(ik))    = zriver(ik)
       tsol_rad(kindex(ik))     = ztsol_rad(ik)
       vevapp(kindex(ik))       = zvevapp(ik)
       temp_sol_new(kindex(ik)) = ztemp_sol_new(ik)
       qsurf(kindex(ik))        = zqsurf(ik)
       albedo(kindex(ik),1)     = zalbedo(ik,1)
       albedo(kindex(ik),2)     = zalbedo(ik,2)
       fluxsens(kindex(ik))     = zfluxsens(ik)
       fluxlat(kindex(ik))      = zfluxlat(ik)
       emis(kindex(ik))         = zemis(ik)
       cdrag(kindex(ik))        = zcdrag(ik)

    ENDDO
    !
    ! Modified fields for variables scattered during the writing
    !
    dcoastal(:) = zcoastal(:)
    driver(:)   = zriver(:)
    !
    IF ( .NOT. l_first_intersurf) THEN
       
       CALL xios_orchidee_send_field("LandPoints" ,(/ ( REAL(ik), ik=1,kjpindex ) /))
       CALL xios_orchidee_send_field("Areas", area)
       CALL xios_orchidee_send_field("Contfrac",contfrac)
       CALL xios_orchidee_send_field("evap",zvevapp*one_day/dt_sechiba)
       CALL xios_orchidee_send_field("evap_alma",zvevapp/dt_sechiba)
       CALL xios_orchidee_send_field("coastalflow",dcoastal/dt_sechiba)
       CALL xios_orchidee_send_field("riverflow",driver/dt_sechiba)
       CALL xios_orchidee_send_field("temp_sol_C",ztemp_sol_new-ZeroCelsius)
       CALL xios_orchidee_send_field("temp_sol_K",ztemp_sol_new)
       CALL xios_orchidee_send_field("fluxsens",zfluxsens)
       CALL xios_orchidee_send_field("fluxlat",zfluxlat)
       CALL xios_orchidee_send_field("alb_vis",zalbedo(:,1))
       CALL xios_orchidee_send_field("alb_nir",zalbedo(:,2))
       CALL xios_orchidee_send_field("tair",ztemp_air)
       CALL xios_orchidee_send_field("qair",zqair)
       CALL xios_orchidee_send_field("q2m",zqair)
       CALL xios_orchidee_send_field("t2m",ztemp_air)
       CALL xios_orchidee_send_field("swnet",zswnet)
       CALL xios_orchidee_send_field("swdown",zswdown)

       IF ( .NOT. almaoutput ) THEN
          !
          !  scattered during the writing
          !
          CALL histwrite_p (hist_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'coastalflow',itau_sechiba, dcoastal, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'riverflow',itau_sechiba, driver, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'temp_sol', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'tsol_max', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'tsol_min', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'fluxlat',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'swnet',    itau_sechiba, zswnet, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'swdown',   itau_sechiba, zswdown, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
          CALL histwrite_p (hist_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
          CALL histwrite_p (hist_id, 'tair',     itau_sechiba, ztemp_air, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'qair',     itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'q2m',     itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist_id, 't2m',     itau_sechiba, ztemp_air, kjpindex, kindex)
       
          CALL histwrite_p (hist2_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'coastalflow',itau_sechiba, dcoastal, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'riverflow',itau_sechiba, driver, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'temp_sol', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tsol_max', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tsol_min', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'fluxlat',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'swnet',    itau_sechiba, zswnet, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'swdown',   itau_sechiba, zswdown, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'tair',     itau_sechiba, ztemp_air, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'qair',     itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'q2m',     itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 't2m',     itau_sechiba, ztemp_air, kjpindex, kindex)
       
       ELSE
          CALL histwrite_p (hist_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'SWnet',    itau_sechiba, zswnet, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'AvgSurfT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'RadT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Tair', itau_sechiba, ztemp_air, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Qair', itau_sechiba, zqair, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'SWnet',    itau_sechiba, zswnet, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'AvgSurfT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
          CALL histwrite_p (hist2_id, 'RadT', itau_sechiba, ztemp_sol_NEW, kjpindex, kindex)
       ENDIF
       
       ! Syncronize output but only if flag ok_histsync is set to true
       IF (ok_histsync .AND. is_omp_root ) THEN
          IF ( (dw .EQ. xrdt) .AND. hist_id > 0 ) THEN
             CALL histsync(hist_id)
          ENDIF
       END IF
   ENDIF
    !
    ! 6. Transform the water fluxes into Kg/m^2s and m^3/s
    !
    DO ik=1, kjpindex

       vevapp(kindex(ik)) = vevapp(kindex(ik))/xrdt
       coastalflow(kindex(ik)) = coastalflow(kindex(ik))/xrdt
       riverflow(kindex(ik)) = riverflow(kindex(ik))/xrdt

    ENDDO
    !
    IF ( lrestart_write .AND. ok_watchout ) THEN
       CALL watchout_close()
    ENDIF
    !
    IF(l_first_intersurf .AND. is_root_prc) CALL getin_dump
    l_first_intersurf = .FALSE.
    !
    IF (long_print) WRITE (numout,*) ' intersurf_main done '
    !
    CALL ipslnlf(new_number=old_fileout)
    !    
  END SUBROUTINE intersurf_main_1d


  SUBROUTINE Init_intersurf(nbp_l_lon,nbp_l_lat,kjpindex,kindex,orch_offset,orch_omp_size,orch_omp_rank,COMM)
    ! This subroutine is called from LMDZ before first call to ORCHIDEE intersurf_main

    USE mod_orchidee_para
    USE timer
    IMPLICIT NONE

    INTEGER,INTENT(IN)  :: nbp_l_lon
    INTEGER,INTENT(IN)  :: nbp_l_lat
    INTEGER,INTENT(IN)  :: kjpindex
    INTEGER,INTENT(IN)  :: kindex(:)
    INTEGER,INTENT(IN)  :: orch_offset
    INTEGER,INTENT(IN)  :: COMM
    INTEGER,INTENT(IN)  :: orch_omp_size
    INTEGER,INTENT(IN)  :: orch_omp_rank

    INTEGER,DIMENSION(kjpindex)  :: kindex_offset

    IF (orch_omp_rank==0) THEN
      CALL Init_timer
      CALL start_timer(timer_mpi)
      CALL set_grid_glo(nbp_l_lon,nbp_l_lat)
    ENDIF
    CALL barrier2_omp()    
    CALL init_orchidee_data_para(kjpindex,kindex,orch_offset,orch_omp_size,orch_omp_rank,COMM)
    CALL Set_stdout_file('out_orchidee')
    
    IF (is_omp_root) CALL Allocate_grid_glo
    CALL barrier2_omp()
    CALL init_ioipsl_para
          
    kindex_offset(:)=kindex(:)+offset
    CALL gather(kindex_offset,index_g)
    CALL bcast(index_g)  

    WRITE(numout,*) "kjpindex = ",kjpindex
    WRITE(numout,*) "offset for OMP = ",offset_omp
    WRITE(numout,*) "Index num local for continental points = ",kindex
    WRITE(numout,*) "Index num global for continental points = ",kindex_offset
    IF (is_omp_root) THEN
       WRITE(numout,*) "ROOT OMP, Index global MPI : ",kindex_mpi(:)
    ENDIF
    IF (is_root_prc) THEN
       WRITE(numout,*) "ROOT global, Index global : ",index_g(:)
    ENDIF
    
  END SUBROUTINE Init_intersurf
!
!
  SUBROUTINE intersurf_gathered (kjit, iim_glo, jjm_glo, kjpindex, kindex, xrdt, &
     & lrestart_read, lrestart_write, latlon, zcontfrac, zneighbours, zresolution, date0, &
! First level conditions
     & zlev,  u, v, qair, temp_air, epot_air, ccanopy, &
! Variables for the implicit coupling
     & cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
! Rain, snow, radiation and surface pressure
     & precip_rain, precip_snow, lwdown, swnet, swdown, pb, &
! Output : Fluxes
     & vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
! Surface temperatures and surface properties
     & tsol_rad, temp_sol_new, qsurf, albedo, emis, z0, lon_scat_g, lat_scat_g, q2m, t2m, &
! Add emission/deposit fields
     & field_out_names, fields_out, field_in_names, fields_in, &
! For VOC
     & sinang)  
    ! routines called : sechiba_main
    !
    USE mod_orchidee_para
    IMPLICIT NONE
    !    
    ! interface description for dummy arguments
    ! input scalar 
    INTEGER(i_std),INTENT (in)                            :: kjit          !! Time step number
    INTEGER(i_std),INTENT (in)                            :: iim_glo, jjm_glo  !! Dimension of global fields
    INTEGER(i_std),INTENT (in)                            :: kjpindex      !! Number of continental points
    REAL(r_std),INTENT (in)                               :: xrdt          !! Time step in seconds
    LOGICAL, INTENT (in)                                 :: lrestart_read !! Logical for _restart_ file to read
    LOGICAL, INTENT (in)                                 :: lrestart_write!! Logical for _restart_ file to write'
    REAL(r_std), INTENT (in)                              :: date0         !! Date at which kjit = 0
    ! input fields
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: kindex        !! Index for continental points
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: u             !! Lowest level wind speed
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: v             !! Lowest level wind speed 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: zlev          !! Height of first layer
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: qair          !! Lowest level specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: precip_rain   !! Rain precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: precip_snow   !! Snow precipitation
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: lwdown        !! Down-welling long-wave flux 
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: swnet         !! Net surface short-wave flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: swdown        !! Downwelling surface short-wave flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: temp_air      !! Air temperature in Kelvin
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: epot_air      !! Air potential energy
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: ccanopy       !! CO2 concentration in the canopy
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: petAcoef      !! Coeficients A from the PBL resolution
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: peqAcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: petBcoef      !! Coeficients B from the PBL resolution
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: peqBcoef      !! One for T and another for q
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: cdrag         !! Cdrag
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: pb            !! Surface pressure
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)        :: latlon        !! Geographical coordinates
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: zcontfrac     !! Fraction of continent
    INTEGER(i_std),DIMENSION (kjpindex,8), INTENT(in)    :: zneighbours   !! neighbours
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)        :: zresolution   !! size of the grid box
    REAL(r_std),DIMENSION (kjpindex), INTENT(in), OPTIONAL          :: q2m          !! Surface specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT(in), OPTIONAL          :: t2m          !! Surface air temperature
    ! output fields
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: z0            !! Surface roughness
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: coastalflow   !! Diffuse flow of water into the ocean (m^3/dt)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: riverflow     !! Largest rivers flowing into the ocean (m^3/dt)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: tsol_rad      !! Radiative surface temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: vevapp        !! Total of evaporation
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: temp_sol_new  !! New soil temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: qsurf         !! Surface specific humidity
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(out)       :: albedo        !! Albedo
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: fluxsens      !! Sensible chaleur flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: fluxlat       !! Latent chaleur flux
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)         :: emis          !! Emissivity


    !
    ! Optional arguments
    !
    ! Names and fields for emission variables : to be transport by GCM to chemistry model.
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN) :: field_out_names
    REAL(r_std),DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: fields_out
    !
    ! Names and fields for deposit variables : to be transport from chemistry model by GCM to ORCHIDEE.
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN) :: field_in_names
    REAL(r_std),DIMENSION(:,:), OPTIONAL, INTENT(IN) :: fields_in
    !
    ! For VOC 
    REAL(r_std), DIMENSION(kjpindex), OPTIONAL, INTENT(in) :: sinang  !! 


    ! LOCAL declaration
    ! work arrays to scatter and/or gather information just before/after sechiba_main call's
    ! and to keep output value for next call
    REAL(r_std),DIMENSION (kjpindex)                      :: zccanopy      !! Work array to keep ccanopy
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_rain  !! Work array to keep precip_rain
    REAL(r_std),DIMENSION (kjpindex)                      :: zprecip_snow  !! Work array to keep precip_snow
    REAL(r_std),DIMENSION (kjpindex)                      :: zz0           !! Work array to keep z0
    REAL(r_std),DIMENSION (kjpindex)                      :: zcdrag        !! Work array for surface drag
    REAL(r_std),DIMENSION (kjpindex)                      :: zcoastal      !! Work array to keep coastal flow
    REAL(r_std),DIMENSION (kjpindex)                      :: zriver        !! Work array to keep river out flow
    REAL(r_std),DIMENSION (kjpindex)                      :: dcoastal      !! Work array to keep coastal flow
    REAL(r_std),DIMENSION (kjpindex)                      :: driver        !! Work array to keep river out flow
    REAL(r_std),DIMENSION (kjpindex)                      :: znetco2       !! Work array to keep netco2flux
    REAL(r_std),DIMENSION (kjpindex)                      :: zcarblu       !! Work array to keep fco2_land_use
    REAL(r_std),DIMENSION (kjpindex)                      :: ztsol_rad     !! Work array to keep tsol_rad
    REAL(r_std),DIMENSION (kjpindex)                      :: zvevapp       !! Work array to keep vevapp
    REAL(r_std),DIMENSION (kjpindex)                      :: ztemp_sol_new !! Work array to keep temp_sol_new
    REAL(r_std),DIMENSION (kjpindex)                      :: zqsurf        !! Work array to keep qsurf
    REAL(r_std),DIMENSION (kjpindex,2)                    :: zalbedo       !! Work array to keep albedo
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxsens     !! Work array to keep fluxsens
    REAL(r_std),DIMENSION (kjpindex)                      :: zfluxlat      !! Work array to keep fluxlat
    REAL(r_std),DIMENSION (kjpindex)                      :: zemis         !! Work array to keep emis
    REAL(r_std),DIMENSION (kjpindex)                      :: zsinang       !! Work array to keep sinang
    REAL(r_std),ALLOCATABLE, DIMENSION (:)                :: soilth_lev    !! Vertical soil axis for thermal scheme (m)

    !
    ! Optional arguments
    !
    REAL(r_std),DIMENSION (iim_glo,jjm_glo), INTENT(IN) :: lon_scat_g, lat_scat_g !! The scattered values for longitude 
    !
    REAL(r_std),DIMENSION (:,:),ALLOCATABLE :: lon_scat, lat_scat !! The scattered values for longitude 
    !                                                             !! and latitude.
    !
    !
!    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dvevapp       !! Diagnostic array for evaporation
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dtemp_sol     !! for surface temperature
!$OMP THREADPRIVATE(dtemp_sol)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dfluxsens     !! for sensible heat flux
!$OMP THREADPRIVATE(dfluxsens)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dfluxlat      !! for latent heat flux
!$OMP THREADPRIVATE(dfluxlat)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dswnet        !! net solar radiation
!$OMP THREADPRIVATE(dswnet)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dswdown       !! Incident solar radiation
!$OMP THREADPRIVATE(dswdown)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:,:)                     :: dalbedo       !! albedo
!$OMP THREADPRIVATE(dalbedo)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dtair         !! air temperature
!$OMP THREADPRIVATE(dtair)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dqair         !! specific air humidity
!$OMP THREADPRIVATE(dqair)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dq2m          !! Surface specific humidity
!$OMP THREADPRIVATE(dq2m)
    REAL(r_std),ALLOCATABLE,SAVE,DIMENSION (:)                       :: dt2m          !! Surface air temperature
!$OMP THREADPRIVATE(dt2m)
    !                                                                          !! and latitude.
    !
    INTEGER(i_std)                                        :: i, j, ik
    INTEGER(i_std)                                        :: ier
    INTEGER(i_std)                                        :: itau_sechiba
    REAL(r_std)                                           :: mx, zlev_mean
    REAL(r_std), ALLOCATABLE, DIMENSION(:,:)              :: tmp_lon, tmp_lat, tmp_lev
    LOGICAL                                               :: do_watch      !! if it's time, write watchout
    INTEGER                                               :: old_fileout   !! old Logical Int for std IO output
    LOGICAL :: check = .FALSE.
    !
    LOGICAL, SAVE                                         :: fatmco2       !! Flag to force the value of atmospheric CO2 for vegetation.
!$OMP THREADPRIVATE(fatmco2)
    REAL(r_std), SAVE                                     :: atmco2        !! atmospheric CO2 
!$OMP THREADPRIVATE(atmco2)
    LOGICAL, SAVE                                     :: ok_q2m_t2m
!$OMP THREADPRIVATE(ok_q2m_t2m)

    !
    ! Number of fields to give (nb_fields_out) or get from (nb_fields_in) GCM :
    INTEGER(i_std), SAVE                                  :: nb_fields_out, nb_fields_in
!$OMP THREADPRIVATE(nb_fields_out, nb_fields_in)
    ! Id of fields to give (nb_fields_out) or get from (nb_fields_in) GCM :
    INTEGER(i_std)                                        :: i_fields_out, i_fields_in

    REAL,ALLOCATABLE,DIMENSION(:,:)                  :: lalo_mpi
    REAL,ALLOCATABLE,DIMENSION(:)                  :: zlev_glo1D
     
    !! Local landpoints vector
    REAL(r_std),DIMENSION (kjpindex)                       :: landpoints
    !
    CALL ipslnlf_p(old_number=old_fileout)
    CALL ipslnlf_p(new_number=numout)
    

    !
    IF (l_first_intersurf) THEN
       ! Test sur la presence de q2m et t2m
       IF (PRESENT(q2m) .AND. PRESENT(t2m)) THEN
         ok_q2m_t2m=.TRUE.
       ELSE
         ok_q2m_t2m=.FALSE.
       ENDIF

       CALL intsurf_time( kjit, date0, xrdt )

       IF ( check ) WRITE(numout,*) 'Initialisation of intersurf'
       WRITE(numout,*) 'Initialisation of intersurf gathered'       !

       CALL ioget_calendar (one_year, one_day)

       IF (is_omp_root) THEN
          ALLOCATE(lon_scat(iim_g,jj_nb))
          ALLOCATE(lat_scat(iim_g,jj_nb))
       ELSE
          ALLOCATE(lon_scat(1,1))
          ALLOCATE(lat_scat(1,1))
       ENDIF

       ALLOCATE(dtemp_sol(kjpindex))
       ALLOCATE(dfluxsens(kjpindex))
       ALLOCATE(dfluxlat(kjpindex))
       ALLOCATE(dswnet(kjpindex))
       ALLOCATE(dswdown(kjpindex))
       ALLOCATE(dalbedo(kjpindex,2))
       ALLOCATE(dtair(kjpindex))
       ALLOCATE(dqair(kjpindex)) 

       IF (ok_q2m_t2m) THEN
          ALLOCATE(dq2m(kjpindex))
          ALLOCATE(dt2m(kjpindex))
       ENDIF

       CALL init_WriteField_p
       !
       ! Allocation of grid variables
       !
       CALL init_grid ( kjpindex )
       !
       !  Create the internal coordinate table
       !
       lalo(:,:) = latlon(:,:)
       CALL gather(lalo,lalo_g)
       !
       !-
       !- Store variable to help describe the grid
       !- once the points are gathered.
       !-
       neighbours(:,:) = zneighbours(:,:)
       CALL gather(neighbours,neighbours_g)
       !
       resolution(:,:) = zresolution(:,:)
       CALL gather(resolution,resolution_g)
       !
       area(:) = resolution(:,1)*resolution(:,2)
       CALL gather(area,area_g)
       !
       !- Store the fraction of the continents only once so that the user
       !- does not change them afterwards.
       !
       contfrac(:) = zcontfrac(:)
       CALL gather(contfrac,contfrac_g)
       !
       !
       !  Create the internal coordinate table
       !
       IF ( (.NOT.ALLOCATED(tmp_lon))) THEN
          ALLOCATE(tmp_lon(iim_g,jj_nb))
       ENDIF
       IF ( (.NOT.ALLOCATED(tmp_lat))) THEN
          ALLOCATE(tmp_lat(iim_g,jj_nb))
       ENDIF
       IF ( (.NOT.ALLOCATED(tmp_lev))) THEN
          ALLOCATE(tmp_lev(iim_g,jj_nb))
       ENDIF
       !
       !  Either we have the scattered coordinates as arguments or
       !  we have to do the work here.
       !
       IF (is_omp_root) THEN
          lon_scat(:,:)=zero
          lat_scat(:,:)=zero 
          CALL scatter2D_mpi(lon_scat_g,lon_scat)
          CALL scatter2D_mpi(lat_scat_g,lat_scat)
          lon_scat(:,1)=lon_scat(:,2)
          lon_scat(:,jj_nb)=lon_scat(:,2)
          lat_scat(:,1)=lat_scat(iim_g,1)
          lat_scat(:,jj_nb)=lat_scat(1,jj_nb)
          
          tmp_lon(:,:) = lon_scat(:,:)
          tmp_lat(:,:) = lat_scat(:,:)
          
          IF (is_mpi_root) THEN
             lon_g(:,:) = lon_scat_g(:,:)
             lat_g(:,:) = lat_scat_g(:,:)
          ENDIF
       ENDIF
       !

!ANNE ici possibilite de calculer la grille du modele ? (lignes 965 ??
       

!ANNE calcul zlev, conserve la version de Martial
       IF (is_root_prc) ALLOCATE(zlev_glo1D(nbp_glo))
       CALL gather(zlev,zlev_glo1D)
       
       IF (is_root_prc) THEN
          DO ik=1, nbp_glo
             j = INT( (index_g(ik)-1) / iim_g ) + 1
             i = index_g(ik) - (j-1) * iim_g
             zlev_g(i,j) = zlev_glo1D(ik)
          ENDDO
       ENDIF  
       !
       !
       !  Configuration of SSL specific parameters
       !
       CALL intsurf_config(control_flags,xrdt)
       !
       !Config  Key  = FORCE_CO2_VEG
       !Config  Desc = Flag to force the value of atmospheric CO2 for vegetation.
       !Config  Def  = FALSE
       !Config  Help = If this flag is set to true, the ATM_CO2 parameter is used
       !Config         to prescribe the atmospheric CO2.
       !Config         This Flag is only use in couple mode.
       !
       fatmco2=.FALSE.
       CALL getin_p('FORCE_CO2_VEG',fatmco2)
       !
       ! Next flag is only use in couple mode with a gcm in intersurf.
       ! In forced mode, it has already been read and set in driver.
       IF ( fatmco2 ) THEN
          !Config  Key  = ATM_CO2
          !Config  IF   = FORCE_CO2_VEG (in not forced mode)
          !Config  Desc = Value for atm CO2 
          !Config  Def  = 350.
          !Config  Help = Value to prescribe the atm CO2.
          !Config         For pre-industrial simulations, the value is 286.2 .
          !Config         348. for 1990 year.
          !
          atmco2=350.
          CALL getin_p('ATM_CO2',atmco2)
          WRITE(numout,*) 'atmco2 ',atmco2
       ENDIF
       
       !
       CALL intsurf_restart(kjit, date0, xrdt, control_flags, rest_id, rest_id_stom, itau_offset)
       itau_sechiba = kjit + itau_offset
       !

       !!- Initialize module for output with XIOS
       !
       ! Get the vertical soil levels for the thermal scheme, to be used in xios_orchidee_init
       ALLOCATE(soilth_lev(ngrnd), stat=ier)
       IF (ier /= 0) THEN
          CALL ipslerr_p(3,'intersurf_gathered', 'Error in allocation of soilth_lev','','')
       END IF
       soilth_lev(1:ngrnd) = thermosoil_levels()
       
       CALL xios_orchidee_init( MPI_COMM_ORCH,               &
            date0,    year,    month,          day,          &
            tmp_lon,  tmp_lat, soilth_lev,     control_flags )
       
       !- Initialize IOIPSL output files
       CALL intsurf_history(iim_g, jj_nb, tmp_lon, tmp_lat,  kindex, kjpindex, itau_sechiba, &
 &                          date0_shifted, xrdt, control_flags, hist_id, hist2_id, hist_id_stom, hist_id_stom_IPCC)
       !

       CALL bcast_omp(hist_id)
       CALL bcast_omp(hist2_id)
       CALL bcast_omp(hist_id_stom)
       CALL bcast_omp(hist_id_stom_IPCC)
       !
       IF ( ok_watchout ) THEN
          IF (is_root_prc) THEN
             zlev_mean = zero
             DO ik=1, nbp_glo
                j = ((index_g(ik)-1)/iim_g) + 1
                i = (index_g(ik) - (j-1)*iim_g)
                
                zlev_mean = zlev_mean + zlev_g(i,j)
             ENDDO
             zlev_mean = zlev_mean / REAL(nbp_glo,r_std)
          ENDIF

          last_action_watch = itau_sechiba
          last_check_watch =  last_action_watch

          ! Only root proc write watchout file
          IF (is_root_prc) &
               CALL watchout_init(iim_g, jjm_g, kjpindex, nbp_glo, &
               & date0_shifted, last_action_watch, dt_watch, index_g, lon_g, lat_g, zlev_mean)
       ENDIF

       ! Prepare fieds out/in for interface with GCM.
       IF (PRESENT(field_out_names)) THEN
          nb_fields_out=SIZE(field_out_names)
       ELSE
          nb_fields_out=0
       ENDIF
       IF (PRESENT(field_in_names)) THEN
          nb_fields_in=SIZE(field_in_names)
       ELSE
          nb_fields_in=0
       ENDIF

!>> VOC in coupled mode 
       IF ( PRESENT(sinang) )  THEN 
          zsinang(:) = sinang(:)
       ELSE
          zsinang(:) = zero
       ENDIF

       !
       IF ( check ) WRITE(numout,*) 'End of Initialisation of intersurf'
       !
    ENDIF

    !
    !! Change to be in the orchidee context for XIOS
    !
    CALL xios_orchidee_change_context("orchidee")
    
    !
    !  Shift the time step to phase the two models
    !
    itau_sechiba = kjit + itau_offset
    
    ! Update the calendar in xios by sending the new time step
    IF (l_first_intersurf) THEN
       ! Special case : the model is only in initialization phase and the current itau_sechiba is not a real time step. 
       ! Therefor give itau_sechiba+1 to xios to have a correct time axis in output files. 
       CALL xios_orchidee_update_calendar(itau_sechiba+1)
    ELSE
       CALL xios_orchidee_update_calendar(itau_sechiba)
    ENDIF

    CALL intsurf_time( itau_sechiba, date0_shifted, xrdt )

    ! Test LMDz inputs
    IF (check_INPUTS) THEN
       IF ( sec <= 7200 ) THEN
          WRITE(numout,*) "sec = ",sec
          CALL WriteFieldI_p("Wu",u)
          WRITE(numout,*) "Lowest level wind speed :",u
          CALL WriteFieldI_p("Wv",v)
          WRITE(numout,*) "Lowest level wind speed  :",v
          CALL WriteFieldI_p("Wzlev",zlev)
          WRITE(numout,*) "Height of first layer :",zlev
          CALL WriteFieldI_p("Wqair",qair)
          WRITE(numout,*) "Lowest level specific humidity :",qair
          CALL WriteFieldI_p("Wprecip_rain",precip_rain)
          WRITE(numout,*) "Rain precipitation :",precip_rain
          CALL WriteFieldI_p("Wprecip_snow",precip_snow)
          WRITE(numout,*) "Snow precipitation :",precip_snow
          CALL WriteFieldI_p("Wlwdown",lwdown)
          WRITE(numout,*) "Down-welling long-wave flux  :",lwdown
          CALL WriteFieldI_p("Wswnet",swnet)
          WRITE(numout,*) "Net surface short-wave flux :",swnet
          CALL WriteFieldI_p("Wswdown",swdown)
          WRITE(numout,*) "Downwelling surface short-wave flux :",swdown
          CALL WriteFieldI_p("Wtemp_air",temp_air)
          WRITE(numout,*) "Air temperature in Kelvin :",temp_air
          CALL WriteFieldI_p("Wepot_air",epot_air)
          WRITE(numout,*) "Air potential energy :",epot_air
          CALL WriteFieldI_p("Wccanopy",ccanopy)
          WRITE(numout,*) "CO2 concentration in the canopy :",ccanopy
          CALL WriteFieldI_p("WpetAcoef",petAcoef)
          WRITE(numout,*) "Coeficients A from the PBL resolution :",petAcoef
          CALL WriteFieldI_p("WpeqAcoef",peqAcoef)
          WRITE(numout,*) "One for T and another for q :",peqAcoef
          CALL WriteFieldI_p("WpetBcoef",petBcoef)
          WRITE(numout,*) "Coeficients B from the PBL resolution :",petBcoef
          CALL WriteFieldI_p("WpeqBcoef",peqBcoef)
          WRITE(numout,*) "One for T and another for q :",peqBcoef
          CALL WriteFieldI_p("Wcdrag",cdrag)
          WRITE(numout,*) "Cdrag :",cdrag
          CALL WriteFieldI_p("Wpb",pb)
          WRITE(numout,*) "Surface pressure :",pb
          CALL WriteFieldI_p("Wzcontfrac",zcontfrac)
          WRITE(numout,*) "Fraction of continent :",zcontfrac
          IF ( ok_q2m_t2m) THEN
             CALL WriteFieldI_p("Wq2m",q2m)
             WRITE(numout,*) "Surface specific humidity :",q2m
             CALL WriteFieldI_p("Wt2m",t2m)
             WRITE(numout,*) "Surface air temperature :",t2m
          ENDIF
       ENDIF
    ENDIF

    !
    ! 1. Just change the units of some input fields
    !
    DO ik=1, kjpindex
       
       zprecip_rain(ik) = precip_rain(ik)*xrdt
       zprecip_snow(ik) = precip_snow(ik)*xrdt
       zcdrag(ik)       = cdrag(ik)
       
    ENDDO
    !
    IF (check_INPUTS) THEN
       WRITE(numout,*) "Intersurf_main_gathered :"
       WRITE(numout,*) "Time step number = ",kjit
       WRITE(numout,*) "Dimension of input fields for local mpi process = ",iim_g, jj_nb
       WRITE(numout,*) "Number of continental points = ",kjpindex
       WRITE(numout,*) "Time step in seconds = ",xrdt
       WRITE(numout,*) "Logical for _restart_ file to read, write = ",lrestart_read,lrestart_write
       WRITE(numout,*) "Date at which kjit = 0  =  ",date0
       WRITE(numout,*) "offset for OMP = ",offset_omp
       WRITE(numout,*) "Index for continental points = ",kindex
       IF (is_omp_root) THEN
          WRITE(numout,*) "ROOT OMP, Index global MPI : ",kindex_mpi(:)
       ENDIF
       IF (is_root_prc) THEN
          WRITE(numout,*) "ROOT global, Index global : ",index_g(:)
       ENDIF
       WRITE(numout,*) "Lowest level wind speed North = ",u
       WRITE(numout,*) "Lowest level wind speed East = ",v
       WRITE(numout,*) "Height of first layer = ",zlev
       WRITE(numout,*) "Lowest level specific humidity = ",qair
       WRITE(numout,*) "Rain precipitation = ",zprecip_rain
       WRITE(numout,*) "Snow precipitation = ",zprecip_snow
       WRITE(numout,*) "Down-welling long-wave flux = ",lwdown
       WRITE(numout,*) "Net surface short-wave flux = ",swnet
       WRITE(numout,*) "Downwelling surface short-wave flux = ",swdown
       WRITE(numout,*) "Air temperature in Kelvin = ",temp_air
       WRITE(numout,*) "Air potential energy = ",epot_air
       WRITE(numout,*) "CO2 concentration in the canopy = ",ccanopy
       WRITE(numout,*) "Coeficients A from the PBL resolution = ",petAcoef
       WRITE(numout,*) "One for T and another for q = ",peqAcoef
       WRITE(numout,*) "Coeficients B from the PBL resolution = ",petBcoef
       WRITE(numout,*) "One for T and another for q = ",peqBcoef
       WRITE(numout,*) "Cdrag = ",zcdrag
       WRITE(numout,*) "Surface pressure = ",pb
! Variables utilisees uniquement dans l_first et donc plus definies ensuite
!       IF (is_omp_root) THEN
!          WRITE(numout,*) "Geographical coordinates lon = "
!          WRITE(numout,*) lon_scat
!          WRITE(numout,*) "Geographical coordinates lat = "
!          WRITE(numout,*) lat_scat 
!       ENDIF
       WRITE(numout,*) "Geographical land coordinates lon = ", lalo(:,2)
       WRITE(numout,*) "Geographical land coordinates lat = ", lalo(:,1)
       WRITE(numout,*) "Fraction of continent in the grid = ",zcontfrac
    ENDIF
    ! Fields for deposit variables : to be transport from chemistry model by GCM to ORCHIDEE.
    DO i_fields_in=1,nb_fields_in
       WRITE(numout,*) i_fields_in," Champ = ",TRIM(field_in_names(i_fields_in)) 
       SELECT CASE(TRIM(field_in_names(i_fields_in)))
       CASE DEFAULT 
          CALL ipslerr_p (3,'intsurf_gathered', &
            &          'You ask in GCM an unknown field '//TRIM(field_in_names(i_fields_in))//&
            &          ' to give to ORCHIDEE for this specific version.',&
            &          'This model won''t be able to continue.', &
            &          '(check your tracer parameters in GCM)')
       END SELECT
    ENDDO


    !
    ! 2. modification of co2
    !
    IF ( fatmco2 ) THEN
       zccanopy(:) = atmco2
       WRITE (numout,*) 'Modification of the ccanopy value. CO2 = ',atmco2
    ELSE
       zccanopy(:) = ccanopy(:)
    ENDIF
    !
    ! 3. save the grid
    !
    IF ( check ) WRITE(numout,*) 'Save the grid'
    !
    IF (l_first_intersurf) THEN
       landpoints(:)=(/ ( REAL(ik), ik=1,kjpindex ) /)
       CALL histwrite_p(hist_id, 'LandPoints',  itau_sechiba+1, landpoints, kjpindex, kindex)
       CALL histwrite_p(hist_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       IF ( control_flags%ok_stomate ) THEN
          CALL histwrite_p(hist_id_stom, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
          IF ( hist_id_stom_ipcc > 0 ) &
               CALL histwrite_p(hist_id_stom_IPCC, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
       ENDIF
       CALL histwrite_p(hist_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)

       ! Syncronize output but only if flag ok_histsync is set to true       
       IF (ok_histsync) THEN
          IF (is_omp_root .AND. hist_id > 0) THEN
             CALL histsync(hist_id)
          END IF
       END IF
       
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'LandPoints',  itau_sechiba+1, landpoints, kjpindex, kindex)
          CALL histwrite_p(hist2_id, 'Areas',  itau_sechiba+1, area, kjpindex, kindex)
          CALL histwrite_p(hist2_id, 'Contfrac',  itau_sechiba+1, contfrac, kjpindex, kindex)

          ! Syncronize output but only if flag ok_histsync is set to true
          IF (ok_histsync .AND. is_omp_root) THEN
             CALL histsync(hist2_id)
          ENDIF
       ENDIF
       !
    ENDIF
    !
    ! 4. call sechiba for continental points only
    !
    IF ( check ) WRITE(numout,*) 'Calling sechiba'
    !
    IF (.NOT. ok_q2m_t2m) THEN
       CALL sechiba_main (itau_sechiba, iim_g*jj_nb, kjpindex, kindex, xrdt, date0_shifted, &
            & lrestart_read, lrestart_write, control_flags, &
            & lalo, contfrac, neighbours, resolution, &
! First level conditions
! Ajout Nathalie - Juin 2006 - passage q2m/t2m pour calcul rveget
!       & zlev, u, v, qair, temp_air, epot_air, ccanopy, &
            & zlev, u, v, qair, qair, temp_air, temp_air, epot_air, zccanopy, &
! Variables for the implicit coupling
            & zcdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
! Rain, snow, radiation and surface pressure
            & zprecip_rain ,zprecip_snow,  lwdown, swnet, swdown, zsinang, pb, &
! Output : Fluxes
            & zvevapp, zfluxsens, zfluxlat, zcoastal, zriver, znetco2, zcarblu, &
! Surface temperatures and surface properties
            & ztsol_rad, ztemp_sol_new, zqsurf, zalbedo, zemis, zz0, &
! File ids
            & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC ) 
    ELSE
       CALL sechiba_main (itau_sechiba, iim_g*jj_nb, kjpindex, kindex, xrdt, date0_shifted, &
            & lrestart_read, lrestart_write, control_flags, &
            & lalo, contfrac, neighbours, resolution, &
! First level conditions
! Ajout Nathalie - Juin 2006 - passage q2m/t2m pour calcul rveget
!       & zlev, u, v, qair, temp_air, epot_air, ccanopy, &
            & zlev, u, v, qair, q2m, t2m, temp_air, epot_air, zccanopy, &
! Variables for the implicit coupling
            & zcdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
! Rain, snow, radiation and surface pressure
            & zprecip_rain ,zprecip_snow,  lwdown, swnet, swdown, zsinang, pb, &
! Output : Fluxes
            & zvevapp, zfluxsens, zfluxlat, zcoastal, zriver, znetco2, zcarblu, &
! Surface temperatures and surface properties
            & ztsol_rad, ztemp_sol_new, zqsurf, zalbedo, zemis, zz0, &
! File ids
            & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC ) 
    ENDIF
    !
    IF ( check ) WRITE(numout,*) 'out of SECHIBA'
    !
    ! 5. save watchout 
    !
    IF ( ok_watchout .AND. .NOT. l_first_intersurf ) THEN
       ! Accumulate last time step
       sum_zlev(:) = sum_zlev(:) + zlev(:)
       sum_u(:) = sum_u(:) + u(:)
       sum_v(:) = sum_v(:) + v(:)
       sum_qair(:) = sum_qair(:) + qair(:) 
       sum_temp_air(:) = sum_temp_air(:) + temp_air(:)
       sum_epot_air(:) = sum_epot_air(:) + epot_air(:)
       sum_ccanopy(:) = sum_ccanopy(:) + ccanopy(:)
       sum_cdrag(:) = sum_cdrag(:) + zcdrag(:)
       sum_petAcoef(:) = sum_petAcoef(:) + petAcoef(:)
       sum_peqAcoef(:) = sum_peqAcoef(:) + peqAcoef(:)
       sum_petBcoef(:) = sum_petBcoef(:) + petBcoef(:)
       sum_peqBcoef(:) = sum_peqBcoef(:) + peqBcoef(:)
       sum_rain(:) = sum_rain(:) + zprecip_rain(:)
       sum_snow(:) = sum_snow(:) + zprecip_snow(:)
       sum_lwdown(:) = sum_lwdown(:) + lwdown(:)
       sum_pb(:) = sum_pb(:) + pb(:)

!!$       IF ( dt_watch > 3600 ) THEN 
!!$          julian_watch = date0_shifted+((itau_sechiba-0.5)/dt_split_watch)*dt_watch/one_day
!!$          CALL solarang (julian_watch, julian0, iim_g, jj_nb, tmp_lon, tmp_lat, sinang)
!!$          WHERE ( sinang(:,:) .LT. EPSILON(un) )
!!$             isinang(:,:) = isinang(:,:) - 1
!!$          ENDWHERE
!!$          mean_sinang(:,:) = mean_sinang(:,:)+sinang(:,:)
!!$          !
!!$          DO ik=1,kjpindex          
!!$             j = ((kindex(ik)-1)/iim_g) + 1
!!$             i = (kindex(ik) - (j-1)*iim_g)
!!$             
!!$             sum_swnet(ik) = sum_swnet(ik) + sinang(i,j)*swnet(ik)
!!$             sum_swdown(ik) = sum_swdown(ik) + sinang(i,j)*swdown(ik)
!!$          ENDDO
!!$       ELSE
          sum_swnet(:) = sum_swnet(:) + swnet(:)
          sum_swdown(:) = sum_swdown(:) + swdown(:)
!!$       ENDIF
          
       do_watch = .FALSE.
       call isittime &
            &  (itau_sechiba,date0_shifted,xrdt,dt_watch,&
            &   last_action_watch,last_check_watch,do_watch)
       last_check_watch = itau_sechiba
       IF (do_watch) THEN
          !
          IF ( check ) WRITE(numout,*) 'save watchout'
          !
          IF (long_print) THEN
             WRITE(numout,*) "intersurf : write watchout for date ",date0,date0_shifted,itau_sechiba, &
                  & last_action_watch,last_check_watch
          ENDIF
          last_action_watch = itau_sechiba

          sum_zlev(:) = sum_zlev(:) / dt_split_watch
          sum_u(:) = sum_u(:) / dt_split_watch
          sum_v(:) = sum_v(:) / dt_split_watch
          sum_qair(:) = sum_qair(:) / dt_split_watch
          sum_temp_air(:) = sum_temp_air(:) / dt_split_watch
          sum_epot_air(:) = sum_epot_air(:) / dt_split_watch
          sum_ccanopy(:) = sum_ccanopy(:) / dt_split_watch
          sum_cdrag(:) = sum_cdrag(:) / dt_split_watch
          sum_petAcoef(:) = sum_petAcoef(:) / dt_split_watch
          sum_peqAcoef(:) = sum_peqAcoef(:) / dt_split_watch
          sum_petBcoef(:) = sum_petBcoef(:) / dt_split_watch
          sum_peqBcoef(:) = sum_peqBcoef(:) / dt_split_watch
          sum_rain(:) = sum_rain(:) / dt_split_watch
          sum_snow(:) = sum_snow(:) / dt_split_watch
          sum_lwdown(:) = sum_lwdown(:) / dt_split_watch
          sum_pb(:) = sum_pb(:) / dt_split_watch

!!$          IF ( dt_watch > 3600 ) THEN 
!!$             WHERE ( isinang(:,:) .GT. 0 )
!!$                mean_sinang(:,:) = mean_sinang(:,:) / isinang(:,:)
!!$             ENDWHERE
!!$             !
!!$             DO ik=1,kjpindex          
!!$                j = ((kindex(ik)-1)/iim_g) + 1
!!$                i = (kindex(ik) - (j-1)*iim_g)
!!$                IF (mean_sinang(i,j) > zero) THEN
!!$                   sum_swdown(ik) = sum_swdown(ik)/mean_sinang(i,j)
!!$                   sum_swnet(ik) =  sum_swnet(ik)/mean_sinang(i,j)
!!$                ELSE
!!$                   sum_swdown(ik) = zero
!!$                   sum_swnet(ik) =  zero
!!$                ENDIF
!!$             ENDDO
!!$          ELSE
             sum_swnet(:) = sum_swnet(:) / dt_split_watch
             sum_swdown(:) = sum_swdown(:) / dt_split_watch
!!$          ENDIF

          CALL watchout_write_p(kjpindex, itau_sechiba, xrdt, sum_zlev, sum_swdown, sum_rain, &
               &   sum_snow, sum_lwdown, sum_pb, sum_temp_air, sum_epot_air, sum_qair, &
               &   sum_u, sum_v, sum_swnet, sum_petAcoef, sum_peqAcoef, sum_petBcoef, sum_peqBcoef, &
               &   sum_cdrag, sum_ccanopy )
       ENDIF       
    ENDIF
    !
    ! 6. scatter output fields
    !
    z0(:)           = undef_sechiba
    coastalflow(:)  = undef_sechiba
    riverflow(:)    = undef_sechiba
    tsol_rad(:)     = undef_sechiba
    vevapp(:)       = undef_sechiba
    temp_sol_new(:) = undef_sechiba
    qsurf(:)        = undef_sechiba
    albedo(:,1)     = undef_sechiba
    albedo(:,2)     = undef_sechiba
    fluxsens(:)     = undef_sechiba
    fluxlat(:)      = undef_sechiba
    emis(:)         = undef_sechiba
    cdrag(:)        = undef_sechiba
    !    

    DO ik=1, kjpindex
       
       z0(ik)           = zz0(ik)
       coastalflow(ik)  = zcoastal(ik)
       riverflow(ik)    = zriver(ik)
       tsol_rad(ik)     = ztsol_rad(ik)
       vevapp(ik)       = zvevapp(ik)
       temp_sol_new(ik) = ztemp_sol_new(ik)
       qsurf(ik)        = zqsurf(ik)
       albedo(ik,1)     = zalbedo(ik,1)
       albedo(ik,2)     = zalbedo(ik,2)
       fluxsens(ik)     = zfluxsens(ik)
       fluxlat(ik)      = zfluxlat(ik)
       emis(ik)         = zemis(ik)
       cdrag(ik)        = zcdrag(ik)
       
      

       
    ENDDO

    !
    ! Modified fields for variables scattered during the writing
    !
    dcoastal(:) = zcoastal(:)
    driver(:)   = zriver(:)
    !
    IF ( .NOT. l_first_intersurf) THEN
       
       CALL xios_orchidee_send_field("LandPoints" ,(/ ( REAL(ik), ik=1,kjpindex ) /))
       CALL xios_orchidee_send_field("Areas", area)
       CALL xios_orchidee_send_field("Contfrac",contfrac)
       CALL xios_orchidee_send_field("evap",zvevapp*one_day/dt_sechiba)
       CALL xios_orchidee_send_field("evap_alma",zvevapp/dt_sechiba)
       CALL xios_orchidee_send_field("coastalflow",dcoastal/dt_sechiba)
       CALL xios_orchidee_send_field("riverflow",driver/dt_sechiba)
       CALL xios_orchidee_send_field("temp_sol_C",ztemp_sol_new-ZeroCelsius)
       CALL xios_orchidee_send_field("temp_sol_K",ztemp_sol_new)
       CALL xios_orchidee_send_field("fluxsens",zfluxsens)
       CALL xios_orchidee_send_field("fluxlat",zfluxlat)
       CALL xios_orchidee_send_field("alb_vis",zalbedo(:,1))
       CALL xios_orchidee_send_field("alb_nir",zalbedo(:,2))
       CALL xios_orchidee_send_field("tair",temp_air)
       CALL xios_orchidee_send_field("qair",qair)
       CALL xios_orchidee_send_field("swnet",swnet)
       CALL xios_orchidee_send_field("swdown",swdown)
          
       IF (ok_q2m_t2m) THEN
          CALL xios_orchidee_send_field("t2m",t2m)
          CALL xios_orchidee_send_field("q2m",q2m)
       ELSE
          CALL xios_orchidee_send_field("t2m",temp_air)
          CALL xios_orchidee_send_field("q2m",qair)
       ENDIF

       IF ( .NOT. almaoutput ) THEN
          !
          !  scattered during the writing
          !           
          CALL histwrite_p (hist_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'coastalflow',itau_sechiba, dcoastal, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'riverflow',itau_sechiba, driver, kjpindex, kindex)
          ! 
          CALL histwrite_p (hist_id, 'temp_sol', itau_sechiba, ztemp_sol_new, kjpindex, kindex)

          CALL histwrite_p (hist_id, 'tsol_max', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'tsol_min', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'fluxlat',  itau_sechiba, zfluxlat,  kjpindex, kindex)
          CALL histwrite_p (hist_id, 'swnet',    itau_sechiba, swnet,    kjpindex, kindex)
          CALL histwrite_p (hist_id, 'swdown',   itau_sechiba, swdown,   kjpindex, kindex)
          CALL histwrite_p (hist_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
          CALL histwrite_p (hist_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
          CALL histwrite_p (hist_id, 'tair',     itau_sechiba, temp_air, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'qair',     itau_sechiba, qair, kjpindex, kindex)
          IF (ok_q2m_t2m) THEN
             CALL histwrite_p (hist_id, 't2m',      itau_sechiba, t2m, kjpindex, kindex)
             CALL histwrite_p (hist_id, 'q2m',      itau_sechiba, q2m, kjpindex, kindex)
          ELSE
            CALL histwrite_p (hist_id, 't2m',      itau_sechiba, temp_air, kjpindex, kindex)
            CALL histwrite_p (hist_id, 'q2m',      itau_sechiba, qair, kjpindex, kindex)
          ENDIF
          !
          IF ( hist2_id > 0 ) THEN
             CALL histwrite_p (hist2_id, 'evap',     itau_sechiba, zvevapp, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'coastalflow',itau_sechiba, dcoastal, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'riverflow',itau_sechiba, driver, kjpindex, kindex)
             ! 
             CALL histwrite_p (hist2_id, 'temp_sol', itau_sechiba, temp_sol_new, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'tsol_max', itau_sechiba, temp_sol_new, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'tsol_min', itau_sechiba, temp_sol_new, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'fluxsens', itau_sechiba, zfluxsens, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'fluxlat',  itau_sechiba, zfluxlat,  kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'swnet',    itau_sechiba, swnet,    kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'swdown',   itau_sechiba, swdown,   kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'alb_vis',  itau_sechiba, zalbedo(:,1), kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'alb_nir',  itau_sechiba, zalbedo(:,2), kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'tair',     itau_sechiba, temp_air, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'qair',     itau_sechiba, qair, kjpindex, kindex)
             IF (ok_q2m_t2m) THEN
                CALL histwrite_p (hist2_id, 't2m',      itau_sechiba, t2m, kjpindex, kindex)
                CALL histwrite_p (hist2_id, 'q2m',      itau_sechiba, q2m, kjpindex, kindex)
             ELSE
               CALL histwrite_p (hist2_id, 't2m',      itau_sechiba, temp_air, kjpindex, kindex)
               CALL histwrite_p (hist2_id, 'q2m',      itau_sechiba, qair, kjpindex, kindex)
            ENDIF
          ENDIF
       ELSE
          CALL histwrite_p (hist_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'SWnet',    itau_sechiba, swnet, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'AvgSurfT', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
          CALL histwrite_p (hist_id, 'RadT', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
          !
          IF ( hist2_id > 0 ) THEN
             CALL histwrite_p (hist2_id, 'Evap', itau_sechiba, zvevapp, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'SWnet',    itau_sechiba, swnet, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'Qh', itau_sechiba, zfluxsens, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'Qle',  itau_sechiba, zfluxlat, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'AvgSurfT', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
             CALL histwrite_p (hist2_id, 'RadT', itau_sechiba, ztemp_sol_new, kjpindex, kindex)
          ENDIF
       ENDIF

       ! Syncronize output but only if flag ok_histsync is set to true
       IF (ok_histsync .AND. is_omp_root) THEN
          IF ( (dw .EQ. xrdt) .AND. hist_id > 0 ) THEN
             CALL histsync(hist_id)
          ENDIF
       ENDIF
       
    ENDIF
    !
    ! 7. Transform the water fluxes into Kg/m^2s and m^3/s
    !
    DO ik=1, kjpindex

       vevapp(ik) = vevapp(ik)/xrdt
       coastalflow(ik) = coastalflow(ik)/xrdt
       riverflow(ik) = riverflow(ik)/xrdt

    ENDDO
    
    ! Fields for emission variables : to be transport by GCM to chemistry model.
    DO i_fields_out=1,nb_fields_out
       SELECT CASE(TRIM(field_out_names(i_fields_out)))
       CASE("fCO2_land") 
          fields_out(:,i_fields_out)=znetco2(:)
       CASE("fCO2_land_use")
          fields_out(:,i_fields_out)=zcarblu(:)
       CASE DEFAULT 
          CALL ipslerr_p (3,'intsurf_gathered', &
            &          'You ask from GCM an unknown field '//TRIM(field_out_names(i_fields_out))//&
            &          ' to ORCHIDEE for this specific version.',&
            &          'This model won''t be able to continue.', &
            &          '(check your tracer parameters in GCM)')
       END SELECT
    ENDDO
    !
    IF ( lrestart_write .AND. ok_watchout .AND. is_root_prc ) THEN
       CALL watchout_close()
    ENDIF
    !
    IF(l_first_intersurf .AND. is_root_prc) CALL getin_dump
    l_first_intersurf = .FALSE.
    !
    IF (long_print) WRITE (numout,*) ' intersurf_main done '

!Test ORCHIDEE Outputs
    IF (check_INPUTS) THEN
       IF ( sec <= 7200 ) THEN
          CALL WriteFieldI_p("Wcdrag",cdrag)
          WRITE(numout,*) "Cdrag :",cdrag
          CALL WriteFieldI_p("Wz0",z0)
          WRITE(numout,*) "Surface roughness :",z0
          CALL WriteFieldI_p("Wcoastalflow",coastalflow)
          WRITE(numout,*) "Diffuse flow of water into the ocean (m^3/dt) :",coastalflow
          CALL WriteFieldI_p("Wriverflow",riverflow)
          WRITE(numout,*) "Largest rivers flowing into the ocean (m^3/dt) :",riverflow
          CALL WriteFieldI_p("Wtsol_rad",tsol_rad)
          WRITE(numout,*) "Radiative surface temperature :",tsol_rad
          CALL WriteFieldI_p("Wvevapp",vevapp)
          WRITE(numout,*) "Total of evaporation :",vevapp
          CALL WriteFieldI_p("Wtemp_sol_new",temp_sol_new)
          WRITE(numout,*) "New soil temperature :",temp_sol_new
          CALL WriteFieldI_p("Wqsurf",qsurf)
          WRITE(numout,*) "Surface specific humidity :",qsurf
          CALL WriteFieldI_p("Walbedo_nir",albedo(:,1))
          WRITE(numout,*) "Albedo nir:",albedo(:,1)
          CALL WriteFieldI_p("Walbedo_vis",albedo(:,2))
          WRITE(numout,*) "Albedo vir :",albedo(:,2)
          CALL WriteFieldI_p("Wfluxsens",fluxsens)
          WRITE(numout,*) "Sensible chaleur flux :",fluxsens
          CALL WriteFieldI_p("Wfluxlat",fluxlat)
          WRITE(numout,*) "Latent chaleur flux :",fluxlat
          CALL WriteFieldI_p("Wemis",emis)
          WRITE(numout,*) "Emissivity :",emis
       ENDIF
    ENDIF
    !
    CALL ipslnlf_p(new_number=old_fileout)
    !        

    !
    !! Finalize the XIOS orchidee context if it is the last call
    !
    IF (lrestart_write) THEN
       CALL xios_orchidee_context_finalize
    END IF
    !
    !! Change back to be in the LMDZ context for XIOS
    !
    CALL xios_orchidee_change_context("LMDZ")

  END SUBROUTINE intersurf_gathered
!


  SUBROUTINE intsurf_time(istp, date0, dt)
    !
    !  This subroutine initialized the time global variables in grid module.
    !
    IMPLICIT NONE
    !
    INTEGER(i_std), INTENT(in)                  :: istp      !! Time step of the restart file
    REAL(r_std), INTENT(in)                     :: date0     !! The date at which itau = 0
    REAL(r_std), INTENT(in)                     :: dt        !! Time step
    !

    IF (l_first_intersurf) THEN

       ! Initialization de dt_sechiba in module constants_var
       dt_sechiba=dt

       CALL ioget_calendar(calendar_str)
       CALL ioget_calendar(one_year, one_day)
       CALL tlen2itau('1Y',dt,date0,year_length)
       IF ( TRIM(calendar_str) .EQ. 'gregorian' ) THEN  
          year_spread=un
       ELSE
          year_spread = one_year/365.2425
       ENDIF

       IF (check_time) THEN
          write(numout,*) "calendar_str =",calendar_str
          write(numout,*) "one_year=",one_year,", one_day=",one_day
          write(numout,*) "dt=",dt,", date0=",date0,", year_length=",year_length,", year_spread=",year_spread
       ENDIF
    ENDIF

    !
    IF (check_time) &
         WRITE(numout,*) "---" 
    ! Dans diffuco (ie date0 == date0_shift !!) 

    IF ( TRIM(calendar_str) .EQ. 'gregorian' ) THEN  
       !
       ! Get Julian date
       in_julian = itau2date(istp, date0, dt)

       ! Real date
       CALL ju2ymds (in_julian, year, month, day, sec)
!!$       jur=zero
!!$       julian_diff = in_julian
!!$       month_len = ioget_mon_len (year,month)
!!$       IF (check_time) THEN
!!$          write(numout,*) "in_julian, jur, julian_diff=",in_julian, jur, julian_diff
!!$          write(numout,*) 'DATE ymds', year, month,'(',month_len,'d)', day, sec, '-- stp --', istp
!!$       ENDIF

       ! julian number for january, the first.
       CALL ymds2ju (year,1,1,zero, julian0)
       julian_diff = in_julian-julian0
       ! real number of seconds
!       sec = (julian_diff-REAL(INT(julian_diff)))*one_day
       sec = NINT((julian_diff-REAL(INT(julian_diff)))*one_day)
       month_len = ioget_mon_len (year,month)
       IF (check_time) THEN
          write(numout,*) "2 in_julian, julian0, julian_diff=",in_julian, julian0, julian_diff
          write(numout,*) '2 DATE ymds', year, month,'(',month_len,'d)', day, sec, '-- stp --', istp
       ENDIF
    ELSE 
!!$       in_julian = itau2date(istp-1, zero, dt)
!!$       CALL ju2ymds (in_julian, year, month, day, sec)
!!$       jur=zero
!!$       julian_diff = in_julian
!!$       month_len = ioget_mon_len (year,month)
!!$       IF (check_time) THEN
!!$          write(numout,*) "in_julian=",in_julian, jur, julian_diff
!!$          write(numout,*) 'DATE ymds', year, month,'(',month_len,'d)', day, sec, '-- stp --', istp
!!$       ENDIF
!!$
!!$
!!$       CALL ymds2ju (year,1,1,zero, jur)
!!$       julian_diff = in_julian-jur
!!$       CALL ju2ymds (julian_diff, year, month, day, sec)
!!$!       sec = (julian_diff-REAL(INT(julian_diff)))*one_day
!!$       sec = NINT((julian_diff-REAL(INT(julian_diff)))*one_day)
!!$       month_len = ioget_mon_len (year,month)
!!$       IF (check_time) THEN
!!$          write(numout,*) "2 in_julian, jur, julian_diff=",in_julian, jur, julian_diff
!!$          write(numout,*) '2 DATE ymds', year, month,'(',month_len,'d)', day, sec, '-- stp --', istp
!!$       ENDIF


!!$       IF (check_time) &
!!$            WRITE(numout,*) "-"

!MM
!PB date0 = celui de Soenke (à tester avec un autre date0)
!       in_julian = itau2date(istp, 153116., dt)
       in_julian = itau2date(istp, date0, dt)
       CALL itau2ymds(istp, dt, year, month, day, sec)
       CALL ymds2ju (year,1,1,zero, julian0)
       julian_diff = in_julian - julian0
       month_len = ioget_mon_len (year,month)
       IF (check_time) THEN
          write(numout,*) "in_julian=",in_julian, julian0, julian_diff
          write(numout,*) 'DATE ymds', year, month,'(',month_len,'d)', day, sec, '-- stp --', istp
       ENDIF
    ENDIF
!!$    IF (check_time) &
!!$         WRITE(numout,*) "---" 

  END SUBROUTINE intsurf_time
!

!-------------------------------------------------------------------------------------
!
  SUBROUTINE intsurf_config(control_flags,dt)

    !
    !  This subroutine reads all the configuration flags which control the behaviour of the model
    !
    IMPLICIT NONE
    !
    REAL, INTENT(in)                           :: dt            !! Time step in seconds
    !
    TYPE(control_type), INTENT(out)            :: control_flags !! Flags that (de)activate parts of the model
    ! LOCAL
    INTEGER(i_std)                         :: jv

    !
    !Config Key   = LONGPRINT
    !Config Desc  = ORCHIDEE will print more messages
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = This flag permits to print more debug messages in the run.
    !Config Units = [FLAG]
    !
    long_print = .FALSE.
    CALL getin_p('LONGPRINT',long_print)
    !
    !Config Key   = CHECKTIME
    !Config Desc  = ORCHIDEE will print messages on time
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = This flag permits to print debug messages on the time.
    !Config Units = [FLAG]
    !
    check_time = .FALSE.
    CALL getin_p('CHECKTIME',check_time)
    !
    !Config Key   = SOILTYPE_CLASSIF
    !Config Desc  = Type of classification used for the map of soil types 
    !Config Def   = zobler
    !Config If    = !IMPOSE_VEG
    !Config Help  = The classification used in the file that we use here 
    !Config         There are three classification supported:  
    !Config         FAO (3 soil types), Zobler (7 converted to 3) and USDA (12) 
    !Config Units = [-]
    !
    !-tdo- Suivant le type de classification utilisee pour le sol, on adapte nscm 
    soil_classif = 'zobler'
    CALL getin_p('SOILTYPE_CLASSIF',soil_classif)
    SELECTCASE (soil_classif)
    CASE ('zobler', 'fao','none')
       nscm = nscm_fao
    CASE ('usda')
       nscm = nscm_usda
    CASE DEFAULT
       WRITE(numout,*) "Unsupported soil type classification. Choose between zobler, fao and usda according to the map"
       STOP 'intsurf_config'
    ENDSELECT
    !
    !Config Key   = ORCHIDEE_WATCHOUT
    !Config Desc  = ORCHIDEE will write out its forcing to a file
    !Config If    =
    !Config Def   = n
    !Config Help  = This flag allows to write to a file all the variables
    !Config         which are used to force the land-surface. The file 
    !Config         has exactly the same format than a normal off-line forcing
    !Config         and thus this forcing can be used for forcing ORCHIDEE.
    !Config Units = [FLAG]
    !
    ok_watchout = .FALSE.
    CALL getin_p('ORCHIDEE_WATCHOUT',ok_watchout)
    !
    IF (ok_watchout) THEN
       !Config Key   = DT_WATCHOUT
       !Config Desc  = ORCHIDEE will write out with this frequency
       !Config If    = ORCHIDEE_WATCHOUT
       !Config Def   = dt
       !Config Help  = This flag indicates the frequency of the write of the variables. 
       !Config Units = [seconds]
       !
       dt_watch = dt
       CALL getin_p('DT_WATCHOUT',dt_watch)
       dt_split_watch = dt_watch / dt
       !
       !Config Key   = WATCHOUT_FILE
       !Config Desc  = Filenane for the ORCHIDEE forcing file
       !Config If    = ORCHIDEE_WATCHOUT
       !Config Def   = orchidee_watchout.nc
       !Config Help  = This is the name of the file in which the
       !Config         forcing used here will be written for later use. 
       !Config Units = [FILE]
       !
       watchout_file = "orchidee_watchout.nc"
       CALL getin_p('WATCHOUT_FILE',watchout_file)
       
       WRITE(numout,*) 'WATCHOUT flag :', ok_watchout
       WRITE(numout,*) 'WATCHOUT file :', watchout_file
    ENDIF
    !
    !Config Key   = RIVER_ROUTING
    !Config Desc  = Decides if we route the water or not
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = This flag allows the user to decide if the runoff
    !Config         and drainage should be routed to the ocean
    !Config         and to downstream grid boxes.
    !Config Units = [FLAG]
    !
    control_flags%river_routing = .FALSE.
    CALL getin_p('RIVER_ROUTING', control_flags%river_routing)
    WRITE(numout,*) "RIVER routing is activated : ",control_flags%river_routing
    !
    !Config key   = HYDROL_CWRR
    !Config Desc  = Allows to switch on the multilayer hydrology of CWRR
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = This flag allows the user to decide if the vertical
    !Config         hydrology should be treated using the multi-layer 
    !Config         diffusion scheme adapted from CWRR by Patricia de Rosnay.
    !Config         by default the Choisnel hydrology is used.
    !Config Units = [FLAG]
    !
    control_flags%hydrol_cwrr = .FALSE.
    CALL getin_p('HYDROL_CWRR', control_flags%hydrol_cwrr)
    WRITE(numout,*) "CWRR hydrology is activated : ",control_flags%hydrol_cwrr
    !
    IF ( control_flags%hydrol_cwrr ) THEN
       dpu_max = 2.0
    ELSE
       dpu_max = 4.0
    END IF

    !Config Key   = HYDROL_SOIL_DEPTH
    !Config Desc  = Total depth of soil reservoir
    !Config If    = OK_SECHIBA 
    !Config Def   = 4./ 2. (if HYDROL_CWRR)
    !Config Help  = By default, ORCHIDEE uses the AR5 configuration (Choisnel-4m).
    !Config Units = [m]
    CALL getin_p ("HYDROL_SOIL_DEPTH", dpu_max)
    WRITE(numout,*) "Soil depth value read in configuration file run.def : ", dpu_max
    !
    DO jv = 1, nbdl-1
       diaglev(jv) = dpu_max/(2**(nbdl-1) -1) * ( ( 2**(jv-1) -1) + ( 2**(jv) -1) ) / deux
    ENDDO
    diaglev(nbdl) = dpu_max
    !
    !Config Key   = DO_IRRIGATION
    !Config Desc  = Should we compute an irrigation flux 
    !Config If    = RIVER_ROUTING 
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         to compute an irigation flux. This performed for the
    !Config         on very simple hypothesis. The idea is to have a good
    !Config         map of irrigated areas and a simple function which estimates
    !Config         the need to irrigate.
    !Config Units = [FLAG]
    !
    control_flags%do_irrigation = .FALSE.
    CALL getin_p('DO_IRRIGATION', control_flags%do_irrigation)
    !
    !Config Key   = DO_FLOODPLAINS
    !Config Desc  = Should we include floodplains 
    !Config If    = RIVER_ROUTING 
    !Config Def   = n
    !Config Help  = This parameters allows the user to ask the model
    !Config         to take into account the flood plains and return 
    !Config         the water into the soil moisture. It then can go 
    !Config         back to the atmopshere. This tried to simulate 
    !Config         internal deltas of rivers.
    !Config Units = [FLAG]  
    !
    control_flags%do_floodplains = .FALSE.
    CALL getin_p('DO_FLOODPLAINS', control_flags%do_floodplains)
    !
    !Config Key   = CHECK_WATERBAL
    !Config Desc  = Should we check the global water balance 
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = This parameters allows the user to check
    !Config         the integrated water balance at the end
    !Config         of each time step
    !Config Units = [FLAG]  
    !
    check_waterbal = .FALSE.
    CALL getin_p('CHECK_WATERBAL', check_waterbal)

    !Config Key   = OK_EXPLICITSNOW
    !Config Desc  = Activate explict snow scheme
    !Config If    = OK_SECHIBA
    !Config Def   = FALSE
    !Config Help  = Activate explicit snow scheme instead of default snow scheme
    !Config Units = [FLAG]
    ok_explicitsnow = .FALSE.
    CALL getin_p('OK_EXPLICITSNOW', ok_explicitsnow)

    !Config Key   = OK_EXPLICITSNOW
    !Config Desc  = Activate explict snow scheme
    !Config If    = OK_SECHIBA
    !Config Def   = FALSE
    !Config Help  = Activate explicit snow scheme instead of default snow scheme
    !Config Units = [FLAG]
    ok_pc = .FALSE.
    CALL getin_p('OK_PC', ok_pc)

    !
    !Config Key   = STOMATE_OK_CO2
    !Config Desc  = Activate CO2?
    !Config If    = OK_SECHIBA 
    !Config Def   = n
    !Config Help  = set to TRUE if photosynthesis is to be activated
    !Config Units = [FLAG]
    !
    control_flags%ok_co2 = .FALSE.
    CALL getin_p('STOMATE_OK_CO2', control_flags%ok_co2)
    WRITE(numout,*) 'photosynthesis: ', control_flags%ok_co2
    !
    !Config Key   = STOMATE_OK_STOMATE
    !Config Desc  = Activate STOMATE?
    !Config If    = OK_SECHIBA and OK_CO2
    !Config Def   = n
    !Config Help  = set to TRUE if STOMATE is to be activated
    !Config Units = [FLAG]
    !
    control_flags%ok_stomate = .FALSE.
    CALL getin_p('STOMATE_OK_STOMATE',control_flags%ok_stomate)
    WRITE(numout,*) 'STOMATE is activated: ',control_flags%ok_stomate
    !
    !Config Key   = STOMATE_OK_DGVM
    !Config Desc  = Activate DGVM?
    !Config If    = OK_STOMATE
    !Config Def   = n
    !Config Help  = set to TRUE if DGVM is to be activated
    !Config Units = [FLAG]
    !
    control_flags%ok_dgvm = .FALSE.
    CALL getin_p('STOMATE_OK_DGVM',control_flags%ok_dgvm)
    !
    !Config Key   = DIFFUCO_OK_INCA
    !Config Desc  = Activate DIFFUCO_INCA?
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = set to TRUE if biogenic emissions calculation is to be activated
    !Config Units = [FLAG]
    !
    control_flags%ok_inca = .FALSE.
    CALL getin_p('DIFFUCO_OK_INCA', control_flags%ok_inca)
    WRITE(numout,*) 'Biogenic emissions: ', control_flags%ok_inca
    !
    !Config Key   = LEAFAGE_OK_INCA
    !Config Desc  = Activate LEAFAGE?
    !Config If    = DIFFUCO_OK_INCA
    !Config Def   = n
    !Config Help  = set to TRUE if biogenic emissions calculation takes leaf age into account
    !Config Units = [FLAG]
    !
    control_flags%ok_leafage = .FALSE.
    CALL getin_p('LEAFAGE_OK_INCA', control_flags%ok_leafage)
    WRITE(numout,*) 'Leaf Age: ', control_flags%ok_leafage
    !
    !Config Key   = CANOPY_EXTINCTION 
    !Config Desc  = Use canopy radiative transfer model?
    !Config If    = DIFFUCO_OK_INCA 
    !Config Def   = n
    !Config Help  = set to TRUE if canopy radiative transfer model is used for biogenic emissions 
    !Config Units = [FLAG]
    !
    control_flags%ok_radcanopy = .FALSE.
    CALL getin_p('CANOPY_EXTINCTION', control_flags%ok_radcanopy)
    WRITE(numout,*) 'Canopy radiative transfer model: ', control_flags%ok_radcanopy
    !
    !Config Key   = CANOPY_MULTILAYER
    !Config Desc  = Use canopy radiative transfer model with multi-layers
    !Config If    = CANOPY_EXTINCTION 
    !Config Def   = n
    !Config Help  = set to TRUE if canopy radiative transfer model is with multiple layers 
    !Config Units = [FLAG]
    !
    control_flags%ok_multilayer = .FALSE.
    CALL getin_p('CANOPY_MULTILAYER', control_flags%ok_multilayer)
    WRITE(numout,*) 'Multi-layer Canopy model: ', control_flags%ok_multilayer
    !
    !Config Key   = NOx_RAIN_PULSE
    !Config Desc  = Calculate NOx emissions with pulse?
    !Config If    = DIFFUCO_OK_INCA 
    !Config Def   = n
    !Config Help  = set to TRUE if NOx rain pulse is taken into account
    !Config Units = [FLAG]
    !
    control_flags%ok_pulse_NOx = .FALSE.
    CALL getin_p('NOx_RAIN_PULSE', control_flags%ok_pulse_NOx)
    WRITE(numout,*) 'Rain NOx pulsing: ', control_flags%ok_pulse_NOx
    !
    !Config Key   = NOx_BBG_FERTIL
    !Config Desc  = Calculate NOx emissions with bbg fertilizing effect?
    !Config If    = DIFFUCO_OK_INCA 
    !Config Def   = n
    !Config Help  = set to TRUE if NOx emissions are calculated with bbg effect 
    !Config         Fertil effect of bbg on NOx soil emissions 
    !Config Units = [FLAG]
    !
    control_flags%ok_bbgfertil_NOx = .FALSE.
    CALL getin_p('NOx_BBG_FERTIL', control_flags%ok_bbgfertil_NOx)
    WRITE(numout,*) 'NOx bbg fertil effect: ', control_flags%ok_bbgfertil_NOx
    !
    !Config Key   = NOx_FERTILIZERS_USE
    !Config Desc  = Calculate NOx emissions with fertilizers use?
    !Config If    = DIFFUCO_OK_INCA 
    !Config Def   = n
    !Config Help  = set to TRUE if NOx emissions are calculated with fertilizers use
    !Config         Fertilizers use effect on NOx soil emissions  
    !Config Units = [FLAG] 
    !
    control_flags%ok_cropsfertil_NOx = .FALSE.
    CALL getin_p('NOx_FERTILIZERS_USE', control_flags%ok_cropsfertil_NOx)
    WRITE(numout,*) 'NOx Fertilizers use: ', control_flags%ok_cropsfertil_NOx
 

    !
    ! control initialisation with sechiba
    !
    control_flags%ok_sechiba = .TRUE.
    !
    !
    ! Ensure consistency
    !
    IF ( control_flags%ok_dgvm ) control_flags%ok_stomate = .TRUE.
    IF ( control_flags%ok_stomate ) control_flags%ok_co2 = .TRUE.
    IF ( control_flags%ok_multilayer .AND. .NOT.(control_flags%ok_radcanopy) ) THEN
       control_flags%ok_radcanopy  = .TRUE.
       WRITE(numout,*) 'You want to use the multilayer model without activating the flag CANOPY_EXTINCTION'
       WRITE(numout,*) 'We set CANOPY_EXTINCTION to TRUE to ensure consistency'
    ENDIF
    !
    !Config Key   = STOMATE_WATCHOUT
    !Config Desc  = STOMATE does minimum service
    !Config If    = OK_SECHIBA 
    !Config Def   = n
    !Config Help  = set to TRUE if you want STOMATE to read
    !Config         and write its start files and keep track
    !Config         of longer-term biometeorological variables.
    !Config         This is useful if OK_STOMATE is not set,
    !Config         but if you intend to activate STOMATE later.
    !Config         In that case, this run can serve as a 
    !Config         spinup for longer-term biometeorological
    !Config         variables.
    !Config Units = [FLAG]
    !
    control_flags%stomate_watchout = .FALSE.
    CALL getin_p('STOMATE_WATCHOUT',control_flags%stomate_watchout)
    WRITE(numout,*) 'STOMATE keeps an eye open: ',control_flags%stomate_watchout
    !
    ! Here we need the same initialisation as above
    !
    control_flags%ok_pheno = .TRUE.

    !
    ! Configuration : number of PFTs and parameters
    !

    ! 1. Number of PFTs defined by the user

    !Config Key   = NVM
    !Config Desc  = number of PFTs  
    !Config If    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 13
    !Config Help  = The number of vegetation types define by the user
    !Config Units = [-]
    !
    CALL getin_p('NVM',nvm)
    WRITE(numout,*)'the number of pfts used by the model is : ', nvm

    ! 2. Should we read the parameters in the run.def file ?

    !Config Key   = IMPOSE_PARAM
    !Config Desc  = Do you impose the values of the parameters?
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = y
    !Config Help  = This flag can deactivate the reading of some parameters.
    !               Useful if you want to use the standard values without commenting the run.def
    !Config Units = [FLAG]
    !
    CALL getin_p('IMPOSE_PARAM',impose_param)

    ! 3. Allocate and intialize the pft parameters

    CALL pft_parameters_main(control_flags)

    ! 4. Activation sub-models of ORCHIDEE

    CALL activate_sub_models(control_flags)

    ! 5. Vegetation configuration (impose_veg, land_use, lcchange...previously in slowproc)

    CALL veget_config

    ! 6. Read the parameters in the run.def file  according the flags

    IF (impose_param ) THEN 
       CALL config_pft_parameters
    ENDIF

    IF ( control_flags%ok_sechiba ) THEN
       IF (impose_param ) THEN
          CALL config_sechiba_parameters
          CALL config_sechiba_pft_parameters(control_flags)
          WRITE(numout,*)'    some sechiba parameters have been imposed '
       ENDIF
    ENDIF
    
    ! Initialize variables in constantes_soil
    CALL config_soil_parameters(control_flags, impose_param)
        

    IF ( control_flags%ok_co2 ) THEN
       IF ( impose_param ) THEN
          CALL config_co2_parameters
          WRITE(numout,*)'    some co2 parameters have been imposed '         
       ENDIF
    ENDIF
    
    IF ( control_flags%ok_stomate ) THEN
       IF ( impose_param ) THEN
          CALL config_stomate_parameters
          CALL config_stomate_pft_parameters
          WRITE(numout,*)'    some stomate parameters have been imposed '
       ENDIF
    ENDIF
    
    IF ( control_flags%ok_dgvm ) THEN
       IF ( impose_param ) THEN
          CALL config_dgvm_parameters
          WRITE(numout,*)'    some dgvm parameters have been imposed '         
       ENDIF
    ENDIF    

    !
    !
    RETURN
    !
  END SUBROUTINE intsurf_config
  !
  !
  !
  SUBROUTINE intsurf_restart(istp, date0, dt, control_flags, rest_id, rest_id_stom, itau_offset)

    USE mod_orchidee_para
    !
    !  This subroutine initialized the restart file for the land-surface scheme
    !
    IMPLICIT NONE
    !
    INTEGER(i_std), INTENT(in)                  :: istp      !! Time step of the restart file
    REAL(r_std)                                 :: date0     !! The date at which itau = 0
    REAL(r_std)                                 :: dt        !! Time step
    INTEGER(i_std), INTENT(out)                 :: rest_id, rest_id_stom   !! ID of the restart file
    INTEGER(i_std), INTENT(out)                 :: itau_offset
    !
    TYPE(control_type), INTENT(in)             :: control_flags !! Flags that (de)activate parts of the model
    !
    !  LOCAL
    !
    REAL(r_std)                 :: dt_rest, date0_rest
    INTEGER(i_std)              :: itau_dep
    INTEGER(i_std),PARAMETER    :: llm=1
    REAL(r_std), DIMENSION(llm) :: lev
    LOGICAL                    :: overwrite_time
    REAL(r_std)                 :: in_julian, rest_julian
    INTEGER(i_std)              :: yy, mm, dd
    REAL(r_std)                 :: ss
    !
    !Config Key   = SECHIBA_restart_in
    !Config Desc  = Name of restart to READ for initial conditions
    !Config If    = OK_SECHIBA 
    !Config Def   = NONE
    !Config Help  = This is the name of the file which will be opened
    !Config         to extract the initial values of all prognostic
    !Config         values of the model. This has to be a netCDF file.
    !Config         Not truly COADS compliant. NONE will mean that
    !Config         no restart file is to be expected.
    !Config Units = [FILE]
!-
    CALL getin_p('SECHIBA_restart_in',restname_in)
    WRITE(numout,*) 'INPUT RESTART_FILE', restname_in
    !-
    !Config Key   = SECHIBA_rest_out
    !Config Desc  = Name of restart files to be created by SECHIBA
    !Config If    = OK_SECHIBA
    !Config Def   = sechiba_rest_out.nc
    !Config Help  = This variable give the name for
    !Config         the restart files. The restart software within
    !Config         IOIPSL will add .nc if needed.
    !Config Units = [FILE]
    !
    CALL getin_p('SECHIBA_rest_out', restname_out)
    !
    !Config Key   = SECHIBA_reset_time
    !Config Desc  = Override the time from the restart files for SECHIBA and STOMATE
    !Config If    = OK_SECHIBA
    !Config Def   = y
    !Config Help  = This option allows the model to override the time
    !Config         found in the restart files of SECHIBA and STOMATE with the time
    !Config         given by the GCM or the forcing file.
    !Config Units = [FLAG]
    !
    overwrite_time = .TRUE.
    CALL getin_p('SECHIBA_reset_time', overwrite_time)
    !
    lev(:) = zero
    itau_dep = istp
    in_julian = itau2date(istp, date0, dt)
    date0_rest = date0
    dt_rest = dt
    !
    IF (is_root_prc) THEN
      CALL restini( restname_in, iim_g, jjm_g, lon_g, lat_g, llm, lev, &
         &  restname_out, itau_dep, date0_rest, dt_rest, rest_id, overwrite_time)
    ELSE
       rest_id=0
    ENDIF
    CALL bcast (itau_dep)
    CALL bcast (date0_rest)
    CALL bcast (dt_rest)
    !
    !  itau_dep of SECHIBA is phased with the GCM if needed
    !
    rest_julian = itau2date(itau_dep, date0_rest, dt_rest)
    !
    IF ( ABS( in_julian - rest_julian) .GT. dt/one_day .AND. .NOT. OFF_LINE_MODE ) THEN
       IF ( overwrite_time ) THEN
          WRITE(numout,*) 'The SECHIBA restart is not for the same timestep as the GCM,'
          WRITE(numout,*) 'the two are synchronized. The land-surface conditions can not impose'
          WRITE(numout,*) 'the chronology of the simulation.'
          WRITE(numout,*) 'Time step of the GCM :', istp, 'Julian day : ', in_julian
          CALL ju2ymds(in_julian, yy, mm, dd, ss)
          WRITE(numout,*) 'In other word (yy,mm,dd,ss) : ', yy, mm, dd, ss
          WRITE(numout,*) 'Time step of SECHIBA :', itau_dep, 'Julian day : ', rest_julian
          CALL ju2ymds(rest_julian, yy, mm, dd, ss)
          WRITE(numout,*) 'In other word (yy,mm,dd,ss) : ', yy, mm, dd, ss

          itau_offset = itau_dep - istp
          date0_shifted = date0 - itau_offset*dt/one_day
!MM_ A VOIR : dans le TAG 1.4 :
!         date0_shifted = in_julian - itau_dep*dt/one_day
!MM_ Bon calcul ?

          WRITE(numout,*) 'The new starting date is :', date0_shifted
          CALL ju2ymds(date0_shifted, yy, mm, dd, ss)
          WRITE(numout,*) 'In other word (yy,mm,dd,ss) : ', yy, mm, dd, ss
       ELSE
          WRITE(numout,*) 'IN -> OUT :', istp, '->', itau_dep
          WRITE(numout,*) 'IN -> OUT :', in_julian, '->', rest_julian
          WRITE(numout,*) 'SECHIBA''s restart file is not consistent with the one of the GCM'
          WRITE(numout,*) 'Correct the time axis of the restart file or force the code to change it.'
          STOP
       ENDIF
    ELSE
       itau_offset = 0
       date0_shifted = date0
    ENDIF
    !
!!!    CALL ioconf_startdate(date0_shifted)
    !
    !=====================================================================
    !- 1.5 Restart file for STOMATE
    !=====================================================================
    IF ( control_flags%ok_stomate .OR. control_flags%stomate_watchout ) THEN 
       !-
       ! STOMATE IS ACTIVATED
       !-
       !Config Key   = STOMATE_RESTART_FILEIN
       !Config Desc  = Name of restart to READ for initial conditions of STOMATE
       !Config If    = STOMATE_OK_STOMATE or STOMATE_WATCHOUT
       !Config Def   = NONE
       !Config Help  = This is the name of the file which will be opened
       !Config         to extract the initial values of all prognostic
       !Config         values of STOMATE.
       !Config Units = [FILE]
       !-
       CALL getin_p('STOMATE_RESTART_FILEIN',stom_restname_in)
       WRITE(numout,*) 'STOMATE INPUT RESTART_FILE', stom_restname_in
       !-
       !Config Key   = STOMATE_RESTART_FILEOUT
       !Config Desc  = Name of restart files to be created by STOMATE
       !Config If    = STOMATE_OK_STOMATE or STOMATE_WATCHOUT
       !Config Def   = stomate_restart.nc
       !Config Help  = This is the name of the file which will be opened
       !Config         to write the final values of all prognostic values
       !Config         of STOMATE.
       !Config Units = [FILE]
       !-
       CALL getin_p('STOMATE_RESTART_FILEOUT', stom_restname_out)
       WRITE(numout,*) 'STOMATE OUTPUT RESTART_FILE', stom_restname_out
       !-
       IF (is_root_prc) THEN
         CALL restini (stom_restname_in, iim_g, jjm_g, lon_g, lat_g, llm, lev, &
            &  stom_restname_out, itau_dep, date0_rest, dt_rest, rest_id_stom, overwrite_time)
       ELSE
         rest_id_stom=0
       ENDIF
       CALL bcast (itau_dep)
       CALL bcast (date0_rest)
       CALL bcast (dt_rest)
       !-
    ENDIF
    !
  END SUBROUTINE intsurf_restart
  
  SUBROUTINE intsurf_history(iim, jjm, lon, lat, kindex, kjpindex, istp_old, date0, dt, control_flags, hist_id, hist2_id, &
       & hist_id_stom, hist_id_stom_IPCC)

    USE mod_orchidee_para

    !
    !   
    !  This subroutine initialized the history files for the land-surface scheme
    !
    IMPLICIT NONE
    !
    INTEGER(i_std), INTENT(in)                  :: iim, jjm  !! Size in x and y of the data to be handeled
    REAL(r_std),DIMENSION (iim,jjm), INTENT(in) :: lon, lat  !! Longitude and latitude of the data points
    INTEGER(i_std),INTENT (in)                            :: kjpindex
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: kindex
    
    INTEGER(i_std), INTENT(in)                  :: istp_old  !! Time step counter
    REAL(r_std), INTENT(in)                     :: date0     !! Julian day at which istp=0
    REAL(r_std), INTENT(in)                     :: dt        !! Time step of the counter in seconds
    !
    TYPE(control_type), INTENT(in)             :: control_flags !! Flags that (de)activate parts of the model
    !
    INTEGER(i_std), INTENT(out)                 :: hist_id !! History file identification for SECHIBA
    INTEGER(i_std), INTENT(out)                 :: hist2_id !! History file 2 identification for SECHIBA (Hi-frequency ?)
    !! History file identification for STOMATE and IPCC
    INTEGER(i_std), INTENT(out)                 :: hist_id_stom, hist_id_stom_IPCC 
    !
    !  LOCAL
    !
    CHARACTER(LEN=80) :: histname,histname2                    !! Name of history files for SECHIBA
    CHARACTER(LEN=80) :: stom_histname, stom_ipcc_histname     !! Name of history files for STOMATE
    LOGICAL           :: ok_histfile2                 !! Flag to switch on histfile 2 for SECHIBA
    REAL(r_std)       :: dw2                          !! frequency of history write (sec.)
    CHARACTER(LEN=30)   :: flux_op                    !! Operations to be performed on fluxes
!    CHARACTER(LEN=30)   :: flux_sc                    !! Operations which do not include a scatter
    CHARACTER(LEN=40)   :: flux_insec, flux_scinsec   !! Operation in seconds
    INTEGER(i_std)     :: hist_level, hist2_level     !! history output level (default is 10 => maximum output)
    CHARACTER(LEN=40),DIMENSION(max_hist_level) :: &
         & ave, avecels, avescatter, fluxop, &
         & fluxop_scinsec, tmincels, tmaxcels, once, sumscatter  !! The various operation to be performed
!!, tmax, fluxop_sc, fluxop_insec, &
    CHARACTER(LEN=40),DIMENSION(max_hist_level) :: &
         & ave2, avecels2, avescatter2, fluxop2, &
         & fluxop_scinsec2, tmincels2, tmaxcels2, once2, sumscatter2  !! The various operation to be performed
!!, tmax2, fluxop_sc2, fluxop_insec2, &
    CHARACTER(LEN=80) :: global_attribute              !! for writing attributes in the output files
    INTEGER(i_std)     :: i, jst
    ! SECHIBA AXIS
    INTEGER(i_std)     :: hori_id                      !! ID of the default horizontal longitude and latitude map.
    INTEGER(i_std)     :: vegax_id, laiax_id, solax_id, soltax_id, nobioax_id !! ID's for two vertical coordinates
    INTEGER(i_std)     :: soildiagax_id                !! ID for diagnostic soil levels
    INTEGER(i_std)     :: solayax_id                   !! ID for the vertical axis of the CWRR hydrology 
    INTEGER(i_std)     :: hori_id2                      !! ID of the default horizontal longitude and latitude map.
    INTEGER(i_std)     :: vegax_id2, laiax_id2, solax_id2, soltax_id2, nobioax_id2, albax_id2 !! ID's for two vertical coordinates
    INTEGER(i_std)     :: solayax_id2                   !! ID for the vertical axis of the CWRR hydrology 
    INTEGER(i_std)     :: snowax_id                     !! ID for snow level axis

    ! STOMATE AXIS
    INTEGER(i_std)     :: hist_PFTaxis_id
! deforestation
    INTEGER(i_std)     :: hist_pool_10axis_id
    INTEGER(i_std)     :: hist_pool_100axis_id
    INTEGER(i_std)     :: hist_pool_11axis_id
    INTEGER(i_std)     :: hist_pool_101axis_id
    ! STOMATE IPCC AXIS
    INTEGER(i_std)     :: hist_IPCC_PFTaxis_id
    !
    ! permafrost vert axis for stomate
    INTEGER(i_std)     :: hist_stomate_deepsoil
    INTEGER(i_std)     :: hist_stomate_snow
    CHARACTER(LEN=10)  :: part_str                      !! string suffix indicating an index
    REAL(r_std),DIMENSION(nsnow)  :: snowlev            !! snow axis
    REAL(r_std),DIMENSION(ngrnd) :: sol_coef

    LOGICAL                               :: rectilinear
    INTEGER(i_std)                         :: ier,jv
    REAL(r_std), ALLOCATABLE, DIMENSION(:) :: lon_rect, lat_rect
    !
    REAL(r_std),DIMENSION(nvm)   :: veg
    REAL(r_std),DIMENSION(nlai+1):: lai
    REAL(r_std),DIMENSION(ngrnd) :: sol
    REAL(r_std),DIMENSION(nstm)  :: soltyp
    REAL(r_std),DIMENSION(nnobio):: nobiotyp
    REAL(r_std),DIMENSION(2)     :: albtyp
    REAL(r_std),DIMENSION(nslm)  :: solay
    !
    CHARACTER(LEN=80)           :: var_name           !! To store variables names
    !
    ! STOMATE history file
    REAL(r_std)                  :: hist_days_stom     !!- GK time step in days for this history file
    REAL(r_std)                  :: hist_dt_stom       !!- GK time step in seconds for this history file
    REAL(r_std)                  :: dt_slow_           !!  for test : time step of slow processes and STOMATE
    REAL(r_std),DIMENSION(nvm)   :: hist_PFTaxis       !!- GK An axis we need for the history files
!
    REAL(r_std),DIMENSION(10)  :: hist_pool_10axis     !! Deforestation axis
    REAL(r_std),DIMENSION(100)  :: hist_pool_100axis     !! Deforestation axis
    REAL(r_std),DIMENSION(11)  :: hist_pool_11axis     !! Deforestation axis
    REAL(r_std),DIMENSION(101)  :: hist_pool_101axis     !! Deforestation axis
    !
    ! IPCC history file
    REAL(r_std)                  :: hist_days_stom_ipcc     !!- GK time step in days for this history file
    REAL(r_std)                  :: hist_dt_stom_ipcc       !!- GK time step in seconds for this history file
!
    !
    !=====================================================================
    !- 3.0 Setting up the history files
    !=====================================================================
    !- 3.1 SECHIBA
    !=====================================================================
    !Config Key   = ALMA_OUTPUT
    !Config Desc  = Should the output follow the ALMA convention
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = If this logical flag is set to true the model
    !Config         will output all its data according to the ALMA 
    !Config         convention. It is the recommended way to write
    !Config         data out of ORCHIDEE.
    !Config Units = [FLAG]
    !-
    almaoutput = .FALSE.
    CALL getin_p('ALMA_OUTPUT', almaoutput)    
    WRITE(numout,*) 'ALMA_OUTPUT', almaoutput
    !-
    !Config Key   = OUTPUT_FILE
    !Config Desc  = Name of file in which the output is going to be written
    !Config If    = OK_SECHIBA
    !Config Def   = sechiba_history.nc
    !Config Help  = This file is going to be created by the model
    !Config         and will contain the output from the model.
    !Config         This file is a truly COADS compliant netCDF file.
    !Config         It will be generated by the hist software from
    !Config         the IOIPSL package.
    !Config Units = [FILE]
    !-
    histname='sechiba_history.nc'
    CALL getin_p('OUTPUT_FILE', histname)
    WRITE(numout,*) 'OUTPUT_FILE', histname
    !-
    !Config Key   = WRITE_STEP
    !Config Desc  = Frequency in seconds at which to WRITE output
    !Config If    = OK_SECHIBA
    !Config Def   = one_day
    !Config Help  = This variables gives the frequency the output of
    !Config         the model should be written into the netCDF file.
    !Config         It does not affect the frequency at which the
    !Config         operations such as averaging are done.
    !Config         That is IF the coding of the calls to histdef
    !Config         are correct !
    !Config Units = [seconds]
    !-
    dw = one_day
    CALL getin_p('WRITE_STEP', dw)
    IF ( dw == 0 ) WRITE(numout,*) ' sechiba_history file will not be created' 
    !
    veg(1:nvm)   = (/ (REAL(i,r_std),i=1,nvm) /)
    lai(1:nlai+1) = (/ (REAL(i,r_std),i=1,nlai+1) /)
    soltyp(1:nstm) = (/ (REAL(i,r_std),i=1,nstm) /)
    nobiotyp(1:nnobio) = (/ (REAL(i,r_std),i=1,nnobio) /)
    albtyp(1:2) = (/ (REAL(i,r_std),i=1,2) /)
    solay(1:nslm) = (/ (REAL(i,r_std),i=1,nslm) /)
    snowlev =  (/ (REAL(i,r_std),i=1,nsnow) /)

    ! Get the vertical soil levels for the thermal scheme
    sol(1:ngrnd) = thermosoil_levels()

    !
    !- We need to flux averaging operation as when the data is written
    !- from within SECHIBA a scatter is needed. In the driver on the other
    !- hand the data is 2D and can be written is it is.
    !-
    WRITE(flux_op,'("ave(scatter(X*",F8.1,"))")') one_day/dt
    ! WRITE(flux_op,'("(ave(scatter(X))*",F8.1,")")') one_day/dt
!    WRITE(flux_sc,'("ave(X*",F8.1,")")') one_day/dt
!    WRITE(flux_insec,'("ave(X*",F8.6,")")') un/dt
!    WRITE(flux_insec,'("ave(X*",F12.10,")")') un/dt
    WRITE(flux_scinsec,'("ave(scatter(X*",F13.10,"))")') un/dt
    WRITE(numout,*) 'flux_op=',flux_op,' one_day/dt=', one_day/dt, ' dt=',dt,' dw=', dw
    !-
    !Config Key   = SECHIBA_HISTLEVEL
    !Config Desc  = SECHIBA history output level (0..10)
    !Config If    = OK_SECHIBA and HF
    !Config Def   = 5
    !Config Help  = Chooses the list of variables in the history file. 
    !Config         Values between 0: nothing is written; 10: everything is 
    !Config         written are available More details can be found on the web under documentation.
    !Config Units = [-]
    !-
    hist_level = 5
    CALL getin_p('SECHIBA_HISTLEVEL', hist_level)
    !-
    WRITE(numout,*) 'SECHIBA history level: ',hist_level
    IF ( (hist_level > max_hist_level).OR.(hist_level < 0) ) THEN
       STOP 'This history level is not allowed'
    ENDIF
    !-
    !- define operations as a function of history level.
    !- Above hist_level, operation='never'
    !-
    ave(1:max_hist_level) = 'ave(scatter(X))'
    IF (hist_level < max_hist_level) THEN
       ave(hist_level+1:max_hist_level) = 'never'
    ENDIF
    sumscatter(1:max_hist_level) = 't_sum(scatter(X))'
    IF (hist_level < max_hist_level) THEN
       sumscatter(hist_level+1:max_hist_level) = 'never'
    ENDIF

!    avecels(1:max_hist_level) = 'ave(cels(X))'
!    IF (hist_level < max_hist_level) THEN
!       avecels(hist_level+1:max_hist_level) = 'never'
!    ENDIF

    avecels(1:max_hist_level) = 'ave(cels(scatter(X)))'
    IF (hist_level < max_hist_level) THEN
       avecels(hist_level+1:max_hist_level) = 'never'
    ENDIF

    avescatter(1:max_hist_level) = 'ave(scatter(X))'
    IF (hist_level < max_hist_level) THEN
       avescatter(hist_level+1:max_hist_level) = 'never'
    ENDIF
    tmincels(1:max_hist_level) = 't_min(cels(scatter(X)))'
    IF (hist_level < max_hist_level) THEN
       tmincels(hist_level+1:max_hist_level) = 'never'
    ENDIF
    tmaxcels(1:max_hist_level) = 't_max(cels(scatter(X)))'
    IF (hist_level < max_hist_level) THEN
       tmaxcels(hist_level+1:max_hist_level) = 'never'
    ENDIF
!!$    tmax(1:max_hist_level) = 't_max(X)'
!!$    IF (hist_level < max_hist_level) THEN
!!$       tmax(hist_level+1:max_hist_level) = 'never'
!!$    ENDIF
    fluxop(1:max_hist_level) = flux_op
    IF (hist_level < max_hist_level) THEN
       fluxop(hist_level+1:max_hist_level) = 'never'
    ENDIF
!!$    fluxop_sc(1:max_hist_level) = flux_sc
!!$    IF (hist_level < max_hist_level) THEN
!!$       fluxop_sc(hist_level+1:max_hist_level) = 'never'
!!$    ENDIF
!!$    fluxop_insec(1:max_hist_level) = flux_insec
!!$    IF (hist_level < max_hist_level) THEN
!!$       fluxop_insec(hist_level+1:max_hist_level) = 'never'
!!$    ENDIF
    fluxop_scinsec(1:max_hist_level) = flux_scinsec
    IF (hist_level < max_hist_level) THEN
       fluxop_scinsec(hist_level+1:max_hist_level) = 'never'
    ENDIF
    once(1:max_hist_level) = 'once(scatter(X))'
    IF (hist_level < max_hist_level) THEN
       once(hist_level+1:max_hist_level) = 'never'
    ENDIF


    !- Initialize sechiba_history output file
    !- 
    IF ( dw == 0 ) THEN
       ! sechiba_history file will not be created.
       hist_id = -1

    ELSE
       ! sechiba_history file will be created

       !- Calculation necessary for initialization of sechiba_history file
       !- Check if we have by any change a rectilinear grid. This would allow us to 
       !- simplify the output files.
    IF (is_omp_root) THEN

       rectilinear = .FALSE.
       IF ( ALL(lon(:,:) == SPREAD(lon(:,1), 2, SIZE(lon,2))) .AND. &
            & ALL(lat(:,:) == SPREAD(lat(1,:), 1, SIZE(lat,1))) ) THEN
          rectilinear = .TRUE.
          ALLOCATE(lon_rect(iim),stat=ier)
          IF (ier .NE. 0) THEN
             WRITE (numout,*) ' error in lon_rect allocation. We stop. We need iim words = ',iim
             STOP 'intersurf_history'
          ENDIF
          ALLOCATE(lat_rect(jjm),stat=ier)
          IF (ier .NE. 0) THEN
             WRITE (numout,*) ' error in lat_rect allocation. We stop. We need jjm words = ',jjm
             STOP 'intersurf_history'
          ENDIF
          lon_rect(:) = lon(:,1)
          lat_rect(:) = lat(1,:)
       ENDIF
       !-
       !-
       !-
       ! Initialize sechiba_history file
       IF ( .NOT. almaoutput ) THEN
          !- 
          IF ( rectilinear ) THEN
#ifdef CPP_PARA
             CALL histbeg(histname, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id,orch_domain_id)
#else
             CALL histbeg(histname, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id)
#endif
             WRITE(numout,*)  'HISTBEG --->',istp_old,date0,dt,dw,hist_id
          ELSE
#ifdef CPP_PARA
             CALL histbeg(histname, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id,domain_id=orch_domain_id)
#else
             CALL histbeg(histname, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id)
#endif
          ENDIF
          !-
          CALL histvert(hist_id, 'veget', 'Vegetation types', '1', &
               &    nvm,   veg, vegax_id)
          CALL histvert(hist_id, 'laiax', 'Nb LAI', '1', &
               &    nlai+1,   lai, laiax_id)
          CALL histvert(hist_id, 'solth', 'Soil levels',      'm', &
               &    ngrnd, sol, solax_id)
          CALL histvert(hist_id, 'soiltyp', 'Soil types',      '1', &
               &    nstm, soltyp, soltax_id)
          CALL histvert(hist_id, 'nobio', 'Other surface types',      '1', &
               &    nnobio, nobiotyp, nobioax_id)
          IF (  control_flags%hydrol_cwrr ) THEN
             CALL histvert(hist_id, 'solay', 'Hydrol soil levels',      'm', &
                  &    nslm, diaglev(1:nslm), solayax_id)
             CALL histvert(hist_id, 'soildiag', 'Diagnostic soil levels', 'm', &
                  &    nbdl, diaglev(1:nbdl), soildiagax_id)
          ENDIF

          CALL histvert(hist_id, 'snowlev', 'Snow levels',      'm', &
               &    nsnow, snowlev, snowax_id)

          !-
          !- SECHIBA_HISTLEVEL = 1
          !-
          CALL histdef(hist_id, 'evap', 'Evaporation', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(1), dt,dw)
          CALL histdef(hist_id, 'coastalflow', 'Diffuse coastal flow', 'm^3/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'riverflow', 'River flow to the oceans', 'm^3/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw) 
          CALL histdef(hist_id, 'temp_sol', 'Surface Temperature', 'C', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avecels(1), dt,dw)
          CALL histdef(hist_id, 'rain', 'Rainfall', 'mm/d',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(1), dt,dw)
          CALL histdef(hist_id, 'snowf', 'Snowfall', 'mm/d',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(1), dt,dw)
          CALL histdef(hist_id, 'netrad', 'Net radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'lai', 'Leaf Area Index', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(1), dt,dw)
          !
          IF (  control_flags%hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'reinf_slope', 'Slope index for each grid box', '1', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1),  dt,dw)
             CALL histdef(hist_id, 'soilindex', 'Soil index', '1', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, once(1),  dt,dw)
          ENDIF
          !
          IF ( control_flags%river_routing ) THEN
             CALL histdef(hist_id, 'basinmap', 'Aproximate map of the river basins', ' ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw) 
             CALL histdef(hist_id, 'nbrivers', 'Number or rivers in the outflow grid box', ' ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)  
          ENDIF
          !-
          !- SECHIBA_HISTLEVEL = 2
          !-
          CALL histdef(hist_id, 'subli', 'Sublimation', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
          CALL histdef(hist_id, 'evapnu', 'Bare soil evaporation', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
          CALL histdef(hist_id, 'runoff', 'Surface runoff', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
          CALL histdef(hist_id, 'drainage', 'Deep drainage', 'mm/d', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
          IF ( control_flags%river_routing ) THEN
             CALL histdef(hist_id, 'riversret', 'Return from endorheic rivers', 'mm/d', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw)
             CALL histdef(hist_id, 'hydrographs', 'Hydrographs of gridbox outflow', 'm^3/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(2), dt,dw)
          ENDIF
          IF ( control_flags%hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'evapnu_soil', 'Bare soil evap for soil type', 'mm/d', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, fluxop(2), dt,dw)
             CALL histdef(hist_id, 'drainage_soil', 'Drainage for soil type', 'mm/d', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, fluxop(2), dt,dw)
             CALL histdef(hist_id, 'transpir_soil', 'Transpir for soil type', 'mm/d', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, fluxop(2), dt,dw)
             CALL histdef(hist_id, 'runoff_soil', 'Runoff for soil type', 'mm/d', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, fluxop(2), dt,dw)
          ENDIF
          !
          CALL histdef(hist_id, 'tair', 'Air Temperature', 'K',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'qair', 'Air humidity', 'g/g',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          ! Ajouts Nathalie - Juillet 2006
          CALL histdef(hist_id, 'q2m', '2m Air humidity', 'g/g',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 't2m', '2m Air Temperature', 'K',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          ! Fin ajouts Nathalie
          CALL histdef(hist_id, 'alb_vis', 'Albedo visible', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'alb_nir', 'Albedo near infrared', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          ! Ajouts Nathalie - Septembre 2008
          CALL histdef(hist_id, 'soilalb_vis', 'Soil Albedo visible', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'soilalb_nir', 'Soil Albedo near infrared', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'vegalb_vis', 'Vegetation Albedo visible', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'vegalb_nir', 'Vegetation Albedo near infrared', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          ! Fin ajouts Nathalie - Septembre 2008
          CALL histdef(hist_id, 'z0', 'Surface roughness', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'roughheight', 'Effective roughness height', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(2), dt,dw)
          CALL histdef(hist_id, 'transpir', 'Transpiration', 'mm/d', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop(2), dt,dw)
          CALL histdef(hist_id, 'inter', 'Interception loss', 'mm/d', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop(2), dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 3
          !-
          CALL histdef(hist_id, 'tsol_max', 'Maximum Surface Temperature',&
               & 'C', iim,jjm, hori_id, 1,1,1, -99, 32, tmaxcels(3), dt,dw)
          CALL histdef(hist_id, 'tsol_min', 'Minimum Surface Temperature',&
               & 'C', iim,jjm, hori_id, 1,1,1, -99, 32, tmincels(3), dt,dw)
          CALL histdef(hist_id, 'fluxsens', 'Sensible Heat Flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(3), dt,dw)
          CALL histdef(hist_id, 'fluxlat', 'Latent Heat Flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(3), dt,dw)
          CALL histdef(hist_id, 'snow', 'Snow mass', 'kg/m^2', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'snowage', 'Snow age', '?', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'snownobio', 'Snow on other surfaces', 'kg/m^2', &
               & iim,jjm, hori_id, nnobio,1, nnobio, nobioax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'snownobioage', 'Snow age on other surfaces', 'd', &
               & iim,jjm, hori_id, nnobio,1, nnobio, nobioax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'vegetfrac', 'Fraction of vegetation', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'maxvegetfrac', 'Maximum fraction of vegetation', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'tot_bare_soil', "Total Bare Soil Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'nobiofrac', 'Fraction of other surface types', '1', &
               & iim,jjm, hori_id, nnobio, 1, nnobio, nobioax_id, 32, avescatter(3), dt,dw)
          IF ( control_flags%do_floodplains ) THEN
             CALL histdef(hist_id, 'flood_frac', 'Flooded fraction', '1', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(3), dt,dw)
             CALL histdef(hist_id, 'reinfiltration', 'Reinfiltration from floodplains', 'mm/d', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(3), dt,dw)
          ENDIF
          IF ( control_flags%hydrol_cwrr ) THEN
             DO jst=1,nstm
             
                ! var_name= "mc_1" ... "mc_3"
                WRITE (var_name,"('moistc_',i1)") jst
                CALL histdef(hist_id, var_name, 'Soil Moisture profile for soil type', 'm3/m3', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(3),  dt,dw)
                
                ! var_name= "vegetsoil_1" ... "vegetsoil_3"
                WRITE (var_name,"('vegetsoil_',i1)") jst
                CALL histdef(hist_id, var_name, 'Fraction of vegetation on soil types', '%', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3),  dt,dw)
                
                ! var_name= "kfact_root_1" ... "kfact_root_3"
                WRITE (var_name,"('kfactroot_',i1)") jst
                CALL histdef(hist_id, var_name, 'Root fraction profile for soil type', '%', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(3), dt,dw)
                
             ENDDO

             IF (ok_freeze_cwrr) THEN
                CALL histdef(hist_id, 'profil_froz_hydro', 'Frozen fraction for each hydrological soil layer', '-', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1),  dt,dw)
                DO jst=1,nstm
                   WRITE (var_name,"('profil_froz_hydro_',i1)") jst
                   CALL histdef(hist_id, trim(var_name), 'Frozen fraction for each hydrological soil layer and soiltile', '-', &
                        & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1),  dt,dw)
                ENDDO
                CALL histdef(hist_id, 'temp_hydro', 'Soil temperature interpolated on hydrological layers', 'K', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1), dt,dw)
                CALL histdef(hist_id, 'kk_moy', 'Mean hydrological conductivity', 'mm/d', &
                     & iim,jjm,hori_id, nslm,1,nslm, solayax_id, 32, avescatter(1),  dt,dw)
             END IF
             
             CALL histdef(hist_id, 'ptn_beg', 'Soil temperature from previous time step', 'K', &
                  & iim,jjm,hori_id, ngrnd,1,ngrnd, solax_id, 32, avescatter(1),  dt,dw)
             CALL histdef(hist_id, 'pcappa_supp', 'Additional heat capacity due to soil freezing for each soil layer', 'J/K', &
                  & iim,jjm,hori_id, ngrnd,1,ngrnd, solax_id, 32, avescatter(1),  dt,dw)
            DO jv = 1, nvm
               WRITE(part_str,'(I2)') jv
               IF (jv < 10) part_str(1:1) = '0'
               CALL histdef(hist_id,'shum_ngrnd_perma_'//part_str(1:LEN_TRIM(part_str)), 'Saturation degree on thethermal axes', '-', &
                    & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(6),  dt,dw)
            END DO

            DO jv = 1, nvm
               WRITE(part_str,'(I2)') jv
               IF (jv < 10) part_str(1:1) = '0'
               CALL histdef(hist_id,'shum_perma_long_'//part_str(1:LEN_TRIM(part_str)), &
                    & 'Long-term Saturation degree on the thermal axes', '-', &
                    & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32,avescatter(6),  dt,dw)
            END DO

            DO jv = 1, nvm
               WRITE(part_str,'(I2)') jv
               IF (jv < 10) part_str(1:1) = '0'
               CALL histdef(hist_id, 'wetdiag_'//part_str(1:LEN_TRIM(part_str)), 'Deep ground moisture', 'fraction', &
                    & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(6),  dt,dw)
            END DO
            !          CALL histdef(hist_id, 'wetdiag', 'Deep ground moisture',
            !          'fraction', &
            !               & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32,
            !               avescatter(6),  dt,dw)
            DO jv = 1, nvm
               WRITE(part_str,'(I2)') jv
               IF (jv < 10) part_str(1:1) = '0'
               CALL histdef(hist_id, 'wetdiaglong_'//part_str(1:LEN_TRIM(part_str)), 'Long-term deep ground moisture', 'fraction', & 
                    & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(6),  dt,dw)
            END DO


             CALL histdef(hist_id, 'shumdiag_perma', 'Saturation degree of the soil', '-', &
                  & iim,jjm,hori_id,nbdl,1,nbdl, soildiagax_id, 32, avescatter(1),  dt,dw)
             CALL histdef(hist_id, 'stempdiag', 'Soil temperature on the diagnostic axis', 'K', &
                  & iim,jjm,hori_id,nbdl,1,nbdl, soildiagax_id, 32, avescatter(1),  dt,dw)
          ENDIF
          !
          CALL histdef(hist_id, 'frac_bare', 'Bare soil fraction for each tile', '-', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'soiltile', 'Fraction of soil tiles', '%', &
               & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, avescatter(4),  dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 4
          !-
          IF ( .NOT. control_flags%hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'dss', 'Up-reservoir Height', 'm',  &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'gqsb', 'Upper Soil Moisture', 'Kg/m^2',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'bqsb', 'Lower Soil Moisture', 'Kg/m^2',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
          ELSE
             CALL histdef(hist_id, 'humtot', 'Total Soil Moisture', 'Kg/m^2', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'humtot_soil', 'Soil Moisture for soil type', 'Kg/m^2', &
                  & iim,jjm, hori_id, nstm, 1, nstm, soltax_id, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'SWI', 'Soil wetness index','-',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'njsc', 'Soil class used for hydrology', '-', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, once(4), dt,dw)

!pss:+
 !         IF ( TOPM_calcul ) THEN
             CALL histdef(hist_id, 'fsat', 'Fraction of soil saturated', '-',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'fwet', 'Fraction of wetland', '-',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'fwt1', 'Fraction of soil with wt [0,xcm]  & ', '-',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'fwt2', 'Fraction of soil with wt [xcm,ycm] ', '-',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'fwt3', 'Fraction of soil with wt [ycm,zcm] ', '-',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
             CALL histdef(hist_id, 'fwt4', 'Fraction of soil with wt [ucm,vcm] ', '-',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(4), dt,dw)
   
             CALL histdef(hist_id, 'ZMIN', 'MIN INDICE TOPO', '-', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
             CALL histdef(hist_id, 'ZMAX', 'MAX INDICE TOPO', '-', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
             CALL histdef(hist_id, 'ZMEAN', 'MEAN INDICE TOPO', '-', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
! 	         CALL histdef(hist_id, 'NB_PIXE', 'NB PIXELS AVC VALEUR', '-', &
!	              & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
             CALL histdef(hist_id, 'ZSTDT', 'STD INDICE TOPO', '-', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      
             CALL histdef(hist_id, 'ZSKEW', 'SKEWNESS INDICE TOPO', '-', &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(3), dt,dw)      

             CALL histdef(hist_id, 'drunoff_tot', 'Surface drunoff', 'mm/d', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(2), dt,dw) 
 !         ENDIF
!pss:-
          ENDIF
          CALL histdef(hist_id, 'qsintveg', 'Water on canopy', 'Kg/m^2', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(4), dt,dw)
          CALL histdef(hist_id, 'rstruct', 'Structural resistance', 's/m', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(4), dt,dw)
          IF ( control_flags%ok_co2 ) THEN
             CALL histdef(hist_id, 'gpp', 'Net assimilation of carbon by the vegetation', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
          ENDIF
          IF ( control_flags%ok_stomate ) THEN
             CALL histdef(hist_id, 'nee', 'Net Ecosystem Exchange', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
             CALL histdef(hist_id, 'maint_resp', 'Maintenance respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
             CALL histdef(hist_id, 'hetero_resp', 'Heterotrophic respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
             CALL histdef(hist_id, 'growth_resp', 'Growth respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(4), dt,dw)
             CALL histdef(hist_id, 'npp', 'Net Primary Production', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt, dw)
          ENDIF
          CALL histdef(hist_id, 'precisol', 'Throughfall', 'mm/d',  &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop(4), dt,dw)
          CALL histdef(hist_id, 'drysoil_frac', 'Fraction of visibly dry soil', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(4), dt,dw)
          CALL histdef(hist_id, 'evapot', 'Potential evaporation', 'mm/d',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(4), dt,dw)
          CALL histdef(hist_id, 'evapot_corr', 'Potential evaporation', 'mm/d',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(4), dt,dw)



          !-
          !- SECHIBA_HISTLEVEL = 5
          !-
          CALL histdef(hist_id, 'swnet', 'Net solar radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(5), dt,dw)
          CALL histdef(hist_id, 'swdown', 'Incident solar radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(5), dt,dw)
          CALL histdef(hist_id, 'lwdown', 'Absorbed downward longwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(5), dt,dw)
          CALL histdef(hist_id, 'lwnet', 'Net surface longwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(5), dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 6
          !-
           call histdef(hist_id, 'ptn_pftmean', 'Soil temperature, PFT-mean','K', &
                & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32,avescatter(6), dt,dw)
           call histdef(hist_id, 'ptn_snow_pftmean', 'Snow temperature,PFT-mean', 'K', &
                & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6), dt,dw)

           DO jv = 1, nvm
              IF (permafrost_veg_exists(jv)) THEN
                 WRITE(part_str,'(I2)') jv
                 IF (jv < 10) part_str(1:1) = '0'
                 CALL histdef(hist_id, 'ptn_'//part_str(1:LEN_TRIM(part_str)),'Soil temperature', 'K', &
                      & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32,avescatter(6),  dt,dw)
              END IF
           ENDDO

          CALL histdef(hist_id, 'snowmelt', 'snow melt', 'kg/m2', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(6), dt,dw)
          CALL histdef(hist_id, 'frac_snow_veg', 'snow fraction on vegeted area','-', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(6), dt,dw)
          CALL histdef(hist_id, 'frac_snow_nobio', 'snow fraction on non-vegeted area', '-', &
               & iim,jjm, hori_id, nnobio, 1,nnobio, nobioax_id, 32, avescatter(6), dt,dw)
          CALL histdef(hist_id, 'pgflux', 'extra energy used for melting top snow layer', '-', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(6), dt,dw)

          IF (ok_explicitsnow) THEN
             CALL histdef(hist_id, 'grndflux', 'ground heat flux', 'W/m2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(6), dt,dw)
             CALL histdef(hist_id, 'snowrho', 'Snow density profile', 'kg/m3', & 
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6), dt,dw)
             CALL histdef(hist_id, 'snowtemp','Snow temperature profile','K', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'soiltemp','Soil temperature profile','K', &
                  & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowdz','Snow depth profile','m', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowliq','Snow liquid content profile','m', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowgrain','Snow grain profile','m', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowheat','Snow Heat profile','J/m2', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(6),dt,dw)
             CALL histdef(hist_id, 'snowflx','Snow flux','W/m2', &
                  & iim,jjm, hori_id, 1, 1, 1, snowax_id, 32,avescatter(1),dt,dw)
            CALL histdef(hist_id, 'snowcap','Snow capacity','W/m2', &
                  & iim,jjm, hori_id, 1, 1, 1, snowax_id, 32,avescatter(1),dt,dw)
            CALL histdef(hist_id, 'temp_sol_add','surface temperature from fluxes','K', &
                  & iim,jjm, hori_id, 1, 1, 1, snowax_id, 32,avescatter(1),dt,dw)
            CALL histdef(hist_id, 'soilflx','Soil flux','W/m2', &
                  & iim,jjm, hori_id, 1, 1, 1, snowax_id, 32,avescatter(1),dt,dw)
            CALL histdef(hist_id, 'soilcap','Soil capacity','W/m2', &
                  & iim,jjm, hori_id, 1, 1, 1, snowax_id, 32,avescatter(1),dt,dw)
            CALL histdef(hist_id, 'cgrnd_snow','cgrnd for snow','-', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(1),dt,dw)
            CALL histdef(hist_id, 'dgrnd_snow','dgrnd for snow','-', &
                  & iim,jjm, hori_id, nsnow, 1, nsnow, snowax_id, 32,avescatter(1),dt,dw)

          END IF
          CALL histdef(hist_id, 'Qg', 'Ground heat flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)

         IF (control_flags%hydrol_cwrr) THEN
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef(hist_id,'profil_froz_'//part_str(1:LEN_TRIM(part_str)), 'Frozen fraction of the soil','-', &
                     & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32,avescatter(6),  dt,dw)
             END IF
          END DO

          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef(hist_id, 'pcapa_'//part_str(1:LEN_TRIM(part_str)),'Apparent heat capacity', 'J/m3/K', &
                     & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32,avescatter(6),  dt,dw)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef(hist_id, 'pkappa_'//part_str(1:LEN_TRIM(part_str)),'Soil thermal conductivity', 'W/m/K', &
                     & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32,avescatter(6),  dt,dw)
             END IF
          END DO

         ENDIF


          !-
          !- SECHIBA_HISTLEVEL = 7
          !-
          IF ( control_flags%river_routing ) THEN
             CALL histdef(hist_id, 'fastr', 'Fast flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(7), dt,dw)
             CALL histdef(hist_id, 'slowr', 'Slow flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(7), dt,dw)
             CALL histdef(hist_id, 'streamr', 'Stream flow reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(7), dt,dw)
             CALL histdef(hist_id, 'lakevol', 'Volume in lake reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(7), dt,dw)
             
             !-
             !- SECHIBA_HISTLEVEL = 8
             !-
             CALL histdef(hist_id, 'pondr', 'Ponds reservoir', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
             CALL histdef(hist_id, 'swampmap', 'Map of swamps', 'm^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(8), dt,dw)
             !
             IF ( control_flags%do_irrigation ) THEN
                CALL histdef(hist_id, 'irrigation', 'Net irrigation', 'mm/d', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(8), dt,dw)
                CALL histdef(hist_id, 'netirrig', 'Net irrigation requirement', 'mm/d', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(8), dt,dw)
                CALL histdef(hist_id, 'irrigmap', 'Map of irrigated surfaces', 'm^2', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, once(8), dt,dw)
             ENDIF
             IF ( control_flags%do_floodplains ) THEN
                CALL histdef(hist_id, 'floodmap', 'Map of floodplains', 'm^2', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, once(8), dt,dw)
                CALL histdef(hist_id, 'floodh', 'Floodplains height', 'mm', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
                CALL histdef(hist_id, 'floodr', 'Floodplains reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
                CALL histdef(hist_id, 'floodout', 'Flow out of floodplains', 'mm/d', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(8), dt,dw)
                CALL histdef(hist_id, 'evapflo', 'Floodplains evaporation', 'mm/d', &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop(8), dt,dw)
             ENDIF
             !
          ENDIF

          IF ( control_flags%hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'k_litt', 'Litter cond', 'mm/d', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          ENDIF
          CALL histdef(hist_id, 'beta', 'Beta Function', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'raero', 'Aerodynamic resistance', 's/m',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          ! Ajouts Nathalie - Novembre 2006
          CALL histdef(hist_id, 'cdrag', 'Drag coefficient for LE and SH', '?',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'Wind', 'Wind speed', 'm/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          ! Fin ajouts Nathalie
          !
          CALL histdef(hist_id, 'qsatt' , 'Surface saturated humidity', 'g/g', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'vbeta1', 'Beta for sublimation', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'vbeta4', 'Beta for bare soil', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          CALL histdef(hist_id, 'vbeta5', 'Beta for floodplains', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(8), dt,dw)
          IF ( control_flags%ok_co2 ) THEN
             CALL histdef(hist_id, 'gsmean', 'mean stomatal conductance', 'umol/m2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(8), dt,dw)
          ENDIF
          CALL histdef(hist_id, 'humrel', 'Soil moisture stress', '-',  &
               & iim,jjm, hori_id, nvm,1,nvm, vegax_id, 32, avescatter(8), dt,dw)
          !-
          !- SECHIBA_HISTLEVEL = 9
          !-
          !-
          !- SECHIBA_HISTLEVEL = 10
          !-
          IF ( control_flags%ok_co2 ) THEN
             CALL histdef(hist_id, 'cimean', 'Stomatal CO2 concentation', 'mmole/m2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
          ENDIF
          CALL histdef(hist_id, 'vbeta3', 'Beta for Transpiration', 'mm/d', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
          CALL histdef(hist_id, 'rveget', 'Canopy resistance', 's/m', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
          IF ( .NOT. control_flags%hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'rsol', 'Soil resistance', 's/m',  &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(10), dt,dw)
          ENDIF
          CALL histdef(hist_id,'vbeta2','Beta for Interception loss','mm/d', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)
          CALL histdef(hist_id, 'qsintmax', 'Maximum Interception capacity', 'Kg/m^2', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(10), dt,dw)

          !- SECHIBA_HISTLEVEL = 11
          !-

          IF ( .NOT. control_flags%hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'mrsos', "Moisture in Upper 0.1 m of Soil Column", "kg m-2", &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
             
             CALL histdef(hist_id, 'mrso', "Total Soil Moisture Content", "kg m-2", &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
             
             CALL histdef(hist_id, 'mrros', "Surface Runoff", "kg m-2 s-1", &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
             
             CALL histdef(hist_id, 'mrro', "Total Runoff", "kg m-2 s-1", &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
             
             CALL histdef(hist_id, 'prveg', "Precipitation onto Canopy", "kg m-2 s-1", &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)

          ENDIF


          CALL histdef(hist_id, 'evspsblveg', "Evaporation from Canopy", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
          
          CALL histdef(hist_id, 'evspsblsoi', "Water Evaporation from Soil", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
          
          CALL histdef(hist_id, 'tran', "Transpiration", "kg m-2 s-1", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(11), dt,dw)
          
          CALL histdef(hist_id, 'treeFrac', "Tree Cover Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'grassFrac', "Natural Grass Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'cropFrac', "Crop Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'baresoilFrac', "Bare Soil Fraction", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          CALL histdef(hist_id, 'residualFrac', &
               & "Fraction of Grid Cell that is Land but Neither Vegetation-Covered nor Bare Soil", "%", &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(11), dt,dw)
          
          IF ( control_flags%ok_inca ) THEN
             CALL histdef(hist_id, 'PAR', 'PAR', 'umol phot/m^2/s',  &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
             IF ( control_flags%ok_radcanopy ) THEN
                CALL histdef(hist_id, 'PARsun', 'Sunlit Leaf PAR', 'umol phot/m^2/s', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'PARsh', 'Shaded Leaf Area PAR', 'umol phot/m^2/s', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'laisun', 'Sunlit Leaf Area Index', '1', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'laish', 'Shaded Leaf Area Index', '1', &
                     & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'Fdf', 'Fdf', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                IF ( control_flags%ok_multilayer ) then 
                   CALL histdef(hist_id, 'PARsuntab', 'Sunlit Leaf PAR', 'umol phot/m^2/s', &
                        & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(11), dt,dw)
                   CALL histdef(hist_id, 'PARshtab', 'Shaded Leaf Area PAR', 'umol phot/m^2/s', &
                        & iim,jjm, hori_id, nlai+1, 1, nlai+1, laiax_id, 32, avescatter(11), dt,dw)
                ENDIF
                CALL histdef(hist_id, 'Sinang', 'Sinang', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'PARdf', 'PARdf', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'PARdr', 'PARdr', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'Trans', 'Trans', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'Day', 'Day', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
                CALL histdef(hist_id, 'Year_length', 'Year_length', '1',  &
                     & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
             END IF
             
             CALL histdef(hist_id, 'flx_fertil_no', 'flx_fertil_no', 'ngN/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'CRF', 'CRF', '1', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_co2_bbg_year', 'flx_co2_bbg_year', 'kgC/m^2/yr ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(11), dt,dw)  
             CALL histdef(hist_id, 'N_qt_WRICE_year', 'N_qt_WRICE_year', 'kgN/yr ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(11), dt,dw)  
             CALL histdef(hist_id, 'N_qt_OTHER_year', 'N_qt_OTHER_year', 'kgN/yr ', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(11), dt,dw)  
             CALL histdef(hist_id, 'ptnlev1', 'ptnlev1', 'K',  &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_iso', 'flx_iso', 'kgC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_mono', 'flx_mono', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_ORVOC', 'flx_ORVOC', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_MBO', 'flx_MBO', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_methanol', 'flx_methanol', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_acetone', 'flx_acetone', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_acetal', 'flx_acetal', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_formal', 'flx_formal', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_acetic', 'flx_acetic', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_formic', 'flx_formic', 'kgC/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_no_soil', 'flx_no_soil', 'ngN/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
             CALL histdef(hist_id, 'flx_no', 'flx_no', 'ngN/m^2/s',&
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(11), dt,dw)
          ENDIF

       ELSE 
          !-
          !- This is the ALMA convention output now
          !-
          !- 
          IF ( rectilinear ) THEN
#ifdef CPP_PARA
             CALL histbeg(histname, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id,orch_domain_id)
#else
             CALL histbeg(histname, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id)
#endif
          ELSE
#ifdef CPP_PARA
             CALL histbeg(histname, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id,domain_id=orch_domain_id)
#else
             CALL histbeg(histname, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id)
#endif
          ENDIF
          !-
          CALL histvert(hist_id, 'veget', 'Vegetation types', '1', &
               &    nvm,   veg, vegax_id)
          CALL histvert(hist_id, 'solth', 'Soil levels',      'm', &
               &    ngrnd, sol, solax_id)
          CALL histvert(hist_id, 'soiltyp', 'Soil types',      '1', &
               &    nstm, soltyp, soltax_id)
          CALL histvert(hist_id, 'nobio', 'Other surface types',      '1', &
               &    nnobio, nobiotyp, nobioax_id)
          IF (  control_flags%hydrol_cwrr ) THEN
             CALL histvert(hist_id, 'solay', 'Hydrol soil levels',      'm', &
                  &    nslm, diaglev(1:nslm), solayax_id)
          ENDIF
          !-
          !-  Vegetation
          !-
          CALL histdef(hist_id, 'vegetfrac', 'Fraction of vegetation', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'maxvegetfrac', 'Maximum fraction of vegetation', '1', &
               & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, avescatter(3), dt,dw)
          CALL histdef(hist_id, 'nobiofrac', 'Fraction of other surface types', '1', &
               & iim,jjm, hori_id, nnobio, 1, nnobio, nobioax_id, 32, avescatter(3), dt,dw)
          !- 
          !-  General energy balance
          !-
          CALL histdef(hist_id, 'Tair', 'Near surface air temperature at forcing level', 'K',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'Qair', 'Near surface specific humidity at forcing level', 'g/g',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'SWnet', 'Net shortwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(1), dt,dw)
          CALL histdef(hist_id, 'LWnet', 'Net longwave radiation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'Qh', 'Sensible heat flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(1), dt,dw)
          CALL histdef(hist_id, 'Qle', 'Latent heat flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(1), dt,dw)
          CALL histdef(hist_id, 'Qg', 'Ground heat flux', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'Qf', 'Energy of fusion', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(2), dt,dw)
          CALL histdef(hist_id, 'Qv', 'Energy of sublimation', 'W/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'DelSurfHeat', 'Change in surface layer heat', 'J/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, sumscatter(1), dt,dw)
          CALL histdef(hist_id, 'DelColdCont', 'Change in snow surface layer cold content', 'J/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, sumscatter(1), dt,dw)
          !-
          !- General water balance
          !-
          CALL histdef(hist_id, 'Snowf', 'Snowfall rate', 'kg/m^2/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Rainf', 'Rainfall rate', 'kg/m^2/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Evap', 'Total Evapotranspiration', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Qs', 'Surface runoff', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Qsb', 'Sub-surface runoff', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'Qsm', 'Snowmelt', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'DelSoilMoist', 'Change in soil moisture', 'kg/m^2',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, sumscatter(1), dt,dw)
          CALL histdef(hist_id, 'DelSurfStor', 'Change in Surface Water Storage','kg/m^2',&
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, sumscatter(1), dt,dw)
          CALL histdef(hist_id, 'DelIntercept', 'Change in interception storage', 'kg/m^2',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, sumscatter(1), dt,dw)
          CALL histdef(hist_id, 'DelSWE', 'Change in Snow Water Equivalent', 'kg/m^2',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, sumscatter(1), dt,dw)
          IF ( control_flags%do_irrigation ) THEN
             CALL histdef(hist_id, 'Qirrig', 'Irrigation', 'kg/m^2/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'Qirrig_req', 'Irrigation requirement', 'kg/m^2/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          ENDIF
          !-
          !- Surface state
          !-
          CALL histdef(hist_id, 'AvgSurfT', 'Average surface temperature', 'K', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(1), dt,dw)
          CALL histdef(hist_id, 'PotSurfT', 'Potential (Unstressed) surface temperature', 'K', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(1), dt,dw)
          CALL histdef(hist_id, 'RadT', 'Surface radiative temperature', 'K', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, ave(1), dt,dw)
          CALL histdef(hist_id, 'Albedo', 'Albedo', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SWI', 'Soil wetness index','1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SurfStor', 'Surface Water Storage','kg/m^2',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SWE', 'Snow Water Equivalent', 'kg/m^2', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          !!-
          !-  Sub-surface state
          !- 
          IF ( .NOT. control_flags%hydrol_cwrr ) THEN
             CALL histdef(hist_id, 'SoilMoist', '3D average layer soil moisture', 'kg/m^2',  &
                  & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(1), dt,dw)
          ELSE
             CALL histdef(hist_id, 'SoilMoist', '3D average layer soil moisture', 'kg/m^2',  &
                  & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1), dt,dw)

             IF (ok_freeze_cwrr) THEN
                CALL histdef(hist_id, 'profil_froz_hydro', 'Frozen fraction for each hydrological soil layer', '-', &
                     & iim,jjm, hori_id, nslm, 1, nslm,solayax_id, 32, avescatter(1),  dt,dw)
                DO jst=1,nstm
                   WRITE (var_name,"('profil_froz_hydro_',i1)") jst
                   CALL histdef(hist_id, trim(var_name), 'Frozen fraction for each hydrological soil layer and soiltile', '-', &
                        & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1),  dt,dw)
                ENDDO

                CALL histdef(hist_id, 'temp_hydro', 'Soil temperature interpolated on hydrological layers', 'K', &
                     & iim,jjm, hori_id, nslm, 1, nslm, solayax_id, 32, avescatter(1), dt,dw)
                CALL histdef(hist_id, 'kk_moy', 'Mean hydrological conductivity', 'mm/d', &
                     & iim,jjm,hori_id, nslm,1,nslm, solayax_id, 32, avescatter(1),  dt,dw)
             ENDIF
          END IF
          CALL histdef(hist_id, 'SoilWet', 'Total soil wetness', '-',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SoilTemp', 'Soil temperature profile', 'K', &
               & iim,jjm, hori_id, ngrnd, 1, ngrnd, solax_id, 32, avescatter(1),  dt,dw)
          !- 
          !-  Evaporation components
          !-
          CALL histdef(hist_id, 'PotEvap', 'Potential evapotranspiration', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'PotEvapOld', 'Potential evapotranspiration old method', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'ECanop', 'Interception evaporation', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'TVeg', 'Transpiration', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'ESoil', 'Bare soil evaporation', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'EWater', 'Open water evaporation', 'kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'RootMoist','Root zone soil water', 'kg/m^2',  &
               & iim,jjm, hori_id, 1, 1, 1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SubSnow','Snow sublimation','kg/m^2/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          CALL histdef(hist_id, 'ACond', 'Aerodynamic conductance', 'm/s',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          IF ( control_flags%do_floodplains ) THEN
             CALL histdef(hist_id, 'Qflood', 'Floodplain Evaporation', 'kg/m^2/s', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(1), dt,dw)
          ENDIF
          !-
          !- Surface turbulence
          !-
          CALL histdef(hist_id, 'Z0', 'Roughness height', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'EffectHeight', 'Effective displacement height (h-d)', 'm',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          !-
          !-
          !-  Cold Season Processes
          !-
          CALL histdef(hist_id, 'SnowFrac', 'Snow cover fraction', '1',  &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SAlbedo', 'Snow albedo', '1', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          CALL histdef(hist_id, 'SnowDepth', '3D snow depth', 'm', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          !-
          !- Hydrologic variables
          !-
          CALL histdef(hist_id, 'SwampMap', 'Map of swamp areas', 'm^2', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)
          CALL histdef(hist_id, 'Dis', 'Simulated River Discharge', 'm^3/s', &
               & iim,jjm, hori_id, 1,1,1, -99, 32, fluxop_scinsec(2), dt,dw)
          IF ( control_flags%do_irrigation ) THEN
             CALL histdef(hist_id, 'IrrigationMap', 'Map of irrigated areas', 'm^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)
          ENDIF
          IF ( control_flags%do_floodplains ) THEN
             CALL histdef(hist_id, 'FloodplainsMap', 'Map of flooded areas', 'm^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)
             CALL histdef(hist_id, 'FloodFrac', 'Floodplain Fraction', '-', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
          ENDIF
          CALL histdef(hist_id, 'humrel', 'Soil moisture stress', '-',  &
               & iim,jjm, hori_id, nvm,1,nvm, vegax_id, 32, avescatter(8), dt,dw)
          !-
          !-  The carbon budget
          !-
          IF ( control_flags%ok_co2 ) THEN
             CALL histdef(hist_id, 'GPP', 'Net assimilation of carbon by the vegetation', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
          ENDIF
          IF ( control_flags%ok_stomate ) THEN
             CALL histdef(hist_id, 'NEE', 'Net Ecosystem Exchange', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'maint_resp', 'Maintenance respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'hetero_resp', 'Heterotrophic respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'growth_resp', 'Growth respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
             CALL histdef(hist_id, 'npp', 'Net Primary Production', 'gC/m^2/s', &
                  & iim,jjm, hori_id, nvm, 1, nvm, vegax_id, 32, fluxop_scinsec(1), dt,dw)
          ENDIF
          !
      ENDIF
       !-
       !- Forcing and grid information
       !-
       CALL histdef(hist_id, 'LandPoints', 'Land Points', '1', &
            & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt,dw)  
       CALL histdef(hist_id, 'Areas', 'Mesh areas', 'm2', &
            & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt, dw)
       CALL histdef(hist_id, 'Contfrac', 'Continental fraction', '1', &
            & iim,jjm, hori_id, 1,1,1, -99, 32, once(1), dt, dw)
       !-
       ! Write the names of the pfts in the history files
       global_attribute="PFT_name"
       DO i=1,nvm
          WRITE(global_attribute(9:10),"(I2.2)") i
          CALL histglobal_attr(hist_id, global_attribute, PFT_name(i))
       ENDDO
       !-
       CALL histend(hist_id)
    ENDIF ! IF (is_omp_root)
 
    END IF !IF ( dw == 0 )
    !
    !
    ! Second SECHIBA hist file
    !
    !-
    !Config Key   = SECHIBA_HISTFILE2
    !Config Desc  = Flag to switch on histfile 2 for SECHIBA (hi-frequency ?)
    !Config If    = OK_SECHIBA
    !Config Def   = n
    !Config Help  = This Flag switch on the second SECHIBA writing for hi (or low) 
    !Config         frequency writing. This second output is optional and not written
    !Config         by default.
    !Config Units = [FLAG]
    !-
    ok_histfile2=.FALSE.
    CALL getin_p('SECHIBA_HISTFILE2', ok_histfile2)
    WRITE(numout,*) 'SECHIBA_HISTFILE2 ', ok_histfile2
    !
    !-
    !Config Key   = WRITE_STEP2
    !Config Desc  = Frequency in seconds at which to WRITE output
    !Config If    = SECHIBA_HISTFILE2
    !Config Def   = 1800.0
    !Config Help  = This variables gives the frequency the output 2 of
    !Config         the model should be written into the netCDF file.
    !Config         It does not affect the frequency at which the
    !Config         operations such as averaging are done.
    !Config         That is IF the coding of the calls to histdef
    !Config         are correct !
    !Config Units = [seconds]
    !-
    dw2 = 1800.0
    CALL getin_p('WRITE_STEP2', dw2)
    
    ! Deactivate sechiba_histfile2 if the frequency is set to zero
    IF ( dw2 == 0 ) THEN
       ok_histfile2=.FALSE.
       WRITE(numout,*) 'WRITE_STEP2 was set to zero and therfore SECHIBA_HISTFILE2 is deactivated.'
    ELSE IF ( hist_id < 0 ) THEN
       ! Deactivate all history files if sechiba_history file is deactivated
       ok_histfile2=.FALSE.
       WRITE(numout,*) 'SECHIBA_HISTFILE2 will not be created because sechiba_history file is deactivated.'
    END IF

    hist2_id = -1
    !
    IF (ok_histfile2) THEN
       !-
       !Config Key   = SECHIBA_OUTPUT_FILE2
       !Config Desc  = Name of file in which the output number 2 is going to be written
       !Config If    = SECHIBA_HISTFILE2
       !Config Def   = sechiba_out_2.nc
       !Config Help  = This file is going to be created by the model
       !Config         and will contain the output 2 from the model.
       !Config Units = [FILE]
       !-
       histname2='sechiba_out_2.nc'
       CALL getin_p('SECHIBA_OUTPUT_FILE2', histname2)
       WRITE(numout,*) 'SECHIBA_OUTPUT_FILE2 ', histname2
       !-
       !Config Key   = SECHIBA_HISTLEVEL2
       !Config Desc  = SECHIBA history 2 output level (0..10)
       !Config If    = SECHIBA_HISTFILE2
       !Config Def   = 1
       !Config Help  = Chooses the list of variables in the history file. 
       !Config         Values between 0: nothing is written; 10: everything is 
       !Config         written are available More details can be found on the web under documentation.
       !Config         web under documentation.
       !Config         First level contains all ORCHIDEE outputs.
       !Config Units = [-] 
       !-
       hist2_level = 1
       CALL getin_p('SECHIBA_HISTLEVEL2', hist2_level)
       !-
       WRITE(numout,*) 'SECHIBA history level 2 : ',hist2_level
       IF ( (hist2_level > max_hist_level).OR.(hist2_level < 0) ) THEN
          STOP 'This history level 2 is not allowed'
       ENDIF
       !
       !-
       !- define operations as a function of history level.
       !- Above hist2_level, operation='never'
       !-
       ave2(1:max_hist_level) = 'ave(scatter(X))'
       IF (hist2_level < max_hist_level) THEN
          ave2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       sumscatter2(1:max_hist_level) = 't_sum(scatter(X))'
       IF (hist2_level < max_hist_level) THEN
          sumscatter2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       avecels2(1:max_hist_level) = 'ave(cels(scatter(X)))'
       IF (hist2_level < max_hist_level) THEN
          avecels2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       avescatter2(1:max_hist_level) = 'ave(scatter(X))'
       IF (hist2_level < max_hist_level) THEN
          avescatter2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       tmincels2(1:max_hist_level) = 't_min(cels(scatter(X)))'
       IF (hist2_level < max_hist_level) THEN
          tmincels2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       tmaxcels2(1:max_hist_level) = 't_max(cels(scatter(X)))'
       IF (hist2_level < max_hist_level) THEN
          tmaxcels2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
!!$       tmax2(1:max_hist_level) = 't_max(X)'
!!$       IF (hist2_level < max_hist_level) THEN
!!$          tmax2(hist2_level+1:max_hist_level) = 'never'
!!$       ENDIF
       fluxop2(1:max_hist_level) = flux_op
       IF (hist2_level < max_hist_level) THEN
          fluxop2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
!!$       fluxop_sc2(1:max_hist_level) = flux_sc
!!$       IF (hist2_level < max_hist_level) THEN
!!$          fluxop_sc2(hist2_level+1:max_hist_level) = 'never'
!!$       ENDIF
!!$       fluxop_insec2(1:max_hist_level) = flux_insec
!!$       IF (hist2_level < max_hist_level) THEN
!!$          fluxop_insec2(hist2_level+1:max_hist_level) = 'never'
!!$       ENDIF
       fluxop_scinsec2(1:max_hist_level) = flux_scinsec
       IF (hist2_level < max_hist_level) THEN
          fluxop_scinsec2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       once2(1:max_hist_level) = 'once(scatter(X))'
       IF (hist2_level < max_hist_level) THEN
          once2(hist2_level+1:max_hist_level) = 'never'
       ENDIF
       ! 
       IF (is_omp_root) THEN
          IF ( .NOT. almaoutput ) THEN
             !- 
             IF ( rectilinear ) THEN
#ifdef CPP_PARA
                CALL histbeg(histname2, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id,orch_domain_id)
#else
                CALL histbeg(histname2, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id)
#endif
                WRITE(numout,*)  'HISTBEG2 --->',istp_old,date0,dt,dw2,hist2_id
             ELSE
#ifdef CPP_PARA
                CALL histbeg(histname2, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id,domain_id=orch_domain_id)
#else
                CALL histbeg(histname2, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id)
#endif
             ENDIF
             !-
             CALL histvert(hist2_id, 'veget', 'Vegetation types', '1', &
                  &    nvm,   veg, vegax_id2)
             CALL histvert(hist2_id, 'laiax', 'Nb LAI', '1', &
                  &    nlai+1,   lai, laiax_id2)
             CALL histvert(hist2_id, 'solth', 'Soil levels',      'm', &
                  &    ngrnd, sol, solax_id2)
             CALL histvert(hist2_id, 'soiltyp', 'Soil types',      '1', &
                  &    nstm, soltyp, soltax_id2)
             CALL histvert(hist2_id, 'nobio', 'Other surface types',      '1', &
                  &    nnobio, nobiotyp, nobioax_id2)
             CALL histvert(hist2_id, 'albtyp', 'Albedo Types',     '1', &
                  &    2, albtyp, albax_id2)
             IF (  control_flags%hydrol_cwrr ) THEN
                CALL histvert(hist2_id, 'solay', 'Hydrol soil levels',      'm', &
                     &    nslm, solay, solayax_id2)
             ENDIF
             !-
             !- SECHIBA_HISTLEVEL2 = 1
             !-
             CALL histdef(hist2_id, 'ptn', 'Deep ground temperature', 'K', &
                  & iim,jjm, hori_id2, ngrnd, 1, ngrnd, solax_id2, 32, avescatter2(1),  dt, dw2)
             IF ( .NOT. control_flags%hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'mrsos', "Moisture in Upper 0.1 m of Soil Column", "kg m-2", &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(1), dt,dw2)
                
                CALL histdef(hist2_id, 'mrro', "Total Runoff", "kg m-2 s-1", &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(1), dt,dw2)
             ENDIF
             !-
             !- SECHIBA_HISTLEVEL2 = 2
             !-
             CALL histdef(hist2_id, 'cdrag', 'Drag coefficient for LE and SH', '?',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt, dw2)
             ! Ajouts Nathalie - Septembre 2008
             CALL histdef(hist2_id, 'soilalb_vis', 'Soil Albedo visible', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
             CALL histdef(hist2_id, 'soilalb_nir', 'Soil Albedo near infrared', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
             CALL histdef(hist2_id, 'vegalb_vis', 'Vegetation Albedo visible', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
             CALL histdef(hist2_id, 'vegalb_nir', 'Vegetation Albedo near infrared', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
             ! Fin ajouts Nathalie - Septembre 2008
             CALL histdef(hist2_id, 'z0', 'Surface roughness', 'm',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt, dw2)
             CALL histdef(hist2_id, 'coastalflow', 'Diffuse coastal flow', 'm^3/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(2), dt, dw2)
             CALL histdef(hist2_id, 'riverflow', 'River flow to the oceans', 'm^3/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(2), dt, dw2) 
             CALL histdef(hist2_id, 'tsol_rad', 'Radiative surface temperature', 'C', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avecels2(2), dt, dw2)
             IF ( control_flags%do_floodplains ) THEN
                CALL histdef(hist2_id, 'floodout', 'Flow out of floodplains', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(2), dt,dw2)
                CALL histdef(hist2_id, 'vevapflo', 'Floodplains evaporation', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(2), dt, dw2)
                CALL histdef(hist2_id, 'flood_frac', 'Flooded fraction', '1', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(2), dt,dw2)
                CALL histdef(hist2_id, 'reinfiltration', 'Reinfiltration from floodplains', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(2), dt,dw2)
             ENDIF
             CALL histdef(hist2_id, 'vevapnu', 'Bare soil evaporation', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(2), dt, dw2)
             CALL histdef(hist2_id, 'temp_sol', 'New Surface Temperature', 'C', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avecels2(2), dt, dw2)
             CALL histdef(hist2_id, 'qsurf', 'Near surface specific humidity', 'g/g',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(2), dt, dw2)
             CALL histdef(hist2_id, 'albedo', 'Albedo', '1', &
                  & iim,jjm, hori_id2, 2,1,2, albax_id2, 32, avescatter2(2), dt, dw2)
             CALL histdef(hist2_id, 'fluxsens', 'Sensible Heat Flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(2), dt, dw2)
             CALL histdef(hist2_id, 'fluxlat', 'Latent Heat Flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(2), dt, dw2)
             CALL histdef(hist2_id, 'emis', 'Surface emissivity', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(2), dt, dw2)
             !-
             !- SECHIBA_HISTLEVEL2 = 3
             !-
             CALL histdef(hist2_id, 'evap', 'Evaporation', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(3), dt, dw2)
             CALL histdef(hist2_id, 'rain', 'Rainfall', 'mm/d',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(3), dt, dw2)
             CALL histdef(hist2_id, 'snowf', 'Snowfall', 'mm/d',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(3), dt, dw2)
             CALL histdef(hist2_id, 'netrad', 'Net radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(3), dt, dw2)
             CALL histdef(hist2_id, 'lai', 'Leaf Area Index', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(3), dt, dw2)
             IF ( control_flags%river_routing ) THEN
                CALL histdef(hist2_id, 'basinmap', 'Aproximate map of the river basins', ' ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(3), dt, dw2) 
                CALL histdef(hist2_id, 'nbrivers', 'Number or rivers in the outflow grid box', ' ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(3), dt, dw2)  
             ENDIF
             IF (check_waterbal) THEN
                CALL histdef(hist2_id, 'TotWater', 'Total amount of water at end of time step', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(3), dt, dw2)
                CALL histdef(hist2_id, 'TotWaterFlux', 'Total water flux', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(3), dt, dw2)
             ENDIF
             !-
             !- SECHIBA_HISTLEVEL2 = 4
             !-
             CALL histdef(hist2_id, 'subli', 'Sublimation', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(4), dt, dw2)
             CALL histdef(hist2_id, 'runoff', 'Surface runoff', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(4), dt, dw2)
             CALL histdef(hist2_id, 'drainage', 'Deep drainage', 'mm/d', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(4), dt, dw2)
             IF ( control_flags%river_routing ) THEN
                CALL histdef(hist2_id, 'riversret', 'Return from endorheic rivers', 'mm/d', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(4), dt, dw2)
                CALL histdef(hist2_id, 'hydrographs', 'Hydrographs of gridbox outflow', 'm^3/s', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(4), dt, dw2)
             ENDIF
             IF ( control_flags%hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'evapnu_soil', 'Bare soil evap for soil type', 'mm/d', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, fluxop2(4), dt, dw2)
                CALL histdef(hist2_id, 'drainage_soil', 'Drainage for soil type', 'mm/d', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, fluxop2(4), dt, dw2)
                CALL histdef(hist2_id, 'transpir_soil', 'Transpir for soil type', 'mm/d', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, fluxop2(4), dt, dw2)
                CALL histdef(hist2_id, 'runoff_soil', 'Runoff for soil type', 'mm/d', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, fluxop2(4), dt, dw2)
             ENDIF
             !
             CALL histdef(hist2_id, 'tair', 'Air Temperature', 'K',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 'qair', 'Air humidity', 'g/g',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             ! Ajouts Nathalie - Juillet 2006
             CALL histdef(hist2_id, 'q2m', '2m Air humidity', 'g/g',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 't2m', '2m Air Temperature', 'K',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             ! Fin ajouts Nathalie
             CALL histdef(hist2_id, 'alb_vis', 'Albedo visible', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 'alb_nir', 'Albedo near infrared', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(4), dt, dw2)
             CALL histdef(hist2_id, 'roughheight', 'Effective roughness height', 'm',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(4), dt, dw2)
             CALL histdef(hist2_id, 'transpir', 'Transpiration', 'mm/d', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop2(4), dt, dw2)
             CALL histdef(hist2_id, 'inter', 'Interception loss', 'mm/d', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop2(4), dt, dw2)
             !-
             !- SECHIBA_HISTLEVEL2 = 5
             !-
             CALL histdef(hist2_id, 'tsol_max', 'Maximum Surface Temperature',&
                  & 'C', iim,jjm, hori_id2, 1,1,1, -99, 32, tmaxcels2(5), dt, dw2)
             CALL histdef(hist2_id, 'tsol_min', 'Minimum Surface Temperature',&
                  & 'C', iim,jjm, hori_id2, 1,1,1, -99, 32, tmincels2(5), dt, dw2)
             CALL histdef(hist2_id, 'snow', 'Snow mass', 'kg/m^2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'snowage', 'Snow age', '?', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'snownobio', 'Snow on other surfaces', 'kg/m^2', &
                  & iim,jjm, hori_id2, nnobio,1, nnobio, nobioax_id2, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'snownobioage', 'Snow age on other surfaces', 'd', &
                  & iim,jjm, hori_id2, nnobio,1, nnobio, nobioax_id2, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'vegetfrac', 'Fraction of vegetation', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'maxvegetfrac', 'Maximum fraction of vegetation', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'nobiofrac', 'Fraction of other surface types', '1', &
                  & iim,jjm, hori_id2, nnobio, 1, nnobio, nobioax_id2, 32, avescatter2(5), dt, dw2)
             IF ( control_flags%hydrol_cwrr ) THEN
                DO jst=1,nstm
                   
                   ! var_name= "mc_1" ... "mc_3"
                   WRITE (var_name,"('moistc_',i1)") jst
                   CALL histdef(hist2_id, var_name, 'Soil Moisture profile for soil type', 'm3/m3', &
                        & iim,jjm, hori_id2, nslm, 1, nslm, solayax_id2, 32, avescatter2(5), dt, dw2)
                   
                   ! var_name= "vegetsoil_1" ... "vegetsoil_3"
                   WRITE (var_name,"('vegetsoil_',i1)") jst
                   CALL histdef(hist2_id, var_name, 'Fraction of vegetation on soil types', '%', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(5), dt, dw2)
                   
                   ! var_name= "kfact_root_1" ... "kfact_root_3"
                   WRITE (var_name,"('kfactroot_',i1)") jst
                   CALL histdef(hist2_id, var_name, 'Root fraction profile for soil type', '%', &
                        & iim,jjm, hori_id2, nslm, 1, nslm, solayax_id2, 32, avescatter2(5), dt,dw2)
                ENDDO

             ENDIF
             !-
             !- SECHIBA_HISTLEVEL2 = 6
             !-
             IF ( .NOT. control_flags%hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'dss', 'Up-reservoir Height', 'm',  &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id, 32, avescatter2(6), dt,dw2)
                CALL histdef(hist2_id, 'gqsb', 'Upper Soil Moisture', 'Kg/m^2',  &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(6), dt, dw2)
                CALL histdef(hist2_id, 'bqsb', 'Lower Soil Moisture', 'Kg/m^2',  &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(6), dt, dw2)
             ELSE
                CALL histdef(hist2_id, 'humtot', 'Total Soil Moisture', 'Kg/m^2', &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(6), dt, dw2)
                CALL histdef(hist2_id, 'humtot_soil', 'Soil Moisture for soil type', 'Kg/m^2', &
                     & iim,jjm, hori_id2, nstm, 1, nstm, soltax_id2, 32, avescatter2(6), dt, dw2)
                CALL histdef(hist2_id, 'SWI', 'Soil wetness index','-',  &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(6), dt,dw2)
             ENDIF
             CALL histdef(hist2_id, 'qsintveg', 'Water on canopy', 'Kg/m^2', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(6), dt, dw2)
             CALL histdef(hist2_id, 'rstruct', 'Structural resistance', 's/m', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(6), dt, dw2)
             IF ( control_flags%ok_co2 ) THEN
                CALL histdef(hist2_id, 'gpp', 'Net assimilation of carbon by the vegetation', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt,dw2)
             ENDIF
             IF ( control_flags%ok_stomate ) THEN
                CALL histdef(hist2_id, 'nee', 'Net Ecosystem Exchange', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt,dw2)
                CALL histdef(hist2_id, 'maint_resp', 'Maintenance respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt,dw2)
                CALL histdef(hist2_id, 'hetero_resp', 'Heterotrophic respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt,dw2)
                CALL histdef(hist2_id, 'growth_resp', 'Growth respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt, dw2)
                CALL histdef(hist2_id, 'npp', 'Net Primary Production', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(6), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'precisol', 'Throughfall', 'mm/d',  &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop2(6), dt, dw2)
             CALL histdef(hist2_id, 'drysoil_frac', 'Fraction of visibly dry soil', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(6), dt, dw2)
             CALL histdef(hist2_id, 'evapot', 'Potential evaporation', 'mm/d',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(6), dt, dw2)
             CALL histdef(hist2_id, 'evapot_corr', 'Potential evaporation', 'mm/d',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(6), dt, dw2)
             CALL histdef(hist2_id, 'snowmelt', 'snow melt', 'kg/m2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter(6), dt, dw2)

             !-
             !- SECHIBA_HISTLEVEL2 = 7
             !-
             CALL histdef(hist2_id, 'swnet', 'Net solar radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(7), dt, dw2)
             CALL histdef(hist2_id, 'swdown', 'Incident solar radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(7), dt, dw2)
             CALL histdef(hist2_id, 'lwdown', 'Absorbed downward longwave radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'lwnet', 'Net surface longwave radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             !-
             !- SECHIBA_HISTLEVEL2 = 8
             !-
             IF ( control_flags%river_routing ) THEN
                CALL histdef(hist2_id, 'fastr', 'Fast flow reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt, dw2)
                CALL histdef(hist2_id, 'slowr', 'Slow flow reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt, dw2)
                CALL histdef(hist2_id, 'streamr', 'Stream flow reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt, dw2)
                CALL histdef(hist2_id, 'floodr', 'Floodplains reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt,dw2)
                CALL histdef(hist2_id, 'floodh', 'Floodplains height', 'mm', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt,dw2)
                CALL histdef(hist2_id, 'pondr', 'Ponds reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt,dw2)
                CALL histdef(hist2_id, 'lakevol', 'Volume in lake reservoir', 'kg/m^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(8), dt, dw2)
                IF ( control_flags%do_irrigation ) THEN
                   CALL histdef(hist2_id, 'irrigation', 'Net irrigation', 'mm/d', &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(8), dt, dw2)
                   CALL histdef(hist2_id, 'netirrig', 'Net irrigation requirement', 'mm/d', &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop2(8), dt, dw2)
                   CALL histdef(hist2_id, 'irrigmap', 'Map of irrigated areas', 'm^2', &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(8), dt, dw2)
                ENDIF
                CALL histdef(hist2_id, 'floodmap', 'Map of floodplains', 'm^2', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(8), dt,dw2)
                CALL histdef(hist2_id, 'swampmap', 'Map of swamps', 'm^2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(8), dt,dw2)
             ENDIF
             !-
             !- SECHIBA_HISTLEVEL2 = 9
             !-
             CALL histdef(hist2_id, 'beta', 'Beta Function', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             CALL histdef(hist2_id, 'raero', 'Aerodynamic resistance', 's/m',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             ! Ajouts Nathalie - Novembre 2006
             CALL histdef(hist2_id, 'Wind', 'Wind speed', 'm/s',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             ! Fin ajouts Nathalie
             CALL histdef(hist2_id, 'qsatt' , 'Surface saturated humidity', 'g/g', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             CALL histdef(hist2_id, 'vbeta1', 'Beta for sublimation', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             CALL histdef(hist2_id, 'vbeta4', 'Beta for bare soil', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             IF ( control_flags%ok_co2 ) THEN
                CALL histdef(hist2_id, 'gsmean', 'mean stomatal conductance', 'umol/m2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(9), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'vbeta5', 'Beta for floodplains', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(9), dt, dw2)
             IF (  control_flags%hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'reinf_slope', 'Slope index for each grid box', '1', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(9),  dt,dw2)
                CALL histdef(hist2_id, 'soilindex', 'Soil index', '1', &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, once2(9),  dt,dw2)
             ENDIF
             CALL histdef(hist2_id, 'humrel', 'Soil moisture stress', '-',  &
                  & iim,jjm, hori_id2, nvm,1,nvm, vegax_id2, 32, avescatter2(9), dt, dw2)
             !-
             !- SECHIBA_HISTLEVEL2 = 10
             !-
             IF ( control_flags%ok_co2 ) THEN
                CALL histdef(hist2_id, 'cimean', 'Stomatal CO2 concentation', 'mmole/m2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'vbeta3', 'Beta for Transpiration', 'mm/d', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             CALL histdef(hist2_id, 'rveget', 'Canopy resistance', 's/m', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             IF ( .NOT. control_flags%hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'rsol', 'Soil resistance', 's/m',  &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt, dw2)
             ENDIF
             CALL histdef(hist2_id,'vbeta2','Beta for Interception loss','mm/d', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             CALL histdef(hist2_id, 'qsintmax', 'Maximum Interception capacity', 'Kg/m^2', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt, dw2)
             ! 
             IF ( control_flags%ok_inca ) THEN
                CALL histdef(hist2_id, 'PAR', 'PAR', 'umol phot/m^2/s',  &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                IF ( control_flags%ok_radcanopy ) THEN
                   CALL histdef(hist2_id, 'PARsun', 'Sunlit Leaf PAR', 'umol phot/m^2/s', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'PARsh', 'Shaded Leaf Area PAR', 'umol phot/m^2/s', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'laisun', 'Sunlit Leaf Area Index', '1', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'laish', 'Shaded Leaf Area Index', '1', &
                        & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'Fdf', 'Fdf', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   IF ( control_flags%ok_multilayer ) then 
                      CALL histdef(hist2_id, 'PARsuntab', 'Sunlit Leaf PAR', 'umol phot/m^2/s', &
                           & iim,jjm, hori_id2, nlai+1, 1, nlai+1, laiax_id2, 32, avescatter2(10), dt,dw2)
                      CALL histdef(hist2_id, 'PARshtab', 'Shaded Leaf Area PAR', 'umol phot/m^2/s', &
                           & iim,jjm, hori_id2, nlai+1, 1, nlai+1, laiax_id2, 32, avescatter2(10), dt,dw2)
                   ENDIF
                   CALL histdef(hist2_id, 'Sinang', 'Sinang', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'PARdf', 'PARdf', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'PARdr', 'PARdr', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'Trans', 'Trans', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'Day', 'Day', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                   CALL histdef(hist2_id, 'Year_length', 'Year_length', '1',  &
                        & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                END IF
                
                CALL histdef(hist2_id, 'flx_fertil_no', 'flx_fertil_no', 'ngN/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'CRF', 'CRF', '1', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_co2_bbg_year', 'flx_co2_bbg_year', 'kgC/m^2/yr ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(10), dt,dw2)  
                CALL histdef(hist2_id, 'N_qt_WRICE_year', 'N_qt_WRICE_year', 'kgN/yr ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(10), dt,dw2)  
                CALL histdef(hist2_id, 'N_qt_OTHER_year', 'N_qt_OTHER_year', 'kgN/yr ', &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(10), dt,dw2)  
                CALL histdef(hist2_id, 'ptnlev1', 'ptnlev1', 'K',  &
                     & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_iso', 'flx_iso', 'kgC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_mono', 'flx_mono', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_ORVOC', 'flx_ORVOC', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_MBO', 'flx_MBO', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_methanol', 'flx_methanol', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_acetone', 'flx_acetone', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_acetal', 'flx_acetal', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_formal', 'flx_formal', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_acetic', 'flx_acetic', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_formic', 'flx_formic', 'kgC/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_no_soil', 'flx_no_soil', 'ngN/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
                CALL histdef(hist2_id, 'flx_no', 'flx_no', 'ngN/m^2/s',&
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(10), dt,dw2)
             ENDIF
         ELSE 
             !-
             !- This is the ALMA convention output now
             !-
             !- 
             IF ( rectilinear ) THEN
#ifdef CPP_PARA
                CALL histbeg(histname2, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id,orch_domain_id)
#else
                CALL histbeg(histname2, iim, lon_rect, jjm, lat_rect, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id)
#endif
                WRITE(numout,*)  'HISTBEG2 --->',istp_old,date0,dt,dw2,hist2_id
             ELSE
#ifdef CPP_PARA
                CALL histbeg(histname2, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id,domain_id=orch_domain_id)
#else
                CALL histbeg(histname2, iim, lon, jjm, lat, 1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id2, hist2_id)
#endif 
             ENDIF
             !-
             CALL histvert(hist2_id, 'veget', 'Vegetation types', '1', &
                  &    nvm,   veg, vegax_id2)
             CALL histvert(hist2_id, 'solth', 'Soil levels',      'm', &
                  &    ngrnd, sol, solax_id2)
             CALL histvert(hist2_id, 'soiltyp', 'Soil types',      '1', &
                  &    nstm, soltyp, soltax_id2)
             CALL histvert(hist2_id, 'nobio', 'Other surface types',      '1', &
                  &    nnobio, nobiotyp, nobioax_id2)
             IF (  control_flags%hydrol_cwrr ) THEN
                CALL histvert(hist2_id, 'solay', 'Hydrol soil levels',      'm', &
                     &    nslm, diaglev(1:nslm), solayax_id2)
             ENDIF
             !-
             !-  Vegetation
             !-
             CALL histdef(hist2_id, 'vegetfrac', 'Fraction of vegetation', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(3), dt, dw2)
             CALL histdef(hist2_id, 'maxvegetfrac', 'Maximum fraction of vegetation', '1', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, avescatter2(3), dt, dw2)
             CALL histdef(hist2_id, 'nobiofrac', 'Fraction of other surface types', '1', &
                  & iim,jjm, hori_id2, nnobio, 1, nnobio, nobioax_id2, 32, avescatter2(3), dt, dw2)
             !- 
             !-  General energy balance
             !-
             CALL histdef(hist2_id, 'SWnet', 'Net shortwave radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(1), dt, dw2)
             CALL histdef(hist2_id, 'LWnet', 'Net longwave radiation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qh', 'Sensible heat flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qle', 'Latent heat flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qg', 'Ground heat flux', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qf', 'Energy of fusion', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(5), dt, dw2)
             CALL histdef(hist2_id, 'Qv', 'Energy of sublimation', 'W/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             CALL histdef(hist2_id, 'DelSurfHeat', 'Change in surface layer heat', 'J/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, sumscatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'DelColdCont', 'Change in snow surface layer cold content', 'J/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, sumscatter2(7), dt, dw2)
             !-
             !- General water balance
             !-
             CALL histdef(hist2_id, 'Snowf', 'Snowfall rate', 'kg/m^2/s',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Rainf', 'Rainfall rate', 'kg/m^2/s',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Evap', 'Total Evapotranspiration', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qs', 'Surface runoff', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qsb', 'Sub-surface runoff', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'Qsm', 'Snowmelt', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'DelSoilMoist', 'Change in soil moisture', 'kg/m^2',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, sumscatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'DelSurfStor', 'Change in Surface Water Storage','kg/m^2',&
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, sumscatter2(7), dt,dw2)
             CALL histdef(hist2_id, 'DelIntercept', 'Change in interception storage', 'kg/m^2',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, sumscatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'DelSWE', 'Change in interception storage Snow Water Equivalent', 'kg/m^2',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, sumscatter2(7), dt, dw2)
             !-
             !- Surface state
             !-
             CALL histdef(hist2_id, 'AvgSurfT', 'Average surface temperature', 'K', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(1), dt, dw2)
             CALL histdef(hist2_id, 'RadT', 'Surface radiative temperature', 'K', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, ave2(5), dt, dw2)
             CALL histdef(hist2_id, 'Albedo', 'Albedo', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'SWI', 'Soil wetness index','1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'SurfStor', 'Surface Water Storage','kg/m^2',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'SWE', 'Snow Water Equivalent', 'kg/m^2', &
                  & iim,jjm, hori_id, 1,1,1, -99, 32, avescatter(1), dt,dw)
             !!-
             !-  Sub-surface state
             !- 
             IF ( .NOT. control_flags%hydrol_cwrr ) THEN
                CALL histdef(hist2_id, 'SoilMoist', '3D average layer soil moisture', 'kg/m^2',  &
                     & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(7), dt, dw2)
             ELSE
                CALL histdef(hist2_id, 'SoilMoist', '3D average layer soil moisture', 'kg/m^2',  &
                     & iim,jjm, hori_id2, nslm, 1, nslm, solayax_id2, 32, avescatter2(7), dt, dw2)
             ENDIF
             CALL histdef(hist2_id, 'SoilWet', 'Total soil wetness', '-',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'SoilTemp', 'Soil temperature profile', 'K', &
                  & iim,jjm, hori_id2, ngrnd, 1, ngrnd, solax_id2, 32, avescatter2(7), dt, dw2)
             !- 
             !-  Evaporation components
             !-
             CALL histdef(hist2_id, 'PotEvap', 'Potential evapotranspiration', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'ECanop', 'Interception evaporation', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(5), dt, dw2)
             CALL histdef(hist2_id, 'TVeg', 'Transpiration', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, fluxop_scinsec2(5), dt, dw2)
             CALL histdef(hist2_id, 'ESoil', 'Bare soil evaporation', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(5), dt, dw2)
             CALL histdef(hist2_id, 'EWater', 'Open water evaporation', 'kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(5), dt, dw2)
             CALL histdef(hist2_id, 'RootMoist','Root zone soil water', 'kg/m^2',  &
                  & iim,jjm, hori_id2, 1, 1, 1, -99, 32, avescatter2(1), dt, dw2)
             CALL histdef(hist2_id, 'SubSnow','Snow sublimation','kg/m^2/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt, dw2)
             CALL histdef(hist2_id, 'ACond', 'Aerodynamic conductance', 'm/s',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(5), dt, dw2)
             !-
             !-
             !-  Cold Season Processes
             !-
             CALL histdef(hist2_id, 'SnowFrac', 'Snow cover fraction', '1',  &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'SAlbedo', 'Snow albedo', '1', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             CALL histdef(hist2_id, 'SnowDepth', '3D snow depth', 'm', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, avescatter2(7), dt, dw2)
             !-
             !- Hydrologic variables
             !-
             CALL histdef(hist2_id, 'IrrigationMap', 'Map of irrigated areas', 'm^2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt, dw2)
             CALL histdef(hist2_id, 'FloodplainsMap', 'Map of flooded areas', 'm^2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt,dw2)
             CALL histdef(hist2_id, 'SwampMap', 'Map of swamp areas', 'm^2', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt,dw2)
             CALL histdef(hist2_id, 'Dis', 'Simulated River Discharge', 'm^3/s', &
                  & iim,jjm, hori_id2, 1,1,1, -99, 32, fluxop_scinsec2(1), dt,dw2)
             CALL histdef(hist2_id, 'humrel', 'Soil moisture stress', '-',  &
                  & iim,jjm, hori_id2, nvm,1,nvm, vegax_id2, 32, avescatter2(9), dt, dw2)
             !-
             !-  The carbon budget
             !-
             IF ( control_flags%ok_co2 ) THEN
                CALL histdef(hist2_id, 'GPP', 'Net assimilation of carbon by the vegetation', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt, dw2)
             ENDIF
             IF ( control_flags%ok_stomate ) THEN
                CALL histdef(hist2_id, 'NEE', 'Net Ecosystem Exchange', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt,dw2)
                CALL histdef(hist2_id, 'maint_resp', 'Maintenance respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt,dw2)
                CALL histdef(hist2_id, 'hetero_resp', 'Heterotrophic respiration', 'gC/m^2/s', &
                  & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt,dw2)
                CALL histdef(hist2_id, 'growth_resp', 'Growth respiration', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt, dw2)
                CALL histdef(hist2_id, 'npp', 'Net Primary Production', 'gC/m^2/s', &
                     & iim,jjm, hori_id2, nvm, 1, nvm, vegax_id2, 32, fluxop_scinsec2(1), dt, dw2)
             ENDIF
             !
          ENDIF
          !-
          CALL histdef(hist2_id, 'LandPoints', 'Land Points', '1', &
               & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt, dw2)  
          CALL histdef(hist2_id, 'Areas', 'Mesh areas', 'm2', &
               & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt, dw2)
          CALL histdef(hist2_id, 'Contfrac', 'Continental fraction', '1', &
               & iim,jjm, hori_id2, 1,1,1, -99, 32, once2(1), dt, dw2)
          !-
          ! Write the names of the pfts in the high frequency sechiba history files
          global_attribute="PFT_name"
          DO i=1,nvm
             WRITE(global_attribute(9:10),"(I2.2)") i
             CALL histglobal_attr(hist2_id, global_attribute, PFT_name(i))
          ENDDO
          !-
          CALL histend(hist2_id)
      ENDIF
  ENDIF

    !-
    !=====================================================================
    !- 3.2 STOMATE's history file
    !=====================================================================
    IF ( control_flags%ok_stomate ) THEN
       !-
       ! STOMATE IS ACTIVATED
       !-
       !Config Key   = STOMATE_OUTPUT_FILE
       !Config Desc  = Name of file in which STOMATE's output is going to be written
       !Config If    = OK_STOMATE
       !Config Def   = stomate_history.nc
       !Config Help  = This file is going to be created by the model
       !Config         and will contain the output from the model.
       !Config         This file is a truly COADS compliant netCDF file.
       !Config         It will be generated by the hist software from
       !Config         the IOIPSL package.
       !Config Units = [FILE]
       !-
       stom_histname='stomate_history.nc'
       CALL getin_p('STOMATE_OUTPUT_FILE', stom_histname)       
       WRITE(numout,*) 'STOMATE_OUTPUT_FILE', TRIM(stom_histname)
       !-
       !Config Key   = STOMATE_HIST_DT
       !Config Desc  = STOMATE history time step
       !Config If    = OK_STOMATE
       !Config Def   = 10.
       !Config Help  = Time step of the STOMATE history file
       !Config Units = [days]
       !-
       hist_days_stom = 10.
       CALL getin_p('STOMATE_HIST_DT', hist_days_stom)       

       IF ( hist_id < 0 ) THEN
          ! Deactivate all history files if sechiba_history file is deactivated
          hist_dt_stom=0
          WRITE(numout,*) 'STOMATE history file will not be created because sechiba_history file is deactivated.'
       ELSE IF ( hist_days_stom == moins_un ) THEN
          hist_dt_stom = moins_un
          WRITE(numout,*) 'output frequency for STOMATE history file (d): one month.'
       ELSE IF ( hist_days_stom == 0 ) THEN
          ! Deactivate this file
          hist_dt_stom=0
          WRITE(numout,*) 'STOMATE history file will not be created'
       ELSE
          hist_dt_stom = NINT( hist_days_stom ) * one_day
          WRITE(numout,*) 'output frequency for STOMATE history file (d): ', &
               hist_dt_stom/one_day
       ENDIF

       ! test consistency between STOMATE_HIST_DT and DT_SLOW parameters
       dt_slow_ = one_day
       CALL getin_p('DT_SLOW', dt_slow_)
       IF ( hist_days_stom /= moins_un .AND. hist_dt_stom/=0) THEN
          IF (dt_slow_ > hist_dt_stom) THEN
             WRITE(numout,*) "DT_SLOW = ",dt_slow_,"  , STOMATE_HIST_DT = ",hist_dt_stom
             CALL ipslerr_p (3,'intsurf_history', &
                  &          'Problem with DT_SLOW > STOMATE_HIST_DT','', &
                  &          '(must be less or equal)')
          ENDIF
       ENDIF
       !-
       !- Initialize stomate_history file
       IF ( hist_dt_stom == 0 ) THEN
          ! Case hist_dt_stom=0 : No creation of stomate_history.nc file
          ! Nothing will be done.
          hist_id_stom=-1
       ELSE
          ! Initialise stomate_history file
       IF (is_omp_root) THEN
          IF ( rectilinear ) THEN
#ifdef CPP_PARA
             CALL histbeg(stom_histname, iim, lon_rect, jjm, lat_rect,  1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id_stom,domain_id=orch_domain_id)
#else
             CALL histbeg(stom_histname, iim, lon_rect, jjm, lat_rect,  1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id_stom)
#endif
          ELSE
#ifdef CPP_PARA
             CALL histbeg(stom_histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id_stom,domain_id=orch_domain_id)
#else
             CALL histbeg(stom_histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
                  &     istp_old, date0, dt, hori_id, hist_id_stom)
#endif
          ENDIF
          !- define PFT axis
          hist_PFTaxis = (/ ( REAL(i,r_std), i=1,nvm ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'PFT', 'Plant functional type', &
               & '1', nvm, hist_PFTaxis, hist_PFTaxis_id)
          ! deforestation
          !- define Pool_10 axis
          hist_pool_10axis = (/ ( REAL(i,r_std), i=1,10 ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'P10', 'Pool 10 years', &
               & '1', 10, hist_pool_10axis, hist_pool_10axis_id)
          
          !- define Pool_100 axis
          hist_pool_100axis = (/ ( REAL(i,r_std), i=1,100 ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'P100', 'Pool 100 years', &
               & '1', 100, hist_pool_100axis, hist_pool_100axis_id)
          
          !- define Pool_11 axis
          hist_pool_11axis = (/ ( REAL(i,r_std), i=1,11 ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'P11', 'Pool 10 years + 1', &
               & '1', 11, hist_pool_11axis, hist_pool_11axis_id)
          
          !- define Pool_101 axis
          hist_pool_101axis = (/ ( REAL(i,r_std), i=1,101 ) /)
          !- declare this axis
          CALL histvert (hist_id_stom, 'P101', 'Pool 100 years + 1', &
               & '1', 101, hist_pool_101axis, hist_pool_101axis_id)
          !- define deep permafrost axis for stomate variables
          CALL histvert(hist_id_stom, 'solth', 'deep soil levels',      'm', &
               &    ngrnd, sol, hist_stomate_deepsoil)

          snowlev = (/ ( REAL(i,r_std), i=1,nsnow ) /)
          CALL histvert(hist_id_stom, 'snowlev', 'snow levels',      'index', &
               &    nsnow, snowlev, hist_stomate_snow)

       ENDIF
       !- define STOMATE history file
       CALL stom_define_history (hist_id_stom, nvm, iim, jjm, &
            & dt, hist_dt_stom, hori_id, hist_PFTaxis_id, &
            & hist_pool_10axis_id, hist_pool_100axis_id, &
            & hist_pool_11axis_id, hist_pool_101axis_id, &
            & hist_stomate_deepsoil, hist_stomate_snow)
       
       !- Write the names of the pfts in the stomate history files 
       IF (is_omp_root) THEN
          global_attribute="PFT_name"
          DO i=1,nvm
             WRITE(global_attribute(9:10),"(I2.2)") i
             CALL histglobal_attr(hist_id_stom, global_attribute, PFT_name(i))
          ENDDO

       !- end definition
          CALL histend(hist_id_stom)
       ENDIF
    END IF ! IF ( hist_dt_stom == 0 )

       !-
       !-
       !-
       ! STOMATE IPCC OUTPUTS IS ACTIVATED
       !-
       !Config Key   = STOMATE_IPCC_OUTPUT_FILE
       !Config Desc  = Name of file in which STOMATE's output is going to be written
       !Config If    = OK_STOMATE
       !Config Def   = stomate_ipcc_history.nc
       !Config Help  = This file is going to be created by the model
       !Config         and will contain the output from the model.
       !Config         This file is a truly COADS compliant netCDF file.
       !Config         It will be generated by the hist software from
       !Config         the IOIPSL package.
       !Config Units = [FILE]
       !-
       stom_ipcc_histname='stomate_ipcc_history.nc'
       CALL getin_p('STOMATE_IPCC_OUTPUT_FILE', stom_ipcc_histname)       
       WRITE(numout,*) 'STOMATE_IPCC_OUTPUT_FILE', TRIM(stom_ipcc_histname)
       !-
       !Config Key   = STOMATE_IPCC_HIST_DT
       !Config Desc  = STOMATE IPCC history time step
       !Config If    = OK_STOMATE
       !Config Def   = 0.
       !Config Help  = Time step of the STOMATE IPCC history file
       !Config Units = [days]
       !-
       hist_days_stom_ipcc = zero
       CALL getin_p('STOMATE_IPCC_HIST_DT', hist_days_stom_ipcc)       
       IF ( hist_days_stom_ipcc == moins_un ) THEN
          hist_dt_stom_ipcc = moins_un
          WRITE(numout,*) 'output frequency for STOMATE IPCC history file (d): one month.'
       ELSE
          hist_dt_stom_ipcc = NINT( hist_days_stom_ipcc ) * one_day
          WRITE(numout,*) 'output frequency for STOMATE IPCC history file (d): ', &
            hist_dt_stom_ipcc/one_day
       ENDIF
       
       IF ( hist_dt_stom_ipcc /= 0 .AND. hist_id < 0 ) THEN
          ! sechiba_history file is not created therefore STOMATE IPCC history file will be deactivated
          hist_dt_stom_ipcc=0
          hist_days_stom_ipcc=0
          WRITE(numout,*) 'STOMATE IPCC history file is not created.'
       END IF

       ! test consistency between STOMATE_IPCC_HIST_DT and DT_SLOW parameters
       dt_slow_ = one_day
       CALL getin_p('DT_SLOW', dt_slow_)
       IF ( hist_days_stom_ipcc > zero ) THEN
          IF (dt_slow_ > hist_dt_stom_ipcc) THEN
             WRITE(numout,*) "DT_SLOW = ",dt_slow_,"  , STOMATE_IPCC_HIST_DT = ",hist_dt_stom_ipcc
             CALL ipslerr_p (3,'intsurf_history', &
                  &          'Problem with DT_SLOW > STOMATE_IPCC_HIST_DT','', &
                  &          '(must be less or equal)')
          ENDIF
       ENDIF

       !Config Key   = OK_HISTSYNC
       !Config Desc  = Syncronize and write at each time step the IOIPSL output files
       !Config If    = 
       !Config Def   = FALSE
       !Config Units = [FLAG]
       ok_histsync=.FALSE.
       CALL getin_p('OK_HISTSYNC', ok_histsync)       



       IF ( hist_dt_stom_ipcc == 0 ) THEN
          hist_id_stom_ipcc = -1
       ELSE
          !-
          !- initialize
          IF (is_omp_root) THEN
             IF ( rectilinear ) THEN
#ifdef CPP_PARA
                CALL histbeg(stom_ipcc_histname, iim, lon_rect, jjm, lat_rect,  1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id, hist_id_stom_ipcc,domain_id=orch_domain_id)
#else
                CALL histbeg(stom_ipcc_histname, iim, lon_rect, jjm, lat_rect,  1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id, hist_id_stom_ipcc)
#endif
             ELSE
#ifdef CPP_PARA
                CALL histbeg(stom_ipcc_histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id, hist_id_stom_ipcc,domain_id=orch_domain_id)
#else
                CALL histbeg(stom_ipcc_histname, iim, lon, jjm, lat,  1, iim, 1, jjm, &
                     &     istp_old, date0, dt, hori_id, hist_id_stom_ipcc)
#endif
             ENDIF
             !- declare this axis
             CALL histvert (hist_id_stom_IPCC, 'PFT', 'Plant functional type', &
                  & '1', nvm, hist_PFTaxis, hist_IPCC_PFTaxis_id)
             
             !- define STOMATE history file
             CALL stom_IPCC_define_history (hist_id_stom_IPCC, nvm, iim, jjm, &
                  & dt, hist_dt_stom_ipcc, hori_id, hist_IPCC_PFTaxis_id)
             
             !- Write the names of the pfts in the stomate history files 
             global_attribute="PFT_name"
             DO i=1,nvm
                WRITE(global_attribute(9:10),"(I2.2)") i
                CALL histglobal_attr(hist_id_stom_IPCC, global_attribute, PFT_name(i))
             ENDDO

             !- end definition
             CALL histend(hist_id_stom_IPCC)
          ENDIF
      ENDIF
   ENDIF


    RETURN

  END SUBROUTINE intsurf_history
  !  
  SUBROUTINE stom_define_history &
       & (hist_id_stom, nvm, iim, jjm, dt, &
       &  hist_dt, hist_hori_id, hist_PFTaxis_id, &
       & hist_pool_10axis_id, hist_pool_100axis_id, &
       & hist_pool_11axis_id, hist_pool_101axis_id, &
       & hist_stomate_deepsoil, hist_stomate_snow)
    ! deforestation axis added as arguments

    !---------------------------------------------------------------------
    !- Tell ioipsl which variables are to be written
    !- and on which grid they are defined
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !-
    !- Input
    !-
    !- File id
    INTEGER(i_std),INTENT(in) :: hist_id_stom
    !- number of PFTs
    INTEGER(i_std),INTENT(in) :: nvm
    !- Domain size
    INTEGER(i_std),INTENT(in) :: iim, jjm
    !- Time step of STOMATE (seconds)
    REAL(r_std),INTENT(in)    :: dt
    !- Time step of history file (s)
    REAL(r_std),INTENT(in)    :: hist_dt
    !- id horizontal grid
    INTEGER(i_std),INTENT(in) :: hist_hori_id
    !- id of PFT axis
    INTEGER(i_std),INTENT(in) :: hist_PFTaxis_id
    !- id of Deforestation axis
    INTEGER(i_std),INTENT(in) :: hist_pool_10axis_id,hist_pool_100axis_id
    INTEGER(i_std),INTENT(in) :: hist_pool_11axis_id,hist_pool_101axis_id
    !-  id of permafrost axis
    INTEGER(i_std),INTENT(in) :: hist_stomate_deepsoil
    INTEGER(i_std),INTENT(in)     :: hist_stomate_snow

    !- 1 local
    !-
    !- maximum history level
    INTEGER(i_std), PARAMETER  :: max_hist_level = 10
    !- output level (between 0 and 10)
    !-  ( 0:nothing is written, 10:everything is written)
    INTEGER(i_std)             :: hist_level
    !- Character strings to define operations for histdef
    CHARACTER(LEN=40),DIMENSION(max_hist_level) :: ave
    !- for looping over PFT dimension for permafrost soil variables
    INTEGER(i_std)     :: jv, m
    CHARACTER(LEN=10)  :: part_str    ! string suffix indicating an index

    !---------------------------------------------------------------------
    !=====================================================================
    !- 1 history level
    !=====================================================================
    !- 1.1 define history levelx
    !=====================================================================
    !Config Key   = STOMATE_HISTLEVEL
    !Config Desc  = STOMATE history output level (0..10)
    !Config If    = OK_STOMATE
    !Config Def   = 10
    !Config Help  = 0: nothing is written; 10: everything is written
    !Config Units = [-]
    !-
    hist_level = 10
    CALL getin_p('STOMATE_HISTLEVEL', hist_level)
    !-
    WRITE(numout,*) 'STOMATE history level: ',hist_level
    IF ( (hist_level > max_hist_level).OR.(hist_level < 0) ) THEN
       STOP 'This history level is not allowed'
    ENDIF
    !=====================================================================
    !- 1.2 define operations according to output level
    !=====================================================================
    ave(1:hist_level) =  'ave(scatter(X))'
    ave(hist_level+1:max_hist_level) =  'never          '
    !=====================================================================
    !- 2 surface fields (2d)
    !- 3 PFT: 3rd dimension
    !=====================================================================


    ! structural litter above ground
    IF (is_omp_root) THEN
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTER_STR_AB       "), &
            &               TRIM("structural litter above ground                    "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! metabolic litter above ground                     
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTER_MET_AB       "), &
            &               TRIM("metabolic litter above ground                     "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! structural litter below ground               
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTER_STR_BE       "), &
            &               TRIM("structural litter below ground                    "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! metabolic litter below ground                
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTER_MET_BE       "), &
            &               TRIM("metabolic litter below ground                     "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! fraction of soil covered by dead leaves           
       CALL histdef (hist_id_stom, &
            &               TRIM("DEADLEAF_COVER      "), &
            &               TRIM("fraction of soil covered by dead leaves           "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)
       
       ! total soil and litter carbon
       CALL histdef (hist_id_stom, &
            &               TRIM("TOTAL_SOIL_CARB     "), &
            &               TRIM("total soil and litter carbon                      "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! active soil carbon in ground                 
       CALL histdef (hist_id_stom, &
            &               TRIM("CARBON_ACTIVE       "), &
            &               TRIM("active soil carbon in ground                      "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! slow soil carbon in ground                   
       CALL histdef (hist_id_stom, &
            &               TRIM("CARBON_SLOW         "), &
            &               TRIM("slow soil carbon in ground                        "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! passive soil carbon in ground                
       CALL histdef (hist_id_stom, &
            &               TRIM("CARBON_PASSIVE      "), &
            &               TRIM("passive soil carbon in ground                     "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
      

       ! active soil carbon in ground                 
       CALL histdef (hist_id_stom, &
           &               TRIM("CARBON_ACTIVE_SURF  "), &
           &               TRIM("active soil carbon in ground over surface soils"), &
           &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id,&
           &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

      ! slow soil carbon in ground                   
      CALL histdef (hist_id_stom, &
           &               TRIM("CARBON_SLOW_SURF    "), &
           &               TRIM("slow soil carbon in ground over surface soils "), &
           &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
           &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

      ! passive soil carbon in ground                
      CALL histdef (hist_id_stom, &
           &               TRIM("CARBON_PASSIVE_SURF "), &
           &               TRIM("passive soil carbon in ground over surface soils"), &
           &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
           &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
 
       ! Long term 2 m temperature                           
       CALL histdef (hist_id_stom, &
            &               TRIM("T2M_LONGTERM        "), &
            &               TRIM("Longterm 2 m temperature                          "), &
            &               TRIM("K                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(9), dt, hist_dt)
       
       ! Monthly 2 m temperature                           
       CALL histdef (hist_id_stom, &
            &               TRIM("T2M_MONTH           "), &
            &               TRIM("Monthly 2 m temperature                           "), &
            &               TRIM("K                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(1), dt, hist_dt)
      
       ! "seasonal" 2 m temperature                           
       CALL histdef (hist_id_stom, &
         &               TRIM("TSEASON             "), &
         &               TRIM("Seasonal 2 m temperature                             "), &
         &               TRIM("K                   "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(10), dt, hist_dt)

       ! Tmin in spring (after onset)                           
       CALL histdef (hist_id_stom, &
         &               TRIM("TMIN_SPRING         "), &
         &               TRIM("Tmin in spring (after onset)                         "), &
         &               TRIM("K                   "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

       ! how many days after onset                           
       CALL histdef (hist_id_stom, &
         &               TRIM("TMIN_SPRING_TIME    "), &
         &               TRIM("how many days after onset                            "), &
         &               TRIM("days                "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

       !                           
       CALL histdef (hist_id_stom, &
         &               TRIM("ONSET_DATE          "), &
         &               TRIM("onset date                                           "), &
         &               TRIM("day                 "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

       ! minimum 2 m temperature                           
       CALL histdef (hist_id_stom, &
         &               TRIM("T2M_MIN_DAILY       "), &
         &               TRIM("minimum 2 m temperature                              "), &
         &               TRIM("K                   "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(10), dt, hist_dt) 
       ! Weekly 2 m temperature                            
       CALL histdef (hist_id_stom, &
            &               TRIM("T2M_WEEK            "), &
            &               TRIM("Weekly 2 m temperature                            "), &
            &               TRIM("K                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(1), dt, hist_dt)
       
       ! heterotr. resp. from ground                  
       CALL histdef (hist_id_stom, &
            &               TRIM("HET_RESP            "), &
            &               TRIM("heterotr. resp. from ground                       "), &
            &               TRIM("gC/m^2 tot/pft/day  "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)
       
       ! black carbon on average total ground              
       CALL histdef (hist_id_stom, &
            &               TRIM("BLACK_CARBON        "), &
            &               TRIM("black carbon on average total ground              "), &
            &               TRIM("gC/m^2 tot          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(10), dt, hist_dt)
       
       ! Fire fraction on ground
       CALL histdef (hist_id_stom, &
            &               TRIM("FIREFRAC            "), &
            &               TRIM("Fire fraction on ground                           "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

       ! Fire index on ground                      
       CALL histdef (hist_id_stom, &
            &               TRIM("FIREINDEX           "), &
            &               TRIM("Fire index on ground                              "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
       
       ! Litter humidity                                   
       CALL histdef (hist_id_stom, &
            &               TRIM("LITTERHUM           "), &
            &               TRIM("Litter humidity                                   "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)
       
       ! CO2 flux                                  
       CALL histdef (hist_id_stom, &
            &               TRIM("CO2FLUX             "), &
            &               TRIM("CO2 flux                                          "), &
            &               TRIM("gC/m^2/pft/mth      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

!!$    CALL histdef(hist_id_stom, &
!!$         &               TRIM("CO2FLUX_MONTHLY_SUM "), &
!!$         &               TRIM("Monthly CO2 flux Sum                              "), &
!!$         &               TRIM("PgC/m^2/mth         "), iim,jjm, hist_hori_id, &
!!$         &               1,1,1, -99, 32, 'inst(scatter(X))', dt, hist_dt)

       ! Output CO2 flux from fire                         
       CALL histdef (hist_id_stom, &
            &               TRIM("CO2_FIRE            "), &
            &               TRIM("Output CO2 flux from fire                         "), &
            &               TRIM("gC/day/m^2/pft      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)
       
       ! CO2 taken from atmosphere for initiate growth     
       CALL histdef (hist_id_stom, &
            &               TRIM("CO2_TAKEN           "), &
            &               TRIM("CO2 taken from atmosphere for initiate growth     "), &
            &               TRIM("gC/day/m^2/pft      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)
       
       ! Leaf Area Index                                   
       CALL histdef (hist_id_stom, &
            &               TRIM("LAI                 "), &
            &               TRIM("Leaf Area Index                                   "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("FPC_MAX             "), &
            &               TRIM("foliage projective cover                          "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("MAXFPC_LASTYEAR     "), &
            &               TRIM("foliage projective cover of last year             "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
       
       ! Maximum vegetation fraction (LAI -> infinity)     
       CALL histdef (hist_id_stom, &
            &               TRIM("VEGET_MAX           "), &
            &               TRIM("Maximum vegetation fraction (LAI -> infinity)     "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)
       
       ! Net primary productivity                          
       CALL histdef (hist_id_stom, &
            &               TRIM("NPP                 "), &
            &               TRIM("Net primary productivity                          "), &
            &               TRIM("gC/day/(m^2 tot)    "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

       ! Gross primary productivity                        
       CALL histdef (hist_id_stom, &
            &               TRIM("GPP                 "), &
            &               TRIM("Gross primary productivity                        "), &
            &               TRIM("gC/day/m^2          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)

       ! Density of individuals                            
       CALL histdef (hist_id_stom, &
            &               TRIM("IND                 "), &
            &               TRIM("Density of individuals                            "), &
            &               TRIM("1/ m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       ! Adaptation to climate
       CALL histdef (hist_id_stom, &
            &               TRIM("ADAPTATION          "), &
            &               TRIM("Adaptation to climate (DGVM)                      "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)
    
       ! Probability from regenerative
       CALL histdef (hist_id_stom, &
            &               TRIM("REGENERATION        "), &
            &               TRIM("Probability from regenerative (DGVM)               "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)
       
       ! crown area of individuals (m**2)
       CALL histdef (hist_id_stom, &
            &               TRIM("CN_IND              "), &
            &               TRIM("crown area of individuals                         "), &
            &               TRIM("m^2                 "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       ! woodmass of individuals (gC)
       CALL histdef (hist_id_stom, &
            &               TRIM("WOODMASS_IND        "), &
            &               TRIM("Woodmass of individuals                           "), &
            &               TRIM("gC/pft              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       ! total living biomass
       CALL histdef (hist_id_stom, &
            &               TRIM("TOTAL_M             "), &
            &               TRIM("Total living biomass                              "), &
            &               TRIM("gC/m^2/pft          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)
       
       ! Leaf mass                                         
       CALL histdef (hist_id_stom, &
            &               TRIM("LEAF_M              "), &
            &               TRIM("Leaf mass                                         "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! Sap mass above ground                             
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_M_AB            "), &
            &               TRIM("Sap mass above ground                             "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)

       ! Sap mass below ground                             
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_M_BE            "), &
            &               TRIM("Sap mass below ground                             "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! Heartwood mass above ground                       
       CALL histdef (hist_id_stom, &
            &               TRIM("HEART_M_AB          "), &
            &               TRIM("Heartwood mass above ground                       "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)

       ! Heartwood mass below ground                       
       CALL histdef (hist_id_stom, &
            &               TRIM("HEART_M_BE          "), &
            &               TRIM("Heartwood mass below ground                       "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)

       ! Root mass                                         
       CALL histdef (hist_id_stom, &
            &               TRIM("ROOT_M              "), &
            &               TRIM("Root mass                                         "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! Fruit mass                                        
       CALL histdef (hist_id_stom, &
            &               TRIM("FRUIT_M             "), &
            &               TRIM("Fruit mass                                        "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! Carbohydrate reserve mass                         
       CALL histdef (hist_id_stom, &
            &               TRIM("RESERVE_M           "), &
            &               TRIM("Carbohydrate reserve mass                         "), &
            &               TRIM("gC/m^2              "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! total turnover rate
       CALL histdef (hist_id_stom, &
            &               TRIM("TOTAL_TURN          "), &
            &               TRIM("total turnover rate                               "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Leaf turnover                                     
       CALL histdef (hist_id_stom, &
            &               TRIM("LEAF_TURN           "), &
            &               TRIM("Leaf turnover                                     "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Sap turnover above                                
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_AB_TURN         "), &
            &               TRIM("Sap turnover above                                "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Root turnover                                     
       CALL histdef (hist_id_stom, &
            &               TRIM("ROOT_TURN           "), &
            &               TRIM("Root turnover                                     "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Fruit turnover                                    
       CALL histdef (hist_id_stom, &
            &               TRIM("FRUIT_TURN          "), &
            &               TRIM("Fruit turnover                                    "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! total conversion of biomass to litter 
       CALL histdef (hist_id_stom, &
            &               TRIM("TOTAL_BM_LITTER     "), &
            &               TRIM("total conversion of biomass to litter             "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Leaf death                                        
       CALL histdef (hist_id_stom, &
            &               TRIM("LEAF_BM_LITTER      "), &
            &               TRIM("Leaf death                                        "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)
       
       ! Sap death above ground                            
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_AB_BM_LITTER    "), &
            &               TRIM("Sap death above ground                            "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Sap death below ground                            
       CALL histdef (hist_id_stom, &
            &               TRIM("SAP_BE_BM_LITTER    "), &
            &               TRIM("Sap death below ground                            "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Heartwood death above ground                      
       CALL histdef (hist_id_stom, &
            &               TRIM("HEART_AB_BM_LITTER  "), &
            &               TRIM("Heartwood death above ground                      "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Heartwood death below ground                      
       CALL histdef (hist_id_stom, &
            &               TRIM("HEART_BE_BM_LITTER  "), &
            &               TRIM("Heartwood death below ground                      "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Root death                                        
       CALL histdef (hist_id_stom, &
            &               TRIM("ROOT_BM_LITTER      "), &
            &               TRIM("Root death                                        "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)
       
       ! Fruit death                                       
       CALL histdef (hist_id_stom, &
            &               TRIM("FRUIT_BM_LITTER     "), &
            &               TRIM("Fruit death                                       "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Carbohydrate reserve death                        
       CALL histdef (hist_id_stom, &
            &               TRIM("RESERVE_BM_LITTER   "), &
            &               TRIM("Carbohydrate reserve death                        "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(4), dt, hist_dt)

       ! Maintenance respiration                           
       CALL histdef (hist_id_stom, &
            &               TRIM("MAINT_RESP          "), &
            &               TRIM("Maintenance respiration                           "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)

       ! Growth respiration                                
       CALL histdef (hist_id_stom, &
            &               TRIM("GROWTH_RESP         "), &
            &               TRIM("Growth respiration                                "), &
            &               TRIM("gC/m^2/day          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(2), dt, hist_dt)
       
       ! age                                               
       CALL histdef (hist_id_stom, &
            &               TRIM("AGE                 "), &
            &               TRIM("age                                               "), &
            &               TRIM("years               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(7), dt, hist_dt)
       
       ! height                                            
       CALL histdef (hist_id_stom, &
            &               TRIM("HEIGHT              "), &
            &               TRIM("height                                            "), &
            &               TRIM("m                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(7), dt, hist_dt)

       ! weekly moisture stress                            
       CALL histdef (hist_id_stom, &
            &               TRIM("MOISTRESS           "), &
            &               TRIM("weekly moisture stress                            "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       ! Maximum rate of carboxylation                     
       CALL histdef (hist_id_stom, &
            &               TRIM("VCMAX               "), &
            &               TRIM("Maximum rate of carboxylation                     "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! leaf age                                          
       CALL histdef (hist_id_stom, &
            &               TRIM("LEAF_AGE            "), &
            &               TRIM("leaf age                                          "), &
            &               TRIM("days                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)
       
       ! Fraction of trees that dies (gap)                 
       CALL histdef (hist_id_stom, &
            &               TRIM("MORTALITY           "), &
            &               TRIM("Fraction of trees that dies (gap)                 "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! Fraction of plants killed by fire                 
       CALL histdef (hist_id_stom, &
            &               TRIM("FIREDEATH           "), &
            &               TRIM("Fraction of plants killed by fire                 "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! Density of newly established saplings             
       CALL histdef (hist_id_stom, &
            &               TRIM("IND_ESTAB           "), &
            &               TRIM("Density of newly established saplings             "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! Establish tree
       CALL histdef (hist_id_stom, &
            &               TRIM("ESTABTREE           "), &
            &               TRIM("Rate of tree establishement                       "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(10), dt, hist_dt)

       ! Establish grass
       CALL histdef (hist_id_stom, &
            &               TRIM("ESTABGRASS          "), &
            &               TRIM("Rate of grass establishement                      "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(10), dt, hist_dt)

!!DZADD
!       ! leaf_frac[:,:,1] after lpj_establish
!       CALL histdef (hist_id_stom, &
!            &               TRIM("EST_LEAF_FRAC1      "), &
!            &               TRIM("leaf_frac_1 after lpj_establishs                  "), &
!            &               TRIM(" /                  "), iim,jjm, hist_hori_id, &
!            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!       ! leaf_frac[:,:,2] after lpj_establish
!       CALL histdef (hist_id_stom, &
!            &               TRIM("EST_LEAF_FRAC2      "), &
!            &               TRIM("leaf_frac_2 after lpj_establishs                  "), &
!            &               TRIM(" /                  "), iim,jjm, hist_hori_id, &
!            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!       ! leaf_frac[:,:,3] after lpj_establish
!       CALL histdef (hist_id_stom, &
!            &               TRIM("EST_LEAF_FRAC3      "), &
!            &               TRIM("leaf_frac_3 after lpj_establishs                  "), &
!            &               TRIM(" /                  "), iim,jjm, hist_hori_id, &
!            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!       ! leaf_frac[:,:,4] after lpj_establish
!       CALL histdef (hist_id_stom, &
!            &               TRIM("EST_LEAF_FRAC4      "), &
!            &               TRIM("leaf_frac_4 after lpj_establishs                  "), &
!            &               TRIM(" /                  "), iim,jjm, hist_hori_id, &
!            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!       CALL histdef (hist_id_stom, &
!            &               TRIM("LEAF_FRAC1      "), &
!            &               TRIM("leaf_frac_1 after stomate_lpj                     "), &
!            &               TRIM(" /                  "), iim,jjm, hist_hori_id, &
!            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!       CALL histdef (hist_id_stom, &
!            &               TRIM("LEAF_FRAC2      "), &
!            &               TRIM("leaf_frac_2 after stomate_lpj                     "), &
!            &               TRIM(" /                  "), iim,jjm, hist_hori_id, &
!            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!       CALL histdef (hist_id_stom, &
!            &               TRIM("LEAF_FRAC3      "), &
!            &               TRIM("leaf_frac_3 after stomate_lpj                     "), &
!            &               TRIM(" /                  "), iim,jjm, hist_hori_id, &
!            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!       CALL histdef (hist_id_stom, &
!            &               TRIM("LEAF_FRAC4      "), &
!            &               TRIM("leaf_frac_4 after stomate_lpj                     "), &
!            &               TRIM(" /                  "), iim,jjm, hist_hori_id, &
!            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!!ENDDZADD

       ! Fraction of plants that dies (light competition)  
       CALL histdef (hist_id_stom, &
            &               TRIM("LIGHT_DEATH         "), &
            &               TRIM("Fraction of plants that dies (light competition)  "), &
            &               TRIM("1/day               "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(6), dt, hist_dt)

       ! biomass allocated to leaves                       
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_LEAF       "), &
            &               TRIM("biomass allocated to leaves                       "), &
            &               TRIM("gC/m**2/pft/dt      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to sapwood above ground         
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_SAP_AB     "), &
            &               TRIM("biomass allocated to sapwood above ground         "), &
            &               TRIM("gC/m**2/pft/dt      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to sapwood below ground         
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_SAP_BE     "), &
            &               TRIM("biomass allocated to sapwood below ground         "), &
            &               TRIM("gC/m**2/pft/dt      "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to roots                        
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_ROOT       "), &
            &               TRIM("biomass allocated to roots                        "), &
            &               TRIM("gC/m**2/pft/dt          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to fruits                       
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_FRUIT      "), &
            &               TRIM("biomass allocated to fruits                       "), &
            &               TRIM("gC/m**2/pft/dt          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! biomass allocated to carbohydrate reserve         
       CALL histdef (hist_id_stom, &
            &               TRIM("BM_ALLOC_RES        "), &
            &               TRIM("biomass allocated to carbohydrate reserve         "), &
            &               TRIM("gC/m**2/pft/dt          "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       ! time constant of herbivore activity               
       CALL histdef (hist_id_stom, &
            &               TRIM("HERBIVORES          "), &
            &               TRIM("time constant of herbivore activity               "), &
            &               TRIM("days                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

       CALL histdef (hist_id_stom, &
         &               TRIM("SENESCENCE          "), &
         &               TRIM("Signal to senescence                                 "), &
         &               TRIM("-                   "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

       ! turnover time for grass leaves                    
       CALL histdef (hist_id_stom, &
            &               TRIM("TURNOVER_TIME       "), &
            &               TRIM("turnover time for grass leaves                    "), &
            &               TRIM("days                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       
       ! 10 year wood product pool                         
       CALL histdef (hist_id_stom, &
            &               TRIM("PROD10              "), &
            &               TRIM("10 year wood product pool                         "), &
            &               TRIM("gC/m**2             "), iim,jjm, hist_hori_id, &
            &               11,1,11, hist_pool_11axis_id,32, ave(5), dt, hist_dt)
       
       ! annual flux for each 10 year wood product pool    
       CALL histdef (hist_id_stom, &
            &               TRIM("FLUX10              "), &
            &               TRIM("annual flux for each 10 year wood product pool    "), &
            &               TRIM("gC/m**2/yr          "), iim,jjm, hist_hori_id, &
            &               10,1,10, hist_pool_10axis_id,32, ave(5), dt, hist_dt)
       
       ! 100 year wood product pool                        
       CALL histdef (hist_id_stom, &
            &               TRIM("PROD100             "), &
            &               TRIM("100 year wood product pool                        "), &
            &               TRIM("gC/m**2             "), iim,jjm, hist_hori_id, &
            &               101,1,101, hist_pool_101axis_id,32, ave(5), dt, hist_dt)

       ! annual flux for each 100 year wood product pool   
       CALL histdef (hist_id_stom, &
            &               TRIM("FLUX100             "), &
            &               TRIM("annual flux for each 100 year wood product pool   "), &
            &               TRIM("gC/m**2/yr          "), iim,jjm, hist_hori_id, &
            &               100,1,100, hist_pool_100axis_id,32, ave(5), dt, hist_dt)

       ! annual release right after deforestation          
       CALL histdef (hist_id_stom, &
            &               TRIM("CONVFLUX            "), &
            &               TRIM("annual release right after deforestation          "), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       ! annual release from all 10 year wood product pools 
       CALL histdef (hist_id_stom, &
            &               TRIM("CFLUX_PROD10        "), &
            &               TRIM("annual release from all 10 year wood product pools"), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

       ! annual release from all 100year wood product pools
       CALL histdef (hist_id_stom, &
            &               TRIM("CFLUX_PROD100       "), &
            &               TRIM("annual release from all 100year wood product pools"), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)
       ! agriculure product
       CALL histdef (hist_id_stom, &
            &               TRIM("HARVEST_ABOVE       "), &
            &               TRIM("annual release product after harvest              "), &
            &               TRIM("gC/m**2/day          "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)


       CALL histdef(hist_id_stom, 'RESOLUTION_X', 'E-W resolution', 'm', &
            & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
       CALL histdef(hist_id_stom, 'RESOLUTION_Y', 'N-S resolution', 'm', &
            & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
       CALL histdef(hist_id_stom, 'CONTFRAC', 'Continental fraction', '1', &
            & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
       CALL histdef(hist_id_stom, 'Areas', 'Mesh areas', 'm2', &
            & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
       
       !  Special outputs for phenology
       CALL histdef (hist_id_stom, &
            &               TRIM("WHEN_GROWTHINIT     "), &
            &               TRIM("Time elapsed from season beginning                "), &
            &               TRIM("d                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("PFTPRESENT          "), &
            &               TRIM("PFT exists                                        "), &
            &               TRIM("d                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("GDD_MIDWINTER       "), &
            &               TRIM("Growing degree days, since midwinter              "), &
            &               TRIM("degK                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &               TRIM("GDD_M5_DORMANCE     "), &
            &               TRIM("Growing degree days, since dormance               "), &
            &               TRIM("degK                "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("NCD_DORMANCE        "), &
            &               TRIM("Number of chilling days, since leaves were lost   "), &
            &               TRIM("d                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("ALLOW_INITPHENO     "), &
            &               TRIM("Allow to declare beginning of the growing season  "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)
       
       CALL histdef (hist_id_stom, &
            &               TRIM("BEGIN_LEAVES        "), &
            &               TRIM("Signal to start putting leaves on                 "), &
            &               TRIM("-                   "), iim,jjm, hist_hori_id, &
            &               nvm,1,nvm, hist_PFTaxis_id,32, ave(9), dt, hist_dt)

!JCADD
!GM0
    CALL histdef (hist_id_stom, &
         &               TRIM("GRAZINGC "), &
         &               TRIM("Grazing C "), &
         &               TRIM("kgC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
!GM1
    CALL histdef (hist_id_stom, &
         &               TRIM("GRAZINGCSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NANIMALTOT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("INTAKE_ANIMAL "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("INTAKE "), &
         &               TRIM("grazing animal intake "), &
         &               TRIM("kgDM/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("INTAKESUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("TRAMPLING "), &
         &               TRIM("litter from trample by animals "), &
         &               TRIM("kgC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MILK "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MILKSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MILKCSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MILKC "), &
         &               TRIM("C export by milk production during animal grazing "), &
         &               TRIM("kgC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
!GM11
    CALL histdef (hist_id_stom, &
         &               TRIM("MILKN "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MILKANIMAL "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("METHANE "), &
         &               TRIM("Methane emission by grazing animal "), &
         &               TRIM("kgC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("METHANE_ANI "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("RANIMALSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("METHANESUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("RANIMAL "), &
         &               TRIM("C loss through grazing animal respiration "), &
         &               TRIM("kgC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FAECESNSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FAECESCSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("URINECSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM21
    CALL histdef (hist_id_stom, &
         &               TRIM("URINENSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NEL "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("URINEN "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("URINEC "), &
         &               TRIM("C in urine "), &
         &               TRIM("kgC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FAECESC "), &
         &               TRIM("C in faeces "), &
         &               TRIM("kgC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FAECESN "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GRAZED_FRAC "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NB_ANI "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("IMPORT_YIELD "), &
         &               TRIM("potential harvest yield of last year "), &
         &               TRIM("kgDM/m^2/yr "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("EXTRA_FEED "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM31
    CALL histdef (hist_id_stom, &
         &               TRIM("COMPT_UGB "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NB_GRAZINGDAYS "), &
         &               TRIM("number of grazing days of last year "), &
         &               TRIM("days "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("AMOUNT_YIELD "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("CONSUMP "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("OUTSIDE_FOOD "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("ADD_NB_ANI "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
! repeated
!    CALL histdef (hist_id_stom, &
!         &               TRIM("GRAZINGCSUM "), &
!         &               TRIM("- "), &
!         &               TRIM("- "), iim,jjm, hist_hori_id, &
!         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!    CALL histdef (hist_id_stom, &
!         &               TRIM("NANIMALTOT "), &
!         &               TRIM("- "), &
!         &               TRIM("- "), iim,jjm, hist_hori_id, &
!         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!    CALL histdef (hist_id_stom, &
!         &               TRIM("INTAKE_ANIMAL "), &
!         &               TRIM("- "), &
!         &               TRIM("- "), iim,jjm, hist_hori_id, &
!         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!    CALL histdef (hist_id_stom, &
!         &               TRIM("INTAKE "), &
!         &               TRIM("- "), &
!         &               TRIM("- "), iim,jjm, hist_hori_id, &
!         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!    CALL histdef (hist_id_stom, &
!         &               TRIM("INTAKESUM "), &
!         &               TRIM("- "), &
!         &               TRIM("- "), iim,jjm, hist_hori_id, &
!         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!
!    CALL histdef (hist_id_stom, &
!         &               TRIM("TRAMPLING "), &
!         &               TRIM("- "), &
!         &               TRIM("- "), iim,jjm, hist_hori_id, &
!         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("BCSyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("BCSmature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("Weightyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("Weightmature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM41
    CALL histdef (hist_id_stom, &
         &               TRIM("Weightcalf "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MPyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MPmature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MPwyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MPwmature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MPposyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("MPposmature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NEByoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NEBmature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NEIyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM51
    CALL histdef (hist_id_stom, &
         &               TRIM("NEImature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DMIcyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DMIcmature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DMIfyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DMIfmature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DMIyoung "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DMImature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DMIcalf "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("OMD "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("Weightcows "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM61
    CALL histdef (hist_id_stom, &
         &               TRIM("BCScows "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("CH4young "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("CH4mature "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("TSOILCUMM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("YIELD_RETURN "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("REGCOUNT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FERTCOUNT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN1 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN2 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN3 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM71
    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN4 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN5 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN6 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN7 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN8 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN9 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GMEAN0 "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("WSH "), &
         &               TRIM("shoot structure mass "), &
         &               TRIM("kgDM/m^2 "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("WSHTOT "), &
         &               TRIM("total shoot structure mass "), &
         &               TRIM("kgDM/m^2 "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("WR "), &
         &               TRIM("root structure mass "), &
         &               TRIM("kgDM/m^2 "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
!GM81
    CALL histdef (hist_id_stom, &
         &               TRIM("WRTOT "), &
         &               TRIM("total root structure mass "), &
         &               TRIM("kgDM/m^2 "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("WSHTOTSUM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("SR_UGB "), &
         &               TRIM("instantaneous stocking rate "), &
         &               TRIM("HeadorLSU/m^2 "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FCORGANICFERTMETABOLIC "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FCORGANICFERTSTRUCT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FNORGANICFERTURINE "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FNORGANICFERTSTRUCT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FNORGANICFERTMETABOLIC "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NFERTNITTOT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NFERTAMMTOT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM91
    CALL histdef (hist_id_stom, &
         &               TRIM("LOSS "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("LOSSC "), &
         &               TRIM("Carbon loss as litter during cutting "), &
         &               TRIM("kg C/m**2 "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("LOSSN "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DM_CUTYEARLY "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("C_CUTYEARLY "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NFERT_TOTAL "), &
         &               TRIM("Total Nitrogen input "), &
         &               TRIM("kg N/ha "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("Ndep "), &
         &               TRIM("Nitrogen deposition from input "), &
         &               TRIM("kg N/ha "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("LEGUME_FRACTION "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("SOIL_FERTILITY "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("C "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM101
    CALL histdef (hist_id_stom, &
         &               TRIM("N "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("FN "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NTOT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NAPO "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NSYM "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("DEVSTAGE "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("TGROWTH "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GRAZINGCSTRUCT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GRAZINGNSTRUCT "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("GRAZINGWN "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)
!GM111
    CALL histdef (hist_id_stom, &
         &               TRIM("GRAZINGWC "), &
         &               TRIM("- "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(10), dt, hist_dt)

       ! 14days 2 m temperature
       CALL histdef (hist_id_stom, &
            &               TRIM("T2M_14            "), &
            &               TRIM("14days 2 m temperature"), &
            &               TRIM("K                   "), iim,jjm, hist_hori_id, &
            &               1,1,1, -99,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("LITTER_RESP "), &
         &               TRIM("heterotr. resp. from litter pool "), &
         &               TRIM("gC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("ACTIVE_RESP "), &
         &               TRIM("heterotr. resp. from active carbon pool "), &
         &               TRIM("gC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("SLOW_RESP "), &
         &               TRIM("heterotr. resp. from slow carbon pool "), &
         &               TRIM("gC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("PASSIVE_RESP "), &
         &               TRIM("heterotr. resp. from passive carbon pool "), &
         &               TRIM("gC/m^2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("N_LIMFERT "), &
         &               TRIM("Nitrogen limitation factor on vcmax "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("SLA_CALC "), &
         &               TRIM("sla calculated by leaf age "), &
         &               TRIM("m**2/gC "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NPP_ABOVE "), &
         &               TRIM("Net above primary productivity "), &
         &               TRIM("gC/day/(m^2 (n/a)) "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)

    CALL histdef (hist_id_stom, &
         &               TRIM("NPP_BELOW "), &
         &               TRIM("Net below primary productivity "), &
         &               TRIM("gC/day/(m^2 (n/a)) "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
!GMtotal120
    CALL histdef (hist_id_stom, &
         &               TRIM("LITTER_STR_AVAIL "), &
         &               TRIM("Structural litter available for grazing "), &
         &               TRIM("gC/day/(m^2 (n/a)) "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
    CALL histdef (hist_id_stom, &
         &               TRIM("LITTER_MET_AVAIL "), &
         &               TRIM("Metabolic litter available for grazing "), &
         &               TRIM("gC/day/(m^2 (n/a)) "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
    CALL histdef (hist_id_stom, &
         &               TRIM("LITTER_STR_NAVAIL "), &
         &               TRIM("Structural litter not available for grazing "), &
         &               TRIM("gC/day/(m^2 (n/a)) "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
    CALL histdef (hist_id_stom, &
         &               TRIM("LITTER_MET_NAVAIL "), &
         &               TRIM("Metabolic litter not available for grazing "), &
         &               TRIM("gC/day/(m^2 (n/a)) "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
    CALL histdef (hist_id_stom, &
         &               TRIM("LITTER_STR_AVAILF "), &
         &               TRIM("Structural litter available fraction for grazing "), &
         &               TRIM("% "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
    CALL histdef (hist_id_stom, &
         &               TRIM("LITTER_MET_AVAILF "), &
         &               TRIM("Metabolic litter available fraction for grazing "), &
         &               TRIM("% "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
    CALL histdef (hist_id_stom, &
         &               TRIM("INTAKE_ANIMAL_LITTER "), &
         &               TRIM("Litter intake per animal "), &
         &               TRIM("kg DM/animal/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
    CALL histdef (hist_id_stom, &
         &               TRIM("INTAKE_LITTER "), &
         &               TRIM("Litter intake per m**2 "), &
         &               TRIM("kg DM/m**2/day "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
    CALL histdef (hist_id_stom, &
         &               TRIM("GRAZING_LITTER "), &
         &               TRIM("Flag of grazing litter 0 AGB 1 Litter 2 none "), &
         &               TRIM("- "), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)



!ENDJCADD
!
!variables for CH4 flux density from wetlands 
!
!pss:+
    CALL histdef (hist_id_stom, &
         &               TRIM("CH4_FLUX_TOT_0      "), &
         &               TRIM("flux density tot of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_dt)
    
    CALL histdef (hist_id_stom, &
         &               TRIM("CH4_FLUX_DIF_0      "), &
         &               TRIM("flux density dif of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_dt)
    
    CALL histdef (hist_id_stom, &
         &               TRIM("CH4_FLUX_BUB_0      "), &
         &               TRIM("flux density bub of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_dt)
    
    CALL histdef (hist_id_stom, &
         &               TRIM("CH4_FLUX_PLA_0      "), &
         &               TRIM("flux density pla of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_dt)
    
    !!pour wetland avec WTD = -x1
    CALL histdef (hist_id_stom, &
         &               TRIM("CH4_FLUX_TOT_wet1    "), &
         &               TRIM("flux density tot of CH4 by wetlands               "), &
         &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(5), dt, hist_dt)
    
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_DIF_wet1    "), &
        &               TRIM("flux density dif of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_BUB_wet1    "), &
        &               TRIM("flux density bub of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_PLA_wet1    "), &
        &               TRIM("flux density pla of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   !pour wetland avc WTD = -x2
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_TOT_wet2    "), &
        &               TRIM("flux density tot of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_DIF_wet2    "), &
        &               TRIM("flux density dif of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_BUB_wet2    "), &
        &               TRIM("flux density bub of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_PLA_wet2    "), &
        &               TRIM("flux density pla of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)

   !pour wetland avec WTD = -x3
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_TOT_wet3    "), &
        &               TRIM("flux density tot of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_DIF_wet3    "), &
        &               TRIM("flux density dif of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_BUB_wet3    "), &
        &               TRIM("flux density bub of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_PLA_wet3    "), &
        &               TRIM("flux density pla of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   !wetland avc WTD = -x4
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_TOT_wet4    "), &
        &               TRIM("flux density tot of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_DIF_wet4    "), &
        &               TRIM("flux density dif of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_BUB_wet4    "), &
        &               TRIM("flux density bub of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)
   
   CALL histdef (hist_id_stom, &
        &               TRIM("CH4_FLUX_PLA_wet4    "), &
        &               TRIM("flux density pla of CH4 by wetlands               "), &
        &               TRIM("mgCH4/d/m**2        "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(5), dt, hist_dt)

   !tsurf_year
   CALL histdef (hist_id_stom, &
        &               TRIM("TSURF_YEAR    "), &
        &               TRIM("Annual surface temperature                      "), &
        &               TRIM("K              "), iim,jjm, hist_hori_id, &
        &               1,1,1, -99,32, ave(1), dt, hist_dt)
   !pss:-

       ! permafrost variables
       ! first read logic on which variables to write to hist file.  (variables
       ! are
       ! stored in constantes_soil.f90)

       CALL getin_p ('writehist_deepC',writehist_deepC)
       CALL getin_p ('writehist_soilgases',writehist_soilgases)
       CALL getin_p ('writehist_deltaC',writehist_deltaC)
       CALL getin_p ('writehist_zimovheat',writehist_zimovheat)
       CALL getin_p ('writehist_deltaC_litter',writehist_deltaC_litter)
       CALL getin_p ('writehist_gascoeff',writehist_gascoeff)

       ! heterotr. resp. from ground                  
       CALL histdef (hist_id_stom, &
            &    TRIM("resp_hetero_litter   "), &
            &    TRIM("heterotr. resp. from litter                      "), &
            &    TRIM("gC/m^2 tot/day      "), iim,jjm, hist_hori_id, &
            &    nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &    TRIM("resp_hetero_soil     "), &
            &    TRIM("heterotr. resp. from standard stomate soil       "), &
            &    TRIM("gC/m^2 tot/day      "), iim,jjm, hist_hori_id, &
            &    nvm,1,nvm, hist_PFTaxis_id,32, ave(3), dt, hist_dt)
!++cdk: end of variables with implicit PFT dimension
       IF (writehist_deepC) THEN
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deepC_a_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("active pool deep soil (permafrost) carbon,PFT:"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("gC/m**3   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deepC_s_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("slow pool deep soil (permafrost) carbon   "), &
                     & TRIM("gC/m**3   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deepC_p_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("passive pool deep soil (permafrost) carbon   "), &
                     & TRIM("gC/m**3   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
       ENDIF
       IF (writehist_soilgases) THEN
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("O2_soil_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("deep soil (permafrost) oxygen   "), &
                     & TRIM("gC/m**3   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("CH4_soil_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("deep soil (permafrost) methane   "), &
                     & TRIM("gC/m**3   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("O2_snow_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("snow oxygen   "), &
                     & TRIM("gC/m**3   "), iim,jjm, hist_hori_id, &
                     & nsnow, 1, nsnow, hist_stomate_snow,32, ave(5), dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("CH4_snow_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("snow methane   "), &
                     & TRIM("gC/m**3   "), iim,jjm, hist_hori_id, &
                     & nsnow, 1, nsnow, hist_stomate_snow,32, ave(5), dt,hist_dt)
             END IF
          END DO
       ENDIF

       IF (writehist_deltaC) THEN
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deltaCH4g_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("methanogenesis   "), &
                     & TRIM("gCH4/m**3 air/s   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deltaCH4_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("methanotrophy   "), &
                     & TRIM("gCH4/m**3 air/s   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deltaC1_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("oxic decomposition   "), &
                     & TRIM("gC/m**3/s   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deltaC2_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("methanogenesis   "), &
                     & TRIM("gC/m**3 soil/s   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deltaC3_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("methanotrophy   "), &
                     & TRIM("gC/m**3 soil/s   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
       ENDIF
       IF (writehist_zimovheat) THEN
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("heat_Zimov_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("heating due to decomposition   "), &
                     & TRIM("W/m**3   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
       ENDIF
       IF (writehist_deltaC_litter) THEN
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deltaC_litter_act_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("litter C input to soil active C pool   "), &
                     & TRIM("gC/m**3 soil/s   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
          DO jv = 1, nvm
             IF (permafrost_veg_exists(jv)) THEN
                WRITE(part_str,'(I2)') jv
                IF (jv < 10) part_str(1:1) = '0'
                CALL histdef (hist_id_stom, &
                     & TRIM("deltaC_litter_slo_"//part_str(1:LEN_TRIM(part_str))), &
                     & TRIM("litter C input to soil slow C pool   "), &
                     & TRIM("gC/m**3 soil/s   "), iim,jjm, hist_hori_id, &
                     & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5),dt,hist_dt)
             END IF
          END DO
       ENDIF
       IF (writehist_gascoeff) THEN
          CALL histdef (hist_id_stom, &
               & TRIM("deltaC_litter_pas_"//part_str(1:LEN_TRIM(part_str))), &
               &               TRIM("litter C input to soil passive C pool   "),&
               &               TRIM("gC/m**3 soil/s   "), iim,jjm, hist_hori_id,&
               &               ndeep, 1, ndeep, hist_stomate_deepsoil,32,ave(5),dt, hist_dt)
          DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histdef (hist_id_stom, &
                  & TRIM("totporO2_soil_"//part_str(1:LEN_TRIM(part_str))), &
                  & TRIM("    "), &
                  & TRIM("    "), iim,jjm, hist_hori_id, &
                  & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt,hist_dt)
          END DO

          DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histdef (hist_id_stom, &
                  & TRIM("diffO2_soil_"//part_str(1:LEN_TRIM(part_str))), &
                  & TRIM("    "), &
                  & TRIM("    "), iim,jjm, hist_hori_id, &
                  & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt,hist_dt)
          END DO

          DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histdef (hist_id_stom, &
                  & TRIM("alphaO2_soil_"//part_str(1:LEN_TRIM(part_str))), &
                  & TRIM("    "), &
                  & TRIM("    "), iim,jjm, hist_hori_id, &
                  & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt,hist_dt)
          END DO

          DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histdef (hist_id_stom, &
                  & TRIM("betaO2_soil_"//part_str(1:LEN_TRIM(part_str))), &
                  & TRIM("    "), &
                  & TRIM("    "), iim,jjm, hist_hori_id, &
                  & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt,hist_dt)
          END DO
          DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histdef (hist_id_stom, &
                  & TRIM("totporCH4_soil_"//part_str(1:LEN_TRIM(part_str))), &
                  & TRIM("    "), &
                  & TRIM("    "), iim,jjm, hist_hori_id, &
                  & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt,hist_dt)
          END DO

          DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histdef (hist_id_stom, &
                  & TRIM("diffCH4_soil_"//part_str(1:LEN_TRIM(part_str))), &
                  & TRIM("    "), &
                  & TRIM("    "), iim,jjm, hist_hori_id, &
                  & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt,hist_dt)
          END DO

          DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histdef (hist_id_stom, &
                  & TRIM("alphaCH4_soil_"//part_str(1:LEN_TRIM(part_str))), &
                  & TRIM("    "), &
                  & TRIM("    "), iim,jjm, hist_hori_id, &
                  & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt,hist_dt)
          END DO
          DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histdef (hist_id_stom, &
                  & TRIM("betaCH4_soil_"//part_str(1:LEN_TRIM(part_str))), &
                  & TRIM("    "), &
                  & TRIM("    "), iim,jjm, hist_hori_id, &
                  & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt, hist_dt)
          END DO
       ENDIF

       call histdef (hist_id_stom, &
            & trim("deepC_a_pftmean"), &
            & trim("active pool deep soil (permafrost) carbon, mean of all PFTs"), &
            & trim("gC/m**3   "), iim,jjm, hist_hori_id, &
            & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt, hist_dt)
       call histdef (hist_id_stom, &
            & trim("deepC_s_pftmean"), &
            & trim("slow pool deep soil (permafrost) carbon, mean of all PFTs"), &
            & trim("gC/m**3   "), iim,jjm, hist_hori_id, &
            & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt, hist_dt)
       call histdef (hist_id_stom, &
            & trim("deepC_p_pftmean"), &
            & trim("passive pool deep soil (permafrost) carbon, mean of all PFTs"),&
            & trim("gC/m**3   "), iim,jjm, hist_hori_id, &
            & ndeep, 1, ndeep, hist_stomate_deepsoil,32, ave(5), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &               TRIM("fluxCH4           "), &
            &               TRIM("   "), &
            &               TRIM("gCH4/m**2/day   "), iim,jjm, hist_hori_id, &
            &               nvm, 1, nvm, hist_PFTaxis_id,32, ave(5), dt,hist_dt)
       CALL histdef (hist_id_stom, &
            &               TRIM("febul           "), &
            &               TRIM("   "), &
            &               TRIM("gCH4/m**2/day   "), iim,jjm, hist_hori_id, &
            &               nvm, 1, nvm, hist_PFTaxis_id,32, ave(5), dt,hist_dt)
       CALL histdef (hist_id_stom, &
            &               TRIM("flupmt           "), &
            &               TRIM("   "), &
            &               TRIM("gCH4/m**2/day   "), iim,jjm, hist_hori_id, &
            &               nvm, 1, nvm, hist_PFTaxis_id,32, ave(5), dt, hist_dt)
       CALL histdef (hist_id_stom, &
            &               TRIM("alt           "), &
            &               TRIM("active layer thickness   "), &
            &               TRIM("m   "), iim,jjm, hist_hori_id, &
            &               nvm, 1, nvm, hist_PFTaxis_id,32, ave(5), dt,hist_dt)
       CALL histdef (hist_id_stom, &
            &               TRIM("altmax           "), &
            &               TRIM("max annual alt   "), &
            &               TRIM("m   "), iim,jjm, hist_hori_id, &
            &               nvm, 1, nvm, hist_PFTaxis_id,32, ave(5), dt,hist_dt)
       CALL histdef (hist_id_stom, &
            &               TRIM("sfluxCH4_deep           "), &
            &               TRIM("total surface CH4 flux   "), &
            &               TRIM("gCH4/m**2/sec   "), iim,jjm, hist_hori_id, &
            &               1, 1, 1, -99,32, ave(5), dt, hist_dt)
       CALL histdef (hist_id_stom, &
            &               TRIM("sfluxCO2_deep           "), &
            &               TRIM("total surface CO2 C flux   "), &
            &               TRIM("gC/m**2/sec   "), iim,jjm, hist_hori_id, &
            &               1, 1, 1, -99,32, ave(5), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &               TRIM("z_organic           "), &
            &               TRIM("depth of organic soil   "), &
            &               TRIM("m   "), iim,jjm, hist_hori_id, &
            &               1, 1, 1, -99,32, 'once(scatter(X))', dt, hist_dt)
       CALL histdef (hist_id_stom, &
            &               TRIM("tsurf          "), &
            &               TRIM("surface temp  "), &
            &               TRIM("K  "), iim,jjm, hist_hori_id, &
            &               1, 1, 1, -99,32, ave(5), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &               TRIM("pb          "), &
            &               TRIM("surface pressure  "), &
            &               TRIM("pa   "), iim,jjm, hist_hori_id, &
            &               1, 1, 1, -99,32, ave(5), dt, hist_dt)

       CALL histdef (hist_id_stom, &
            &               TRIM("mu_soil          "), &
            &               TRIM("mu_soil  "), &
            &               TRIM("   "), iim,jjm, hist_hori_id, &
            &               1, 1, 1, -99,32, ave(5), dt, hist_dt)


    ENDIF
    !---------------------------------
  END SUBROUTINE stom_define_history
  !
  SUBROUTINE stom_IPCC_define_history &
       & (hist_id_stom_IPCC, nvm, iim, jjm, dt, &
       &  hist_dt, hist_hori_id, hist_PFTaxis_id)
    ! deforestation axis added as arguments

    !---------------------------------------------------------------------
    !- Tell ioipsl which variables are to be written
    !- and on which grid they are defined
    !---------------------------------------------------------------------
    IMPLICIT NONE
    !-
    !- Input
    !-
    !- File id
    INTEGER(i_std),INTENT(in) :: hist_id_stom_IPCC
    !- number of PFTs
    INTEGER(i_std),INTENT(in) :: nvm
    !- Domain size
    INTEGER(i_std),INTENT(in) :: iim, jjm
    !- Time step of STOMATE (seconds)
    REAL(r_std),INTENT(in)    :: dt
    !- Time step of history file (s)
    REAL(r_std),INTENT(in)    :: hist_dt
    !- id horizontal grid
    INTEGER(i_std),INTENT(in) :: hist_hori_id
    !- id of PFT axis
    INTEGER(i_std),INTENT(in) :: hist_PFTaxis_id
    !-
    !- 1 local
    !-
    !- Character strings to define operations for histdef
    CHARACTER(LEN=40),DIMENSION(max_hist_level) :: ave

    !=====================================================================
    !- 1 define operations
    !=====================================================================
    ave(1) =  'ave(scatter(X))'
    !=====================================================================
    !- 2 surface fields (2d)
    !=====================================================================
    ! Carbon in Vegetation
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cVeg"), &
         &               TRIM("Carbon in Vegetation"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Litter Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cLitter"), &
         &               TRIM("Carbon in Litter Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Soil Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cSoil"), &
         &               TRIM("Carbon in Soil Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Products of Land Use Change
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cProduct"), &
         &               TRIM("Carbon in Products of Land Use Change"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon Mass Variation
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cMassVariation"), &
         &               TRIM("Terrestrial Carbon Mass Variation"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Leaf Area Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("lai"), &
         &               TRIM("Leaf Area Fraction"), &
         &               TRIM("1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Gross Primary Production
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("gpp"), &
         &               TRIM("Gross Primary Production"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Autotrophic Respiration
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("ra"), &
         &               TRIM("Autotrophic Respiration"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Net Primary Production
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("npp"), &
         &               TRIM("Net Primary Production"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Heterotrophic Respiration
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("rh"), &
         &               TRIM("Heterotrophic Respiration"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Emission from Fire
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fFire"), &
         &               TRIM("CO2 Emission from Fire"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! CO2 Flux to Atmosphere from Crop Harvesting
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fHarvest"), &
         &               TRIM("CO2 Flux to Atmosphere from Crop Harvesting"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux to Atmosphere from Land Use Change
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fLuc"), &
         &               TRIM("CO2 Flux to Atmosphere from Land Use Change"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Net Biospheric Production
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nbp"), &
         &               TRIM("Net Biospheric Production"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Total Carbon Flux from Vegetation to Litter
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fVegLitter"), &
         &               TRIM("Total Carbon Flux from Vegetation to Litter"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Total Carbon Flux from Litter to Soil
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("fLitterSoil"), &
         &               TRIM("Total Carbon Flux from Litter to Soil"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Carbon in Leaves
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cLeaf"), &
         &               TRIM("Carbon in Leaves"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Wood
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cWood"), &
         &               TRIM("Carbon in Wood"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Roots
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cRoot"), &
         &               TRIM("Carbon in Roots"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Other Living Compartments
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cMisc"), &
         &               TRIM("Carbon in Other Living Compartments"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Carbon in Above-Ground Litter
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cLitterAbove"), &
         &               TRIM("Carbon in Above-Ground Litter"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Below-Ground Litter
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cLitterBelow"), &
         &               TRIM("Carbon in Below-Ground Litter"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Fast Soil Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cSoilFast"), &
         &               TRIM("Carbon in Fast Soil Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Medium Soil Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cSoilMedium"), &
         &               TRIM("Carbon in Medium Soil Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Carbon in Slow Soil Pool
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("cSoilSlow"), &
         &               TRIM("Carbon in Slow Soil Pool"), &
         &               TRIM("kg C m-2"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    !- 3 PFT: 3rd dimension
    ! Fractional Land Cover of PFT
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("landCoverFrac"), &
         &               TRIM("Fractional Land Cover of PFT"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               nvm,1,nvm, hist_PFTaxis_id,32, ave(1), dt, hist_dt)


    ! Total Primary Deciduous Tree Cover Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("treeFracPrimDec"), &
         &               TRIM("Total Primary Deciduous Tree Cover Fraction"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Total Primary Evergreen Tree Cover Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("treeFracPrimEver"), &
         &               TRIM("Total Primary Evergreen Tree Cover Fraction"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    ! Total C3 PFT Cover Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("c3PftFrac"), &
         &               TRIM("Total C3 PFT Cover Fraction"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Total C4 PFT Cover Fraction
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("c4PftFrac"), &
         &               TRIM("Total C4 PFT Cover Fraction"), &
         &               TRIM("%"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Growth Autotrophic Respiration
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("rGrowth"), &
         &               TRIM("Growth Autotrophic Respiration"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Maintenance Autotrophic Respiration
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("rMaint"), &
         &               TRIM("Maintenance Autotrophic Respiration"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux from Atmosphere due to NPP Allocation to Leaf
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nppLeaf"), &
         &               TRIM("CO2 Flux from Atmosphere due to NPP Allocation to Leaf"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux from Atmosphere due to NPP Allocation to Wood
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nppWood"), &
         &               TRIM("CO2 Flux from Atmosphere due to NPP Allocation to Wood"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! CO2 Flux from Atmosphere due to NPP Allocation to Root
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nppRoot"), &
         &               TRIM("CO2 Flux from Atmosphere due to NPP Allocation to Root"), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)
    ! Net Carbon Mass Flux out of Atmophere due to Net Ecosystem Productivity on Land.
    CALL histdef (hist_id_stom_IPCC, &
         &               TRIM("nep"), &
         &               TRIM("Net Carbon Mass Flux out of Atmophere due to Net Ecosystem Productivity."), &
         &               TRIM("kg C m-2 s-1"), iim,jjm, hist_hori_id, &
         &               1,1,1, -99,32, ave(1), dt, hist_dt)

    CALL histdef(hist_id_stom_IPCC, 'RESOLUTION_X', 'E-W resolution', 'm', &
         & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
    CALL histdef(hist_id_stom_IPCC, 'RESOLUTION_Y', 'N-S resolution', 'm', &
         & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
    CALL histdef(hist_id_stom_IPCC, 'CONTFRAC', 'Continental fraction', '1', &
         & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)
    CALL histdef(hist_id_stom_IPCC, 'Areas', 'Mesh areas', 'm2', &
         & iim,jjm, hist_hori_id, 1,1,1, -99, 32, 'once(scatter(X))', dt, hist_dt)

    !---------------------------------
  END SUBROUTINE stom_IPCC_define_history
END MODULE intersurf
