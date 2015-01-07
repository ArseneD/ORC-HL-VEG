MODULE watchout

  USE defprec
  USE mod_orchidee_para 
  USE constantes
  USE netcdf

  PRIVATE
  PUBLIC :: watchout_init, watchout_write_p, watchout_close
!watchout_write ??

  LOGICAL,SAVE,PUBLIC             :: ok_watchout = .FALSE.
!$OMP THREADPRIVATE(ok_watchout)
  REAL, SAVE,PUBLIC               :: dt_watch = zero
!$OMP THREADPRIVATE(dt_watch)
  INTEGER, SAVE,PUBLIC            :: last_action_watch = 0, &
       & last_check_watch = 0
!$OMP THREADPRIVATE(last_action_watch, last_check_watch)
  CHARACTER(LEN=80),SAVE, PUBLIC   :: watchout_file
!$OMP THREADPRIVATE(watchout_file)
  ! At module level we need the ids of the variables for the ORCHIDEE_WATCH. They will be 
  ! shared by the watchout_init and watchout_write routines.
  ! The flag which will control all this is watchout
  !
  INTEGER(i_std),SAVE      :: time_id, timestp_id
!$OMP THREADPRIVATE(time_id, timestp_id)
  INTEGER(i_std),SAVE      :: watchfid, zlevid, soldownid, rainfid, snowfid, lwradid, &
       & psolid, tairid, eairid, qairid, uid, vid, &
       & solnetid, petAcoefid, peqAcoefid, petBcoefid, peqBcoefid, cdragid, ccanopyid
!$OMP THREADPRIVATE(watchfid, zlevid, soldownid, rainfid)
!$OMP THREADPRIVATE(snowfid, lwradid, psolid, tairid) 
!$OMP THREADPRIVATE( eairid, qairid, uid, vid) 
!$OMP THREADPRIVATE(solnetid, petAcoefid, peqAcoefid, petBcoefid) 
!$OMP THREADPRIVATE(peqBcoefid, cdragid, ccanopyid) 
  INTEGER(i_std),SAVE      :: watchoffset
!$OMP THREADPRIVATE(watchoffset) 
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_zlev
!$OMP THREADPRIVATE(sum_zlev) 
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_u, sum_v
!$OMP THREADPRIVATE(sum_u, sum_v)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_qair
!$OMP THREADPRIVATE(sum_qair) 
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_temp_air
!$OMP THREADPRIVATE(sum_temp_air) 
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_epot_air
!$OMP THREADPRIVATE(sum_epot_air) 
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_ccanopy
!$OMP THREADPRIVATE(sum_ccanopy) 
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_cdrag
!$OMP THREADPRIVATE(sum_cdrag) 
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_petAcoef, sum_peqAcoef, sum_petBcoef, sum_peqBcoef
!$OMP THREADPRIVATE(sum_petAcoef, sum_peqAcoef) 
!$OMP THREADPRIVATE(sum_petBcoef, sum_peqBcoef)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_rain, sum_snow
!$OMP THREADPRIVATE(sum_rain, sum_snow)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_lwdown
!$OMP THREADPRIVATE(sum_lwdown)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_swnet
!$OMP THREADPRIVATE(sum_swnet)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_swdown
!$OMP THREADPRIVATE(sum_swdown)
  REAL(r_std), ALLOCATABLE, DIMENSION(:), SAVE, PUBLIC  :: sum_pb
!$OMP THREADPRIVATE(sum_pb)
  !! Short wave mean : compute with solar angle, as in readdim2
!!$  REAL(r_std), ALLOCATABLE, DIMENSION(:,:), SAVE, PUBLIC  :: sinang, mean_sinang
!!$  INTEGER(i_std), ALLOCATABLE, DIMENSION(:,:), SAVE, PUBLIC :: isinang

  REAL(r_std), SAVE, PUBLIC :: dt_split_watch
!$OMP THREADPRIVATE(dt_split_watch) 
  !! mean julian time for each time step
  REAL(r_std), SAVE, PUBLIC :: julian_watch
!$OMP THREADPRIVATE(julian_watch)
CONTAINS 
SUBROUTINE watchout_init(iim, jjm, kjpindex, igmax, date0, itau, dt, kindex, lon, lat, lev0)
  !
  IMPLICIT NONE
  !
  ! This routine will allow to set up forcing files for ORCHIDEE. The idea is that
  ! during a coupled simulation one write's out all the forcing so that ORCHIDEE can
  ! can be re-run (to equilibrium or for sensitivity) afterwards.
  !
  ! INPUT
  INTEGER(i_std), INTENT(in)   :: iim, jjm, igmax, kjpindex
  REAL(r_std), INTENT(in)      :: date0, dt
  INTEGER(i_std), INTENT(in)   :: itau, kindex(igmax)
  REAL(r_std), INTENT(in)      :: lon(iim,jjm), lat(iim,jjm), lev0
  !
  ! OUTPUT
  !
  !
  ! LOCAL
  !
  INTEGER(i_std)   :: iret, nlonid1, nlatid1, nlevid1, fid, nlandid1, tdimid1
  INTEGER(i_std)   :: dims(3)
  INTEGER(i_std)   :: nlonid, nlatid, nlevid, nlandid, varid, contid, resolxid, resolyid
  INTEGER(i_std), DIMENSION(8)  ::  neighid
  REAL(r_std)      :: lon_min, lon_max, lat_min, lat_max, lev_min, lev_max
  INTEGER(i_std)   :: yy, mm, dd, hh, mn, i, j, ig, direction
  REAL(r_std)      :: ss
  REAL(r_std),ALLOCATABLE :: tmpdata(:,:)
  CHARACTER(LEN=3)  :: cal(12)
  CHARACTER(LEN=10) :: today, att, axx
  CHARACTER(LEN=30) :: str30
  CHARACTER(LEN=70) :: str70, var, unit, titre, assoc
  CHARACTER(LEN=80) :: stamp, lon_name, lat_name, land_name,time_name
  !    
  INTEGER,PARAMETER :: kind_r_watch=nf90_real8
  !
  ! Only root proc write watchout file
  IF (is_root_prc) THEN

     cal(1) = 'JAN'
     cal(2) = 'FEB'
     cal(3) = 'MAR'
     cal(4) = 'APR'
     cal(5) = 'MAY'
     cal(6) = 'JUN'
     cal(7) = 'JUL'
     cal(8) = 'AUG'
     cal(9) = 'SEP'
     cal(10) = 'OCT'
     cal(11) = 'NOV'
     cal(12) = 'DEC'
     !
     iret = NF90_CREATE (TRIM(watchout_file), NF90_CLOBBER, fid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not create file :',TRIM(watchout_file), &
             &          '(Problem with disk place or filename ?)')
     ENDIF
     !
     !   Dimensions
     !
     iret = NF90_DEF_DIM(fid, 'x', iim, nlonid1)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Dimension "x" can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_DEF_DIM(fid, 'y', jjm, nlatid1)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Dimension "y" can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_DEF_DIM(fid, 'z', 1, nlevid1)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Dimension "z" can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_DEF_DIM(fid, 'land', igmax, nlandid1)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Dimension "land" can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_DEF_DIM(fid, 'tstep', NF90_UNLIMITED, tdimid1)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Dimension "tstep" can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     !   Coordinate  VARIABLES
     !
     dims(1) = nlonid1
     dims(2) = nlatid1
     !
     lon_name = 'nav_lon'
     iret = NF90_DEF_VAR(fid, lon_name, kind_r_watch, dims(1:2), nlonid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//lon_name//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlonid, 'units', "degrees_east") 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lon_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     lon_min = -180.
     lon_max = 180.
     !
     iret = NF90_PUT_ATT(fid, nlonid, 'valid_min', lon_min)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lon_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlonid, 'valid_max', lon_max)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lon_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_ATT(fid, nlonid, 'long_name', "Longitude")
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lon_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     lat_name = 'nav_lat'
     iret = NF90_DEF_VAR(fid, lat_name, kind_r_watch, dims(1:2), nlatid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//lat_name//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlatid, 'units', "degrees_north")
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lat_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     lat_max = 90.
     lat_min = -90.
     !
     iret = NF90_PUT_ATT(fid, nlatid, 'valid_min', lat_min)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lat_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlatid, 'valid_max', lat_max)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lat_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlatid, 'long_name', "Latitude")
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lat_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     lat_name = 'level'
     iret = NF90_DEF_VAR(fid, lat_name, kind_r_watch,(/ nlevid1 /), nlevid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//lat_name//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlevid, 'units', "m")
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lat_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     lev_max = lev0
     lev_min = lev0
     !
     iret = NF90_PUT_ATT(fid, nlevid, 'valid_min', lev_min)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lat_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlevid, 'valid_max', lev_max)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lat_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlevid, 'long_name', "Vertical levels")
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//lat_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     !
     land_name = 'land'
     iret = NF90_DEF_VAR(fid, land_name, NF90_INT, (/ nlandid1 /), nlandid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//land_name//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, nlandid, 'compress', "y x") 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//land_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     !   Time in real days !
     !
     time_name = 'time'
     iret = NF90_DEF_VAR(fid, time_name, kind_r_watch, (/ tdimid1 /), time_id)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//time_name//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     ! Compute an itau offset so that we can relate the itau of the model
     ! to the position in the file
     !
     watchoffset = itau
     !
    CALL ju2ymds(date0, yy, mm, dd, ss)
    hh = INT(ss/3600.)
    ss = ss - hh*3600.
     mn = INT(ss/60.) 
     ss = ss - mn*60.
    WRITE(str70,7000) yy, mm, dd, hh, mn, INT(ss)
!!MM : Time axis by month : 
!!$     hh = INT(sec/3600.)
!!$     ss = sec - hh*3600.
!!$     mn = INT(ss/60.) 
!!$     ss = ss - mn*60.
!!$     WRITE(str70,7000) year, month, day, hh, mn, INT(ss)
     !
     iret = NF90_PUT_ATT(fid, time_id, 'units', TRIM(str70))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     CALL ioget_calendar(str30)
     iret = NF90_PUT_ATT(fid, time_id, 'calendar', TRIM(str30))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_ATT(fid, time_id, 'title', 'Time')
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_ATT(fid, time_id, 'long_name', 'Time axis')
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
    WRITE(str70,7001) yy, cal(mm), dd, hh, mn, INT(ss)
!!MM : Time axis by month : 
!!!$     WRITE(str70,7001) year, CAL(month), day, hh, mn, INT(ss)
     iret = NF90_PUT_ATT(fid, time_id, 'time_origin', TRIM(str70))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     !   Time steps
     !
     time_name = 'timestp'
     iret = NF90_DEF_VAR(fid, time_name, NF90_INT, (/ tdimid1 /), timestp_id)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//time_name//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
    WRITE(str70,7002) yy, mm, dd, hh, mn, INT(ss)
!!MM : Time axis by month : 
!!!$     WRITE(str70,7002) year, month, day, hh, mn, INT(ss)
     iret = NF90_PUT_ATT(fid, timestp_id, 'units', TRIM(str70))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_ATT(fid, timestp_id, 'title', 'Time steps')
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_ATT(fid, timestp_id, 'tstep_sec', dt)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_ATT(fid, timestp_id, 'long_name', 'Time step axis')
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
    WRITE(str70,7001) yy, cal(mm), dd, hh, mn, INT(ss)
!!MM : Time axis by month : 
!!!$     WRITE(str70,7001) year, CAL(month), day, hh, mn, INT(ss)
     iret = NF90_PUT_ATT(fid, timestp_id, 'time_origin', TRIM(str70))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//time_name//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
7000 FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
7001 FORMAT(' ', I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
7002 FORMAT('timesteps since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
     !

     dims(1) = nlandid1
     dims(2) = tdimid1

     assoc = 'time (nav_lat nav_lon)'
     axx='TYX'
     !
     var = 'SWdown'
     unit = 'W/m^2'
     titre = 'Surface incident shortwave radiation'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     soldownid = varid
     !
     var = 'SWnet'
     unit = 'W/m^2'
     titre = 'Net surface short-wave flux'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     solnetid = varid
     !
     var = 'Rainf'
     unit = 'Kg/m^2s'
     titre = 'Rainfall rate'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     rainfid = varid
     !
     var = 'Snowf'
     unit = 'Kg/m^2s'
     titre = 'Snowfall rate'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     snowfid = varid
     !
     var = 'LWdown'
     unit = 'W/m^2'
     titre = 'Surface incident longwave radiation'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     lwradid = varid
     !
     var = 'PSurf'
     unit = 'Pa'
     titre = 'Surface pressure'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     psolid = varid
     !
     !
     !  3D Variables to be written
     !
     dims(1) = nlandid1
     dims(2) = nlevid1
     dims(3) = tdimid1
     !
     assoc = 'time level (nav_lat nav_lon)'
     axx='TZYX'
     !
     lat_name = 'levels'
     iret = NF90_DEF_VAR(fid, lat_name, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', "m")
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', "Vertical levels")
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     lev_min = 2.
     lev_max = 100.
     iret = NF90_PUT_ATT(fid, varid, 'valid_min', lev_min)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'valid_max', lev_max)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     zlevid = varid
     !
     !
     var = 'Tair'
     unit = 'K'
     titre = 'Near surface air temperature'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     tairid = varid
     !
     var = 'Eair'
     unit = 'J/m^2'
     titre = 'Air potential energy'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     eairid = varid
     !
     var = 'Qair'
     unit = 'Kg/Kg'
     titre = 'Near surface specific humidity'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     qairid = varid
     !
     var = 'Wind_N'
     unit = 'm/s'
     titre = 'Near surface northward wind component'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     uid = varid
     !
     var = 'Wind_E'
     unit = 'm/s'
     titre = 'Near surface eastward wind component'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     vid = varid
     !
     var = 'petAcoef'
     unit = '-'
     titre = 'Coeficients A from the PBL resolution for T'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     petAcoefid = varid
     !
     var = 'peqAcoef'
     unit = '-'
     titre = 'Coeficients A from the PBL resolution for q'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     peqAcoefid = varid
     !
     var = 'petBcoef'
     unit = '-'
     titre = 'Coeficients B from the PBL resolution for T'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     petBcoefid = varid
     !
     var = 'peqBcoef'
     unit = '-'
     titre = 'Coeficients B from the PBL resolution for q'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     peqBcoefid = varid
     !
     var = 'cdrag'
     unit = '-'
     titre = 'Surface drag'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     cdragid = varid
     !
     !
     var = 'ccanopy'
     unit = '-'
     titre = 'CO2 concentration in the canopy'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:3), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     ccanopyid = varid
     !
     !
     ! Time fixed variable
     !
     dims(1) = nlonid1
     dims(2) = nlatid1
     !
     var = 'contfrac'
     unit = '-'
     titre = 'Fraction of continent'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     contid=varid
     !
     !
     var = 'neighboursNN'
     unit = '-'
     titre = 'indices of North neighbours of each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     neighid(1)=varid
     !
     var = 'neighboursNE'
     unit = '-'
     titre = 'indices of North-East neighbours of each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     neighid(2)=varid
     !
     var = 'neighboursEE'
     unit = '-'
     titre = 'indices of East neighbours of each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     neighid(3)=varid
     !
     var = 'neighboursSE'
     unit = '-'
     titre = 'indices of South-East neighbours of each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     neighid(4)=varid
     !
     var = 'neighboursSS'
     unit = '-'
     titre = 'indices of South neighbours of each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     neighid(5)=varid
     !
     var = 'neighboursSW'
     unit = '-'
     titre = 'indices of South-West neighbours of each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     neighid(6)=varid
     !
     var = 'neighboursWW'
     unit = '-'
     titre = 'indices of West neighbours of each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     neighid(7)=varid
     !
     var = 'neighboursNW'
     unit = '-'
     titre = 'indices of North-West neighbours of each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     neighid(8)=varid
     !
     !
     var = 'resolutionX'
     unit = 'm'
     titre = 'resolution in x at each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     resolxid=varid
     !
     var = 'resolutionY'
     unit = 'm'
     titre = 'resolution in y at each grid point'
     assoc = 'nav_lat nav_lon'
     axx='YX'
     iret = NF90_DEF_VAR(fid,var, kind_r_watch, dims(1:2), varid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &         'Variable '//var//' can not be defined for the file : ', &
             &         TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'axis',  TRIM(axx) )
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'units', TRIM(unit)) 
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'long_name', TRIM(titre))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'associate', TRIM(assoc))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, varid, 'missing_value', undef_sechiba)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add attribut to variable '//var//' for the file :', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     resolyid=varid
     !
     !
     !  Global attributes
     !
     CALL DATE_AND_TIME(today, att)
     stamp = "Forcing generated by intersurf in a previous run "//today(1:LEN_TRIM(today))//" at "//att(1:LEN_TRIM(att))
     iret = NF90_PUT_ATT(fid, NF90_GLOBAL, 'Conventions', "GDT 1.2")
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add global attribut to the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, NF90_GLOBAL, 'file_name', TRIM(watchout_file))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add global attribut to the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     iret = NF90_PUT_ATT(fid, NF90_GLOBAL, 'production', TRIM(stamp))
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not add global attribut to the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_ENDDEF(fid)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not end definitions in the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     !    Write coordinates
     !
     iret = NF90_PUT_VAR(fid, nlonid, lon)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not put variable nav_lon  in the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_VAR(fid, nlatid, lat)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not put variable nav_lat  in the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_VAR(fid, nlevid, lev0)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not put variable level  in the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     iret = NF90_PUT_VAR(fid, nlandid, kindex)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not put variable land  in the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     IF ( .NOT. ALLOCATED(tmpdata)) THEN
        ALLOCATE(tmpdata(iim,jjm))
     ENDIF
     !
     tmpdata(:,:) = undef_sechiba
     DO ig=1,igmax

        j = ((kindex(ig)-1)/iim) + 1
        i = (kindex(ig) - (j-1)*iim)

        tmpdata(i,j) = contfrac_g(ig)

     ENDDO
     iret = NF90_PUT_VAR(fid, contid, tmpdata)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not put variable contfrac  in the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !    
     DO direction=1,8
        tmpdata(:,:) = undef_sechiba
        DO ig=1,igmax

           j = ((kindex(ig)-1)/iim) + 1
           i = (kindex(ig) - (j-1)*iim)

           tmpdata(i,j) = REAL( neighbours_g(ig,direction) )

        ENDDO
        iret = NF90_PUT_VAR(fid, neighid(direction), tmpdata)
        IF (iret /= NF90_NOERR) THEN
           CALL ipslerr_p (3,'watchout_init', &
                &             'Could not put variable neighbours  in the file : ', &
                &             TRIM(watchout_file),'(Solution ?)')
        ENDIF
     ENDDO
     !
     tmpdata(:,:) = undef_sechiba
     DO ig=1,igmax

        j = ((kindex(ig)-1)/iim) + 1
        i = (kindex(ig) - (j-1)*iim)

        tmpdata(i,j) = resolution_g(ig,1)

     ENDDO
     iret = NF90_PUT_VAR(fid, resolxid, tmpdata)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not put variable resolutionx  in the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     tmpdata(:,:) = undef_sechiba
     DO ig=1,igmax

        j = ((kindex(ig)-1)/iim) + 1
        i = (kindex(ig) - (j-1)*iim)

        tmpdata(i,j) = resolution_g(ig,2)

     ENDDO
     iret = NF90_PUT_VAR(fid, resolyid, tmpdata)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr_p (3,'watchout_init', &
             &          'Could not put variable resolutiony  in the file : ', &
             &          TRIM(watchout_file),'(Solution ?)')
     ENDIF
     !
     DEALLOCATE(tmpdata)
     !
     watchfid = fid
     !
  ENDIF
  !
  ALLOCATE(sum_zlev(kjpindex))
  ALLOCATE(sum_u(kjpindex), sum_v(kjpindex))
  ALLOCATE(sum_qair(kjpindex))
  ALLOCATE(sum_temp_air(kjpindex))
  ALLOCATE(sum_epot_air(kjpindex))
  ALLOCATE(sum_ccanopy(kjpindex))
  ALLOCATE(sum_cdrag(kjpindex))
  ALLOCATE(sum_petAcoef(kjpindex), sum_peqAcoef(kjpindex), sum_petBcoef(kjpindex), sum_peqBcoef(kjpindex))
  ALLOCATE(sum_rain(kjpindex), sum_snow(kjpindex))
  ALLOCATE(sum_lwdown(kjpindex))
  ALLOCATE(sum_swnet(kjpindex))
  ALLOCATE(sum_swdown(kjpindex))
  ALLOCATE(sum_pb(kjpindex))
!!$  ALLOCATE(mean_sinang(iim,jjm))
!!$  ALLOCATE(sinang(iim,jjm))
!!$  ALLOCATE(isinang(iim,jjm))

  sum_zlev(:) = zero
  sum_u(:) = zero
  sum_v(:) = zero
  sum_qair(:) = zero
  sum_temp_air(:) = zero
  sum_epot_air(:) = zero
  sum_ccanopy(:) = zero
  sum_cdrag(:) = zero
  sum_petAcoef(:) = zero
  sum_peqAcoef(:) = zero
  sum_petBcoef(:) = zero
  sum_peqBcoef(:) = zero
  sum_rain(:) = zero
  sum_snow(:) = zero
  sum_lwdown(:) = zero
  sum_swnet(:) = zero
  sum_swdown(:) = zero
  sum_pb(:) = zero

!!$  mean_sinang(:,:) = zero
!!$  sinang(:,:) = zero
!!$  isinang(:,:) = dt_split_watch

END SUBROUTINE watchout_init
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE watchout_write_p(igmax, itau, dt, levels, &
       &                        soldown, rain, snow, lwdown, psurf, temp, eair, qair, u, v, &
       &                        solnet, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy )
    !
    ! 
    IMPLICIT NONE
    !
    ! INPUT
    !
    INTEGER(i_std), INTENT(in) :: igmax, itau
    REAL(r_std), INTENT(inout) :: levels(igmax)
    REAL(r_std), DIMENSION(igmax), INTENT(inout) :: soldown, rain, snow, lwdown, psurf, &
         &                                      solnet, petAcoef, peqAcoef, petBcoef, peqBcoef, &
         &                                      cdrag, ccanopy
    REAL(r_std), INTENT(inout) :: temp(igmax), eair(igmax), qair(igmax), u(igmax), v(igmax)
    REAL(r_std), INTENT(in) :: dt
    !
    ! LOCAL
    !
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: levels_g
!$OMP THREADPRIVATE(levels_g)
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:) :: soldown_g, rain_g, snow_g, lwdown_g, psurf_g, &
         &                             solnet_g, petAcoef_g, peqAcoef_g, petBcoef_g, peqBcoef_g, cdrag_g, &
         &                             temp_g, eair_g, qair_g, u_g, v_g, ccanopy_g
!$OMP THREADPRIVATE(soldown_g, rain_g, snow_g, lwdown_g)
!$OMP THREADPRIVATE(psurf_g,solnet_g, petAcoef_g, peqAcoef_g)
!$OMP THREADPRIVATE( petBcoef_g, peqBcoef_g, cdrag_g,temp_g)
!$OMP THREADPRIVATE(eair_g, qair_g, u_g, v_g, ccanopy_g)
    !
    LOGICAL, SAVE                          :: is_first_time=.TRUE.
!$OMP THREADPRIVATE(is_first_time)
    INTEGER(i_std)                         :: ier

    IF (is_first_time .AND. is_root_prc) THEN
       ALLOCATE(levels_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in levels_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(soldown_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in soldown_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(rain_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in rain_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(snow_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in snow_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(lwdown_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in lwdown_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(psurf_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in psurf_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(solnet_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in solnet_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(petAcoef_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in petAcoef_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(peqAcoef_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in peqAcoef_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(petBcoef_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in petBcoef_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(peqBcoef_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in peqBcoef_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(cdrag_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in cdrag_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(temp_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in temp_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(eair_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in eair_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(qair_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in qair_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(u_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in u_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(v_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in v_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF

       ALLOCATE(ccanopy_g(nbp_glo),stat=ier)
       IF (ier .NE. 0) THEN
          WRITE (numout,*) ' error in ccanopy_g allocation. We stop. We need iim words = ',nbp_glo
          STOP 'watchout_write_p'
       ENDIF
    ENDIF
    is_first_time=.FALSE.

    CALL gather(levels,levels_g)
    CALL gather(soldown,soldown_g)
    CALL gather(rain,rain_g)
    CALL gather(snow,snow_g)
    CALL gather(lwdown,lwdown_g)
    CALL gather(psurf,psurf_g)
    CALL gather(solnet,solnet_g)
    CALL gather(petAcoef,petAcoef_g)
    CALL gather(peqAcoef,peqAcoef_g)
    CALL gather(petBcoef,petBcoef_g)
    CALL gather(peqBcoef,peqBcoef_g)
    CALL gather(cdrag,cdrag_g)
    CALL gather(temp,temp_g)
    CALL gather(eair,eair_g)
    CALL gather(qair,qair_g)
    CALL gather(u,u_g)
    CALL gather(v,v_g)
    CALL gather(ccanopy,ccanopy_g)

    IF (is_root_prc) THEN
      CALL watchout_write(nbp_glo, itau, dt, levels_g, &
       &                      soldown_g, rain_g, snow_g, lwdown_g, psurf_g, temp_g, eair_g, qair_g, u_g, v_g, &
       &                      solnet_g, petAcoef_g, peqAcoef_g, petBcoef_g, peqBcoef_g, cdrag_g, ccanopy_g )
    ENDIF

    levels(:) = zero
    soldown(:) = zero
    rain(:) = zero
    snow(:) = zero
    lwdown(:) = zero
    psurf(:) = zero
    solnet(:) = zero
    petAcoef(:) = zero
    peqAcoef(:) = zero
    petBcoef(:) = zero
    peqBcoef(:) = zero
    cdrag(:) = zero
    temp(:) = zero
    eair(:) = zero
    qair(:) = zero
    u(:) = zero
    v(:) = zero
    ccanopy(:) = zero

!!$    mean_sinang(:,:) = zero
!!$    isinang(:,:) = dt_split_watch

  END SUBROUTINE watchout_write_p

  SUBROUTINE watchout_write(igmax, itau, dt, levels, &
       &                        soldown, rain, snow, lwdown, psurf, temp, eair, qair, u, v, &
       &                        solnet, petAcoef, peqAcoef, petBcoef, peqBcoef, cdrag, ccanopy )
    !
    ! 
    IMPLICIT NONE
    !
    ! This subroutine will write to the file the current fields which force the
    ! land-surface scheme. It will be in exactly the same format as the other forcing
    ! files, i.e. ALMA convention !
    !
    !
    ! INPUT
    !
    INTEGER(i_std) :: igmax, itau
    REAL(r_std) :: levels(igmax)
    REAL(r_std), DIMENSION(igmax), INTENT(in) :: soldown, rain, snow, lwdown, psurf, &
         &                                      solnet, petAcoef, peqAcoef, petBcoef, peqBcoef, &
         &                                      cdrag, ccanopy
    REAL(r_std) :: temp(igmax), eair(igmax), qair(igmax), u(igmax), v(igmax)
    REAL(r_std) :: dt
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iret
    INTEGER(i_std) :: corner(3), edges(3)
    REAL(r_std) :: timestp
    LOGICAL    :: check=.FALSE.
    REAL(r_std),ALLOCATABLE :: tmpdata(:)
    INTEGER(i_std) :: corner_tstp
    !
    ! For dt_watch non equal to dt :
    !
    corner_tstp = NINT((itau - watchoffset)/dt_split_watch)
    !
    corner(1) = corner_tstp
    edges(1) = 1
    IF ( check ) &
         WRITE(numout,*) 'watchout_write corners, edges : ', corner(1), edges(1)
    !
    timestp = itau/dt_split_watch
!!MM : Time axis by month : 
!!!$    timestp = (itau - watchoffset)/dt_split_watch
    IF ( check ) &
         WRITE(numout,*) "watchout_write : timestp = ",timestp
    iret = NF90_PUT_VAR(watchfid, timestp_id, (/ timestp /), &
         &              start=(/ corner(1) /), count=(/ edges(1) /))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable timestp  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    timestp=timestp*dt_watch
    IF ( check ) &
         WRITE(numout,*) "watchout_write : time = ",timestp
    iret = NF90_PUT_VAR(watchfid, time_id, (/ timestp /), &
         &              start=(/ corner(1) /), count=(/ edges(1) /))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable time  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    corner(1) = 1
    edges(1) = igmax
    corner(2) = corner_tstp
    edges(2) = 1
    !
    IF ( .NOT. ALLOCATED(tmpdata)) THEN
       ALLOCATE(tmpdata(igmax))
    ENDIF
    !
    ! 2D
    IF ( check ) THEN
       WRITE(numout,*) '--',itau, ' SOLDOWN : ', MINVAL(soldown), MAXVAL(soldown)
    ENDIF
    iret = NF90_PUT_VAR(watchfid, soldownid, soldown, start=corner(1:2), count=edges(1:2))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable SWdown  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    iret = NF90_PUT_VAR(watchfid, solnetid, solnet, start=corner(1:2), count=edges(1:2))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable SWnet  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    ! Bring back to kg/m^2/s
    !
    tmpdata = rain/dt
    IF ( check ) THEN
       WRITE(numout,*) '--',itau, ' RAIN : ', MINVAL(tmpdata), MAXVAL(tmpdata)
    ENDIF
    iret = NF90_PUT_VAR(watchfid, rainfid, tmpdata, start=corner(1:2), count=edges(1:2))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable Rainf  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    tmpdata = snow/dt
    iret = NF90_PUT_VAR(watchfid, snowfid, tmpdata, start=corner(1:2), count=edges(1:2))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable Snowf  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    iret = NF90_PUT_VAR(watchfid, lwradid, lwdown, start=corner(1:2), count=edges(1:2)) 
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable LWdown  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    !  Bring back to Pa
    !
    tmpdata = psurf*100.
    iret = NF90_PUT_VAR(watchfid, psolid, tmpdata, start=corner(1:2), count=edges(1:2))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable PSurf  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    ! 3D
    corner(2) = 1
    edges(2) = 1
    corner(3) = corner_tstp
    edges(3) = 1
    !
    iret = NF90_PUT_VAR(watchfid, zlevid, levels, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable levels  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    iret = NF90_PUT_VAR(watchfid, tairid, temp, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable Tair  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, eairid, eair, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable Eair  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, qairid, qair, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable Qair  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, uid, u, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable Wind_N  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, vid, v, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable Wind_E  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
    iret = NF90_PUT_VAR(watchfid, petAcoefid, petAcoef, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable petAcoef  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, peqAcoefid, peqAcoef, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable peqAcoef  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, petBcoefid, petBcoef, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable petBcoef  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, peqBcoefid, peqBcoef, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable peqBcoef  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, cdragid, cdrag, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable cdrag  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    iret = NF90_PUT_VAR(watchfid, ccanopyid, ccanopy, start=corner(1:3), count=edges(1:3))
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_write', &
 &          'Could not put variable ccanopy  in the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF

    DEALLOCATE(tmpdata)
    !
  END SUBROUTINE watchout_write
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  SUBROUTINE watchout_close()
    !
    !  Close the watch files
    !
    IMPLICIT NONE
    !
    ! LOCAL
    !
    INTEGER(i_std) :: iret
    LOGICAL       :: check = .FALSE.
    !
    !
    IF ( check )  THEN
       WRITE(numout,*) 'watchout_close : closing file : ', watchfid
    ENDIF
    iret = NF90_CLOSE(watchfid)
    IF (iret /= NF90_NOERR) THEN
       CALL ipslerr_p (3,'watchout_close','Could not close the file : ', &
 &          TRIM(watchout_file),'(Solution ?)')
    ENDIF
    !
  END SUBROUTINE watchout_close

END MODULE watchout
