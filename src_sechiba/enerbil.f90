!  ==============================================================================================================================
!  MODULE		 			: enerbil
!
!  CONTACT		 			: orchidee-help _at_ ipsl.jussieu.fr
!
!  LICENCE	 			        : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF					This module computes the energy balance on 
!! continental surfaces. The module contains the following subroutines: enerbil_main, enerbil_init,
!! enerbil_clear, enerbil_var_init, enerbil_begin, enerbil_surftemp, enerbil_flux, enerbil_evapveg,
!! enerbil_fusion and enerbil_t2mdiag
!!
!!\n DESCRIPTION                                : 
!! \n
!! \latexonly 
!!     \input{enerbil_intro2.tex}
!! \endlatexonly
!!
!! IMPORTANT NOTE: The coefficients A and B are defined differently than those in the referenced 
!! literature and from those in the code and documentation of the atmospheric model LMDZ. For the
!! avoidance of doubt, the coefficients as described here always refer to the ORCHIDEE coefficients, 
!! and are denoted as such (with the marker: ORC). The re-definition of the coefficients takes place 
!! within LMDZ before they are passed to ORCHIDEE. The following sequence of expressions is to be 
!! found within the LMDZ module 'surf_land_orchidee':\n
!!
!! \latexonly 
!!     \input{surflandLMDZ1.tex}
!!     \input{surflandLMDZ2.tex}
!!     \input{surflandLMDZ3.tex}
!!     \input{surflandLMDZ4.tex}
!! \endlatexonly
!! \n
!!
!! \latexonly 
!!     \input{enerbil_symbols.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S)                             : None
!!
!! REFERENCE(S)	                                : None 
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_sechiba/enerbil.f90 $
!! $Date: 2014-07-08 07:51:19 +0200 (Tue, 08 Jul 2014) $
!! $Revision: 2222 $
!! \n
!_ ================================================================================================================================

MODULE enerbil

  ! routines called : restput, restget
  !
  ! modules used
  USE ioipsl
  USE xios_orchidee
  USE ioipsl_para 
  USE constantes
  USE pft_parameters
  USE qsat_moisture
  USE sechiba_io

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: enerbil_main, enerbil_fusion,enerbil_clear

  ! variables used inside enerbil module : declaration and initialisation
 
  LOGICAL, SAVE                              :: l_first_enerbil=.TRUE.  !! Initialisation has to be done one time
!$OMP THREADPRIVATE(l_first_enerbil)
  ! one dimension array allocated, computed and used in enerbil module exclusively

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: psold         !! Old surface dry static energy (J kg^{-1})
!$OMP THREADPRIVATE(psold)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: qsol_sat      !! Saturated specific humudity for old temperature (kg kg^{-1})
!$OMP THREADPRIVATE(qsol_sat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: pdqsold       !! Deriv. of saturated specific humidity at old temp
                                                                 !! (kg (kg s)^{-1})
!$OMP THREADPRIVATE(pdqsold)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: psnew         !! New surface static energy (J kg^{-1})
!$OMP THREADPRIVATE(psnew)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: qsol_sat_new  !! New saturated surface air moisture (kg kg^{-1})
!$OMP THREADPRIVATE(qsol_sat_new)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: netrad        !! Net radiation (W m^{-2})
!$OMP THREADPRIVATE(netrad)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: lwabs         !! LW radiation absorbed by the surface (W m^{-2})
!$OMP THREADPRIVATE(lwabs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: lwup          !! Long-wave up radiation (W m^{-2})
!$OMP THREADPRIVATE(lwup)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: lwnet         !! Net Long-wave radiation (W m^{-2})
!$OMP THREADPRIVATE(lwnet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: fluxsubli     !! Energy of sublimation (mm day^{-1})
!$OMP THREADPRIVATE(fluxsubli)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: qsat_air      !! Air saturated specific humidity (kg kg^{-1})
!$OMP THREADPRIVATE(qsat_air)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:) :: tair          !! Air temperature (K)
!$OMP THREADPRIVATE(tair)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: q_sol_pot               !! Potential surface humidity
!$OMP THREADPRIVATE(q_sol_pot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: temp_sol_pot            !! Potential surface temperature (unstressed)
!$OMP THREADPRIVATE(temp_sol_pot)
  
  CONTAINS  
    
   
  !!  ===========================================================================================================================
  !! SUBROUTINE		 		    : enerbil_main
  !!
  !!
  !>\BRIEF				    Calls each part of the energy budget calculation in sequence 
  !!
  !! DESCRIPTION			    : 
  !! This is the main routine for the 'enerbil' module. It is called
  !! once during the initialisation of ORCHIDEE, and then once for each time step. It is called a final
  !! time at the culmination of the last time step, for the writing of a restart file.\n
  !!
  !! The algorithm first calls 'enerbil_begin' for initialisation, followed by 'enerbil_surftemp' to
  !! calculate the surface static energy and the saturated surface humidity for the new time step.
  !! Next is the module 'enerbil_flux' which calculates the new surface temperature, the net radiation in
  !! the surface layer, the total evaporation and the latent and sensible heat fluxes. Finally comes
  !! 'enerbil_evapveg', which calculates the evaporation and transpiration from the vegetation.\n

  !! \n
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S)	: evapot, evapot_corr, temp_sol, qsurf, fluxsens, fluxlat, tsol_rad,
  !! vevapp, temp_sol_new, vevapnu, vevapsno, transpir, vevapwet, t2mdiag
  !!
  !! REFERENCE(S)		:
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ, Technical note, available (22/12/11):
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
  !! sur le cycle de l'eau, PhD Thesis, available (22/12/11):
  !! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf 
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers\n
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663.
  !!
  !! FLOWCHART			:
  !! \latexonly 
  !!     \includegraphics[scale=0.5]{enerbil_main_flowchart.png}
  !! \endlatexonly
  !! \n
  !_ ==============================================================================================================================
   
  SUBROUTINE enerbil_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, &
       & index, indexveg, zlev, lwdown, swnet, epot_air, temp_air, u, v, petAcoef, petBcoef, &
       & qair, peqAcoef, peqBcoef, pb, rau, vbeta, valpha, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, &
       & emis, soilflx, soilcap, q_cdrag, humrel, fluxsens, fluxlat, &
       & vevapp, transpir, transpot, vevapnu, vevapwet, vevapsno, vevapflo, t2mdiag, temp_sol, tsol_rad, &
       & temp_sol_new, qsurf, evapot, evapot_corr, rest_id, hist_id, hist2_id, &
       & precip_rain,snow,snowdz,snowrho,pgflux,snowflx,snowcap,temp_sol_add)
 

    !! 0 Variable and parameter description

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                           :: kjit             !! Time step number (-)
    INTEGER(i_std), INTENT(in)                           :: kjpindex         !! Domain size (-)
    INTEGER(i_std),INTENT (in)                           :: rest_id,hist_id  !! _Restart_ file and _history_ file identifier (-)
    INTEGER(i_std),INTENT (in)                           :: hist2_id         !! _history_ file 2 identifier (-)
    REAL(r_std), INTENT (in)                             :: dtradia          !! Time (s)
    LOGICAL, INTENT(in)                                  :: ldrestart_read   !! Logical for _restart_ file to read (-)
    LOGICAL, INTENT(in)                                  :: ldrestart_write  !! Logical for _restart_ file to write (-)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: index            !! Indeces of the points on the map (-)
    INTEGER(i_std),DIMENSION(kjpindex*nvm), INTENT(in)   :: indexveg         !! Indeces of the points on the 3D map
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: zlev             !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: lwdown           !! Down-welling long-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: swnet            !! Net surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: epot_air         !! Air potential energy (?? J)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_air         !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: u                !! Eastward Lowest level wind speed  (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: v                !! Northward Lowest level wind speed (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: petAcoef         !! PetAcoef (see note)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: petBcoef         !! PetBcoef (see note)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: qair             !! Lowest level specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: peqAcoef         !! PeqAcoef (see note)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: peqBcoef         !! PeqBcoef (see note)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: pb               !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: rau              !! Air density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: vbeta            !! Resistance coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: valpha           !! Resistance coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: vbeta1           !! Snow resistance (-) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: vbeta4           !! Bare soil resistance (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: vbeta5           !! Floodplains resistance
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: emis             !! Emissivity (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: soilflx          !! Soil flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: soilcap          !! Soil calorific capacity (J K^{-1])
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: q_cdrag          !! Surface drag (m s^{-1})

    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (in)   :: humrel           !! Soil moisture stress (within range 0 to 1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: vbeta2           !! Interception resistance (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: vbeta3           !! Vegetation resistance (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: vbeta3pot        !! Vegetation resistance for potential transpiration
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snow             !! Snow mass [kg/m2]
    REAL(r_std),DIMENSION (kjpindex),INTENT (in)         :: precip_rain      !! Rainfall
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)    :: snowrho          !! Snow density at each snow layer
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)    :: snowdz           !! Snow depth at each snow layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snowflx          !! Soil flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: snowcap          !! Soil flux (J m^{-2}^{-1})

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: vevapnu          !! Bare soil evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: vevapsno         !! Snow evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: vevapflo         !! Floodplains evaporation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: transpir         !! Transpiration (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: transpot         !! Potential transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)   :: vevapwet         !! Interception (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: t2mdiag          !! Diagnostic temperature at 2m above ground (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: temp_sol_new     !! New soil temperature (K)
    
    !! 0.3 Modified variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: evapot           !! Soil Potential Evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: evapot_corr      !! Soil Potential Evaporation Correction (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: temp_sol         !! Soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: qsurf            !! Surface specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: fluxsens         !! Sensible heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: fluxlat          !! Latent heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: tsol_rad         !! Tsol_rad (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: vevapp           !! Total of evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: pgflux           !! Net energy into snowpack(W/m^2)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)     :: temp_sol_add 


    !! 0.4 Local variables
    REAL(r_std),DIMENSION (kjpindex)                     :: epot_air_new, qair_new
    CHARACTER(LEN=80)                                    :: var_name         !! To store variables names for I/O
    INTEGER(i_std)                                       :: ji,ii,iv


!_ ================================================================================================================================
    
    !! 1. Run the initialisation routines 

    IF (l_first_enerbil) THEN

        IF (long_print) WRITE (numout,*) ' l_first_enerbil : call enerbil_init '
        CALL enerbil_init (kjit, ldrestart_read, kjpindex, index, rest_id, qair, temp_sol, temp_sol_new,  &
             &             qsurf, tsol_rad, vevapp, fluxsens, fluxlat, evapot, evapot_corr)

        CALL enerbil_var_init (kjpindex, temp_air, t2mdiag)

        RETURN 

    ENDIF

    
    !! 2. Prepare a restart file for the next simulation
    
    IF (ldrestart_write) THEN

        IF (long_print) WRITE (numout,*) ' we have to complete restart file with ENERBIL variables '

        var_name= 'temp_sol'  
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  temp_sol, 'scatter',  nbp_glo, index_g)

        var_name= 'qsurf'  
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  qsurf, 'scatter',  nbp_glo, index_g)

        var_name= 'evapot'  
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  evapot, 'scatter',  nbp_glo, index_g)

        var_name= 'evapot_corr'  
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  evapot_corr, 'scatter',  nbp_glo, index_g)

        var_name= 'tsolrad'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  tsol_rad, 'scatter',  nbp_glo, index_g)

        var_name= 'evapora'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  vevapp, 'scatter',  nbp_glo, index_g)

        var_name= 'fluxlat'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  fluxlat, 'scatter',  nbp_glo, index_g)

        var_name= 'fluxsens'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit,  fluxsens, 'scatter',  nbp_glo, index_g)

        var_name= 'tempsolpot'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, temp_sol_pot, 'scatter',  nbp_glo, index_g)

        var_name= 'qsolpot'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, q_sol_pot, 'scatter',  nbp_glo, index_g)


        RETURN 

    END IF

    
    !! 3. Computes some initialisation variables
 
    !  Computes some initialisation variables psold (the old surface static energy), qsol_sat
    ! (the saturated surface humidity) and pdqsold (the derivative of the saturated surface humidity,
    ! calculated with respect to the surface temperature at the 'old' timestep)
    CALL enerbil_begin (kjpindex, temp_sol, lwdown, swnet, pb, psold, qsol_sat, pdqsold, netrad, emis)


    !! 4. Computes psnew (the surface static energy at the 'new' timestep) 

    ! Computes psnew (the surface static energy at the 'new' timestep), qsol_sat_new (the surface
    ! saturated humidity at the 'new' timestep), temp_sol_new (the surface temperature at the 'new'
    ! timestep), qair_new (the lowest atmospheric humidity at the 'new' timestep) and epot_air_new
    ! (the lowest atmospheric evaporation potential at the 'new' timestep)
    CALL enerbil_surftemp (kjpindex, dtradia, zlev, emis, &
       & epot_air, petAcoef, petBcoef, qair, peqAcoef, peqBcoef, soilflx, rau, u, v, q_cdrag, vbeta,&
       & valpha, vbeta1, vbeta5, soilcap, lwdown, swnet, psnew, qsol_sat_new, temp_sol_new, &
       & qair_new, epot_air_new, snowdz, snowflx,snowcap)


    ! 2.1 computes surface temperature and humidity for a saturated surface
    !
    CALL enerbil_pottemp (kjpindex, dtradia, zlev, emis, &
       & epot_air, petAcoef, petBcoef, qair, peqAcoef, peqBcoef, soilflx, rau, u, v, q_cdrag, vbeta,&
       & valpha, vbeta1, vbeta5, soilcap, lwdown, swnet, q_sol_pot, temp_sol_pot)


    !! 5. Diagnoses components of the energy budget
    
    ! Diagnoses lwup (upwards long wave radiation), lwnet (net long wave radiation), tsol_rad (radiative
    ! ground temperature), netrad (net radiation), qsurf (surface humidity), vevapp (total evaporation),
    ! evapot (evaporation potential), evapot_corr (evaporation potential correction factor), 
    ! fluxlat (latent heat flux), fluxsubli (sublimination heat flux) and fluxsens (sensible heat flux).

    CALL enerbil_flux (kjpindex, dtradia, emis, temp_sol, rau, u, v, q_cdrag, vbeta, valpha, vbeta1, vbeta5, &
       & qair_new, epot_air_new, psnew, qsurf, &
       & fluxsens , fluxlat , fluxsubli, vevapp, temp_sol_new, lwdown, swnet, &
       & lwup, lwnet, pb, tsol_rad, netrad, evapot, evapot_corr, &
       & precip_rain,temp_air,pgflux,snowdz,snowcap,temp_sol_add)


    !! 6. Diagnoses the values for evaporation and transpiration 
    
    ! Diagnoses the values for evaporation and transpiration: vevapsno (snow evaporation), vevapnu 
    ! (bare soil evaporation), transpir (transpiration) and vevapwet (interception)
    CALL enerbil_evapveg (kjpindex, dtradia, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5,  &
       & rau, u, v, q_cdrag, qair_new, humrel, vevapsno, vevapnu , vevapflo, vevapwet, &
       & transpir, transpot, evapot, &
       & snowdz,snowrho,evapot_corr)

    !! 7. Assigns a temperature 2 metres above the ground 

    ! Assigns a temperature 2 metres above the ground. The output variable - t2mdiag, is used just for
    ! diagnostic purposes, for now.
    CALL enerbil_t2mdiag (kjpindex, temp_air, t2mdiag)

    CALL xios_orchidee_send_field("netrad",netrad)
    CALL xios_orchidee_send_field("evapot",evapot*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("evapot_alma",evapot/dt_sechiba)
    CALL xios_orchidee_send_field("evapot_corr",evapot_corr*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("evapot_corr_alma",evapot_corr/dt_sechiba)
    CALL xios_orchidee_send_field("lwdown",lwabs)
    CALL xios_orchidee_send_field("lwnet",lwnet)
    CALL xios_orchidee_send_field("Qv",fluxsubli)
    CALL xios_orchidee_send_field("PotSurfT",temp_sol_pot)
          

    IF ( .NOT. almaoutput ) THEN
       CALL histwrite_p(hist_id, 'netrad', kjit, netrad, kjpindex, index)
       CALL histwrite_p(hist_id, 'evapot', kjit, evapot, kjpindex, index)
       CALL histwrite_p(hist_id, 'evapot_corr', kjit, evapot_corr, kjpindex, index)
       CALL histwrite_p(hist_id, 'lwdown', kjit, lwabs,  kjpindex, index)
       CALL histwrite_p(hist_id, 'lwnet',  kjit, lwnet,  kjpindex, index)
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'netrad', kjit, netrad, kjpindex, index)
          CALL histwrite_p(hist2_id, 'evapot', kjit, evapot, kjpindex, index)
          CALL histwrite_p(hist2_id, 'evapot_corr', kjit, evapot_corr, kjpindex, index)
          CALL histwrite_p(hist2_id, 'lwdown', kjit, lwabs,  kjpindex, index)
          CALL histwrite_p(hist2_id, 'lwnet',  kjit, lwnet,  kjpindex, index)
       ENDIF
    ELSE
       CALL histwrite_p(hist_id, 'LWnet', kjit, lwnet, kjpindex, index)
       CALL histwrite_p(hist_id, 'Qv', kjit, fluxsubli, kjpindex, index)
       CALL histwrite_p(hist_id, 'PotEvap', kjit, evapot_corr, kjpindex, index)
       CALL histwrite_p(hist_id, 'PotEvapOld', kjit, evapot, kjpindex, index)
       CALL histwrite_p(hist_id, 'PotSurfT', kjit, temp_sol_pot, kjpindex, index)
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'LWnet', kjit, lwnet, kjpindex, index)
          CALL histwrite_p(hist2_id, 'Qv', kjit, fluxsubli, kjpindex, index)
          CALL histwrite_p(hist2_id, 'PotEvap', kjit, evapot_corr, kjpindex, index)
       ENDIF
    ENDIF

    IF (long_print) WRITE (numout,*) ' enerbil_main Done '

  END SUBROUTINE enerbil_main


  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : enerbil_init
  !!
  !>\BRIEF			            Performs dynamic allocation for local arrays, and retrieves
  !! previous data from restart file (if the switch is set to do so)
  !!
  !! DESCRIPTION			    : The purpose of this module is, firstly, to allocate space
  !! in memory for key variables within the 'enerbil' module. The second task is to retrieve previous data
  !! from the restart file, if the 'ldrestart-read' flag is set.
  !!
  !! RECENT CHANGE(S)			    : None
  !!
  !! MAIN OUTPUT VARIABLE(S)	            : temp_sol, temp_sol_new, qsurf, tsol_rad, vevapp, fluxsens,
  !! fluxlat, evapot, evapot_corr
  !!
  !! REFERENCE(S)			    : None
  !! 
  !! FLOWCHART                              : None
  !! \n
  !_ ==============================================================================================================================
  
  SUBROUTINE enerbil_init (kjit, ldrestart_read, kjpindex, index, rest_id, qair, temp_sol, temp_sol_new, &
       &                   qsurf, tsol_rad, vevapp, fluxsens, fluxlat, evapot, evapot_corr)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables
  
    INTEGER(i_std), INTENT (in)                              :: kjit               !! Time step number (-) 
    LOGICAL ,INTENT (in)                                     :: ldrestart_read     !! Logical for _restart_ file to read (-)
    INTEGER(i_std), INTENT (in)                              :: kjpindex           !! Domain size (-)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: index              !! Indeces of the points on the map (-)
    INTEGER(i_std), INTENT (in)                              :: rest_id            !! _Restart_ file identifier (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair               !! Lowest level specific humidity (kg kg^{-1})
    
    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: temp_sol           !! Ground temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: temp_sol_new       !! New ground temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsurf              !! Near surface specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: tsol_rad           !! Tsol_rad (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vevapp             !! Total of evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsens           !! Sensible heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxlat            !! Latent heat flux (W m^{-2})

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: evapot             !! Soil Potential Evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: evapot_corr        !! Soil Potential Evaporation Correction 
                                                                                   !! (mm day^{-1})

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ier
    CHARACTER(LEN=80)                                        :: var_name            !! To store variables names for I/O

!_ ================================================================================================================================

  !! 1 Initialisation of variables
 
    !! Initialisation of variables psold (old surface static energy), qsol_sat (saturated surface air moisture), pdqsold
    !! (derivative of satured specific humidity at the old temperature), psnew (new surface static energy), qsol_sat_new
    !! (new saturated surface air moisture), netrad (net radiation), lwabs (long wave absolute radiation), lwup (long wave
    !! (upwelling radiation), lwnet (long wave net radiation), fluxsubli (sublimination flux), qsat_air (saturated air humidity)
    !! and tair (air temperature).
    IF (l_first_enerbil) THEN 
        l_first_enerbil=.FALSE.
    ELSE 
        WRITE (numout,*) ' l_first_enerbil false . we stop '
        STOP 'enerbil_init'
    ENDIF

    ALLOCATE (psold(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in psold allocation. We stop.We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (qsol_sat(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in qsol_sat allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (pdqsold(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pdqsold allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (psnew(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in psnew allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (qsol_sat_new(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in qsol_sat_new allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (netrad(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in netrad allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (lwabs(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in lwabs allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (lwup(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in lwup allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (lwnet(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in lwnet allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (fluxsubli(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in fluxsubli allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (qsat_air(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in qsat_air allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (tair(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in tair allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF


    ALLOCATE (q_sol_pot(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in q_sol_pot allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF

    ALLOCATE (temp_sol_pot(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in temp_sol_pot allocation. We stop. We need kjpindex words = ',kjpindex
        STOP 'enerbil_init'
    END IF


  !! 2. Open restart input file and read data from the restart input file

    !! Open restart input file done by enerbil_init and read data from the restart input file
    !! that corresponds to processes within the 'enerbil' module (if the 'ldrestart_read' flag is set).\n
    !! The variables read are: temp_sol (surface temperature), qsurf (near surface specific humidity),
    !! evapot (soil potential evaporation), evapot_corr (corrected soil potential evaporation), tsolrad
    !! (radiative surface temperature), evapora (evaporation), fluxlat (latent heat flux), fluxsens
    !! (sensible heat flux) and temp_sol_new (new surface temperature).
    IF (ldrestart_read) THEN

        IF (long_print) WRITE (numout,*) ' we have to read a restart file for ENERBIL variables'

        var_name='temp_sol'
        CALL ioconf_setatt_p('UNITS', 'K')
        CALL ioconf_setatt_p('LONG_NAME','Surface temperature')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp_sol, "gather", nbp_glo, index_g)

        IF ( MINVAL(temp_sol) < 0. .OR. MAXVAL(temp_sol) > 500.) THEN
           WRITE(numout, *) "Big problem with surface temperature at restart : ",  MINVAL(temp_sol), MAXVAL(temp_sol)
           WHERE ( temp_sol < 0.) temp_sol = 280.
        ENDIF

        !
        !Config Key   = ENERBIL_TSURF
        !Config Desc  = Initial temperature if not found in restart
        !Config If    = OK_SECHIBA
        !Config Def   = 280.
        !Config Help  = The initial value of surface temperature if its value is not found
        !Config         in the restart file. This should only be used if the model is 
        !Config         started without a restart file.
        !Config Units = Kelvin [K]
        CALL setvar_p (temp_sol, val_exp,'ENERBIL_TSURF', 280._r_std)
        
        var_name= 'qsurf'
        CALL ioconf_setatt_p('UNITS', 'g/g')
        CALL ioconf_setatt_p('LONG_NAME','near surface specific humidity')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., qsurf, "gather", nbp_glo, index_g)
        IF ( ALL( qsurf(:) .EQ. val_exp ) ) THEN 
           qsurf(:) = qair(:)
        ENDIF
        
        var_name= 'evapot'
        CALL ioconf_setatt_p('UNITS', 'mm day^{-1}')
        CALL ioconf_setatt_p('LONG_NAME','Soil Potential Evaporation')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., evapot, "gather", nbp_glo, index_g)
        CALL setvar_p (evapot, val_exp, 'ENERBIL_EVAPOT', zero)
        
        var_name= 'evapot_corr'
        CALL ioconf_setatt_p('UNITS', 'mm day^{-1}')
        CALL ioconf_setatt_p('LONG_NAME','Corrected Soil Potential Evaporation')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., evapot_corr, "gather", nbp_glo, index_g)
        
        !Config Key   = ENERBIL_EVAPOT
        !Config Desc  = Initial Soil Potential Evaporation
        !Config If    = OK_SECHIBA       
        !Config Def   = 0.0
        !Config Help  = The initial value of soil potential evaporation if its value 
        !Config         is not found in the restart file. This should only be used if
        !Config         the model is started without a restart file. 
        !Config Units = 
        CALL setvar_p (evapot_corr, val_exp, 'ENERBIL_EVAPOT', zero)
                
        var_name= 'tsolrad'
        CALL ioconf_setatt_p('UNITS', 'K')
        CALL ioconf_setatt_p('LONG_NAME','Radiative surface temperature')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., tsol_rad, "gather", nbp_glo, index_g)
        IF ( ALL( tsol_rad(:) .EQ. val_exp ) ) THEN 
           tsol_rad(:) = temp_sol(:)
        ENDIF
        
        !! 1.3 Set the fluxes so that we have something reasonable and not NaN on some machines
        var_name= 'evapora'
        CALL ioconf_setatt_p('UNITS', 'Kg/m^2/dt')
        CALL ioconf_setatt_p('LONG_NAME','Evaporation')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., vevapp, "gather", nbp_glo, index_g)
        IF ( ALL( vevapp(:) .EQ. val_exp ) ) THEN 
           vevapp(:) = zero
        ENDIF
        
        var_name= 'fluxlat'
        CALL ioconf_setatt_p('UNITS', 'W/m^2')
        CALL ioconf_setatt_p('LONG_NAME','Latent heat flux')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., fluxlat, "gather", nbp_glo, index_g)
        IF ( ALL( fluxlat(:) .EQ. val_exp ) ) THEN 
           fluxlat(:) = zero
        ENDIF
        
        var_name= 'fluxsens'
        CALL ioconf_setatt_p('UNITS', 'W/m^2')
        CALL ioconf_setatt_p('LONG_NAME','Sensible heat flux')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., fluxsens, "gather", nbp_glo, index_g)
        IF ( ALL( fluxsens(:) .EQ. val_exp ) ) THEN 
           fluxsens(:) = zero
        ENDIF
        
        var_name= 'tempsolpot'
        CALL ioconf_setatt_p('UNITS', 'K')
        CALL ioconf_setatt_p('LONG_NAME','Potential surface temperature')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp_sol_pot, "gather", nbp_glo, index_g)
        IF ( ALL( temp_sol_pot(:) .EQ. val_exp ) ) THEN 
           temp_sol_pot = temp_sol
        ENDIF
        !
        var_name= 'qsolpot'
        CALL ioconf_setatt_p('UNITS', 'kg/m^2')
        CALL ioconf_setatt_p('LONG_NAME','Potential saturated surface humidity')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., q_sol_pot, "gather", nbp_glo, index_g)
        IF ( ALL( q_sol_pot(:) .EQ. val_exp ) ) THEN 
           q_sol_pot = qsurf
        ENDIF


        ! initialises temp_sol_new saved in thermosoil
        var_name='temp_sol_new'
        CALL ioconf_setatt_p('UNITS', 'K')
        CALL ioconf_setatt_p('LONG_NAME','New Surface temperature')
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp_sol_new, "gather", nbp_glo, index_g)
        IF ( ALL( temp_sol_new(:) .EQ. val_exp ) ) THEN
           ! Variables not found in restart file.
           temp_sol_new(:) = temp_sol(:)
        ENDIF
        
    ELSE

       ! initialises temp_sol_new 
       temp_sol_new(:) = temp_sol(:)

    ENDIF

    iF (long_print) WRITE (numout,*) ' enerbil_init done '

  END SUBROUTINE enerbil_init


  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : enerbil_clear
  !!
  !>\BRIEF				    Routine deallocates clear output variables if already allocated.
  !!
  !! DESCRIPTION			    : This is a 'housekeeping' routine that deallocates the key output
  !! variables, if they have already been allocated. The variables that are deallocated are psold,
  !! qsol_sat, pdqsold, psnew, qsol_sat_new, netrad, lwabs, lwup, lwnet, fluxsubli, qsat_air, tair
  !!
  !! RECENT CHANGE(S)			    : None
  !!
  !! MAIN OUTPUT VARIABLE(S)	            : None
  !!
  !! REFERENCES				    : None
  !! 
  !! FLOWCHART                              : None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE enerbil_clear ()
    l_first_enerbil=.TRUE.
    IF ( ALLOCATED (psold)) DEALLOCATE (psold)
    IF ( ALLOCATED (qsol_sat)) DEALLOCATE (qsol_sat)
    IF ( ALLOCATED (pdqsold)) DEALLOCATE (pdqsold)
    IF ( ALLOCATED (psnew)) DEALLOCATE (psnew)
    IF ( ALLOCATED (qsol_sat_new)) DEALLOCATE (qsol_sat_new)
    IF ( ALLOCATED (netrad)) DEALLOCATE (netrad)
    IF ( ALLOCATED (lwabs)) DEALLOCATE (lwabs)
    IF ( ALLOCATED (lwup)) DEALLOCATE (lwup)
    IF ( ALLOCATED (lwnet)) DEALLOCATE (lwnet)
    IF ( ALLOCATED (fluxsubli)) DEALLOCATE (fluxsubli)
    IF ( ALLOCATED (qsat_air)) DEALLOCATE (qsat_air)
    IF ( ALLOCATED (tair)) DEALLOCATE (tair)
    IF ( ALLOCATED (q_sol_pot)) DEALLOCATE (q_sol_pot)
    IF ( ALLOCATED (temp_sol_pot)) DEALLOCATE (temp_sol_pot)
   
   ! ??open restart input file done by enerbil_init
   ! and read data from restart input file for ENERBIL process??
 
  END SUBROUTINE enerbil_clear


  !!  =============================================================================================================================
  !! SUBROUTINE		 			: enerbil_var_init
  !!
  !>\BRIEF					Housekeeping routine that calls enerbil_t2mdiag.
  !!
  !! DESCRIPTION				: Housekeeping routine that calls the 2 metre temperature
  !! diagnostic subroutine 'enerbil_2tmdiag'.
  !!
  !! MAIN OUTPUT VARIABLE(S)	                : t2mdiag
  !!
  !! REFERENCE(S)				: None
  !!
  !! FLOWCHART                                  : None
  !! \n
  !_ ==============================================================================================================================
 
  SUBROUTINE enerbil_var_init (kjpindex, temp_air, t2mdiag)

    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                         :: kjpindex       !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_air       !! Air temperature (K)

    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Output fields
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: t2mdiag        !! Temperature at 2 metres above ground(K)
!_ ================================================================================================================================

    CALL enerbil_t2mdiag (kjpindex, temp_air, t2mdiag)

    IF (long_print) WRITE (numout,*) ' enerbil_var_init done '

  END SUBROUTINE enerbil_var_init


  !!  =============================================================================================================================
  !! SUBROUTINE		 			: enerbil_begin
  !!
  !>\BRIEF					Preliminary variables required for the calculation of
  !! the energy budget are derived.
  !!
  !! DESCRIPTION				: This routines computes preliminary variables required
  !! for the calculation of the energy budget: the old surface static energy (psold), the surface saturation
  !! humidity (qsol_sat), the derivative of satured specific humidity at the old temperature (pdqsold) 
  !! and the net radiation (netrad).
  !!
  !! MAIN OUTPUT VARIABLE(S)	                : psold, qsol_sat, pdqsold, netrad
  !!
  !! REFERENCE(S)				:
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ, Technical note, available (22/12/11):
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers\n
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663.
  !!
  !! FLOWCHART  : None                     
  !! \n
  !_ ==============================================================================================================================
  
  SUBROUTINE enerbil_begin (kjpindex, temp_sol, lwdown, swnet, pb, psold, qsol_sat, pdqsold, netrad, emis)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: temp_sol         !! Soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: lwdown           !! Down-welling long-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: swnet            !! Net surface short-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: pb               !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      :: emis             !! Emissivity (-)

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: psold            !! Old surface dry static energy (J kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: qsol_sat	   !! Saturated specific humudity for old temperature 
                                                                           !! (kg kg^{-1})    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: pdqsold	   !! Derivative of satured specific humidity at the old 
                                                                           !! temperature (kg (kg s)^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)     :: netrad           !! Net radiation (W m^{-2})

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                     :: ji
    REAL(r_std), DIMENSION(kjpindex)                   :: dev_qsol
    REAL(r_std), PARAMETER                             :: missing = 999998.


!_ ================================================================================================================================
  
  !! 1. Computes psold (the surface static energy for the old timestep)
   
    !! We here define the surface static energy for the 'old' timestep, in terms of the surface
    !! temperature and heat capacity.
    !! \latexonly 
    !!     \input{enerbilbegin1.tex}
    !! \endlatexonly
    psold(:) = temp_sol(:)*cp_air
  !! 2. Computes qsol_sat (the surface saturated humidity).
 
    !! We call the routine 'qsatcalc' from within the module 'src_parameters/constantes_veg'.
    CALL qsatcalc (kjpindex, temp_sol, pb, qsol_sat)
    IF ( diag_qsat ) THEN
      IF ( ANY(ABS(qsol_sat(:)) .GT. missing) ) THEN
        DO ji = 1, kjpindex
          IF ( ABS(qsol_sat(ji)) .GT. missing) THEN
            WRITE(numout,*) 'ERROR on ji = ', ji
            WRITE(numout,*) 'temp_sol(ji),  pb(ji) :', temp_sol(ji),  pb(ji)
            CALL ipslerr_p (3,'enerbil_begin', &
 &           'qsol too high ','','')
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    
  !! 3. Computes pdqsold 
    
    !! Computes pdqsold (the derivative of the saturated humidity with respect to temperature
    !! analysed at the surface temperature at the 'old' timestep.
    !! We call the routine 'dev_qsatcalc' from qsat_moisture module.
    CALL dev_qsatcalc (kjpindex, temp_sol, pb, dev_qsol)
    
    !! \latexonly 
    !!     \input{enerbilbegin2.tex}
    !! \endlatexonly
    pdqsold(:) = dev_qsol(:)
    

    IF ( diag_qsat ) THEN
      IF ( ANY(ABS( pdqsold(:)) .GT. missing) ) THEN
        DO ji = 1, kjpindex
          IF ( ABS( pdqsold(ji)) .GT. missing ) THEN
            WRITE(numout,*) 'ERROR on ji = ', ji
            WRITE(numout,*) 'temp_sol(ji),  pb(ji) :', temp_sol(ji),  pb(ji)
            CALL ipslerr_p (3,'enerbil_begin', &
 &           'pdqsold too high ','','')
          ENDIF
        ENDDO
      ENDIF
    ENDIF

  !! 4. Computes the net radiation and the absorbed LW radiation absorbed at the surface. 

    !! Long wave radiation absorbed by the surface is the product of emissivity and downwelling LW radiation    
    !! \latexonly 
    !!     \input{enerbilbegin3.tex}
    !! \endlatexonly
    lwabs(:) = emis(:) * lwdown(:)

    !! Net radiation is calculated as:
    !! \latexonly 
    !!     \input{enerbilbegin4.tex}
    !! \endlatexonly
    netrad(:) = lwdown(:) + swnet (:) - (emis(:) * c_stefan * temp_sol(:)**4 + (un - emis(:)) * lwdown(:)) 

    IF (long_print) WRITE (numout,*) ' enerbil_begin done '

  END SUBROUTINE enerbil_begin


  !!  =============================================================================================================================
  !! SUBROUTINE		 			: enerbil_surftemp
  !!
  !>\BRIEF					This routine computes the new surface static energy
  !! (psnew) and the saturated humidity at the surface (qsol_sat_new). 
  !!
  !! DESCRIPTION				: This is the key part of the enerbil module, for which
  !! the energy budget in the surface layer is solved and changes over the model timestep 'dtradia' are
  !! quantified for the surface static energy, surface temperature, surface humidity, the 'air' (or lowest
  !! level atmospheric model) temperature, the 'air' (or lowest level atmospheric model) humidity, 
  !! according to the method that is laid out by Dufresne \& Ghattas (2009) and Shultz et al. (2001).
  !!
  !! It computes the energy balance at the surface with an implicit scheme, that is connected to the 
  !! Richtmyer and Morton algorithm of the PBL. By computing the surface temperature and surface humidity 
  !! the routine also implicitly estimates the various fluxes which balance in order to give the new 
  !! temperature. Thus once the surface temperature has been obtained all the fluxes need to be diagnosed. 
  !!
  !! If ok_explicitsnow is used, the notion of skin layer in ORCHIDEE is abandoned and
  !! the first thin snow layer is adopted to solve the surface energy budget.
  !! Implicit method is still used for coupling to the atmosphere.
  !! In this new scheme, the snow temperature profile is simulatenously solved
  !! along with the surface snow temperature, and the details are referenced to Boone et al. (2010)
  !!
  !! MAIN OUTPUT VARIABLE(S)	: psnew, qsol_sat_new, temp_sol_new, qair_new, epot_air_new
  !!
  !! REFERENCE(S)					:
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ, Technical note, available (22/12/11):
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
  !! sur le cycle de l'eau, PhD Thesis, available (22/12/11):
  !! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pd 
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663.
  !!
  !! FLOWCHART    : None
  !!
  !_ ==============================================================================================================================  
  

  SUBROUTINE enerbil_surftemp (kjpindex, dtradia, zlev, emis, epot_air, &
     & petAcoef, petBcoef, qair, peqAcoef, peqBcoef, soilflx, rau, u, v, q_cdrag, vbeta,&
     & valpha, vbeta1, vbeta5, soilcap, lwdown, swnet, psnew, qsol_sat_new, temp_sol_new, &
     & qair_new, epot_air_new, &
     & snowdz,snowflx,snowcap)


    !! 0. Variable and parameter declaration 

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex      !! Domain size (-)
    REAL(r_std), INTENT(in)                                  :: dtradia       !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: zlev          !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: emis          !! Emissivity (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: epot_air      !! Air potential energy (?? J)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petAcoef      !! PetAcoef (see note)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petBcoef      !! PetBcoef (see note)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair          !! Lowest level specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqAcoef      !! PeqAcoef (see note)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqBcoef      !! PeqBcoef (see note)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: soilflx       !! Soil flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau           !! Air density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u, v          !! Wind velocity by directional components 
                                                                              !! u (Eastwards) and v (Northwards) (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag       !! Surface drag (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta         !! Resistance coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: valpha        !! Resistance coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta1        !! Snow resistance (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta5        !! Floodplains resistance
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: soilcap       !! Soil calorific capacity (J K^{-1])
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown        !! Down-welling long-wave flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet         !! Net surface short-wave flux (W m^{-2})
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)      :: snowdz        !! Snow depth
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: snowflx       !! Soil flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: snowcap       !! Soil flux (J m^{-2}^{-1})


    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: psnew         !! New surface static energy (J kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsol_sat_new  !! New saturated surface air moisture (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: temp_sol_new  !! New soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qair_new      !! New air moisture (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: epot_air_new  !! New air temperature (K)

    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                   :: ji
    REAL(r_std),DIMENSION (kjpindex) :: zicp
    REAL(r_std)                      :: fevap
    REAL(r_std)                      :: sensfl_old, larsub_old, lareva_old, dtheta, sum_old, sum_sns
    REAL(r_std)                      :: zikt, zikq, netrad_sns, sensfl_sns, larsub_sns, lareva_sns
    REAL(r_std)                      :: speed

!_ ================================================================================================================================
 
    zicp = un / cp_air
    
    DO ji=1,kjpindex

    !! 1. Derivation of auxiliary variables
      
      !! \latexonly 
      !!     \input{enerbilsurftemp1.tex}
      !! \endlatexonly
      speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
      !
      !! \latexonly 
      !!     \input{enerbilsurftemp2.tex}
      !! \endlatexonly
      zikt = 1/(rau(ji) * speed * q_cdrag(ji))
      zikq = 1/(rau(ji) * speed * q_cdrag(ji))
      
    !! 2. Computation of fluxes for the old surface conditions
     
      !! As a reference, we first calculate the sensible and latent heat for the 'old' timestep.
 
      !! 2.1 Sensible heat for the 'old' timestep
      !! This is equation (64) of (Dufresne & Ghattas, 2009), rewritten in terms of H_old. We make the 
      !! approximation that (P_r/P)^{kappa}=1 and convert the A and B coefficients to the ORCHIDEE
      !! format (see introductory note, above). Also, the equation here is in terms of surface static
      !! energy, which is defined as per section 3.1, above.
      !! \latexonly 
      !!     \input{enerbilsurftemp3.tex}
      !! \endlatexonly
      sensfl_old = (petBcoef(ji) -  psold(ji)) / (zikt -  petAcoef(ji))
      !! 2.2 Latent heat for the 'old' timestep
      !! This is equation (70) of (Dufresne & Ghattas, 2009), rewritten in terms of {\lambda E}_{old}.
      !! Again we convert the A and B coefficients from the LMDZ format to the ORCHIDEE format.\n
      !! There are two forms of the equation - first for the latent heat as a result of sublimination
      !! processes and then as a result of evaporative processes:
      !! \latexonly 
      !!     \input{enerbilsurftemp4.tex}
      !! \endlatexonly
      larsub_old = chalsu0 * vbeta1(ji) * (un - vbeta5(ji)) * (peqBcoef(ji) -  qsol_sat(ji)) / (zikq - peqAcoef(ji))
      !! \latexonly 
      !!     \input{enerbilsurftemp5.tex}
      !! \endlatexonly
      lareva_old = chalev0 * (un - vbeta1(ji)) * (un - vbeta5(ji)) * vbeta(ji) * &
           & (peqBcoef(ji) -  valpha(ji) * qsol_sat(ji)) / (zikq - peqAcoef(ji)) &
           & + chalev0 * vbeta5(ji) * (peqBcoef(ji) -  qsol_sat(ji)) / (zikq - peqAcoef(ji))
      
    !! 3. Calculation of sensitivity terms
     
      !! We here calculate the sensitivity for the different fluxes to changes in dtheta, which is the
      !! change in the surface static energy over the model time step (dtradia).
      
      !! 3.1 Net radiation sensitivity
      !! This is an explicit-implicit representation of the Stefan-Boltzmann law - the explicit terms 
      !! are ${ps}_{old}^3$, and it is completed when multiplied by dtheta (defined above).
      !! \latexonly 
      !!     \input{enerbilsurftemp6.tex}
      !! \endlatexonly
      netrad_sns = zicp(ji) * quatre * emis(ji) * c_stefan * ((zicp(ji) * psold(ji))**3)
      
      !! 3.2 Sensible heat flux sensitivity
      !! This is the temperature sensitivity term N_1^h derived in equation (66) of (Dufresne & Ghattas, 2009),
      !! where we again assume that (P_r/P)^{kappa}=1, convert the A and B coefficients to the ORCHIDEE format,
      !! and rewrite in terms of surface static energy, rather than temperature (see section 3.1, above).
      !! \latexonly 
      !!     \input{enerbilsurftemp7.tex}
      !! \endlatexonly
      sensfl_sns = un / (zikt -  petAcoef(ji))
      !! 3.3 Latent heat flux sensitivity 
      !! This is the humidity sensitivity term N_1^q derived in equation (72) of (Dufresne & Ghattas, 2009), where
      !! the A and B coefficients are written in the ORCHIDEE format.
      !! larsub_sns is the latent heat sensitivity for sublimination. The coefficient vbeta1 is the evaporation
      !! coefficient. It represents the relationship between the real and potential evaporation. It is derived
      !! in the module 'src_sechiba/diffuco_snow'.
      !! \latexonly 
      !!     \input{enerbilsurftemp8.tex}
      !! \endlatexonly
      larsub_sns = chalsu0 * vbeta1(ji) * (un - vbeta5(ji)) * zicp(ji) * pdqsold(ji) / (zikq - peqAcoef(ji))

      !! lareva_sns is the latent heat sensitivity for evaporation. vbeta1 (src_sechiba/diffuco_snow), 
      !! vbeta (src_sechiba/diffuco_comb) and valpha (src_sechiba/diffuco_comb) are the evaporation
      !! coefficients.
      !! \latexonly 
      !!     \input{enerbilsurftemp9.tex}
      !! \endlatexonly
      lareva_sns = chalev0 * ((un - vbeta1(ji))*(un - vbeta5(ji)) * vbeta(ji) * valpha(ji) + vbeta5(ji)) * &
           & zicp(ji) * pdqsold(ji) / (zikq - peqAcoef(ji))

    !! 4. Solution of the energy balance
      
      !! 4.1 We calculate the total flux for the 'old' timestep.
      !! \latexonly 
      !!     \input{enerbilsurftemp10.tex}
      !! \endlatexonly
      IF (SUM(snowdz(ji,:)) .GT. 0.01 .AND. ok_explicitsnow) THEN
         sum_old = netrad(ji) + sensfl_old + larsub_old + lareva_old + snowflx(ji)
      ELSE
         sum_old = netrad(ji) + sensfl_old + larsub_old + lareva_old + soilflx(ji)
      ENDIF

      !! 4.2 We calculate the total sensitivity dtheta (the change in the
      !! surface static energy over the timestep.
      !! \latexonly 
      !!     \input{enerbilsurftemp11.tex}
      !! \endlatexonly
      sum_sns = netrad_sns + sensfl_sns + larsub_sns + lareva_sns
      !! 4.3 Calculation of dtheta (change in surface static energy over the
      !! timestep.
      !! \latexonly 
      !!     \input{enerbilsurftemp12.tex}
      !! \endlatexonly
      IF (SUM(snowdz(ji,:)) .GT. 0.01 .AND. ok_explicitsnow) THEN
         dtheta = dtradia * sum_old / (zicp(ji) * snowcap(ji) + dtradia * sum_sns)
      ELSE
         dtheta = dtradia * sum_old / (zicp(ji) * soilcap(ji) + dtradia * sum_sns)
      ENDIF

      !! 4.4 Determination of state variables for the 'new' timestep
      !! No we have dtheta, we can determine the surface static energy that
      !! corresponds to the 'new' timestep.
      !! \latexonly 
      !!     \input{enerbilsurftemp13.tex}
      !! \endlatexonly
      psnew(ji) =  psold(ji) + dtheta

      
      !! The new surface saturated humidity can be calculated by equation (69)
      !! of (Dufresne & Ghattas, 2009), in which we substitute dtheta for the
      !! change between old and new temperature using the relationship from 3.1,
      !! above.
      !! \latexonly 
      !!     \input{enerbilsurftemp14.tex}
      !! \endlatexonly
      qsol_sat_new(ji) = qsol_sat(ji) + zicp(ji) * pdqsold(ji) * dtheta
      !! The new surface temperature is determined from the new surface static temperature,
      !! again by using the relationship in section 3.1.
      !! \latexonly 
      !!     \input{enerbilsurftemp15.tex}
      !! \endlatexonly

      temp_sol_new(ji) = psnew(ji) / cp_air

      !! 4.5 Calculation of new evaporation potential and new evaporation latent heat
      !! flux (???)
      !! \latexonly 
      !!     \input{enerbilsurftemp16.tex}
      !! \endlatexonly
      epot_air_new(ji) = zikt * (sensfl_old - sensfl_sns * dtheta) + psnew(ji)
      
      !! \latexonly 
      !!     \input{enerbilsurftemp17.tex}
      !! \endlatexonly
      fevap = (lareva_old - lareva_sns * dtheta) + (larsub_old - larsub_sns * dtheta)

      IF ( ABS(fevap) < EPSILON(un) ) THEN
        
        !! \latexonly 
        !!     \input{enerbilsurftemp18.tex}
        !! \endlatexonly
        qair_new(ji) = qair(ji)
      ELSE
        !! \latexonly 
        !!     \input{enerbilsurftemp19.tex}
        !! \endlatexonly     
        qair_new(ji) = zikq * un / ( chalsu0 *  vbeta1(ji) * (un - vbeta5(ji)) + &
           & chalev0 * ((un - vbeta1(ji))*(un - vbeta5(ji)) * vbeta(ji) * valpha(ji) + vbeta5(ji)) ) &
           & * fevap + qsol_sat_new(ji)
      ENDIF
    ENDDO

    !write(2,*) 'tao debug in enerbil_surftemp,temp_sol,snowflx,soilflx',temp_sol_new(:),snowflx(:),soilflx(:)
    !write(2,*) 'tao debug in enerbil_surftemp,snowdz',sum(snowdz,2)
    IF (long_print) WRITE (numout,*) ' enerbil_surftemp done '

  END SUBROUTINE enerbil_surftemp


 !! =============================================================================================================================
 !! SUBROUTINE		 			: enerbil_pottemp
 !!
 !>\BRIEF		This subroutine computes the surface temperature and humidity should the surface been unstressed	
 !!
 !! DESCRIPTION				: 
 !!
 !! MAIN OUTPUT VARIABLE(S)	: 
 !!
 !! REFERENCE(S)		:
 !!
 !! FLOWCHART    : None
 !!
 !_ ==============================================================================================================================  

  SUBROUTINE enerbil_pottemp (kjpindex, dtradia, zlev, emis, epot_air, &
     & petAcoef, petBcoef, qair, peqAcoef, peqBcoef, soilflx, rau, u, v, q_cdrag, vbeta,&
     & valpha, vbeta1, vbeta5, soilcap, lwdown, swnet, q_sol_pot, temp_sol_pot) 

    ! interface 
    ! input scalar 
    INTEGER(i_std), INTENT(in)                               :: kjpindex      !! Domain size
    REAL(r_std), INTENT(in)                                  :: dtradia       !! Time step in seconds
    ! input fields
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: zlev          !! Height of first layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: emis          !! Emissivity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: epot_air      !! Air potential energy
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petAcoef      !! PetAcoef
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petBcoef      !! PetBcoef
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair          !! Lowest level specific humidity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqAcoef      !! PeqAcoef
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqBcoef      !! PeqBcoef
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: soilflx       !! Soil flux
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau           !! Density
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u, v          !! Wind
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag       !!
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta         !! Resistance coefficient
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: valpha        !! Resistance coefficient
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta1        !! Snow resistance
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta5        !! Floodplains resistance
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: soilcap       !! Soil calorific capacity
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown        !! Down-welling long-wave flux
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet         !! Net surface short-wave flux
    ! output fields
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: q_sol_pot     !! Potential surface air moisture
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: temp_sol_pot  !! Potential soil temperature


    ! local declaration
    INTEGER(i_std)                  :: ji
    REAL(r_std),DIMENSION (kjpindex) :: zicp
    REAL(r_std)                      :: fevap
    REAL(r_std)                      :: sensfl_old, larsub_old, lareva_old, dtheta, sum_old, sum_sns
    REAL(r_std)                      :: zikt, zikq, netrad_sns, sensfl_sns, larsub_sns, lareva_sns
    REAL(r_std)                      :: speed

!_ ==============================================================================================================================

    zicp = un / cp_air
    !
    DO ji=1,kjpindex

       dtheta = zero
       fevap = zero

       temp_sol_pot(ji) = temp_sol_pot(ji) + dtheta

       q_sol_pot(ji) = q_sol_pot(ji) + fevap

    ENDDO

  END SUBROUTINE enerbil_pottemp


  !!  =============================================================================================================================
  !! SUBROUTINE		 			: enerbil_flux
  !!
  !>\BRIEF					Computes the new soil temperature, net radiation and
  !! latent and sensible heat flux for the new time step.
  !!
  !! DESCRIPTION				: This routine diagnoses, based on the new soil temperature, 
  !! the net radiation, the total evaporation, the latent heat flux, the sensible heat flux and the 
  !! sublimination flux. It also diagnoses the potential evaporation used for the fluxes (evapot) and the potential
  !! as defined by Penman & Monteith (Monteith, 1965) based on the correction term developed by Chris 
  !! Milly (1992). This Penman-Monteith formulation is required for the estimation of bare soil evaporation 
  !! when the 11 layer CWRR moisture scheme is used for the soil hydrology.
  !!
  !! MAIN OUTPUT VARIABLE(S)	: qsurf, fluxsens, fluxlat, fluxsubli, vevapp, lwup, lwnet,
  !! tsol_rad, netrad, evapot, evapot_corr
  !!
  !! REFERENCE(S)					:
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ, Technical note, available (22/12/11):
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
  !! sur le cycle de l'eau, PhD Thesis, available (22/12/11):
  !! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
  !! - Monteith, JL, 1965. Evaporation and Environment, paper presented at Symposium of the Society
  !! for Experimental Biology
  !! - Monteith & Unsworth, 2008. Principles of Environmental Physics (third edition), published Elsevier
  !! ISBN 978-0-12-505103-3
  !! - Milly, P. C. D., 1992: Potential Evaporation and Soil Moisture in General Circulation Models. 
  !! Journal of Climate, 5, pp. 209–226.
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663.
  !!
  !! FLOWCHART					: None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE enerbil_flux (kjpindex, dtradia, emis, temp_sol, rau, u, v, q_cdrag, vbeta, valpha, vbeta1, vbeta5, &
       & qair, epot_air, psnew, qsurf, fluxsens, fluxlat, fluxsubli, vevapp, temp_sol_new, &
       & lwdown, swnet, lwup, lwnet, pb, tsol_rad, netrad, evapot, evapot_corr,&
       & precip_rain,temp_air,pgflux,snowdz,snowcap,temp_sol_add)

    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex      !! Domain size (-)
    REAL(r_std), INTENT(in)                                  :: dtradia       !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: emis          !! Emissivity (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol      !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: rau           !! Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u,v           !! Wind velocity in components u and v (m s^{-1}) 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q_cdrag       !! Surface drag (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta         !! Resistance coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: valpha        !! Resistance coefficient (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta1        !! Snow resistance  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: vbeta5        !! Flood resistance 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair          !! Lowest level specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: epot_air      !! Air potential energy (J)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: psnew         !! New surface static energy (J kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new  !! New soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb            !! Lowest level pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown        !! Downward long wave radiation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet         !! Net short wave radiation (W m^{-2})
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)      :: snowdz        !! Snow depth
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)             :: precip_rain   !! Rainfall
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)             :: temp_air      !! Air temperature
    REAL(r_std), DIMENSION (kjpindex),INTENT(in)             :: snowcap
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsurf         !! Surface specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsens      !! Sensible heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxlat       !! Latent heat flux (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsubli     !! Energy of sublimation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vevapp        !! Total of evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: lwup          !! Long-wave up radiation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: lwnet         !! Long-wave net radiation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: tsol_rad      !! Radiative surface temperature (W m^{-2})

    
    !! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: netrad        !! Net radiation (W m^{-2})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: evapot        !! Soil Potential Evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: evapot_corr   !! Soil Potential Evaporation Correction 
                                                                              !! (?? mm day^{-1})
    REAL(r_std), DIMENSION (kjpindex),INTENT(inout)          :: temp_sol_add
    REAL(r_std), DIMENSION (kjpindex),INTENT(inout)          :: pgflux        !! net energy into snowpack
    !! 0.4 Local variables

    INTEGER(i_std)                                               :: ji
    REAL(r_std),DIMENSION (kjpindex)        		         :: grad_qsat
    REAL(r_std)                               		         :: correction
    REAL(r_std)                                			 :: speed              !! Speed (m s^{-1}), 
    REAL(r_std)                                			 :: qc          !! Surface drag coefficient (??)
    LOGICAL,DIMENSION (kjpindex)            		         :: warning_correction
    REAL(r_std),DIMENSION (kjpindex)                             :: fluxsens_tmp
    REAL(r_std),DIMENSION (kjpindex)                             :: fluxlat_tmp
    REAL(r_std),DIMENSION (kjpindex)                             :: netrad_tmp,lwup_tmp,tsol_rad_tmp
    REAL(r_std), DIMENSION (kjpindex)                            :: zgflux
    REAL(r_std)                                                  :: allowed_err
    REAL(r_std), DIMENSION (kjpindex)                            :: qsol_sat_tmp
    REAL(r_std), DIMENSION (kjpindex)                            :: zerocelcius
    REAL(r_std), DIMENSION (kjpindex)                            :: PHPSNOW
!_ ================================================================================================================================
    
    zerocelcius(:) = tp_00
    CALL qsatcalc (kjpindex, zerocelcius, pb, qsol_sat_tmp)

    DO ji=1,kjpindex
 
    !! 1. Determination of 'housekeeping' variables
      
      !! The horizontal wind speed is calculated from the velocities along each axis using Pythagorus' theorem. 
      !! \latexonly 
      !!     \input{enerbilflux1.tex}
      !! \endlatexonly
      speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
      
      !! From definition, the surface drag coefficient is defined:
      !! \latexonly 
      !!     \input{enerbilflux2.tex}
      !! \endlatexonly
      qc = speed * q_cdrag(ji)
    
    !! 2. Calculation of the net upward flux of longwave radiation
      
      !! We first of all calculate the radiation as a result of the Stefan-Boltzmann equation,
      !! which is the sum of the calculated values at the surface temperature  at the 'old' 
      !! temperature and the value that corresponds to the difference between the old temperature
      !! and the temperature at the 'new' timestep.
      !! \latexonly 
      !!     \input{enerbilflux3.tex}
      !! \endlatexonly
      lwup(ji) = emis(ji) * c_stefan * temp_sol(ji)**4 + &
           &     quatre * emis(ji) * c_stefan * temp_sol(ji)**3 * &
           &     (temp_sol_new(ji) - temp_sol(ji))      
      
      !! We then add to this value that from the reflected longwave radiation:
      !! \latexonly 
      !!     \input{enerbilflux4.tex}
      !! \endlatexonly
      lwup(ji) = lwup(ji)  +  (un - emis(ji)) * lwdown(ji)
      
      !! The radiative surface temperature is calculated by inverting the Stefan-Boltzmann relation:
      !! \latexonly 
      !!     \input{enerbilflux5.tex}
      !! \endlatexonly
      !! The implicit solution computes an emitted long-wave flux which is a limited Taylor expansion
      !! of the future surface temperature around the current values. Thus the long-wave flux does not 
      !! correspond to the new surface temperature but some intermediate value. So we need to deduce the
      !! radiative surface temperature to which this flux corresponds.
      !!
      tsol_rad(ji) = (lwup(ji)/ (emis(ji) * c_stefan)) **(1./quatre)
     
      !! qsurf (the surface specific humidity) is a simple diagnostic which will be used by the 
      !! GCM to compute the dependence of of the surface layer stability on moisture.
      !! \latexonly 
      !!     \input{enerbilflux6.tex}
      !! \endlatexonly
      qsurf(ji) = (vbeta1(ji) * (un - vbeta5(ji)) + vbeta5(ji)) * qsol_sat_new(ji) + &
           & (un - vbeta1(ji))*(un - vbeta5(ji)) * vbeta(ji) * valpha(ji) * qsol_sat_new(ji)
      
      !! \latexonly 
      !!     \input{enerbilflux7.tex}
      !! \endlatexonly
      qsurf(ji) = MAX(qsurf(ji), qair(ji))
      
      !! Net downward radiation is the sum of the down-welling less the up-welling long wave flux, plus
      !! the short wave radiation.
      !! \latexonly 
      !!     \emissivity absorbed lw radiationinput{enerbilflux8.tex}
      !! \endlatexonly
      netrad(ji) = lwdown(ji) + swnet(ji) - lwup(ji) 
      
      !! 'vevapp' is the sum of the total evaporative processes (snow plus non-snow processes).
      !! \latexonly 
      !!     \input{enerbilflux9.tex}
      !! \endlatexonly
      vevapp(ji) = dtradia * rau(ji) * qc * (vbeta1(ji) * (un - vbeta5(ji)) + vbeta5(ji)) * &
           & (qsol_sat_new(ji) - qair(ji)) + &
           &  dtradia * rau(ji) * qc * (un - vbeta1(ji))*(un-vbeta5(ji)) * vbeta(ji) * &
           & (valpha(ji) * qsol_sat_new(ji) - qair(ji))
      
      !! The total latent heat flux is the sum of the snow plus non-snow processes.
      !! \latexonly 
      !!     \input{enerbilflux10.tex}
      !! \endlatexonly
      fluxlat(ji) = chalsu0 * rau(ji) * qc * vbeta1(ji) * (un - vbeta5(ji)) * &
           & (qsol_sat_new(ji) - qair(ji)) + &
           &  chalev0 * rau(ji) * qc * vbeta5(ji) *&
           & (qsol_sat_new(ji) - qair(ji)) + &
           &  chalev0 * rau(ji) * qc * (un - vbeta1(ji)) * (un - vbeta5(ji)) * vbeta(ji) * &
           & (valpha(ji) * qsol_sat_new(ji) - qair(ji))
      
      !! The sublimination flux concerns is calculated using vbeta1, the snow resistance.
      !! \latexonly 
      !!     \input{enerbilflux11.tex}
      !! \endlatexonly
      fluxsubli(ji) = chalsu0 * rau(ji) * qc * vbeta1(ji) * (un - vbeta5(ji)) * &
           & (qsol_sat_new(ji) - qair(ji)) 
      
      !! The sensible heat flux is a factor of the difference between the new surface static energy
      !! and the potential energy of air.
      !! \latexonly 
      !!     \input{enerbilflux12.tex}
      !! \endlatexonly
      fluxsens(ji) =  rau(ji) * qc * (psnew(ji) - epot_air(ji))
      
      !! This is the net longwave downwards radiation.
      !! \latexonly 
      !!     \input{enerbilflux13.tex}
      !! \endlatexonly
      lwnet(ji) = lwdown(ji) - lwup(ji)
      
      !! Diagnoses the potential evaporation used for the fluxes (evapot)
      !! \latexonly 
      !!     \input{enerbilflux14.tex}
      !! \endlatexonly  
      evapot(ji) =  MAX(zero, dtradia * rau(ji) * qc * (qsol_sat_new(ji) - qair(ji)))
     
      !! From definition we can say:
      !! \latexonly 
      !!     \input{enerbilflux15.tex}
      !! \endlatexonly 
      tair(ji)  =  epot_air(ji) / cp_air
      !! To calculate net energy flux into the snowpack
      IF (ok_explicitsnow) THEN
          PHPSNOW(ji) = precip_rain(ji)*(4.218E+3)*(MAX(tp_00,temp_air(ji))-tp_00)/dtradia ! (w/m2)
          pgflux(ji)  = netrad(ji) - fluxsens(ji) - fluxlat(ji) + PHPSNOW(ji)
      ENDIF


      !! To get the extra energy used to melt the snowpack
      IF (ok_explicitsnow .AND. temp_sol_new (ji) > tp_00 .AND. SUM(snowdz(ji,:)) .GT. zero .AND. snowcap(ji) .GT. zero) THEN

         lwup_tmp(ji) = emis(ji) * c_stefan * temp_sol(ji)**4 + &
           &     quatre * emis(ji) * c_stefan * temp_sol(ji)**3 * &
           &     (tp_00 - temp_sol(ji))

         lwup_tmp(ji) = lwup_tmp(ji)  +  (un - emis(ji)) * lwdown(ji)
         tsol_rad_tmp(ji) = emis(ji) * c_stefan * temp_sol(ji)**4 + lwup_tmp(ji)
         netrad_tmp(ji) = lwdown(ji) + swnet(ji) - lwup_tmp(ji)
         fluxsens_tmp(ji) =  rau(ji) * qc * cp_air * (tp_00 - epot_air(ji)/cp_air)
         fluxlat_tmp(ji) = chalsu0 * rau(ji) * qc * vbeta1(ji) * (un-vbeta5(ji)) * &
                        & (qsol_sat_tmp(ji) - qair(ji)) + &
                        &  chalev0 * rau(ji) * qc * vbeta5(ji) *&
                        & (qsol_sat_tmp(ji) - qair(ji)) + &
                        &  chalev0 * rau(ji) * qc * (un - vbeta1(ji)) * (un-vbeta5(ji))* vbeta(ji) * &
                        & (valpha(ji) * qsol_sat_tmp(ji) - qair(ji))

         zgflux(ji)  = netrad_tmp(ji) - fluxsens_tmp(ji) - fluxlat_tmp(ji)+PHPSNOW(ji)

         temp_sol_add(ji) = -(pgflux(ji) - zgflux(ji))*dtradia/snowcap(ji)

         pgflux(ji) = zgflux(ji)

      ELSE

         temp_sol_add(ji) = zero

      ENDIF

     ENDDO 

  !! 3. Define qsat_air with the subroutine src_parameter:

    CALL qsatcalc(kjpindex, tair, pb, qsat_air)

    CALL dev_qsatcalc(kjpindex, tair, pb, grad_qsat)
 
  ! grad_qsat(:)= (qsol_sat_new(:)- qsat_air(:)) / ((psnew(:) - epot_air(:)) / cp_air) ! * dtradia
    !- Penser a sortir evapot en meme temps qu'evapot_corr tdo.

    warning_correction(:)=.FALSE.
    DO ji=1,kjpindex
      
      !! \latexonly 
      !!     \input{enerbilflux16.tex}
      !! \endlatexonly 
      speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
      
      !! \latexonly 
      !!     \input{enerbilflux17.tex}
      !! \endlatexonly 
      qc = speed * q_cdrag(ji)
       
      !! Derive the potential as defined by Penman & Monteith (Monteith, 1965) based on the correction 
      !! term developed by Chris Milly (1992).
       IF ((evapot(ji) .GT. zero) .AND. ((psnew(ji) - epot_air(ji)) .NE. zero )) THEN
          !
          !! \latexonly 
          !!     \input{enerbilflux18.tex}
          !! \endlatexonly 
          correction =  (quatre * emis(ji) * c_stefan * tair(ji)**3 + rau(ji) * qc * cp_air + &
               &                  chalev0 * rau(ji) * qc * grad_qsat(ji) * vevapp(ji) / evapot(ji) )
          
          !! \latexonly 
          !!     \input{enerbilflux19.tex}
          !! \endlatexonly 
          IF (ABS(correction) .GT. min_sechiba) THEN
             correction = chalev0 * rau(ji) * qc * grad_qsat(ji) * (un - vevapp(ji)/evapot(ji)) / correction
          ELSE
             warning_correction(ji)=.TRUE.
          ENDIF
       ELSE
          correction = zero
       ENDIF
       correction = MAX (zero, correction)
      
      !! \latexonly 
      !!     \input{enerbilflux20.tex}
      !! \endlatexonly 
       evapot_corr(ji) = evapot(ji) / (un + correction)
       
    ENDDO
    IF ( ANY(warning_correction) ) THEN
       DO ji=1,kjpindex
          IF ( warning_correction(ji) ) THEN
             WRITE(numout,*) ji,"Denominateur de la correction de milly nul! Aucune correction appliquee"
          ENDIF
       ENDDO
    ENDIF
    IF (long_print) WRITE (numout,*) ' enerbil_flux done '

  END SUBROUTINE enerbil_flux


  !!  ===========================================================================================================================
  !! SUBROUTINE		 			: enerbil_evapveg
  !!
  !>\BRIEF					Splits total evaporation loss into individual process
  !! components and calculates the assimilation.
  !!
  !! DESCRIPTION				: Based on the estimation of the fluxes in enerbil_flux, 
  !! this routine splits the total evaporation into transpiration interception loss from vegetation, 
  !! bare soil evaporation and snow sublimation. It then calculates the assimilation.
  !!
  !! MAIN OUTPUT VARIABLE(S)	                : vevapsno, vevapnu, transpir, vevapwet
  !!
  !! REFERENCE(S)					:
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ, Technical note, available (22/12/11):
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
  !! sur le cycle de l'eau, PhD Thesis, available (22/12/11): 
  !! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663. 
  !!
  !! FLOWCHART   : None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE enerbil_evapveg (kjpindex, dtradia, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, &
     & rau, u, v, q_cdrag, qair, humrel, vevapsno, vevapnu , vevapflo, vevapwet, &
     & transpir, transpot, evapot, &
     & snowdz,snowrho,evapot_corr)

    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjpindex          !! Domain size (-)
    REAL(r_std), INTENT(in)                             :: dtradia           !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: vbeta1            !! Snow resistance (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: vbeta4            !! Bare soil resistance (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: vbeta5            !! Floodplains resistance
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: rau               !! Density (kg m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: u, v              !! Wind velocity in directions u and v (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: q_cdrag           !! Surface drag (m s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: qair              !! Lowest level specific humidity (kg kg^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: evapot            !! Soil Potential Evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex, nvm), INTENT (in)  :: humrel            !! Soil moisture stress (within range 0 to 1)
!!$ DS 15022011 humrel was used in a previous version of Orchidee, developped by Nathalie. Need to be discussed if it should be introduce again             
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: vbeta2            !! Interception resistance (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: vbeta3            !! Vegetation resistance (-)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)   :: vbeta3pot         !! Vegetation resistance for potential transpiration
    REAL(r_std), DIMENSION(kjpindex,nsnow),INTENT(in)   :: snowrho           !! Snow density
    REAL(r_std), DIMENSION(kjpindex,nsnow),INTENT(in):: snowdz               !! Snow depth 
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: vevapsno          !! Snow evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: vevapnu           !! Bare soil evaporation (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: vevapflo          !! Floodplains evaporation
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: transpir          !! Transpiration (mm day^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: transpot          !! Potential Transpiration
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)  :: vevapwet          !! Interception (mm day^{-1})

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION(kjpindex),INTENT(inout)      :: evapot_corr       !! Correction term in evapotranspiration
    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jv
    REAL(r_std), DIMENSION(kjpindex)                    :: xx
    REAL(r_std), DIMENSION(kjpindex)                    :: vbeta2sum, vbeta3sum
    REAL(r_std)                                    	:: speed	     !! Speed (m s^{-1})
    REAL(r_std), DIMENSION(kjpindex)                    :: ZSNOWEVAPS,ZEVAPCOR,ZSNOWDZ

!_ ==============================================================================================================================
    ! initialisation: utile pour calculer l'evaporation des floodplains dans lesquelles il y a de la vegetation
    vbeta2sum(:) = 0.
    vbeta3sum(:) = 0.
    DO jv = 1, nvm
      vbeta2sum(:) = vbeta2sum(:) + vbeta2(:,jv)
      vbeta3sum(:) = vbeta3sum(:) + vbeta3(:,jv)
    ENDDO 
    
    !! 1. Compute vevapsno (evaporation from snow) and vevapnu (bare soil evaporation)
        
    DO ji=1,kjpindex 

      !! \latexonly 
      !!     \input{enerbilevapveg1.tex}
      !! \endlatexonly
      speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))

      
      !! 1.1 Snow sublimation
      !! \latexonly 
      !!     \input{enerbilevapveg2.tex}
      !! \endlatexonly
      vevapsno(ji) = (un - vbeta5(ji)) * vbeta1(ji) * dtradia * rau(ji) * speed * q_cdrag(ji) * (qsol_sat_new(ji) - qair(ji))


      ! For the case of ok_explicitsnow activated, remove mass from uppermost snow layer in 
      ! response to sublimation. Do not take into account evaporation. The evaporation on the
      ! snow surface could be possible in this case since the new snow model allows liquid water to be existed. 
      IF ( ok_explicitsnow ) THEN
         ! Sublimation of snow ice:
         ! Reduce snow layer thickness and total snow depth if sublimation. 
         ! Add to correction term if potential sublimation exceeds available snow cover.
         
         ZSNOWEVAPS(ji)  = vevapsno(ji)/(chalsu0*snowrho(ji,1))
         ZSNOWDZ(ji)     = snowdz(ji,1) - ZSNOWEVAPS(ji)

         ZEVAPCOR(ji)    = MAX(0.0,-ZSNOWDZ(ji))*snowrho(ji,1)/dtradia

         ! Update evaporation correction term:
          evapot_corr(ji) = evapot_corr(ji) + ZEVAPCOR(ji)
      ENDIF

      !! 1.2 Bare soil evaporation
      !! \latexonly 
      !!     \input{enerbilevapveg3.tex}
      !! \endlatexonly
      vevapnu(ji) = (un - vbeta1(ji)) * (un-vbeta5(ji)) * vbeta4(ji) * dtradia * rau(ji) * speed * q_cdrag(ji) &
         & * (qsol_sat_new(ji) - qair(ji))
      !
      ! 1.3 floodplains evaporation - transpiration et interception prioritaires dans les floodplains
      !
      vevapflo(ji) = vbeta5(ji) * (1 - vbeta2sum(ji) - vbeta3sum(ji)) &
           & * dtradia * rau(ji) * speed * q_cdrag(ji) * (qsol_sat_new(ji) - qair(ji))

    END DO

    !! 2. Compute transpir (transpiration) and vevapwet (interception)
  
    !! Preliminaries
    DO ji = 1, kjpindex
       
       !! \latexonly 
       !!     \input{enerbilevapveg4.tex}
       !! \endlatexonly
       speed = MAX(min_wind, SQRT (u(ji)*u(ji) + v(ji)*v(ji)))
       
       !! \latexonly 
       !!     \input{enerbilevapveg5.tex}
       !! \endlatexonly
       xx(ji) = dtradia * (un-vbeta1(ji)) * (qsol_sat_new(ji)-qair(ji)) * rau(ji) * speed * q_cdrag(ji)
     
    ENDDO
    
    DO jv=1,nvm 
      DO ji=1,kjpindex 
        
        !! 2.1 Calculate interception loss
        !! \latexonly 
        !!     \input{enerbilevapveg6.tex}
        !! \endlatexonly
        vevapwet(ji,jv) = xx(ji) * vbeta2(ji,jv)
        ! 
        !! 2.2 Calculate transpiration
        !! \latexonly 
        !!     \input{enerbilevapveg7.tex}
        !! \endlatexonly
        transpir (ji,jv) = xx(ji) * vbeta3(ji,jv)

        transpot(ji,jv) = xx(ji) * vbeta3pot(ji,jv)

      END DO
    END DO

    

    IF (long_print) WRITE (numout,*) ' enerbil_evapveg done '

  END SUBROUTINE enerbil_evapveg

 
  !!  =============================================================================================================================
  !! SUBROUTINE		 			: enerbil_fusion
  !!
  !>\BRIEF					Computes new soil temperature due to
  !! ice and snow melt. 
  !!
  !! DESCRIPTION				: This routine computes new soil temperature due to
  !! ice and snow melt. It is the second part of main routine for enerbil module, and is called every
  !! time step. 
  !!
  !! MAIN OUTPUT VARIABLE(S)			: temp_sol_new, fusion
  !!
  !! REFERENCE(S)				:
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ. Technical note, available (22/12/11):
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
  !! sur le cycle de l'eau. PhD Thesis, available (22/12/11):
  !! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663.
  !!
  !! FLOWCHART  : None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE enerbil_fusion (kjpindex, dtradia, tot_melt, soilcap, snowdz, &
                             temp_sol_new, fusion)


    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex       !! Domain size (-)
    REAL(r_std), INTENT(in)                                  :: dtradia        !! Time step (s)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: tot_melt       !! Total melt (??)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: soilcap        !! Soil heat capacity (J K^{-1])
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in)        :: snowdz         !! Snow layer thickness (m)
    
    !! 0.2 Output variables
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fusion         !! Fusion (??)

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: temp_sol_new   !! New soil temperature (K)

    !! 0.4 Local variables
    INTEGER(i_std)                                	     :: ji	       !! (-)

  !_ ==============================================================================================================================
    
    !! Initialisation
    IF (long_print) WRITE (numout,*) ' enerbil_fusion start ', MINVAL(soilcap), MINLOC(soilcap),&
         & MAXVAL(soilcap), MAXLOC(soilcap)
    
    !! 1. Calculate new soil temperature due to ice and snow melt
   
    IF  (ok_explicitsnow) THEN
       !! Case with explicit snow activated :
       !! Surface temperature is forced to zero celcius if its value is larger than melting point 

       DO ji=1,kjpindex 
          IF  (SUM(snowdz(ji,:)) .GT. 0.0) THEN
             IF (temp_sol_new(ji) .GE. tp_00) THEN
                temp_sol_new(ji) = tp_00
             ENDIF
          END IF
       END DO

    ELSE
       !! Default case :
       !! Surface temperature is reduced if there is snow melt

       DO ji=1,kjpindex 
          !! \latexonly 
          !!     \input{enerbilfusion1.tex}
          !! \endlatexonly
          fusion(ji) = tot_melt(ji) * chalfu0 / dtradia
          
          !! \latexonly 
          !!     \input{enerbilfusion2.tex}
          !! \endlatexonly

          temp_sol_new(ji) = temp_sol_new(ji) - ((tot_melt(ji) * chalfu0) / soilcap(ji))

       END DO
    ENDIF
           
    IF (long_print) WRITE (numout,*) ' enerbil_fusion done '

  END SUBROUTINE enerbil_fusion


  !!  =============================================================================================================================
  !! SUBROUTINE		 			: enerbil_t2mdiag
  !!
  !> BRIEF					This routine diagnoses the air temperature at 2
  !! metres above ground.
  !!
  !! DESCRIPTION				: The 2 metre above ground air temperature is currently
  !! set to be the same as the general air temperaturem, as the surface model currently simulates the 
  !! above ground atmosphere as a uniform box (and the lowest layer of the LMDZ atmospheric model).
  !!
  !! MAIN OUTPUT VARIABLE(S)			: t2mdiag
  !!
  !! REFERENCE(S)				:
  !! - Best, MJ, Beljaars, A, Polcher, J & Viterbo, P, 2004. A proposed structure for coupling tiled
  !! surfaces with the planetary boundary layer. Journal of Hydrometeorology, 5, pp.1271-1278
  !! - de Noblet-Ducoudré, N, Laval, K & Perrier, A, 1993. SECHIBA, a new set of parameterisations
  !! of the hydrologic exchanges at the land-atmosphere interface within the LMD Atmospheric General
  !! Circulation Model. Journal of Climate, 6, pp.248-273
  !! - Dufresne, J-L & Ghattas, J, 2009. Description du schéma de la couche limite turbulente et la
  !! interface avec la surface planetaire dans LMDZ, Technical note, available from:
  !! http://lmdz.lmd.jussieu.fr/developpeurs/notes-techniques/ressources/pbl_surface.pdf
  !! - Guimberteau, M, 2010. Modélisation de l'hydrologie continentale et influences de l'irrigation
  !! sur le cycle de l'eau, PhD Thesis, available from:
  !! http://www.sisyphe.upmc.fr/~guimberteau/docs/manuscrit_these.pdf
  !! - Polcher, J. McAvaney, B, Viterbo, P, Gaertner, MA, Hahmann, A, Mahfouf, J-F, Noilhan, J
  !! Phillips, TJ, Pitman, AJ, Schlosser, CA, Schulz, J-P, Timbal, B, Verseghy, D &
  !! Xue, Y, 1998. A proposal for a general interface between land surface schemes and
  !! general circulation models. Global and Planetary Change, 19, pp.261-276
  !! - Richtmeyer, RD, Morton, KW, 1967. Difference Methods for Initial-Value Problems.
  !! Interscience Publishers\n
  !! - Schulz, Jan-Peter, Lydia Dümenil, Jan Polcher, 2001: On the Land Surface–Atmosphere 
  !! Coupling and Its Impact in a Single-Column Atmospheric Model. J. Appl. Meteor., 40, 642–663.
  !!
  !! FLOWCHART  : None
  !! \n
  !_ ==============================================================================================================================

  SUBROUTINE enerbil_t2mdiag (kjpindex, temp_air, t2mdiag)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                         		:: kjpindex       !! Domain size (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      		:: temp_air       !! Air temperature (K)
    
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)    		:: t2mdiag        !! 2-metre temperature (K)

    !! 0.3 Modified variables

    !! 0.4 Local variables

  !_ ==============================================================================================================================

    !! The 2 metre air temperature is assumed to be the same of the value of temp_air - the atmospheric air
    !! temperature, which also taken as the temperature at the lowest level of the LMDZ. 't2mdiag' is just
    !! a value used for diagnostic purposes, admd it is no applied amywhere else.
    !! \latexonly 
 	!!     \input{enerbilt2mdiag1.tex}
    !! \endlatexonly
    t2mdiag(:) = temp_air(:)

    IF (long_print) WRITE (numout,*) ' enerbil_t2mdiag done '

  END SUBROUTINE enerbil_t2mdiag



END MODULE enerbil
