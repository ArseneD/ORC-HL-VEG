! =================================================================================================================================
! MODULE       : thermosoil
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Calculates the soil temperatures by solving the heat
!! diffusion equation within the soil
!!
!!\n DESCRIPTION : General important informations about the numerical scheme and
!!                 the soil vertical discretization:\n
!!               - the soil is divided into "ngrnd" (=7 by default) layers, reaching to as
!!                 deep as 5.5m down within the soil, with thiscknesses
!!                 following a geometric series of ration 2.\n
!!               - "jg" is usually used as the index going from 1 to ngrnd to describe the
!!                  layers, from top (jg=1) to bottom (jg=ngrnd)\n
!!               - the thermal numerical scheme is implicit finite differences.\n
!!                 -- When it is resolved in thermosoil_profile at the present timestep t, the
!!                 dependancy from the previous timestep (t-1) is hidden in the
!!                 integration coefficients cgrnd and dgrnd, which are therefore
!!                 calculated at the very end of thermosoil_main (call to
!!                 thermosoil_coef) for use in the next timestep.\n
!!                 -- At timestep t, the system becomes :\n 
!!
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!
!!                 (the bottom boundary condition has been used to obtained this equation).\n
!!                 To solve it, the uppermost soil temperature T(1) is required.
!!                 It is obtained from the surface temperature Ts, which is
!!                 considered a linear extrapolation of T(1) and T(2)\n
!!
!!                           Ts=(1-lambda)*T(1) -lambda*T(2) \n 
!!                                      -- EQ2--\n
!!
!!                 -- caveat 1 : Ts is called 'temp_soil_new' in this routine,
!!                 don' t act.\n
!!                 -- caveat 2 : actually, the surface temperature at time t Ts
!!                 depends on the soil temperature at time t through the
!!                 ground heat flux. This is again implicitly solved, with Ts(t)
!!                 expressed as :\n
!!
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflux+otherfluxes(Ts(t))\n 
!!                                      -- EQ3 --\n
!!
!!                 and the dependency from the previous timestep is hidden in
!!                 soilcap and soilflux (apparent surface heat capacity and heat
!!                 flux respectively). Soilcap and soilflux are therefore
!!                 calculated at the previsou timestep, at the very end of thermosoil
!!                 (final call to thermosoil_coef) and stored to be used at the next time step.
!!                 At timestep t, EQ3 is solved for Ts in enerbil, and Ts
!!                 is used in thermosoil to get T(1) and solve EQ1.\n
!!
!! - lambda is the @tex $\mu$ @endtex of F. Hourdin' s PhD thesis, equation (A28); ie the
!! coefficient of the linear extrapolation of Ts (surface temperature) from T1 and T2 (ptn(jg=1) and ptn(jg=2)), so that:\n
!! Ts= (1+lambda)*T(1)-lambda*T(2) --EQ2-- \n
!! lambda = (zz_coeff(1))/((zz_coef(2)-zz_coef(1))) \n
!!
!! - cstgrnd is the attenuation depth of the diurnal temperature signal
!! (period : one_day) as a result of the heat conduction equation
!! with no coefficients :
!!\latexonly
!!\input{thermosoil_var_init0.tex}
!!\endlatexonly
!!  -- EQ4 --\n
!! This equation results from the change of variables :
!! z' =z*sqrt(Cp/K) where z' is the new depth (homogeneous
!! to sqrt(time) ), z the real depth (in m), Cp and K the soil heat
!! capacity and conductivity respectively.\n
!!
!! the attenuation depth of a diurnal thermal signal for EQ4 is therefore homogeneous to sqrt(time) and
!! equals : \n
!! cstgrnd = sqrt(oneday/Pi)
!!
!! - lskin is the attenuation depth of the diurnal temperature signal
!! (period : one_day) within the soil for the complete heat conduction equation
!! (ie : with coefficients)
!!\latexonly
!!\input{thermosoil_var_init00.tex}
!!\endlatexonly
!! -- EQ5 --  \n
!! it can be retrieved from cstgrnd using the change of variable  z' =z*sqrt(Cp/K):\n
!! lskin = sqrt(K/Cp)*cstgrnd =  sqrt(K/Cp)*sqrt(oneday//Pi)\n
!! 
!! In thermosoil, the ratio lskin/cstgrnd is frequently used as the
!! multiplicative factor to go from
!!'adimensional' depths (like z' ) to real depths (z). z' is not really
!! adimensional but is reffered to like this in the code.
!!
!!---------------------------------------------------------------------------------------
!!   Modified by Dmitry Khvorostyanov and Gerhard Krinner 12-14/12/06 to account
!for permafrost
!!
!!    - new subroutine 'thermosoil_getdiff' that computes soil heat conductivity
!!      and heat capacity with account for liquid and frozen phases
!!
!!    - new subroutine 'thermosoil_wlupdate' that computes long-term soil
!!      humidity ensuring energy conservation when soil freezes
!!
!!    - in 'thermosoil_coef' and 'thermosoil_var_init' the part computing the
!!      soil capa and kappa has been rewritten in terms of the new routine 'thermosoil_getdiff'
!!
!!    - in the call to 'thermosoil_var_init' the variable 'snow' is now passed
!!      as an input argument in order to be able to use 'thermosoil_getdiff'
!!
!!    -  'thermosoil_wlupdate' is called in 'thermosoil_main', just after
!!       'thermosoil_humlev', to update the long-term humidity
!!
!!    - new module constants related to permafrost have been added
!!
!!    - modifications related to the thermosoil autonomy (optional output using
!!      flio, no use of restart files, initial and boundary conditions specified in the call to
!!      thermosoil_main)
!!---------------------------------------------------------------------------------------
!!  Modified by Charlie Koven 2008-2010 and Tao Wang 2014 to:
!!  
!!  - added PFT dimension to soil thermal calculations
!!
!!  - three-layer snow module from ISBA-ES
!!  
!!  - take into account soil carbon in the thermal properties of soil 
!!    Two possible modes, with separate organic layer on top of soil, or mixed mieral/organic soils
!!  
!!  - take into account exothermic heat of decomposition in heat budget of soil layers
!!  
!!  - output some diagnostic properties (soil moisture, etc) on soil temperature grid for use in 
!!    calculations within the permafrost soil carbon code.
!!
!!---------------------------------------------------------------------------------------
!!
!!
!! RECENT CHANGE(S) : None
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE-MICT/ORCHIDEE/src_sechiba/thermosoil.f90 $
!! $Date: 2013-07-29 16:10:47 +0200 (Mon, 29 Jul 2013) $
!! $Revision: 1397 $
!! \n
!_ ================================================================================================================================

MODULE thermosoil

  USE ioipsl_para
  USE xios_orchidee
  USE constantes
  USE constantes_soil
  USE sechiba_io
  USE grid
  USE pft_parameters_var
  USE constantes_var
  USE mod_orchidee_para

  IMPLICIT NONE

  !private and public routines :
  PRIVATE
  PUBLIC :: thermosoil_main,thermosoil_clear,thermosoil_vert_axes,thermosoil_levels 

  LOGICAL, SAVE                                   :: l_first_thermosoil=.TRUE.!! does the initialisation of the routine 
                                                                              !! (true/false)
!$OMP THREADPRIVATE(l_first_thermosoil)
  CHARACTER(LEN=80) , SAVE                        :: var_name                 !! To store variables names for the 
                                                                              !! input-outputs dealt with by IOIPSL
!$OMP THREADPRIVATE(var_name)
  REAL(r_std), SAVE                               :: lambda, cstgrnd, lskin   !! See Module description
!$OMP THREADPRIVATE(lambda, cstgrnd, lskin)
  REAL(r_std), SAVE                               :: fz1               !! usefull constants for diverse use
!$OMP THREADPRIVATE(fz1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: ptn                      !! vertically discretized 
                                                                              !! soil temperatures @tex ($K$) @endtex. 
!$OMP THREADPRIVATE(ptn)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: ptn_pftmean             !! Different levels soil temperature, mean across all pfts
!$OMP THREADPRIVATE(ptn_pftmean)
  REAL(r_std),  ALLOCATABLE, SAVE, DIMENSION (:)            :: zz                       !! depths of the soil thermal numerical nodes. 
                                                                              !! Caveats: they are not exactly the centers of the
                                                                              !! thermal layers, see the calculation in 
                                                                              !! ::thermosoil_var_init  @tex ($m$) @endtex.
!$OMP THREADPRIVATE(zz)
  REAL(r_std),  ALLOCATABLE,SAVE, DIMENSION (:)            :: zz_coef                  !! depths of the boundaries of the thermal layers,
                                                                              !! see the calculation in 
                                                                              !! thermosoil_var_init  @tex ($m$) @endtex.
!$OMP THREADPRIVATE(zz_coef)
  REAL(r_std),  ALLOCATABLE,SAVE, DIMENSION (:)            :: dz1                      !! numerical constant used in the thermal numerical
                                                                              !! scheme  @tex ($m^{-1}$) @endtex. ; it corresponds
                                                                              !! to the coefficient  @tex $d_k$ @endtex of equation
                                                                              !! (A.12) in F. Hourdin PhD thesis.
!$OMP THREADPRIVATE(dz1)
  REAL(r_std),  ALLOCATABLE,SAVE, DIMENSION (:)            :: dz2                      !! thicknesses of the thermal layers  @tex ($m$)
                                                                              !! @endtex; typically: 
                                                                              !! dz2(jg)=zz_coef(jg+1)-zz_coef(jg); calculated once 
                                                                              !! and for all in thermosoil_var_init
!$OMP THREADPRIVATE(dz2)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: z1                       !! constant of the numerical scheme; it is an 
                                                                              !! intermediate buffer for the calculation of the 
                                                                              !! integration coefficients cgrnd and dgrnd.
!$OMP THREADPRIVATE(z1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: cgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
!$OMP THREADPRIVATE(cgrnd)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: dgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
!$OMP THREADPRIVATE(dgrnd)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: pcapa                    !! volumetric vertically discretized soil heat 
                                                                              !! capacity  @tex ($J K^{-1} m^{-3}$) @endtex. 
                                                                              !! It depends on the soil
                                                                              !! moisture content (wetdiag) and is calculated at 
                                                                              !! each time step in thermosoil_coef.
!$OMP THREADPRIVATE(pcapa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: pkappa                   !! vertically discretized soil thermal conductivity 
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex. Same as pcapa.
!$OMP THREADPRIVATE(pkappa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: zdz1                     !! numerical constant of the numerical scheme; it is
                                                                              !! an intermediate buffer for the calculation of the 
                                                                              !! integration coefficients cgrnd and dgrnd 
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex 
!$OMP THREADPRIVATE(zdz1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: zdz2                     !! numerical constant of the numerical scheme; it is 
                                                                              !! an intermediate buffer for the calculation of the 
                                                                              !! integration coefficients cgrnd and dgrnd
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex
!$OMP THREADPRIVATE(zdz2)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: pcapa_en                 !! heat capacity used for surfheat_incr and 
                                                                              !! coldcont_incr 
!$OMP THREADPRIVATE(pcapa_en)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: ptn_beg                  !! vertically discretized temperature at the 
                                                                              !! beginning of the time step  @tex ($K$) @endtex; 
                                                                              !! is used in 
                                                                              !! thermosoil_energy for energy-related diagnostic of
                                                                              !! the routine.
!$OMP THREADPRIVATE(ptn_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: temp_sol_beg             !! Surface temperature at the beginning of the 
                                                                              !! timestep  @tex ($K$) @endtex
!$OMP THREADPRIVATE(temp_sol_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: surfheat_incr            !! Change in soil heat content during the timestep 
                                                                              !!  @tex ($J$) @endtex.
!$OMP THREADPRIVATE(surfheat_incr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: coldcont_incr            !! Change in snow heat content  @tex ($J$) @endtex.
!$OMP THREADPRIVATE(coldcont_incr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:) :: wetdiag                  !! Soil wetness on the thermodynamical levels 
                                                                              !! (1, ngrnd) (0-1, dimensionless). corresponds to the
                                                                              !! relative soil humidity to the wilting point when 
                                                                              !! the 11-layers hydrology (hydrol) is used, see more
                                                                              !! precisions in thermosoil_humlev.
!$OMP THREADPRIVATE(wetdiag)
 !Isa------------------------------------------------------------------------
  !  Permafrost-related constants: 
  !
!  ! temperature range over which soil freezing takes place (K)
!  REAL(r_std), PARAMETER        	 :: fr_dT =2
!
!  ! Isa : the new porosity is a correction, old value was 0.15
!  REAL(r_std), PARAMETER                :: poros = .41 ! (moyenne des mcs sur classif USDA)

! heat capacity of ice-filled saturated soil (J m**-3 K**-1, at porosity of 0.41)
!Isa : the new heat capacity for saturated frozen soil is a correction ; old
!value was 1.81e6 (corresponding to poros = 0.15), the new value is computed
!based on : so_capa_ice = so_capa_dry + poros*capa_ice*rho_ice, with
!capa_ice=2.06 J/Kg/K at 0°C.
! this is a divergence from Gouttevin et al., 2012, where so_capa_ice = 2.3e6
! was considered (mistake...)
 
  INTEGER, PARAMETER	   	 :: bavard = 4
  INTEGER, SAVE   	         :: flioid	  !! flio output file ID
REAL(r_std),SAVE :: so_cond = 1.5396
    REAL(r_std),SAVE :: so_capa = 2.0514e+6

!Isa 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: profil_froz
!$OMP THREADPRIVATE(profil_froz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: wetdiaglong          !! Long-term soil humidity (for permafrost) if ok_freeze_thermix ; wetdiag sinon.
!$OMP THREADPRIVATE(wetdiaglong)
!  LOGICAL, SAVE    				   :: ok_freeze_thermix
!$OMP THREADPRIVATE(ok_freeze_thermix)
!Isa E
! CR: option supl pour converg avec trunc
!        LOGICAL, SAVE    :: ok_thermix_trunc
!$OMP THREADPRIVATE(ok_thermix_trunc)
        LOGICAL, SAVE    :: ok_wetdiaglong
!$OMP THREADPRIVATE(ok_wetdiaglong)
!        LOGICAL, SAVE    :: ok_converge_isaorig ! on le met dans constantes

    REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:,:,:)    :: pcappa_supp
!$OMP THREADPRIVATE(pcappa_supp)
    REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:,:,:)    :: E_sol_lat_couche
!$OMP THREADPRIVATE(E_sol_lat_couche)
    REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:,:)      :: pcappa_supp_pftmean
!$OMP THREADPRIVATE(pcappa_supp_pftmean)
    REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:,:)      :: E_sol_lat_couche_pftmean
!$OMP THREADPRIVATE(E_sol_lat_couche_pftmean)

!    LOGICAL, SAVE     				   :: ok_Ecorr
!$OMP THREADPRIVATE(ok_Ecorr)
!Isa permafrost_map
    REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:) :: overburden
!$OMP THREADPRIVATE(overburden)
    REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:) :: excess_ice
!$OMP THREADPRIVATE(excess_ice)
    REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:) :: permafrost
!$OMP THREADPRIVATE(permafrost)
!    LOGICAL, SAVE     				   :: read_permafrost_map
!$OMP THREADPRIVATE(read_permafrost_map)
!Isa reftemp
    REAL(r_std), ALLOCATABLE, SAVE,DIMENSION(:,:) :: reftemp
!$OMP THREADPRIVATE(reftemp)
!    LOGICAL, SAVE     				   :: read_reftemp
!$OMP THREADPRIVATE(read_reftemp)
    ! end isa

!Tao
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: snow_thick
!$OMP THREADPRIVATE(snow_thick)
!Vertical Permafrost Carbon
  LOGICAL, SAVE                                :: use_toporganiclayer_tempdiff = .FALSE. 
  LOGICAL, SAVE                                :: use_soilc_tempdiff = .TRUE. 
  LOGICAL, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: veget_mask_2d
  LOGICAL, SAVE                                :: satsoil = .FALSE.
!  REAL(r_std), PARAMETER                :: so_capa_ice = 2.11e6 !!Add from Tao.. but not running. On constantes_var.f90


CONTAINS

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_main
!!
!>\BRIEF        Thermosoil_main computes the soil thermal properties and dynamics, ie solves
!! the heat diffusion equation within the soil. The soil temperature profile is
!! then interpolated onto the diagnostic axis.
!!
!! DESCRIPTION : The resolution of the soil heat diffusion equation 
!! relies on a numerical finite-difference implicit scheme
!! fully described in the reference and in the header of the thermosoil module.
!! - The dependency of the previous timestep hidden in the 
!! integration coefficients cgrnd and dgrnd (EQ1), calculated in thermosoil_coef, and 
!! called at the end of the routine to prepare for the next timestep.
!! - The effective computation of the new soil temperatures is performed in thermosoil_profile. 
!!
!! - The calling sequence of thermosoil_main is summarized in the flowchart below.
!! - Thermosoil_init and thermosoil_var_init initialize the variables from
!! restart files or with default values; they also set up
!! the vertical discretization for the numerical scheme.
!! - thermosoil_coef calculates the coefficients for the numerical scheme for the very first iteration of thermosoil;
!! after that, thermosoil_coef is called only at the end of the module to calculate the coefficients for the next timestep.
!! - thermosoil_profile solves the numerical scheme.\n
!!
!! - Flags : one unique flag : THERMOSOIL_TPRO (to be set to the desired initial soil in-depth temperature in K; by default 280K)
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): vertically discretized soil temperatures ptn, soil
!! thermal properties (pcapa, pkappa), apparent surface heat capacity (soilcap)
!! and heat flux (soilflux) to be used in enerbil at the next timestep to solve
!! the surface energy balance.
!!
!! REFERENCE(S) : 
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!!  Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin' s PhD thesis relative to the thermal
!!  integration scheme has been scanned and is provided along with the documentation, with name : 
!!  Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{thermosoil_flowchart.png}
!! \endlatexonly
!! 
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, lalo,indexgrnd, &
!Isa
       & indexnbdl, control_in, &
       ! GK
        & temp_sol_new, snow, soilcap, soilflx,  &
        & shumdiag_perma, stempdiag, ptnlev1, rest_id, hist_id, hist2_id, &
       & soiltemp,pb,grndflux,snowrho,snowdz,snowtemp,gthick,gtemp,gpkappa,&
       & pkappa_snow,cgrnd_snow,dgrnd_snow,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil,&
       & thawed_humidity, organic_layer_thick, heat_Zimov, deeptemp_prof, deephum_prof,&
       & soilc_total, veget_max)

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id,hist_id  !! Restart_ file and history file identifier 
                                                                              !! (unitless)
    INTEGER(i_std),INTENT (in)                            :: hist2_id         !! history file 2 identifier (unitless)
    REAL(r_std), INTENT (in)                              :: dtradia          !! model iteration time step in seconds (s)
    LOGICAL, INTENT(in)                                   :: ldrestart_read   !! Logical for restart files to be read 
                                                                              !! (true/false)
    LOGICAL, INTENT(in)                                   :: ldrestart_write  !! Logical for restart files to be writen 
                                                                              !! (true/false)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: index            !! Indeces of the points on the map (unitless)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)         :: lalo            !! coordinates
    INTEGER(i_std),DIMENSION (kjpindex*ngrnd), INTENT (in):: indexgrnd        !! Indeces of the points on the 3D map (vertical 
                                                                              !! dimension towards the ground) (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new     !! Surface temperature at the present time-step,
                                                                              !! Ts @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: snow             !! Snow mass @tex ($kg$) @endtex.
                                                                              !! Caveat: when there is snow on the
                                                                              !! ground, the snow is integrated into the soil for
                                                                              !! the calculation of the thermal dynamics. It means
                                                                              !! that the uppermost soil layers can completely or 
                                                                              !! partially consist in snow. In the second case, zx1
                                                                              !! and zx2 are the fraction of the soil layer 
                                                                              !! consisting in snow and 'normal' soil, respectively
                                                                              !! This is calculated in thermosoil_coef.
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)    :: shumdiag_perma         !! Relative soil humidity on the diagnostic axis 
                                                                              !! (0-1, unitless). Caveats: when "hydrol" 
                                                                              !! (the 11-layers hydrology) 
                                                                              !! is used, this humidity is 
                                                                              !! calculated with respect to the wilting point: 
                                                                              !! shumdiag_perma= (mc-mcw)/(mcs-mcw), with mc : moisture 
                                                                              !! content; mcs : saturated soil moisture content; 
                                                                              !! mcw: soil moisture content at the wilting point. 
                                                                              !! When the 2-layers hydrology "hydrolc" is used, 
                                                                              !! shumdiag_perma is just a soil wetness index, from 0 to 1
                                                                              !! but cannot direcly be linked to a soil moisture 
                                                                              !! content.
    ! specific Isa:
    INTEGER(i_std),DIMENSION (kjpindex*nbdl), INTENT (in) :: indexnbdl       !! Indeces of the points on the 3D map
    ! CR: peut-être à enlever, je ne suis pas sure de bien comprendre la gestion
    ! des flags
    TYPE(control_type), INTENT (in)                    :: control_in       !! Flags that (de)activate parts of the model
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: soilcap          !! apparent surface heat capacity
                                                                              !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: soilflx          !! apparent soil heat flux @tex ($W m^{-2}$) @endtex
                                                                              !! , positive 
                                                                              !! towards the soil, writen as Qg (ground heat flux) 
                                                                              !! in the history files, and computed at the end of 
                                                                              !! thermosoil for the calculation of Ts in enerbil, 
                                                                              !! see EQ3.
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (inout) :: stempdiag        !! diagnostic temperature profile @tex ($K$) @endtex
                                                                              !! , eg on the 
                                                                              !! diagnostic axis (levels:1:nbdl). The soil 
                                                                              !! temperature is put on this diagnostic axis to be
                                                                              !! used by other modules (slowproc.f90; routing.f90;
                                                                              !! hydrol or hydrolc when a frozen soil 
                                                                              !! parametrization is used..)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)       :: ptnlev1           !! 1st level soil temperature   

    REAL(r_std), DIMENSION(kjpindex),   INTENT (in)  :: thawed_humidity       !! specified humidity of thawed soil
    REAL(r_std), DIMENSION(kjpindex),   INTENT (inout) :: organic_layer_thick !! how deep is the organic soil?
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm), INTENT (in)   :: heat_Zimov   !! heating associated with decomposition

    REAL(r_std), DIMENSION (kjpindex,ndeep,nvm), INTENT (out) :: deephum_prof !! moisture on a deep thermodynamic profile for permafrost calcs
    REAL(r_std), DIMENSION (kjpindex,ndeep,nvm), INTENT (out):: deeptemp_prof !! temp on a deep thermodynamic profile for permafrost calcs
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (in) :: soilc_total  !! total soil carbon for use in thermal calcs
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in):: veget_max            !! Fraction of vegetation type 
    REAL(r_std), DIMENSION (kjpindex,nvm)             ::  veget_max_bg,veget_mask_real        !! Fraction of vegetation type 
    LOGICAL, SAVE                          :: ok_zimov
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: cgrnd_soil
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: dgrnd_soil
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: zdz1_soil
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: zdz2_soil
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)    :: cgrnd_snow
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)    :: dgrnd_snow

    !! 0.3 Modified variables

    !! 0.4 Local variables

    REAL(r_std),DIMENSION (kjpindex,ngrnd)                :: temp             !! buffer
    REAL(r_std),DIMENSION (kjpindex,ngrnd-1)              :: temp1            !! buffer
    REAL(r_std),DIMENSION (kjpindex)                      :: temp2            !! buffer
    CHARACTER(LEN=80)                                     :: var_name         !! To store variables names for I/O

    !##Tao
    !in thermosoil module, there is no distinguish between vegetation types.
    !Before passing the 3D snow variables to this module, snow variables should
    !first be
    !converted to 2D based on veget_max. This is done given the snowpack.f90 is
    !written for 3D snow variables.  
    REAL(r_std),DIMENSION (kjpindex,ngrnd),INTENT(inout)  :: soiltemp ! soil temperature profile
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)               :: pb !surface pressure 
    REAL(r_std),DIMENSION (kjpindex),INTENT(inout)               :: grndflux 
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowrho,snowdz,snowtemp
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)           :: gthick !the first soil layer thickness
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)           :: gtemp  ! the first soil layer temperature
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)           :: gpkappa ! the first soil layer thermal conductivity
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT(in)       :: pkappa_snow
    REAL(r_std),DIMENSION (kjpindex,ngrnd) :: pkappa_pftmean           
    INTEGER(i_std)   :: jv,ji,ii,m,jg
    CHARACTER(LEN=10) :: part_str    ! string suffix indicating an index
    !!Tao
    !
!_ ================================================================================================================================

  !! 1. do initialisation

  IF (l_first_thermosoil) THEN

        veget_max_bg(:,2:nvm) = veget_max(:,2:nvm)
        veget_max_bg(:,1) = MAX((un - SUM(veget_max(:,2:nvm), 2)), zero)
        IF (long_print) WRITE (numout,*) ' l_first_thermosoil : call thermosoil_init '
        CALL getin_p('OK_ZIMOV',ok_zimov)

        
        !! 1.1. Allocate and initialize soil temperatures variables
        !! by reading restart files or using default values.
        CALL thermosoil_init (kjit, ldrestart_read, kjpindex, index, lalo, rest_id, &
                             & snowdz)

        !! 1.2.Computes physical constants and arrays; initializes soil thermal properties; produces the first stempdiag
        !!  Computes some physical constants and arrays depending on the soil vertical discretization 
        !! (lskin, cstgrnd, zz, zz_coef, dz1, dz2); get the vertical humidity onto the thermal levels, and 
        !! initializes soil thermal properties (pkappa, pcapa); produces the first temperature diagnostic stempdiag.

        
        CALL thermosoil_var_init (kjpindex, zz, zz_coef, dz1, dz2, pkappa, pcapa, pcapa_en, &
        &        shumdiag_perma, stempdiag, snow, &
!Isa
        & profil_froz,pb,snowtemp,snowrho,snowdz,&
        & thawed_humidity,organic_layer_thick, soilc_total, veget_max_bg)
        !
        !! 1.3. Computes cgrd, dgrd, soilflx and soilcap coefficients from restart values or initialisation values.
        ! computes cgrd and dgrd coefficient from previous time step (restart)
        !
        !WRITE(numout,*) 'zd _coef IF(l_first_) 1 ','ptn(1,:,10)',ptn(1,:,10)

        CALL thermosoil_coef (kjpindex, dtradia, temp_sol_new, snow, ptn, soilcap, soilflx, zz, dz1, dz2, z1, zdz1,&
           & zdz2, cgrnd, dgrnd, zdz1_soil, zdz2_soil, cgrnd_soil, dgrnd_soil, pcapa, pcapa_en, pkappa,&
!Isa
           & profil_froz, pcappa_supp, stempdiag, &
           & organic_layer_thick, soilc_total,veget_max_bg,snowdz)
        !WRITE(numout,*) 'zd _coef IF(l_first_) 2 ','ptn(1,:,10)',ptn(1,:,10)
        
        
        !! 1.4. call to thermosoil_energy, if you wish to perform some checks (?)
        !!?? the usefulness of this routine seems questionable.
        CALL thermosoil_energy (kjpindex, temp_sol_new, soilcap, .TRUE., veget_max_bg) 

    
        ! ajout CR from trunk
        !! 1.5. read restart files for other variables than ptn.
        !!?? mind the use of ok_var here.
        !!?? ok_var is a function of sechiba_io_p.f90, documented as follows :
        !!!! pour déclancher les restarts rajoutés avec un paramètre externe
           !!FUNCTION ok_var ( varname )
           !!CHARACTER(LEN=*), INTENT(IN) :: varname
           !!LOGICAL ok_var
           !!ok_var=.FALSE.
           !!CALL getin_p(varname, ok_var)
           !!END FUNCTION ok_var
         !!
         !! from what we understand, it looks for the chain varname in
         !!run.def; if absent, returns .FALSE., and the variable named
         !!'varname' is not searched for in the restart. This looks like a
         !!trick to read variables in restart files when they are not read
         !!there by default. For all variables in the following sequence, ok_var
         !!is by default false, so don' t bother about this.
         !! this is also logical as those variables have been initialized
         !!above.
         !!?? so maybe this part of the code could be deleted to add clarity.
        
        
!    if (control%ok_converge_isaorig) then
!             ! il n'y avait pas d'init???
!         else ! if (ok_converge_isaorig) then
         IF (ldrestart_read) THEN
           IF (long_print) WRITE (numout,*) ' we have to READ a restart file for THERMOSOIL variables'

           var_name= 'cgrnd'
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','Cgrnd coefficient.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, .TRUE., temp1, "gather", nbp_glo, index_g)
              IF (MINVAL(temp1) < MAXVAL(temp1) .OR. MAXVAL(temp1) .NE. val_exp) THEN
                 DO m=1,nvm
                    cgrnd(:,:,m)=temp1(:,:)
                 END DO
              ENDIF
           ENDIF

           var_name= 'dgrnd'
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','Dgrnd coefficient.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, .TRUE., temp1, "gather", nbp_glo, index_g)
              IF (MINVAL(temp1) < MAXVAL(temp1) .OR. MAXVAL(temp1) .NE. val_exp) THEN
                DO m=1,nvm
                    dgrnd(:,:,m)=temp1(:,:)
                END DO
              ENDIF
           ENDIF

           var_name= 'z1'
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp2, "gather", nbp_glo, index_g)
              IF (MINVAL(temp2) < MAXVAL(temp2) .OR. MAXVAL(temp2) .NE. val_exp) THEN
                 z1(:)=temp2(:)
              ENDIF
           ENDIF

           var_name= 'pcapa'
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 DO m=1,nvm
                    pcapa(:,:,m)=temp(:,:)
                 END DO
              ENDIF
           ENDIF

           var_name= 'pcapa_en'
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 DO m=1,nvm
                    pcapa_en(:,:,m)=temp(:,:)
                 END DO
              ENDIF
           ENDIF

           var_name= 'pkappa'
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 DO m=1,nvm
                    pkappa(:,:,m)=temp(:,:)
                 END DO
              ENDIF
           ENDIF

           var_name= 'zdz1'
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, .TRUE., temp1, "gather", nbp_glo, index_g)
              IF (MINVAL(temp1) < MAXVAL(temp1) .OR. MAXVAL(temp1) .NE. val_exp) THEN
                 DO m=1,nvm
                    zdz1(:,:,m)=temp1(:,:)
                 END DO
              ENDIF
           ENDIF

           var_name= 'zdz2'
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','?.')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., temp, "gather", nbp_glo, index_g)
              IF (MINVAL(temp) < MAXVAL(temp) .OR. MAXVAL(temp) .NE. val_exp) THEN
                 DO m=1,nvm
                    zdz2(:,:,m)=temp(:,:)
                 END DO
              ENDIF
           ENDIF

           var_name='temp_sol_beg'
           CALL ioconf_setatt_p('UNITS', 'K')
           CALL ioconf_setatt_p('LONG_NAME','Old Surface temperature')
           IF ( ok_var(var_name) ) THEN
              CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., temp2, "gather", nbp_glo, index_g)
              IF (MINVAL(temp2) < MAXVAL(temp2) .OR. MAXVAL(temp2) .NE. val_exp) THEN
                 temp_sol_beg(:) = temp2(:)
              ENDIF
           ENDIF

        ENDIF !ldrestart_read

        RETURN 

    ENDIF !l_first_thermosoil

    !     make vegetation masks so that we don't bother to calculated pfts on
    !     gridcells where they don's exist
    veget_max_bg(:,2:nvm) = veget_max(:,2:nvm)
    veget_max_bg(:,1) = MAX((un - SUM(veget_max(:,2:nvm), 2)), zero)
    veget_mask_2d(:,:) = .TRUE.
    !veget_mask_2d(:,:) = ( veget_max_bg(:,:) .GT. EPSILON(0.))
    WHERE ( veget_mask_2d(:,:) )
       veget_mask_real(:,:) = un
    ELSEWHERE
       veget_mask_real(:,:) = zero
    END WHERE

    
  !! 2. Prepares the restart files for the next simulation

    !!?? do all the coefficients (cgrnd, dgrnd...) be put in the restart file
    !! as they are by default not read there, but calculated in
    !!thermosoil_var_init from the restart or initial temperature ?
    !! exceptions are soilcap and soilflx, used in enerbil, and of course ptn.
    IF (ldrestart_write) THEN

        IF (long_print) WRITE (numout,*) ' we have to complete restart file with THERMOSOIL variables'

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'ptn_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, ptn(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        if (ok_wetdiaglong) then
            DO m=1,nvm
               WRITE(part_str,'(I2)') m
               IF (m < 10) part_str(1:1) = '0'
               var_name = 'wetdiaglong_'//part_str(1:LEN_TRIM(part_str))
               CALL restput_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit,wetdiaglong(:,:,m), 'scatter', nbp_glo, index_g) !need to add veg dim    
            END DO
        end if

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'wetdiag_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, wetdiag(:,:,m), 'scatter', nbp_glo, index_g)      !need to add veg dim
        END DO

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'E_lat_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, E_sol_lat_couche(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        !!var_name= 'E_sol_lat_couche'
        !!CALL restput_p(rest_id, var_name, nbp_glo,ngrnd , 1, kjit, E_sol_lat_couche, 'scatter', nbp_glo, index_g)

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'cgrnd_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, cgrnd(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'dgrnd_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, dgrnd(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        var_name= 'z1'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, z1, 'scatter', nbp_glo, index_g)

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'pcapa_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, pcapa(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'pcapa_en_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, pcapa_en(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'pkappa_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, pkappa(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'zdz1_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd-1, 1, kjit, zdz1(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF (m < 10) part_str(1:1) = '0'
           var_name = 'zdz2_'//part_str(1:LEN_TRIM(part_str))
           CALL restput_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, zdz2(:,:,m), 'scatter', nbp_glo, index_g)
        END DO

        var_name= 'temp_sol_beg'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, temp_sol_beg, 'scatter', nbp_glo, index_g)

        var_name= 'soilcap'  
        CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  soilcap, 'scatter',  nbp_glo, index_g)
        
        var_name= 'soilflx'  
        CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  soilflx, 'scatter',  nbp_glo, index_g)

        ! read in enerbil
        var_name= 'temp_sol_new'
        CALL restput_p(rest_id, var_name, nbp_glo, 1, 1, kjit, temp_sol_new, 'scatter', nbp_glo, index_g)

        RETURN 

    END IF !ldrestart_write


  !! 3. Put the soil wetness diagnostic on the levels of the soil temperature

    !!?? this could logically be put just before the last call to
    !!thermosoil_coef, as the results are used there...    !
    !##Tao
    snow_thick(:) = 0.0
    DO ji=1,kjpindex
       snow_thick(ji) = SUM(snowdz(ji,:))
    ENDDO
    !!Tao

    CALL thermosoil_humlev(kjpindex, shumdiag_perma, thawed_humidity)
    
    ! Compute long-term soil humidity (for permafrost)
    !    

    if (ok_wetdiaglong) then
        CALL thermosoil_wlupdate( kjpindex, dtradia, ptn, wetdiag, wetdiaglong )
    else
        wetdiaglong(:,:,:)=wetdiag(:,:,:)
    endif


  !! 4. Effective computation of the soil temperatures profile, using the cgrd and !dgrd coefficients from previsou tstep.
    !WRITE(numout,*) 'zd _profile 1 ','ptn(1,:,10)',ptn(1,:,10),'cgrnd_snow(1,:)',cgrnd_snow(1,:),'dgrnd_snow(1,:)',dgrnd_snow(1,:),'snowtemp(1,:)',snowtemp(1,:),'temp_sol_new(1)',temp_sol_new(1),'cgrnd(1,:,10)',cgrnd(1,:,10),'dgrnd(1,:,10)',dgrnd(1,:,10)
    CALL thermosoil_profile (kjpindex, temp_sol_new, ptn, &
                &stempdiag,pkappa_snow,snowdz,snowtemp,grndflux,dtradia,veget_max,cgrnd_snow,dgrnd_snow)
    !WRITE(numout,*) 'zd _profile 2 ','ptn(1,:,10)',ptn(1,:,10)
! on appelle thermosoil_profile comme dans le trunk, car ça ne change rien

  !! 5. Call to thermosoil_energy, still to be clarified..

    CALL thermosoil_energy (kjpindex, temp_sol_new, soilcap, .FALSE.,veget_max_bg)
    ptn_pftmean(:,:) = zero
    do m=1,nvm
       do jg = 1, ngrnd
          ptn_pftmean(:,jg) = ptn_pftmean(:,jg) + ptn(:,jg,m) * veget_max_bg(:,m)
       end do
    end do
    soiltemp=ptn_pftmean
   
    !! make pft-mean pcappa_supp and E_sol_lat_couche
    pcappa_supp_pftmean(:,:) = zero
    E_sol_lat_couche_pftmean(:,:) = zero
    do m=1,nvm
       do jg = 1, ngrnd
          pcappa_supp_pftmean(:,jg) = pcappa_supp_pftmean(:,jg) + pcappa_supp(:,jg,m) * veget_max_bg(:,m)
       end do
    end do
    do m=1,nvm
       do jg = 1, ngrnd
          E_sol_lat_couche_pftmean(:,jg) = E_sol_lat_couche_pftmean(:,jg) + E_sol_lat_couche(:,jg,m) * veget_max_bg(:,m)
       end do
    end do
 
    IF ( .NOT. almaoutput ) THEN
       !!need to write with PFT dimension
       DO jv = 1, nvm
          WRITE(part_str,'(I2)') jv
          IF (jv < 10) part_str(1:1) = '0'
          CALL histwrite_p(hist_id, 'ptn_'//part_str(1:LEN_TRIM(part_str)), &
               kjit, ptn(:,:,jv), kjpindex*ngrnd, indexgrnd)
       END DO
       CALL histwrite_p(hist_id, 'ptn_pftmean', kjit, ptn_pftmean, kjpindex*ngrnd, indexgrnd)
       IF (control_in%hydrol_cwrr) THEN
         DO jv = 1, nvm
             WRITE(part_str,'(I2)') jv
             IF (jv < 10) part_str(1:1) = '0'
             CALL histwrite_p(hist_id, 'pcapa_'//part_str(1:LEN_TRIM(part_str)), &
                  kjit, pcapa(:,:,jv), kjpindex*ngrnd, indexgrnd)
             !CALL histwrite_p(hist_id, 'pcappa_supp_'//part_str(1:LEN_TRIM(part_str)), &
             !     kjit, pcappa_supp(:,:,jv), kjpindex*ngrnd, indexgrnd)
             CALL histwrite_p(hist_id, 'pkappa_'//part_str(1:LEN_TRIM(part_str)), &
                  kjit, pkappa(:,:,jv), kjpindex*ngrnd, indexgrnd)
             CALL histwrite_p(hist_id, 'wetdiag_'//part_str(1:LEN_TRIM(part_str)), &
                  kjit, wetdiag(:,:,jv), kjpindex*ngrnd, indexgrnd)
             CALL histwrite_p(hist_id,'wetdiaglong_'//part_str(1:LEN_TRIM(part_str)), &
                  kjit, wetdiaglong(:,:,jv), kjpindex*ngrnd, indexgrnd)
             !CALL histwrite_p(hist_id,'ptn_beg_'//part_str(1:LEN_TRIM(part_str)), &
             !     kjit, ptn_beg(:,:,jv), kjpindex*ngrnd, indexgrnd)
             !CALL histwrite_p(hist_id,'profil_froz_'//part_str(1:LEN_TRIM(part_str)), &
             !     kjit, profil_froz(:,:,jv), kjpindex*ngrnd, indexgrnd)
         END DO
         !CALL histwrite_p(hist_id, 'shumdiag_perma', kjit, shumdiag_perma, kjpindex*nbdl, indexnbdl)
         CALL histwrite_p(hist_id, 'stempdiag', kjit, stempdiag, kjpindex*nbdl,indexnbdl)
       END IF
      CALL histwrite_p(hist_id, 'Qg', kjit, soilflx, kjpindex, index)

    ELSE !IF ( .NOT. almaoutput ) THEN
      CALL histwrite_p(hist_id, 'SoilTemp', kjit, ptn_pftmean, kjpindex*ngrnd, indexgrnd)
      !CALL histwrite_p(hist_id, 'SoilTemp', kjit, ptn, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'Qg', kjit, soilflx, kjpindex, index)
      CALL histwrite_p(hist_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
      CALL histwrite_p(hist_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
    ENDIF  !IF ( .NOT. almaoutput ) THEN
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
          CALL histwrite_p(hist_id, 'SoilTemp', kjit, ptn_pftmean, kjpindex*ngrnd, indexgrnd)
          !CALL histwrite_p(hist2_id, 'ptn', kjit, ptn, kjpindex*ngrnd, indexgrnd)
       ELSE
          CALL histwrite_p(hist2_id, 'SoilTemp', kjit, ptn_pftmean, kjpindex*ngrnd, indexgrnd)
          !CALL histwrite_p(hist2_id, 'SoilTemp', kjit, ptn, kjpindex*ngrnd, indexgrnd)
          CALL histwrite_p(hist2_id, 'Qg', kjit, soilflx, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
       ENDIF
    ENDIF
   
  !! 7. Considering the heat released by microbial respiration
    IF (ok_zimov) THEN
       CALL add_heat_Zimov(kjpindex, dtradia, ptn, heat_zimov)
    END IF

  !! 8. A last final call to thermosoil_coef
 
    !! A last final call to thermosoil_coef, which calculates the different
    !!coefficients (cgrnd, dgrnd, dz1, z1, zdz2, soilcap, soilflx) from this time step to be
    !!used at the next time step, either in the surface temperature calculation
    !!(soilcap, soilflx) or in the soil thermal numerical scheme.
    !
    !Isa
! computes cgrd and dgrd coefficient
    !WRITE(numout,*) 'zd _coef 1 ','ptn(1,:,10)',ptn(1,:,10)
    CALL thermosoil_coef (kjpindex, dtradia, temp_sol_new, snow, ptn,soilcap, soilflx, zz,&
    & dz1, dz2, z1, zdz1,zdz2, cgrnd, dgrnd, zdz1_soil, zdz2_soil, cgrnd_soil, dgrnd_soil,&
    & pcapa, pcapa_en, pkappa, profil_froz, pcappa_supp, stempdiag,&
    & organic_layer_thick, soilc_total,veget_max_bg,snowdz)
    !WRITE(numout,*) 'zd _coef 2 ','ptn(1,:,10)',ptn(1,:,10)

    !save some useful variables for new snow model
    ptn_pftmean(:,:) = zero
    pkappa_pftmean(:,:) = zero
    do m=1,nvm
       do jg = 1, ngrnd
          ptn_pftmean(:,jg) = ptn_pftmean(:,jg) + ptn(:,jg,m) * veget_max_bg(:,m)
          pkappa_pftmean(:,jg) = pkappa_pftmean(:,jg) + pkappa(:,jg,m) * veget_max_bg(:,m)
       end do
    end do

    DO ji=1,kjpindex
       gthick(ji) = zz_coef(1)
       gtemp(ji) = ptn_pftmean(ji,1)
       gpkappa(ji)=pkappa_pftmean(ji,1)
       grndflux(ji)=0.0
    ENDDO
    !!Tao
           ! CR: ajout du trunk
    ptnlev1(:) = ptn_pftmean(:,1)

    !++cdk prep updated temp and moisture fields so they can be sent to stomate
    !permafrost calcs
    deephum_prof = wetdiaglong
    deeptemp_prof = ptn
    !--cdk

 !   ptn_beg(:,:)      = ptn(:,:) ! seulement dans isa: 
 ! pas necessaire dans trunk car deja fait dans thermosoil_energy
    IF (long_print) WRITE (numout,*) ' thermosoil_main done '

  END SUBROUTINE thermosoil_main

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_init
!!
!>\BRIEF        Allocates local and global arrays; initializes soil temperatures using either restart files
!! or a fixed value set by the flag THERMOSOIL_TPRO.
!!		  
!! DESCRIPTION  : flag : THERMOSOIL_TPRO (to be set to the desired initial temperature in K; by default 280K).
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

  SUBROUTINE thermosoil_init(kjit, ldrestart_read, kjpindex, index, lalo, rest_id, &
                             & snowdz)
    !! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number (unitless) 
    LOGICAL,INTENT (in)                                 :: ldrestart_read     !! Logical for restart file to read (true/false)
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map (unitless)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT(in)      :: lalo               !! coordinates
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! Restart file identifier (unitless)
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    REAL(r_std),DIMENSION (kjpindex,ngrnd,nvm)         :: reftemp_3d
    INTEGER(i_std)                                     :: ier, i, m, jv
    !##Tao
    INTEGER(i_std)  :: ji
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowdz
    !!Tao
    CHARACTER(LEN=80)                                  :: var_name            !! To store variables names for I/O
    CHARACTER(LEN=10) :: part_str    ! string suffix indicating an index

!_ ================================================================================================================================

  !! 1. Initialisation

    !! Initialisation has to be done only one time, so the logical
    !! logical l_first_thermosoil has to be set to .FALSE. now..
    IF (l_first_thermosoil) THEN 
        l_first_thermosoil=.FALSE.
    ELSE 
        WRITE (numout,*) ' l_first_thermosoil false . we stop '
        STOP 'thermosoil_init'
    ENDIF

     ok_wetdiaglong = .FALSE.
     CALL getin_p ('OK_WETDIAGLONG',ok_wetdiaglong)
    if (ok_freeze_thermix .AND. ok_pc) then
        ok_wetdiaglong = .TRUE.
    endif

    CALL getin_p('satsoil', satsoil)
    if (ok_freeze_thermix .AND. ok_pc) then
        use_toporganiclayer_tempdiff = .false.
        CALL getin_p('USE_TOPORGANICLAYER_TEMPDIFF',use_toporganiclayer_tempdiff)

        use_soilc_tempdiff = .false.
        CALL getin_p('USE_SOILC_TEMPDIFF', use_soilc_tempdiff)
        IF (use_toporganiclayer_tempdiff .AND. use_soilc_tempdiff) THEN
           WRITE(*,*) 'warning: thermosoil_getdiff: cant have both use_toporganiclayer_tempdiff and'
           WRITE(*,*) 'use_soilc_tempdiff set to .true.. using only use_soilc_tempdiff.'
           use_toporganiclayer_tempdiff = .FALSE.
        ENDIF
    endif

  !! 2. Arrays allocations

    ALLOCATE (ptn(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ptn allocation. We stop. We need ',kjpindex,' fois ',ngrnd,' words = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (ptn_pftmean(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ptn_pftmean allocation. We stop. We need',kjpindex,' fois ',ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (zz(ngrnd),stat=ier)
    IF (ier /= 0) THEN
       CALL ipslerr_p(3,'thermosoil_init', 'Error in allocation of zz','','')
    END IF

    ALLOCATE (zz_coef(ngrnd),stat=ier)
    IF (ier /= 0) THEN
       CALL ipslerr_p(3,'thermosoil_init', 'Error in allocation of zz_coef','','')
    END IF

    ALLOCATE (dz1(ngrnd),stat=ier)
    IF (ier /= 0) THEN
       CALL ipslerr_p(3,'thermosoil_init', 'Error in allocation of dz1','','')
    END IF

    ALLOCATE (dz2(ngrnd),stat=ier)
    IF (ier /= 0) THEN
       CALL ipslerr_p(3,'thermosoil_init', 'Error in allocation of dz2','','')
    END IF


    ALLOCATE (z1(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in z1 allocation. We STOP. We need ',kjpindex,' words '
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (cgrnd(kjpindex,ngrnd-1,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in cgrnd allocation. We STOP. We need ',kjpindex,' fois ',ngrnd-1 ,' words  = '&
           & , kjpindex*(ngrnd-1)
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (dgrnd(kjpindex,ngrnd-1,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in dgrnd allocation. We STOP. We need ',kjpindex,' fois ',ngrnd-1 ,' words  = '&
           & , kjpindex*(ngrnd-1)
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (pcapa(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pcapa allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (pkappa(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pkappa allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (zdz1(kjpindex,ngrnd-1,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zdz1 allocation. We STOP. We need ',kjpindex,' fois ',ngrnd-1 ,' words  = '&
           & , kjpindex*(ngrnd-1)
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (zdz2(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in zdz2 allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (surfheat_incr(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in surfheat_incr allocation. We STOP. We need ',kjpindex,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (coldcont_incr(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in coldcont_incr allocation. We STOP. We need ',kjpindex,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (pcapa_en(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pcapa_en allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (ptn_beg(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ptn_beg allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (temp_sol_beg(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in temp_sol_beg allocation. We STOP. We need ',kjpindex,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (wetdiag(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in wetdiag allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    wetdiag(:,:,:)=val_exp
    ALLOCATE (wetdiaglong(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in wetdiaglong allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

!##Tao
    ALLOCATE (snow_thick(kjpindex),stat=ier);
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in snow_thick allocation. We STOP. We need',kjpindex,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF

!!Tao

!Isa ++
    ALLOCATE (profil_froz(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in profil_froz allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    
     ALLOCATE (pcappa_supp(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pcapa_supp allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF
    
     ALLOCATE (E_sol_lat_couche(kjpindex,ngrnd,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in E_sol_sat allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

     ALLOCATE (pcappa_supp_pftmean(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in pcapa_supp allocation. We STOP. We need',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

     ALLOCATE (E_sol_lat_couche_pftmean(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in E_sol_sat allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex*ngrnd
        STOP 'thermosoil_init'
    END IF

    ALLOCATE (permafrost(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in permafrost allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF
    ALLOCATE (excess_ice(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in excess_ice allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF
    ALLOCATE (overburden(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in overburden allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF
    ALLOCATE (reftemp(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in reftemp allocation. We STOP. We need ',kjpindex,' fois ',ngrnd ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF
    ALLOCATE (veget_mask_2d(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in reftemp allocation. We STOP. We need',kjpindex,' fois ',nvm ,' words  = '&
           & , kjpindex
        STOP 'thermosoil_init'
    END IF

    !Isa --

    !
  !! 3. Reads restart files for soil temperatures only 
    
    !! Reads restart files for soil temperatures only. If no restart file is
    !! found,  the initial soil temperature is by default set to 280K at all depths. The user
    !! can decide to initialize soil temperatures at an other value, in which case he should set the flag THERMOSOIL_TPRO
    !! to this specific value in the run.def.

    IF (ldrestart_read) THEN
        IF (long_print) WRITE (numout,*) ' we have to READ a restart file for THERMOSOIL variables'

        ptn(:,:,:) = val_exp
        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF ( m < 10 ) part_str(1:1) = '0'
           var_name = 'ptn_'//part_str(1:LEN_TRIM(part_str))
           CALL ioconf_setatt_p('UNITS', 'K')
           CALL ioconf_setatt_p('LONG_NAME','Soil Temperature profile')
           CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE., ptn(:,:,m), "gather", nbp_glo, index_g) !need to add veg dim
        END DO

        !Config Key   = THERMOSOIL_TPRO
        !Config Desc  = Initial soil temperature profile if not found in restart
        !Config Def   = 280.
        !Config If    = OK_SECHIBA
        !Config Help  = The initial value of the temperature profile in the soil if 
        !Config         its value is not found in the restart file. This should only 
        !Config         be used if the model is started without a restart file. Here
        !Config         we only require one value as we will assume a constant 
        !Config         throughout the column.
        !Config Units = Kelvin [K]
!
	if (read_permafrost_map) then 
		CALL read_permafrostmap(kjpindex,lalo,overburden,excess_ice,permafrost)
		DO jv=1,nvm
                  DO i=1, ngrnd
                     WHERE (veget_mask_2d(:,jv))
		       where (ptn(:,i,jv).eq.val_exp.AND.permafrost(:).eq.1)
		             ptn(:,i,jv)=272._r_std
		       elsewhere 
		             ptn(:,i,jv)=280._r_std
		       endwhere
                     ENDWHERE
		  ENDDO
                ENDDO
	else if (read_reftemp) then
		CALL read_reftempfile(kjpindex,lalo,reftemp)
                DO jv = 1,nvm
                  reftemp_3d(:,:,jv)=reftemp(:,:)
                ENDDO
                ptn(:,:,:) = reftemp_3d(:,:,:)
       		!CALL setvar_p (ptn, val_exp, 'NO_KEYWORD' ,reftemp_3d)
	else
       		CALL setvar_p (ptn, val_exp,'THERMOSOIL_TPRO',272._r_std)
	endif
 
        ptn_beg(:,:,:) = ptn(:,:,:)
     


!! Arsene 15-06-2015 Remove (see below)
!        wetdiaglong(:,:,:) = val_exp
!        DO m=1,nvm
!           WRITE(part_str,'(I2)') m
!           IF ( m < 10 ) part_str(1:1) = '0'
!           var_name = 'wetdiaglong_'//part_str(1:LEN_TRIM(part_str))
!           CALL ioconf_setatt_p('UNITS', '-')
!           CALL ioconf_setatt_p('LONG_NAME','Long-term soil humidity')
!           CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE.,wetdiaglong(:,:,m), "gather", nbp_glo, index_g) !need to add veg dim
!        END DO
!! Arsene 15-06-2015 Remove (see below)


        wetdiag(:,:,:) = val_exp
        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF ( m < 10 ) part_str(1:1) = '0'
           var_name = 'wetdiag_'//part_str(1:LEN_TRIM(part_str))
           CALL ioconf_setatt_p('UNITS', '-')
           CALL ioconf_setatt_p('LONG_NAME','soil humidity')
           CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE.,wetdiag(:,:,m), "gather", nbp_glo, index_g) !need to add veg dim 
        END DO


        IF (ok_wetdiaglong) THEN   !! Arsene15-06-2015 Add... because if is not wetdiaglong is not def (and  wetdiaglong(:,:,:) = val_exp)
           wetdiaglong(:,:,:) = val_exp
           DO m=1,nvm
              WRITE(part_str,'(I2)') m
              IF ( m < 10 ) part_str(1:1) = '0'
              var_name = 'wetdiaglong_'//part_str(1:LEN_TRIM(part_str))
              CALL ioconf_setatt_p('UNITS', '-')
              CALL ioconf_setatt_p('LONG_NAME','Long-term soil humidity')
              CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit, .TRUE.,wetdiaglong(:,:,m), "gather", nbp_glo, index_g) !need to add veg dim
           END DO
        ELSE                                     !! Arsene 15-06-2015 Add because the variable is not in the restart file if .NOT.
            wetdiaglong(:,:,:)=wetdiag(:,:,:)    !! Arsene 15-06-2015 Add because the variable is not in the restart file if .NOT.
        ENDIF                                    !! Arsene 15-06-2015 Add because the variable is not in the restart file if .NOT.


        IF ( ALL(ABS(wetdiag(:,:,:)-val_exp).LT.EPSILON(val_exp)) ) THEN
           wetdiag = 1.
        ENDIF
        IF ( ALL(ABS(wetdiaglong(:,:,:)-val_exp).LT.EPSILON(val_exp)) ) THEN
           wetdiaglong = 1.
        ENDIF

        E_sol_lat_couche(:,:,:) = val_exp
        DO m=1,nvm
           WRITE(part_str,'(I2)') m
           IF ( m < 10 ) part_str(1:1) = '0'
           var_name = 'E_lat_'//part_str(1:LEN_TRIM(part_str))
           CALL ioconf_setatt_p('UNITS', 'J/m2/layer')
           CALL ioconf_setatt_p('LONG_NAME','Latent Energy')
           CALL restget_p (rest_id, var_name, nbp_glo, ngrnd, 1, kjit,.TRUE.,E_sol_lat_couche(:,:,m), "gather", nbp_glo, index_g)
       ENDDO

       CALL setvar_p (E_sol_lat_couche, val_exp,'NO_KEYWORD',zero)

   ENDIF

    !##Tao 
    ! init all the other snow variables to zero so that we can see that they
    ! aren't used
    snow_thick(:) = 0.0
    DO ji=1,kjpindex
        snow_thick(ji) =  SUM(snowdz(ji,:))
    ENDDO
    !!Tao

    IF (long_print) WRITE (numout,*) ' thermosoil_init done '


  END SUBROUTINE thermosoil_init

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_clear
!!
!>\BRIEF        Sets the flag l_first_thermosoil to true and desallocates the allocated arrays.
!! ??!! the call of thermosoil_clear originates from sechiba_clear but the calling sequence and 
!! its purpose require further investigation.
!!
!! DESCRIPTION  : None
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
 SUBROUTINE thermosoil_clear()

        l_first_thermosoil=.TRUE.
 
        IF ( ALLOCATED (ptn)) DEALLOCATE (ptn)
        IF ( ALLOCATED (ptn_pftmean)) DEALLOCATE (ptn_pftmean)
        IF ( ALLOCATED (z1)) DEALLOCATE (z1) 
        IF ( ALLOCATED (cgrnd)) DEALLOCATE (cgrnd) 
        IF ( ALLOCATED (dgrnd)) DEALLOCATE (dgrnd) 
        IF ( ALLOCATED (pcapa)) DEALLOCATE (pcapa)
        IF ( ALLOCATED (pkappa))  DEALLOCATE (pkappa)
        IF ( ALLOCATED (zdz1)) DEALLOCATE (zdz1)
        IF ( ALLOCATED (zdz2)) DEALLOCATE (zdz2)
        IF ( ALLOCATED (pcapa_en)) DEALLOCATE (pcapa_en)
        IF ( ALLOCATED (ptn_beg)) DEALLOCATE (ptn_beg)
        IF ( ALLOCATED (temp_sol_beg)) DEALLOCATE (temp_sol_beg)
        IF ( ALLOCATED (surfheat_incr)) DEALLOCATE (surfheat_incr)
        IF ( ALLOCATED (coldcont_incr)) DEALLOCATE (coldcont_incr)
        IF ( ALLOCATED (wetdiag)) DEALLOCATE (wetdiag)
!Isa
        IF ( ALLOCATED (profil_froz)) DEALLOCATE (profil_froz)
        IF ( ALLOCATED (wetdiaglong)) DEALLOCATE (wetdiaglong)
        !##Tao
        IF ( ALLOCATED (snow_thick)) DEALLOCATE (snow_thick)
        !!Tao

  END SUBROUTINE thermosoil_clear
  !!
  !!
  FUNCTION fz(rk) RESULT (fz_result)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    REAL(r_std), INTENT(in)                        :: rk
    
    !! 0.2 Output variables

    REAL(r_std)                                    :: fz_result
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

!_ ================================================================================================================================

    fz_result = fz1 * (zalph ** rk - un) / (zalph - un)

  END FUNCTION fz


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_var_init
!!
!>\BRIEF        Define and initializes the soil thermal parameters
!!		  
!! DESCRIPTION	: This routine\n
!! 1. Defines the parameters ruling the vertical grid of the thermal scheme (fz1, zalpha).\n
!! 2. Defines the scaling coefficients for adimensional depths (lskin, cstgrnd, see explanation in the 
!!    variables description of thermosoil_main). \n
!! 3. Calculates the vertical discretization of the soil (zz, zz_coef, dz2) and the constants used
!!    in the numerical scheme and which depend only on the discretization (dz1, lambda).\n
!! 4. Initializes the soil thermal parameters (capacity, conductivity) based on initial soil moisture content.\n
!! 5. Produces a first temperature diagnostic based on temperature initialization.\n
!!
!! The scheme comprizes ngrnd=7 layers by default.
!! The layer' s boundaries depths (zz_coef) follow a geometric series of ratio zalph=2 and first term fz1.\n
!! zz_coef(jg)=fz1.(1-zalph^jg)/(1-zalph) \n
!! The layers' boudaries depths are first calculated 'adimensionally', ie with a
!! discretization adapted to EQ5. This discretization is chosen for its ability at
!! reproducing a thermal signal with periods ranging from days to centuries. (see
!! Hourdin, 1992). Typically, fz1 is chosen as : fz1=0.3*cstgrnd (with cstgrnd the
!! adimensional attenuation depth). \n
!! The factor lskin/cstgrnd is then used to go from adimensional depths to
!! depths in m.\n
!! zz(real)=lskin/cstgrnd*zz(adimensional)\n
!! Similarly, the depths of the numerical nodes are first calculated
!! adimensionally, then the conversion factor is applied.\n
!! the numerical nodes (zz) are not exactly the layers' centers : their depths are calculated as follows:\n
!! zz(jg)=fz1.(1-zalph^(jg-1/2))/(1-zalph)\n
!! The values of zz and zz_coef used in the default thermal discretization are in the following table.
!! \latexonly
!! \includegraphics{thermosoil_var_init1.jpg}
!! \endlatexonly\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : None
!!
!! REFERENCE(S)	:
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of 
!! planetary atmospheres, Ph.D. thesis, Paris VII University.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_var_init(kjpindex, zz, zz_coef, dz1, dz2, pkappa, pcapa, pcapa_en, &
  &     shumdiag_perma, stempdiag,&
 !Isa
  & snow, profil_froz,pb,snowtemp,snowrho,snowdz, &
  & thawed_humidity,organic_layer_thick, soilc_total, veget_max)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size (unitless)
    REAL(r_std), DIMENSION (kjpindex,nbdl), INTENT (in)      :: shumdiag_perma    !! Relative soil humidity on the diagnostic axis 
                                                                                  !! (unitless), [0,1]. (see description of the 
                                                                                  !! variables of thermosoil_main for more 
                                                                                  !! explanations) 
    
    !! 0.2 Output variables

    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: zz                !! depths of the layers'numerical nodes 
                                                                                  !! @tex ($m$)@endtex 
    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: zz_coef		  !! depths of the layers'boundaries 
                                                                                  !! @tex ($m$)@endtex 
    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: dz1               !! numerical constant depending on the vertical
                                                                                  !! thermal grid only @tex  ($m^{-1}$) @endtex. 
                                                                                  !! (see description
                                                                                  !! of the variables of thermosoil_main for more
                                                                                  !! explanations)
    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: dz2               !! thicknesses of the soil thermal layers 
                                                                                  !! @tex ($m$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(out)     :: pcapa             !! volumetric vertically discretized soil heat 
                                                                                  !! capacity @tex ($J K^{-1} m^{-3}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(out)     :: pcapa_en	  !! volumetric vertically discretized heat 
                                                                                  !! capacity used in thermosoil_energy
                                                                                  !! @tex ($J K^{-1} m^{-3}$) @endtex ;
                                                                                  !! usefulness still to be clarified.
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(out)     :: pkappa        !! vertically discretized soil thermal 
                                                                                  !! conductivity @tex ($W m^{-1} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,nbdl), INTENT (out)     :: stempdiag         !! Diagnostic temperature profile @tex ($K$)
                                                                                  !! @endtex                                                                                  
!Isa++
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(out)     :: profil_froz
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)      	    :: snow               !! Snow quantity   
!Isa --
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)           :: pb                !! surface pressure
    REAL(r_std), DIMENSION (kjpindex,nsnow),INTENT(in)       :: snowtemp,snowrho,snowdz
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)    :: organic_layer_thick      !! how deep is the organic soil?    
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm), INTENT (in)  :: soilc_total       !! total soil carbon for use in thermal calcs
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT (in)        :: veget_max         !!Fraction of vegetation type 
    REAL(r_std), DIMENSION(kjpindex),   INTENT (in)          :: thawed_humidity   !!specified humidity of thawed soil
    ! local declaration
    INTEGER(i_std)                                :: ier, ji, jg, jv
    REAL(r_std)                                    :: sum
    !
    !
  !! 1. Initialization of the parameters of the vertical discretization and of the attenuation depths
    cstgrnd=SQRT(one_day / pi)
    lskin = SQRT(so_cond / so_capa * one_day / pi)
    fz1 = 0.3_r_std * cstgrnd
    !zalph = deux !this value has been changed to 1.18 in the src_parameter
    !directory if 32 levels have been
    !used
    
  !! 2.  Computing the depth of the thermal levels (numerical nodes) and the layers boundaries
   
    !! Computing the depth of the thermal levels (numerical nodes) and 
    !! the layers boundariesusing the so-called
    !! adimentional variable z' = z/lskin*cstgrnd (with z in m)
    
    !! 2.1 adimensional thicknesses of the layers
    DO jg=1,ngrnd

    !!?? code simplification hopefully possible here with up-to-date compilers !
    !!! This needs to be solved soon. Either we allow CPP options in SECHIBA or the VPP
    !!! fixes its compiler 
    !!!#ifdef VPP5000
      dz2(jg) = fz(REAL(jg,r_std)-undemi+undemi) - fz(REAL(jg-1,r_std)-undemi+undemi)
    !!!#else
    !!!      dz2(jg) = fz(REAL(jg,r_std)) - fz(REAL(jg-1,r_std))
    !!!#endif
    ENDDO
    
    !! 2.2 adimentional depth of the numerical nodes and layers' boudaries
    DO jg=1,ngrnd
      zz(jg)      = fz(REAL(jg,r_std) - undemi)
      zz_coef(jg) = fz(REAL(jg,r_std)-undemi+undemi)
    ENDDO

    !! 2.3 Converting to meters
    DO jg=1,ngrnd
      zz(jg)      = zz(jg) /  cstgrnd * lskin
      zz_coef(jg) = zz_coef(jg) / cstgrnd * lskin 
      dz2(jg)     = dz2(jg) /  cstgrnd * lskin
    ENDDO

    !! 2.4 Computing some usefull constants for the numerical scheme
    DO jg=1,ngrnd-1
      dz1(jg)  = un / (zz(jg+1) - zz(jg))
    ENDDO
    lambda = zz(1) * dz1(1)
    !! 2.5 Making vegetation masks so that we don't bother to calculated pfts on
    !!     gridcells where they don't exist
    veget_mask_2d(:,:) = .TRUE.
    !veget_mask_2d(:,:) = ( veget_max(:,:) .GT. EPSILON(0.))
    WHERE( ALL((.NOT. veget_mask_2d(:,:)), dim=2) )
       veget_mask_2d(:,1) = .TRUE.
    END WHERE

    !! 2.6 Get the wetness profile on the thermal vertical grid from the diagnostic axis
    CALL thermosoil_humlev(kjpindex, shumdiag_perma, thawed_humidity)
    !
    ! Compute long-term soil humidity (for permafrost)
    !CALL setvar_p (wetdiaglong, val_exp,'NO_KEYWORD',wetdiag(:,:)) !has already
    !been considered in thermosoil_init
    ! cette routine veut dire que wetdiaglong=wetdiag si wetdiaglong=val_exp

    !! 2.7 Thermal conductivity at all levels
    if (ok_explicitsnow) then
       CALL thermosoil_getdiff( kjpindex, ptn, wetdiaglong, pkappa, pcapa,&
                pcapa_en,profil_froz, pcappa_supp, organic_layer_thick, soilc_total)
       ! this is for the thin snow in order to prevent the warm surface
       CALL thermosoil_getdiff_thinsnow (kjpindex, ptn, wetdiaglong, snowdz, pkappa, pcapa, pcapa_en,profil_froz)
    else
       !if (ok_thermix_trunc) then
       !    ! pour convergence avec le trunc
       !    CALL thermosoil_getdiff_old_thermix_trunc2( kjpindex, pkappa, pcapa, pcapa_en )
       !else
       !    CALL thermosoil_getdiff_old_thermix_with_snow( kjpindex, ptn, wetdiaglong, snow, pkappa, pcapa, pcapa_en,profil_froz )
       !endif 
    endif

  !! 3. Diagnostics : consistency checks on the vertical grid.
    sum = zero
    DO jg=1,ngrnd
      sum = sum + dz2(jg)
      WRITE (numout,*) zz(jg),sum
    ENDDO

  !! 4. Compute a first diagnostic temperature profile

    CALL thermosoil_diaglev(kjpindex, stempdiag, veget_max)

    IF (long_print) WRITE (numout,*) ' thermosoil_var_init done '

  END SUBROUTINE thermosoil_var_init


  

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_coef
!!
!>\BRIEF        Calculate soil thermal properties, integration coefficients, apparent heat flux,
!! surface heat capacity,  
!!
!! DESCRIPTION	: This routine computes : \n
!!		1. the soil thermal properties. \n 
!!		2. the integration coefficients of the thermal numerical scheme, cgrnd and dgrnd,
!!              which depend on the vertical grid and on soil properties, and are used at the next 
!!              timestep.\n
!!              3. the soil apparent heat flux and surface heat capacity soilflux
!!              and soilcap, used by enerbil to compute the surface temperature at the next
!!              timestep.\n
!!             -  The soil thermal properties depend on water content (wetdiag) and on the presence 
!!              of snow : snow is integrated into the soil for the thermal calculations, ie if there 
!!              is snow on the ground, the first thermal layer(s) consist in snow, depending on the 
!!              snow-depth. If a layer consists out of snow and soil, wheighed soil properties are 
!!              calculated\n
!!             - The coefficients cgrnd and dgrnd are the integration
!!              coefficients for the thermal scheme \n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!              They correspond respectively to $\beta$ and $\alpha$ from F. Hourdin\'s thesis and 
!!              their expression can be found in this document (eq A19 and A20)
!!             - soilcap and soilflux are the apparent surface heat capacity and flux
!!               used in enerbil at the next timestep to solve the surface
!!               balance for Ts (EQ3); they correspond to $C_s$ and $F_s$ in F.
!!               Hourdin\'s PhD thesis and are expressed in eq. A30 and A31. \n
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflux+otherfluxes(Ts(t)) \n
!!                                      -- EQ3 --\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cgrnd, dgrnd, pcapa, pkappa, soilcap, soilflx
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_coef (kjpindex, dtradia, temp_sol_new, snow, ptn, soilcap, soilflx, zz, dz1, dz2, z1, zdz1,&
           & zdz2, cgrnd, dgrnd, zdz1_soil, zdz2_soil, cgrnd_soil, dgrnd_soil, pcapa, pcapa_en, pkappa, &
!Isa
	& profil_froz, pcappa_supp, stempdiag, &
        & organic_layer_thick, soilc_total,veget_max,snowdz)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                             :: kjpindex     !! Domain size (unitless)
    REAL(r_std), INTENT (in)                               :: dtradia      !! Time step in seconds @tex ($s$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new !! soil surface temperature @tex ($K$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: snow         !! snow mass @tex ($Kg$) @endtex
    REAL(r_std), DIMENSION (ngrnd), INTENT(in)             :: zz           !! depths of the soil thermal numerical nodes 
                                                                           !! @tex ($m$) @endtex 
    REAL(r_std), DIMENSION (ngrnd), INTENT(in)             :: dz1          !! numerical constant depending on the vertical 
                                                                           !! thermal grid only @tex ($m^{-1}$) @endtex 
    REAL(r_std), DIMENSION (ngrnd), INTENT(in)             :: dz2          !! thicknesses of the soil thermal layers
                                                                           !! @tex ($m$) @endtex 
    ! Isa: ptn devient inout alors qu'il était in dans le trunk                                                                       
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT (inout)   :: ptn          !! vertically discretized soil temperatures  
                                                                           !! @tex ($K$) @endtex
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)            :: veget_max         !!Fraction of vegetation type
    REAL(r_std), DIMENSION(kjpindex),   INTENT (in)               :: organic_layer_thick !! how deep is the organic soil?
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm),   INTENT (in)     :: soilc_total !! total soil carbon for use in thermal calcs
    REAL(r_std), DIMENSION(kjpindex,nsnow),   INTENT (in)     :: snowdz !! snow depth
    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: soilcap      !! surface heat capacity
                                                                           !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: soilflx      !! surface heat flux @tex ($W m^{-2}$) @endtex,
                                                                           !! positive towards the 
                                                                           !! soil, writen as Qg (ground heat flux) in the history 
                                                                           !! files.
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: z1           !! numerical constant @tex ($W m^{-1} K^{-1}$) @endtex

    REAL(r_std), DIMENSION (kjpindex,ngrnd-1,nvm), INTENT(out) :: cgrnd        !! matrix coefficient for the computation of soil 
                                                                           !! temperatures (beta in F. Hourdin thesis)
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1,nvm), INTENT(out) :: dgrnd        !! matrix coefficient for the computation of soil 
                                                                           !! temperatures (alpha in F. Hourdin thesis)
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1,nvm), INTENT(out) :: zdz1         !! numerical (buffer) constant 
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex

    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(out)   :: zdz2         !! numerical (buffer) constant  
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex
!Isa
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(out)     :: profil_froz
!Isa E
     REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pcappa_supp
     REAL(r_std),DIMENSION(kjpindex,nbdl),INTENT(inout)    :: stempdiag
     REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: cgrnd_soil
     REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: dgrnd_soil
     REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: zdz1_soil
     REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: zdz2_soil

    !! 0.3 Modified variable

    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(inout) :: pcapa        !! volumetric vertically discretized soil heat capacity
                                                                           !! @tex ($J K^{-1} m^{-3}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(inout) :: pcapa_en     !! volumetric vertically discretized heat capacity used 
                                                                           !! to calculate surfheat_incr
                                                                           !! @tex ($J K^{-1} m^{-3}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd,nvm), INTENT(inout) :: pkappa       !! vertically discretized soil thermal conductivity 
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex

    !! 0.4 Local variables
    REAL(r_std), DIMENSION (kjpindex,nvm)              :: soilcap_pft !!
    REAL(r_std), DIMENSION (kjpindex,nvm)              :: soilflx_pft !!

    INTEGER(i_std)                                         :: ji, jg,jv
    REAL(r_std), DIMENSION(kjpindex)                       :: snow_h       !! snow_h is the snow height @tex ($m$) @endtex 
    REAL(r_std), DIMENSION(kjpindex)                       :: zx1, zx2     !! zx1 and zx2 are the layer fraction consisting in snow
                                                                           !! and soil respectively.
!_ ================================================================================================================================

  !! 1. Computation of the soil thermal properties
    ! Computation of the soil thermal properties; snow properties are also accounted for


    ! Isa
    if (ok_explicitsnow) then
       CALL thermosoil_getdiff( kjpindex, ptn, wetdiaglong, pkappa, pcapa,&
                pcapa_en,profil_froz, pcappa_supp, organic_layer_thick, soilc_total)
    else
       !if (ok_thermix_trunc) then
       !    ! pour convergence avec le trunc
       !    CALL thermosoil_getdiff_old_thermix_trunc2( kjpindex, pkappa, pcapa, pcapa_en )
       !else
       !    CALL thermosoil_getdiff_old_thermix_with_snow( kjpindex, ptn, wetdiaglong, snow, pkappa, pcapa, pcapa_en,profil_froz )
       !endif
    endif

    !added by Tao for very thin snow layer
!Isa
if (ok_Ecorr) then
    CALL thermosoil_readjust(kjpindex, ptn)
endif

!if (control%ok_converge_isaorig) then
!if (1.eq.0) then
!    write(*,*) 'thermosoil 1042: call thermosoil_diaglev'
!    CALL thermosoil_diaglev(kjpindex, stempdiag)
!endif !if (control%ok_converge_isaorig) then
! sinon, on le fait plutôt à la fin de thermosoil_profile
! en fait, ça ne change rien, donc même quand on teste convergence avec isa, on
! appelle thermosoil_diaglev a la fin de thermosoil_profile

    !! 2. computation of the coefficients of the numerical integration scheme

    ! cgrnd, dgrnd

    !! 2.1.  some "buffer" values
     DO jv=1,nvm
       DO jg=1,ngrnd
          WHERE (veget_mask_2d(:,jv))
           zdz2(:,jg,jv)=pcapa(:,jg,jv) * dz2(jg)/dtradia
          ENDWHERE
       ENDDO
       DO jg=1,ngrnd-1
          WHERE (veget_mask_2d(:,jv))
           zdz1(:,jg,jv) = dz1(jg) * pkappa(:,jg,jv)
          ENDWHERE
       ENDDO !DO jg=1,ngrnd-1
    
    !! 2.2.  the coefficients ! 
       WHERE (veget_mask_2d(:,jv))
        z1(:) = zdz2(:,ngrnd,jv) + zdz1(:,ngrnd-1,jv)
        cgrnd(:,ngrnd-1,jv) = (phigeoth + zdz2(:,ngrnd,jv) * ptn(:,ngrnd,jv)) / z1(:)
        dgrnd(:,ngrnd-1,jv) = zdz1(:,ngrnd-1,jv) / z1(:)
       ENDWHERE
       DO jg = ngrnd-1,2,-1
         WHERE (veget_mask_2d(:,jv))
          z1(:) = un / (zdz2(:,jg,jv) + zdz1(:,jg-1,jv) + zdz1(:,jg,jv) * (un - dgrnd(:,jg,jv)))
          cgrnd(:,jg-1,jv) = (ptn(:,jg,jv) * zdz2(:,jg,jv) + zdz1(:,jg,jv) * cgrnd(:,jg,jv)) * z1(:)
          dgrnd(:,jg-1,jv) = zdz1(:,jg-1,jv) * z1(:)
         ENDWHERE
       ENDDO

     !! 3. Computation of the apparent ground heat flux 
       
       !! Computation of the apparent ground heat flux (> towards the soil) and
       !! apparent surface heat capacity, used at the next timestep by enerbil to
       !! compute the surface temperature.
       WHERE (veget_mask_2d(:,jv)) !soil
         soilflx_pft(:,jv) = zdz1(:,1,jv) * (cgrnd(:,1,jv) + (dgrnd(:,1,jv)-1.) * ptn(:,1,jv))
         soilcap_pft(:,jv) = (zdz2(:,1,jv) * dtradia + dtradia * (un - dgrnd(:,1,jv)) * zdz1(:,1,jv))
         z1(:) = lambda * (un - dgrnd(:,1,jv)) + un
         soilcap_pft(:,jv) = soilcap_pft(:,jv) / z1(:)
         soilflx_pft(:,jv) = soilflx_pft(:,jv) + &
            & soilcap_pft(:,jv) * (ptn(:,1,jv) * z1(:) - lambda * cgrnd(:,1,jv) - temp_sol_new(:)) / dtradia 
       ENDWHERE
    ENDDO

    ! 4 here is where I normalize to take the weighted means of each of the
    ! PFTs for surface energy fluxes
    soilflx(:) = zero
    soilcap(:) = zero
    cgrnd_soil(:) = zero
    dgrnd_soil(:) = zero
    zdz1_soil(:) = zero
    zdz2_soil(:) = zero

    DO ji = 1,kjpindex
          DO jv=1,nvm !pft
             IF ( veget_mask_2d(ji,jv) .AND. sum(snowdz(ji,:)) .LE. 0.01) THEN
                soilflx(ji) = soilflx(ji) + (soilflx_pft(ji,jv)*veget_max(ji,jv))
                soilcap(ji) = soilcap(ji) + (soilcap_pft(ji,jv)*veget_max(ji,jv))

                cgrnd_soil(ji) = cgrnd_soil(ji) + (cgrnd(ji,1,jv)*veget_max(ji,jv))
                dgrnd_soil(ji) = dgrnd_soil(ji) + (dgrnd(ji,1,jv)*veget_max(ji,jv))
                zdz1_soil(ji)  = zdz1_soil(ji)  + (zdz1(ji,1,jv)*veget_max(ji,jv))
                zdz2_soil(ji)  = zdz2_soil(ji)  + (zdz2(ji,1,jv)*veget_max(ji,jv))

             END IF
          END DO
    END DO

    IF (long_print) WRITE (numout,*) ' thermosoil_coef done '

  END SUBROUTINE thermosoil_coef
  

  !! ================================================================================================================================
!! SUBROUTINE   : thermosoil_profile
!!
!>\BRIEF        In this routine solves the numerical soil thermal scheme, ie calculates the new soil temperature profile; 
!! This profile is then exported onto the diagnostic axis (call thermosoil_diaglev)
!!
!! DESCRIPTION	: The calculation of the new soil temperature profile is based on
!! the cgrnd and dgrnd values from the previous timestep and the surface temperature Ts aka temp_sol_new. (see detailed
!! explanation in the header of the thermosoil module or in the reference).\n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k)\n
!!                                      -- EQ1 --\n
!!                           Ts=(1-lambda)*T(1) -lambda*T(2)\n 
!!                                      -- EQ2--\n
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): ptn (soil temperature profile on the thermal axis), 
!!                          stempdiag (soil temperature profile on the diagnostic axis)
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

 SUBROUTINE thermosoil_profile (kjpindex, temp_sol_new, ptn, stempdiag,&
                                pkappa_snow,snowdz,snowtemp,grndflux,dtradia,veget_max,cgrnd_snow,dgrnd_snow)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex       !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new   !! Surface temperature at the present time-step 
                                                                               !! @tex ($K$) @endtex
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in):: veget_max             !! Fraction of vegetation type 
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out)      :: stempdiag      !! diagnostic temperature profile 
                                                                               !! @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: cgrnd_snow
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: dgrnd_snow
    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex,ngrnd,nvm), INTENT (inout)   :: ptn            !! vertically discretized soil temperatures 
                                                                               !! @tex ($K$) @endtex


    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jg, jv
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in) :: snowdz
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in) :: snowtemp
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in) :: pkappa_snow
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)       :: grndflux
    REAL(r_std),INTENT(in) :: dtradia

!_ ================================================================================================================================
    
  !! 1. Computes the soil temperatures ptn.

    !! 1.1. ptn(jg=1) using EQ1 and EQ2
    DO jv = 1,nvm 
      DO ji = 1,kjpindex
         IF ( veget_mask_2d(ji,jv) ) THEN
           IF (ok_explicitsnow .AND. SUM(snowdz(ji,:)) .GT. 0.01) THEN
                 ptn(ji,1,jv) = cgrnd_snow(ji,nsnow) + dgrnd_snow(ji,nsnow) * snowtemp(ji,nsnow)
           ELSE
                 ptn(ji,1,jv) = (lambda * cgrnd(ji,1,jv) + temp_sol_new(ji)) / (lambda *(un - dgrnd(ji,1,jv)) + un)
           ENDIF
         ENDIF 
      ENDDO
     
      !! 1.2. ptn(jg=2:ngrnd) using EQ1.
      DO jg = 1,ngrnd-1
        DO ji = 1,kjpindex
         IF ( veget_mask_2d(ji,jv) ) THEN 
          ptn(ji,jg+1,jv) = cgrnd(ji,jg,jv) + dgrnd(ji,jg,jv) * ptn(ji,jg,jv)
         ENDIF
        ENDDO
      ENDDO
    ENDDO
  !! 2. Put the soil temperatures onto the diagnostic axis 
  
    !! Put the soil temperatures onto the diagnostic axis for convenient
    !! use in other routines (stomate..)
    CALL thermosoil_diaglev(kjpindex, stempdiag, veget_max)

    IF (long_print) WRITE (numout,*) ' thermosoil_profile done '

  END SUBROUTINE thermosoil_profile
!!
!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_diaglev
!!
!>\BRIEF        Interpolation of the soil in-depth temperatures onto the diagnostic profile.
!!
!! DESCRIPTION  : This is a very easy linear interpolation, with intfact(jd, jg) the fraction
!! the thermal layer jg comprised within the diagnostic layer jd. The depths of
!! the diagnostic levels are diaglev(1:nbdl), computed in slowproc.f90.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): stempdiag (soil temperature profile on the diagnostic axis)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_diaglev(kjpindex, stempdiag, veget_max)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                          :: kjpindex       !! Domain size (unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in):: veget_max        !! Fraction of vegetation type 
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (out) :: stempdiag      !! Diagnostoc soil temperature profile @tex ($K$) @endtex
    
    !! 0.3 Modified variables

    !! 0.4 Local variables

    INTEGER(i_std)                                      :: ji, jd, jg,jv
    REAL(r_std)                                         :: lev_diag, prev_diag, lev_prog, prev_prog
    REAL(r_std), SAVE, ALLOCATABLE, DIMENSION(:,:)      :: intfact
    REAL(r_std),DIMENSION (kjpindex,ngrnd)              :: ptnmoy
    LOGICAL, PARAMETER                                  :: check=.FALSE.
!_ ================================================================================================================================
    
  !! 1. Computes intfact(jd, jg)

    !! Computes intfact(jd, jg), the fraction
    !! the thermal layer jg comprised within the diagnostic layer jd.

    IF ( .NOT. ALLOCATED(intfact)) THEN
        
        ALLOCATE(intfact(nbdl, ngrnd))
        
        prev_diag = zero
        DO jd = 1, nbdl
          lev_diag = diaglev(jd)
          prev_prog = zero
          DO jg = 1, ngrnd
             IF ( jg == ngrnd .AND. (prev_prog + dz2(jg)) < lev_diag ) THEN
                lev_prog = lev_diag
             ELSE
                lev_prog = prev_prog + dz2(jg)
             ENDIF
            intfact(jd,jg) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog),&
                        & zero)/(lev_diag-prev_diag)
            prev_prog = lev_prog
          ENDDO
           prev_diag = lev_diag
        ENDDO

        IF ( check ) THEN
           WRITE(numout,*) 'thermosoil_diagev -- thermosoil_diaglev -- thermosoil_diaglev --' 
           DO jd = 1, nbdl
              WRITE(numout,*) jd, '-', intfact(jd,1:ngrnd)
           ENDDO
           WRITE(numout,*) "SUM -- SUM -- SUM SUM -- SUM -- SUM"
           DO jd = 1, nbdl
              WRITE(numout,*) jd, '-', SUM(intfact(jd,1:ngrnd))
           ENDDO
           WRITE(numout,*) 'thermosoil_diaglev -- thermosoil_diaglev -- thermosoil_diaglev --' 
        ENDIF
        
    ENDIF

 !! 2. does the interpolation
    ptnmoy(:,:) = 0.
    DO jv = 1, nvm
      DO jg = 1, ngrnd
        ptnmoy(:,jg) = ptnmoy(:,jg) + ptn(:,jg,jv)*veget_max(:,jv)
      ENDDO
    ENDDO

    stempdiag(:,:) = zero
    DO jg = 1, ngrnd
      DO jd = 1, nbdl
        DO ji = 1, kjpindex
           stempdiag(ji,jd) = stempdiag(ji,jd) + ptnmoy(ji,jg)*intfact(jd,jg)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE thermosoil_diaglev

!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_humlev
!!
!>\BRIEF        Interpolates the diagnostic soil humidity profile shumdiag_perma(nbdl, diagnostic axis) onto 
!!              the thermal axis, which gives wetdiag(ngrnd, thermal axis).
!!
!! DESCRIPTION  : Same as in thermosoil_diaglev : This is a very easy linear interpolation, with intfactw(jd, jg) the fraction
!! the thermal layer jd comprised within the diagnostic layer jg. 
!!?? I would think wise to change the indeces here, to keep jD for Diagnostic
!!?? and jG for Ground thermal levels...
!! 
!! The depths of the diagnostic levels are diaglev(1:nbdl), computed in slowproc.f90.
!! Recall that when the 11-layer hydrology is used,
!! wetdiag and shumdiag_perma are with reference to the moisture content (mc)
!! at the wilting point mcw : wetdiag=(mc-mcw)/(mcs-mcw).
!! with mcs the saturated soil moisture content.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): wetdiag (soil humidity profile on the thermal axis)
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
  SUBROUTINE thermosoil_humlev(kjpindex, shumdiag_perma, thawed_humidity)
  
  !! 0. Variables and parameter declaration

    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                            :: kjpindex    !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)    :: shumdiag_perma    !! Relative soil humidity on the diagnostic axis. 
                                                                         !! (0-1, unitless). Caveats : when "hydrol" (the 11-layers
                                                                         !! hydrology) is used, this humidity is calculated with 
                                                                         !! respect to the wilting point : 
                                                                         !! shumdiag_perma= (mc-mcw)/(mcs-mcw), with mc : moisture 
                                                                         !! content; mcs : saturated soil moisture content; mcw: 
                                                                         !! soil moisture content at the wilting point. when the 2-layers
                                                                         !! hydrology "hydrolc" is used, shumdiag_perma is just
                                                                         !! a diagnostic humidity index, with no real physical 
                                                                         !! meaning.
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                       :: ji, jd, jg, jv
    REAL(r_std)                                          :: lev_diag, prev_diag, lev_prog, prev_prog
    REAL(r_std), DIMENSION(ngrnd,nbdl)                   :: intfactw     !! fraction of each diagnostic layer (jd) comprized within
                                                                         !! a given thermal layer (jg)(0-1, unitless) 
    INTEGER(i_std), SAVE  :: proglevel_bottomdiaglev !! for keeping track of where the base of the diagnostic level meets the prognostic level
    INTEGER(i_std), SAVE  :: proglevel_zdeep         !! for keeping track of where the prognostic levels meet zdeep
    LOGICAL               :: at_zdeep=.FALSE.
    LOGICAL               :: at_bottomdiaglev=.FALSE.
    REAL(r_std), DIMENSION(kjpindex), INTENT (in)  :: thawed_humidity !! specified humidity of thawed soil

    LOGICAL, PARAMETER :: check=.FALSE.

!_ ================================================================================================================================
    
  !! 1. computes intfactw(jd,jg), the fraction of each diagnostic layer (jg) comprized within a given thermal layer (jd)
    IF ( check ) &
         WRITE(numout,*) 'thermosoil_humlev --' 

    wetdiag(:,:,:) = zero
    prev_diag = zero
    DO jd = 1, ngrnd
       lev_diag = prev_diag + dz2(jd)
       prev_prog = zero
       DO jg = 1, nbdl
          IF ( jg == nbdl .AND. diaglev(jg) < lev_diag ) THEN
             lev_prog = lev_diag
          ELSE
             lev_prog = diaglev(jg)
          ENDIF
          intfactw(jd,jg) = MAX(MIN(lev_diag,lev_prog)-MAX(prev_diag, prev_prog), zero)/(lev_diag-prev_diag)
          prev_prog = lev_prog
       ENDDO
       prev_diag = lev_diag
    ENDDO

    !!calculate the indices where the thermodynamic levels meet the base of the
    !!moisture levels and zdeep
    jd = 1
    DO WHILE (jd .LT. ngrnd .AND. (.not. at_zdeep ) )
       IF (zz(jd) .GE. z_deepsoil) THEN
          at_zdeep = .TRUE.
          proglevel_zdeep = jd
       END IF
       jd = jd + 1
    END DO
    !
    jd = 1
    DO WHILE (jd .LT. ngrnd .AND. ( .not. at_bottomdiaglev ) )
       IF (zz(jd) .GE. diaglev(nbdl)) THEN
          at_bottomdiaglev = .TRUE.
          proglevel_bottomdiaglev = jd
       END IF
       jd = jd + 1
    END DO

    IF ( check ) THEN
       WRITE(*,*) 'cdk: proglevel_zdeep = ', proglevel_zdeep
       WRITE(*,*) 'cdk: proglevel_bottomdiaglev = ', proglevel_bottomdiaglev
    END IF

    IF (.NOT. satsoil ) THEN
       !++cdk separate permafrost and non-permafrost
       ! only to z_deep for the permafrost
       DO jd = 1, proglevel_zdeep
          WHERE ( veget_mask_2d(:,:) )
             wetdiag(:,jd,:) = 0.0
          END WHERE
       ENDDO


       DO jv = 1, nvm
          DO ji = 1, kjpindex
             IF ( veget_mask_2d(ji,jv) ) THEN
                DO jg = 1, nbdl
                   DO jd = 1, proglevel_zdeep
                      wetdiag(ji,jd,jv) = wetdiag(ji,jd,jv) + shumdiag_perma(ji,jg)*intfactw(jd,jg)
                   END DO
                ENDDO
             END IF
          END DO
       END DO


       ! now update the deep permafrost soil moisture separately
       CALL update_deep_soil_moisture(kjpindex, intfactw, shumdiag_perma,proglevel_bottomdiaglev, proglevel_zdeep, &
            thawed_humidity)

    ELSE
       wetdiag(:,:,:) = 1.
    ENDIF


    !IF ( check ) &
!         WRITE(numout,*) 'thermosoil_humlev --' 
!        write(*,*) 'thermosoil_humlev 1205: wetdiag=',wetdiag(1,1,1)

  END SUBROUTINE thermosoil_humlev

    !!-------------------------------------------------------------------------------

    SUBROUTINE update_deep_soil_moisture (kjpindex, intfactw, shumdiag_perma, proglevel_bottomdiaglev, &
         proglevel_zdeep, thawed_humidity)
    INTEGER(i_std), INTENT(in)                            :: kjpindex    !! Domain size
    REAL(r_std),DIMENSION (kjpindex,nbdl), INTENT (in)    :: shumdiag_perma      !! Diagnostoc profile
    REAL(r_std), DIMENSION(ngrnd, nbdl), INTENT (in)      :: intfactw
    REAL(r_std), DIMENSION(kjpindex),   INTENT (in)       :: thawed_humidity !specified humidity of thawed soil
    INTEGER(i_std), INTENT (in)  :: proglevel_bottomdiaglev !! for keeping track of where the base of the diagnostic level meets the prognostic level
    INTEGER(i_std), INTENT (in)  :: proglevel_zdeep         !! for keeping track of where the prognostic levels meet zdeep
    INTEGER(i_std) :: ji, jd, zi_deep, jv                   !!
    IF (long_print) WRITE (numout,*) 'entering update_deep_soil_misture'


    DO ji = 1, kjpindex
       DO jv = 1,nvm
          IF ( veget_mask_2d(ji,jv) ) THEN
             DO jd = proglevel_zdeep, ngrnd
                IF ( (ptn(ji,jd,jv) .GT. (ZeroCelsius+fr_dT/2.)) ) THEN
                   wetdiag(ji,jd,jv) = thawed_humidity(ji)
                END IF
             END DO
          END IF
       END DO
    END DO

    DO jd =  proglevel_bottomdiaglev, proglevel_zdeep-1
       DO ji = 1, kjpindex
          DO jv = 1,nvm
             IF (veget_mask_2d(ji,jv)) THEN
                CALL lint (diaglev(nbdl), shumdiag_perma(ji,nbdl), z_deepsoil,wetdiag(ji,proglevel_zdeep,jv), &
                     zz(jd), wetdiag(ji,jd,jv), 1)
             END IF
          END DO
       END DO
    END DO

    IF (long_print) WRITE (numout,*) ' update_deep_soil_misture done'
    
    END SUBROUTINE update_deep_soil_moisture

  SUBROUTINE lint(x1,y1,x2,y2,x,y,NY)
    ! Interpolation linéaire entre des points (x1,y1) et (x2,y2)
    ! Ces commentaires en mauvais français permettent savoir qui a
    ! ecrit la subroutine :-) - DK
    !
    IMPLICIT NONE
    REAL, PARAMETER                    :: EPSILON = 1.E-10
    REAL, INTENT(in)                   ::  x1,x2,y1,y2,x
    INTEGER, INTENT(in)                ::  NY
    REAL, DIMENSION(NY), INTENT(inout) :: y
    INTEGER i
    
    IF (ABS(x1 - x2) .LT. EPSILON) THEN
       PRINT *, 'ERROR IN lint(x1,y1,x2,y2,y,NY) : x1==x2!'
       PRINT *, 'x1=',x1,'  x2=',x2
       PRINT *, 'y1=',y1,'  y2=',y2
       STOP
    END IF
    
    IF (x1 .LE. x .AND. x .LE. x2) THEN
       y = x*(y2-y1)/(x2-x1) + (y1*x2 - y2*x1)/(x2-x1)
       !      ELSE
       !        y = UNDEF
    END IF
    
  END SUBROUTINE lint


!!
!================================================================================================================================
!! SUBROUTINE   : thermosoil_energy
!!
!>\BRIEF         Energy check-up.
!!
!! DESCRIPTION  : I didn\'t comment this routine since at do not understand its
!use, please
!! ask initial designers (Jan ? Nathalie ?).
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ??
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_
!===============================================================================================================
  SUBROUTINE thermosoil_energy(kjpindex, temp_sol_new, soilcap, first_call, veget_max)
    ! interface description
    ! input scalar
    INTEGER(i_std), INTENT(in)                            :: kjpindex    !! Domain size
    LOGICAL, INTENT (in)                                 :: first_call  !!
    ! input fields
    !
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: temp_sol_new !! New soil temperature
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)        :: soilcap      !! Soil capacity
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in):: veget_max        !! Fraction of vegetation type 
    !
    ! modified fields
    !
    ! output fields
    !
    ! local variable
    !
    INTEGER(i_std)  :: ji, jg, jv
    !
    !
    IF (long_print) WRITE (numout,*) 'entering thermosoil_energy'
    !
    IF (first_call) THEN

     DO ji = 1, kjpindex
      surfheat_incr(ji) = zero
      coldcont_incr(ji) = zero
      temp_sol_beg(ji)  = temp_sol_new(ji)
      !
      DO jg = 1, ngrnd
       ptn_beg(ji,jg,:)   = ptn(ji,jg,:)
      ENDDO
      !
     ENDDO
    
     RETURN

    ENDIF

     DO ji = 1, kjpindex
      surfheat_incr(ji) = zero
      coldcont_incr(ji) = zero
     ENDDO
    !
    !  Sum up the energy content of all layers in the soil.
    !
    DO ji = 1, kjpindex
    !
       IF (SUM(pcapa_en(ji,1,:)*veget_max(ji,:)) .LE. sn_capa) THEN
          !
          ! Verify the energy conservation in the surface layer
          !
          coldcont_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          surfheat_incr(ji) = zero
       ELSE
          !
          ! Verify the energy conservation in the surface layer
          !
          surfheat_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          coldcont_incr(ji) = zero
       ENDIF
    ENDDO
    
    ptn_beg(:,:,:)      = ptn(:,:,:)
    temp_sol_beg(:)   = temp_sol_new(:)

  END SUBROUTINE thermosoil_energy

!Isa++
  SUBROUTINE thermosoil_readjust(kjpindex, ptn)
    ! interface description
    ! input scalar
    INTEGER(i_std), INTENT(in)                             :: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(inout)    :: ptn
    INTEGER(i_std)  :: ji, jg, jv
    REAL(r_std) :: ptn_tmp
    !WARNING : pcapa est en J/K/m; pcappa_supp est en J/K
    !local
    REAL(r_std), DIMENSION(kjpindex, ngrnd, ngrnd) :: intfact_soil_th
    do jv = 1,nvm
      do jg=1, ngrnd
          do ji=1, kjpindex
            if (veget_mask_2d(ji,jv)) then
               !Isa : here, all soil latent energy is put into E_sol_lat_couche(ji, 1)
               !because the variable soil layers make it difficult to keep track of all
               !layers in this version
               E_sol_lat_couche(ji, 1,jv)=E_sol_lat_couche(ji, 1, jv)+pcappa_supp(ji,jg,jv)*(ptn(ji,jg,jv)-ptn_beg(ji,jg,jv))
            endif
          enddo
      enddo
    enddo

!Isa test
    do jv = 1,nvm
      do ji=1, kjpindex
       if (veget_mask_2d(ji,jv)) then
        if (E_sol_lat_couche(ji,1,jv).GT.min_sechiba.AND.MINVAL(ptn(ji,:,jv)).GT.ZeroCelsius+fr_dT/2.) then
                !1. normalement, à cette température, il n'y a plus de neige
                !donc pas la peine d'utiliser x1 et x2
                !2. répartition de l'excess d'énergie sur 2.7m = 6 premiers
                !niveaux..
         !plus d'énergie au dégel qu'au gel; cette énergie aurait dû chauffer au lieu de dégeler 
                !=> on monte les températures en se limitant à 0.5°C à chaque pas de temps..
                do jg=1,6
                ptn_tmp=ptn(ji,jg,jv)

                ptn(ji,jg,jv)=ptn(ji,jg,jv)+min(E_sol_lat_couche(ji,1,jv)/pcapa(ji,jg,jv)/zz_coef(6), 0.5)
                E_sol_lat_couche(ji,1,jv)=E_sol_lat_couche(ji,1,jv)-(ptn(ji,jg,jv)-ptn_tmp)*pcapa(ji,jg,jv)*dz2(jg)
                enddo
                else if (E_sol_lat_couche(ji,1,jv).LT.-min_sechiba.AND.MINVAL(ptn(ji,:,jv)).GT.ZeroCelsius+fr_dT/2.) then
                !pas assez d'énergie lors du dégel. Cette énergie aurait dû dégeler au lieu de chauffer;
                !=> on rabaisse les températures en conséquence
                do jg=1,6
                ptn_tmp=ptn(ji,jg,jv)

                ptn(ji,jg,jv)=max(ZeroCelsius+fr_dT/2., ptn_tmp+E_sol_lat_couche(ji,1,jv)/pcapa(ji,jg,jv)/zz_coef(6))
                E_sol_lat_couche(ji,1,jv)=E_sol_lat_couche(ji,1,jv)+(ptn_tmp-ptn(ji,jg,jv))*pcapa(ji,jg,jv)*dz2(jg)
               enddo
        endif 
      endif
     enddo
    enddo

  END SUBROUTINE thermosoil_readjust
! Isa --
  !Isa ajout Bruno
  !-------------------------------------------------------------------

  SUBROUTINE thermosoil_wlupdate( kjpindex, dt, ptn, hsd, hsdlong )
  ! Updates the long-term soil humidity
    INTEGER(i_std),INTENT(in)				:: kjpindex
    REAL(r_std),INTENT(in)			 	:: dt
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)	:: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)   	:: hsd
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(inout) 	:: hsdlong
    ! CR: debugage pour voir où c'est actif
    integer(i_std) ::  il,ji,jg,nactif
    REAL(r_std), PARAMETER               :: tau_freezesoil = 30.*86400. 

    !
    do il = 1, ndeep
       WHERE ( ( ptn(:,il,:) .GT. ZeroCelsius + fr_dT/2. ) .AND. veget_mask_2d(:,:))
          hsdlong(:,il,:) = ( hsd(:,il,:) * dt + hsdlong(:,il,:) *(tau_freezesoil-dt) ) / tau_freezesoil
       ENDWHERE
    end do

    IF (long_print) WRITE (numout,*) 'entering thermosoil_wlupdate'

   END SUBROUTINE thermosoil_wlupdate
   
!-------------------------------------------------------------------

!Isa ajout from thermosoil_bruno
  SUBROUTINE thermosoil_getdiff( kjpindex, ptn, wetdiaglong, pkappa, pcapa, pcapa_en,&
                      & profil_froz, pcappa_supp, organic_layer_thick, soilc_total)
   !
   ! Computes soil heat capacity and conductivity   
   !
    INTEGER(i_std),INTENT(in)				:: kjpindex
!Isa E_corr : in -> inout
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(inout)	:: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)   	:: wetdiaglong
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pkappa
!Isa
     REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pcappa_supp
!Isa modif bruno
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: profil_froz
    !    
    REAL						:: x, p, E_supp, Epaiss !Isa
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm) 		:: zx1, zx2    
    ! Heat capacity of ice/water mixture
    REAL			   			:: cap_iw
    ! Thermal conductivity for saturated soil
    REAL						:: csat
    REAL(r_std)                               :: so_capa_dry_net
    REAL(r_std), DIMENSION(kjpindex,ngrnd,nvm):: poros_net
    REAL(r_std)                               :: cond_solid_net
    REAL(r_std)                               :: so_cond_dry_net
    INTEGER						:: ji,jg,jv
    REAL(r_std), DIMENSION(kjpindex), INTENT (in)           :: organic_layer_thick ! how deep is the organic soil?
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm), INTENT (in) :: soilc_total !! total soil carbon for use in thermal calcs

     ! Organic and anorgaic layer fraction
     !
     ! Default: organic layer not taken into account
     zx1(:,:,:) = 0.
     !
     IF ( use_toporganiclayer_tempdiff ) THEN
       !
       ! level 1
       !
       DO jv = 1,nvm
         DO ji = 1,kjpindex
           IF ( organic_layer_thick(ji) .GT. zz_coef(1) ) THEN
             !! the 1st level is in the organic => the 1st layer is entirely organic
             zx1(ji,1,jv) = 1. !!zx1 being the fraction of each level that is organic, zx2 is the remainder
           ELSE IF ( organic_layer_thick(ji) .GT. zero ) THEN
             !! the 1st level is beyond the organic and the organic is present
             zx1(ji,1,jv) = organic_layer_thick(ji) / zz_coef(1)
           ELSE
             ! there is no organic at all
             zx1(ji,1,jv) = 0.
           ENDIF
         ENDDO
       ENDDO
       !
       ! other levels
       !
       DO jg = 2, ngrnd !- 2
         DO ji = 1,kjpindex
           IF ( organic_layer_thick(ji) .GT. zz_coef(jg) ) THEN
             ! the current level is in the organic => the current layer is
             ! entirely organic
             zx1(ji,jg,1) = 1.
           ELSE IF ( organic_layer_thick(ji) .GT. zz_coef(jg-1) ) THEN
             ! the current layer is partially organic
             zx1(ji,jg,1) = (organic_layer_thick(ji) - zz_coef(jg-1)) / (zz_coef(jg) - zz_coef(jg-1))
           ELSE
             ! both levels are out of organic => the current layer is entirely
             ! mineral soil       
             zx1(ji,jg,1) = 0.
           ENDIF
         ENDDO
       ENDDO
       DO jv = 2, nvm
         zx1(ji,jg,jv) = zx1(ji,jg,1)
       ENDDO
       !
     ELSEIF ( use_soilc_tempdiff ) THEN
       !
       DO jv = 1,nvm
         DO jg = 1, ngrnd
           DO ji = 1,kjpindex
             zx1(ji,jg,jv) = MIN((soilc_total(ji,jg,jv)/soilc_max),1.)   !after lawrence and slater
           ENDDO
         ENDDO
       ENDDO
       !
     ENDIF 
     !
     zx2(:,:,:) = 1.-zx1(:,:,:)

     DO jv = 1,nvm
       DO jg = 1, ngrnd
         DO ji = 1,kjpindex
           IF (veget_mask_2d(ji,jv)) THEN
             !
             ! 1. Calculate dry heat capacity and conductivity, taking
             ! into account the organic and mineral fractions in the layer
             !
             so_capa_dry_net = zx1(ji,jg,jv) * so_capa_dry_org + zx2(ji,jg,jv) * so_capa_dry
             cond_solid_net  = un / ( zx1(ji,jg,jv) / cond_solid_org  + zx2(ji,jg,jv) / cond_solid  )
             poros_net(ji,jg,jv) = zx1(ji,jg,jv) * poros_org + zx2(ji,jg,jv) * poros
             !
             so_cond_dry_net = un / ( zx1(ji,jg,jv) / so_cond_dry_org + zx2(ji,jg,jv) / so_cond_dry )
             !
             ! 2. Calculate heat capacity with allowance for permafrost
             !    For the moment we don't take into account porosity changes related
             !    mx_eau_eau changes
             !    which would be reflected in so_capa_wet. 
             !    so_capa_ice implies porosity of 0.15 (mx_eau_eau=150) -> not
             !    anymore :
             !    Isa : changed to account for poros = 0.4
             !    Isa : changed wetdiag to have the real fraction of porosity filled
             !    with water

             IF (ok_freeze_thermix) THEN
                ! 2.1. soil heat capacity depending on temperature and humidity
                IF (ptn(ji,jg,jv) .LT. ZeroCelsius-fr_dT/2.) THEN
                  ! frozen soil
                  !! this is from Koven's version: pcapa(ji,jg,jv) = so_capa_dry_net + wetdiaglong(ji,jg,jv)*poros_net(ji,jg,jv)*capa_ice*rho_ice
                  pcapa(ji,jg,jv) = so_capa_dry_net + wetdiaglong(ji,jg,jv)*(so_capa_ice - so_capa_dry_net)!Isa : old version, proved to be correct
                  pcappa_supp(ji,jg, jv)= 0.
                  profil_froz(ji,jg,jv) = 1.
                ELSEIF (ptn(ji,jg,jv) .GT. ZeroCelsius+fr_dT/2.) THEN
                  ! unfrozen soil         
                  !! this is from Koven's version: pcapa(ji,jg,jv) =  so_capa_dry_net + wetdiaglong(ji,jg,jv)*poros_net(ji,jg,jv)*capa_water*rho_water
                  pcapa(ji,jg,jv) = so_capa_dry_net + wetdiaglong(ji,jg,jv)*(so_capa_wet - so_capa_dry_net) 
                  pcappa_supp(ji,jg,jv)= 0.
                  profil_froz(ji,jg,jv) = 0.
                ELSE
   
                ! x is the unfrozen fraction of soil water              
                x = (ptn(ji,jg,jv)-(ZeroCelsius-fr_dT/2.)) / fr_dT
                profil_froz(ji,jg,jv) = (1. - x)
                ! net heat capacity of the ice/water mixture
                cap_iw = x * so_capa_wet + (1.-x) * so_capa_ice
                ! cap_iw = x * 4.E6 + (1.-x) * 2.E6 !DKtest - compar. w/ theor. sol. 
                pcapa(ji,jg,jv) = so_capa_dry_net + wetdiaglong(ji,jg,jv)*(cap_iw-so_capa_dry_net) + &
                                  wetdiaglong(ji,jg,jv)*poros_net(ji,jg,jv)*lhf*rho_water/fr_dT
                pcappa_supp(ji,jg,jv)= wetdiaglong(ji,jg,jv)*poros_net(ji,jg,jv)*lhf*rho_water/fr_dT*dz2(jg)

                ENDIF
             ELSE !++cdk this is physically wrong and only to be used to test the influence of latent heat
                pcapa(ji,jg,jv) = so_capa_dry_net + wetdiag(ji,jg,jv)*(so_capa_wet - so_capa_dry_net)
                profil_froz(ji,jg,jv) = 0.
             ENDIF
             !
             ! 3. Calculate the heat conductivity with allowance for permafrost (Farouki,
             ! 1981, Cold Reg. Sci. Technol.)
             !
             ! 3.1. unfrozen fraction
             p = poros_net(ji,jg,jv)
             x = (ptn(ji,jg,jv)-(ZeroCelsius-fr_dT/2.)) / fr_dT * p
             x = MIN( p, MAX( 0., x ) )
             !++cdk: DKorig: x = (ptn(ji,jg)-(ZeroCelsius-fr_dT/2.)) / fr_dT * poros
             !++cdk: DKorig: x = MIN( poros, MAX( 0., x ) )

             ! 3.2. saturated conductivity
             csat = cond_solid_net**(1.-p) * cond_ice**(p-x) * cond_water**x
             !++cdk: DKorig: csat = cond_solid**(1.-poros) * cond_ice**(poros-x)
             !* cond_water**x

             ! 3.3. unsaturated conductivity
             pkappa(ji,jg,jv) = (csat - so_cond_dry_net)*wetdiaglong(ji,jg,jv) + so_cond_dry_net
             !++cdk: DKorig: pkappa(ji,jg) = (csat - so_cond_dry)*humdiag(ji,jg)
             !+ so_cond_dry
             !
           ENDIF
          ENDDO
         ENDDO
        ENDDO

        pcapa_en(:,:,:) = pcapa(:,:,:)


   END SUBROUTINE thermosoil_getdiff


!
     SUBROUTINE thermosoil_getdiff_thinsnow (kjpindex, ptn, wetdiaglong, snowdz, pkappa, pcapa, pcapa_en,profil_froz)
   !
   ! Computes soil heat capacity and conductivity   
   !
    INTEGER(i_std),INTENT(in)                           :: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)        :: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)        :: wetdiaglong
    REAL(r_std),DIMENSION(kjpindex,nsnow),INTENT (in)           :: snowdz
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)       :: pcapa
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pkappa
!modif bruno
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: profil_froz
    !    
    REAL                                                :: x
    REAL(r_std), DIMENSION(kjpindex)                    :: snow_h
    REAL(r_std), DIMENSION(kjpindex,ngrnd)              :: zx1, zx2
    INTEGER                                             :: ji,jg,jv


    DO ji = 1,kjpindex

      ! 1. Determine the fractions of snow and soil

      snow_h(ji) = SUM(snowdz(ji,:))

      IF (snow_h(ji) .LE. 0.01) THEN

         !
         !  1.1. The first level
         !
         IF ( snow_h(ji) .GT. zz_coef(1) ) THEN

             ! the 1st level is in the snow => the 1st layer is entirely snow
             zx1(ji,1) = 1.
             zx2(ji,1) = 0.
                
         ELSE IF ( snow_h(ji) .GT. zero ) THEN

             ! the 1st level is beyond the snow and the snow is present
             zx1(ji,1) = snow_h(ji) / zz_coef(1)
             zx2(ji,1) = ( zz_coef(1) - snow_h(ji)) / zz_coef(1)        
         ENDIF

         !
         DO jv = 1,nvm
          DO jg = 1, 1
           IF (veget_mask_2d(ji,jv)) THEN
            !
            ! 2. Calculate frozen profile for hydrolc.f90
        !
            IF (ptn(ji,jg,jv) .LT. ZeroCelsius-fr_dT/2.) THEN
                profil_froz(ji,jg,jv) = 1.

                 ELSEIF (ptn(ji,jg,jv) .GT. ZeroCelsius+fr_dT/2.) THEN
                profil_froz(ji,jg,jv) = 0.
                 ELSE

                   ! x is the unfrozen fraction of soil water              
                   x = (ptn(ji,jg,jv)-(ZeroCelsius-fr_dT/2.)) / fr_dT   
              profil_froz(ji,jg,jv) = (1. - x)

            ENDIF

            ! 3. heat capacity calculation
        !
            ! 3.0 old heat capacity calculation
            pcapa(ji,jg,jv) = so_capa_dry + wetdiaglong(ji,jg,jv)*(so_capa_wet - so_capa_dry)

        ! 3.1. Still some improvement from the old_version : Take into account the snow and soil fractions in the layer

            pcapa(ji,jg,jv) = zx1(ji,jg) * sn_capa + zx2(ji,jg) * pcapa(ji,jg,jv)

        ! 3.2. Calculate the heat capacity for energy conservation check 
        !        (??, does not influence other results, just written to history file)

        IF ( zx1(ji,jg).GT.0. ) THEN
               pcapa_en(ji,jg,jv) = sn_capa
        ELSE
               pcapa_en(ji,jg,jv) = pcapa(ji,jg,jv)
        ENDIF
            !
            !4. heat conductivity calculation
        !
            !4.0 old heat conductivity calculation
            pkappa(ji,jg,jv) = so_cond_dry + wetdiaglong(ji,jg,jv)*(so_cond_wet - so_cond_dry)

            !4.0 Still some improvement from the old_version : Take into account the snow and soil fractions in the layer

            pkappa(ji,jg,jv) = un / ( zx1(ji,jg) / sn_cond + zx2(ji,jg) / pkappa(ji,jg,jv) )



          ENDIF
         END DO
        END DO
      ENDIF
    ENDDO


   END SUBROUTINE thermosoil_getdiff_thinsnow

!Isa
     SUBROUTINE thermosoil_getdiff_old_thermix_with_snow( kjpindex, ptn, wetdiag, snow, pkappa, pcapa, pcapa_en,profil_froz)
   !
   ! Computes soil heat capacity and conductivity   
   !
    INTEGER(i_std),INTENT(in)				:: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)	:: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)   	:: wetdiag
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pkappa
!modif bruno
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: profil_froz
    !    
    REAL						:: x
    REAL(r_std), DIMENSION(kjpindex)             	:: snow_h
    REAL(r_std), DIMENSION(kjpindex,ngrnd) 		:: zx1, zx2    
    INTEGER						:: ji,jg,jv


    DO ji = 1,kjpindex
      
      ! 1. Determine the fractions of snow and soil
      
      snow_h(ji) = snow(ji) / sn_dens
      
      !
      !  1.1. The first level
      !
      IF ( snow_h(ji) .GT. zz_coef(1) ) THEN

          ! the 1st level is in the snow => the 1st layer is entirely snow
          zx1(ji,1) = 1.
	  zx2(ji,1) = 0.
	  	  
      ELSE IF ( snow_h(ji) .GT. zero ) THEN      
      
          ! the 1st level is beyond the snow and the snow is present
          zx1(ji,1) = snow_h(ji) / zz_coef(1)
          zx2(ji,1) = ( zz_coef(1) - snow_h(ji)) / zz_coef(1)	
	 
      ELSE
      
          ! there is no snow at all, quoi ;-)
          zx1(ji,1) = 0.
	  zx2(ji,1) = 1.       
	  
      ENDIF
      
      !
      !  1.2. The other levels except the two last (too deep to be accounted for by the current hydrology??)
      !

      DO jg = 2, ngrnd !- 2
        IF ( snow_h(ji) .GT. zz_coef(jg) ) THEN

            ! the current level is in the snow => the current layer is entirely snow
            zx1(ji,jg) = 1.
	    zx2(ji,jg) = 0.
	  	  
        ELSE IF ( snow_h(ji) .GT. zz_coef(jg-1) ) THEN
	
   	    ! the current layer is partially snow
            zx1(ji,jg) = (snow_h(ji) - zz_coef(jg-1)) / (zz_coef(jg) - zz_coef(jg-1))
            zx2(ji,jg) = ( zz_coef(jg) - snow_h(ji)) / (zz_coef(jg) - zz_coef(jg-1))

        ELSE
	
	    ! both levels are out of snow => the current layer is entirely soil	  
            zx1(ji,jg) = 0.
	    zx2(ji,jg) = 1.       
	  	
        ENDIF
      ENDDO
            
      DO jv = 1,nvm
       DO jg = 1, ngrnd !-2 
        IF (veget_mask_2d(ji,jv)) THEN
         !
         ! 2. Calculate frozen profile for hydrolc.f90
	 !
	 IF (bavard.GE.8) write(*,'(A24,I3,F5.2,G14.7)')  'lev, z, Soil_temp: ', jg, zz(jg), ptn(ji,jg,jv)
	 
         IF (ptn(ji,jg,jv) .LT. ZeroCelsius-fr_dT/2.) THEN
	    profil_froz(ji,jg,jv) = 1.

	    IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, Soil_temp: ', jg, zz(jg), ptn(ji,jg,jv)	 
	    
     	 ELSEIF (ptn(ji,jg,jv) .GT. ZeroCelsius+fr_dT/2.) THEN
	    profil_froz(ji,jg,jv) = 0.

	    IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, : Soil_temp', jg, zz(jg), ptn(ji,jg,jv)	 
	   
     	 ELSE
	 
     	   ! x is the unfrozen fraction of soil water	   	   
     	   x = (ptn(ji,jg,jv)-(ZeroCelsius-fr_dT/2.)) / fr_dT	   
           profil_froz(ji,jg,jv) = (1. - x)
	   
         ENDIF

         ! 3. heat capacity calculation
	 !
         ! 3.0 old heat capacity calculation
         pcapa(ji,jg,jv) = so_capa_dry + wetdiag(ji,jg,jv)*(so_capa_wet - so_capa_dry)

	 ! 3.1. Still some improvement from the old_version : Take into account the snow and soil fractions in the layer

         pcapa(ji,jg,jv) = zx1(ji,jg) * sn_capa + zx2(ji,jg) * pcapa(ji,jg,jv)

	 ! 3.2. Calculate the heat capacity for energy conservation check 
	 !        (??, does not influence other results, just written to history file)
	 
	 IF ( zx1(ji,jg).GT.0. ) THEN
            pcapa_en(ji,jg,jv) = sn_capa
	 ELSE
            pcapa_en(ji,jg,jv) = pcapa(ji,jg,jv)
	 ENDIF
         !
         !4. heat conductivity calculation
	 !
         !4.0 old heat conductivity calculation
         pkappa(ji,jg,jv) = so_cond_dry + wetdiag(ji,jg,jv)*(so_cond_wet - so_cond_dry)

         !4.0 Still some improvement from the old_version : Take into account the snow and soil fractions in the layer

         pkappa(ji,jg,jv) = un / ( zx1(ji,jg) / sn_cond + zx2(ji,jg) / pkappa(ji,jg,jv) )


	 
       ENDIF
      END DO
     END DO            
      !
    ENDDO   
    
   
   END SUBROUTINE thermosoil_getdiff_old_thermix_with_snow

!Isa
     SUBROUTINE thermosoil_getdiff_old_thermix( kjpindex, ptn, wetdiag, snow, pkappa, pcapa, pcapa_en,profil_froz)
   !
   ! Computes soil heat capacity and conductivity   
   !
    INTEGER(i_std),INTENT(in)				:: kjpindex
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)	:: ptn
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(in)   	:: wetdiag
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pkappa
!modif bruno
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: profil_froz
    !    
    REAL						:: x 
    INTEGER						:: ji,jg,jv

    
    DO ji = 1,kjpindex            
      DO jg = 1, ngrnd
       DO jv = 1, nvm
        IF (veget_mask_2d(ji,jv)) THEN 
         !
         ! 2. Calculate frozen profile for hydrolc.f90
	 !
	 IF (bavard.GE.8) write(*,'(A24,I3,F5.2,G14.7)')  'lev, z, Soil_temp: ', jg, zz(jg), ptn(ji,jg,jv)
	 
         IF (ptn(ji,jg,jv) .LT. ZeroCelsius-fr_dT/2.) THEN
	    profil_froz(ji,jg,jv) = 1.

	    IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, Soil_temp: ', jg, zz(jg), ptn(ji,jg,jv)	 
	    
     	 ELSEIF (ptn(ji,jg,jv) .GT. ZeroCelsius+fr_dT/2.) THEN
	    profil_froz(ji,jg,jv) = 0.

	    IF (bavard.GE.7) write(*,'(A24,I3,F5.2,G14.7 )')  'lev, z, : Soil_temp', jg, zz(jg), ptn(ji,jg,jv)	 
	   
     	 ELSE
	 
     	   ! x is the unfrozen fraction of soil water	   	   
     	   x = (ptn(ji,jg,jv)-(ZeroCelsius-fr_dT/2.)) / fr_dT	   
           profil_froz(ji,jg,jv) = (1. - x)
	   
         ENDIF

         ! 3. heat capacity calculation
	 !

         ! 3.1 old heat capacity calculation
         pcapa(ji,jg,jv) = so_capa_dry + wetdiag(ji,jg,jv)*(so_capa_wet - so_capa_dry)

	 ! 3.2. Calculate the heat capacity for energy conservation check 
	 !        (??, does not influence other results, just written to history file)
	 pcapa_en(ji,1,jv) = so_capa_dry + wetdiag(ji,1,jv)*(so_capa_wet - so_capa_dry)

         !4. heat conductivity calculation
	 !
         !4.0 old heat conductivity calculation
         pkappa(ji,jg,jv) = so_cond_dry + wetdiag(ji,jg,jv)*(so_cond_wet - so_cond_dry)

	 
        END IF
       END DO
      END DO            
      !
    ENDDO   
   
   END SUBROUTINE thermosoil_getdiff_old_thermix

   SUBROUTINE thermosoil_getdiff_old_thermix_trunc( kjpindex, snow, pkappa, pcapa, pcapa_en )

    INTEGER(i_std), INTENT(in) :: kjpindex
    
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pkappa
    INTEGER						:: ji,jg,jv
    REAL(r_std), DIMENSION(kjpindex)                       :: snow_h       !! snow_h is the snow height @tex ($m$) @endtex 
    REAL(r_std), DIMENSION(kjpindex)                       :: zx1, zx2     !! zx1 and zx2 
                             !! are the layer fraction consisting in snowand soil respectively.

     
    ! Computation of the soil thermal properties; snow properties are also accounted for

    DO jv = 1, nvm
       DO ji = 1,kjpindex
        IF (veget_mask_2d(ji,jv)) THEN
          snow_h(ji) = snow(ji) / sn_dens
       
          IF ( snow_h(ji) .GT. zz_coef(1) ) THEN
              pcapa(ji,1,jv) = sn_capa
              pcapa_en(ji,1,jv) = sn_capa
              pkappa(ji,1,jv) = sn_cond
          ELSE IF ( snow_h(ji) .GT. zero ) THEN
              pcapa_en(ji,1,jv) = sn_capa
              zx1(ji) = snow_h(ji) / zz_coef(1)
              zx2(ji) = ( zz_coef(1) - snow_h(ji)) / zz_coef(1)
              pcapa(ji,1,jv) = zx1(ji) * sn_capa + zx2(ji) * so_capa_wet
              pkappa(ji,1,jv) = un / ( zx1(ji) / sn_cond + zx2(ji) / so_cond_wet )
          ELSE
              pcapa(ji,1,jv) = so_capa_dry + wetdiag(ji,1,jv)*(so_capa_wet - so_capa_dry)
              pkappa(ji,1,jv) = so_cond_dry + wetdiag(ji,1,jv)*(so_cond_wet - so_cond_dry)
              pcapa_en(ji,1,jv) = so_capa_dry + wetdiag(ji,1,jv)*(so_capa_wet - so_capa_dry)
          ENDIF
          !
          DO jg = 2, ngrnd - 2
            IF ( snow_h(ji) .GT. zz_coef(jg) ) THEN
                pcapa(ji,jg,jv) = sn_capa
                pkappa(ji,jg,jv) = sn_cond
                pcapa_en(ji,jg,jv) = sn_capa
            ELSE IF ( snow_h(ji) .GT. zz_coef(jg-1) ) THEN
                zx1(ji) = (snow_h(ji) - zz_coef(jg-1)) / (zz_coef(jg) - zz_coef(jg-1))
                zx2(ji) = ( zz_coef(jg) - snow_h(ji)) / (zz_coef(jg) - zz_coef(jg-1))
                pcapa(ji,jg,jv) = zx1(ji) * sn_capa + zx2(ji) * so_capa_wet
                pkappa(ji,jg,jv) = un / ( zx1(ji) / sn_cond + zx2(ji) / so_cond_wet )
                pcapa_en(ji,jg,jv) = sn_capa
            ELSE
                pcapa(ji,jg,jv) = so_capa_dry + wetdiag(ji,jg,jv)*(so_capa_wet - so_capa_dry)
                pkappa(ji,jg,jv) = so_cond_dry + wetdiag(ji,jg,jv)*(so_cond_wet - so_cond_dry)
                pcapa_en(ji,jg,jv) = so_capa_dry + wetdiag(ji,jg,jv)*(so_capa_wet - so_capa_dry)
            ENDIF
          ENDDO
        END IF
       ENDDO ! DO ji = 1,kjpindex
     ENDDO 
    END SUBROUTINE thermosoil_getdiff_old_thermix_trunc

    SUBROUTINE thermosoil_getdiff_old_thermix_trunc2( kjpindex,  pkappa, pcapa, pcapa_en )

    INTEGER(i_std), INTENT(in) :: kjpindex
    
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)  	:: pcapa   
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pcapa_en
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(out)    :: pkappa
    INTEGER						:: ji,jg,jv

    DO jv = 1,nvm  
      DO jg = 1,ngrnd
        DO ji = 1,kjpindex
         IF (veget_mask_2d(ji,jv)) THEN
          pkappa(ji,jg,jv) = so_cond_dry + wetdiag(ji,jg,jv)*(so_cond_wet - so_cond_dry)
          pcapa(ji,jg,jv) = so_capa_dry + wetdiag(ji,jg,jv)*(so_capa_wet - so_capa_dry)
          pcapa_en(ji,jg,jv) = so_capa_dry + wetdiag(ji,jg,jv)*(so_capa_wet - so_capa_dry)
         END IF
        ENDDO
      ENDDO
    ENDDO

    END SUBROUTINE thermosoil_getdiff_old_thermix_trunc2

   !-------------------------------------------------------------------
   SUBROUTINE add_heat_Zimov(kjpindex, dtradia, ptn, heat_Zimov)
     INTEGER(i_std),INTENT(in)                          :: kjpindex
    REAL(r_std), INTENT (in)                           :: dtradia          !! Time step in seconds
    REAL(r_std),DIMENSION(kjpindex,ngrnd,nvm),INTENT(inout)     :: ptn
    REAL(r_std), DIMENSION(kjpindex,ndeep,nvm), INTENT (in)   :: heat_Zimov !! heating associated with decomposition
    INTEGER :: ji, jg, jv

    IF (long_print) WRITE (numout,*) 'entering add_heat_Zimov'

    DO ji = 1, kjpindex
       DO jv = 1,nvm
          IF ( veget_mask_2d(ji,jv) ) THEN
             DO jg = 1, ngrnd
                ptn(ji,jg,jv) = ptn(ji,jg,jv) + heat_zimov(ji,jg,jv) * dtradia / ( pcapa(ji,jg,jv) * dz2(jg) )
             END DO
          END IF
       END DO
    END DO

    IF (long_print) WRITE (numout,*) ' add_heat_Zimov done'

  END SUBROUTINE add_heat_Zimov


  SUBROUTINE read_permafrostmap(kjpindex,lalo,overburden,excess_ice,permafrost)
    
    INTEGER(i_std), INTENT(in) :: kjpindex
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in) :: lalo
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout) :: overburden
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout) :: excess_ice
    REAL(r_std), DIMENSION(kjpindex), INTENT(inout) :: permafrost
    
    INTEGER(i_std) :: il, ils, ip, ix, iy, imin, jmin, ier
    REAL(r_std) :: dlon, dlonmin, dlat, dlatmin
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: xx,yy, permafrost_file, continuous_file, discontinuous_file
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: sporadic_file, isolated_file, overburden_file, excess_ice_file
    REAL(r_std),ALLOCATABLE,DIMENSION(:) :: x,y
    REAL(r_std) :: lev(1), date, dt
    INTEGER(i_std) :: itau(1)
    REAL(r_std),DIMENSION(kjpindex) :: tref
    
    ! plus bas, on prend la temperature lue dans un fichier climato si celui-ci existe
    filename = 'NONE'
    CALL getin('PERMAFROST_MAP_FILE',filename)
    IF ( filename .EQ. "NONE" ) THEN
    ELSE
       CALL flininfo(filename,iml, jml, lml, tml, fid)
       ALLOCATE (yy(iml,jml), stat=ier)
       ALLOCATE (xx(iml,jml), stat=ier)
       ALLOCATE (x(iml),y(jml), stat=ier)
       ALLOCATE (continuous_file(iml,jml), stat=ier)
       ALLOCATE (discontinuous_file(iml,jml), stat=ier)
       ALLOCATE (sporadic_file(iml,jml), stat=ier)
       ALLOCATE (isolated_file(iml,jml), stat=ier)
       ALLOCATE (overburden_file(iml,jml), stat=ier)
       ALLOCATE (excess_ice_file(iml,jml), stat=ier)
       ALLOCATE (permafrost_file(iml,jml), stat=ier)
       CALL flinopen (filename, .FALSE., iml, jml, lml, &
            xx, yy, lev, tml, itau, date, dt, fid)
       CALL flinget (fid, 'continuous_permafrost', iml, jml, lml, tml, &
            1, 1, continuous_file)
       CALL flinget (fid, 'discontinuous_permafrost', iml, jml, lml, tml, &
            1, 1, discontinuous_file)
       CALL flinget (fid, 'sporadic_permafrost', iml, jml, lml, tml, &
            1, 1, sporadic_file)
       CALL flinget (fid, 'isolated_permafrost', iml, jml, lml, tml, &
            1, 1, isolated_file)
       CALL flinget (fid, 'thick_overburden', iml, jml, lml, tml, &
            1, 1, overburden_file)
       CALL flinget (fid, 'high_ground_ice_content', iml, jml, lml, tml, &
            1, 1, excess_ice_file)
       CALL flinclo (fid)
       ! On suppose que le fichier est regulier.
       ! Si ce n'est pas le cas, tant pis. Les temperatures seront mal
       ! initialisees et puis voila. De toute maniere, il faut avoir
       ! l'esprit mal tourne pour avoir l'idee de faire un fichier de
       ! climatologie avec une grille non reguliere.
       permafrost_file(:,:) = continuous_file + discontinuous_file + sporadic_file + isolated_file
       x(:) = xx(:,1)
       y(:) = yy(1,:)
       ! prendre la valeur la plus proche
       DO ip = 1, kjpindex
          dlonmin = HUGE(1.)
          DO ix = 1,iml
             dlon = MIN( ABS(lalo(ip,2)-x(ix)), ABS(lalo(ip,2)+360.-x(ix)), ABS(lalo(ip,2)-360.-x(ix)) )
             IF ( dlon .LT. dlonmin ) THEN
                imin = ix
                dlonmin = dlon
             ENDIF
          ENDDO
          dlatmin = HUGE(1.)
          DO iy = 1,jml
             dlat = ABS(lalo(ip,1)-y(iy))
             IF ( dlat .LT. dlatmin ) THEN
                jmin = iy
                dlatmin = dlat
             ENDIF
          ENDDO
          permafrost(ip) = permafrost_file(imin,jmin)
          overburden(ip) = overburden_file(imin,jmin)
          excess_ice(ip) = excess_ice_file(imin,jmin)
       ENDDO
       DEALLOCATE (yy)
       DEALLOCATE (xx)
       DEALLOCATE (x)
       DEALLOCATE (continuous_file)
       DEALLOCATE (discontinuous_file)
       DEALLOCATE (sporadic_file)
       DEALLOCATE (isolated_file)
       DEALLOCATE (overburden_file)
       DEALLOCATE (excess_ice_file)
       DEALLOCATE (permafrost_file)
    ENDIF
    WRITE(*,*) 'cdk: #points permafrost', sum(permafrost)
    WRITE(*,*) 'cdk: #points overburden', sum(overburden)
    WRITE(*,*) 'cdk: #points excess_ice', sum(excess_ice)
    
  END SUBROUTINE read_permafrostmap

 SUBROUTINE read_reftempfile(kjpindex,lalo,reftemp)
    
    INTEGER(i_std), INTENT(in) :: kjpindex
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in) :: lalo
    REAL(r_std), DIMENSION(kjpindex, ngrnd), INTENT(inout) :: reftemp

    
    INTEGER(i_std) :: il, ils, ip, ix, iy, imin, jmin, ier
    REAL(r_std) :: dlon, dlonmin, dlat, dlatmin
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: iml, jml, lml, tml, fid
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: xx,yy
    REAL(r_std),ALLOCATABLE,DIMENSION(:,:) :: reftemp_file
    REAL(r_std),ALLOCATABLE,DIMENSION(:) :: x,y
    REAL(r_std) :: lev(1), date, dt
    INTEGER(i_std) :: itau(1)
    REAL(r_std),DIMENSION(kjpindex) :: tref
    
    ! plus bas, on prend la temperature lue dans un fichier climato si celui-ci existe
    filename = 'reftemp.nc'
    CALL getin('REFTEMP_FILE',filename)

       CALL flininfo(filename,iml, jml, lml, tml, fid)
       ALLOCATE (yy(iml,jml), stat=ier)
       ALLOCATE (xx(iml,jml), stat=ier)
       ALLOCATE (x(iml),y(jml), stat=ier)
       ALLOCATE (reftemp_file(iml,jml), stat=ier)

       CALL flinopen (filename, .FALSE., iml, jml, lml, &
            xx, yy, lev, tml, itau, date, dt, fid)
       CALL flinget (fid, 'temperature', iml, jml, lml, tml, &
            1, 1, reftemp_file)
       CALL flinclo (fid)
       ! On suppose que le fichier est regulier.
       ! Si ce n'est pas le cas, tant pis. Les temperatures seront mal
       ! initialisees et puis voila. De toute maniere, il faut avoir
       ! l'esprit mal tourne pour avoir l'idee de faire un fichier de
       ! climatologie avec une grille non reguliere.
       x(:) = xx(:,1)
       y(:) = yy(1,:)
       ! prendre la valeur la plus proche
       DO ip = 1, kjpindex
          dlonmin = HUGE(1.)
          DO ix = 1,iml
             dlon = MIN( ABS(lalo(ip,2)-x(ix)), ABS(lalo(ip,2)+360.-x(ix)), ABS(lalo(ip,2)-360.-x(ix)) )
             IF ( dlon .LT. dlonmin ) THEN
                imin = ix
                dlonmin = dlon
             ENDIF
          ENDDO
          dlatmin = HUGE(1.)
          DO iy = 1,jml
             dlat = ABS(lalo(ip,1)-y(iy))
             IF ( dlat .LT. dlatmin ) THEN
                jmin = iy
                dlatmin = dlat
             ENDIF
          ENDDO
          reftemp(ip, :) = reftemp_file(imin,jmin)+273.15
       ENDDO
       DEALLOCATE (yy)
       DEALLOCATE (xx)
       DEALLOCATE (x)
       DEALLOCATE (reftemp_file)


    
  END SUBROUTINE read_reftempfile


  SUBROUTINE thermosoil_vert_axes( zz, zz_coef)

    ! interface description
    ! output fields
    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: zz
!!
    REAL(r_std), DIMENSION (ngrnd), INTENT(out)              :: zz_coef
    ! local declaration
    INTEGER(i_std)                                          ::  jg

    IF (long_print) WRITE (numout,*) 'entering thermosoil_vert_axes'

    !
    !     0. initialisation
    !
    cstgrnd=SQRT(one_day / pi)
    lskin = SQRT(so_cond / so_capa * one_day / pi)
    fz1 = 0.3_r_std * cstgrnd

    !
    !     1.  Computing the depth of the Temperature level, using a
    !         non dimentional variable x = z/lskin, lskin beeing
    !         the skin depth
    !

    !
    !     1.2 The undimentional depth is computed.
    !         ------------------------------------
    DO jg=1,ngrnd
      zz(jg)      = fz(REAL(jg,r_std) - undemi)
      zz_coef(jg) = fz(REAL(jg,r_std)-undemi+undemi)
    ENDDO
    !
    !     1.3 Converting to meters.
    !         --------------------
    DO jg=1,ngrnd
      zz(jg)      = zz(jg) /  cstgrnd * lskin
      zz_coef(jg) = zz_coef(jg) / cstgrnd * lskin
    ENDDO

    IF (long_print) WRITE (numout,*) ' thermosoil_vert_axes done'

  END SUBROUTINE thermosoil_vert_axes


!!
!================================================================================================================================
!! FUNCTION     : thermosoil_levels
!!
!>\BRIEF          Depth of nodes for the thermal layers in meters.
!!
!! DESCRIPTION  : Calculate and return the depth in meters of the nodes of the
!soil layers. This calculation is the same
!!                as done in thermosoil_var_init for zz. See thermosoil_var_init
!for more details.
!!
!! RECENT CHANGE(S) : None
!!
!! RETURN VALUE : Vector of soil depth for the nodes in meters
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  FUNCTION thermosoil_levels() RESULT (zz_out)

    !! 0.1 Return variable

    REAL(r_std), DIMENSION (ngrnd)  :: zz_out      !! Depth of soil layers in meters

    !! 0.2 Local variables
    INTEGER(i_std)                  :: jg
    REAL(r_std)                     :: so_capa
    REAL(r_std)                     :: so_cond

!_
!================================================================================================================================

    !! 1. Define some parameters
    so_capa = (so_capa_dry + so_capa_wet)/deux
    so_cond = (so_cond_dry + so_cond_wet)/deux

    cstgrnd=SQRT(one_day / pi)
    lskin = SQRT(so_cond / so_capa * one_day / pi)

    !! Parameters needed by fz function
    fz1 = 0.3_r_std * cstgrnd
    !zalph = deux

    !!  2. Get adimentional depth of the numerical nodes
    DO jg=1,ngrnd
       zz_out(jg) = fz(REAL(jg,r_std) - undemi)
    ENDDO

    !! 3. Convert to meters
    DO jg=1,ngrnd
       zz_out(jg) = zz_out(jg) /  cstgrnd * lskin
    END DO

  END FUNCTION thermosoil_levels

END MODULE thermosoil
