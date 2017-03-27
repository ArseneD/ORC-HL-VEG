! =================================================================================================================================
! MODULE       : constantes_var
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        constantes_var module contains most constantes like pi, Earth radius, etc...
!!              and all externalized parameters except pft-dependent constants.
!!
!!\n DESCRIPTION: This module contains most constantes and the externalized parameters of ORCHIDEE which 
!!                are not pft-dependent.\n
!!                In this module, you can set the flag diag_qsat in order to detect the pixel where the
!!                temperature is out of range (look qsatcalc and dev_qsatcalc in qsat_moisture.f90).\n
!!                The Earth radius is approximated by the Equatorial radius.The Earth's equatorial radius a,
!!                or semi-major axis, is the distance from its center to the equator and equals 6,378.1370 km.
!!                The equatorial radius is often used to compare Earth with other planets.\n
!!                The meridional mean is well approximated by the semicubic mean of the two axe yielding 
!!                6367.4491 km or less accurately by the quadratic mean of the two axes about 6,367.454 km
!!                or even just the mean of the two axes about 6,367.445 km.\n
!!                This module is already USE in module constantes. Therefor no need to USE it seperatly except
!!                if the subroutines in module constantes are not needed.\n
!!                
!! RECENT CHANGE(S):
!!
!! REFERENCE(S)	: 
!! - Louis, Jean-Francois (1979), A parametric model of vertical eddy fluxes in the atmosphere. 
!! Boundary Layer Meteorology, 187-202.\n
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2014-09-04 14:46:14 +0200 (Thu, 04 Sep 2014) $
!! $Revision: 2282 $
!! \n
!_ ================================================================================================================================

MODULE constantes_var

  USE defprec

  IMPLICIT NONE
!-

                         !-----------------------!
                         !  ORCHIDEE CONSTANTS   !
                         !-----------------------!

  !
  ! FLAGS 
  !
  TYPE control_type
    LOGICAL :: river_routing      !! activate river routing (true/false)
    LOGICAL :: hydrol_cwrr        !! activate 11 layers hydrolgy model (true/false)
    LOGICAL :: do_floodplains
    LOGICAL :: do_irrigation
    LOGICAL :: ok_sechiba         !! activate physic of the model (true/false)
    LOGICAL :: ok_co2             !! activate photosynthesis (true/false)
    LOGICAL :: ok_stomate         !! activate carbon cycle (true/false)
    LOGICAL :: ok_dgvm            !! activate dynamic vegetation (true/false)
    LOGICAL :: stomate_watchout   !! activate the creation of restart files for STOMATE even if STOMATE is not activated (true/false)
    LOGICAL :: ok_pheno           !! activate the calculation of lai using stomate rather than a prescription (true/false)
    LOGICAL :: do_land_use
    LOGICAL :: ok_inca            !! activate biogenic volatile organic coumpounds ? (true/false)
    LOGICAL :: ok_leafage         !! activate leafage? (true/false)
    LOGICAL :: ok_radcanopy       !! use canopy radiative transfer model (true/false)
    LOGICAL :: ok_multilayer      !! use canopy radiative transfer model with multi-layers (true/false)
    LOGICAL :: ok_pulse_NOx       !! calculate NOx emissions with pulse (true/false)
    LOGICAL :: ok_bbgfertil_NOx   !! calculate NOx emissions with bbg fertilizing effect (true/false)
    LOGICAL :: ok_cropsfertil_NOx !! calculate NOx emissions with fertilizers use (true/false)
  END TYPE control_type

  !-
  TYPE(control_type), SAVE :: control  !! Flags that (de)activate parts of the model
!$OMP THREADPRIVATE(control)
  LOGICAL, SAVE :: OFF_LINE_MODE = .FALSE.  !! ORCHIDEE detects if it is coupled with a GCM or 
                                            !! just use with one driver in OFF-LINE. (true/false)
!$OMP THREADPRIVATE(OFF_LINE_MODE)
  CHARACTER(LEN=80), SAVE     :: restname_in       = 'NONE'                 !! Input Restart files name for Sechiba component  
!$OMP THREADPRIVATE(restname_in)
  CHARACTER(LEN=80), SAVE     :: restname_out      = 'sechiba_rest_out.nc'  !! Output Restart files name for Sechiba component
!$OMP THREADPRIVATE(restname_out)
  CHARACTER(LEN=80), SAVE     :: stom_restname_in  = 'NONE'                 !! Input Restart files name for Stomate component
!$OMP THREADPRIVATE(stom_restname_in)
  CHARACTER(LEN=80), SAVE     :: stom_restname_out = 'stomate_rest_out.nc'  !! Output Restart files name for Stomate component
!$OMP THREADPRIVATE(stom_restname_out)

  !
  ! TIME
  !
  REAL(r_std), SAVE :: one_day  !! One day in seconds (s)
!$OMP THREADPRIVATE(one_day)
  REAL(r_std), SAVE :: one_year !! One year in seconds (s)
!$OMP THREADPRIVATE(one_year)
  REAL(r_std), PARAMETER :: one_hour = 3600.0  !! One hour in seconds (s)

  ! TIME STEP
  REAL(r_std)            :: dt_sechiba         !! Time step for in sechiba
!$OMP THREADPRIVATE(dt_sechiba)

  !
  ! SPECIAL VALUES 
  !
  INTEGER(i_std), PARAMETER :: undef_int = 999999999     !! undef integer for integer arrays (unitless)
  !-
  REAL(r_std), SAVE :: val_exp = 999999.                 !! Specific value if no restart value  (unitless)
!$OMP THREADPRIVATE(val_exp)
  REAL(r_std), PARAMETER :: undef = -9999.               !! Special value for stomate (unitless)
  !-
  REAL(r_std), PARAMETER :: min_sechiba = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL(r_std), PARAMETER :: undef_sechiba = 1.E+20_r_std !! The undef value used in SECHIBA (unitless)
  !-
  REAL(r_std), PARAMETER :: min_stomate = 1.E-8_r_std    !! Epsilon to detect a near zero floating point (unitless)
  REAL(r_std), PARAMETER :: large_value = 1.E33_r_std    !! some large value (for stomate) (unitless)


  !
  !  DIMENSIONING AND INDICES PARAMETERS  
  !
  INTEGER(i_std), PARAMETER :: ibare_sechiba = 1 !! Index for bare soil in Sechiba (unitless)
  INTEGER(i_std), PARAMETER :: ivis = 1          !! index for albedo in visible range (unitless)
  INTEGER(i_std), PARAMETER :: inir = 2          !! index for albeod i near-infrared range (unitless) 
  INTEGER(i_std), PARAMETER :: nnobio = 1        !! Number of other surface types: land ice (lakes,cities, ...) (unitless)
  INTEGER(i_std), PARAMETER :: iice = 1          !! Index for land ice (see nnobio) (unitless)
  !-
  !! Soil
  INTEGER(i_std), PARAMETER :: classnb = 9       !! Levels of soil colour classification (unitless)
  !-
  INTEGER(i_std), PARAMETER :: nleafages = 4     !! leaf age discretisation ( 1 = no discretisation )(unitless)
  !-
  !! litter fractions: indices (unitless)
  INTEGER(i_std), PARAMETER :: ileaf = 1         !! Index for leaf compartment (unitless)
  INTEGER(i_std), PARAMETER :: isapabove = 2     !! Index for sapwood above compartment (unitless)
  INTEGER(i_std), PARAMETER :: isapbelow = 3     !! Index for sapwood below compartment (unitless)
  INTEGER(i_std), PARAMETER :: iheartabove = 4   !! Index for heartwood above compartment (unitless)
  INTEGER(i_std), PARAMETER :: iheartbelow = 5   !! Index for heartwood below compartment (unitless)
  INTEGER(i_std), PARAMETER :: iroot = 6         !! Index for roots compartment (unitless)
  INTEGER(i_std), PARAMETER :: ifruit = 7        !! Index for fruits compartment (unitless)
  INTEGER(i_std), PARAMETER :: icarbres = 8      !! Index for reserve compartment (unitless)
  INTEGER(i_std), PARAMETER :: nparts = 8        !! Number of biomass compartments (unitless)
  !-
  !! indices for assimilation parameters 
  INTEGER(i_std), PARAMETER :: ivcmax = 1        !! Index for vcmax (assimilation parameters) (unitless)
  INTEGER(i_std), PARAMETER :: npco2 = 1         !! Number of assimilation parameters (unitless)
  !-
  !! trees and litter: indices for the parts of heart-
  !! and sapwood above and below the ground 
  INTEGER(i_std), PARAMETER :: iabove = 1       !! Index for above part (unitless)
  INTEGER(i_std), PARAMETER :: ibelow = 2       !! Index for below part (unitless)
  INTEGER(i_std), PARAMETER :: nlevs = 2        !! Number of levels for trees and litter (unitless)
  !-
  !! litter: indices for metabolic and structural part
  INTEGER(i_std), PARAMETER :: imetabolic = 1   !! Index for metabolic litter (unitless)
  INTEGER(i_std), PARAMETER :: istructural = 2  !! Index for structural litter (unitless)
  INTEGER(i_std), PARAMETER :: nlitt = 2        !! Number of levels for litter compartments (unitless)
  !-
  !! carbon pools: indices
  INTEGER(i_std), PARAMETER :: iactive = 1      !! Index for active carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: islow = 2        !! Index for slow carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ipassive = 3     !! Index for passive carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ncarb = 3        !! Number of soil carbon pools (unitless)
  !-
  !! For isotopes and nitrogen
  INTEGER(i_std), PARAMETER :: nelements = 1    !! Number of isotopes considered
  INTEGER(i_std), PARAMETER :: icarbon = 1      !! Index for carbon 
  !
  !! Indices used for analytical spin-up
  INTEGER(i_std), PARAMETER :: nbpools = 7              !! Total number of carbon pools (unitless)
  INTEGER(i_std), PARAMETER :: istructural_above = 1    !! Index for structural litter above (unitless)
  INTEGER(i_std), PARAMETER :: istructural_below = 2    !! Index for structural litter below (unitless)
  INTEGER(i_std), PARAMETER :: imetabolic_above = 3     !! Index for metabolic litter above (unitless)
  INTEGER(i_std), PARAMETER :: imetabolic_below = 4     !! Index for metabolic litter below (unitless)
  INTEGER(i_std), PARAMETER :: iactive_pool = 5         !! Index for active carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: islow_pool   = 6         !! Index for slow carbon pool (unitless)
  INTEGER(i_std), PARAMETER :: ipassive_pool = 7        !! Index for passive carbon pool (unitless)


  !
  ! NUMERICAL AND PHYSICS CONSTANTS
  !
  !

  !-
  ! 1. Mathematical and numerical constants
  !-
  REAL(r_std), PARAMETER :: pi = 3.141592653589793238   !! pi souce : http://mathworld.wolfram.com/Pi.html (unitless)
  REAL(r_std), PARAMETER :: euler = 2.71828182845904523 !! e source : http://mathworld.wolfram.com/e.html (unitless)
  REAL(r_std), PARAMETER :: zero = 0._r_std             !! Numerical constant set to 0 (unitless)
  REAL(r_std), PARAMETER :: undemi = 0.5_r_std          !! Numerical constant set to 1/2 (unitless)
  REAL(r_std), PARAMETER :: un = 1._r_std               !! Numerical constant set to 1 (unitless)
  REAL(r_std), PARAMETER :: moins_un = -1._r_std        !! Numerical constant set to -1 (unitless)
  REAL(r_std), PARAMETER :: deux = 2._r_std             !! Numerical constant set to 2 (unitless)
  REAL(r_std), PARAMETER :: trois = 3._r_std            !! Numerical constant set to 3 (unitless)
  REAL(r_std), PARAMETER :: quatre = 4._r_std           !! Numerical constant set to 4 (unitless)
  REAL(r_std), PARAMETER :: cinq = 5._r_std             !![DISPENSABLE] Numerical constant set to 5 (unitless)
  REAL(r_std), PARAMETER :: six = 6._r_std              !![DISPENSABLE] Numerical constant set to 6 (unitless)
  REAL(r_std), PARAMETER :: huit = 8._r_std             !! Numerical constant set to 8 (unitless)
  REAL(r_std), PARAMETER :: mille = 1000._r_std         !! Numerical constant set to 1000 (unitless)

  !-
  ! 2 . Physics
  !-
  REAL(r_std), PARAMETER :: R_Earth = 6378000.              !! radius of the Earth : Earth radius ~= Equatorial radius (m)
  REAL(r_std), PARAMETER :: mincos  = 0.0001                !! Minimum cosine value used for interpolation (unitless) 
  REAL(r_std), PARAMETER :: pb_std = 1013.                  !! standard pressure (hPa)
  REAL(r_std), PARAMETER :: ZeroCelsius = 273.15            !! 0 degre Celsius in degre Kelvin (K)
  REAL(r_std), PARAMETER :: tp_00 = 273.15                  !! 0 degre Celsius in degre Kelvin (K)
  REAL(r_std), PARAMETER :: chalsu0 = 2.8345E06             !! Latent heat of sublimation (J.kg^{-1})
  REAL(r_std), PARAMETER :: chalev0 = 2.5008E06             !! Latent heat of evaporation (J.kg^{-1}) 
  REAL(r_std), PARAMETER :: chalfu0 = chalsu0-chalev0       !! Latent heat of fusion (J.kg^{-1}) 
  REAL(r_std), PARAMETER :: c_stefan = 5.6697E-8            !! Stefan-Boltzman constant (W.m^{-2}.K^{-4})
  REAL(r_std), PARAMETER :: cp_air = 1004.675               !! Specific heat of dry air (J.kg^{-1}.K^{-1}) 
  REAL(r_std), PARAMETER :: cte_molr = 287.05               !! Specific constant of dry air (kg.mol^{-1}) 
  REAL(r_std), PARAMETER :: kappa = cte_molr/cp_air         !! Kappa : ratio between specific constant and specific heat 
                                                            !! of dry air (unitless)
  REAL(r_std), PARAMETER :: msmlr_air = 28.964E-03          !! Molecular weight of dry air (kg.mol^{-1})
  REAL(r_std), PARAMETER :: msmlr_h2o = 18.02E-03           !! Molecular weight of water vapor (kg.mol^{-1}) 
  REAL(r_std), PARAMETER :: cp_h2o = &                      !! Specific heat of water vapor (J.kg^{-1}.K^{-1}) 
       & cp_air*(quatre*msmlr_air)/( 3.5_r_std*msmlr_h2o) 
  REAL(r_std), PARAMETER :: cte_molr_h2o = cte_molr/quatre  !! Specific constant of water vapor (J.kg^{-1}.K^{-1}) 
  REAL(r_std), PARAMETER :: retv = msmlr_air/msmlr_h2o-un   !! Ratio between molecular weight of dry air and water 
                                                            !! vapor minus 1(unitless)  
  REAL(r_std), PARAMETER :: rvtmp2 = cp_h2o/cp_air-un       !! Ratio between specific heat of water vapor and dry air
                                                            !! minus 1 (unitless)
  REAL(r_std), PARAMETER :: cepdu2 = (0.1_r_std)**2         !! Squared wind shear (m^2.s^{-2}) 
  REAL(r_std), PARAMETER :: ct_karman = 0.35_r_std          !! Van Karmann Constant (unitless)
  REAL(r_std), PARAMETER :: cte_grav = 9.80665_r_std        !! Acceleration of the gravity (m.s^{-2})
  REAL(r_std), PARAMETER :: pa_par_hpa = 100._r_std         !! Transform pascal into hectopascal (unitless)
  REAL(r_std), PARAMETER :: RR = 8.314                      !! Ideal gas constant (J.mol^{-1}.K^{-1})
  REAL(r_std), PARAMETER :: Sct = 1370.                     !! Solar constant (W.m^{-2}) 


  !-
  ! 3. Climatic constants
  !-
  !! Constantes of the Louis scheme 
  REAL(r_std), SAVE :: cb = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
  REAL(r_std), SAVE :: cc = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
  REAL(r_std), SAVE :: cd = 5._r_std              !! Constant of the Louis scheme (unitless);
                                                  !! reference to Louis (1979)
  !-
  REAL(r_std), SAVE :: rayt_cste = 125.           !! Constant in the computation of surface resistance (W.m^{-2})
  REAL(r_std), SAVE :: defc_plus = 23.E-3         !! Constant in the computation of surface resistance (K.W^{-1})
  REAL(r_std), SAVE :: defc_mult = 1.5            !! Constant in the computation of surface resistance (K.W^{-1})

  !-
  ! 4. Soil thermodynamics constants
  !-
  ! Look at constantes_soil.f90


  !
  ! OPTIONAL PARTS OF THE MODEL
  !
  LOGICAL, SAVE     :: long_print = .false.       !! To set for more printing
!$OMP THREADPRIVATE(long_print)
  LOGICAL,PARAMETER :: diag_qsat = .TRUE.         !! One of the most frequent problems is a temperature out of range
                                                  !! we provide here a way to catch that in the calling procedure. 
                                                  !! (from Jan Polcher)(true/false) 
  LOGICAL, SAVE     :: almaoutput                 !! Selects the type of output for the model.(true/false)
                                                  !! Value is read from run.def in intersurf_history
!$OMP THREADPRIVATE(almaoutput)

  !
  ! DIVERSE
  !
  CHARACTER(LEN=100), SAVE :: stomate_forcing_name='NONE'  !! NV080800 Name of STOMATE forcing file (unitless)
                                                           ! Compatibility with Nicolas Viovy driver.
!$OMP THREADPRIVATE(stomate_forcing_name)
  CHARACTER(LEN=100), SAVE :: stomate_Cforcing_name='NONE' !! NV080800 Name of soil forcing file (unitless)
                                                           ! Compatibility with Nicolas Viovy driver.
!$OMP THREADPRIVATE(stomate_Cforcing_name)
  INTEGER(i_std), SAVE :: forcing_id                 !! Index of the forcing file (unitless)
!$OMP THREADPRIVATE(forcing_id)




                         !------------------------!
                         !  SECHIBA PARAMETERS    !
                         !------------------------!
 

  !
  ! GLOBAL PARAMETERS   
  !
  REAL(r_std), SAVE :: min_wind = 0.1      !! The minimum wind (m.s^{-1})
!$OMP THREADPRIVATE(min_wind)
  REAL(r_std), SAVE :: snowcri = 1.5       !! Sets the amount above which only sublimation occures (kg.m^{-2})
!$OMP THREADPRIVATE(snowcri)


  !
  ! FLAGS ACTIVATING SUB-MODELS
  !
  LOGICAL, SAVE :: treat_expansion = .FALSE.   !! Do we treat PFT expansion across a grid point after introduction? (true/false)
!$OMP THREADPRIVATE(treat_expansion)
  LOGICAL, SAVE :: ok_herbivores = .FALSE.     !! flag to activate herbivores (true/false)
!$OMP THREADPRIVATE(ok_herbivores)
  LOGICAL, SAVE :: harvest_agri = .TRUE.       !! flag to harvest aboveground biomass from agricultural PFTs)(true/false)
!$OMP THREADPRIVATE(harvest_agri)
  LOGICAL, SAVE :: lpj_gap_const_mort = .TRUE. !! constant moratlity (true/false)
!$OMP THREADPRIVATE(lpj_gap_const_mort)
  LOGICAL, SAVE :: disable_fire = .FALSE.      !! flag that disable fire (true/false)
!$OMP THREADPRIVATE(disable_fire)
  LOGICAL, SAVE :: spinup_analytic = .FALSE.   !! Flag to activate analytical resolution for spinup (true/false)
!$OMP THREADPRIVATE(spinup_analytic)
  LOGICAL, SAVE :: ok_explicitsnow             !! Flag to activate explicit snow scheme instead of default snow scheme
!$OMP THREADPRIVATE(ok_explicitsnow)
  LOGICAL, SAVE :: ok_pc                       !! Flag to activate permafrost carbon (vertical carbon and soil carbon thermal insulation)
!$OMP THREADPRIVATE(ok_pc)

  !
  ! CONFIGURATION VEGETATION
  !
  LOGICAL, SAVE :: agriculture = .TRUE.    !! allow agricultural PFTs (true/false)
!$OMP THREADPRIVATE(agriculture)
  LOGICAL, SAVE :: impveg = .FALSE.        !! Impose vegetation ? (true/false)
!$OMP THREADPRIVATE(impveg)
  LOGICAL, SAVE :: impsoilt = .FALSE.      !! Impose soil ? (true/false)
!$OMP THREADPRIVATE(impsoilt)
  LOGICAL, SAVE :: lcchange = .FALSE.      !! Land cover change flag (true/false)
!$OMP THREADPRIVATE(lcchange)
  LOGICAL, SAVE :: read_lai = .FALSE.      !! Flag to read a map of LAI if STOMATE is not activated (true/false)
!$OMP THREADPRIVATE(read_lai)
  LOGICAL, SAVE :: old_lai = .FALSE.       !! Flag for the old LAI map interpolation (SHOULD BE DROPED ??)(true/false)
!$OMP THREADPRIVATE(old_lai)
  LOGICAL, SAVE :: old_veget = .FALSE.     !! Flag to use the old vegetation Map interpolation (SHOULD BE DROPED ?)(true/false)
!$OMP THREADPRIVATE(old_veget)
  LOGICAL, SAVE :: land_use = .TRUE.       !! flag to account or not for Land Use  (true/false)
!$OMP THREADPRIVATE(land_use)
  LOGICAL, SAVE :: veget_reinit = .TRUE.   !! To change LAND USE file in a run. (true/false)
!$OMP THREADPRIVATE(veget_reinit)
  LOGICAL, SAVE :: shrubs_like_trees = .FALSE.  !! Arsene 03-08-2015 - Add shrubs_like_trees !! Flag to use equation closer trees for shrubs  
!$OMP THREADPRIVATE(shrubs_like_trees)          !! Arsene 03-08-2015 - Add shrubs_like_trees
  LOGICAL, SAVE :: new_moist_func = .TRUE.  !! Arsene 29-12-2015 - ADD for LUT: new litter and soil_carbon moist dependence decomposition 
!$OMP THREADPRIVATE(new_moist_func)         !! Arsene 29-12-2015
  LOGICAL, SAVE :: new_litter_discret = .TRUE.     !! Arsene 29-12-2015 - ADD for new soil discretisation for stomate_litter  
!$OMP THREADPRIVATE(new_litter_discret)            !! Arsene 29-12-2015

  !
  ! PARAMETERS USED BY BOTH HYDROLOGY MODELS
  !
  REAL(r_std), SAVE :: max_snow_age = 50._r_std !! Maximum period of snow aging (days)
!$OMP THREADPRIVATE(max_snow_age)
  REAL(r_std), SAVE :: snow_trans = 0.3_r_std   !! Transformation time constant for snow (m)
!$OMP THREADPRIVATE(snow_trans)
  REAL(r_std), SAVE :: sneige                   !! Lower limit of snow amount (kg.m^{-2})
!$OMP THREADPRIVATE(sneige)
  REAL(r_std), SAVE :: maxmass_snow = 3000.     !! The maximum mass of snow (kg.m^{-2})
!$OMP THREADPRIVATE(maxmass_snow)

  !! Heat capacity
  REAL(r_std), PARAMETER :: capa_ice = 2.228*1.E3       !! Heat capacity of ice (J/kg/K)
  REAL(r_std), SAVE      :: so_capa_ice = 2.11e6        !! Heat capacity of saturated frozen soil (J/K/m3)  !! Arsene 19-12-2014 Value add from Tao
!$OMP THREADPRIVATE(so_capa_ice)
  REAL(r_std), PARAMETER :: rho_water = 1000.           !! Density of water (kg/m3)
  REAL(r_std), PARAMETER :: rho_ice = 920.              !! Density of ice (kg/m3)

  !! Thermal conductivities
  REAL(r_std), PARAMETER :: cond_water = 0.6            !! Thermal conductivity of liquid water (W/m/K)
  REAL(r_std), PARAMETER :: cond_ice = 2.2              !! Thermal conductivity of ice (W/m/K)
  REAL(r_std), PARAMETER :: cond_solid = 2.32           !! Thermal conductivity of mineral soil particles (W/m/K)

  !! Time constant of long-term soil humidity (s) 
  REAL(r_std), PARAMETER :: lhf = 0.3336*1.E6           !! Latent heat of fusion (J/kg)

  INTEGER(i_std), PARAMETER :: nsnow=3                  !! Number of levels in the snow for explicit snow scheme   
  REAL(r_std), PARAMETER    :: XMD    = 28.9644E-3 
  REAL(r_std), PARAMETER    :: XBOLTZ      = 1.380658E-23 
  REAL(r_std), PARAMETER    :: XAVOGADRO   = 6.0221367E+23 
  REAL(r_std), PARAMETER    :: XRD    = XAVOGADRO * XBOLTZ / XMD 
  REAL(r_std), PARAMETER    :: XCPD   = 7.* XRD /2. 
  REAL(r_std), PARAMETER    :: phigeoth = 0.057 ! 0. DKtest 
  REAL(r_std), PARAMETER    :: thick_min_snow = .01 

  !! The maximum snow density and water holding characterisicts 
  REAL(r_std), SAVE         :: xrhosmax = 750.  ! (kg m-3) 
  REAL(r_std), SAVE         :: xwsnowholdmax1   = 0.03  ! (-) 
  REAL(r_std), SAVE         :: xwsnowholdmax2   = 0.10  ! (-) 
  REAL(r_std), SAVE         :: xsnowrhohold     = 200.0 ! (kg/m3) 
  REAL(r_std), SAVE         :: xrhosmin = 50. 
  REAL(r_std), PARAMETER    :: xci = 2.106e+3 
  REAL(r_std), PARAMETER    :: xrv = 6.0221367e+23 * 1.380658e-23 /18.0153e-3 

  !! ISBA-ES Critical snow depth at which snow grid thicknesses constant 
  REAL(r_std), PARAMETER    :: xsnowcritd = 0.03  ! (m) 
  
  !! ISBA-ES CROCUS (Pahaut 1976): snowfall density coefficients: 
  REAL(r_std), PARAMETER       :: snowfall_a_sn = 109.0  !! (kg/m3) 
  REAL(r_std), PARAMETER       :: snowfall_b_sn =   6.0  !! (kg/m3/K) 
  REAL(r_std), PARAMETER       :: snowfall_c_sn =  26.0  !! [kg/(m7/2 s1/2)] 

  REAL(r_std), PARAMETER       :: dgrain_new_max=  2.0e-4!! (m) : Maximum grain size of new snowfall 
  
  !! Minimum snow layer thickness for thermal calculations. Used to prevent 
  !! numerical problems as snow becomes vanishingly thin. 
  REAL(r_std), PARAMETER                :: psnowdzmin = .0001   ! m 
  REAL(r_std), PARAMETER                :: xsnowdmin = .000001  ! m 

  REAL(r_std), PARAMETER                :: ph2o = 1000.         !! Water density [kg/m3] 
  
  ! ISBA-ES Thermal conductivity coefficients from Anderson (1976): 
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002) 
  REAL(r_std), SAVE                     :: ZSNOWTHRMCOND1 = 0.02    ! [W/m/K] 
  REAL(r_std), SAVE                     :: ZSNOWTHRMCOND2 = 2.5E-6  ! [W m5/(kg2 K)] 
  
  ! ISBA-ES Thermal conductivity: Implicit vapor diffn effects 
  ! (sig only for new snow OR high altitudes) 
  ! from Sun et al. (1999): based on data from Jordan (1991) 
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002) 
  ! 
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_AVAP  = -0.06023 ! (W/m/K) 
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_BVAP  = -2.5425  ! (W/m) 
  REAL(r_std), SAVE                       :: ZSNOWTHRMCOND_CVAP  = -289.99  ! (K) 
  
  REAL(r_std),SAVE :: xansmax = 0.85      !! Maxmimum snow albedo
  REAL(r_std),SAVE :: xansmin = 0.50      !! Miniumum snow albedo
  REAL(r_std),SAVE :: xans_todry = 0.008  !! Albedo decay rate for dry snow
  REAL(r_std),SAVE :: xans_t = 0.240      !! Albedo decay rate for wet snow

  ! ISBA-ES Thermal conductivity coefficients from Anderson (1976):
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
  REAL(r_std), PARAMETER                  :: XP00 = 1.E5

  ! ISBA-ES Thermal conductivity: Implicit vapor diffn effects
  ! (sig only for new snow OR high altitudes)
  ! from Sun et al. (1999): based on data from Jordan (1991)
  ! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
  !

!! Arsene 23-09-2014 Change for snowpack - START
!  REAL(r_std), SAVE          :: ZSNOWCMPCT_RHOD  = 150.0        !! (kg/m3)
!  REAL(r_std), SAVE          :: ZSNOWCMPCT_ACM   = 2.8e-6       !! (1/s)
!  REAL(r_std), SAVE          :: ZSNOWCMPCT_BCM   = 0.04         !! (1/K) 
!  REAL(r_std), SAVE          :: ZSNOWCMPCT_CCM   = 460.         !! (m3/kg)
!  REAL(r_std), SAVE          :: ZSNOWCMPCT_V0    = 3.7e7        !! (Pa/s) 
!  REAL(r_std), SAVE          :: ZSNOWCMPCT_VT    = 0.081        !! (1/K)
!  REAL(r_std), SAVE          :: ZSNOWCMPCT_VR    = 0.018        !! (m3/kg)

!! Dimension = 3, for strate herbacée, arbustive et arborescente 
  REAL(r_std), SAVE, DIMENSION(3)   :: ZSNOWCMPCT_RHOD  = (/  150.0,  150.0,  150.0 /)        ! (kg/m3)

  REAL(r_std), SAVE, DIMENSION(3)   :: ZSNOWCMPCT_ACM   = (/ 1.4e-6, 4.2e-6, 4.2e-6 /)        ! (1/s)
  REAL(r_std), SAVE, DIMENSION(3)   :: ZSNOWCMPCT_BCM   = (/   0.02,   0.06,   0.06 /)        ! (1/K) 
  REAL(r_std), SAVE, DIMENSION(3)   :: ZSNOWCMPCT_CCM   = (/   230.,   690.,   690. /)        ! (m3/kg)

  REAL(r_std), SAVE, DIMENSION(3)   :: ZSNOWCMPCT_V0    = (/ 1.85e7, 5.55e7, 5.55e7 /)        ! (Pa/s)
  REAL(r_std), SAVE, DIMENSION(3)   :: ZSNOWCMPCT_VT    = (/ 0.0405,   0.12,   0.12 /)        ! (1/K)
  REAL(r_std), SAVE, DIMENSION(3)   :: ZSNOWCMPCT_VR    = (/  0.009,  0.027,  0.027 /)        ! (m3/kg)

!! Arsene 23-09-2014 Change for snowpack - END

  !
  ! BVOC : Biogenic activity  for each age class
  !
  REAL(r_std), SAVE, DIMENSION(nleafages) :: iso_activity = (/0.5, 1.5, 1.5, 0.5/)     !! Biogenic activity for each 
                                                                                       !! age class : isoprene (unitless)
!$OMP THREADPRIVATE(iso_activity)
  REAL(r_std), SAVE, DIMENSION(nleafages) :: methanol_activity = (/1., 1., 0.5, 0.5/)  !! Biogenic activity for each
                                                                                       !! age class : methanol (unnitless)
!$OMP THREADPRIVATE(methanol_activity)

  !
  ! condveg.f90
  !

  ! 1. Scalar

  ! 1.1 Flags used inside the module

  LOGICAL, SAVE :: alb_bare_model = .FALSE. !! Switch for choosing values of bare soil 
                                            !! albedo (see header of subroutine)
                                            !! (true/false)
!$OMP THREADPRIVATE(alb_bare_model)
!! Arsene 16-07-2016 - Add new Albedo - START
  LOGICAL, SAVE :: alb_bg_modis = .FALSE.   !! Switch for choosing values of bare soil 
                                            !! albedo read from file
                                            !! (true/false)
!$OMP THREADPRIVATE(alb_bg_modis)
!! Arsene 16-07-2016 - Add new Albedo - END
  LOGICAL, SAVE :: impaze = .FALSE.         !! Switch for choosing surface parameters
                                            !! (see header of subroutine).  
                                            !! (true/false)
!$OMP THREADPRIVATE(impaze)
  LOGICAL, SAVE :: z0cdrag_ave = .TRUE.     !! Chooses between two methods to calculate the 
                                            !! grid average of the roughness (see header of subroutine)   
                                            !! (true/false)
!$OMP THREADPRIVATE(z0cdrag_ave)
  ! 1.2 Others 

  REAL(r_std), SAVE :: z0_over_height = un/16.           !! Factor to calculate roughness height from 
                                                         !! vegetation height (unitless)   
!$OMP THREADPRIVATE(z0_over_height)
  REAL(r_std), SAVE :: height_displacement = 0.75        !! Factor to calculate the zero-plane displacement
                                                         !! height from vegetation height (m)
!$OMP THREADPRIVATE(height_displacement)
  REAL(r_std), SAVE :: z0_bare = 0.01                    !! bare soil roughness length (m)
!$OMP THREADPRIVATE(z0_bare)
  REAL(r_std), SAVE :: z0_ice = 0.001                    !! ice roughness length (m)
!$OMP THREADPRIVATE(z0_ice)
  REAL(r_std), SAVE :: tcst_snowa = 5.0                  !! Time constant of the albedo decay of snow (days)
!$OMP THREADPRIVATE(tcst_snowa)
  REAL(r_std), SAVE :: snowcri_alb = 10.                 !! Critical value for computation of snow albedo (cm)
!$OMP THREADPRIVATE(snowcri_alb)
!! Arsene 22-11-2016 - Add for frac_snow_veg such as Boone 2002 - START
  REAL(r_std), SAVE :: snowcrit_z0alb1 = 5.              !! Critical parameter value for computation of snow fraction on vegetation (-)
!$OMP THREADPRIVATE(snowcrit_z0alb1)
  REAL(r_std), SAVE :: snowcrit_z0alb2 = 1.              !! Critical parameter value for computation of snow fraction on vegetation (-)
!$OMP THREADPRIVATE(snowcrit_z0alb2)
!! Arsene 22-11-2016 - Add for frac_snow_veg such as Boone 2002 - END 
  REAL(r_std), SAVE :: fixed_snow_albedo = undef_sechiba !! To choose a fixed snow albedo value (unitless)
!$OMP THREADPRIVATE(fixed_snow_albedo)
  REAL(r_std), SAVE :: z0_scal = 0.15                    !! Surface roughness height imposed (m)
!$OMP THREADPRIVATE(z0_scal)
  REAL(r_std), SAVE :: roughheight_scal = zero           !! Effective roughness Height depending on zero-plane 
                                                         !! displacement height (m) (imposed)
!$OMP THREADPRIVATE(roughheight_scal)
  REAL(r_std), SAVE :: emis_scal = 1.0                   !! Surface emissivity imposed (unitless)
!$OMP THREADPRIVATE(emis_scal)
  ! 2. Arrays

  REAL(r_std), SAVE, DIMENSION(2) :: alb_deadleaf = (/ .12, .35/)    !! albedo of dead leaves, VIS+NIR (unitless)
!$OMP THREADPRIVATE(alb_deadleaf)
  REAL(r_std), SAVE, DIMENSION(2) :: alb_ice = (/ .60, .20/)         !! albedo of ice, VIS+NIR (unitless)
!$OMP THREADPRIVATE(alb_ice)
  REAL(r_std), SAVE, DIMENSION(2) :: albedo_scal = (/ 0.25, 0.25 /)  !! Albedo values for visible and near-infrared 
                                                                     !! used imposed (unitless) 
!$OMP THREADPRIVATE(albedo_scal)
  REAL(r_std) , SAVE, DIMENSION(classnb) :: vis_dry = (/0.24,&
       &0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.27/)  !! Soil albedo values to soil colour classification:
                                                          !! dry soil albedo values in visible range
!$OMP THREADPRIVATE(vis_dry)
  REAL(r_std), SAVE, DIMENSION(classnb) :: nir_dry = (/0.48,&
       &0.44, 0.40, 0.36, 0.32, 0.28, 0.24, 0.20, 0.55/)  !! Soil albedo values to soil colour classification:
                                                          !! dry soil albedo values in near-infrared range 
!$OMP THREADPRIVATE(nir_dry)
  REAL(r_std), SAVE, DIMENSION(classnb) :: vis_wet = (/0.12,&
       &0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15/)  !! Soil albedo values to soil colour classification:
                                                          !! wet soil albedo values in visible range 
!$OMP THREADPRIVATE(vis_wet)
  REAL(r_std), SAVE, DIMENSION(classnb) :: nir_wet = (/0.24,&
       &0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.31/)  !! Soil albedo values to soil colour classification:
                                                          !! wet soil albedo values in near-infrared range
!$OMP THREADPRIVATE(nir_wet)
  REAL(r_std), SAVE, DIMENSION(classnb) :: albsoil_vis = (/ &
       &0.18, 0.16, 0.16, 0.15, 0.12, 0.105, 0.09, 0.075, 0.25/)   !! Soil albedo values to soil colour classification:
                                                                   !! Averaged of wet and dry soil albedo values
                                                                   !! in visible and near-infrared range
!$OMP THREADPRIVATE(albsoil_vis) 
  REAL(r_std), SAVE, DIMENSION(classnb) :: albsoil_nir = (/ &
       &0.36, 0.34, 0.34, 0.33, 0.30, 0.25, 0.20, 0.15, 0.45/)  !! Soil albedo values to soil colour classification:
                                                                !! Averaged of wet and dry soil albedo values
                                                                !! in visible and near-infrared range
!$OMP THREADPRIVATE(albsoil_nir)

  LOGICAL, SAVE :: new_frac_snow_veg = .TRUE.     !! Arsene 29-12-2015 - Change the calculation of frac_snow_veg (conveg)  
!$OMP THREADPRIVATE(new_frac_snow_veg)            

  !
  ! diffuco.f90
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: Tetens_1 = 0.622         !! Ratio between molecular weight of water vapor and molecular weight  
                                                     !! of dry air (unitless)
  REAL(r_std), PARAMETER :: Tetens_2 = 0.378         !!
  REAL(r_std), PARAMETER :: ratio_H2O_to_CO2 = 1.6   !! Ratio of water vapor diffusivity to the CO2 diffusivity (unitless)
  REAL(r_std), PARAMETER :: mmol_to_m_1 = 0.0244     !!
  REAL(r_std), PARAMETER :: RG_to_PAR = 0.5          !!
  REAL(r_std), PARAMETER :: W_to_mmol = 4.6          !! W_to_mmol * RG_to_PAR = 2.3

  ! 1. Scalar

  INTEGER(i_std), SAVE :: nlai = 20             !! Number of LAI levels (unitless)
!$OMP THREADPRIVATE(nlai)
  LOGICAL, SAVE :: ldq_cdrag_from_gcm = .FALSE. !! Set to .TRUE. if you want q_cdrag coming from GCM
!$OMP THREADPRIVATE(ldq_cdrag_from_gcm)
  REAL(r_std), SAVE :: laimax = 12.             !! Maximal LAI used for splitting LAI into N layers (m^2.m^{-2})
!$OMP THREADPRIVATE(laimax)
  LOGICAL, SAVE :: downregulation_co2 = .FALSE.            !! Set to .TRUE. if you want CO2 downregulation.
!$OMP THREADPRIVATE(downregulation_co2)
  REAL(r_std), SAVE :: downregulation_co2_baselevel = 280. !! CO2 base level (ppm)
!$OMP THREADPRIVATE(downregulation_co2_baselevel)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: lai_level_depth = 0.15  !!
!$OMP THREADPRIVATE(lai_level_depth)
!
  REAL(r_std), SAVE, DIMENSION(6) :: dew_veg_poly_coeff = &            !! coefficients of the 5 degree polynomomial used
  & (/ 0.887773, 0.205673, 0.110112, 0.014843,  0.000824,  0.000017 /) !! in the equation of coeff_dew_veg
!$OMP THREADPRIVATE(dew_veg_poly_coeff)
!
  REAL(r_std), SAVE               :: Oi=210000.    !! Intercellular oxygen partial pressure (ubar)
!$OMP THREADPRIVATE(Oi)
  !
  ! slowproc.f90 
  !

  ! 1. Scalar

  INTEGER(i_std), SAVE :: veget_year_orig = 0        !!  first year for landuse (number)
!$OMP THREADPRIVATE(veget_year_orig)
  REAL(r_std), SAVE :: clayfraction_default = 0.2    !! Default value for clay fraction (0-1, unitless)
!$OMP THREADPRIVATE(clayfraction_default)
  REAL(r_std), SAVE :: min_vegfrac = 0.001           !! Minimal fraction of mesh a vegetation type can occupy (0-1, unitless)
!$OMP THREADPRIVATE(min_vegfrac)
  REAL(r_std), SAVE :: frac_nobio_fixed_test_1 = 0.0 !! Value for frac_nobio for tests in 0-dim simulations (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_fixed_test_1)
  
  REAL(r_std), SAVE :: stempdiag_bid = 280.          !! only needed for an initial LAI if there is no restart file
!$OMP THREADPRIVATE(stempdiag_bid)


                           !-----------------------------!
                           !  STOMATE AND LPJ PARAMETERS !
                           !-----------------------------!


  !
  ! lpj_constraints.f90
  !
  
  ! 1. Scalar

  REAL(r_std), SAVE  :: too_long = 5.      !! longest sustainable time without 
                                           !! regeneration (vernalization) (years)
!$OMP THREADPRIVATE(too_long)


  !
  ! lpj_establish.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: estab_max_tree = 0.12   !! Maximum tree establishment rate (0-1, unitless)
!$OMP THREADPRIVATE(estab_max_tree)
  REAL(r_std), SAVE :: estab_max_grass = 0.12  !! Maximum grass establishment rate (0-1, unitless)
!$OMP THREADPRIVATE(estab_max_grass)
  
  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: establish_scal_fact = 5.  !!
!$OMP THREADPRIVATE(establish_scal_fact)
  REAL(r_std), SAVE :: max_tree_coverage = 0.98  !! (0-1, unitless)
!$OMP THREADPRIVATE(max_tree_coverage)
  REAL(r_std), SAVE :: ind_0_estab = 0.2         !! = ind_0 * 10.
!$OMP THREADPRIVATE(ind_0_estab)


  !
  ! lpj_fire.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: tau_fire = 30.           !! Time scale for memory of the fire index (days).
!$OMP THREADPRIVATE(tau_fire)
  REAL(r_std), SAVE :: litter_crit = 200.       !! Critical litter quantity for fire
                                                !! below which iginitions extinguish 
                                                !! @tex $(gC m^{-2})$ @endtex
!$OMP THREADPRIVATE(litter_crit)
  REAL(r_std), SAVE :: fire_resist_struct = 0.5 !!
!$OMP THREADPRIVATE(fire_resist_struct)
  ! 2. Arrays

  REAL(r_std), SAVE, DIMENSION(nparts) :: co2frac = &    !! The fraction of the different biomass 
       & (/ .95, .95, 0., 0.3, 0., 0., .95, .95 /)       !! compartments emitted to the atmosphere 
!$OMP THREADPRIVATE(co2frac)                                                         !! when burned (unitless, 0-1)  

  ! 3. Coefficients of equations

  REAL(r_std), SAVE, DIMENSION(3) :: bcfrac_coeff = (/ .3,  1.3,  88.2 /)         !! (unitless)
!$OMP THREADPRIVATE(bcfrac_coeff)
  REAL(r_std), SAVE, DIMENSION(4) :: firefrac_coeff = (/ 0.45, 0.8, 0.6, 0.13 /)  !! (unitless)
!$OMP THREADPRIVATE(firefrac_coeff)

  !
  ! lpj_gap.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: ref_greff = 0.035         !! Asymptotic maximum mortality rate
                                                 !! @tex $(year^{-1})$ @endtex
!$OMP THREADPRIVATE(ref_greff)
  ! 3. Coefficients of equations

!  REAL(r_std), SAVE :: availability_fact = 0.1   !!
!$OMP THREADPRIVATE(availability_fact)

  !               
  ! lpj_light.f90 
  !              

  ! 1. Scalar
  
  LOGICAL, SAVE :: annual_increase = .TRUE. !! for diagnosis of fpc increase, compare today's fpc to last year's maximum (T) or
                                            !! to fpc of last time step (F)? (true/false)
!$OMP THREADPRIVATE(annual_increase)
  REAL(r_std), SAVE :: min_cover = 0.05     !! For trees, minimum fraction of crown area occupied
                                            !! (due to its branches etc.) (0-1, unitless)
                                            !! This means that only a small fraction of its crown area
                                            !! can be invaded by other trees.
!$OMP THREADPRIVATE(min_cover)
  !
  ! lpj_pftinout.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: min_avail = 0.01         !! minimum availability
!$OMP THREADPRIVATE(min_avail)
  REAL(r_std), SAVE :: ind_0 = 0.02             !! initial density of individuals
!$OMP THREADPRIVATE(ind_0)
  ! 3. Coefficients of equations
  
  REAL(r_std), SAVE :: RIP_time_min = 1.25      !! test whether the PFT has been eliminated lately (years)
!$OMP THREADPRIVATE(RIP_time_min)
  REAL(r_std), SAVE :: npp_longterm_init = 10.  !! Initialisation value for npp_longterm (gC.m^{-2}.year^{-1})
!$OMP THREADPRIVATE(npp_longterm_init)
  REAL(r_std), SAVE :: everywhere_init = 0.05   !!
!$OMP THREADPRIVATE(everywhere_init)


  !
  ! stomate_alloc.f90
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: max_possible_lai = 10. !! (m^2.m^{-2})
  REAL(r_std), PARAMETER :: Nlim_Q10 = 10.         !!
!JCADD
   REAL(r_std), SAVE  ::  reserve_time_cut = 20.
   REAL(r_std), SAVE  ::  lai_happy_cut = 0.25
   REAL(r_std), SAVE  ::  tau_leafinit_cut = 10
   REAL(r_std), SAVE  ::  tau_t2m_14 = 14.
!ENDJCADD
  ! 1. Scalar

  LOGICAL, SAVE :: ok_minres = .TRUE.              !! [DISPENSABLE] Do we try to reach a minimum reservoir even if
                                                   !! we are severely stressed? (true/false)
!$OMP THREADPRIVATE(ok_minres)
  REAL(r_std), SAVE :: reserve_time_tree = 30.     !! Maximum number of days during which
                                                   !! carbohydrate reserve may be used for 
                                                   !! trees (days)
!$OMP THREADPRIVATE(reserve_time_tree)
  REAL(r_std), SAVE :: reserve_time_shrub = 25.    !! Maximum number of days during which   !! Arsene 31-07-2014 modifications
                                                   !! carbohydrate reserve may be used for  !! Arsene 31-07-2014 modifications
                                                   !! shrub (days)                          !! Arsene 31-07-2014 modifications
!$OMP THREADPRIVATE(reserve_time_tree)
  REAL(r_std), SAVE :: reserve_time_grass = 20.    !! Maximum number of days during which
                                                   !! carbohydrate reserve may be used for 
                                                   !! grasses (days)
!$OMP THREADPRIVATE(reserve_time_grass)

  REAL(r_std), SAVE :: f_fruit = 0.1               !! Default fruit allocation (0-1, unitless)
!$OMP THREADPRIVATE(f_fruit)
  REAL(r_std), SAVE :: alloc_sap_above_grass = 1.0 !! fraction of sapwood allocation above ground
                                                   !! for grass (0-1, unitless)
!$OMP THREADPRIVATE(alloc_sap_above_grass)
  REAL(r_std), SAVE :: min_LtoLSR = 0.2            !! Prescribed lower bounds for leaf 
                                                   !! allocation (0-1, unitless)
!$OMP THREADPRIVATE(min_LtoLSR)
  REAL(r_std), SAVE :: max_LtoLSR = 0.5            !! Prescribed upper bounds for leaf 
                                                   !! allocation (0-1, unitless)
!$OMP THREADPRIVATE(max_LtoLSR)
  REAL(r_std), SAVE :: z_nitrogen = 0.2            !! Curvature of the root profile (m)
!$OMP THREADPRIVATE(z_nitrogen)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: Nlim_tref = 25.             !! (C)
!$OMP THREADPRIVATE(Nlim_tref)


  !
  ! stomate_data.f90 
  !

  ! 1. Scalar 

  ! 1.1 Parameters for the pipe model

  REAL(r_std), SAVE :: pipe_tune1 = 100.0        !! crown area = pipe_tune1. stem diameter**(1.6) (Reinicke's theory) (unitless)
!$OMP THREADPRIVATE(pipe_tune1)
  REAL(r_std), SAVE :: pipe_tune2 = 40.0         !! height=pipe_tune2 * diameter**pipe_tune3 (unitless)
!$OMP THREADPRIVATE(pipe_tune2)
  REAL(r_std), SAVE :: pipe_tune3 = 0.5          !! height=pipe_tune2 * diameter**pipe_tune3 (unitless)
!$OMP THREADPRIVATE(pipe_tune3)
!!  REAL(r_std), SAVE :: pipe_tune4 = 0.3          !! needed for stem diameter (unitless)  !! Arsene 11-08-2015 - New maxdia: pipe_tune4 not use
!!!$OMP THREADPRIVATE(pipe_tune4)                  !! Arsene 11-08-2015 - New maxdia: pipe_tune4 not use
  REAL(r_std), SAVE :: pipe_density = 2.e5       !! Density
!$OMP THREADPRIVATE(pipe_density)
  REAL(r_std), SAVE :: pipe_k1 = 8.e3            !! one more SAVE
!$OMP THREADPRIVATE(pipe_k1)
  REAL(r_std), SAVE :: pipe_tune_exp_coeff = 1.6 !! pipe tune exponential coeff (unitless)
!$OMP THREADPRIVATE(pipe_tune_exp_coeff)

!! Arsene 03-08-2015 - New parametrisation for shrub... START
  REAL(r_std), SAVE :: pipe_density_shrub = 2.e5       !! Density for shrub
!$OMP THREADPRIVATE(pipe_density)
  REAL(r_std), SAVE :: pipe_k1_shrub = 13.e3           !! one more SAVE for shrub [TESTS: between 12 to 16]
!$OMP THREADPRIVATE(pipe_k1)
!! Arsene 03-08-2015 - If shrubs_like_trees
  REAL(r_std), SAVE :: pipe_tune1_for_shrub = 217.        !! crown area = pipe_tune1. stem diameter**(1.6) (Reinicke's theory) (unitless)
!$OMP THREADPRIVATE(pipe_tune1)
  REAL(r_std), SAVE :: pipe_tune2_for_shrub = 8.          !! height=pipe_tune2 * diameter**pipe_tune3 (unitless)
!$OMP THREADPRIVATE(pipe_tune2)
  REAL(r_std), SAVE :: pipe_tune3_for_shrub = 0.55        !! height=pipe_tune2 * diameter**pipe_tune3 (unitless)
!$OMP THREADPRIVATE(pipe_tune3)
  REAL(r_std), SAVE :: pipe_tune_exp_coeff_for_shrub = 1.6!! pipe tune exponential coeff (unitless)
!$OMP THREADPRIVATE(pipe_tune_exp_coeff)
!! Arsene 03-08-2015 - If .NOT. shrubs_like_trees
  REAL(r_std), SAVE :: pipe_tune_shrub1 = 10**2.42    !! total crown area = pipe_tune_shrub1. tatal basal area**(pipe_tune_shrub_exp_coeff) (Aiba & Kohyama (1996) (unitless)
!$OMP THREADPRIVATE(pipe_tune1)
  REAL(r_std), SAVE :: pipe_tune_shrub2 = 0.75        !! 1/height=1/(pipe_tune_shrub2 * (100*diameter)**pipe_tune_shrub3) + 1/height_presc (unitless)
!$OMP THREADPRIVATE(pipe_tune2)
  REAL(r_std), SAVE :: pipe_tune_shrub3 = 1.15        !! 1/height=1/(pipe_tune_shrub2 * (100*diameter)**pipe_tune_shrub3) + 1/height_presc (unitless)
!$OMP THREADPRIVATE(pipe_tune3)
  REAL(r_std), SAVE :: pipe_tune_shrub_exp_coeff = 0.8!! pipe tune exponential coeff (unitless)
!$OMP THREADPRIVATE(pipe_tune_exp_coeff)
!! Arsene 03-08-2015 - New parametrisation for shrub... END
!! Arsene 03-09-2015 - For iteration (new allomerty for shrubs)
  LOGICAL, SAVE :: shrub_it_ok = .FALSE.              !! Active shrub iteration for allometry (with DGVM) - else Look-Up Table 
!$OMP THREADPRIVATE(shrub_it_ok)
  REAL(r_std), SAVE :: accept_sigma_it = 0.01         !! Define the precision wanted for the iteration result 
!$OMP THREADPRIVATE(accept_sigma_it)
  REAL(r_std), SAVE :: factor_div_it = 5.             !! Define the factor division bewtween iteration
!$OMP THREADPRIVATE(factor_div_it)
!! Arsene 03-09-2015 - For iteration (new allomerty for shrubs)

!! Arsene 16-10-2015 - For shrub allometry array
  INTEGER(i_std), SAVE :: shrub_allom_lig =200       !! Define the number of line for shrub_allom_array
!$OMP THREADPRIVATE(shrub_allom_lig)
  REAL(r_std), SAVE :: shrub_lim_maxdia = 0.96        !! Define the proportion of real maxdia (= maxdia * shrub_lim_maxdia )
!$OMP THREADPRIVATE(shrub_lim_maxdia)
!! Arsene 16-10-2015 - For shrub allometry array

!! Arsene 13-01-2015 - Roughness for shrubs and snow
  REAL(r_std), SAVE :: z0_sensib = 0.3                !! Define the variation of rougness sensibity for shrub at snow depth limit
!$OMP THREADPRIVATE(z0_sensib)
!! Arsene 13-01-2015 - Roughness for shrubs and snow


  ! 1.2 climatic parameters 

  REAL(r_std), SAVE :: precip_crit = 100.        !! minimum precip, in (mm/year)
!$OMP THREADPRIVATE(precip_crit)
  REAL(r_std), SAVE :: gdd_crit_estab = 150.     !! minimum gdd for establishment of saplings
!$OMP THREADPRIVATE(gdd_crit_estab)
  REAL(r_std), SAVE :: fpc_crit = 0.95           !! critical fpc, needed for light competition and establishment (0-1, unitless)
!$OMP THREADPRIVATE(fpc_crit)

  ! 1.3 sapling characteristics

  REAL(r_std), SAVE :: alpha_grass = 0.5         !! alpha coefficient for grasses (unitless)
!$OMP THREADPRIVATE(alpha_grass)
  REAL(r_std), SAVE :: alpha_tree = 1.           !! alpha coefficient for trees (unitless)
!$OMP THREADPRIVATE(alpha_tree)
  REAL(r_std), SAVE :: alpha_shrub = 0.8         !! alpha coefficient for shrubs (unitless)    !! Arsene 31-07-2014 modifications
!$OMP THREADPRIVATE(alpha_shrub)                                                               !! Arsene 31-07-2014 modifications
  REAL(r_std), SAVE :: mass_ratio_heart_sap = 3. !! mass ratio (heartwood+sapwood)/sapwood (unitless)
!$OMP THREADPRIVATE(mass_ratio_heart_sap)
!    REAL(r_std), SAVE :: shrub_ind_frac = 2000.   !! ratio between number of branches (ind) and real individual  (unitless)     !! Arsene 22-05-2015 modifications
!!$OMP THREADPRIVATE(shrub_ind_frac)               !! Arsene 22-05-2015 modifications

  ! 1.4  time scales for phenology and other processes (in days)

  REAL(r_std), SAVE :: tau_hum_month = 20.        !! (days)       
!$OMP THREADPRIVATE(tau_hum_month)
  REAL(r_std), SAVE :: tau_hum_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_hum_week)
  REAL(r_std), SAVE :: tau_t2m_month = 20.        !! (days)      
!$OMP THREADPRIVATE(tau_t2m_month)
  REAL(r_std), SAVE :: tau_t2m_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_t2m_week)
  REAL(r_std), SAVE :: tau_tsoil_month = 20.      !! (days)     
!$OMP THREADPRIVATE(tau_tsoil_month)
  REAL(r_std), SAVE :: tau_soilhum_month = 20.    !! (days)     
!$OMP THREADPRIVATE(tau_soilhum_month)
  REAL(r_std), SAVE :: tau_gpp_week = 7.          !! (days)  
!$OMP THREADPRIVATE(tau_gpp_week)
  REAL(r_std), SAVE :: tau_gdd = 40.              !! (days)  
!$OMP THREADPRIVATE(tau_gdd)
  REAL(r_std), SAVE :: tau_ngd = 50.              !! (days)  
!$OMP THREADPRIVATE(tau_ngd)
  REAL(r_std), SAVE :: coeff_tau_longterm = 3.    !! (unitless)
!$OMP THREADPRIVATE(coeff_tau_longterm)
  REAL(r_std), SAVE :: tau_longterm               !! (days)  
!$OMP THREADPRIVATE(tau_longterm)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: bm_sapl_carbres = 5.             !!
!$OMP THREADPRIVATE(bm_sapl_carbres)
  REAL(r_std), SAVE :: bm_sapl_sapabove = 0.5           !!
!$OMP THREADPRIVATE(bm_sapl_sapabove)
  REAL(r_std), SAVE :: bm_sapl_heartabove = 2.          !!
!$OMP THREADPRIVATE(bm_sapl_heartabove)
  REAL(r_std), SAVE :: bm_sapl_heartbelow = 2.          !!
!$OMP THREADPRIVATE(bm_sapl_heartbelow)
  REAL(r_std), SAVE :: init_sapl_mass_leaf_nat = 0.1    !!
!$OMP THREADPRIVATE(init_sapl_mass_leaf_nat)
  REAL(r_std), SAVE :: init_sapl_mass_leaf_agri = 1.    !!
!$OMP THREADPRIVATE(init_sapl_mass_leaf_agri)
  REAL(r_std), SAVE :: init_sapl_mass_leaf_novasc = 0.01!! Arsene 18-08-2015 - Add vor non vascular plant. If 1 ==> NaN error (via humtress from leaf biomass= -Infinity)
!$OMP THREADPRIVATE(init_sapl_mass_leaf_novasc)         !! Arsene 18-08-2015 - Add vor non vascular plant. If 0.1 to 0.03 ==>  Negative (>10^5) biomass at day 2
  REAL(r_std), SAVE :: init_sapl_mass_carbres = 5.      !!
!$OMP THREADPRIVATE(init_sapl_mass_carbres)
  REAL(r_std), SAVE :: init_sapl_mass_root = 0.1        !!
!$OMP THREADPRIVATE(init_sapl_mass_root)
  REAL(r_std), SAVE :: init_sapl_mass_fruit = 0.3       !!  
!$OMP THREADPRIVATE(init_sapl_mass_fruit)
  REAL(r_std), SAVE :: cn_sapl_init = 0.5               !!
!$OMP THREADPRIVATE(cn_sapl_init)
  REAL(r_std), SAVE :: migrate_tree = 10.*1.E3          !!
!$OMP THREADPRIVATE(migrate_tree)
  REAL(r_std), SAVE :: migrate_shrub = 10.*1.E3         !! Arsene 31-07-2014 modifications
!$OMP THREADPRIVATE(migrate_shrub)                      !! Arsene 31-07-2014 modifications
  REAL(r_std), SAVE :: migrate_grass = 10.*1.E3         !!
!$OMP THREADPRIVATE(migrate_grass)
  REAL(r_std), SAVE :: lai_initmin_tree = 0.3           !!
!$OMP THREADPRIVATE(lai_initmin_tree)
  REAL(r_std), SAVE :: lai_initmin_shrub = 0.2          !! Arsene 31-07-2014 modifications
!$OMP THREADPRIVATE(lai_initmin_shrub
  REAL(r_std), SAVE :: lai_initmin_grass = 0.1          !!
!$OMP THREADPRIVATE(lai_initmin_grass)
  REAL(r_std), SAVE, DIMENSION(2) :: dia_coeff = (/ 4., 0.5 /)            !!
!$OMP THREADPRIVATE(dia_coeff)
!!  REAL(r_std), SAVE, DIMENSION(2) :: maxdia_coeff =(/ 100., 0.01/)       !! Arsene 11-08-2015 - New maxdia ==> maxdia_coeff not use 
!!!$OMP THREADPRIVATE(maxdia_coeff)                                        !! Arsene 11-08-2015 - New maxdia ==> maxdia_coeff not use
!!  REAL(r_std), SAVE, DIMENSION(4) :: bm_sapl_leaf = (/ 4., 4., 0.8, 5./) !! Arsene 11-08-2015 - New calcul of bm_sapl ==> not use
!!!$OMP THREADPRIVATE(bm_sapl_leaf)                                        !! Arsene 11-08-2015 - New calcul of bm_sapl ==> not use



  !
  ! stomate_litter.f90 
  !

  ! 0. Constants

  REAL(r_std), PARAMETER :: Q10 = 10.               !!

  ! 1. Scalar

  REAL(r_std), SAVE :: z_decomp = 0.2               !!  Maximum depth for soil decomposer's activity (m)
!$OMP THREADPRIVATE(z_decomp)

  ! 2. Arrays

  REAL(r_std), SAVE :: frac_soil_struct_aa = 0.55   !! corresponding to frac_soil(istructural,iactive,iabove) 
!$OMP THREADPRIVATE(frac_soil_struct_aa)
  REAL(r_std), SAVE :: frac_soil_struct_ab = 0.45   !! corresponding to frac_soil(istructural,iactive,ibelow)
!$OMP THREADPRIVATE(frac_soil_struct_ab)
  REAL(r_std), SAVE :: frac_soil_struct_sa = 0.7    !! corresponding to frac_soil(istructural,islow,iabove)
!$OMP THREADPRIVATE(frac_soil_struct_sa)
  REAL(r_std), SAVE :: frac_soil_struct_sb = 0.7    !! corresponding to frac_soil(istructural,islow,ibelow)
!$OMP THREADPRIVATE(frac_soil_struct_sb)
  REAL(r_std), SAVE :: frac_soil_metab_aa = 0.45    !! corresponding to frac_soil(imetabolic,iactive,iabove)
!$OMP THREADPRIVATE(frac_soil_metab_aa)
  REAL(r_std), SAVE :: frac_soil_metab_ab = 0.45    !! corresponding to frac_soil(imetabolic,iactive,ibelow)
!$OMP THREADPRIVATE(frac_soil_metab_ab)
  REAL(r_std), SAVE, DIMENSION(nparts) :: CN = &    !! C/N ratio of each plant pool (0-100, unitless)
       & (/ 40., 40., 40., 40., 40., 40., 40., 40. /) 
!$OMP THREADPRIVATE(CN)
  REAL(r_std), SAVE, DIMENSION(nparts) :: LC = &    !! Lignin/C ratio of different plant parts (0,22-0,35, unitless)
       & (/ 0.22, 0.35, 0.35, 0.35, 0.35, 0.22, 0.22, 0.22 /)
!$OMP THREADPRIVATE(LC)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: metabolic_ref_frac = 0.85    !! used by litter and soilcarbon (0-1, unitless)
!$OMP THREADPRIVATE(metabolic_ref_frac)
  REAL(r_std), SAVE :: metabolic_LN_ratio = 0.018   !! (0-1, unitless)   
!$OMP THREADPRIVATE(metabolic_LN_ratio)
  REAL(r_std), SAVE :: tau_metabolic = 0.066        !!
!$OMP THREADPRIVATE(tau_metabolic)
  REAL(r_std), SAVE :: tau_struct = 0.245           !!
!$OMP THREADPRIVATE(tau_struct)
  REAL(r_std), SAVE :: soil_Q10 = 0.69              !!= ln 2
!$OMP THREADPRIVATE(soil_Q10)
  REAL(r_std), SAVE :: tsoil_ref = 30.              !!
!$OMP THREADPRIVATE(tsoil_ref)
  REAL(r_std), SAVE :: litter_struct_coef = 3.      !! 
!$OMP THREADPRIVATE(litter_struct_coef)
  REAL(r_std), SAVE, DIMENSION(3) :: moist_coeff = (/ 1.1,  2.4,  0.29 /) !!
!$OMP THREADPRIVATE(moist_coeff)
  REAL(r_std), SAVE :: moistcont_min = 0.25  !! minimum soil wetness to limit the heterotrophic respiration
!$OMP THREADPRIVATE(moistcont_min)
!! Arsene 29-12-2015 - START - ADD for LUT: new litter moist dependence. IF new_moist_func
  REAL(r_std), SAVE, DIMENSION(4) :: moist_coeff_new = (/ -1.4,  2.22,  -1.12, 1.178 /) !!
!$OMP THREADPRIVATE(moist_coeff_new)
  REAL(r_std), SAVE :: moist_interval = 0.01  !! Interval for Array. Take care: 1./moist_interval have to be an interger !
!$OMP THREADPRIVATE(moist_interval)           !! TAKE CARE: on stomate litter, valid only if = 0.01
!! Arsene 29-12-2015 - END - ADD for LUT: new litter moist dependence

  !
  ! stomate_lpj.f90
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: frac_turnover_daily = 0.55  !! (0-1, unitless)
!$OMP THREADPRIVATE(frac_turnover_daily)

  REAL(r_std), SAVE :: fact_min_height = 10.  !! (1-10, unitless)  !! Arsene 20-05-2015 Add - Note: impoortante pour calc of bm_sap
!$OMP THREADPRIVATE(fact_min_height)                              !! Arsene 20-05-2015 Add

  !
  ! stomate_npp.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: tax_max = 0.8 !! Maximum fraction of allocatable biomass used 
                                     !! for maintenance respiration (0-1, unitless)
!$OMP THREADPRIVATE(tax_max)


  !
  ! stomate_phenology.f90
  !

  ! 1. Scalar

  LOGICAL, SAVE :: always_init = .FALSE.           !! take carbon from atmosphere if carbohydrate reserve too small? (true/false)
!$OMP THREADPRIVATE(always_init)
  REAL(r_std), SAVE :: min_growthinit_time = 300.  !! minimum time since last beginning of a growing season (days)
!$OMP THREADPRIVATE(min_growthinit_time)
  REAL(r_std), SAVE :: moiavail_always_tree = 1.0  !! moisture monthly availability above which moisture tendency doesn't matter
                                                   !!  - for trees (0-1, unitless)
!$OMP THREADPRIVATE(moiavail_always_tree)
  REAL(r_std), SAVE :: moiavail_always_shrub = 0.9 !! moisture monthly availability above which moisture tendency doesn't matter  !! Arsene 31-07-2014 modifications
                                                   !!  - for shrubs (0-1, unitless)                                               !! Arsene 31-07-2014 modifications
!$OMP THREADPRIVATE(moiavail_always_shrub)                                                                                        !! Arsene 31-07-2014 modifications
  REAL(r_std), SAVE :: moiavail_always_grass = 0.6 !! moisture monthly availability above which moisture tendency doesn't matter
                                                   !! - for grass (0-1, unitless)
!$OMP THREADPRIVATE(moiavail_always_grass)
  REAL(r_std), SAVE :: t_always                    !! monthly temp. above which temp. tendency doesn't matter
!$OMP THREADPRIVATE(t_always)
  REAL(r_std), SAVE :: t_always_add = 10.          !! monthly temp. above which temp. tendency doesn't matter (C)
!$OMP THREADPRIVATE(t_always_add)

  ! 3. Coefficients of equations
  
  REAL(r_std), SAVE :: gddncd_ref = 603.           !!
!$OMP THREADPRIVATE(gddncd_ref)
  REAL(r_std), SAVE :: gddncd_curve = 0.0091       !!
!$OMP THREADPRIVATE(gddncd_curve)
  REAL(r_std), SAVE :: gddncd_offset = 64.         !!
!$OMP THREADPRIVATE(gddncd_offset)


  !
  ! stomate_prescribe.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: bm_sapl_rescale = 40.       !!
!$OMP THREADPRIVATE(bm_sapl_rescale)


  !
  ! stomate_resp.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: maint_resp_min_vmax = 0.3   !!
!$OMP THREADPRIVATE(maint_resp_min_vmax)
  REAL(r_std), SAVE :: maint_resp_coeff = 1.4      !!
!$OMP THREADPRIVATE(maint_resp_coeff)


  !
  ! stomate_soilcarbon.f90 
  !

  ! 2. Arrays 

  ! 2.1 frac_carb_coefficients

  REAL(r_std), SAVE :: frac_carb_ap = 0.004  !! from active pool: depends on clay content  (0-1, unitless)
                                             !! corresponding to frac_carb(:,iactive,ipassive)
!$OMP THREADPRIVATE(frac_carb_ap)
  REAL(r_std), SAVE :: frac_carb_sa = 0.42   !! from slow pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,islow,iactive)
!$OMP THREADPRIVATE(frac_carb_sa)
  REAL(r_std), SAVE :: frac_carb_sp = 0.03   !! from slow pool (0-1, unitless) 
                                             !! corresponding to frac_carb(:,islow,ipassive)
!$OMP THREADPRIVATE(frac_carb_sp)
  REAL(r_std), SAVE :: frac_carb_pa = 0.45   !! from passive pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,ipassive,iactive)
!$OMP THREADPRIVATE(frac_carb_pa)
  REAL(r_std), SAVE :: frac_carb_ps = 0.0    !! from passive pool (0-1, unitless)
                                             !! corresponding to frac_carb(:,ipassive,islow)
!$OMP THREADPRIVATE(frac_carb_ps)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: active_to_pass_clay_frac = 0.68  !! (0-1, unitless)
!$OMP THREADPRIVATE(active_to_pass_clay_frac)
  !! residence times in carbon pools (days)
  REAL(r_std), SAVE :: carbon_tau_iactive = 0.149   !! residence times in active pool (days)
!$OMP THREADPRIVATE(carbon_tau_iactive)
  REAL(r_std), SAVE :: carbon_tau_islow = 5.48      !! residence times in slow pool (days)
!$OMP THREADPRIVATE(carbon_tau_islow)
  REAL(r_std), SAVE :: carbon_tau_ipassive = 241.   !! residence times in passive pool (days)
!$OMP THREADPRIVATE(carbon_tau_ipassive)
  REAL(r_std), SAVE, DIMENSION(3) :: flux_tot_coeff = (/ 1.2, 1.4, .75/)
!$OMP THREADPRIVATE(flux_tot_coeff)

  !
  ! stomate_turnover.f90
  !

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: new_turnover_time_ref = 20. !!(days)
!$OMP THREADPRIVATE(new_turnover_time_ref)
  REAL(r_std), SAVE :: dt_turnover_time = 10.      !!(days)
!$OMP THREADPRIVATE(dt_turnover_time)
  REAL(r_std), SAVE :: leaf_age_crit_tref = 20.    !! (C)
!$OMP THREADPRIVATE(leaf_age_crit_tref)
  REAL(r_std), SAVE, DIMENSION(3) :: leaf_age_crit_coeff = (/ 1.5, 0.75, 10./) !! (unitless)
!$OMP THREADPRIVATE(leaf_age_crit_coeff)

  REAL(r_std), SAVE, DIMENSION(4) :: npp0_c = (/20., 60., 130., 0.01/) !!(/days,days,days,unitless) !! Arsene 25-06-2014 NPPcumul
!$OMP THREADPRIVATE(npp0_c)                                                                         !! Arsene 25-06-2014 NPPcumul
  REAL(r_std), SAVE :: llai_coef = 0.06                           !! (unitless)       !! Arsene 29-04-2015 LLAI_COEF for NVP turnover 
!$OMP THREADPRIVATE(llai_coef)                                                                         !! Arsene 29-04-2015 LLAI_COEF for NVP turnover
!! Arsene 25-06-2014 NPPcumul. 1: nb of day before impact of npp0_cumu. 2:max of impact. 3: end of impact. 4: max dturnover(%) by day


  !
  ! stomate_vmax.f90
  !
 
  ! 1. Scalar

  REAL(r_std), SAVE :: vmax_offset = 0.3        !! minimum leaf efficiency (unitless)
!$OMP THREADPRIVATE(vmax_offset)
  REAL(r_std), SAVE :: leafage_firstmax = 0.03  !! relative leaf age at which efficiency
                                                !! reaches 1 (unitless)
!$OMP THREADPRIVATE(leafage_firstmax)
  REAL(r_std), SAVE :: leafage_lastmax = 0.5    !! relative leaf age at which efficiency
                                                !! falls below 1 (unitless)
!$OMP THREADPRIVATE(leafage_lastmax)
  REAL(r_std), SAVE :: leafage_old = 1.         !! relative leaf age at which efficiency
                                                !! reaches its minimum (vmax_offset) 
                                                !! (unitless)
!$OMP THREADPRIVATE(leafage_old)

  REAL(r_std), SAVE :: vcmax_offset = 0.3       !! offset of vcmax reduce by humrel_month (unitless) [0-1]   !! Arsene 30-04-2015 Add
!$OMP THREADPRIVATE(vcmax_offset)               !! offset of vcmax reduce by humrel_month (unitless) [0-1]   !! Arsene 30-04-2015 Add 
  REAL(r_std), SAVE :: humrel_mmin = 0.8        !! Below this value, impact of humrel_month (unitless) [0-1] !! Arsene 30-04-2015 Add
!$OMP THREADPRIVATE(humrel_mmin)                !! Below this value, impact of humrel_month (unitless) [0-1] !! Arsene 30-04-2015 Add 
  REAL(r_std), SAVE :: vcmax_offset2 = 0.2      !! offset of vcmax reduce by humrel_month (unitless) [0-1]   !! Arsene 04-11-2015 Add
!$OMP THREADPRIVATE(vcmax_offset2)              !! offset of vcmax reduce by humrel_month (unitless) [0-1]   !! Arsene 04-11-2015 Add 
  REAL(r_std), SAVE :: humrel_mmax = 0.97       !! Below this value, impact of humrel_month (unitless) [0-1] !! Arsene 04-11-2015 Add
!$OMP THREADPRIVATE(humrel_mmax)                !! Below this value, impact of humrel_month (unitless) [0-1] !! Arsene 04-11-2015 Add 

  !
  ! stomate_season.f90 
  !

  ! 1. Scalar

  REAL(r_std), SAVE :: gppfrac_dormance = 0.2  !! report maximal GPP/GGP_max for dormance (0-1, unitless)
!$OMP THREADPRIVATE(gppfrac_dormance)
  REAL(r_std), SAVE :: tau_climatology = 20.   !! tau for "climatologic variables (years)
!$OMP THREADPRIVATE(tau_climatology)
  REAL(r_std), SAVE :: hvc1 = 0.019            !! parameters for herbivore activity (unitless)
!$OMP THREADPRIVATE(hvc1)
  REAL(r_std), SAVE :: hvc2 = 1.38             !! parameters for herbivore activity (unitless)
!$OMP THREADPRIVATE(hvc2)
  REAL(r_std), SAVE :: leaf_frac_hvc = 0.33    !! leaf fraction (0-1, unitless)
!$OMP THREADPRIVATE(leaf_frac_hvc)
  REAL(r_std), SAVE :: tlong_ref_max = 303.1   !! maximum reference long term temperature (K)
!$OMP THREADPRIVATE(tlong_ref_max)
  REAL(r_std), SAVE :: tlong_ref_min = 253.1   !! minimum reference long term temperature (K)
!$OMP THREADPRIVATE(tlong_ref_min)

  ! 3. Coefficients of equations

  REAL(r_std), SAVE :: ncd_max_year = 3.
!$OMP THREADPRIVATE(ncd_max_year)
  REAL(r_std), SAVE :: gdd_threshold = 5.
!$OMP THREADPRIVATE(gdd_threshold)
  REAL(r_std), SAVE :: green_age_ever = 2.
!$OMP THREADPRIVATE(green_age_ever)
  REAL(r_std), SAVE :: green_age_dec = 0.5
!$OMP THREADPRIVATE(green_age_dec)

  ! permafrost carbon related
  REAL(r_std), parameter :: O2_init_conc = 298.3                  !! gO2/m**3 mean for Cherskii
!$OMP THREADPRIVATE(O2_init_conc)
  REAL(r_std), parameter :: CH4_init_conc = 0.001267              !! gCH4/m**3 mean for Cherskii
!$OMP THREADPRIVATE(CH4_init_conc)
  REAL, PARAMETER          :: z_root_max = 0.5                    !! Depth at which litter carbon input decays e-fold (root depth); 0.5 for compar w/WH
!$OMP THREADPRIVATE(z_root_max)
  REAL, PARAMETER          :: diffO2_air = 1.596E-5               !! oxygen diffusivity in air (m**2/s)
!$OMP THREADPRIVATE(diffO2_air)
  REAL, PARAMETER          :: diffO2_w = 1.596E-9                 !! oxygen diffusivity in water (m**2/s)
!$OMP THREADPRIVATE(diffO2_w)
  REAL, PARAMETER          :: O2_surf = 0.209                     !! oxygen concentration in surface air (molar fraction)
!$OMP THREADPRIVATE(O2_surf)
  REAL, PARAMETER          :: diffCH4_air = 1.702E-5              !! methane diffusivity in air (m**2/s)
!$OMP THREADPRIVATE(diffCH4_air)
  REAL, PARAMETER          :: diffCH4_w = 2.0E-9                  !! methane diffusivity in water (m**2/s)
!$OMP THREADPRIVATE(diffCH4_w)
  REAL, PARAMETER          :: CH4_surf = 1700.E-9                 !! methane concentration in surface air (molar fraction)
!$OMP THREADPRIVATE(CH4_surf)
  REAL, SAVE               :: tetasat =  .5                       !! volumetric water content at saturation (porosity)
!$OMP THREADPRIVATE(tetasat)
  REAL, SAVE               :: tetamoss =  0.98                    !! porosity of moss !! Arsene 14-01-2016 - New value (don't use before), from O'donnell et al, 2009
!$OMP THREADPRIVATE(tetamoss)
  REAL, SAVE               :: rho_moss = 0.5E4                    !! density of moss (gC/m**3) !! Arsene 14-01-2016 - Add new value for depth (use in thermosoil)  
!$OMP THREADPRIVATE(rho_moss)
  REAL, PARAMETER          :: zmoss = 0.2                         !! thickness of moss layer (in permafrost regions,m) 0. ! 0.001 DKtest for compar w/WH
!$OMP THREADPRIVATE(zmoss)
  REAL, PARAMETER          :: h_snowmoss = 0.2                    !! snow height above which we consider the moss layer to be to compressed to be effective (m)
!$OMP THREADPRIVATE(h_snowmoss)
  REAL, PARAMETER          :: BunsenO2 = 0.038                    !! Bunsen coefficient for O2 (10C, 1bar)
!$OMP THREADPRIVATE(BunsenO2)
  REAL, PARAMETER          :: BunsenCH4 = 0.043                   !! Bunsen coefficient for CH4 (10C, Wiesenburg et Guinasso, Jr., 1979)
!$OMP THREADPRIVATE(BunsenCH4)
  REAL, PARAMETER          :: ebuthr = 0.9                        !! Soil humidity threshold for ebullition
!$OMP THREADPRIVATE(ebuthr)
  REAL, PARAMETER          :: wCH4 = 16.                          !! molar weight of CH4 (g/mol)
!$OMP THREADPRIVATE(wCH4)
  REAL, PARAMETER          :: wO2 = 32.                           !! molar weight of O2 (g/mol)
!$OMP THREADPRIVATE(wO2)
  REAL, PARAMETER          :: wC = 12.                            !! molar weight of C (g/mol)
!$OMP THREADPRIVATE(wC)
  REAL, PARAMETER          :: avm = .01                           !! minimum air volume (m**3 air/m**3 soil)
!$OMP THREADPRIVATE(avm)
  REAL(R_STD), PARAMETER   :: hmin_tcalc = .001                   !! minimum total snow layer thickness below which we ignore diffusion across the snow layer
!$OMP THREADPRIVATE(hmin_tcalc)

  ! are we running the soil carbon spinup routine?
  LOGICAL, SAVE           :: soilc_isspinup = .false.

  ! which variables to write to history tapes?
  LOGICAL, SAVE           :: writehist_deepC = .true.
  LOGICAL, SAVE           :: writehist_soilgases = .true.
  LOGICAL, SAVE           :: writehist_deltaC = .false.
  LOGICAL, SAVE           :: writehist_zimovheat = .false.
  LOGICAL, SAVE           :: writehist_deltaC_litter = .false.
  LOGICAL, SAVE           :: writehist_gascoeff = .false.

! Arsene 17-03-2016 - Remove and place in constante_mtc to be PFT adaptable
!  LOGICAL,PARAMETER,DIMENSION(16) :: permafrost_veg_exists = &  !!! set all as true to avoid masking     !! Arsene 05-03-2015 remove 13 to nvm
!   & (/ .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.,  .TRUE., &
!   &    .TRUE., .TRUE.,  .TRUE.,  .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE. /)                       !! Arsene 05-03-2015 Add 2 PFT


!pss:+
!valeurs des bornes des differentes classes de WTD pour TOPMODEL
  REAL(r_std),SAVE :: WTD1_borne=0.03
  REAL(r_std),SAVE :: WTD2_borne=0.09
  REAL(r_std),SAVE :: WTD3_borne=0.15
  REAL(r_std),SAVE :: WTD4_borne=0.21
!valeurs du shift de la distribution topo pour passer de fsat a fwet
  REAL(r_std),SAVE :: SHIFT_fsat_fwet=5.
!pss:-

!
! stomate cste WETLAND 
!
!pss+
!ancien para.h sans ngrid et ntime
!  INTERGER(i_std),SAVE  :: nvert = 371
  INTEGER(i_std),SAVE  :: nvert = 171
  INTEGER(i_std),SAVE  :: ns = 151
  INTEGER(i_std),SAVE  :: nday = 24
  REAL(r_std),SAVE  :: h = 0.1
  REAL(r_std),SAVE  :: rk = 1
  REAL(r_std),SAVE  :: rkh = 100 !rk/h**2
  REAL(r_std),SAVE  :: diffair = 7.2
  REAL(r_std),SAVE  :: pox = 0.5
  REAL(r_std),SAVE  :: dveg = 0.001
  REAL(r_std),SAVE  :: rkm = 5.0
  REAL(r_std),SAVE  :: xvmax = 20.0
  REAL(r_std),SAVE  :: oxq10 = 2.0
!  REAL(r_std),SAVE  :: catm = 0.0033
  REAL(r_std),SAVE  :: funit = 3.84 !3.84/rk
  REAL(r_std),SAVE  :: scmax = 500.
  REAL(r_std),SAVE  :: sr0pl = 600.
!!fin de l ancien para.h

!valeur de WTD pour les routines de calcul de densite de flux de CH4
  REAL(r_std),SAVE  :: pwater_wet1=-3
  REAL(r_std),SAVE  :: pwater_wet2=-9
  REAL(r_std),SAVE  :: pwater_wet3=-15
  REAL(r_std),SAVE  :: pwater_wet4=-21
!

  REAL(r_std),SAVE  :: rpv = 0.5
  REAL(r_std),SAVE  :: iother = -1.0
  !!!!plus necessaire maintenant que je mets 2 subroutines differents pour des wt differentes
  !REAL(r_std),SAVE  :: pwater = 0.0

!!rq10 et alpha
!!pour l instant je les mets constantes pour toutes les latitudes
  REAL(r_std),SAVE  :: rq10 = 3.0

!  REAL(r_std),SAVE  :: alpha = 0.010
!  REAL(r_std),SAVE ,DIMENSION(3) :: alpha = (/0.009,0.004,0.021/)
!  REAL(r_std),SAVE , DIMENSION(3) :: alpha_CH4 = (/0.004,0.003,0.018/)
!pss modify parameters
  REAL(r_std),SAVE , DIMENSION(3) :: alpha_CH4 = (/0.006,0.004,0.028/)
!pss-


END MODULE constantes_var
