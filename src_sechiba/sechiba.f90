!  ==============================================================================================================================\n
!  MODULE 	: sechiba
! 
!  CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Structures the calculation of atmospheric and hydrological 
!! variables by calling diffuco_main, enerbil_main, hydrolc_main (or hydrol_main),
!! enerbil_fusion, condveg_main and thermosoil_main. Note that sechiba_main
!! calls slowproc_main and thus indirectly calculates the biogeochemical
!! processes as well.
!!
!!\n DESCRIPTION  : :: shumdiag, :: litterhumdiag and :: stempdiag are not 
!! saved in the restart file because at the first time step because they 
!! are recalculated. However, they must be saved as they are in slowproc 
!! which is called before the modules which calculate them.
!! 
!! RECENT CHANGE(S): None 
!! 
!! REFERENCE(S) : None
!!   
!! SVN     :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_sechiba/sechiba.f90 $ 
!! $Date: 2014-07-17 17:02:21 +0200 (Thu, 17 Jul 2014) $
!! $Revision: 2247 $
!! \n
!_ ================================================================================================================================
 
MODULE sechiba
 
  USE ioipsl
  USE xios_orchidee
  
  ! modules used :
  USE constantes
  USE constantes_soil
  USE pft_parameters
  USE diffuco
  USE condveg
  USE enerbil
  USE hydrol                                                                     !! 
  USE hydrolc
  USE thermosoil
  USE sechiba_io
  USE slowproc
  USE routing
  use ioipsl_para


  IMPLICIT NONE

  ! Private and public routines

  PRIVATE
  PUBLIC sechiba_main,sechiba_clear

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexveg       !! indexing array for the 3D fields of vegetation
!$OMP THREADPRIVATE(indexveg)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexlai       !! indexing array for the 3D fields of vegetation
!$OMP THREADPRIVATE(indexlai)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexnobio     !! indexing array for the 3D fields of other surfaces (ice,
                                                                     !! lakes, ...)
!$OMP THREADPRIVATE(indexnobio)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexsoil      !! indexing array for the 3D fields of soil types (kjpindex*nstm)
!$OMP THREADPRIVATE(indexsoil)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexgrnd      !! indexing array for the 3D ground heat profiles (kjpindex*ngrnd)
!$OMP THREADPRIVATE(indexgrnd)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexlayer     !! indexing array for the 3D fields of soil layers in CWRR (kjpindex*nslm)
!$OMP THREADPRIVATE(indexlayer)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexnbdl      !! indexing array for the 3D fields of diagnostic soil layers (kjpindex*nbdl)
!$OMP THREADPRIVATE(indexnbdl)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexalb       !! indexing array for the 2 fields of albedo
!$OMP THREADPRIVATE(indexalb)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexsnow      !! indexing array for the 3D fields snow layers
!$OMP THREADPRIVATE(indexsnow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: veget          !! Fraction of vegetation type (unitless, 0-1)       
!$OMP THREADPRIVATE(veget)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: veget_max      !! Max. fraction of vegetation type (LAI -> infty, unitless)
!$OMP THREADPRIVATE(veget_max)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: height         !! Vegetation Height (m)
!$OMP THREADPRIVATE(height)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: totfrac_nobio  !! Total fraction of continental ice+lakes+cities+...
                                                                     !! (unitless, 0-1)
!$OMP THREADPRIVATE(totfrac_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: floodout       !! Flow out of floodplains from hydrol
!$OMP THREADPRIVATE(floodout)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: runoff         !! Surface runoff calculated by hydrol or hydrolc 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(runoff)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: drainage       !! Deep drainage calculatedd by hydrol or hydrolc 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(drainage)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: returnflow     !! Water flow from lakes and swamps which returns to 
                                                                     !! the grid box @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(returnflow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: reinfiltration !! Routed water which returns into the soil
!$OMP THREADPRIVATE(reinfiltration)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: irrigation     !! Irrigation flux taken from the routing reservoirs and 
                                                                     !! being put into the upper layers of the soil 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(irrigation)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: emis           !! Surface emissivity (unitless)
!$OMP THREADPRIVATE(emis)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: z0             !! Surface roughness (m)
!$OMP THREADPRIVATE(z0)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: roughheight    !! Effective height for roughness (m)
!$OMP THREADPRIVATE(roughheight)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: reinf_slope    !! slope coefficient (reinfiltration) 
!$OMP THREADPRIVATE(reinf_slope)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag       !! Mean relative soil moisture in the different levels used 
                                                                     !! by thermosoil.f90 (unitless, 0-1)
!$OMP THREADPRIVATE(shumdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag_perma !! Saturation degree of the soil 
!$OMP THREADPRIVATE(shumdiag_perma)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: k_litt         !! litter cond.
!$OMP THREADPRIVATE(k_litt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: litterhumdiag  !! Litter dryness factor (unitless, 0-1)
!$OMP THREADPRIVATE(litterhumdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: stempdiag      !! Temperature which controls canopy evolution (K)
!$OMP THREADPRIVATE(stempdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: qsintveg       !! Water on vegetation due to interception 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(qsintveg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vbeta2         !! Interception resistance (unitless,0-1)
!$OMP THREADPRIVATE(vbeta2)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vbeta3         !! Vegetation resistance (unitless,0-1)
!$OMP THREADPRIVATE(vbeta3)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vbeta3pot      !! Potential vegetation resistance
!$OMP THREADPRIVATE(vbeta3pot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: gsmean         !! Mean stomatal conductance for CO2 (umol m-2 s-1) 
!$OMP THREADPRIVATE(gsmean) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: cimean         !! STOMATE: mean intercellular CO2 concentration (ppm)
!$OMP THREADPRIVATE(cimean)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vevapwet       !! Interception loss over each PFT 
                                                                     !! @tex $(kg m^{-2} days^{-1})$ @endtex
!$OMP THREADPRIVATE(vevapwet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: transpir       !! Transpiration @tex $(kg m^{-2} days^{-1})$ @endtex
!$OMP THREADPRIVATE(transpir)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: transpot       !! Potential Transpiration (needed for irrigation)
!$OMP THREADPRIVATE(transpot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: qsintmax       !! Maximum amount of water in the canopy interception 
                                                                     !! reservoir @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(qsintmax)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: rveget         !! Surface resistance for the vegetation 
                                                                     !! @tex $(s m^{-1})$ @endtex
!$OMP THREADPRIVATE(rveget)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: rstruct        !! Vegetation structural resistance
!$OMP THREADPRIVATE(rstruct)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: snow_nobio     !! Snow mass of non-vegetative surfaces 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(snow_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: snow_nobio_age !! Snow age on non-vegetative surfaces (days)
!$OMP THREADPRIVATE(snow_nobio_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: frac_nobio     !! Fraction of non-vegetative surfaces (continental ice, 
                                                                     !! lakes, ...) (unitless, 0-1)
!$OMP THREADPRIVATE(frac_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: albedo         !! Surface albedo for visible and near-infrared 
                                                                     !! (unitless, 0-1)
!$OMP THREADPRIVATE(albedo)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: assim_param    !! min+max+opt temps, vcmax, vjmax for photosynthesis
!$OMP THREADPRIVATE(assim_param)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: lai            !! Surface foliaire
!$OMP THREADPRIVATE(lai)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: gpp            !! STOMATE: GPP. gC/m**2 of total area
!$OMP THREADPRIVATE(gpp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: temp_growth       !! Growth temperature (Â°C) - Is equal to t2m_month 
!$OMP THREADPRIVATE(temp_growth) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: humrel         !! Relative humidity
!$OMP THREADPRIVATE(humrel)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vegstress      !! Vegetation moisture stress (only for vegetation growth)
!$OMP THREADPRIVATE(vegstress)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: frac_age       !! Age efficacity from STOMATE for isoprene 
!$OMP THREADPRIVATE(frac_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: soiltile       !! Fraction of each soil tile (0-1, unitless)
!$OMP THREADPRIVATE(soiltile)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
!$OMP THREADPRIVATE(njsc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vbeta1         !! Snow resistance 
!$OMP THREADPRIVATE(vbeta1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vbeta4         !! Bare soil resistance
!$OMP THREADPRIVATE(vbeta4)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vbeta5         !! Floodplains resistance
!$OMP THREADPRIVATE(vbeta5)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: soilcap        !!
!$OMP THREADPRIVATE(soilcap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: soilflx        !!
!$OMP THREADPRIVATE(soilflx)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: temp_sol       !! Soil temperature
!$OMP THREADPRIVATE(temp_sol)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: qsurf          !! near soil air moisture
!$OMP THREADPRIVATE(qsurf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: flood_res      !! flood reservoir estimate
!$OMP THREADPRIVATE(flood_res)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: flood_frac     !! flooded fraction
!$OMP THREADPRIVATE(flood_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: snow           !! Snow mass [Kg/m^2]
!$OMP THREADPRIVATE(snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: snow_age       !! Snow age
!$OMP THREADPRIVATE(snow_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: drysoil_frac   !! Fraction of visibly (albedo) Dry soil (Between 0 and 1)
!$OMP THREADPRIVATE(drysoil_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: rsol           !! resistance to bare soil evaporation
!$OMP THREADPRIVATE(rsol)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evap_bare_lim  !! Bare soil stress
!$OMP THREADPRIVATE(evap_bare_lim)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: co2_flux       !! CO2 flux (gC/m**2 of average ground/s)
!$OMP THREADPRIVATE(co2_flux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evapot         !! Soil Potential Evaporation
!$OMP THREADPRIVATE(evapot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evapot_corr    !! Soil Potential Evaporation Correction (Milly 1992)
!$OMP THREADPRIVATE(evapot_corr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapflo       !! Floodplains evaporation
!$OMP THREADPRIVATE(vevapflo)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapsno       !! Snow evaporation
!$OMP THREADPRIVATE(vevapsno)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapnu        !! Bare soil evaporation
!$OMP THREADPRIVATE(vevapnu)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: t2mdiag        !! 2 meter temperature
!$OMP THREADPRIVATE(t2mdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: tot_melt       !! Total melt
!$OMP THREADPRIVATE(tot_melt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vbeta          !! Resistance coefficient
!$OMP THREADPRIVATE(vbeta)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: valpha         !! Resistance coefficient
!$OMP THREADPRIVATE(valpha)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: fusion         !!
!$OMP THREADPRIVATE(fusion)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: rau            !! Density
!$OMP THREADPRIVATE(rau)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: deadleaf_cover !! Fraction of soil covered by dead leaves
!$OMP THREADPRIVATE(deadleaf_cover)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: ptnlev1        !! 1st level Different levels soil temperature
!$OMP THREADPRIVATE(ptnlev1)


  LOGICAL, SAVE                                    :: l_first_sechiba = .TRUE. !! Flag controlling the intialisation (true/false)
!$OMP THREADPRIVATE(l_first_sechiba)
  LOGICAL, SAVE                                    :: river_routing            !! Flag that decides if we route surface runoff 
                                                                               !! and deep drainage to the ocean through the 
                                                                               !! rivers
!$OMP THREADPRIVATE(river_routing)
  LOGICAL, SAVE                                    :: hydrol_cwrr              !! Flag that decides if we use the 11-layer 
                                                                               !! hydrology module. If not, Choisnel module is 
                                                                               !! taken.
!$OMP THREADPRIVATE(hydrol_cwrr)
  LOGICAL, SAVE                                    :: myfalse=.FALSE.          !! Local flag (true/false)
!$OMP THREADPRIVATE(myfalse)

  ! Variables related to snow processes calculations  

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowrho      !! snow density for each layer
!$OMP THREADPRIVATE(snowrho)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowheat     !! snow heat content for each layer (J/m2)
!$OMP THREADPRIVATE(snowheat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowliq      !! liquid water content (m)
!$OMP THREADPRIVATE(snowliq)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowgrain    !! snow grain size (m)
!$OMP THREADPRIVATE(snowgrain)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowtemp     !! snow temperature profile (K)
!$OMP THREADPRIVATE(snowtemp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: soiltemp     !! soil temperature profile in thermal layers (K)
!$OMP THREADPRIVATE(soiltemp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowdz       !! snow layer thickness (m)
!$OMP THREADPRIVATE(snowdz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: grndflux     !! net energy into soil (W/m2)
!$OMP THREADPRIVATE(grndflux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: gthick       !! soil surface layer thickness
!$OMP THREADPRIVATE(gthick)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: gtemp        !! soil surface temperature
!$OMP THREADPRIVATE(gtemp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: soilflxresid !! Net flux to the snowpack
!$OMP THREADPRIVATE(soilflxresid)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: gpkappa      !! soil surface conductivity
!$OMP THREADPRIVATE(gpkappa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: pgflux       !! net energy into snow pack
!$OMP THREADPRIVATE(pgflux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa_snow      !! snow thermal conductivity
!$OMP THREADPRIVATE(pkappa_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: cgrnd_soil
!$OMP THREADPRIVATE(cgrnd_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: dgrnd_soil
!$OMP THREADPRIVATE(dgrnd_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: zdz1_soil
!$OMP THREADPRIVATE(zdz1_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: zdz2_soil
!$OMP THREADPRIVATE(zdz2_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: cgrnd_snow
!$OMP THREADPRIVATE(cgrnd_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: dgrnd_snow
!$OMP THREADPRIVATE(dgrnd_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: lambda_snow
!$OMP THREADPRIVATE(lambda_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: snowflx
!$OMP THREADPRIVATE(snowflx)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: snowcap
!$OMP THREADPRIVATE(snowcap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: temp_sol_add
!$OMP THREADPRIVATE(temp_sol_add)
  ! Variables related to deep permafrost calculations
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: tdeep          !! deep temperature profile
!$OMP THREADPRIVATE(tdeep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: hsdeep         !! deep soil humidity profile
!$OMP THREADPRIVATE(hsdeep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: heat_Zimov     !! heating associated with decomposition
!$OMP THREADPRIVATE(heat_Zimov)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: sfluxCH4_deep    !! surface flux of CH4 to atmosphere from permafrost
!$OMP THREADPRIVATE(sfluxCH4_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: sfluxCO2_deep    !! surface flux of CO2 to atmosphere from permafrost
!$OMP THREADPRIVATE(sfluxCO2_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: thawed_humidity  
!$OMP THREADPRIVATE(thawed_humidity)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: depth_organic_soil
!$OMP THREADPRIVATE(depth_organic_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: zz_deep
!$OMP THREADPRIVATE(zz_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: zz_coef_deep
!$OMP THREADPRIVATE(zz_coef_deep)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: soilc_total    !! total  soil carbon for use in thermal calcs
!$OMP THREADPRIVATE(soilc_total)

!pss:+
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: drunoff_tot         !! Surface runoff generated Dune process
!$OMP THREADPRIVATE(drunoff_tot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: fwet_out      !! wetland fraction
!$OMP THREADPRIVATE(fwet_out)

!pss:-

CONTAINS


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_main
!!
!>\BRIEF        Main routine for the sechiba module performing three functions:
!! variable initialization (during the first call only), calculating temporal
!! evolution of all variables (every call including the first call) and preparation
!! of output and restart files (during the last call only)
!!
!!\n DESCRIPTION : Main routine for the sechiba module. This module is called 
!! two times:
!! - a first time to set initial values
!! - a second time to compute the entire algorithm.\n 
!! Every time the module is called, three major if/then loops are checked:
!! - initialization (first call only),
!! - time step evolution of all variables (every call including the first),
!! - preparation of output and storage for the restart arrays (last call only).\n
!! One time step evolution consists of:
!! - call sechiba_var_init to do some initialization,
!! - call slowproc_main to do some daily initialization,
!! - call diffuco_main for diffusion coefficient calculation,
!! - call enerbil_main for energy budget calculation,
!! - call hydrolc_main (or hydrol_main) for hydrologic processes calculation,
!! - call enerbil_fusion : last part with fusion,
!! - call condveg_main for surface conditions such as roughness, albedo, and emmisivity,
!! - call thermosoil_main for soil thermodynamic calculation,
!! - call sechiba_end to swap previous to new fields.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): Hydrological variables (:: coastalflow and :: riverflow),
!! components of the energy budget (:: tsol_rad, :: vevapp, :: fluxsens, 
!! :: temp_sol_new and :: fluxlat), surface characteristics (:: z0_out, :: emis_out, 
!! :: tq_cdrag and :: albedo_out) and land use related CO2 fluxes (:: netco2flux and 
!! :: fco2_lu)            
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale = 0.5]{sechibamainflow.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE sechiba_main (kjit, kjpij, kjpindex, index, dtradia, date0, &
       & ldrestart_read, ldrestart_write, control_in, &
       & lalo, contfrac, neighbours, resolution,&
       ! First level conditions
       & zlev, u, v, qair, q2m, t2m, temp_air, epot_air, ccanopy, &
       ! Variables for the implicit coupling
       & tq_cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       ! Rain, snow, radiation and surface pressure
       & precip_rain, precip_snow, lwdown, swnet, swdown, sinang, pb, &
       ! Output : Fluxes
       & vevapp, fluxsens, fluxlat, coastalflow, riverflow, netco2flux, fco2_lu, &
       ! Surface temperatures and surface properties
       & tsol_rad, temp_sol_new, qsurf_out, albedo_out, emis_out, z0_out, &
       ! File ids
       & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC)

!! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                               :: kjit              !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpij             !! Total size of the un-compressed grid 
                                                                                  !! (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size - terrestrial pixels only 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id           !! _Restart_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id           !! _History_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist2_id          !! _History_ file 2 identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id_stom      !! STOMATE's _Restart_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id_stom      !! STOMATE's _History_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT(in)                                :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file file 
                                                                                  !! identifier (unitless)
    REAL(r_std), INTENT (in)                                 :: dtradia           !! Time step (s)
    REAL(r_std), INTENT (in)                                 :: date0             !! Initial date (??unit??)
    LOGICAL, INTENT(in)                                      :: ldrestart_read    !! Logical for _restart_ file to read 
                                                                                  !! (true/false)
    LOGICAL, INTENT(in)                                      :: ldrestart_write   !! Logical for _restart_ file to write 
                                                                                  !! (true/false)
    TYPE(control_type), INTENT(in)                           :: control_in        !! Flags that (de)activate parts of the model
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)          :: lalo              !! Geographic coordinates (latitude,longitude)
                                                                                  !! for grid cells (degrees)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: contfrac          !! Fraction of continent in the grid 
                                                                                  !! (unitless, 0-1)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: index             !! Indices of the pixels on the map. 
                                                                                  !! Sechiba uses a reduced grid excluding oceans
                                                                                  !! ::index contains the indices of the 
                                                                                  !! terrestrial pixels only! (unitless)
    INTEGER(i_std), DIMENSION (kjpindex,8), INTENT(in)       :: neighbours        !! Neighboring grid points if land!(unitless)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)          :: resolution        !! Size in x and y of the grid (m)
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u                 !! Lowest level wind speed in direction u 
                                                                                  !! @tex $(m.s^{-1})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v                 !! Lowest level wind speed in direction v 
                                                                                  !! @tex $(m.s^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: zlev              !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair              !! Lowest level specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: q2m               !! 2m specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: t2m               !! 2m air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_rain       !! Rain precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_snow       !! Snow precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown            !! Down-welling long-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: sinang            !! Sine of the solar angle (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet             !! Net surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown            !! Down-welling surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_air          !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: epot_air          !! Air potential energy (??J)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: ccanopy           !! CO2 concentration in the canopy (ppm)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petAcoef          !! Coefficients A for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqAcoef          !! Coefficients A for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petBcoef          !! Coefficients B for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqBcoef          !! Coefficients B for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb                !! Surface pressure (hPa)


!! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: coastalflow       !! Outflow on coastal points by small basins.
                                                                                  !! This is the water which flows in a disperse 
                                                                                  !! way into the ocean
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: riverflow         !! Outflow of the major rivers.
                                                                                  !! The flux will be located on the continental 
                                                                                  !! grid but this should be a coastal point  
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: tsol_rad          !! Radiative surface temperature 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vevapp            !! Total of evaporation 
                                                                                  !! @tex $(kg m^{-2} days^{-1})$ @endtex
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsurf_out         !! Surface specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0_out            !! Surface roughness (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out)         :: albedo_out        !! VIS and NIR albedo (output diagnostic, 
                                                                                  !! unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsens          !! Sensible heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxlat           !! Latent heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: emis_out          !! Emissivity (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: netco2flux        !! Sum CO2 flux over PFTs 
                                                                                  !! ??(gC m^{-2} s^{-1})??
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fco2_lu           !! Land Cover Change CO2 flux 
                                                                                  !! ??(gC m^{-2} s^{-1})??

!! 0.3 Modified

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: tq_cdrag          !! Surface drag coefficient 
                                                                                  !! @tex $(m.s^{-1})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: temp_sol_new      !! New ground temperature (K)

!! 0.4 local variables

    INTEGER(i_std)                                           :: ji, jv		  !! Index (unitless)
    REAL(r_std), ALLOCATABLE, DIMENSION (:)                  :: runoff1           !! ??Temporary surface runoff calculated by 
                                                                                  !! hydrol @tex $(kg m^{-2})$ @endtex
    REAL(r_std), ALLOCATABLE, DIMENSION (:)                  :: drainage1         !! ??Temporary deep drainage calculatedd by 
                                                                                  !! hydrol @tex $(kg m^{-2})$ @endtex
    REAL(r_std), ALLOCATABLE, DIMENSION (:)                  :: soilcap1          !! ??Temporary soil heat capacity 
                                                                                  !! @tex $(J K^{-1})$ @endtex
    REAL(r_std), ALLOCATABLE, DIMENSION (:)                  :: soilflx1          !! ??Temporary soil heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std), ALLOCATABLE, DIMENSION (:)                  :: snowcap1          !! ??Temporary soil heat capacity
    REAL(r_std), ALLOCATABLE, DIMENSION (:)                  :: snowflx1

    REAL(r_std), ALLOCATABLE, DIMENSION (:,:)                :: shumdiag1         !! ??Temporary relative soil moisture 
                                                                                  !! (unitless, 0-1)
    REAL(r_std), DIMENSION(kjpindex)                         :: histvar           !! Computations for history files (unitless)
    CHARACTER(LEN=80)                                        :: var_name          !! To store variables names for I/O (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treefrac      !! Total fraction occupied by trees (0-1, uniless) 
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_grassfrac     !! Total fraction occupied by grasses (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_cropfrac      !! Total fraction occcupied by crops (0-1, unitess)


!_ ================================================================================================================================

    IF (long_print) WRITE(numout,*) ' sechiba kjpindex =',kjpindex

!! 1. Initialize variables on first call!

    IF (l_first_sechiba) THEN

       !! 1.1 Initialize most of sechiba's variables
       CALL sechiba_init (kjit, ldrestart_read, kjpij, kjpindex, index, rest_id, control_in, lalo)
     
       ALLOCATE(runoff1 (kjpindex),drainage1 (kjpindex), soilcap1 (kjpindex),soilflx1 (kjpindex),&
              & snowcap1 (kjpindex),snowflx1 (kjpindex))
       ALLOCATE(shumdiag1(kjpindex,nbdl))
        
       !! 1.2 Initialize some variables of energy budget from restart file     
       IF (ldrestart_read) THEN
          
          IF (long_print) WRITE (numout,*) ' we have to read a restart file for SECHIBA variables'
          var_name='soilcap' ;
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Soil calorific capacity')
          soilcap1=val_exp
          IF ( ok_var(var_name) ) THEN
             CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., soilcap1, "gather", nbp_glo, index_g)
             IF (MINVAL(soilcap1) < MAXVAL(soilcap1) .OR. MAXVAL(soilcap1) < val_exp) THEN
                soilcap(:) = soilcap1(:)
             ENDIF
          ENDIF
          
          var_name='soilflx' ;
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Soil flux')
          soilflx1=val_exp
          IF ( ok_var(var_name) ) THEN
             CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., soilflx1, "gather", nbp_glo, index_g)
             IF (MINVAL(soilflx1) < MAXVAL(soilflx1)  .OR. MAXVAL(soilflx1) < val_exp) THEN
                soilflx(:) = soilflx1(:)
             ENDIF
          ENDIF
         
          var_name='snowcap' ;
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Snow surface calorific capacity')
          snowcap1=val_exp
          IF ( ok_var(var_name) ) THEN
             CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE.,snowcap1, "gather", nbp_glo, index_g)
             IF (MINVAL(snowcap1) < MAXVAL(snowcap1) .OR. MAXVAL(snowcap1) < val_exp) THEN
                snowcap(:) = snowcap1(:)
             ENDIF
          ENDIF

          var_name='snowflx' ;
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Snow surface flux')
          snowflx1=val_exp
          IF ( ok_var(var_name) ) THEN
             CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE.,snowflx1, "gather", nbp_glo, index_g)
             IF (MINVAL(snowflx1) < MAXVAL(snowflx1)  .OR. MAXVAL(snowflx1) < val_exp) THEN
                snowflx(:) = snowflx1(:)
             ENDIF
          ENDIF
 
          var_name='shumdiag' ;
          CALL ioconf_setatt_p('UNITS', '-')
          CALL ioconf_setatt_p('LONG_NAME','Relative soil moisture')
          shumdiag1=val_exp
          IF ( ok_var(var_name) ) THEN
             CALL restget_p (rest_id, var_name, nbp_glo, nbdl, 1, kjit, .TRUE., shumdiag1, "gather", nbp_glo, index_g)
             IF (MINVAL(shumdiag1) < MAXVAL(shumdiag1) .OR. MAXVAL(shumdiag1) < val_exp) THEN
                shumdiag(:,:) = shumdiag1(:,:)
             ENDIF
          ENDIF
       ENDIF ! ldrestart_read

       !! 1.3 Initialize stomate's variables
       CALL slowproc_main (kjit, kjpij, kjpindex, dtradia, date0, &
            ldrestart_read, ldrestart_write, control%ok_co2, control%ok_stomate, &
            index, indexveg, lalo, neighbours, resolution, contfrac, soiltile, reinf_slope, &
            t2mdiag, t2mdiag, temp_sol, stempdiag, &
            vegstress, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
            deadleaf_cover, &
            assim_param, &
            lai, frac_age, height, veget, frac_nobio, njsc, veget_max, totfrac_nobio, qsintmax, &
            rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
            co2_flux, fco2_lu, temp_growth,&
            tdeep, hsdeep, snow, heat_Zimov, pb, &
            sfluxCH4_deep, sfluxCO2_deep, &
            thawed_humidity, depth_organic_soil, zz_deep, zz_coef_deep, &
            soilc_total,snowdz,snowrho)
       netco2flux(:) = zero
       DO jv = 2,nvm
          netco2flux(:) = netco2flux(:) + co2_flux(:,jv)*veget_max(:,jv)
       ENDDO
        
       !! 1.4 Initialize diffusion coefficients
       CALL diffuco_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, indexveg, indexlai, u, v, &
            & zlev, z0, roughheight, temp_sol, temp_air, temp_growth, rau, tq_cdrag, qsurf, qair, q2m, t2m, pb ,  &
            & rsol, evap_bare_lim, evapot, evapot_corr, snow, flood_frac, flood_res, frac_nobio, snow_nobio, totfrac_nobio, &
            & swnet, swdown, sinang, ccanopy, humrel, veget, veget_max, lai, qsintveg, qsintmax, assim_param, &
            & vbeta, valpha, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, gsmean, rveget, rstruct, cimean, gpp, &
            & lalo, neighbours, resolution, ptnlev1, precip_rain, frac_age, &
            & rest_id, hist_id, hist2_id)

       !! 1.5 Initialize remaining variables of energy budget
       CALL enerbil_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, &
            & index, indexveg, zlev, lwdown, swnet, epot_air, temp_air, u, v, petAcoef, petBcoef,&
            & qair, peqAcoef, peqBcoef, pb, rau, vbeta, valpha, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, &
            & emis, soilflx, soilcap, tq_cdrag, humrel, fluxsens, fluxlat, &
            & vevapp, transpir, transpot, vevapnu, vevapwet, vevapsno, vevapflo, t2mdiag, temp_sol, tsol_rad, &
            & temp_sol_new, qsurf, evapot, evapot_corr, rest_id, hist_id, hist2_id, &
            & precip_rain,snow,snowdz,snowrho,pgflux,snowflx,snowcap,temp_sol_add)
 

       !! 1.6 Initialize some hydrological variables from restart file
       IF (ldrestart_read) THEN
          
          var_name='runoff' ;
          CALL ioconf_setatt_p('UNITS', 'mm/d')
          CALL ioconf_setatt_p('LONG_NAME','Complete runoff')
          runoff1=val_exp
          IF ( ok_var(var_name) ) THEN
             CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., runoff1, "gather", nbp_glo, index_g)
             IF (MINVAL(runoff1) < MAXVAL(runoff1) .OR. MAXVAL(runoff1) < val_exp) THEN
                runoff(:) = runoff1(:)
             ENDIF
          ENDIF

          var_name='drainage' ;
          CALL ioconf_setatt_p('UNITS', 'mm/d')
          CALL ioconf_setatt_p('LONG_NAME','Deep drainage')
          drainage1=val_exp
          IF ( ok_var(var_name) ) THEN
             CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., drainage1, "gather", nbp_glo, index_g)
             IF (MINVAL(drainage1) < MAXVAL(drainage1) .OR. MAXVAL(drainage1) < val_exp) THEN
                drainage(:) = drainage1(:)
             ENDIF
          ENDIF

          IF ( ok_var("shumdiag") ) THEN
             IF (MINVAL(shumdiag1) < MAXVAL(shumdiag1) .OR. MAXVAL(shumdiag1) < val_exp) THEN
                shumdiag(:,:) = shumdiag1(:,:)
             ENDIF
          ENDIF
       ENDIF ! ldrestart_read
       
       !! 1.7 Initialize remaining hydrological variables
       IF ( .NOT. hydrol_cwrr ) THEN
          ! 1.7.1 Initialize remaining hydrological variables from Choisnel module (2 soil layers)
          CALL hydrolc_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, indexveg, control_in, &
               & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max,&
               & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,&
               & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, humrel, &
               & vegstress, rsol, drysoil_frac, evapot, evapot_corr, flood_frac, flood_res, shumdiag, litterhumdiag, &
               & soilcap, rest_id, hist_id, hist2_id, &
               & temp_air, pb, u, v, pgflux, &
               & snowrho,snowtemp,soiltemp,snowgrain,snowdz,snowheat,snowliq,&
               & grndflux,gtemp,gthick,gpkappa,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil, &
               & soilflxresid,snowflx,snowcap,pkappa_snow,lambda_snow,cgrnd_snow,dgrnd_snow,temp_sol_add)

          evap_bare_lim(:) = -un
          k_litt(:) = huit

          ! No specific calculation for shumdiag_perma. We assume it to shumdiag.
          shumdiag_perma(:,:)=shumdiag(:,:)
       ELSE
          !! 1.7.2 Initialize remaining hydrological variables from CWRR module (11 soil layers)
          !WRITE(numout,*) 'zd hydrol_main 1 ','snowtemp(1,:)',snowtemp(1,:)
          CALL hydrol_main (kjit, kjpindex, & !pss:+             
               & lalo, resolution, & !pss:-
               & dtradia, ldrestart_read, ldrestart_write, &
               & index, indexveg, indexsoil, indexlayer, indexnbdl, control_in, &
               & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max, njsc, &
               & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,&
               & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, humrel, &
               & vegstress, drysoil_frac, evapot, evapot_corr, evap_bare_lim, flood_frac, flood_res, &
               & shumdiag, shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, reinf_slope,  &
               & rest_id, hist_id, hist2_id, &
               & stempdiag, &
               & temp_air, pb, u, v, pgflux, &
               & snowrho,snowtemp,soiltemp,snowgrain,snowdz,snowheat,snowliq,&
               & grndflux,gtemp,gthick,gpkappa,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil, &
               & soilflxresid,snowflx,snowcap,pkappa_snow,lambda_snow,cgrnd_snow,dgrnd_snow,temp_sol_add, & !pss:+
               & drunoff_tot, fwet_out) !pss:-
          !WRITE(numout,*) 'zd hydrol_main 2 ','snowtemp(1,:)',snowtemp(1,:)

       ENDIF
       
       !! 1.8 Initialize some water balance variables from restart file
       IF (ldrestart_read) THEN
          
          IF ( ok_var("runoff") ) THEN
             IF (MINVAL(runoff1) < MAXVAL(runoff1) .OR. MAXVAL(runoff1) < val_exp) THEN
                runoff(:) = runoff1(:)
             ENDIF
          ENDIF

          IF ( ok_var("drainage") ) THEN
             IF (MINVAL(drainage1) < MAXVAL(drainage1) .OR. MAXVAL(drainage1) < val_exp) THEN
                drainage(:) = drainage1(:)
             ENDIF
          ENDIF
          
          IF ( ok_var("shumdiag") ) THEN
             IF (MINVAL(shumdiag1) < MAXVAL(shumdiag1) .OR. MAXVAL(shumdiag1) < val_exp) THEN
                shumdiag(:,:) = shumdiag1(:,:)
             ENDIF
          ENDIF

          IF ( ok_var("snowcap") ) THEN
             IF (MINVAL(snowcap1) < MAXVAL(snowcap1) .OR. MAXVAL(snowcap1) < val_exp) THEN
                snowcap(:) = snowcap1(:)
             ENDIF
          ENDIF

          IF ( ok_var("snowflx") ) THEN
             IF (MINVAL(snowflx1) < MAXVAL(snowflx1)  .OR. MAXVAL(snowflx1) < val_exp) THEN
                snowflx(:) = snowflx1(:)
             ENDIF
          ENDIF

       ENDIF ! ldrestart_read

        
       !! 1.9 Initialize surface parameters (emissivity, albedo and roughness)
       CALL condveg_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, &
            & lalo, neighbours, resolution, contfrac, veget, veget_max, frac_nobio, totfrac_nobio, &
            & zlev, snow, snow_age, snow_nobio, snow_nobio_age, &
            & drysoil_frac, height,  emis, albedo, z0, roughheight, rest_id, hist_id, hist2_id)
       
       !! 1.10 Initialization of soil thermodynamics
       !WRITE(numout,*) 'zd thermosoil_main 1 ','snowtemp(1,:)',snowtemp(1,:)
       CALL thermosoil_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, &
            & index,lalo, indexgrnd,indexnbdl, control_in,temp_sol_new, snow, soilcap, soilflx, &
            & shumdiag_perma, stempdiag, ptnlev1, rest_id, hist_id, hist2_id, &
            & soiltemp,pb,grndflux,snowrho,snowdz,snowtemp,gthick,gtemp,gpkappa,&
            & pkappa_snow,cgrnd_snow,dgrnd_snow,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil,&
            & thawed_humidity, depth_organic_soil, heat_Zimov, tdeep, hsdeep,&
            & soilc_total, veget_max)
       !WRITE(numout,*) 'zd thermosoil_main 2 ','snowtemp(1,:)',snowtemp(1,:)

       !! 1.11 Initialize some soil thermodynamics from restart file
       IF (ldrestart_read) THEN
          IF ( ok_var("soilcap") ) THEN
             IF (MINVAL(soilcap1) < MAXVAL(soilcap1) .OR. MAXVAL(soilcap1) < val_exp) THEN
                soilcap(:) = soilcap1(:)
             ENDIF
          ENDIF
         
          IF ( ok_var("soilflx") ) THEN
             IF (MINVAL(soilflx1) < MAXVAL(soilflx1)  .OR. MAXVAL(soilflx1) < val_exp) THEN
                soilflx(:) = soilflx1(:)
             ENDIF
          ENDIF
       ENDIF ! ldrestart_read

       !! 1.12 Initialize river routing
       IF ( river_routing .AND. nbp_glo .GT. 1) THEN
          !! 1.12.1 Initialize river routing
          CALL routing_main (kjit, kjpindex, dtradia, control_in, ldrestart_read, ldrestart_write, index, &
               & lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
               & drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
               & stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, &
               & rest_id, hist_id, hist2_id)
       ELSE
          !! 1.12.2 No routing, set variables to zero
          riverflow(:) = zero
          coastalflow(:) = zero
          returnflow(:) = zero
          reinfiltration(:) = zero
          irrigation(:) = zero
          flood_frac(:) = zero
          flood_res(:) = zero
       ENDIF
       
       !! 1.13 Write internal variables to output fields
       z0_out(:) = z0(:)
       emis_out(:) = emis(:)
       albedo_out(:,:) = albedo(:,:) 
       qsurf_out(:) = qsurf(:)

       !! 1.14 Deallocate memory
       DEALLOCATE(runoff1,drainage1,soilcap1,soilflx1,snowcap1,snowflx1)
       DEALLOCATE(shumdiag1)
       
       RETURN ! Out of sechiba main (should remain the last line!)
 
    ENDIF ! l_first_sechiba

 
!! 2. Computes SECHIBA's variables
    
    !! 2.1 Initialize variables at each time step
    CALL sechiba_var_init (kjpindex, rau, pb, temp_air) 

    !! 2.2 Compute diffusion coefficients
    CALL diffuco_main (kjit, kjpindex, dtradia, ldrestart_read, myfalse, index, indexveg, indexlai, u, v, &
!        & zlev, z0, roughheight, temp_sol, temp_air, rau, tq_cdrag, qsurf, qair, pb ,  &  !!?? could this line be deleted?
         & zlev, z0, roughheight, temp_sol, temp_air, temp_growth, rau, tq_cdrag, qsurf, qair, q2m, t2m, pb ,  &
         & rsol, evap_bare_lim, evapot, evapot_corr, snow, flood_frac, flood_res, frac_nobio, snow_nobio, totfrac_nobio, &
         & swnet, swdown, sinang, ccanopy, humrel, veget, veget_max, lai, qsintveg, qsintmax, assim_param, &
         & vbeta, valpha, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, gsmean, rveget, rstruct, cimean, gpp, &
         & lalo, neighbours, resolution, ptnlev1, precip_rain, frac_age, &
         & rest_id, hist_id, hist2_id)
   
    !! 2.3 Compute energy balance
    CALL enerbil_main (kjit, kjpindex, dtradia, ldrestart_read, myfalse, &
         & index, indexveg, zlev, lwdown, swnet, epot_air, temp_air, u, v, petAcoef, petBcoef, &
         & qair, peqAcoef, peqBcoef, pb, rau, vbeta, valpha, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, &
         & emis, soilflx, soilcap, tq_cdrag, humrel, fluxsens, fluxlat, &
         & vevapp, transpir, transpot, vevapnu, vevapwet, vevapsno, vevapflo, t2mdiag, temp_sol, tsol_rad, &
         & temp_sol_new, qsurf, evapot, evapot_corr, rest_id, hist_id, hist2_id,&
         & precip_rain,snow,snowdz,snowrho,pgflux,snowflx,snowcap,temp_sol_add)
    
    !! 2.4 Compute hydrology
    IF ( .NOT. hydrol_cwrr ) THEN
       ! 2.4.1 Water balance from Choisnel module (2 soil layers)
       CALL hydrolc_main (kjit, kjpindex, dtradia, ldrestart_read, myfalse, index, indexveg, control_in,&
            & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max,&
            & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,&
            & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, humrel, &
            & vegstress, rsol, drysoil_frac, evapot, evapot_corr, flood_frac, flood_res, shumdiag, litterhumdiag, &
            & soilcap, rest_id, hist_id, hist2_id, &
            & temp_air, pb, u, v, pgflux, &
            & snowrho,snowtemp,soiltemp,snowgrain,snowdz,snowheat,snowliq,&
            & grndflux,gtemp,gthick,gpkappa,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil, &
            & soilflxresid,snowflx,snowcap,pkappa_snow,lambda_snow,cgrnd_snow,dgrnd_snow,temp_sol_add)

       evap_bare_lim(:) = -un
       k_litt(:) = huit
       
       ! No specific calculation for shumdiag_perma. We assume it to shumdiag.
       shumdiag_perma(:,:)=shumdiag(:,:)
    ELSE
       !! 2.4.1 Water balance from CWRR module (11 soil layers)
       !WRITE(numout,*) 'zd hydrol_main 3 ','snowtemp(1,:)',snowtemp(1,:)
       CALL hydrol_main (kjit, kjpindex, & !pss:+             
            & lalo, resolution, & !pss:-
            & dtradia, ldrestart_read, myfalse, &
            & index, indexveg, indexsoil, indexlayer, indexnbdl, control_in, &
            & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max, njsc,&
            & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,&
            & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, humrel, &
            & vegstress, drysoil_frac, evapot, evapot_corr, evap_bare_lim, flood_frac, flood_res, &
            & shumdiag, shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, reinf_slope, &
            & rest_id, hist_id, hist2_id,&
            & stempdiag, &
            & temp_air, pb, u, v, pgflux, &
            & snowrho,snowtemp,soiltemp,snowgrain,snowdz,snowheat,snowliq,&
            & grndflux,gtemp,gthick,gpkappa,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil, &
            & soilflxresid,snowflx,snowcap,pkappa_snow,lambda_snow,cgrnd_snow,dgrnd_snow,temp_sol_add, & !pss:+
            & drunoff_tot, fwet_out) !pss:-
       !WRITE(numout,*) 'zd hydrol_main 4 ','snowtemp(1,:)',snowtemp(1,:)

       rsol(:) = -un

    ENDIF
     
    !! 2.5 Compute remaining components of the energy balance
    CALL enerbil_fusion (kjpindex, dtradia, tot_melt, soilcap, snowdz, &
                         temp_sol_new, fusion)
         
    
    !! 2.6 Compute surface variables (emissivity, albedo and roughness)
    CALL condveg_main (kjit, kjpindex, dtradia, ldrestart_read, myfalse, index,&
         & lalo, neighbours, resolution, contfrac, veget, veget_max, frac_nobio, totfrac_nobio, &
         & zlev, snow, snow_age, snow_nobio, snow_nobio_age, &
         & drysoil_frac, height,  emis, albedo, z0, roughheight, rest_id, hist_id, hist2_id)

    !! 2.7 Compute soil thermodynamics
    !WRITE(numout,*) 'zd thermosoil_main 3 ','snowtemp(1,:)',snowtemp(1,:)
    CALL thermosoil_main (kjit, kjpindex, dtradia, ldrestart_read, &
         & myfalse, index,lalo, indexgrnd,indexnbdl, &
         & control_in,temp_sol_new, snow, soilcap, soilflx, shumdiag_perma, stempdiag, &
         & ptnlev1, rest_id, hist_id, hist2_id, &
         & soiltemp,pb,grndflux,snowrho,snowdz,snowtemp,gthick,gtemp,gpkappa,&
         & pkappa_snow,cgrnd_snow,dgrnd_snow,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil,&
         & thawed_humidity, depth_organic_soil, heat_Zimov, tdeep, hsdeep,&
         & soilc_total, veget_max)
    !WRITE(numout,*) 'zd thermosoil_main 4 ','snowtemp(1,:)',snowtemp(1,:)


    !! 2.8 Compute river routing 
    IF ( river_routing .AND. nbp_glo .GT. 1) THEN
       !! 2.8.1 River routing
       CALL routing_main (kjit, kjpindex, dtradia, control_in, ldrestart_read, myfalse, index, &
            & lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
            & drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
            & stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, rest_id, hist_id, hist2_id)
    ELSE
       !! 2.8.2 No routing, set variables to zero
       riverflow(:) = zero
       coastalflow(:) = zero
       returnflow(:) = zero
       reinfiltration(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero
    ENDIF

    !! 2.9 Compute slow processes (i.e. 'daily' and annual time step)
    ! ::ok_co2 and ::ok_stomate are flags that determine whether the
    ! forcing files are written.
    CALL slowproc_main (kjit, kjpij, kjpindex, dtradia, date0, &
         ldrestart_read, myfalse, control%ok_co2, control%ok_stomate, &
         index, indexveg, lalo, neighbours, resolution, contfrac, soiltile, reinf_slope, &
         t2mdiag, t2mdiag, temp_sol, stempdiag, &
         vegstress, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
         deadleaf_cover, &
         assim_param, &
         lai, frac_age, height, veget, frac_nobio, njsc, veget_max, totfrac_nobio, qsintmax, &
         rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
         co2_flux, fco2_lu, temp_growth,&
         tdeep, hsdeep, snow, heat_Zimov, pb, &
         sfluxCH4_deep, sfluxCO2_deep, &
         thawed_humidity, depth_organic_soil, zz_deep, zz_coef_deep, &
         soilc_total,snowdz,snowrho)
   
    !! 2.9 Compute global CO2 flux
    netco2flux(:) = zero
    DO jv = 2,nvm
       netco2flux(:) = netco2flux(:) + co2_flux(:,jv)*veget_max(:,jv)
    ENDDO
 
    !! 2.10 Update the temperature (temp_sol) with newly computed values
    !WRITE(numout,*) 'zd sechiba_end 1 ','snowtemp(1,:)',snowtemp(1,:)
    CALL sechiba_end (kjpindex, dtradia, temp_sol_new, snowtemp, snowdz, &
                      temp_sol)
    !WRITE(numout,*) 'zd sechiba_end 2 ','snowtemp(1,:)',snowtemp(1,:)

   
    !! 2.11 Write internal variables to output fields
    z0_out(:) = z0(:)
    emis_out(:) = emis(:)
    albedo_out(:,:) = albedo(:,:) 
    qsurf_out(:) = qsurf(:)
 
    !! 2.12 Write global variables to history files
    sum_treefrac(:) = zero
    sum_grassfrac(:) = zero
    sum_cropfrac(:) = zero
    DO jv = 2, nvm 
       IF (is_tree(jv) .AND. natural(jv)) THEN
          sum_treefrac(:) = sum_treefrac(:) + veget_max(:,jv)
       ELSE IF ((.NOT. is_tree(jv))  .AND. natural(jv)) THEN
          sum_grassfrac(:) = sum_grassfrac(:) + veget_max(:,jv)
       ELSE 
          sum_cropfrac = sum_cropfrac(:) + veget_max(:,jv)
       ENDIF
    ENDDO          
    !

    CALL xios_orchidee_send_field("evapnu",vevapnu*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("snow",snow)
    CALL xios_orchidee_send_field("snowage",snow_age)
    CALL xios_orchidee_send_field("snownobio",snow_nobio)
    CALL xios_orchidee_send_field("snownobioage",snow_nobio_age)
    CALL xios_orchidee_send_field("reinf_slope",reinf_slope)
    CALL xios_orchidee_send_field("soilindex",REAL(njsc, r_std))
    CALL xios_orchidee_send_field("vegetfrac",veget)
    CALL xios_orchidee_send_field("maxvegetfrac",veget_max)
    CALL xios_orchidee_send_field("nobiofrac",frac_nobio)
    CALL xios_orchidee_send_field("soiltile",soiltile)
    CALL xios_orchidee_send_field("rstruct",rstruct)
    CALL xios_orchidee_send_field("gpp",gpp/dt_sechiba)
    CALL xios_orchidee_send_field("nee",co2_flux/dt_sechiba)
    CALL xios_orchidee_send_field("drysoil_frac",drysoil_frac)
    CALL xios_orchidee_send_field("evapflo",vevapflo*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("evapflo_alma",vevapflo/dt_sechiba)
    CALL xios_orchidee_send_field("k_litt",k_litt)
    CALL xios_orchidee_send_field("beta",vbeta)
    CALL xios_orchidee_send_field("vbeta1",vbeta1)
    CALL xios_orchidee_send_field("vbeta2",vbeta2)
    CALL xios_orchidee_send_field("vbeta3",vbeta3)
    CALL xios_orchidee_send_field("vbeta4",vbeta4)
    CALL xios_orchidee_send_field("vbeta5",vbeta5)
    CALL xios_orchidee_send_field("gsmean",gsmean)
    CALL xios_orchidee_send_field("cimean",cimean)
    CALL xios_orchidee_send_field("rveget",rveget)
    CALL xios_orchidee_send_field("rsol",rsol)

! Note that 0.0005555556 is used instead of one_day/dt temporary to have the same 
! order of cacluation as using IOIPSL for compairing of resluts. This will be changed in futur update.   
    histvar(:)=SUM(vevapwet(:,:),dim=2)
    CALL xios_orchidee_send_field("evspsblveg",histvar/dt_sechiba)
    histvar(:)= vevapnu(:)+vevapsno(:)
    CALL xios_orchidee_send_field("evspsblsoi",histvar/dt_sechiba)
    histvar(:)=SUM(transpir(:,:),dim=2)
    CALL xios_orchidee_send_field("tran",histvar/dt_sechiba)
    histvar(:)= sum_treefrac(:)*100*contfrac(:)
    CALL xios_orchidee_send_field("treeFrac",histvar)
    histvar(:)= sum_grassfrac(:)*100*contfrac(:)
    CALL xios_orchidee_send_field("grassFrac",histvar)
    histvar(:)= sum_cropfrac(:)*100*contfrac(:)
    CALL xios_orchidee_send_field("cropFrac",histvar)
    histvar(:)=veget_max(:,1)*100*contfrac(:)
    CALL xios_orchidee_send_field("baresoilFrac",histvar)
    histvar(:)=SUM(frac_nobio(:,1:nnobio),dim=2)*100*contfrac(:)
    CALL xios_orchidee_send_field("residualFrac",histvar)

    CALL xios_orchidee_send_field("tsol_rad",tsol_rad-273.15)
    CALL xios_orchidee_send_field("qsurf",qsurf)
    CALL xios_orchidee_send_field("albedo",albedo)
    CALL xios_orchidee_send_field("emis",emis)
    CALL xios_orchidee_send_field("z0",z0)
    CALL xios_orchidee_send_field("roughheight",roughheight)
    CALL xios_orchidee_send_field("lai",lai)
    CALL xios_orchidee_send_field("subli",vevapsno*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("vevapnu",vevapnu*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("vevapnu_alma",vevapnu/dt_sechiba)
    CALL xios_orchidee_send_field("transpir",transpir*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("inter",vevapwet*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("Qf",fusion)
    histvar(:)=zero
    DO jv=1,nvm
      histvar(:) = histvar(:) + vevapwet(:,jv)
    ENDDO
    CALL xios_orchidee_send_field("ECanop",histvar/dt_sechiba)
    histvar(:)=zero
    DO jv=1,nvm
      histvar(:) = histvar(:) + transpir(:,jv)
    ENDDO
    CALL xios_orchidee_send_field("TVeg",histvar/dt_sechiba)
    CALL xios_orchidee_send_field("ACond",tq_cdrag)

    IF ( .NOT. almaoutput ) THEN
       ! Write history file in IPSL-format
       CALL histwrite_p(hist_id, 'beta', kjit, vbeta, kjpindex, index)
       CALL histwrite_p(hist_id, 'z0', kjit, z0, kjpindex, index)
       CALL histwrite_p(hist_id, 'roughheight', kjit, roughheight, kjpindex, index)
       CALL histwrite_p(hist_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
       CALL histwrite_p(hist_id, 'lai', kjit, lai, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'subli', kjit, vevapsno, kjpindex, index)
       CALL histwrite_p(hist_id, 'evapnu', kjit, vevapnu, kjpindex, index)
       CALL histwrite_p(hist_id, 'transpir', kjit, transpir, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'inter', kjit, vevapwet, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'vbeta1', kjit, vbeta1, kjpindex, index)
       CALL histwrite_p(hist_id, 'vbeta2', kjit, vbeta2, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'vbeta3', kjit, vbeta3, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'vbeta4', kjit, vbeta4, kjpindex, index)    
       CALL histwrite_p(hist_id, 'vbeta5', kjit, vbeta5, kjpindex, index)    
       CALL histwrite_p(hist_id, 'drysoil_frac', kjit, drysoil_frac, kjpindex, index)
       CALL histwrite_p(hist_id, 'rveget', kjit, rveget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'rstruct', kjit, rstruct, kjpindex*nvm, indexveg)
       IF ( .NOT. control_in%hydrol_cwrr ) THEN
          CALL histwrite_p(hist_id, 'rsol', kjit, rsol, kjpindex, index)
       ENDIF
       CALL histwrite_p(hist_id, 'snow', kjit, snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'snowage', kjit, snow_age, kjpindex, index)
       CALL histwrite_p(hist_id, 'snownobio', kjit, snow_nobio, kjpindex*nnobio, indexnobio)
       CALL histwrite_p(hist_id, 'snownobioage', kjit, snow_nobio_age, kjpindex*nnobio, indexnobio)

       IF (ok_explicitsnow) THEN
          CALL histwrite_p(hist_id, 'grndflux', kjit, grndflux, kjpindex,index)
          CALL histwrite_p(hist_id, 'snowtemp',kjit,snowtemp,kjpindex*nsnow,indexsnow)
          CALL histwrite_p(hist_id, 'soiltemp',kjit,soiltemp,kjpindex*ngrnd,indexgrnd)
          CALL histwrite_p(hist_id, 'snowliq', kjit,snowliq,kjpindex*nsnow,indexsnow)
          CALL histwrite_p(hist_id, 'snowdz', kjit,snowdz,kjpindex*nsnow,indexsnow)
          CALL histwrite_p(hist_id, 'snowrho', kjit,snowrho,kjpindex*nsnow,indexsnow)
          CALL histwrite_p(hist_id, 'snowgrain',kjit,snowgrain,kjpindex*nsnow,indexsnow)
          CALL histwrite_p(hist_id, 'snowheat',kjit,snowheat,kjpindex*nsnow,indexsnow)
       END IF

       CALL histwrite_p(hist_id, 'pgflux',kjit,pgflux,kjpindex,index)
       CALL histwrite_p(hist_id, 'soiltile',  kjit, soiltile, kjpindex*nstm, indexsoil)
       !
       IF ( control_in%hydrol_cwrr ) THEN
          CALL histwrite_p(hist_id, 'soilindex',  kjit, REAL(njsc, r_std), kjpindex, index)
          CALL histwrite_p(hist_id, 'reinf_slope',  kjit, reinf_slope, kjpindex, index)
          CALL histwrite_p(hist_id, 'k_litt', kjit, k_litt, kjpindex, index)
       ENDIF
       IF ( control_in%do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'evapflo', kjit, vevapflo, kjpindex, index)
          CALL histwrite_p(hist_id, 'flood_frac', kjit, flood_frac, kjpindex, index)
       ENDIF
       IF ( control%ok_co2 ) THEN
          CALL histwrite_p(hist_id, 'gsmean', kjit, gsmean, kjpindex*nvm, indexveg)    
          CALL histwrite_p(hist_id, 'gpp', kjit, gpp, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist_id, 'cimean', kjit, cimean, kjpindex*nvm, indexveg)    
       ENDIF
       IF ( control%ok_stomate ) THEN
          CALL histwrite_p(hist_id, 'nee', kjit, co2_flux, kjpindex*nvm, indexveg)    
       ENDIF

       histvar(:)=SUM(vevapwet(:,:),dim=2)
       CALL histwrite_p(hist_id, 'evspsblveg', kjit, histvar, kjpindex, index)

       histvar(:)= vevapnu(:)+vevapsno(:)
       CALL histwrite_p(hist_id, 'evspsblsoi', kjit, histvar, kjpindex, index)

       histvar(:)=SUM(transpir(:,:),dim=2)
       CALL histwrite_p(hist_id, 'tran', kjit, histvar, kjpindex, index)

       histvar(:)= sum_treefrac(:)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'treeFrac', kjit, histvar, kjpindex, index) 

       histvar(:)= sum_grassfrac(:)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'grassFrac', kjit, histvar, kjpindex, index) 

       histvar(:)= sum_cropfrac(:)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'cropFrac', kjit, histvar, kjpindex, index)

       histvar(:)=veget_max(:,1)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'baresoilFrac', kjit, histvar, kjpindex, index)

       histvar(:)=SUM(frac_nobio(:,1:nnobio),dim=2)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'residualFrac', kjit, histvar, kjpindex, index)
    ELSE
       ! Write history file in ALMA format 
       CALL histwrite_p(hist_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
       CALL histwrite_p(hist_id, 'Qf', kjit, fusion, kjpindex, index)
       CALL histwrite_p(hist_id, 'ESoil', kjit, vevapnu, kjpindex, index)
       CALL histwrite_p(hist_id, 'EWater', kjit, vevapflo, kjpindex, index)
       CALL histwrite_p(hist_id, 'SWE', kjit, snow, kjpindex, index)
       histvar(:)=zero
       DO jv=1,nvm
          histvar(:) = histvar(:) + transpir(:,jv)
       ENDDO
       CALL histwrite_p(hist_id, 'TVeg', kjit, histvar, kjpindex, index)
       histvar(:)=zero
       DO jv=1,nvm
          histvar(:) = histvar(:) + vevapwet(:,jv)
       ENDDO
       CALL histwrite_p(hist_id, 'ECanop', kjit, histvar, kjpindex, index)
       CALL histwrite_p(hist_id, 'ACond', kjit, tq_cdrag, kjpindex, index)
       CALL histwrite_p(hist_id, 'SnowFrac', kjit, vbeta1, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'Z0', kjit, z0, kjpindex, index)
       CALL histwrite_p(hist_id, 'EffectHeight', kjit, roughheight, kjpindex, index)
       !
       IF ( control_in%do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'Qflood', kjit, vevapflo, kjpindex, index)
          CALL histwrite_p(hist_id, 'FloodFrac', kjit, flood_frac, kjpindex, index)
       ENDIF
       !
       IF ( control%ok_co2 ) THEN
          CALL histwrite_p(hist_id, 'GPP', kjit, gpp, kjpindex*nvm, indexveg)
       ENDIF
       IF ( control%ok_stomate ) THEN
             CALL histwrite_p(hist_id, 'NEE', kjit, co2_flux, kjpindex*nvm, indexveg)    
       ENDIF
    ENDIF ! almaoutput
    
    !! 2.13 Write additional output file with higher frequency
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
          ! Write history file in IPSL-format
          CALL histwrite_p(hist2_id, 'tsol_rad', kjit, tsol_rad, kjpindex, index)
          CALL histwrite_p(hist2_id, 'qsurf', kjit, qsurf, kjpindex, index)
          CALL histwrite_p(hist2_id, 'albedo', kjit, albedo, kjpindex*2, indexalb)
          CALL histwrite_p(hist2_id, 'emis', kjit, emis, kjpindex, index)
          CALL histwrite_p(hist2_id, 'beta', kjit, vbeta, kjpindex, index)
          CALL histwrite_p(hist2_id, 'z0', kjit, z0, kjpindex, index)
          CALL histwrite_p(hist2_id, 'roughheight', kjit, roughheight, kjpindex, index)
          CALL histwrite_p(hist2_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'lai', kjit, lai, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'subli', kjit, vevapsno, kjpindex, index)
          IF ( control_in%do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'vevapflo', kjit, vevapflo, kjpindex, index)
             CALL histwrite_p(hist2_id, 'flood_frac', kjit, flood_frac, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'vevapnu', kjit, vevapnu, kjpindex, index)
          CALL histwrite_p(hist2_id, 'transpir', kjit, transpir, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'inter', kjit, vevapwet, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'vbeta1', kjit, vbeta1, kjpindex, index)
          CALL histwrite_p(hist2_id, 'vbeta2', kjit, vbeta2, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'vbeta3', kjit, vbeta3, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'vbeta4', kjit, vbeta4, kjpindex, index)    
          CALL histwrite_p(hist2_id, 'vbeta5', kjit, vbeta5, kjpindex, index)    
          CALL histwrite_p(hist2_id, 'drysoil_frac', kjit, drysoil_frac, kjpindex, index)
          CALL histwrite_p(hist2_id, 'rveget', kjit, rveget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'rstruct', kjit, rstruct, kjpindex*nvm, indexveg)
          IF ( .NOT. control_in%hydrol_cwrr ) THEN
             CALL histwrite_p(hist2_id, 'rsol', kjit, rsol, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'snow', kjit, snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowage', kjit, snow_age, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snownobio', kjit, snow_nobio, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'snownobioage', kjit, snow_nobio_age, kjpindex*nnobio, indexnobio)
          !
          IF (  control_in%hydrol_cwrr ) THEN
             CALL histwrite_p(hist2_id, 'soilindex',  kjit, REAL(njsc, r_std), kjpindex, index)
             CALL histwrite_p(hist2_id, 'reinf_slope',  kjit, reinf_slope, kjpindex, index)
          ENDIF
          !
          IF ( control%ok_co2 ) THEN
             CALL histwrite_p(hist2_id, 'gsmean', kjit, gsmean, kjpindex*nvm, indexveg)    
             CALL histwrite_p(hist2_id, 'gpp', kjit, gpp, kjpindex*nvm, indexveg)
             CALL histwrite_p(hist2_id, 'cimean', kjit, cimean, kjpindex*nvm, indexveg)    
          ENDIF
          IF ( control%ok_stomate ) THEN
             CALL histwrite_p(hist2_id, 'nee', kjit, co2_flux, kjpindex*nvm, indexveg)    
          ENDIF
       ELSE
          ! Write history file in ALMA format
          CALL histwrite_p(hist2_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'Qf', kjit, fusion, kjpindex, index)
          CALL histwrite_p(hist2_id, 'ESoil', kjit, vevapnu, kjpindex, index)
          IF ( control_in%do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'EWater', kjit, vevapflo, kjpindex, index)
             CALL histwrite_p(hist2_id, 'FloodFrac', kjit, flood_frac, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'SWE', kjit, snow, kjpindex, index)
          histvar(:)=zero
          DO jv=1,nvm
             histvar(:) = histvar(:) + transpir(:,jv)
          ENDDO
          CALL histwrite_p(hist2_id, 'TVeg', kjit, histvar, kjpindex, index)
          histvar(:)=zero
          DO jv=1,nvm
             histvar(:) = histvar(:) + vevapwet(:,jv)
          ENDDO
          CALL histwrite_p(hist2_id, 'ECanop', kjit, histvar, kjpindex, index)
          CALL histwrite_p(hist2_id, 'ACond', kjit, tq_cdrag, kjpindex, index)
          CALL histwrite_p(hist2_id, 'SnowFrac', kjit, vbeta1, kjpindex, index)
          IF ( control%ok_co2 ) THEN
             CALL histwrite_p(hist2_id, 'GPP', kjit, gpp, kjpindex*nvm, indexveg)
          ENDIF
          IF ( control%ok_stomate ) THEN
             CALL histwrite_p(hist2_id, 'NEE', kjit, co2_flux, kjpindex*nvm, indexveg)    
          ENDIF
       ENDIF ! almaoutput
    ENDIF ! hist2_id

!! 3. Write restart file for the next simulation from SECHIBA and other modules
    IF (ldrestart_write) THEN

       ! Only called during the last run
       IF (long_print) WRITE (numout,*) ' we have to write a restart file '
       !! 3.1 Call diffuco_main to write restart files
       CALL diffuco_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, indexveg, indexlai, u, v, &
!           & zlev, z0, roughheight, temp_sol, temp_air, rau, tq_cdrag, qsurf, qair, pb ,  & !!?? Could this line be deleted?
            & zlev, z0, roughheight, temp_sol, temp_air, temp_growth, rau, tq_cdrag, qsurf, qair, q2m, t2m, pb ,  &
            & rsol, evap_bare_lim, evapot, evapot_corr, snow, flood_frac, flood_res, frac_nobio, snow_nobio, totfrac_nobio, &
            & swnet, swdown, sinang, ccanopy, humrel, veget, veget_max, lai, qsintveg, qsintmax, assim_param, &
            & vbeta, valpha, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, gsmean, rveget, rstruct, cimean, gpp, &
            & lalo, neighbours, resolution, ptnlev1, precip_rain, frac_age, &
            & rest_id, hist_id, hist2_id)

 
       !! 3.2 Call energy budget to write restart files
       CALL enerbil_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, &
            & index, indexveg, zlev, lwdown, swnet, epot_air, temp_air, u, v, petAcoef, petBcoef,&
            & qair, peqAcoef, peqBcoef, pb, rau, vbeta, valpha, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, &
            & emis, soilflx, soilcap, tq_cdrag, humrel, fluxsens, fluxlat, &
            & vevapp, transpir, transpot, vevapnu, vevapwet, vevapsno, vevapflo, t2mdiag, temp_sol, tsol_rad, &
            & temp_sol_new, qsurf, evapot, evapot_corr, rest_id, hist_id, hist2_id, &
            & precip_rain,snow,snowdz,snowrho,pgflux,snowflx,snowcap,temp_sol_add)

       !! 3.3 Call hydrology to write restart files
       IF ( .NOT. hydrol_cwrr ) THEN
          !! 3.3.1 Call water balance from Choisnel module (2 soil layers) to write restart file
          CALL hydrolc_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, indexveg, control_in, &
               & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max,&
               & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,&
               & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, &
               & humrel, vegstress, rsol, drysoil_frac, evapot, evapot_corr, flood_frac, flood_res, shumdiag, litterhumdiag, &
               & soilcap, rest_id, hist_id, hist2_id, &
               & temp_air, pb, u, v, pgflux, &
               & snowrho,snowtemp,soiltemp,snowgrain,snowdz,snowheat,snowliq,&
               & grndflux,gtemp,gthick,gpkappa,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil, &
               & soilflxresid,snowflx,snowcap,pkappa_snow,lambda_snow,cgrnd_snow,dgrnd_snow,temp_sol_add)

          evap_bare_lim(:) = -un
          k_litt(:) = huit
          shumdiag_perma(:,:)=shumdiag(:,:)
       ELSE
          !! 3.3.2 Call water balance from CWRR module (11 soil layers) to write restart file
          !WRITE(numout,*) 'zd hydrol_main 5 ','snowtemp(1,:)',snowtemp(1,:)
          CALL hydrol_main (kjit, kjpindex, & !pss:+             
               & lalo, resolution, & !pss:-
               & dtradia, ldrestart_read, ldrestart_write, &
               & index, indexveg, indexsoil, indexlayer, indexnbdl, control_in, &
               & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max, njsc, &
               & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,&
               & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, humrel, &
               & vegstress, drysoil_frac, evapot, evapot_corr, evap_bare_lim, flood_frac, flood_res, &
               & shumdiag,shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, reinf_slope,  &
               & rest_id, hist_id, hist2_id, &
               & stempdiag, &
               & temp_air, pb, u, v, pgflux, &
               & snowrho,snowtemp,soiltemp,snowgrain,snowdz,snowheat,snowliq,&
               & grndflux,gtemp,gthick,gpkappa,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil, &
               & soilflxresid,snowflx,snowcap,pkappa_snow,lambda_snow,cgrnd_snow,dgrnd_snow,temp_sol_add, & !pss:+
               & drunoff_tot, fwet_out) !pss:-
          !WRITE(numout,*) 'zd hydrol_main 6 ','snowtemp(1,:)',snowtemp(1,:)

          rsol(:) = -un
       ENDIF ! hydrol_cwrr
   
       !! 3.4 Call condveg to write surface variables to restart files
       CALL condveg_main (kjit, kjpindex, dtradia, ldrestart_read, ldrestart_write, index, &
            & lalo, neighbours, resolution, contfrac, veget, veget_max, frac_nobio, totfrac_nobio, &
            & zlev, snow, snow_age, snow_nobio, snow_nobio_age, &
            & drysoil_frac, height,  emis, albedo, z0, roughheight, rest_id, hist_id, hist2_id)
 
       !! 3.5 Call soil thermodynamic to write restart files
       !WRITE(numout,*) 'zd thermosoil_main 5 ','snowtemp(1,:)',snowtemp(1,:)
       CALL thermosoil_main (kjit, kjpindex, dtradia, &
            & ldrestart_read, ldrestart_write, index,lalo, indexgrnd,indexnbdl, &
            & control_in,temp_sol_new, snow, soilcap, soilflx, shumdiag_perma, &
            & stempdiag, ptnlev1, rest_id, hist_id, hist2_id, &
            & soiltemp,pb,grndflux,snowrho,snowdz,snowtemp,gthick,gtemp,gpkappa,&
            & pkappa_snow,cgrnd_snow,dgrnd_snow,zdz1_soil,zdz2_soil,cgrnd_soil,dgrnd_soil,&
            & thawed_humidity, depth_organic_soil, heat_Zimov, tdeep, hsdeep,&
            & soilc_total, veget_max)
       !WRITE(numout,*) 'zd thermosoil_main 6 ','snowtemp(1,:)',snowtemp(1,:)

       !! 3.6 Add river routing to restart files  
       IF ( river_routing .AND. nbp_glo .GT. 1) THEN
          !! 3.6.1 Call river routing to write restart files 
          CALL routing_main (kjit, kjpindex, dtradia, control_in, ldrestart_read, ldrestart_write, index, &
               & lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
               & drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, &
               & stempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, &
               & rest_id, hist_id, hist2_id)
       ELSE
          !! 3.6.2 No routing, set variables to zero
          riverflow(:) = zero
          coastalflow(:) = zero
          reinfiltration(:) = zero
          returnflow(:) = zero
          irrigation(:) = zero
          flood_frac(:) = zero
          flood_res(:) = zero
       ENDIF ! river_routing

       !! 3.7 Call slowproc_main to add 'daily' and annual variables to restart file
       CALL slowproc_main (kjit, kjpij, kjpindex, dtradia, date0, &
            ldrestart_read, ldrestart_write, control%ok_co2, control%ok_stomate, &
            index, indexveg, lalo, neighbours, resolution, contfrac, soiltile, reinf_slope, &
            t2mdiag, t2mdiag, temp_sol, stempdiag, &
            vegstress, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
            deadleaf_cover, &
            assim_param, &
            lai, frac_age, height, veget, frac_nobio, njsc, veget_max, totfrac_nobio, qsintmax, &
            rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
            co2_flux, fco2_lu, temp_growth,&
            tdeep, hsdeep, snow, heat_Zimov, pb, &
            sfluxCH4_deep, sfluxCO2_deep, &
            thawed_humidity, depth_organic_soil, zz_deep, zz_coef_deep, &
            soilc_total,snowdz,snowrho)
       ! Compute global CO2 flux !*
       netco2flux(:) = zero
       DO jv = 2,nvm
          netco2flux(:) = netco2flux(:) + co2_flux(:,jv)*veget_max(:,jv)
       ENDDO

       !! 3.8 Add some remaining variables to the restart file
       var_name= 'shumdiag'  
       CALL restput_p(rest_id, var_name, nbp_glo,   nbdl, 1, kjit,  shumdiag, 'scatter',  nbp_glo, index_g)
       var_name= 'runoff'  
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  runoff, 'scatter',  nbp_glo, index_g)
       var_name= 'drainage'  
       CALL restput_p(rest_id, var_name, nbp_glo,   1, 1, kjit,  drainage, 'scatter',  nbp_glo, index_g)
    
    END IF ! ldrestart_write 


    IF (long_print) WRITE (numout,*) ' sechiba_main done '

  END SUBROUTINE sechiba_main

  
!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_init
!!
!>\BRIEF        Dynamic allocation of the variables, the dimensions of the 
!! variables are determined by user-specified settings. If ::ldrestart_read is
!! true, initial varaible values are read from the restart file.
!! 
!! DESCRIPTION  : The domain size (:: kjpindex) is used to allocate the correct
!! dimensions to all variables in sechiba. Depending on the variable, its 
!! dimensions are also determined by the number of PFT's (::nvm), number of 
!! soil types (::nstm), number of non-vegetative surface types (::nnobio),
!! number of soil levels (::ngrnd), number of soil layers in the hydrological 
!! model (i.e. cwrr) (::nslm). Values for these variables are set in
!! constantes_soil.f90 and constantes_veg.f90.\n
!!
!! Memory is allocated for all Sechiba variables and new indexing tables
!! are build making use of both (::kjpij) and (::kjpindex). New indexing tables 
!! are needed because a single pixel can contain several PFTs, soil types, etc.
!! The new indexing tables have separate indices for the different
!! PFTs, soil types, etc.\n
!!
!! The routine uses the flag ::ldrestart_read to retrieve initial variable 
!! settings from the restart file.\n
!!  
!! ??Flags are getting swapped by it is not entirly clear why?? 
!!
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): Strictly speaking the subroutine has no output 
!! variables. However, the routine allocates memory and builds new indexing 
!! variables for later use.
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================ 

  SUBROUTINE sechiba_init (kjit, ldrestart_read, kjpij, kjpindex, index, rest_id, control_in, lalo)

!! 0.1 Input variables
 
    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number (unitless)
    INTEGER(i_std), INTENT (in)                         :: kjpij              !! Total size of the un-compressed grid (unitless)
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size - terrestrial pixels only (unitless)
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! _Restart_ file identifier (unitless)
    TYPE(control_type), INTENT(in)                      :: control_in         !! Flags that (de)activate parts of the model
    LOGICAL,INTENT (in)                                 :: ldrestart_read     !! Logical for restart file to read (true/false)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map (unitless)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo               !! Geographical coordinates (latitude,longitude) 
                                                                              !! for pixels (degrees)
!! 0.2 Output variables

!! 0.3 Modified variables

!! 0.4 Local variables

    INTEGER(i_std)                                      :: ier                !! Check errors in memory allocation (unitless)
    INTEGER(i_std)                                      :: ji, jv             !! Indeces (unitless)
!_ ==============================================================================================================================

!! 1. Initialize variables 
    
    ! Dynamic allocation with user-specified dimensions on first call
    IF (l_first_sechiba) THEN 
       l_first_sechiba=.FALSE.
    ELSE 
       WRITE (numout,*) ' l_first_sechiba false . we stop '
       STOP 'sechiba_init'
    ENDIF


    !! 1.1 Initialize 3D vegetation indexation table
    ALLOCATE (indexveg(kjpindex*nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexveg allocation. We stop. We need kjpindex words = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (indexlai(kjpindex*(nlai+1)),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexlai allocation. We stop. We need kjpindex words = ',kjpindex*(nlai+1)
       STOP 'sechiba_init'
    END IF

    ALLOCATE (indexsoil(kjpindex*nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexsoil allocation. We stop. We need kjpindex words = ',kjpindex*nstm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (indexnobio(kjpindex*nnobio),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexnobio allocation. We stop. We need kjpindex words = ',kjpindex*nnobio
       STOP 'sechiba_init'
    END IF

    ALLOCATE (indexgrnd(kjpindex*ngrnd),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexgrnd allocation. We stop. We need kjpindex words = ',kjpindex*ngrnd
       STOP 'sechiba_init'
    END IF

    ALLOCATE (indexsnow(kjpindex*nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexsnow allocation. We stop. We need kjpindex words = ',kjpindex*nsnow
       STOP 'sechiba_init'
    END IF

    ALLOCATE (indexlayer(kjpindex*nslm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexlayer allocation. We stop. We need kjpindex words = ',kjpindex*nslm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (indexnbdl(kjpindex*nbdl),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexgrnd allocation. We stop. We need kjpindex words = ',kjpindex*nbdl
       STOP 'sechiba_init'
    END IF

    ALLOCATE (indexalb(kjpindex*2),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in indexalb allocation. We stop. We need kjpindex words = ',kjpindex*2
       STOP 'sechiba_init'
    END IF

    !! 1.2  Initialize 1D array allocation with restartable value
    IF (long_print) WRITE (numout,*) 'Allocation of 1D variables. We need for each kjpindex words = ',kjpindex
    ALLOCATE (flood_res(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in flood_res allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    flood_res(:) = undef_sechiba

    IF (long_print) WRITE (numout,*) 'Allocation of 1D variables. We need for each kjpindex words = ',kjpindex
    ALLOCATE (flood_frac(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in flood_frac allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    flood_frac(:) = undef_sechiba

    IF (long_print) WRITE (numout,*) 'Allocation of 1D variables. We need for each kjpindex words = ',kjpindex
    ALLOCATE (snow(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snow allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    snow(:) = undef_sechiba

    ALLOCATE (snow_age(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snow_age allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    snow_age(:) = undef_sechiba

    ALLOCATE (drysoil_frac(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in drysoil_frac allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    drysoil_frac(:) = zero

    ALLOCATE (rsol(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in rsol allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (evap_bare_lim(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in evap_bare_lim allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (evapot(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in evapot allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    evapot(:) = undef_sechiba

    ALLOCATE (evapot_corr(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in evapot_corr allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (humrel(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in humrel allocation. We stop. we need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    humrel(:,:) = undef_sechiba

    ALLOCATE (vegstress(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vegstress allocation. We stop. we need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    vegstress(:,:) = undef_sechiba

    ALLOCATE (njsc(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in njsc allocation. We stop. we need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    njsc(:)= undef_int

    ALLOCATE (soiltile(kjpindex,nstm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soiltile allocation. We stop. we need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    soiltile(:,:)=undef_sechiba

    ALLOCATE (reinf_slope(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in reinf_slope allocation. We stop. we need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    reinf_slope(:)=undef_sechiba

    ALLOCATE (vbeta1(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vbeta1 allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (vbeta4(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vbeta4 allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (vbeta5(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vbeta5 allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (soilcap(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soilcap allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (soilflx(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soilflx allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (temp_sol(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in temp_sol allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    temp_sol(:) = undef_sechiba

    ALLOCATE (qsurf(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in qsurf allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    qsurf(:) = undef_sechiba

    !! 1.3 Initialize 2D array allocation with restartable value
    ALLOCATE (qsintveg(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in qsintveg allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm,' = ',kjpindex*nvm 
       STOP 'sechiba_init'
    END IF
    qsintveg(:,:) = undef_sechiba

    ALLOCATE (vbeta2(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vbeta2 allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (vbeta3(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vbeta3 allocation. We stop.We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (vbeta3pot(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vbeta3pot allocation. We stop.We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (gsmean(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gsmean allocation. We stop.We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (cimean(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in cimean allocation. We stop.We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (gpp(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gpp allocation. We stop.We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    gpp(:,:) = undef_sechiba

 
    ALLOCATE (temp_growth(kjpindex),stat=ier) 
    IF (ier.NE.0) THEN 
       WRITE (numout,*) ' error in temp_growth allocation. We stop.We need kjpindex words = ',& 
            & kjpindex,' = ',kjpindex 
       STOP 'sechiba_init' 
    END IF
    temp_growth(:) = undef_sechiba 

    ALLOCATE (veget(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in veget allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    veget(:,:)=undef_sechiba

    ALLOCATE (veget_max(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in veget_max allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    veget_max(:,:)=undef_sechiba

    ALLOCATE (lai(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in lai allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    lai(:,:)=undef_sechiba

    ALLOCATE (frac_age(kjpindex,nvm,nleafages),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in frac_age allocation. We stop. We need kjpindex x nvm words = ',&
             & kjpindex,' x ' ,nvm, 'x',nleafages,' = ',kjpindex*nvm*nleafages
        STOP 'sechiba_init'
    END IF
    frac_age(:,:,:)=undef_sechiba

    ALLOCATE (height(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in height allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    height(:,:)=undef_sechiba

    ALLOCATE (frac_nobio(kjpindex,nnobio),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in frac_nobio allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    frac_nobio(:,:) = undef_sechiba

    ALLOCATE (albedo(kjpindex,2),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in albedo allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (snow_nobio(kjpindex,nnobio),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snow_nobio allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    snow_nobio(:,:) = undef_sechiba

    ALLOCATE (snow_nobio_age(kjpindex,nnobio),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snow_nobio_age allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    snow_nobio_age(:,:) = undef_sechiba

    ALLOCATE (assim_param(kjpindex,nvm,npco2),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in assim_param allocation. We stop. We need kjpindex x nvm x npco2 words = ',&
            & kjpindex,' x ' ,nvm,' x ',npco2, ' = ',kjpindex*nvm*npco2
       STOP 'sechiba_init'
    END IF

    !! 1.4 Initialize 1D array allocation 

!pss:+
    ALLOCATE (fwet_out(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in fwet_out allocation. We stop. we need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    fwet_out(:) = undef_sechiba
!pss:-
!pss:+
    ALLOCATE (drunoff_tot(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in drunoff_tot allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
!pss:-

    ALLOCATE (vevapflo(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vevapflo allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    vevapflo(:)=zero

    ALLOCATE (vevapsno(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vevapsno allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (vevapnu(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vevapnu allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (t2mdiag(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in t2mdiag allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (totfrac_nobio(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in totfrac_nobio allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (floodout(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in floodout allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (runoff(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in runoff allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (drainage(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in drainage allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (returnflow(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in returnflow allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    returnflow(:) = zero

    ALLOCATE (reinfiltration(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in reinfiltration allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    reinfiltration(:) = zero

    ALLOCATE (irrigation(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in irrigation allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    irrigation(:) = zero

    ALLOCATE (z0(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in z0 allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (roughheight(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in roughheight allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (emis(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in emis allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (tot_melt(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tot_melt allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (valpha(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in valpha allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (vbeta(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vbeta allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (fusion(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in fusion allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (rau(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in rau allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (deadleaf_cover(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in deadleaf_cover allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (stempdiag(kjpindex, nbdl),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in stempdiag allocation. We stop. We need kjpindex*nbdl words = ',&
            & kjpindex*nbdl
       STOP 'sechiba_init'
    END IF

    ALLOCATE (co2_flux(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in co2_flux allocation. We stop. We need kjpindex*nvm words = ' ,kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    co2_flux(:,:)=zero

    ALLOCATE (shumdiag(kjpindex,nbdl),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in shumdiag allocation. We stop. We need kjpindex*nbdl words = ',&
            & kjpindex*nbdl
       STOP 'sechiba_init'
    END IF
    
    ALLOCATE (shumdiag_perma(kjpindex,nbdl),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in shumdiag_perma allocation. We stop. We need kjpindex*nbdl words = ',&
            & kjpindex*nbdl
       STOP 'sechiba_init'
    END IF

    ALLOCATE (litterhumdiag(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in litterhumdiag allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    ALLOCATE (ptnlev1(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
        WRITE (numout,*) ' error in ptnlev1 allocation. We stop. We need kjpindex words = ', kjpindex
        STOP 'sechiba_init'
    END IF

    ALLOCATE (k_litt(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in k_litt allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF

    !! 1.5 Initialize 2D array allocation
    ALLOCATE (vevapwet(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in vevapwet allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF
    vevapwet(:,:)=undef_sechiba

    ALLOCATE (transpir(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in transpir allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (transpot(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in transpot allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (qsintmax(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in qsintmax allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm 
       STOP 'sechiba_init'
    END IF

    ALLOCATE (rveget(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in rveget allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (rstruct(kjpindex,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in rstruct allocation. We stop. We need kjpindex x nvm words = ',&
            & kjpindex,' x ' ,nvm, ' = ',kjpindex*nvm
       STOP 'sechiba_init'
    END IF

    ALLOCATE (pgflux(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in pgflux allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    END IF
    pgflux(:)= 0.0

    ALLOCATE (grndflux(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in grndflux allocation. We stop. We need kjpindex words = ',&
                         & kjpindex, ' = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    grndflux(:) = 0.0

    ALLOCATE (pkappa_snow(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in pkappa_snow allocation. We stop. We need kjpindex x nsnow words = ',&
                       & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    pkappa_snow(:,:) = 0.0


    ALLOCATE (gthick(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gthick allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    gthick(:) = 0

    ALLOCATE (gtemp(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gtemp allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    gtemp(:) = 280.0

    ALLOCATE (soilflxresid(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soilflxresid allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    soilflxresid(:) = 0.0

    ALLOCATE (gpkappa(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in gpkappa allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    gpkappa(:) = 0.2

    ALLOCATE (snowrho(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowrho allocation. We stop. We need kjpindex x nsnow words = ',&
            & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    snowrho(:,:) = xrhosmin

    ALLOCATE (snowheat(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowheat allocation. We stop. We need kjpindex x nsnow  words = ',&
            & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    snowheat(:,:) =0.0

    ALLOCATE (snowliq(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowliq allocation. We stop. We need kjpindex x nsnow words = ',&
            & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    snowliq(:,:) = 0.0

    ALLOCATE (snowgrain(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowgrain allocation. We stop. We need kjpindex x nsnow words = ',&
            & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    snowgrain(:,:) = 0.0

    ALLOCATE (snowtemp(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowtemp allocation. We stop. We need kjpindex x nsnow words = ',&
            & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    snowtemp(:,:) = tp_00

    ALLOCATE (soiltemp(kjpindex,ngrnd),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soiltemp allocation. We stop. We need kjpindex x ngrnd words = ',&
            & kjpindex,' x ' ,ngrnd, ' = ',kjpindex*ngrnd
       STOP 'sechiba_init'
    ENDIF
    soiltemp(:,:) = tp_00

    ALLOCATE (snowdz(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowdepth allocation. We stop. We need kjpindex x nsnow words = ',&
            & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    snowdz(:,:)=0.0

    ALLOCATE (cgrnd_soil(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in cgrnd_soil allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    cgrnd_soil(:) = 0

    ALLOCATE (dgrnd_soil(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in dgrnd_soil allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    dgrnd_soil(:) = 0

    ALLOCATE (zdz1_soil(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in zdz1_soil allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    zdz1_soil(:) = 0

    ALLOCATE (zdz2_soil(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in zdz2_soil allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    zdz2_soil(:) = 0

    ALLOCATE (cgrnd_snow(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in cgrnd_snow allocation. We stop. We need kjpindex x nsnow words = ',&
            & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    cgrnd_snow(:,:) = 0

    ALLOCATE (dgrnd_snow(kjpindex,nsnow),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in dgrnd_snow allocation. We stop. We need kjpindex x nsnow words = ',&
            & kjpindex,' x ' ,nsnow, ' = ',kjpindex*nsnow
       STOP 'sechiba_init'
    ENDIF
    dgrnd_snow(:,:) = 0

    ALLOCATE (lambda_snow(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in lambda_snow allocation. We stop. We need kjpindex words = ',&
            & kjpindex, ' = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    lambda_snow(:) = 0
    ALLOCATE (snowflx(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowflx allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
!    snowflx(:) = 0.0

    ALLOCATE (snowcap(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in snowcap allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
!    snowcap(:) = 0.0

    ALLOCATE (temp_sol_add(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in temp_sol_add allocation. We stop. We need kjpindex words = ',kjpindex
       STOP 'sechiba_init'
    ENDIF
    temp_sol_add(:) = 0.0

    !allocate arrays needed for permafrost calculations
    ALLOCATE(tdeep(kjpindex,ndeep,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in tdeep allocation. We stop. We need kjpindex x ndeep x nvm words = ',&
            & kjpindex,' x ', ndeep, ' x ', nvm, ' = ', kjpindex*ndeep*nvm
       STOP 'sechiba init'
    END IF

    ALLOCATE(hsdeep(kjpindex,ndeep,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in hsdeep allocation. We stop. We need kjpindex x ndeep x nvm words = ',&
            & kjpindex,' x ', ndeep, ' x ', nvm, ' = ', kjpindex*ndeep*nvm
       STOP 'sechiba init'
    END IF
    tdeep(:,:,:) = 250.
    hsdeep(:,:,:) = 1.0

    ALLOCATE(heat_Zimov(kjpindex,ndeep,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in heat_Zimov allocation. We stop. We need kjpindex x ndeep x nvm words = ',&
            & kjpindex,' x ', ndeep, ' x ', nvm, ' = ', kjpindex*ndeep*nvm
       STOP 'sechiba init'
    END IF
    heat_Zimov(:,:,:) = zero

  ! 1d arrays (xy)
    ALLOCATE(sfluxCH4_deep(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in sfluxCH4_deep allocation. We stop. We need kjpindex = ',kjpindex
       STOP 'sechiba init'
    END IF
    ALLOCATE(sfluxCO2_deep(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in sfluxCO2_deep allocation. We stop. We need kjpindex = ',kjpindex
       STOP 'sechiba init'
    END IF
    ALLOCATE(thawed_humidity(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in thawed_humidity allocation. We stop. We need kjpindex = ',kjpindex
       STOP 'sechiba init'
    END IF
    ALLOCATE(depth_organic_soil(kjpindex),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in depth_organic_soil allocation. We stop. We need kjpindex = ',kjpindex
       STOP 'sechiba init'
    END IF
    ! 1d arrays (ndeep)
    ALLOCATE(zz_deep(ndeep),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in zz_deep allocation. We stop. We need ndeep = ',ndeep
       STOP 'sechiba init'
    END IF
    ALLOCATE(zz_coef_deep(ndeep),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in zz_coef_deep allocation. We stop. We need ndeep = ',ndeep
       STOP 'sechiba init'
    END IF
    ALLOCATE(soilc_total(kjpindex,ndeep,nvm),stat=ier)
    IF (ier.NE.0) THEN
       WRITE (numout,*) ' error in soilc_total allocation. We stop. We need kjpindex x ndeep x nvm words = ',&
            & kjpindex,' x ', ndeep, ' x ', nvm, ' = ', kjpindex*ndeep*nvm
       STOP 'sechiba init'
    END IF

    !! 1.6 Initialize indexing table for the vegetation fields. 
    ! In SECHIBA we work on reduced grids but to store in the full 3D filed vegetation variable 
    ! we need another index table : indexveg, indexsoil, indexnobio and indexgrnd
    DO ji = 1, kjpindex
       !
       DO jv = 1, nlai+1
          indexlai((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, nvm
          indexveg((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !      
       DO jv = 1, nstm
          indexsoil((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !      
       DO jv = 1, nnobio
          indexnobio((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, ngrnd
          indexgrnd((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, nsnow
          indexsnow((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO

       DO jv = 1, nbdl
          indexnbdl((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO

       DO jv = 1, nslm
          indexlayer((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, 2
          indexalb((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
    ENDDO
    ! define thermodynamic grid for permafrost diffusion calculations here
    CALL thermosoil_vert_axes(zz_deep, zz_coef_deep)
    !

!! 2. Read restart files to set initial variable values
    
    IF (ldrestart_read) THEN

       IF (long_print) WRITE (numout,*) ' we have to read a restart file for SECHIBA variables'
       ! Open restart file and read data
       ! ?? Where is the code to open the restart file?? !!?? I guess, the IF just writes out the comment.
       !!?? The single routines are called in sechiba_main. It would be more logical to put it in sechiba_main
       ! Read the default value that will be put into variable which are not in the restart file
       CALL ioget_expval(val_exp)

    ENDIF

!! 3. Swapping flags
    !??What is going with the flags (also next section??). Why is this done?? Where are the values specified??
    river_routing = control_in%river_routing
    hydrol_cwrr = control_in%hydrol_cwrr
    
!! 4. Run control: store flags in a common variable

    control%river_routing = control_in%river_routing
    control%hydrol_cwrr = control_in%hydrol_cwrr
    control%ok_co2 = control_in%ok_co2
    control%ok_sechiba = control_in%ok_sechiba
    control%ok_stomate = control_in%ok_stomate
    control%ok_dgvm = control_in%ok_dgvm
    control%do_land_use = control_in%do_land_use
    control%ok_pheno = control_in%ok_pheno
    control%stomate_watchout = control_in%stomate_watchout
    control%ok_inca = control_in%ok_inca                    
    control%ok_leafage = control_in%ok_leafage         
    control%ok_radcanopy = control_in%ok_radcanopy
    control%ok_multilayer = control_in%ok_multilayer
    control%ok_pulse_NOx = control_in%ok_pulse_NOx       
    control%ok_bbgfertil_NOx = control_in%ok_bbgfertil_NOx     
    control%ok_cropsfertil_NOx = control_in%ok_cropsfertil_NOx

    IF (long_print) WRITE (numout,*) ' sechiba_init done '

  END SUBROUTINE sechiba_init
  

!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_clear
!!
!>\BRIEF        Deallocate memory of sechiba's variables
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================ 

  SUBROUTINE sechiba_clear (forcing_name,cforcing_name)

    CHARACTER(LEN=100), INTENT(in)           :: forcing_name       !! Name of forcing file (unitless)
    CHARACTER(LEN=100), INTENT(in)           :: cforcing_name      !! Name of forcing file with carbon related variables (unitless)
!_ ================================================================================================================================
    
!! 1. Initialize first run

    l_first_sechiba=.TRUE.

!! 2. Deallocate dynamic variables of sechiba

    IF ( ALLOCATED (indexveg)) DEALLOCATE (indexveg)
    IF ( ALLOCATED (indexlai)) DEALLOCATE (indexlai)
    IF ( ALLOCATED (indexsoil)) DEALLOCATE (indexsoil)
    IF ( ALLOCATED (indexnobio)) DEALLOCATE (indexnobio)
    IF ( ALLOCATED (indexsnow)) DEALLOCATE (indexsnow)
    IF ( ALLOCATED (indexgrnd)) DEALLOCATE (indexgrnd)
    IF ( ALLOCATED (indexlayer)) DEALLOCATE (indexlayer)
    IF ( ALLOCATED (indexnbdl)) DEALLOCATE (indexnbdl)
    IF ( ALLOCATED (indexalb)) DEALLOCATE (indexalb)
    IF ( ALLOCATED (flood_res)) DEALLOCATE (flood_res)
    IF ( ALLOCATED (flood_frac)) DEALLOCATE (flood_frac)
    IF ( ALLOCATED (snow)) DEALLOCATE (snow)
    IF ( ALLOCATED (snow_age)) DEALLOCATE (snow_age)
    IF ( ALLOCATED (drysoil_frac)) DEALLOCATE (drysoil_frac)
    IF ( ALLOCATED (rsol)) DEALLOCATE (rsol)
    IF ( ALLOCATED (evap_bare_lim)) DEALLOCATE (evap_bare_lim)
    IF ( ALLOCATED (evapot)) DEALLOCATE (evapot)
    IF ( ALLOCATED (evapot_corr)) DEALLOCATE (evapot_corr)
    IF ( ALLOCATED (humrel)) DEALLOCATE (humrel)
    IF ( ALLOCATED (vegstress)) DEALLOCATE (vegstress)
    IF ( ALLOCATED (soiltile)) DEALLOCATE (soiltile)
    IF ( ALLOCATED (njsc)) DEALLOCATE (njsc)
    IF ( ALLOCATED (reinf_slope)) DEALLOCATE (reinf_slope)
    IF ( ALLOCATED (vbeta1)) DEALLOCATE (vbeta1)
    IF ( ALLOCATED (vbeta4)) DEALLOCATE (vbeta4)
    IF ( ALLOCATED (vbeta5)) DEALLOCATE (vbeta5)
    IF ( ALLOCATED (soilcap)) DEALLOCATE (soilcap)
    IF ( ALLOCATED (soilflx)) DEALLOCATE (soilflx)
    IF ( ALLOCATED (snowcap)) DEALLOCATE (snowcap)
    IF ( ALLOCATED (snowflx)) DEALLOCATE (snowflx)
    IF ( ALLOCATED (temp_sol)) DEALLOCATE (temp_sol)
    IF ( ALLOCATED (qsurf)) DEALLOCATE (qsurf)
    IF ( ALLOCATED (qsintveg)) DEALLOCATE (qsintveg)
    IF ( ALLOCATED (vbeta2))  DEALLOCATE (vbeta2)
    IF ( ALLOCATED (vbeta3)) DEALLOCATE (vbeta3)
    IF ( ALLOCATED (vbeta3pot)) DEALLOCATE (vbeta3pot)
    IF ( ALLOCATED (gsmean)) DEALLOCATE (gsmean)
    IF ( ALLOCATED (cimean)) DEALLOCATE (cimean)
    IF ( ALLOCATED (gpp)) DEALLOCATE (gpp)
    IF ( ALLOCATED (temp_growth)) DEALLOCATE (temp_growth) 
    IF ( ALLOCATED (veget)) DEALLOCATE (veget)
    IF ( ALLOCATED (veget_max)) DEALLOCATE (veget_max)
    IF ( ALLOCATED (lai)) DEALLOCATE (lai)
    IF ( ALLOCATED (frac_age)) DEALLOCATE (frac_age)
    IF ( ALLOCATED (height)) DEALLOCATE (height)
    IF ( ALLOCATED (roughheight)) DEALLOCATE (roughheight)
    IF ( ALLOCATED (frac_nobio)) DEALLOCATE (frac_nobio)
    IF ( ALLOCATED (snow_nobio)) DEALLOCATE (snow_nobio)
    IF ( ALLOCATED (snow_nobio_age)) DEALLOCATE (snow_nobio_age)
    IF ( ALLOCATED (assim_param)) DEALLOCATE (assim_param)
    IF ( ALLOCATED (vevapflo)) DEALLOCATE (vevapflo)
    IF ( ALLOCATED (vevapsno)) DEALLOCATE (vevapsno)
    IF ( ALLOCATED (vevapnu)) DEALLOCATE (vevapnu)
    IF ( ALLOCATED (t2mdiag)) DEALLOCATE (t2mdiag)
    IF ( ALLOCATED (totfrac_nobio)) DEALLOCATE (totfrac_nobio)
    IF ( ALLOCATED (floodout)) DEALLOCATE (floodout)
    IF ( ALLOCATED (runoff)) DEALLOCATE (runoff)
    IF ( ALLOCATED (drainage)) DEALLOCATE (drainage)
    IF ( ALLOCATED (reinfiltration)) DEALLOCATE (reinfiltration)
    IF ( ALLOCATED (irrigation)) DEALLOCATE (irrigation)
    IF ( ALLOCATED (tot_melt)) DEALLOCATE (tot_melt)
    IF ( ALLOCATED (valpha)) DEALLOCATE (valpha)
    IF ( ALLOCATED (vbeta)) DEALLOCATE (vbeta)
    IF ( ALLOCATED (fusion)) DEALLOCATE (fusion)
    IF ( ALLOCATED (rau)) DEALLOCATE (rau)
    IF ( ALLOCATED (deadleaf_cover)) DEALLOCATE (deadleaf_cover)
    IF ( ALLOCATED (stempdiag)) DEALLOCATE (stempdiag)
    IF ( ALLOCATED (co2_flux)) DEALLOCATE (co2_flux)
    IF ( ALLOCATED (shumdiag)) DEALLOCATE (shumdiag)
    IF ( ALLOCATED (shumdiag_perma)) DEALLOCATE (shumdiag_perma)
    IF ( ALLOCATED (litterhumdiag)) DEALLOCATE (litterhumdiag)
    IF ( ALLOCATED (ptnlev1)) DEALLOCATE (ptnlev1)
    IF ( ALLOCATED (k_litt)) DEALLOCATE (k_litt)
    IF ( ALLOCATED (vevapwet)) DEALLOCATE (vevapwet)
    IF ( ALLOCATED (transpir)) DEALLOCATE (transpir)
    IF ( ALLOCATED (transpot)) DEALLOCATE (transpot)
    IF ( ALLOCATED (qsintmax)) DEALLOCATE (qsintmax)
    IF ( ALLOCATED (rveget)) DEALLOCATE (rveget)
    IF ( ALLOCATED (rstruct)) DEALLOCATE (rstruct)
    IF ( ALLOCATED (snowrho)) DEALLOCATE (snowrho)
    IF ( ALLOCATED (snowgrain)) DEALLOCATE (snowgrain)
    IF ( ALLOCATED (snowtemp)) DEALLOCATE (snowtemp)
    IF ( ALLOCATED (soiltemp)) DEALLOCATE (soiltemp)
    IF ( ALLOCATED (snowdz)) DEALLOCATE (snowdz)
    IF ( ALLOCATED (snowliq)) DEALLOCATE (snowliq)
    IF ( ALLOCATED (snowheat)) DEALLOCATE (snowheat)
    IF ( ALLOCATED (grndflux)) DEALLOCATE (grndflux)
    IF ( ALLOCATED (gtemp)) DEALLOCATE (gtemp)
    IF ( ALLOCATED (soilflxresid)) DEALLOCATE (soilflxresid)
    IF ( ALLOCATED (gpkappa)) DEALLOCATE (gpkappa)
    IF ( ALLOCATED (gthick)) DEALLOCATE (gthick)
    IF ( ALLOCATED (pgflux)) DEALLOCATE (pgflux)
    IF ( ALLOCATED (pkappa_snow)) DEALLOCATE (pkappa_snow)

!pss:+
    IF ( ALLOCATED (fwet_out)) DEALLOCATE (fwet_out)
    IF ( ALLOCATED (drunoff_tot)) DEALLOCATE (drunoff_tot)
!pss:-

!! 3. Clear all allocated memory

    CALL pft_parameters_clear
    CALL slowproc_clear 
    CALL diffuco_clear 
    CALL enerbil_clear  
    IF ( hydrol_cwrr ) THEN
       CALL hydrol_clear 
    ELSE
       CALL hydrolc_clear  
    ENDIF
    CALL condveg_clear 
    CALL thermosoil_clear
    CALL routing_clear

  END SUBROUTINE sechiba_clear


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_var_init
!!
!>\BRIEF        Calculate air density as a function of air temperature and 
!! pressure for each terrestrial pixel.
!! 
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): air density (::rau, kg m^{-3}).
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE sechiba_var_init (kjpindex, rau, pb, temp_air) 

!! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                    :: kjpindex        !! Domain size - terrestrial pixels only (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: pb              !! Surface pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: temp_air        !! Air temperature (K)
    
!! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out) :: rau             !! Air density @tex $(kg m^{-3})$ @endtex

!! 0.3 Modified variables

!! 0.4 Local variables

    INTEGER(i_std)                                 :: ji              !! Indices (unitless)
!_ ================================================================================================================================
    
!! 1. Calculate intial air density (::rau)
   
    DO ji = 1,kjpindex
       rau(ji) = pa_par_hpa * pb(ji) / (cte_molr*temp_air(ji))
    END DO

    IF (long_print) WRITE (numout,*) ' sechiba_var_init done '

  END SUBROUTINE sechiba_var_init


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_end
!!
!>\BRIEF        Swap old for newly calculated soil temperature.
!! 
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): soil temperature (::temp_sol; K)
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!! ================================================================================================================================ 

  SUBROUTINE sechiba_end (kjpindex, dtradia, temp_sol_new, snowtemp, snowdz, &
                          temp_sol)
                         

!! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                       :: kjpindex           !! Domain size - terrestrial pixels only (unitless)
    REAL(r_std),INTENT (in)                           :: dtradia            !! Time step (seconds)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: temp_sol_new       !! New soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in) :: snowtemp           !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT(in) :: snowdz             !! Snow layer thickness (m)
    
    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)    :: temp_sol           !! Soil temperature (K)

    !! 0.3 Local variables
    INTEGER(i_std) :: ji

!_ ================================================================================================================================
    
!! 1. Swap temperature

    IF (ok_explicitsnow) THEN
       DO ji=1,kjpindex
          ! When the snow does not exist on the ground, the old surface temperature
          ! should be taken from the previous snow surface temperature
          IF (SUM(snowdz(ji,:)) .GT. 0.0) THEN
             temp_sol(ji) = temp_sol_new(ji)
          ELSE
             temp_sol(ji) = temp_sol_new(ji)
          END IF
       END DO

    ELSE
       temp_sol(:) = temp_sol_new(:)
    END IF

    IF (long_print) WRITE (numout,*) ' sechiba_end done '

  END SUBROUTINE sechiba_end

END MODULE sechiba
