MODULE Grassland_Management
  ! this module include grassland management from PaSim
  ! graze - cut - fertilisation 
  ! with auto management or user-defined management

  USE constantes_PaSim
  USE constantes
  USE fonctions_PaSim  
  USE Animaux
  USE applic_plant
  USE Fauche
  USE Fertilisation
  USE ioipsl
  USE ioipsl_para
!  USE parallel

  IMPLICIT NONE

  LOGICAL, SAVE :: first_call_grassland_manag = .TRUE. 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intakemax
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake_litter
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake_animal_litter
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: litter_avail_totDM
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: grazing_litter
  ! shoot dry matter afer cutting Kg/m2
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: wshtotcutinit
  ! lai after cutting (m**2 leaf/m**2)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: lcutinit     
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: devstage
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: faecesc
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: faecesn
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: urinen
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: urinec
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nel
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nanimaltot
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: tgrowth               
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wsh
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshtotinit
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshtot
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wr
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wrtot
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wanimal
  ! concentration totale en N (kg n/kg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ntot   
  ! concentration en C du substrat de la plante (kg C/kg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: c      
  ! concentration en N du substrat de la plante (kg N/kg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: n      
  ! n in structral mass kgN/kg
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fn     
  ! n concentration of apoplast (kgN/kg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: napo   
  ! n concentration of symplast (kgN/kg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nsym   
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wnapo
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wnsym
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wn
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nanimal
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tanimal
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: danimal             
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tcut  
  ! day of fertilisation (management) (d)           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert   
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Nliquidmanure    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Nslurry           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Nsolidmanure      
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: legume_fraction
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)    :: soil_fertility
  ! Threshold shoot dry matter, under which animals are moved out (kg/m2)
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: Animalwgrazingmin
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: AnimalkintakeM
  ! parameter for calculation of vegetation compartement selection by animals (-)
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: AnimalDiscremineQualite
  ! Valeurs associ�es � la croissance a�rienne entre 2 �v�nements de fertilisation (autogestion fauches)
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: controle_azote 
  ! Carbon flux from Organic fertilization to metabolic SOM pool (kg C m-2 day-1)
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fcOrganicFertmetabolicsum    
  ! Carbon flux from Organic fertilization to strcutural SOM pool (kg C m-2 day-1)
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fcOrganicFertstructsum       
  ! Nitrogen flux from Organic fertilization to strcutural SOM pool (kg N m-2 day-1)
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFertmetabolicsum    
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFertstructsum
  ! Nitrogen flux coming from slurry and liquid manure (k N.m-2)
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFerturinesum        
  ! Nitrogen deposition (kg N m-2 year-1)
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnatmsum                     
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: controle_azote_sum
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: nfertamm
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: nfertnit
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intakesum
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intakensum
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake_animal
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: intake_animalsum
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: PIYcow
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: PIMcow
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: BCSYcow
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: BCSMcow
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: PICcow
  ! Age of dairy primi cow
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: AGE_cow_P           
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: AGE_cow_M
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Autogestion_out
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: Forage_quantity
  REAL(r_std ), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tcut_modif
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: countschedule
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: mux
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: mugmean
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: sigx
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: sigy
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: gmeanslope
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: gzero
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: gcor      
  INTEGER (i_std)  , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: cuttingend
  LOGICAL   , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tcut_verif
  INTEGER(i_std)   , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: regcount
  INTEGER(i_std)                                       :: tcutmodel
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshcutinit      
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: gmean           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tgmean           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wc_frac              
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wgn             
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: tasum           
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: loss            
  ! perte en C lors de la fauche
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: lossc                 
  ! perte en N lors de la fauche
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: lossn                 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: tlossstart         
  INTEGER(i_std) , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: flag_fertilisation
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fertcount
  ! C/N dans le r�servoir stucturel  = 150
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: c2nratiostruct        
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertammtot
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertnittot
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertammtotyear
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertnittotyear
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertammtotprevyear
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nfertnittotprevyear
  ! metabolic C in slurry and manure (kg C/m**2/d)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fcOrganicFertmetabolic 
  ! structural C in slurry and manure (kg C/m**2/d)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fcOrganicFertstruct   
  ! urine N in slurry and manure (kg N/m**2/d)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFerturine     
  ! metabolic N in slurry and manure (kg N/m**2/d)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fnOrganicFertmetabolic            

  ! variables pour l'auto gestion de nicolas
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nsatur_somerror_temp
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nsatur_somerror
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tfert_modif
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nnonlimit_SOMerror
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: nnonlimit_SOMerrormax
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: controle_azote_sum_mem
  INTEGER(i_std) , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: n_auto
  INTEGER(i_std) , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: stoplimitant
  INTEGER(i_std) , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fertcount_start
  INTEGER(i_std) , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fertcount_current
  LOGICAL   , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: fertil_year
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: toto
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION(:,:)   :: wshtotsumprevyear
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_sr_ugb_C3
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_nb_ani_C3
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_grazed_frac_C3
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_import_yield_C3
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_wshtotsum_C3
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_sr_ugb_C4
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_nb_ani_C4
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_grazed_frac_C4
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_import_yield_C4
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: tmp_wshtotsum_C4

  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION(:,:)   :: DM_cutyearly
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION(:,:)   :: C_cutyearly
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: YIELD_RETURN
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)   :: sr_ugb_init

  INTEGER(i_std)                  , SAVE                 :: cut_year
  INTEGER(i_std)                  , SAVE                 :: compt_fert
  INTEGER(i_std)                  , SAVE                 :: min_fert
  INTEGER(i_std)                  , SAVE                 :: fert_max
  INTEGER(i_std)                  , SAVE                 :: i_compt
  REAL(r_std)                     , SAVE                 :: deltat             ! = 0.02
  ! how may years that managemeng applied
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:)        :: nb_year_management 
  ! couter number of years of simulation     
  INTEGER(i_std)                  , SAVE                 :: count_year            
  INTEGER(i_std)                  , SAVE                 :: year_count1
  ! couter number of years of simulation
  INTEGER(i_std)                  , SAVE                 :: year_count2
  ! couter number of years of simulation

  CHARACTER(len=500), ALLOCATABLE,SAVE,DIMENSION (:)     :: file_management
  CHARACTER(len=500)              , SAVE                 :: file_param_init
  CHARACTER(len=500)              , SAVE                 :: file_import_yield

  INTEGER(i_std)                  , SAVE                 :: Type_animal
  INTEGER(i_std)                  , SAVE                 :: mcut_C3
  INTEGER(i_std)                  , SAVE                 :: mauto_C3
  INTEGER(i_std)                  , SAVE                 :: mcut_C4
  INTEGER(i_std)                  , SAVE                 :: mauto_C4
  ! yearly total azote by fertilization
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: apport_azote
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: trampling

  ! new variables for get map of management
  INTEGER(i_std)                  , SAVE                 :: f_management_map 
  CHARACTER(len=500)              , SAVE                 :: management_map
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:)        :: management_intensity
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:)        :: management_start
  CHARACTER(len=500)              , SAVE                 :: fertility_map

  INTEGER(i_std)                  , SAVE                 :: f_deposition_map
  CHARACTER(len=500)              , SAVE                 :: deposition_map
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:)        :: deposition_start
  INTEGER(i_std)                  , SAVE                 :: f_grazing_map
  CHARACTER(len=500)              , SAVE                 :: grazing_map
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: ndeposition
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)   :: N_fert_total
  REAL(r_std)                     , SAVE                 :: N_effect
CONTAINS

  SUBROUTINE Main_Grassland_Management(&
     npts         , &
     dt           , &
     tjulian      , &
     ta           , &
     t2m_min_daily, &
     t2m_14       , &
     tsoil        , &
     snow         , &
     biomass      , &
     bm_to_litter , &
     litter       , &
     litter_avail , &
     litter_not_avail , &
     new_day      , &
     new_year     , &
     flag_cutting , &
     when_growthinit_cut , &
     lai          , &         
     sla_calc     , &
     leaf_age     , &
     leaf_frac    , &
     wshtotsum    , &
     sr_ugb       , &
     compt_ugb    , &
     nb_ani       , &
     grazed_frac  , &
     import_yield , &
     N_limfert)

    INTEGER(i_std)                                , INTENT(in)   :: npts   
    REAL(r_std)                             , INTENT(in)   :: dt          
    INTEGER(i_std)                             , INTENT(in)   :: tjulian 
    ! julien day
    REAL(r_std), DIMENSION(npts)            , INTENT(in)   :: ta      
    ! air temperature 
    REAL(r_std), DIMENSION(npts)            , INTENT(in)   ::  t2m_min_daily
    ! daily minimum temperature
    REAL(r_std), DIMENSION(npts)            , INTENT(in)   ::  t2m_14   
    ! 14 days mean temperature
    REAL(r_std), DIMENSION(npts)            , INTENT(in)   :: tsoil     
    ! soil surface t (k)
    REAL(r_std), DIMENSION(npts)            , INTENT(in)   :: snow       
    ! snow mass (mm)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass       
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: bm_to_litter 
    ! conv of biomass to litter (gC/(m**2/agri ground)) / day
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm) , INTENT(inout):: litter_not_avail
    LOGICAL                                , INTENT(in)   :: new_day   
    ! flag indicate new day
    LOGICAL                                , INTENT(in)   :: new_year   
    ! flag indicate new year
    INTEGER(i_std)   , DIMENSION(npts,nvm)            , INTENT(out)  :: flag_cutting
    ! flag indicate if there is a cut
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                 :: when_growthinit_cut
    ! how many days ago was the beginning of the last cut
    REAL(r_std), DIMENSION(npts,nvm)       , INTENT(out)  :: lai        
    ! leaf area index OF AN INDIVIDUAL PLANT
    REAL(r_std), DIMENSION(npts,nvm)       , INTENT(in)  :: sla_calc
    REAL(r_std), DIMENSION(npts,nvm,nleafages)       , INTENT(inout)  :: leaf_frac
    REAL(r_std), DIMENSION(npts,nvm,nleafages)       , INTENT(inout)  :: leaf_age
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                 :: wshtotsum
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  sr_ugb
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  compt_ugb
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  nb_ani
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  grazed_frac
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  import_yield
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)  ::  N_limfert

    LOGICAL :: l_error = .FALSE.
    INTEGER(i_std) :: ier, i, j, k,h, m
    REAL(r_std), DIMENSION(npts)        :: xtmp_npts
    REAL(r_std), DIMENSION(npts,ngmean) :: xtmp_npts_3d
    REAL(r_std), DIMENSION(npts,nvm)        :: regcount_real
    REAL(r_std), DIMENSION(npts,nvm)        :: fertcount_real
    INTEGER(i_std) :: fertcount_next
    REAL(r_std) :: intakemax_t
    REAL(r_std) :: wanimal_t
    REAL(r_std), DIMENSION(ncut)        ::wshtotcutinit_t
    REAL(r_std), DIMENSION(npts,nvm)        :: lm_before
    REAL(r_std), DIMENSION(npts,nvm)        :: lm_after
    REAL(r_std), DIMENSION(npts,nvm)        :: bm_cut
!    REAL(r_std), DIMENSION(npts,nvm)        :: N_fert_total
    REAL(r_std) :: fertility_legume_t
    ! initialisations
    ! 

    init_grassland : IF (first_call_grassland_manag) THEN

      first_call_grassland_manag = .FALSE. 

      ALLOCATE (intake                (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (intakemax             (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (intake_litter         (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (intake_animal_litter  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (grazing_litter        (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (litter_avail_totDM    (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wshtotcutinit         (npts,nvm,ncut)     , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (lcutinit              (npts,nvm,ncut)     , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (devstage              (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)   
      ALLOCATE (faecesc               (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (faecesn               (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (urinen                (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (urinec                (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (nel                   (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (nanimaltot            (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tgrowth               (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wsh                   (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wshtot                (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wshtotinit            (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wr                    (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wrtot                 (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wanimal               (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (ntot                  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (c                     (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (n                     (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (fn                    (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (napo                  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (nsym                  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wnapo                 (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wnsym                 (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wn                    (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (nanimal               (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tanimal               (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (danimal               (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tcut                  (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tfert                 (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (Nliquidmanure         (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (Nslurry               (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (Nsolidmanure          (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (legume_fraction       (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (soil_fertility        (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (Animalwgrazingmin     (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (AnimalkintakeM        (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (AnimalDiscremineQualite (npts,nvm)        , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (controle_azote        (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (fcOrganicFertmetabolicsum (npts,nvm)      , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fcOrganicFertstructsum (npts,nvm)         , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fnOrganicFertmetabolicsum (npts,nvm)      , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fnOrganicFertstructsum (npts,nvm)         , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fnOrganicFerturinesum (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fnatmsum              (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (controle_azote_sum    (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (nfertamm              (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0) 
      ALLOCATE (nfertnit              (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0) 
      ALLOCATE (intakesum             (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (intakensum            (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (intake_animal         (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (intake_animalsum      (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (PIYcow                (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (PIMcow                (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (BCSYcow               (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (BCSMcow               (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (PICcow                (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (AGE_cow_P             (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (AGE_cow_M             (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (Autogestion_out       (npts,nvm,n_out)    , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (Forage_quantity       (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tcut_modif            (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (countschedule         (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (mux                   (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (mugmean               (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (sigx                  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (sigy                  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (gmeanslope            (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (gzero                 (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (gcor                  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (cuttingend            (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tcut_verif            (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tfert_verif           (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tfert_verif2          (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tfert_verif3          (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (regcount              (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wshcutinit            (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (gmean                 (npts,nvm,ngmean)   , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tgmean                (npts,nvm,ngmean)   , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wc_frac               (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wgn                   (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tasum                 (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (loss                  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (lossc                 (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (lossn                 (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tlossstart            (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (flag_fertilisation    (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (fertcount             (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0) 
      ALLOCATE (c2nratiostruct        (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (nfertammtot           (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (nfertnittot           (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (nfertammtotyear       (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)   
      ALLOCATE (nfertnittotyear       (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)   
      ALLOCATE (nfertammtotprevyear   (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)   
      ALLOCATE (nfertnittotprevyear   (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)   
      ALLOCATE (fcOrganicFertmetabolic (npts,nvm)         , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fcOrganicFertstruct   (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fnOrganicFerturine    (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fnOrganicFertstruct   (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (fnOrganicFertmetabolic (npts,nvm)         , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (nsatur_somerror_temp  (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (nsatur_somerror       (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tfert_modif           (npts,nvm,nstocking), stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (nnonlimit_SOMerror    (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (nnonlimit_SOMerrormax (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (controle_azote_sum_mem (npts,nvm)         , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (n_auto                (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (stoplimitant          (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (fertcount_start       (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (fertcount_current     (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wshtotsumprev         (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (fertil_year           (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (toto                  (npts)              , stat=ier); l_error=l_error .OR. (ier .NE. 0)  
      ALLOCATE (apport_azote          (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (trampling             (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (wshtotsumprevyear     (npts,nvm)          , stat=ier); l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE(file_management        (nvm)               , stat=ier);l_error =l_error .OR. (ier .NE. 0)
      ALLOCATE(nb_year_management     (nvm)               , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_sr_ugb_C3         (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_nb_ani_C3         (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_grazed_frac_C3    (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_import_yield_C3   (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_wshtotsum_C3      (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_sr_ugb_C4         (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_nb_ani_C4         (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_grazed_frac_C4    (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_import_yield_C4   (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (tmp_wshtotsum_C4      (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (DM_cutyearly          (npts,nvm)          , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (C_cutyearly           (npts,nvm)          , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (YIELD_RETURN          (npts,nvm)          , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (sr_ugb_init           (npts)              , stat=ier);l_error=l_error .OR. (ier .NE. 0)

      !  new variables for get map of management
      ALLOCATE (management_intensity  (nvm)               , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (management_start      (nvm)               , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (deposition_start      (nvm)               , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (N_fert_total          (npts,nvm)          , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      ALLOCATE (ndeposition           (npts,nvm)          , stat=ier);l_error=l_error .OR. (ier .NE. 0)
      IF (l_error) THEN

          STOP 'erreur allocation memory in Grassland_management'

      END IF

      ! concentrations : mean value
      c(:,:)                  = 0.0365122     !  4.22e-02
      n(:,:)                  = 0.00732556    !  8.17e-03
      napo(:,:)               = 0.000542054   !  6.39e-04
      nsym(:,:)               = 0.0108071     !  6.15e-03
      fn(:,:)                 = 0.0316223     !  4.15e-02   ! 2.64e-02
      ntot(:,:)               = 0.03471895    !  2.89e-02  

      ! set flags and variables need to read in Pasim
 
      ! saturant N supply
      f_saturant = 0
      CALL getin('F_SATURANT',f_saturant)
      ! N fertilization without limitation
      f_nonlimitant = 0
      CALL getin('F_NONLIMITANT',f_nonlimitant)
      ! f_autogestion = 1-5
      ! 1: auto cut for PFT m_auto
      ! 2: auto graze for PFT m_auto
      ! 3: auto cut and graze for PFT m_cut and m_grazed with increasing sr_ugb
      ! 4: auto cut and graze for PFT m_cut and m_grazed with constant sr_ugb
      ! 5: auto graze for PFT m_grazed with grazing litter during winter for LGM period
      f_autogestion = 0
      CALL getin('F_AUTOGESTION',f_autogestion)
      ! whether animal is fed by extra feedstuffs
      f_complementation = 0
      CALL getin('F_COMPLEMENTATION',f_complementation)
      ! whether apply fertilizer
      f_fertilization = 1         
      CALL getin('F_FERTILIZATION',f_fertilization)
      ! number of management year cycled
      nb_year_management(:) = 0
      CALL getin_p('NB_YEAR_MANAGEMENT',nb_year_management)
      WRITE(numout,*) 'NB_YEAR_MANAGEMENT',nb_year_management
      ! f_postauto = 0-5
      ! 1: after f_autogestion=2 with varied sr_ugb and nb_ani
      ! 2: after f_postauto=1 with varied sr_ugb and nb_ani
      f_postauto = 0
      CALL getin('F_POSTAUTO',f_postauto)
      ! the maximum impact to vcmax due to N fertilization
      ! N_effect =0.0 - 1.0
      N_effect=0.6
      CALL getin('N_EFFECT',N_effect)
      IF (N_effect .LT. 0.0 .OR. N_effect .GT. 1.0) THEN
        N_effect =0.6
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! READ INITIAL CONDITION FILE FOR OLD/NEW ANIMAL MODULE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      file_param_init='/home/orchidee_ns/lhli/Modele_ORCHIDEE/Management/param_init.txt'

      CALL getin('FILE_PARAM_INIT',file_param_init)
      WRITE (numout,*) 'FILE_PARAM_INIT',file_param_init
      OPEN(unit=61, file = file_param_init)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) (wshtotcutinit_t(h), h=1,ncut)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) tcutmodel
      READ(61, *, iostat = ier) intakemax_t

      READ(61, *, iostat = ier) wanimal_t
      READ(61, *, iostat = ier) Type_animal
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      READ(61, *, iostat = ier) toto(:)
      READ(61, *, iostat = ier) toto(:)

      intakemax(:,:)=intakemax_t
      wanimal(:,:)=wanimal_t
      DO h=1,ncut
        wshtotcutinit(:,:,h)=wshtotcutinit_t(h)
      END DO
      CLOSE (61)
      WRITE(numout,*) 'type',Type_Animal
      ! initialisation the variables lied on animals cattle or sheep ?
      ! Type_Animal = 1,2,3 : Dairy cows, suckler cows, cows in old module 
      ! Type_Animal = 4,5 : Dairy heifers, suckler heifers
      ! Type_Animal = 6 : sheep in old module
      IF (Type_Animal==3)  THEN ! old module
      !090810 AIG changement du seuil de sortie des animaux Animalwgrazingmin trop faible
      ! changement de AnimalkintakeM pour garder que l'ingere garde la meme pente
      ! en fonction de la biomasse disponible et pour eviter un artefact de calcul
      ! Animalwgrazingmin = 0.03 ! Threshold shoot dry matter, under which animals are moved out for old module (kg.m-2) 
  
        Animalwgrazingmin(:,:)        = 0.03  ! N. Vuichard
        !AnimalkintakeM           = 0.18 ! AI Graux
        AnimalkintakeM(:,:)           = 0.1   ! N. Vuichard
        AnimalDiscremineQualite(:,:)  = 2    ! AI Graux  
      ELSEIF (Type_Animal .EQ. 6)THEN ! Sheep
        Animalwgrazingmin(:,:)        = 0.015
        AnimalkintakeM(:,:)           = 0.045
        AnimalDiscremineQualite(:,:)  = 3
      ELSE !new module
        !Animalwgrazingmin        = 0.11 ! AI Graux ! unsued in the new module
        !AnimalkintakeM           = 0.18 ! AI Graux ! unsued in the new module
        AnimalDiscremineQualite(:,:)  = 2    ! AI Graux  
      ENDIF ! Type_Animal

      ! initialisations
      intake(:,:)                = 0.0
      intake_litter(:,:)         = 0.0
      intake_animal_litter(:,:)  = 0.0
      grazing_litter(:,:)        = 2
      litter_avail_totDM(:,:)    = 0.0
      devstage(:,:)              = 0.0
      faecesc(:,:)               = 0.0
      faecesn(:,:)               = 0.0
      urinen(:,:)                = 0.0
      urinec(:,:)                = 0.0
      nel(:,:)                   = 0.0
      nanimaltot(:,:)            = 0.0
      grazingcstruct(:,:)        = 0.0
      grazingnstruct(:,:)        = 0.0
      tgrowth(:,:)               = 0.0
      wshtot(:,:) = (biomass(:,:,ileaf,icarbon) + biomass(:,:,isapabove,icarbon) + &
                    & biomass(:,:,ifruit,icarbon))/(1000*CtoDM) ! Unit: kgDM/m2
      wsh(:,:) = wshtot(:,:) / (1.0 + (mc /12.0)*c(:,:) + (mn /14.0)*n(:,:) )
      wshtotinit(:,:)            = wshtot(:,:)         
      wrtot(:,:) = (biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon))/ &
                   & (1000*CtoDM)   ! Unit: kg/m2
      wr(:,:) = wrtot(:,:) / (1.0 + (mc /12.0)*c(:,:) + (mn /14.0)*n(:,:) )
      wnapo(:,:)                 = 0.0
      wnsym(:,:)                 = 0.0
      mux(:,:)                   = 0.0
      mugmean(:,:)               = 0.0
      sigx(:,:)                  = 0.0
      sigy(:,:)                  = 0.0
      gmeanslope(:,:)            = 0.0
      gzero(:,:)                 = 0.0
      gcor(:,:)                  = 0.0
      countschedule(:,:)         = 0
      cuttingend(:,:)            = 0
      regcount(:,:)              = 1
      gmean(:,:,:)               = 0.0
      tgmean(:,:,:)              = 0.0
      wc_frac(:,:)               = 0.0
      wgn(:,:)                   = 0.0
      tasum(:,:)                 = 0.0
      loss(:,:)                  = 0.0
      lossc(:,:)                 = 0.0
      lossn(:,:)                 = 0.0
      tlossstart(:,:)            = 0.0
      wshcutinit(:,:)            = 0.0
      deltat                     = dt
      fertcount(:,:)             = 0.0
      c2nratiostruct(:,:)        = 150.0
      nfertammtot(:,:)           = 0.0
      nfertnittot(:,:)           = 0.0
      nfertammtotyear(:,:)       = 0.0
      nfertnittotyear(:,:)       = 0.0
      nfertammtotprevyear(:,:)   = 0.0
      nfertnittotprevyear(:,:)   = 0.0
      fcOrganicFertmetabolic(:,:)      = 0.0
      fcOrganicFertstruct(:,:)         = 0.0
      fnOrganicFertmetabolic(:,:)      = 0.0
      fnOrganicFertstruct(:,:)         = 0.0
      fnOrganicFerturine(:,:)          = 0.0
      flag_fertilisation(:,:)    = 0
      fertil_year(:,:)           = .TRUE.
      tcut_verif(:,:,:)          = .FALSE. 
      tfert_verif(:,:,:)         = .FALSE. 
      tfert_verif2(:,:,:)        = .FALSE.
      tfert_verif3(:,:,:)        = .FALSE.

      nsatur_somerror_temp(:,:)          = 0.0
      nsatur_somerror(:,:)               = 0.0
      stoplimitant(:,:)                  = 0
      fertcount_start(:,:)               = 0
      fertcount_current(:,:)             = 0
      nnonlimit_SOMerror(:,:)            = 0.0
      nnonlimit_SOMerrormax(:,:)         = 0.5

      controle_azote_sum_mem(:,:)        = 0.0
      n_auto(:,:)                        = 4
      flag_fertilisation(:,:)            = 0

      YIELD_RETURN(:,:) = 0.0
      sr_ugb_init(:) = 0.0

      ! Define PFT that used for optimization, cutting, and grazing
      DO j=2,nvm
        IF (is_grassland_cut(j).AND.(.NOT.is_grassland_grazed(j)) .AND. &
           (.NOT. is_c4(j)) .AND. (.NOT.is_tree(j)))THEN
          mcut_C3=j
        END IF
        IF ( is_grassland_manag(j) .AND.(.NOT. is_grassland_cut(j)).AND. &
          (.NOT.is_grassland_grazed(j)).AND. (.NOT. is_c4(j)) .AND. (.NOT.is_tree(j)))THEN
          mauto_C3=j
        END IF
        IF ( is_grassland_manag(j) .AND.(is_grassland_grazed(j)).AND. &
          (.NOT.is_grassland_cut(j)) .AND. (.NOT. is_c4(j)) .AND. (.NOT.is_tree(j)))THEN
          mgraze_C3=j
        END IF
        IF (is_grassland_cut(j).AND.(.NOT.is_grassland_grazed(j)) .AND. &
           (is_c4(j)) .AND. (.NOT.is_tree(j)))THEN
          mcut_C4=j
        END IF
        IF ( is_grassland_manag(j) .AND.(.NOT. is_grassland_cut(j)).AND. &
          (.NOT.is_grassland_grazed(j)).AND. (is_c4(j)) .AND. (.NOT.is_tree(j)))THEN
          mauto_C4=j
        END IF
        IF ( is_grassland_manag(j) .AND.(is_grassland_grazed(j)).AND. &
          (.NOT.is_grassland_cut(j)) .AND. (is_c4(j)) .AND. (.NOT.is_tree(j)))THEN
          mgraze_C4=j
        END IF
      END DO ! nvm
      WRITE(numout,*) 'PFT_M',mauto_C3,mcut_C3,mgraze_C3,mauto_C4,mcut_C4,mgraze_C4 
      ! avoid PFT = 0
      IF (mauto_C4 .EQ. 0) THEN 
        mauto_C4=1
      ENDIF
      IF (mcut_C4 .EQ. 0) THEN
        mcut_C4=1
      ENDIF
      IF (mgraze_C4 .EQ. 0) THEN
        mgraze_C4=1
      ENDIF
      IF (mauto_C3 .EQ. 0) THEN
        mauto_C3=1
      ENDIF
      IF (mcut_C3 .EQ. 0) THEN
        mcut_C3=1
      ENDIF
      IF (mgraze_C3 .EQ. 0) THEN
        mgraze_C3=1
      ENDIF
      WRITE(numout,*) 'PFT_M2',mauto_C3,mcut_C3,mgraze_C3,mauto_C4,mcut_C4,mgraze_C4
  
      ! Initialization of management related parameters for each management option
      IF ((f_postauto .EQ. 0) .AND. (f_autogestion .LE. 2))THEN
        sr_ugb         = 1e-5
        compt_ugb      = 0.0
        nb_ani         = 1e-10
        grazed_frac         = 0.50
      ELSE IF (f_postauto .EQ. 5) THEN
        sr_ugb         = 1e-5
        compt_ugb      = 0.0
        nb_ani         = 5e-10
        grazed_frac         = 0.50
      ELSE IF ((f_autogestion .EQ. 4) .OR. (f_autogestion .EQ. 3) .OR. &
              & (f_autogestion .EQ. 5)) THEN
        tmp_sr_ugb_C3(:)=sr_ugb(:,mgraze_C3)
        tmp_sr_ugb_C4(:)=sr_ugb(:,mgraze_C4)
        sr_ugb         = 1e-6
        sr_ugb(:,mgraze_C3)      = tmp_sr_ugb_C3(:)
        sr_ugb(:,mgraze_C4)      = tmp_sr_ugb_C4(:)
        compt_ugb      = 0.0
        grazed_frac         = 0.50
        nb_ani         = 1e-10
      ELSE IF (f_postauto .EQ. 1) THEN
        tmp_sr_ugb_C3(:)=sr_ugb(:,mauto_C3)
        tmp_sr_ugb_C4(:)=sr_ugb(:,mauto_C4)
        sr_ugb         = 1e-5
        sr_ugb(:,mgraze_C3)      = tmp_sr_ugb_C3(:)
        sr_ugb(:,mgraze_C4)      = tmp_sr_ugb_C4(:)
        compt_ugb      = 0.0
        tmp_nb_ani_C3(:)=nb_ani(:,mauto_C3)
        tmp_nb_ani_C4(:)=nb_ani(:,mauto_C4)
        nb_ani         = 1e-10
        nb_ani(:,mgraze_C3)         = tmp_nb_ani_C3(:)
        nb_ani(:,mgraze_C4)         = tmp_nb_ani_C4(:)
        tmp_grazed_frac_C3(:)=grazed_frac(:,mauto_C3)
        tmp_grazed_frac_C4(:)=grazed_frac(:,mauto_C4)
        grazed_frac         = 0.50
        grazed_frac(:,mgraze_C3)         = tmp_grazed_frac_C3(:)
        grazed_frac(:,mgraze_C4)         = tmp_grazed_frac_C4(:)
       WHERE (sr_ugb(:,mgraze_C3) .GT. 0.0)
        grazed_frac(:,mgraze_C3)  = nb_ani(:,mgraze_C3)/sr_ugb(:,mgraze_C3) 
       ELSEWHERE
        grazed_frac(:,mgraze_C3)  = tmp_grazed_frac_C3(:)
       ENDWHERE
       WHERE (sr_ugb(:,mgraze_C4) .GT. 0.0)
        grazed_frac(:,mgraze_C4)  = nb_ani(:,mgraze_C4)/sr_ugb(:,mgraze_C4)  
       ELSEWHERE
        grazed_frac(:,mgraze_C4)  = tmp_grazed_frac_C4(:)
       ENDWHERE
  
      ELSE IF ((f_postauto .EQ. 2) .OR. (f_postauto .EQ. 3) .OR. &
              & (f_postauto .EQ. 4)) THEN
        tmp_sr_ugb_C3(:)=sr_ugb(:,mgraze_C3)
        tmp_sr_ugb_C4(:)=sr_ugb(:,mgraze_C4)
        sr_ugb         = 1e-5
        sr_ugb(:,mgraze_C3)      = tmp_sr_ugb_C3(:)
        sr_ugb(:,mgraze_C4)      = tmp_sr_ugb_C4(:)
        sr_ugb_init(:) = tmp_sr_ugb_C3(:)
        compt_ugb      = 0.0
        tmp_nb_ani_C3(:)=nb_ani(:,mgraze_C3)
        tmp_nb_ani_C4(:)=nb_ani(:,mgraze_C4)
        nb_ani         = 1e-10
        nb_ani(:,mgraze_C3)         = tmp_nb_ani_C3(:)
        nb_ani(:,mgraze_C4)         = tmp_nb_ani_C4(:)
        tmp_grazed_frac_C3(:)=grazed_frac(:,mgraze_C3)
        tmp_grazed_frac_C4(:)=grazed_frac(:,mgraze_C4)
        grazed_frac         = 0.50
        grazed_frac(:,mgraze_C3)         = tmp_grazed_frac_C3(:)
        grazed_frac(:,mgraze_C4)         = tmp_grazed_frac_C4(:)
       WHERE (sr_ugb(:,mgraze_C3) .GT. 0.0)
        grazed_frac(:,mgraze_C3)  = nb_ani(:,mgraze_C3)/sr_ugb(:,mgraze_C3)
       ELSEWHERE
        grazed_frac(:,mgraze_C3)  = tmp_grazed_frac_C3(:)
       ENDWHERE
       WHERE (sr_ugb(:,mgraze_C4) .GT. 0.0)
        grazed_frac(:,mgraze_C4)  = nb_ani(:,mgraze_C4)/sr_ugb(:,mgraze_C4)
       ELSEWHERE
        grazed_frac(:,mgraze_C4)  = tmp_grazed_frac_C4(:)
       ENDWHERE
  
      ENDIF ! f_autogestion or f_postauto
  
      IF ((f_autogestion .NE. 2 ) .AND. (f_postauto .EQ. 0))  THEN
        wshtotsum (:,:) = 0.0
        import_yield (:,:) = 0.0
      ELSE IF (f_postauto .EQ. 5) THEN
        wshtotsum (:,:) = 0.0
        import_yield (:,:) = 0.0
      END IF
  
      IF (f_autogestion .EQ. 2 ) THEN
        CALL getin_p('NB_CUT_YEAR',cut_year)
        DO j=2,nvm
          IF (is_grassland_manag(j) .AND. (.NOT.is_grassland_cut(j)) .AND. &
            (.NOT.is_grassland_grazed(j)))THEN
            WHERE ( wshtotsum(:,j) .GE. 0.0)
              import_yield(:,j) = wshtotsum(:,j)/cut_year
            ENDWHERE
          END IF
        END DO
      END IF
  
      IF (f_postauto .EQ. 1 ) THEN
        tmp_import_yield_C3(:) = import_yield(:,mauto_C3)
        tmp_import_yield_C4(:) = import_yield(:,mauto_C4)
        import_yield = 0.0
        import_yield (:,mgraze_C3) = tmp_import_yield_C3(:) 
        import_yield (:,mgraze_C4) = tmp_import_yield_C4(:)
      END IF
  
      IF ((f_postauto .EQ. 2 ) .OR. (f_postauto .EQ. 3) &
         .OR. (f_postauto .EQ. 4)) THEN
        tmp_import_yield_C3(:) = import_yield(:,mgraze_C3)
        tmp_import_yield_C4(:) = import_yield(:,mgraze_C4)
        import_yield = 0.0
        import_yield (:,mgraze_C3) = tmp_import_yield_C3(:)
        import_yield (:,mgraze_C4) = tmp_import_yield_C4(:)
        tmp_wshtotsum_C3(:) = wshtotsum(:,mcut_C3)
        tmp_wshtotsum_C4(:) = wshtotsum(:,mcut_C4)
        wshtotsumprevyear(:,:) = 0.0
        wshtotsumprevyear(:,mcut_C3) = tmp_wshtotsum_C3(:)
        wshtotsumprevyear(:,mcut_C4) = tmp_wshtotsum_C4(:)
      END IF
  
      DM_cutyearly(:,:)=0.0
      C_cutyearly(:,:) =0.0
      IF (f_postauto .LT. 2 .OR. f_postauto .EQ. 5) THEN
      wshtotsumprevyear(:,:) = 0.0
      ENDIF
      wshtotsumprev   (:,:) = 0.0
      controle_azote(:,:,:)       = 0.0
      controle_azote_sum(:,:)        = 0.0
      trampling(:,:)              = 0.0
      count_year            = 1
      year_count1 = 0
      year_count2 = 0
      tcut(:,:,:) = 500.0
      tfert(:,:,:) = 500.0
      nfertamm(:,:,:) = 0.0
      nfertnit(:,:,:) = 0.0
      nanimal(:,:,:) = 0.0
      tanimal(:,:,:) = 500.0
      danimal(:,:,:) = 0.0
      nliquidmanure(:,:,:) = 0.0
      nslurry(:,:,:) = 0.0
      nsolidmanure(:,:,:) = 0.0
      legume_fraction(:,:) =0.0
      soil_fertility(:,:) = 1.0
      ndeposition(:,:) = 0.0
      PIYcow(:,:,:) = 0.0
      PIMcow(:,:,:) = 0.0
      BCSYcow(:,:,:) = 0.0
      BCSMcow(:,:,:) = 0.0
      PICcow(:,:,:) = 0.0
      AGE_cow_P(:,:,:) = 36.0
      AGE_cow_M(:,:,:) = 54.0
      Forage_quantity(:,:,:) = 0
  
      IF (blabla_pasim) PRINT *, 'PASIM : end memory allocation'
  
      ! get_map of 1 spatial .nc file or 0 old txt/dat file  
      CALL getin ('F_MANAGEMENT_MAP',f_management_map)
      WRITE(numout,*)  'F_MANAGEMENT_MAP',f_management_map
      CALL getin ('F_DEPOSITION_MAP',f_deposition_map)
      WRITE(numout,*)  'F_DEPOSITION_MAP',f_deposition_map  
      CALL getin ('F_GRAZING_MAP',f_grazing_map)
      WRITE(numout,*)  'F_GRAZING_MAP',f_grazing_map
   
      IF (f_management_map .EQ. 1) THEN    !! Arsene 19-12-2014. before: IF (f_management_map) THEN
        management_map='/ccc/work/cont003/dsm/p529chan/data/eur_management_interpolated.nc'
        CALL getin('MANAGEMENT_MAP',management_map)
        WRITE(numout,*) 'MANAGEMENT_MAP',management_map
        ! managemeng_intensity 1->Low 2->Medium 3->High
        CALL getin_p('MANAGEMENT_INTENSITY',management_intensity)
        WRITE(numout,*) 'MANAGEMENT_INTENSITY',management_intensity
        CALL getin_p('MANAGEMENT_START',management_start)
        WRITE(numout,*) 'MANAGEMENT_START',management_start
        fertility_map='/ccc/work/cont003/dsm/p529chan/data/eur_fertility.nc'
        CALL getin('FERTILITY_MAP',fertility_map)
        WRITE(numout,*) 'FERTILITY_MAP',fertility_map
  
        deposition_map='/ccc/work/cont003/dsm/p529chan/data/eur_Ndeposition_NCAR.nc'
        CALL getin('DEPOSITION_MAP',deposition_map)
        WRITE(numout,*) 'DEPOSITION_MAP',deposition_map  
        CALL getin_p('DEPOSITION_START',deposition_start)
        WRITE(numout,*) 'DEPOSITION_START',deposition_start
        grazing_map='/ccc/scratch/cont003/dsm/p529chan/glbdata/glb_sr_ugb_1961_2010_adjusted.nc'
        CALL getin('GRAZING_MAP',grazing_map)
        WRITE(numout,*) 'GRAZING_MAP',grazing_map
  
        ! read management map    
        CALL reading_map_manag(&
               npts, count_year, nb_year_management,& 
               management_intensity,&
               management_start,&
               tcut, tfert, nfertamm, nfertnit,&
               nanimal, tanimal, danimal,&
               nliquidmanure, nslurry, nsolidmanure,&
               legume_fraction,soil_fertility,&
               deposition_start,ndeposition,sr_ugb)
        ! calculate effect of N fertilizer to vcmax
        CALL calc_N_limfert(&
               npts,nfertamm, nfertnit,&
               nliquidmanure, nslurry, nsolidmanure,&
               legume_fraction,soil_fertility,ndeposition,&
               N_fert_total,N_limfert)
  
  
      ELSE
        ! re-initial management variables
        tcut(:,:,:) = 500.0
        tfert(:,:,:) = 500.0
        nfertamm(:,:,:) = 0.0
        nfertnit(:,:,:) = 0.0
        nanimal(:,:,:) = 0.0
        tanimal(:,:,:) = 500.0
        danimal(:,:,:) = 0.0
        nliquidmanure(:,:,:) = 0.0
        nslurry(:,:,:) = 0.0
        nsolidmanure(:,:,:) = 0.0
        ndeposition(:,:) = 0.0
  !!! delete FIRE_MANAGEMENT READ: not used in LGM
  !      CALL getin_p('FILE_MANAGEMENT',file_management)
  !      WRITE(numout,*)  'FILE_MANAGEMENT',file_management
  !      IF (blabla_pasim) PRINT *, 'PASIM : reading management conditions'
  !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      !!!!!!!!! READ NEW MANAGEMENT TXT DAT FILE JCADD
  !      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      CALL reading_new_animal(&
  !           npts           , &
  !           nb_year_management , &
  !           tcutmodel      , &
  !           tcut           , &
  !           tfert          , &
  !           nfertamm       , &
  !           nfertnit       , &
  !           nanimal        , &
  !           tanimal        , &
  !           danimal        , &
  !           nliquidmanure  , &
  !           nslurry        , &
  !           nsolidmanure   , &
  !           PIYcow         , &
  !           PIMcow         , &
  !           BCSYcow        , &
  !           BCSMcow        , &
  !           PICcow         , &
  !           AGE_cow_P      , &
  !           AGE_cow_M      , &
  !           Forage_quantity)
  !      CALL getin('SOIL_FERTILITY',fertility_legume_t)
  !      soil_fertility(:,:)=fertility_legume_t
  !      CALL getin('LEGUME_FRACTION',fertility_legume_t)
  !      legume_fraction(:,:)=fertility_legume_t
  !
  !      CALL calc_N_limfert(&
  !             npts,nfertamm, nfertnit,&
  !             nliquidmanure, nslurry, nsolidmanure,&
  !             legume_fraction,soil_fertility,ndeposition,&
  !             N_fert_total,N_limfert)
  
      ENDIF ! f_management_map
  
      DO k=1,nstocking
        WHERE (tfert(:,:,k) .NE. 500) 
          apport_azote(:,:) = apport_azote(:,:) + nfertamm(:,:,k) + nfertnit(:,:,k)    
        END WHERE  
      END DO
        !************************************************
        !************************************************
        !
        ! AUTO GESTION INITIALISATION VARIABLES
        ! modifs Nico 20/07/2004
        !
        !************************************************
        !************************************************
        ! MODIF INN
      IF (f_nonlimitant .EQ. 1) THEN
        IF (f_autogestion .NE. 2) THEN
          WHERE (tcut(:,:,1) .EQ. 500.0)
            stoplimitant(:,:) = 1
          END WHERE
        ENDIF
        DO j=2,nvm
          DO i=1,npts
            IF (tfert(i,j,1) .EQ. 500.0) THEN
              stoplimitant(i,j) = 1
            ELSE
              compt_fert = 1
              min_fert   = 1
              DO WHILE (tfert(i,j,compt_fert) .NE. 500.0)
                 print *, compt_fert, min_fert
                 print *, controle_azote(i,j,compt_fert)
                 print *, controle_azote(i,j,min_fert)
                IF (controle_azote(i,j,compt_fert) .GT. controle_azote(i,j,min_fert)) THEN
                  min_fert = compt_fert
                ENDIF
                  compt_fert = compt_fert + 1
              END DO
              fert_max = compt_fert - 1
              IF ((min_fert - 1) .EQ. 0) THEN
                fertcount_start(i,j) = fert_max
              ELSE
                fertcount_start(i,j) = min_fert - 1
              ENDIF
                i_compt = min_fert + 1
              DO WHILE ( tfert(i,j,i_compt) .NE. 500.0 )
                controle_azote(i,j,i_compt) = controle_azote(i,j,i_compt - 1)+&
                  controle_azote(i,j,i_compt)
                i_compt = i_compt + 1
              END DO
              IF ( min_fert .NE. 1. ) THEN
                controle_azote(i,j,1) = controle_azote(i,j,1) + controle_azote(i,j,fert_max)
                i_compt = 2
                DO WHILE (i_compt .NE. min_fert)
                  controle_azote(i,j,i_compt) = controle_azote(i,j,i_compt-1)+&
                    controle_azote(i,j,i_compt)
                  i_compt = i_compt + 1
                END DO
              ENDIF
            ENDIF
          END DO ! i
        END DO !j
          fertcount_current(:,:) = fertcount_start(:,:)
      ENDIF
        ! fin initialisation auto gestion nicolas
    END IF init_grassland

    ! update the root/shoot dry matter variables
    wshtot(:,:) = (biomass(:,:,ileaf,icarbon) + biomass(:,:,isapabove,icarbon) + &
                 & biomass(:,:,ifruit,icarbon))/(1000*CtoDM) ! Unit: kgDM/m2
    wsh(:,:) = wshtot(:,:) / (1.0 + (mc /12.0)*c(:,:) + (mn /14.0)*n(:,:) )
    wrtot(:,:) = (biomass(:,:,iroot,icarbon) + biomass(:,:,isapbelow,icarbon))/ &
                 & (1000*CtoDM)   ! Unit: kg/m2
    wr(:,:) = wrtot(:,:) / (1.0 + (mc /12.0)*c(:,:) + (mn /14.0)*n(:,:) )

    n_day : IF (new_day) THEN

      ! GMEAN
      ! Taux de croissance moyen de la repousse
      h  = 1

      DO WHILE (h  .LT. ngmean)
        gmean(:,:,h ) = gmean(:,:,h +1)
        h  = h  + 1
      END DO

      DO j=2,nvm  
        DO i=1,npts
          IF ((tgrowth(i,j) .GT. 0.0) .AND. (devstage(i,j) .GE. 2.0)) THEN
            gmean(i,j,ngmean) = MAX (0.0, (wshtot(i,j) - wshtotcutinit(i,j,regcount(i,j)))/tgrowth(i,j))
          ELSEIF ((tgrowth(i,j) .GT. 0.0) .AND. (devstage(i,j) .GT. 0.0) .AND. &
            & (regcount(i,j) .GT. 1)) THEN
            gmean(i,j,ngmean) = MAX (0.0, (wshtot(i,j) - wshtotcutinit(i,j,regcount(i,j)))/tgrowth(i,j))
          ELSEIF ((tgrowth(i,j) .GT. 0.0) .AND. (devstage(i,j) .GT. 0.0) .AND. &
            & (regcount(i,j) .EQ. 1)) THEN
            gmean(i,j,ngmean) = MAX (0.0, (wshtot(i,j)  - wshtotinit(i,j))/tgrowth(i,j))
          ELSE
            gmean(i,j,ngmean) = 0.0
          ENDIF
        END DO
      ENDDO
 
      h = 1
      DO WHILE (h .LE. ngmean) 
        tgmean(:,:,h) = h
        h = h + 1
      END DO
      
    END IF n_day


    n_year : IF (new_year) THEN

      tcut_verif(:,:,:)         = .FALSE. 
      fertil_year(:,:)          = .TRUE. 
      tasum(:,:)                = 0.0
      regcount(:,:)             = 1
      nfertammtotprevyear(:,:)  = nfertammtot 
      nfertnittotprevyear(:,:)  = nfertnittot 
      fertcount(:,:)            = 0
      nfertammtotyear(:,:)      = 0.0
      nfertnittotyear(:,:)      = 0.0
      fnatmsum(:,:)             = 0.0
      tfert_verif(:,:,:)        = .FALSE.
      tfert_verif2(:,:,:)       = .FALSE.
      tfert_verif3(:,:,:)       = .FALSE.
      fcOrganicFertmetabolicsum(:,:) = 0.0
      fcOrganicFertstructsum(:,:)    = 0.0
      fnOrganicFertmetabolicsum(:,:) = 0.0
      fnOrganicFertstructsum(:,:)    = 0.0
      fnOrganicFerturinesum(:,:)     = 0.0
      devstage(:,:)             = 0.0
      fertcount(:,:)            = 0
      tgrowth (:,:)             = 0.0
      tfert_modif(:,:,:)        = 500.0

      IF (f_saturant .EQ. 1) THEN
         nfertamm(:,:,:)  = 0.025
         nfertnit(:,:,:)  = 0.025
         nsatur_somerror(:,:)      = 0.0
         nsatur_somerror_temp(:,:) = 0.0
      END IF
      IF ((f_postauto .EQ. 1) .OR. (f_autogestion  .EQ. 4) .OR. &
         !!!! JCMODIF 290714 for postaut = 5
         & (f_postauto .GE. 2)) THEN
        import_yield(:,mgraze_C3) = wshtotsum(:,mcut_C3)-wshtotsumprevyear(:,mcut_C3)
        wshtotsumprevyear(:,mcut_C3) = wshtotsum(:,mcut_C3)
        import_yield(:,mgraze_C4) = wshtotsum(:,mcut_C4)-wshtotsumprevyear(:,mcut_C4)
        wshtotsumprevyear(:,mcut_C4) = wshtotsum(:,mcut_C4)
      END IF

      DM_cutyearly(:,:)= wshtotsum(:,:)-wshtotsumprevyear(:,:)
      C_cutyearly(:,:) = DM_cutyearly(:,:) * 1000 * CtoDM
      wshtotsumprevyear(:,:) = wshtotsum(:,:)
      wshtotsumprev(:,:)          = 0.0
      c(:,:)                  = 0.0365122     !  4.22e-02
      n(:,:)                  = 0.00732556    !  8.17e-03
      napo(:,:)               = 0.000542054   !  6.39e-04
      nsym(:,:)               = 0.0108071     !  6.15e-03
      fn(:,:)                 = 0.0316223     !  4.15e-02   ! 2.64e-02
      ntot(:,:)               = 0.03471895    !  2.89e-02  

      count_year = count_year + 1
      IF (count_year .LT. 30) THEN
        year_count1 = count_year-1
        year_count2 = 0
      ELSEIF (count_year .GE. 30) THEN
        year_count1 = 29
        year_count2 = count_year - 29
      ELSE 
        year_count1 = 29
        year_count2 = 21
      ENDIF
      ! get_map of spatial .nc file or old txt/dat file  
      IF (f_management_map .EQ. 1) THEN    !! Arsene 19-12-2014. before: IF (f_management_map) THEN
        ! re-initial management variables
        tcut(:,:,:) = 500.0
        tfert(:,:,:) = 500.0
        nfertamm(:,:,:) = 0.0
        nfertnit(:,:,:) = 0.0
        nanimal(:,:,:) = 0.0
        tanimal(:,:,:) = 500.0
        danimal(:,:,:) = 0.0
        nliquidmanure(:,:,:) = 0.0
        nslurry(:,:,:) = 0.0
        nsolidmanure(:,:,:) = 0.0
        ndeposition(:,:) = 0.0
        CALL reading_map_manag(&  
               npts, count_year, nb_year_management,&
               management_intensity,&
               management_start,& 
               tcut, tfert, nfertamm, nfertnit,&
               nanimal, tanimal, danimal,&
               nliquidmanure, nslurry, nsolidmanure,&
               legume_fraction,soil_fertility,&
               deposition_start,ndeposition,sr_ugb)
  
        CALL calc_N_limfert(&
               npts,nfertamm, nfertnit,&
               nliquidmanure, nslurry, nsolidmanure,&
               legume_fraction,soil_fertility,ndeposition,&
               N_fert_total,N_limfert)
  
      ELSE
        IF (ANY(nb_year_management(:) .GT. 1)) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!! READ NEW MANAGEMENT FILE JCADD
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! delete FIRE_MANAGEMENT READ: not used in LGM
    !        CALL reading_new_animal(&
    !           npts           , &
    !           nb_year_management , &
    !           tcutmodel      , &
    !           tcut           , &
    !           tfert          , &
    !           nfertamm       , &
    !           nfertnit       , &
    !           nanimal        , &
    !           tanimal        , &
    !           danimal        , &
    !           nliquidmanure  , &
    !           nslurry        , &
    !           nsolidmanure   , &
    !           PIYcow         , &
    !           PIMcow         , &
    !           BCSYcow        , &
    !           BCSMcow        , &
    !           PICcow         , &
    !           AGE_cow_P      , &
    !           AGE_cow_M      , &
    !           Forage_quantity)
          ! re-initial management variables
          tcut(:,:,:) = 500.0
          tfert(:,:,:) = 500.0
          nfertamm(:,:,:) = 0.0
          nfertnit(:,:,:) = 0.0
          nanimal(:,:,:) = 0.0
          tanimal(:,:,:) = 500.0
          danimal(:,:,:) = 0.0
          nliquidmanure(:,:,:) = 0.0
          nslurry(:,:,:) = 0.0
          nsolidmanure(:,:,:) = 0.0
          ndeposition(:,:) = 0.0
          CALL calc_N_limfert(&
                 npts,nfertamm, nfertnit,&
                 nliquidmanure, nslurry, nsolidmanure,&
                 legume_fraction,soil_fertility,ndeposition,&
                 N_fert_total,N_limfert)
    
        END IF ! nb_year_management

        DO k=1,nstocking

          WHERE (tfert(:,:,k) .NE. 500) 
                  
            apport_azote(:,:) = apport_azote(:,:)  + nfertamm(:,:,k) + nfertnit(:,:,k)

          END WHERE

        END DO

      END IF ! f_management_map

    END IF n_year

    !------
    ! Fertilisation from PaSim 2011
    !-------
    !****** RUN USERS OR RUN SATURANT *****
    users_or_saturant_fert : IF ((tcutmodel .EQ. 0) .AND. (f_saturant .EQ. 0)) THEN

      ! flag_fertilisation : flag for spatialization of cutting

      DO k=1,nfert
        flag_fertilisation(:,:) = 0

        IF (ANY(tfert_verif(:,:,k) .EQV. .FALSE. )) THEN
          WHERE ((tjulian .GE. tfert(:,:,k)) .AND. (tjulian .LE. tfert(:,:,k)+0.9) .AND. &
            & (tfert_verif(:,:,k) .EQV. .FALSE.))
                tfert_verif(:,:,k) = .TRUE.
                flag_fertilisation(:,:) = 1
          END WHERE
        END IF

        IF (ANY(flag_fertilisation(:,:) .EQ. 1)) THEN
            CALL fertilisation_spa(&
               npts               , &
               flag_fertilisation , &
               fertcount_start    , &
               tjulian            , &
               tfert              , &
               nfertnittotprevyear, &
               nfertammtotprevyear, &
               nfertnit, nfertamm , &
               fertcount       , &
               nfertammtot     , &
               nfertnittot     , &
               nfertammtotyear    , &
               nfertnittotyear    , &
               controle_azote_sum , &
               controle_azote_sum_mem)

        END IF

      END DO
      !*****************************************
      ! MODIFS NICO AUTO MANAGEMENT DE PASIM
      !*****************************************

    ELSE IF (f_saturant .EQ. 1) THEN   !***** RUN SATURANT *******

      flag_fertilisation(:,:) = 0
      DO j=2 ,nvm
        DO i=1,npts

          IF (( tjulian .GE. tfert(i,j,fertcount(i,j) + 1)) .AND. &
             ( tjulian .LT. (tfert(i,j,fertcount(i,j) + 2) - 1))) THEN
             !JCmodif 110523  with problem
             !above means tjulian between two tfert 
             !undaily(i) uptake n daily always=0
             !thetas volumetric water content in soil layer h
             !thetasfc water field capacity
             !!!!! For we did not consider undaily , there will be no point need to fert???                  
             ! IF ((undaily(i) .GT. 0.0) .AND. (thetas(i,1) .LE. thetasfc(i,1))) THEN
                  flag_fertilisation(i,j) = 1
          ELSE
              flag_fertilisation(i,j) = 2
          END IF

        END DO
      END DO

      IF (ANY(flag_fertilisation(:,:) .EQ. 1)) THEN
          CALL fertilisation_spa(&
               npts                , &
               flag_fertilisation  , &
               fertcount_start     , &
               tjulian             , &
               tfert               , &
               nfertnittotprevyear , &
               nfertammtotprevyear , &
               nfertnit            , &
               nfertamm            , &
               fertcount           , &
               nfertammtot         , &
               nfertnittot         , &
               nfertammtotyear     , &
               nfertnittotyear     , &
               controle_azote_sum  , &
               controle_azote_sum_mem)
        DO j=2, nvm
          DO i=1,npts

            IF (flag_fertilisation(i,j) .EQ. 1) THEN
              IF (controle_azote_sum(i,j) .GT. 0.) THEN

                 nsatur_somerror_temp(i,j) = &
                   ABS(controle_azote(i,j,fertcount(i,j)) - controle_azote_sum(i,j)) / &
                   controle_azote_sum(i,j)

              ENDIF

            IF (nsatur_somerror_temp(i,j) .GT. nsatur_somerror(i,j)) THEN

              nsatur_somerror(i,j) = nsatur_somerror_temp(i,j)

            ENDIF
              controle_azote(i,j,fertcount(i,j) ) = controle_azote_sum(i,j)

              tfert_modif(i,j,fertcount(i,j) )    = tjulian
            END IF

          END DO
        END DO
      END IF ! flag_fertilisation(:,:) .EQ. 1

      IF (ANY(flag_fertilisation(:,:) .EQ. 2)) THEN

        WHERE (flag_fertilisation(:,:) .NE. 2)

          flag_fertilisation(:,:) = 0

        END WHERE
        DO j = 2, nvm
          DO i=1,npts
            IF ((tjulian .GE. (tfert(i,j,fertcount(i,j)+2)-1))  .AND. &
              (tjulian .LE. (tfert(i,j,fertcount(i,j)+2)-0.1)) .AND. &
              (.NOT.(tfert_verif2(i,j,fertcount(i,j)+2)) ) .AND. &
              (flag_fertilisation(i,j) .EQ. 2) ) THEN

              flag_fertilisation(i,j) = 1
              tfert_verif2(i,j,fertcount(i,j)+2) = .TRUE.

              nfertamm(i,j,fertcount(i,j) + 1) = 0.
              nfertnit(i,j,fertcount(i,j) + 1) = 0.

              tfert_modif(i,j,fertcount(i,j) + 1) = 500.0
            END IF
          END DO
        END DO

        IF (ANY(flag_fertilisation(:,:) .EQ. 1)) THEN

          CALL fertilisation_spa(&
                   npts                  , &
                   flag_fertilisation    , &
                   fertcount_start       , &
                   tjulian               , &
                   tfert                 , &
                   nfertnittotprevyear   , &
                   nfertammtotprevyear   , &
                   nfertnit              , &
                   nfertamm              , &
                   fertcount          , &
                   nfertammtot        , &
                   nfertnittot        , &
                   nfertammtotyear       , &
                   nfertnittotyear       , &
                   controle_azote_sum    , &
                   controle_azote_sum_mem)

        END IF

      END IF ! flag_fertilisation(:,:) .EQ. 2

    END IF users_or_saturant_fert


    !***** RUN NONLIMITANT *****
    ! recherche des erreurs pour l'�quilibre
    ! recherche de stoplimitant (fin du run)
    run_nonlimitant : IF ((f_nonlimitant .EQ. 1) .AND. (ANY(stoplimitant(:,:) .EQ. 0))) THEN   ! any ?

      DO j=2,nvm
        DO i=1,npts

          ! on recherche le dernier temps de fertilisation
          IF (tfert(i,j,fertcount_current(i,j) + 1) .EQ. 500) THEN
              fertcount_next = 1
          ELSE
              fertcount_next = fertcount_current(i,j) + 1
          ENDIF


          ! si tjulian correspond au prochain temps de fertilisation
          IF ((tjulian .GE. tfert(i,j,fertcount_next)) .AND. &
             (tjulian .LE. tfert(i,j,fertcount_next)+0.9) .AND. &
             (tfert_verif2(i,j,fertcount_next) .EQV. .FALSE.)) THEN

              tfert_verif2(i,j,fertcount_next) = .TRUE.

              ! calcul de somerror
              IF(controle_azote(i,j,fertcount_next) .GT. 0.) THEN
                  nnonlimit_SOMerror(i,j) = &
                     (controle_azote(i,j,fertcount_next) - controle_azote_sum_mem(i,j))/ &
                     controle_azote(i,j,fertcount_next)
              ELSE
                  nnonlimit_SOMerror(i,j) = 0.
              ENDIF
              ! on regarde si on ne d�passe pas l'erreur max voulue
              ! puis on r�ajuste cette erreur max suivant dans quel cas
              ! nous sommes
              IF (nnonlimit_SOMerror(i,j) .GT. nnonlimit_SOMerrormax(i,j)) THEN

                  nfertamm(i,j,fertcount_current(i,j)) = nfertamm(i,j,fertcount_current(i,j)) + 0.00125
                  nfertnit(i,j,fertcount_current(i,j)) = nfertnit(i,j,fertcount_current(i,j)) + 0.00125

                  PRINT *, '!!! apport en azote !!! pour fertcount_current = ', fertcount_current(i,j) &
                              ,' nfertamm= ',nfertamm(i,j,fertcount_current(i,j))

              ELSE

                  fertcount_current(i,j) = fertcount_current(i,j) + 1

                  IF(tfert(i,j,fertcount_current(i,j)) .EQ. 500.) THEN
                      fertcount_current(i,j) = 1
                  ENDIF

                  IF (fertcount_current(i,j) .EQ. fertcount_start(i,j)) THEN

                      nnonlimit_SOMerrormax(i,j) = nnonlimit_SOMerrormax(i,j) - n_auto(i,j)*0.05

                      n_auto(i,j) = n_auto(i,j) - 1.

                      IF ( nnonlimit_SOMerrormax(i,j) .LE. 0.) THEN
                         nnonlimit_SOMerrormax(i,j)=0.025
                      ENDIF

                      IF (n_auto(i,j) .LT. 0.) THEN

                          stoplimitant(i,j) = 1
                          print *,'*********************************'
                          print *,'stoplimitant =1 '
                          print *,'********************************'
                      ENDIF ! n_auto
                  ENDIF ! fertcount_current
              ENDIF ! nnonlimit
          ENDIF ! tjulian
        END DO ! npts
      END DO ! nvm
    END IF run_nonlimitant

    print *,'before pas temps'
    CALL fertilisation_pas_temps(&
       npts                           , &
       fertcount                      , &
       dt                             , &
       tjulian                        , &
       deltat                         , &
       tfert                          , &
       Nliquidmanure                  , &
       Nslurry                        , &
       Nsolidmanure                   , &
       fcOrganicFertmetabolic         , &
       fcOrganicFertstruct            , &
       fnOrganicFerturine             , &
       fnOrganicFertmetabolic         , &
       c2nratiostruct)
    print *,'before animal'

    DO j=2,nvm
      CALL Euler_funct(npts, dt, MAX(0.0,(ta - 278.15)), tasum(:,j))
    END DO

    ! calculate variables that not included in ORCHIDEE
    ! liste : 
    ! * devstage              
    ! * tgrowth               
    CALL Main_appl_pre_animal(&
       npts                  , &
       dt                    , &
       tjulian               , &
       ta                    , &
       tsoil                 , &
       new_day               , &
       new_year              , &
       regcount              , &
       tcut                  , &
       devstage              , &
       tgrowth              )

    !JCADD prepare litter_avail for grazing litter
    ! kg DM/m^2
    litter_avail_totDM(:,:) = (litter_avail(:,istructural,:) + &
       & litter_avail(:,imetabolic,:)) / (1000. * CtoDM) 
    !ENDJCADD

    IF ((Type_animal.EQ.3).OR.(Type_animal.EQ.6)) THEN ! old animal module
      CALL Animaux_main(&
       npts, dt, devstage, wsh, intakemax, &
       snow, wshtot, Animalwgrazingmin, &
       AnimalkintakeM, nel, wanimal, nanimaltot, &
       ntot, intake, urinen, faecesn, urinec, faecesc, &
       tgrowth, new_year, new_day, &
       nanimal, tanimal, danimal, &
       tcutmodel, tjulian, import_yield, &
       intakesum, intakensum, fn, c, n, leaf_frac, &
       intake_animal, intake_animalsum, &
       biomass, trampling, sr_ugb, &
       compt_ugb, nb_ani, grazed_frac,AnimalDiscremineQualite, &
       YIELD_RETURN,sr_ugb_init,year_count1,year_count2, & 
       grazing_litter, litter_avail_totDM, &
       intake_animal_litter, intake_litter)

    ELSE ! new animal module

      CALL Animaux_main_dynamic(&
        npts, dt, devstage                  , &
        intakemax, snow, wshtot, wsh        , &
        nel, nanimaltot                     , &
        intake                              , &
        import_yield                        , &
        new_year, new_day                   , &
        nanimal, tanimal, danimal           , &
        PIYcow, PIMcow, BCSYcow             , &
        BCSMcow, PICcow, AGE_cow_P          , &
        AGE_cow_M, tcutmodel, tjulian       , &
        intakesum                           , &
        intakensum, fn, ntot, c, n,leaf_frac, &
        intake_animal, intake_animalsum     , &
        t2m_min_daily, type_animal          , &
        ta, intakemax, Autogestion_out      , &
        Forage_quantity,t2m_14              , &
        intake_tolerance                    , &
        q_max_complement                    , &
        biomass, urinen, faecesn, urinec,faecesc, &
        file_param_init,trampling,sr_ugb    , &
        compt_ugb, nb_ani,grazed_frac,AnimalDiscremineQualite, &
        grazing_litter) 

    ENDIF ! Type_Animal
    print *,'after animal'

    !!!!!!!!!
    ! CUTTING
    !!!!!!!!!
    ! Cutting Management: auto_fauche and user_fauche

    flag_cutting(:,:) = 0

    ! user defined cut 
    user_fauche : IF ((f_autogestion .EQ. 0) .AND. (f_postauto .EQ. 0)) THEN

      flag_cutting(:,:) = 0
      when_growthinit_cut(:,:) = when_growthinit_cut(:,:) + dt

      lm_before(:,:)= biomass(:,:,ileaf,icarbon)
      DO k=1,nstocking

        flag_cutting(:,:) = 0

        DO j=2,nvm
          IF (is_grassland_manag(j) )THEN

            IF (ANY(tcut_verif(:,j,k) .EQV. .FALSE.)) THEN
              WHERE ((tjulian .GE. tcut(:,j,k)) .AND. (tjulian .LE. tcut(:,j,k)+0.9) .AND. &
                (tcut_verif(:,j,k) .EQV. .FALSE.))

                tcut_verif(:,j,k) = .TRUE. 
                flag_cutting(:,j) = 1
                when_growthinit_cut(:,j) = 0.0                  
              END WHERE
            END IF
          END IF
        END DO              
        IF (ANY (flag_cutting(:,:) .EQ. 1)) THEN 

          IF (blabla_pasim)  PRINT *, 'cutting users', tjulian

          CALL cutting_spa(&
                   npts              , &
                   tjulian           , &
                   flag_cutting      , &
                   wshtotcutinit     , &
                   lcutinit          , &
                   wsh               , &
                   wshtot            , &
                   wr                , &
                   c                 , &
                   n                 , &
                   napo              , &
                   nsym              , &
                   fn                , &
                   tjulian           , &
                   nel               , &
                   biomass           , &
                   devstage          , &
                   regcount          , &
                   wshcutinit        , &
                   gmean             , &
                   wc_frac                , &
                   wnapo             , &
                   wnsym             , &
                   wgn               , &
                   tasum             , &
                   tgrowth           , &
                   loss              , &
                   lossc             , &
                   lossn             , &
                   tlossstart        , &
                   lai               , &
                   tcut              , &
                   tcut_modif        , &
                   wshtotsum         , &
                   controle_azote_sum)

          WHERE ((wsh + wr .GT. 0.0).AND. (flag_cutting .EQ. 1)) 
              
            c = wc_frac / (wsh + wr)
            n = (wnapo + wnsym) / (wsh + wr) 
            fn = wgn / (wr + wsh)
            napo = wnapo / (wsh + wr)
            nsym = wnsym / (wsh + wr)

          END WHERE
               
          WHERE (wshtot + wrtot .GT. 0.0)

            ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)

          END WHERE
        END IF

      END DO !nstocking

    END IF user_fauche

    ! auto cut
    n_day_autofauche :  IF (new_day) THEN
      DO  j=2,nvm
        CALL linreg_pasim (&
           npts          , &
           ngmean        , &
           tgmean(:,j,:)        , &
           gmean(:,j,:)         , &
           ngmean        , &
           misval        , &
           mux(:,j)           , &
           mugmean(:,j)       , &
           sigx(:,j)          , &
           sigy(:,j)          , &
           gmeanslope(:,j)    , &
           gzero(:,j)         , &
           gcor(:,j))
      END DO
      countschedule(:,:)  = 1

      auto_fauche : IF (f_autogestion .EQ. 1) THEN ! for optimalize sr_ugb and nb_ani

        flag_cutting(:,:) = 0
        when_growthinit_cut(:,:) = when_growthinit_cut(:,:) + dt
        DO j=2,nvm
          IF (is_grassland_manag(j) .AND. (.NOT. is_grassland_cut(j)) .AND. &
            (.NOT. is_grassland_grazed(j)))THEN

        ! FIRST test for automanagement > 45 days 
          WHERE((nanimal(:,j,1) .EQ. 0.0) .AND. (cuttingend(:,j) .EQ. 0) .AND. &
            (countschedule(:,j) .EQ. 1) .AND. (((tgrowth(:,j) .GE. tgrowthmin) .AND. &
            (gmean(:,j,ngmean) .GT. 0.0) .AND.(lai(:,j) .GE. 2.5)  .AND. &
            (devstage(:,j) .GT. devstagemin ) .AND. &
            (gmeanslope(:,j) .LT. gmeansloperel * mugmean(:,j)))))

            flag_cutting(:,j)  = 1
            countschedule(:,j) = countschedule(:,j)  + 1             

          END WHERE 
          END IF
        END DO !nvm

        ! If there is at least one point concerned (flag_cutting = 1)
        IF (ANY (flag_cutting(:,:) .EQ. 1)) THEN 
          IF (blabla_pasim) PRINT *, 'FAUCHE AVEC METHODE NV ', tjulian
            ! There will be one fertilization the day after cutting
            ! A COURT-CIRCUITER si couplage autogestion ferti avec INN
            ! AIG 06/10/2009

            IF (f_fertilization.NE.1) THEN
              DO j=2,nvm   
                DO i=1,npts
                  IF (flag_cutting(i,j) .EQ. 1) THEN

                    tfert(i,j,regcount(i,j) + 1) = tjulian + 1
                    print*, 'FERTILISATION AVEC METHODE NV', tjulian

                  END IF
                END DO
              END DO  
            END IF

            CALL cutting_spa(&
               npts              , &
               tjulian           , &
               flag_cutting      , &
               wshtotcutinit     , &
               lcutinit          , &
               wsh               , &
               wshtot            , &
               wr                , &
               c                 , &
               n                 , &
               napo              , &
               nsym              , &
               fn                , &
               tjulian           , &
               nel               , &
               biomass           , &
               devstage          , &
               regcount          , &
               wshcutinit        , &
               gmean             , &
               wc_frac                , &
               wnapo             , &
               wnsym             , &
               wgn               , &
               tasum             , &
               tgrowth           , &
               loss              , &
               lossc             , &
               lossn             , &
               tlossstart        , &
               lai               , &
               tcut              , &
               tcut_modif        , &
               wshtotsum         , &
               controle_azote_sum)

            ! ******************************************************
            ! mise � jour des concentrations
            ! ******************************************************

            WHERE ((wsh + wr .GT. 0.0)  .AND. (flag_cutting .EQ. 1) )
                
                c = wc_frac / (wsh + wr)
                n = (wnapo + wnsym) / (wsh + wr) 
                fn = wgn / (wr + wsh)
                napo = wnapo / (wsh + wr)
                nsym = wnsym / (wsh + wr)
                
            END WHERE

            WHERE (wshtot + wrtot .GT. 0.0)

                ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)

            END WHERE

          END IF

          WHERE (flag_cutting(:,:) .EQ. 1)

            when_growthinit_cut(:,:) = 0.0

          END WHERE

        ! SECOND test for automanagement lai & accumulated temperature over shreshold
        flag_cutting(:,:) = 0

        DO j=2,nvm
          IF (is_grassland_manag(j) .AND. (.NOT. is_grassland_cut(j)) .AND. &
            (.NOT.is_grassland_grazed(j)))THEN
        
            WHERE ((countschedule(:,j) .EQ. 1) .AND. (nanimal(:,j,1) .EQ. 0.0) .AND. &
              (devstage(:,j) .LT. 2.0) .AND. (tasum(:,j) .GE. tasumrep ) .AND. &
              (lai(:,j) .GE. 2.5))

              flag_cutting(:,j) = 1

              countschedule(:,j) = countschedule(:,j)  + 1

            END WHERE
          END IF
        END DO !nvm      

        ! If there is at least one point concerned
        IF (ANY (flag_cutting(:,:) .EQ. 1)) THEN 
            IF (blabla_pasim) PRINT *, 'FAUCHE AVEC METHODE NV', tjulian

            ! There will be one fertilization the day after cutting
            ! MODIF INN
            !courciruiter le calcul de tfert si f_fertiliZation = 0
            IF (f_fertilization.NE.1) THEN
              DO j=2,nvm
                DO i=1,npts
                  IF (flag_cutting(i,j) .EQ. 1) THEN
                       tfert(i,j,regcount(i,j) + 1) = tjulian + 1
                      print*, 'FERTILISATION AVEC METHODE NV', tjulian
                   END IF
                END DO
              END DO
            END IF

            CALL cutting_spa(&
               npts              , &
               tjulian           , &
               flag_cutting      , &
               wshtotcutinit     , &
               lcutinit          , &
               wsh               , &
               wshtot            , &
               wr                , &
               c                 , &
               n                 , &
               napo              , &
               nsym              , &
               fn                , &
               tjulian           , &
               nel               , &
               biomass           , &
               devstage          , &
               regcount          , &
               wshcutinit        , &
               gmean             , &
               wc_frac                , &
               wnapo             , &
               wnsym             , &
               wgn               , &
               tasum             , &
               tgrowth           , &
               loss              , &
               lossc             , &
               lossn             , &
               tlossstart        , &
               lai               , &
               tcut              , &
               tcut_modif        , &
               wshtotsum         , &
               controle_azote_sum)

            ! ******************************************************
            ! mise � jour des concentrations
            ! ******************************************************
            WHERE ((wsh + wr .GT. 0.0) .AND. (flag_cutting .EQ. 1))
                c = wc_frac / (wsh + wr)
                n = (wnapo + wnsym) / (wsh + wr) 
                fn = wgn / (wr + wsh)
                napo = wnapo / (wsh + wr)
                nsym = wnsym / (wsh + wr)
            END WHERE

            WHERE (wshtot + wrtot .GT. 0.0)
                ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)
            END WHERE

        END IF

      WHERE (flag_cutting(:,:) .EQ. 1)

        when_growthinit_cut(:,:) = 0.0

      END WHERE

      !If there are ncut cutting, it's finish
    
      WHERE (regcount(:,:) .EQ. ncut ) 
        cuttingend(:,:) = 1
      END WHERE
     
      !end of the cutting season by snow fall

    ELSE IF ((f_postauto .EQ.1 ) .OR. (f_autogestion .EQ. 3) .OR. &
      (f_autogestion .EQ. 4) &
    !! JCMODIF 290714 for postauto 5
      .OR. (f_postauto .GE. 2)) THEN

        flag_cutting(:,:) = 0
        when_growthinit_cut(:,:) = when_growthinit_cut(:,:) + dt

        ! FIRST test for automanagement
        WHERE((nanimal(:,mcut_C3,1) .EQ. 0.0) .AND. (cuttingend(:,mcut_C3) .EQ. 0) .AND. &
          (countschedule(:,mcut_C3) .EQ. 1) .AND. (((tgrowth(:,mcut_C3).GE.tgrowthmin) .AND. &
          (gmean(:,mcut_C3,ngmean).GT. 0.0) .AND. &
          (lai(:,mcut_C3) .GE. 2.5)  .AND. (devstage(:,mcut_C3) .GT. devstagemin ) .AND. &
          (gmeanslope(:,mcut_C3) .LT.gmeansloperel * mugmean(:,mcut_C3)))))

            flag_cutting(:,mcut_C3)  = 1
            countschedule(:,mcut_C3) = countschedule(:,mcut_C3)  + 1
        END WHERE

        WHERE((nanimal(:,mcut_C4,1) .EQ. 0.0) .AND. (cuttingend(:,mcut_C4) .EQ. 0).AND. &
          (countschedule(:,mcut_C4) .EQ. 1) .AND. (((tgrowth(:,mcut_C4).GE.tgrowthmin).AND. &
          (gmean(:,mcut_C4,ngmean).GT. 0.0) .AND. &
          (lai(:,mcut_C4) .GE. 2.5)  .AND. (devstage(:,mcut_C4) .GT. devstagemin ).AND. &
          (gmeanslope(:,mcut_C4) .LT.gmeansloperel * mugmean(:,mcut_C4)))))

            flag_cutting(:,mcut_C4)  = 1
            countschedule(:,mcut_C4) = countschedule(:,mcut_C4)  + 1
        END WHERE

        ! If there is at least one point concerned (flag_cutting = 1)

        IF ((ANY(flag_cutting(:,mcut_C3) .EQ. 1)) .OR. &
            (ANY(flag_cutting(:,mcut_C4) .EQ. 1))) THEN
            IF (blabla_pasim) PRINT *, 'FAUCHE AVEC METHODE NV ', tjulian
                ! There will be one fertilization the day after cutting
                ! A COURT-CIRCUITER si couplage autogestion ferti avec INN
                ! AIG 06/10/2009

            IF (f_fertilization.NE.1) THEN
             DO j=2,nvm
               DO i=1,npts
                  IF (flag_cutting(i,j) .EQ. 1) THEN

                      tfert(i,j,regcount(i,j) + 1) = tjulian + 1
                      print*, 'FERTILISATION AVEC METHODE NV', tjulian

                  END IF
                END DO
              END DO
            END IF
           CALL cutting_spa(&
               npts              , &
               tjulian           , &
               flag_cutting      , &
               wshtotcutinit     , &
               lcutinit          , &
               wsh               , &
               wshtot            , &
               wr                , &
               c                 , &
               n                 , &
               napo              , &
               nsym              , &
               fn                , &
               tjulian           , &
               nel               , &
               biomass           , &
               devstage          , &
               regcount          , &
               wshcutinit        , &
               gmean             , &
               wc_frac                , &
               wnapo             , &
               wnsym             , &
               wgn               , &
               tasum             , &
               tgrowth           , &
               loss              , &
               lossc             , &
               lossn             , &
               tlossstart        , &
               lai               , &
               tcut              , &
               tcut_modif        , &
               wshtotsum         , &
               controle_azote_sum)

            ! ******************************************************
            ! mise � jour des concentrations
            ! ******************************************************

            WHERE ((wsh + wr .GT. 0.0)  .AND. (flag_cutting .EQ. 1) )

                c = wc_frac / (wsh + wr)
                n = (wnapo + wnsym) / (wsh + wr)
                fn = wgn / (wr + wsh)
                napo = wnapo / (wsh + wr)
                nsym = wnsym / (wsh + wr)

            END WHERE

            WHERE (wshtot + wrtot .GT. 0.0)

                ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)

            END WHERE

        END IF

    WHERE (flag_cutting(:,mcut_C3) .EQ. 1)

        when_growthinit_cut(:,mcut_C3) = 0.0

    END WHERE
    WHERE (flag_cutting(:,mcut_C4) .EQ. 1)

        when_growthinit_cut(:,mcut_C4) = 0.0

    END WHERE


        ! SECOND test for automanagement
        flag_cutting(:,mcut_C3) = 0
        flag_cutting(:,mcut_C4) = 0

        WHERE ((countschedule(:,mcut_C3) .EQ. 1) .AND. (nanimal(:,mcut_C3,1) .EQ.0.0) .AND. &
          (devstage(:,mcut_C3) .LT. 2.0) .AND. (tasum(:,mcut_C3) .GE. tasumrep ) .AND. &
          (lai(:,mcut_C3) .GE. 2.5))

            flag_cutting(:,mcut_C3) = 1

            countschedule(:,mcut_C3) = countschedule(:,mcut_C3)  + 1

        END WHERE

        WHERE ((countschedule(:,mcut_C4) .EQ. 1) .AND. (nanimal(:,mcut_C4,1).EQ.0.0) .AND. &
          (devstage(:,mcut_C4) .LT. 2.0) .AND. (tasum(:,mcut_C4) .GE. tasumrep ).AND. &
          (lai(:,mcut_C4) .GE. 2.5))

            flag_cutting(:,mcut_C4) = 1

            countschedule(:,mcut_C4) = countschedule(:,mcut_C4)  + 1

        END WHERE

       ! If there is at least one point concerned
        IF ((ANY(flag_cutting(:,mcut_C3) .EQ. 1)) .OR. &
            (ANY(flag_cutting(:,mcut_C4) .EQ. 1))) THEN

            IF (blabla_pasim) PRINT *, 'FAUCHE AVEC METHODE NV', tjulian

            ! There will be one fertilization the day after cutting
            ! MODIF INN
            !courciruiter le calcul de tfert si f_fertiliZation = 0
            IF (f_fertilization.NE.1) THEN
              DO j=2,nvm
                DO i=1,npts
                  IF (flag_cutting(i,j) .EQ. 1) THEN
                       tfert(i,j,regcount(i,j) + 1) = tjulian + 1
                      print*, 'FERTILISATION AVEC METHODE NV', tjulian
                   END IF
                END DO
              END DO
            END IF

           CALL cutting_spa(&
               npts              , &
               tjulian           , &
               flag_cutting      , &
               wshtotcutinit     , &
               lcutinit          , &
               wsh               , &
               wshtot            , &
               wr                , &
               c                 , &
               n                 , &
               napo              , &
               nsym              , &
               fn                , &
               tjulian           , &
               nel               , &
               biomass           , &
               devstage          , &
               regcount          , &
               wshcutinit        , &
               gmean             , &
               wc_frac                , &
               wnapo             , &
               wnsym             , &
               wgn               , &
               tasum             , &
               tgrowth           , &
               loss              , &
               lossc             , &
               lossn             , &
               tlossstart        , &
               lai               , &
               tcut              , &
               tcut_modif        , &
               wshtotsum         , &
               controle_azote_sum)

            ! ******************************************************
            ! mise � jour des concentrations
            ! ******************************************************
           WHERE ((wsh + wr .GT. 0.0) .AND. (flag_cutting .EQ. 1))
                c = wc_frac / (wsh + wr)
                n = (wnapo + wnsym) / (wsh + wr)
                fn = wgn / (wr + wsh)
                napo = wnapo / (wsh + wr)
                nsym = wnsym / (wsh + wr)
            END WHERE

            WHERE (wshtot + wrtot .GT. 0.0)
                ntot = (wnapo + wnsym + wgn) / (wshtot + wrtot)
            END WHERE

        END IF

      WHERE (flag_cutting(:,mcut_C3) .EQ. 1)

        when_growthinit_cut(:,mcut_C3) = 0.0


      END WHERE

        !If there are ncut cutting, it's finish

        WHERE (regcount(:,mcut_C3) .EQ. ncut )
            cuttingend(:,mcut_C3) = 1
        END WHERE

    WHERE (flag_cutting(:,mcut_C4) .EQ. 1)

        when_growthinit_cut(:,mcut_C4) = 0.0


      END WHERE

        !If there are ncut cutting, it's finish

        WHERE (regcount(:,mcut_C4) .EQ. ncut )
            cuttingend(:,mcut_C4) = 1
        END WHERE

        !end of the cutting season by snow fall

      END IF auto_fauche
    END IF n_day_autofauche

    ! maintenant nous allons regarder les changements que le management apporte � Orchidee. 
    ! Dans un premier temps uniquement ceux sur le Carbone vu qu'Orchidee n'a pas d'azote.
    ! ******************************************************
    ! 1 - changement sur les r�servoirs du sol biologique. 
    ! ******************************************************
    CALL histwrite_p(hist_id_stomate ,'YIELD_RETURN',itime, YIELD_RETURN,npts*nvm, horipft_index)
    CALL chg_sol_bio(&
       npts                     , &
       tjulian                  , &
       bm_to_litter             , &
       litter                   , &
       litter_avail             , &
       litter_not_avail         , &
       litter_avail_totDM         , &
       intake_litter            , &
       biomass                  , &
       faecesc                  , &
       urinec                   , &
       fcOrganicFertmetabolic   , &
       fcOrganicFertstruct      , &
       fnOrganicFerturine       , &
       fnOrganicFertstruct      , &
       fnOrganicFertmetabolic    , &
       trampling                , &
       YIELD_RETURN)

    lai(:,:) = biomass(:,:,ileaf,icarbon)*sla_calc(:,:)
    ! ******************************************************
    ! 2 - changement sur les r�servoirs de la plante. 
    ! ******************************************************
    ! HISTWRITE
    regcount_real  = regcount
    fertcount_real = fertcount

    CALL histwrite_p(hist_id_stomate ,'REGCOUNT' ,itime ,regcount_real , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'FERTCOUNT',itime ,fertcount_real, npts*nvm, horipft_index)

    CALL histwrite_p(hist_id_stomate ,'GMEAN1',itime ,gmean(:,:,1) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN2',itime ,gmean(:,:,2) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN3',itime ,gmean(:,:,3) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN4',itime ,gmean(:,:,4) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN5',itime ,gmean(:,:,5) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN6',itime ,gmean(:,:,6) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN7',itime ,gmean(:,:,7) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN8',itime ,gmean(:,:,8) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN9',itime ,gmean(:,:,9) ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GMEAN0',itime ,gmean(:,:,10) ,npts*nvm, horipft_index)

    CALL histwrite_p(hist_id_stomate ,'WSH'   ,itime , wsh   , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WSHTOT',itime , wshtot, npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WR',    itime , wr,     npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WRTOT', itime , wrtot,  npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'WSHTOTSUM', itime , wshtotsum,  npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'SR_UGB', itime , sr_ugb,  npts*nvm,horipft_index)
    ! HISTWRITE POUR LA FERTIILSATION
!    CALL histwrite(hist_id_stomate ,'FCORGANICFERTMETABOLIC',itime , fcOrganicFertmetabolic,npts*nvm, horipft_index)
!    CALL histwrite(hist_id_stomate ,'FCORGANICFERTSTRUCT'   ,itime , fcOrganicFertstruct   ,npts*nvm, horipft_index)
!    CALL histwrite(hist_id_stomate ,'FNORGANICFERTURINE'    ,itime , fnOrganicFerturine    ,npts*nvm, horipft_index)
!    CALL histwrite(hist_id_stomate ,'FNORGANICFERTSTRUCT'   ,itime , fnOrganicFertstruct   ,npts*nvm, horipft_index)
!    CALL histwrite(hist_id_stomate ,'FNORGANICFERTMETABOLIC',itime , fnOrganicFertmetabolic,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'NFERTNITTOT'          ,itime , nfertnit(:,:,1)    ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'NFERTAMMTOT'          ,itime , nfertamm(:,:,1)    ,npts*nvm, horipft_index)
    ! HISTWRITE POUR LA FAUCHE
    CALL histwrite_p(hist_id_stomate ,'LOSS' ,itime ,loss  ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'LOSSC',itime ,lossc ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'LOSSN',itime ,lossn ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'DM_CUTYEARLY',itime ,DM_cutyearly ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'C_CUTYEARLY',itime ,C_cutyearly ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'NFERT_TOTAL',itime ,N_fert_total ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'Ndep',itime ,ndeposition ,npts*nvm,horipft_index)
    CALL histwrite_p(hist_id_stomate ,'LEGUME_FRACTION',itime ,legume_fraction ,npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'SOIL_FERTILITY',itime ,soil_fertility ,npts*nvm, horipft_index)

    CALL histwrite_p(hist_id_stomate ,'C'       ,itime, c       , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'N'       ,itime, n       , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'FN'      ,itime, fn      , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'NTOT'    ,itime, ntot    , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'NAPO'    ,itime, napo    , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'NSYM'    ,itime, nsym    , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'DEVSTAGE',itime, devstage, npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'TGROWTH' ,itime, tgrowth , npts*nvm, horipft_index)

    CALL histwrite_p(hist_id_stomate ,'GRAZINGCSTRUCT',itime, grazingcstruct      , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GRAZINGNSTRUCT',itime, grazingnstruct      , npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GRAZINGWN'     ,itime, Substrate_grazingwn, npts*nvm, horipft_index)
    CALL histwrite_p(hist_id_stomate ,'GRAZINGWC'     ,itime, Substrate_grazingwc, npts*nvm, horipft_index)

  END SUBROUTINE Main_Grassland_Management

  SUBROUTINE chg_sol_bio(&
     npts                     , &
     tjulian                  , &
     bm_to_litter             , &
     litter                   , &
     litter_avail             , &
     litter_not_avail         , &
     litter_avail_totDM         , &
     intake_litter            , &
     biomass                  , &
     faecesc                  , &
     urinec                   , &
     fcOrganicFertmetabolic    , &       
     fcOrganicFertstruct       , &
     fnOrganicFerturine        , &
     fnOrganicFertstruct       , &
     fnOrganicFertmetabolic    , &
     trampling                 , &
     YIELD_RETURN)

    INTEGER                                , INTENT(in)   :: npts
    INTEGER(i_std)                             , INTENT(in)   :: tjulian                 ! jour julien    
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: bm_to_litter 
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail 
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_not_avail 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout):: litter_avail_totDM
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in):: intake_litter
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(in)   :: faecesc
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(in)   :: urinec
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(in)   :: biomass           
    ! totalit� de masse s�che du shoot(kg/m**2)
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(in)   :: fcOrganicFertmetabolic
    ! metabolic C in slurry and manure (kg C/m**2/d)
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(in)   :: fcOrganicFertstruct 
    ! structural C in slurry and manure (kg C/m**2/d)
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFerturine    
    ! urine N in slurry and manure (kg N/m**2/d)
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFertstruct   
    ! structural N in slurry and manure (kg N/m**2/d)
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(in)   :: fnOrganicFertmetabolic  
    ! metabolic N in slurry and manure (kg N/m**2/d)           
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(in)   :: trampling
    REAL(r_std), DIMENSION(npts,nvm)            , INTENT(inout)   :: YIELD_RETURN

    REAL(r_std), DIMENSION(npts,nvm) :: litter_avail_totDM_old
    REAL(r_std), DIMENSION(npts,nvm) :: fcloss  
    REAL(r_std), DIMENSION(npts,nvm) :: fnloss
    REAL(r_std), DIMENSION(npts,nvm) :: floss
    REAL(r_std), DIMENSION(npts,nvm) :: fcplantsoil
    REAL(r_std), DIMENSION(npts,nvm) :: fnplantsoil
    REAL(r_std), DIMENSION(npts,nvm) :: fplantsoil
    REAL(r_std), DIMENSION(npts,nvm) :: l2nratio
    REAL(r_std), DIMENSION(npts,nvm) :: fmetabolic
    REAL(r_std), DIMENSION(npts,nvm) :: manure_barn
    INTEGER(i_std) :: j

    IF (blabla_pasim) PRINT *, 'PASIM main grassland : call chg_sol_bio'
    fmetabolic = 0.0
    DO j=2,nvm
      WHERE ((tjulian .GE. tlossstart(:,j)) .AND. (tjulian .LT. (tlossstart(:,j) + deltat/2.0))) 

        fcloss(:,j) = lossc(:,j)/deltat 
        fnloss(:,j) = lossn(:,j)/deltat 
        floss(:,j)  = loss(:,j) /deltat 

      ELSEWHERE

        fcloss(:,j) = 0.0
        fnloss(:,j) = 0.0
        floss(:,j)  = 0.0

      END WHERE

      fcplantsoil(:,j) = fcloss(:,j)
      fnplantsoil(:,j) = fnloss(:,j)
      fplantsoil(:,j)  = floss(:,j) 

      WHERE (fnplantsoil(:,j) .GT. 0.0) 

        l2nratio(:,j) = fligninresidue * fplantsoil(:,j)/ fnplantsoil(:,j)

      ELSEWHERE

        l2nratio(:,j) = 0.0

      END WHERE

      IF (is_grassland_cut(j).AND.(.NOT.is_grassland_grazed(j)))THEN
       
      ! Manure produced at barn
      ! 0.05 yieldloss 0.95 import_yield/harvest 0.85 loss during trasportation 
      ! 0.12 fraction of manure spread to field in total intake dry matter at barn (0.85*0.95harvest) 
      !JCcomment for accounting for Manurefert only not return 
      !        WHERE (YIELD_RETURN(:,:) .GT. 0.0)
      !          manure_barn(:,j) = fcplantsoil(:,j) / 0.05 * 0.95 * 0.85 *0.12 +&
      !            YIELD_RETURN(:,:) * CtoDM 
      !          YIELD_RETURN(:,:) = 0.0
      !        ELSEWHERE
      !          manure_barn(:,j) = fcplantsoil(:,j) / 0.05 * 0.95 * 0.85 *0.12
      !        ENDWHERE  
      manure_barn(:,j) = 0.0

      ELSE 
      manure_barn(:,j) = 0.0
      END IF
 
      fmetabolic(:,j) = MAX(0.625,MIN(0.85 - 0.018 * l2nratio(:,j), 1.0 - fligninresidue))
      bm_to_litter(:,j,ileaf,icarbon) = bm_to_litter(:,j,ileaf,icarbon) + &
         & fmetabolic(:,j) * (fcplantsoil(:,j) * 1000.0 + trampling(:,j))

      bm_to_litter(:,j,isapabove,icarbon) = bm_to_litter(:,j,isapabove,icarbon) + & 
         & (1.0 - fmetabolic(:,j)) * (fcplantsoil(:,j) * 1000.0 + trampling(:,j)) 
      litter_avail_totDM_old(:,j) = litter_avail_totDM(:,j)
      ! new litter available tot DM after intake litter
      litter_avail_totDM(:,j) = litter_avail_totDM(:,j) - intake_litter(:,j)
      IF (ANY(litter_avail_totDM(:,j) .LT. 0.0 ) ) THEN
        WRITE(numout,*) 'zd ','litter avail', j, litter_avail_totDM_old(:,j)
        WRITE(numout,*) 'zd ','intake litter', j, intake_litter(:,j)
        STOP 'available litter is not enough for grazing'

      ENDIF
      ! litter available C left is recalculated 
      ! assuming the same structural and metabolic fraction    
      WHERE (litter_avail_totDM_old(:,j) .GT. 0.0 )
      litter_avail(:,istructural,j) = litter_avail(:,istructural,j) * &
            & (litter_avail_totDM(:,j)/litter_avail_totDM_old(:,j))
      litter_avail(:,imetabolic,j) = litter_avail(:,imetabolic,j) * &
            & (litter_avail_totDM(:,j)/litter_avail_totDM_old(:,j))
      ELSEWHERE
      litter_avail(:,istructural,j) = litter_avail(:,istructural,j)
      litter_avail(:,imetabolic,j) = litter_avail(:,imetabolic,j)
      ENDWHERE
      ! new litter not available after manure/urine

      litter_not_avail(:,istructural,j) = litter_not_avail(:,istructural,j) + &
            & (faecesc(:,j) + urinec(:,j) + manure_barn(:,j) ) * 1000.0 * (1.0 - fmetabolic(:,j)) + &
            &  fcOrganicFertstruct(:,j) * 1000.0

      litter_not_avail(:,imetabolic,j) = litter_not_avail(:,imetabolic,j) + &
            & (faecesc(:,j) + urinec(:,j) + manure_barn(:,j) ) * 1000.0 * fmetabolic(:,j) + &
            &  fcOrganicFertmetabolic(:,j) * 1000.0
      ! update litter
      litter(:,:,j,iabove,icarbon) = litter_avail(:,:,j) + litter_not_avail(:,:,j)

    END DO
  END SUBROUTINE chg_sol_bio

  ! ________________________________________________________________
  ! ________________________________________________________________
  ! 
  ! Fonction permettant de lire les donn�es d'entr�es management.dat
  ! ________________________________________________________________
  ! ________________________________________________________________

  SUBROUTINE reading_new_animal(&
           npts           , &
           nb_year_management, &
           tcutmodel      , &
           tcut           , &
           tfert          , &
           nfertamm       , &
           nfertnit       , &
           nanimal        , &
           tanimal        , &
           danimal        , &
           nliquidmanure  , &
           nslurry        , &
           nsolidmanure   , &
           PIYcow         , &
           PIMcow         , &
           BCSYcow        , &
           BCSMcow        , &
           PICcow         , &
           AGE_cow_P      , &
           AGE_cow_M      , &
           Forage_quantity)

    INTEGER (i_std)                             , INTENT(in)  :: npts
    INTEGER(i_std)                              , INTENT(in) :: tcutmodel
    INTEGER(i_std),DIMENSION(nvm)             , INTENT(in) :: nb_year_management
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: tcut
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: tfert
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nfertamm
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nfertnit
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nanimal
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: tanimal
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: danimal
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nliquidmanure
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nslurry
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nsolidmanure
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: PIYcow
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: PIMcow
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: BCSYcow
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: BCSMcow
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: PICcow
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: AGE_cow_P
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: AGE_cow_M
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: Forage_quantity

    REAL(r_std), DIMENSION(nstocking)          :: nanimal_t
    REAL(r_std), DIMENSION(nstocking)          :: tanimal_t
    REAL(r_std), DIMENSION(nstocking)          :: danimal_t
    REAL(r_std), DIMENSION(nstocking)          :: tcut_t
    REAL(r_std), DIMENSION(nstocking)          :: tfert_t
    REAL(r_std), DIMENSION(nstocking)          :: nfertamm_t
    REAL(r_std), DIMENSION(nstocking)          :: nfertnit_t
    REAL(r_std), DIMENSION(nstocking)          :: nliquidmanure_t
    REAL(r_std), DIMENSION(nstocking)          :: nslurry_t
    REAL(r_std), DIMENSION(nstocking)          :: nsolidmanure_t

    REAL(r_std), DIMENSION(nstocking)          :: PIYcow_t
    REAL(r_std), DIMENSION(nstocking)          :: PIMcow_t
    REAL(r_std), DIMENSION(nstocking)          :: BCSYcow_t
    REAL(r_std), DIMENSION(nstocking)          :: BCSMcow_t
    REAL(r_std), DIMENSION(nstocking)          :: PICcow_t
    REAL(r_std), DIMENSION(nstocking)          :: AGE_cow_P_t
    REAL(r_std), DIMENSION(nstocking)          :: AGE_cow_M_t
    REAL(r_std), DIMENSION(nstocking)          :: Forage_quantity_t

    INTEGER(i_std)            :: ier, i, year, fin,j
    CHARACTER(len=200) :: description

    DO j=2,nvm
      OPEN(unit=60, file = file_management(j))

      READ(60, *   , iostat = ier) description
      read_management : IF (tcutmodel .EQ. 0) THEN
        IF (blabla_pasim) PRINT *, 'USERS MANAGEMENT'

        IF (nb_year_management(j) .LT. 1 ) STOP 'error with the nb_year_management'

        IF (MOD(count_year,nb_year_management(j))  .EQ. 0) THEN
            fin = nb_year_management(j)
        ELSE
            fin = MOD(count_year,nb_year_management(j))
        END IF

        DO year = 1, fin
            READ(60, *, iostat=ier) tcut_t(:)
            READ(60, *, iostat=ier) tfert_t(:)
            READ(60, *, iostat=ier) nfertamm_t(:)
            READ(60, *, iostat=ier) nfertnit_t(:)
            READ(60, *, iostat=ier) nanimal_t(:)
            READ(60, *, iostat=ier) tanimal_t(:)
            READ(60, *, iostat=ier) danimal_t(:)
            READ(60, *, iostat=ier) nliquidmanure_t(:)
            READ(60, *, iostat=ier) nslurry_t(:)
            READ(60, *, iostat=ier) nsolidmanure_t(:)

            READ(60, *, iostat=ier) PIYcow_t(:)
            READ(60, *, iostat=ier) PIMcow_t(:)
            READ(60, *, iostat=ier) BCSYcow_t(:)
            READ(60, *, iostat=ier) BCSMcow_t(:)
            READ(60, *, iostat=ier) PICcow_t(:)
            READ(60, *, iostat=ier) AGE_cow_P_t(:)
            READ(60, *, iostat=ier) AGE_cow_M_t(:)
            READ(60, *, iostat=ier) Forage_quantity_t(:)
          DO i=1,npts
            nanimal(i,j,:)=nanimal_t(:)
            tanimal(i,j,:)=tanimal_t(:)
            danimal(i,j,:)=danimal_t(:)
            tcut(i,j,:)=tcut_t(:)
            tfert(i,j,:)=tfert_t(:)
            nfertamm(i,j,:)=nfertamm_t(:)
            nfertnit(i,j,:)=nfertnit_t(:)
            nliquidmanure(i,j,:)=nliquidmanure_t(:)
            nslurry(i,j,:)=nslurry_t(:)
            nsolidmanure(i,j,:)=nsolidmanure_t(:)

            PIYcow(i,j,:)=PIYcow_t(:)
            PIMcow(i,j,:)=PIMcow_t(:)
            BCSYcow(i,j,:)=BCSYcow_t(:)
            BCSMcow(i,j,:)=BCSMcow_t(:)
            PICcow(i,j,:)=PICcow_t(:)
            AGE_cow_P(i,j,:)=AGE_cow_P_t(:)
            AGE_cow_M(i,j,:)=AGE_cow_M_t(:)
            Forage_quantity(i,j,:)=Forage_quantity_t(:)

          END DO
        END DO

      ELSE IF (tcutmodel .EQ. 1) THEN

        PRINT *, 'AUTO MANAGEMENT'
        READ(61, *, iostat = ier) toto(:)
        READ(61, *, iostat = ier) toto(:)
        READ(61, *, iostat = ier) toto(:)
        READ(61, *, iostat = ier) toto(:)

        READ(60,     *, iostat=ier)    nanimal_t
        DO i=1,npts
          nanimal(i,j,:)=nanimal_t(:)
        END DO
      ELSE

        STOP 'PASIM ERROR :: tcutmodel must be 0 or 1'

      END IF read_management
      CLOSE(60)
    END DO !nvm
  END SUBROUTINE reading_new_animal

  !!!!!!!!!!!!!!!!!!!!!!!
  !!!! 05212013 JCADD subroutine for reading management from map nc file
  !!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE reading_map_manag(&
           npts           , &
           count_year     , &
           nb_year_management, &
           management_intensity, &
           management_start, &
           tcut           , &
           tfert          , &
           nfertamm       , &
           nfertnit       , &
           nanimal        , &
           tanimal        , &
           danimal        , &
           nliquidmanure  , &
           nslurry        , &
           nsolidmanure   , &
           legume_fraction, &
           soil_fertility , &
           deposition_start, &
           ndeposition, &
           sr_ugb)

    INTEGER (i_std)                             , INTENT(in)  :: npts
    INTEGER (i_std)                             , INTENT(in)  :: count_year
    INTEGER(i_std),DIMENSION(nvm)             , INTENT(in) :: nb_year_management
    INTEGER(i_std),DIMENSION(nvm)             , INTENT(in) :: management_intensity
    INTEGER(i_std),DIMENSION(nvm)             , INTENT(in) :: management_start
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: tcut
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: tfert
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nfertamm
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nfertnit
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nanimal
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: tanimal
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: danimal
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nliquidmanure
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nslurry
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(out) :: nsolidmanure
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(out) :: legume_fraction
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(out) :: soil_fertility
    INTEGER(i_std),DIMENSION(nvm)             , INTENT(in) :: deposition_start
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(out) :: ndeposition
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(inout) :: sr_ugb
    ! new variables for get map of management
    REAL(r_std), DIMENSION(npts)        :: nfert_temp
    REAL(r_std), DIMENSION(npts)        :: nmanure_temp
    REAL(r_std), DIMENSION(npts)        :: nanimal_temp
    REAL(r_std), DIMENSION(npts)        :: tcut_temp
    REAL(r_std), DIMENSION(npts)        :: grazing_temp
    INTEGER(i_std)                      :: management_year
    INTEGER(i_std)                      :: deposition_year

    INTEGER(i_std)            :: j

      tcut(:,:,:) = 500.0
      tfert(:,:,:) = 500.0
      nfertamm(:,:,:) = 0.0
      nfertnit(:,:,:) = 0.0
      nanimal(:,:,:) = 0.0
      tanimal(:,:,:) = 500.0
      danimal(:,:,:) = 0.0
      nliquidmanure(:,:,:) = 0.0
      nslurry(:,:,:) = 0.0 
      nsolidmanure(:,:,:) = 0.0

      legume_fraction(:,:) =0.0
      soil_fertility(:,:) = 1.0
      ndeposition(:,:) = 0.0 

      nfert_temp(:) =0.0
      nmanure_temp(:) =0.0
      nanimal_temp(:) = 0.0
      DO j=2,nvm
        IF (nb_year_management(j) .LT. 1 ) STOP 'error with the nb_year_management'
        ! get which year of management should be read 
        IF (MOD(count_year,nb_year_management(j))  .EQ. 0) THEN
            management_year = nb_year_management(j) + management_start(j)-1
            deposition_year = nb_year_management(j) + deposition_start(j)-1
        ELSE
            management_year = MOD(count_year,nb_year_management(j)) + management_start(j)-1
            deposition_year = MOD(count_year,nb_year_management(j)) + deposition_start(j)-1

        END IF
        WRITE(numout,*)  management_year,deposition_year
        !!!! read deposition global file for all grassland including nature
        IF ( (.NOT. is_tree(j)) .AND. natural(j) .AND. (f_deposition_map .EQ. 1)) THEN
          !CALL get_map(npts,TRIM(deposition_map),deposition_year,"Ndep",ndeposition(:,j)) 
          !!print *,'ndep',ndeposition(:,j) 
        ELSE
          ndeposition(:,j)=0.0
        ENDIF
        !!!! read fertilization global file
        IF (management_intensity(j) .EQ. 4) THEN
          !CALL get_map(npts,TRIM(management_map),management_year,"Nmanure",nmanure_temp(:))
          nslurry(:,j,1)=nmanure_temp(:)/10000.
          !CALL get_map(npts,TRIM(management_map),management_year,"Nmineral",nfert_temp(:))
          nfertamm(:,j,1)= 0.5*nfert_temp(:)/10000.
          nfertnit(:,j,1)= 0.5*nfert_temp(:)/10000.
          
        ENDIF
        !!!! read sr_ugb global file
        IF (f_postauto .EQ. 5 .AND. f_grazing_map .EQ. 1) THEN
          !CALL get_map(npts,TRIM(grazing_map),management_year,"sr_ugb",grazing_temp(:))
          sr_ugb(:,mgraze_C3)=grazing_temp(:)/10000.
          sr_ugb(:,mgraze_C4)=grazing_temp(:)/10000.
        print *,'sr_ugb',sr_ugb(:,j)
          WHERE (sr_ugb(:,mgraze_C3) .GT. 10.)
            sr_ugb(:,mgraze_C3) = 10.1
            sr_ugb(:,mgraze_C4) = 10.1
          END WHERE
        print *, count_year,management_year
          if (ANY(sr_ugb(:,mgraze_C3) .EQ. 10.1)) then 
          print *, 'error sr_ugb',sr_ugb(:,mgraze_C3)
          endif
        ENDIF

        IF (management_intensity(j) .EQ. 1) THEN
          !CALL get_map(npts,TRIM(management_map),management_year,"tCut1L",tcut(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"tFert1L",tfert(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"NMinFert1L",nfert_temp(:))
          nfertamm(:,j,1)= 0.5*nfert_temp(:)/10000
          nfertnit(:,j,1)= 0.5*nfert_temp(:)/10000
          !CALL get_map(npts,TRIM(management_map),management_year,"tGrazingL",tanimal(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"dGrazingL",danimal(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"StockingRateL",nanimal_temp(:))
          nanimal(:,j,1)=nanimal_temp(:)/10000.
          print *,'tcutl',tcut(:,j,1)
          print *,'nanimall',nanimal(:,j,1)
          print *,'danimall',danimal(:,j,1)
          !CALL get_map(npts,TRIM(management_map),management_year,"NManureFert1L",nmanure_temp(:))
          nslurry(:,j,1)=nmanure_temp(:)/10000
          print *,'NManureFert1L',nslurry(:,j,1)

          !CALL get_map(npts,TRIM(management_map),management_year,"fracLegumesL",legume_fraction(:,j))
          !CALL get_map(npts,TRIM(fertility_map),management_year,"nutrient",soil_fertility(:,j))

        ELSEIF (management_intensity(j) .EQ. 2) THEN
!          !CALL get_map(npts,TRIM(management_map),management_year,"tCut1M",tcut(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"tFert1M",tfert(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"NMinFert1M",nfert_temp(:))
          nfertamm(:,j,1)= 0.5*nfert_temp(:)/10000
          nfertnit(:,j,1)= 0.5*nfert_temp(:)/10000
!          !CALL get_map(npts,TRIM(management_map),management_year,"tGrazingM",tanimal(:,j,1))
!          !CALL get_map(npts,TRIM(management_map),management_year,"dGrazingM",danimal(:,j,1))
!          !CALL get_map(npts,TRIM(management_map),management_year,"StockingRateM",nanimal_temp(:))
!          nanimal(:,j,1)=nanimal_temp(:)/10000.
          !CALL get_map(npts,TRIM(management_map),management_year,"NManureFert1M",nmanure_temp(:))
          nslurry(:,j,1)=nmanure_temp(:)/10000
          nslurry(:,mcut_C3,1)=nslurry(:,j,1)
          nslurry(:,mgraze_C3,1)=nslurry(:,j,1)
          nslurry(:,mcut_C4,1)=nslurry(:,j,1)
          nslurry(:,mgraze_C4,1)=nslurry(:,j,1)
!          tfert(:,mcut_C3,1)=tfert(:,j,1)
!          tfert(:,mgraze,1)=tfert(:,j,1)
          print *,'tcutm',tcut(:,j,1)
          print *,'nanimalm',nanimal(:,j,1)
          print *,'danimalm',danimal(:,j,1)
          print *,'NManureFert1M',nslurry(:,j,1)

          !CALL get_map(npts,TRIM(management_map),management_year,"fracLegumesM",legume_fraction(:,j))
          !CALL get_map(npts,TRIM(fertility_map),management_year,"nutrient",soil_fertility(:,j))

        ELSEIF (management_intensity(j) .EQ. 3) THEN
          !CALL get_map(npts,TRIM(management_map),management_year,"tCut1H",tcut(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"tFert1H",tfert(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"NMinFert1H",nfert_temp(:))
          nfertamm(:,j,1)= 0.5*nfert_temp(:)/10000
          nfertnit(:,j,1)= 0.5*nfert_temp(:)/10000
          !CALL get_map(npts,TRIM(management_map),management_year,"tGrazingH",tanimal(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"dGrazingH",danimal(:,j,1))
          !CALL get_map(npts,TRIM(management_map),management_year,"StockingRateH",nanimal_temp(:))
          nanimal(:,j,1)=nanimal_temp(:)/10000.
          !CALL get_map(npts,TRIM(management_map),management_year,"NManureFert1H",nmanure_temp(:))
          nslurry(:,j,1)=nmanure_temp(:)/10000
          !second cut and fert from 1981
          !CALL get_map(npts,TRIM(management_map),management_year,"tCut2H",tcut(:,j,2))
          !CALL get_map(npts,TRIM(management_map),management_year,"tFert2H",tfert(:,j,2))
          !CALL get_map(npts,TRIM(management_map),management_year,"NMinFert2H",nfert_temp(:))
          nfertamm(:,j,2)= 0.5*nfert_temp(:)/10000
          nfertnit(:,j,2)= 0.5*nfert_temp(:)/10000
          !CALL get_map(npts,TRIM(management_map),management_start,"NManureFert2H",nmanure_temp(:))
          nslurry(:,j,2)=nmanure_temp(:)/10000
          print *,'tcuth',tcut(:,j,1)
          print *,'nanimalh',nanimal(:,j,1)
          print *,'danimalm',danimal(:,j,1)
          print *,'NManureFert1H',nslurry(:,j,1)

          !CALL get_map(npts,TRIM(management_map),management_year,"fracLegumesH",legume_fraction(:,j))
          !CALL get_map(npts,TRIM(fertility_map),management_year,"nutrient",soil_fertility(:,j))

        ENDIF

      END DO ! nvm
    END SUBROUTINE reading_map_manag

    SUBROUTINE calc_N_limfert(&
             npts,nfertamm, nfertnit,&
             nliquidmanure, nslurry, nsolidmanure,&
             legume_fraction,soil_fertility,ndeposition,&
             N_fert_total,N_limfert)

    INTEGER (i_std)                             , INTENT(in)  :: npts
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nfertamm
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nfertnit
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nliquidmanure
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nslurry
    REAL(r_std), DIMENSION(npts,nvm,nstocking), INTENT(in) :: nsolidmanure
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(in) :: legume_fraction
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(in) :: soil_fertility
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(out) :: N_fert_total
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(out) :: N_limfert
    REAL(r_std), DIMENSION(npts,nvm)          , INTENT(in) :: ndeposition

    INTEGER(i_std) :: k,j

      N_fert_total(:,:) = 0.0 
      DO k=1,nstocking
        N_fert_total(:,:) = (N_fert_total(:,:) + nfertamm(:,:,k) + &
                            nfertnit(:,:,k) + nliquidmanure(:,:,k) + &
                            nslurry(:,:,k) + nsolidmanure(:,:,k))
      ENDDO
      N_fert_total(:,:) = N_fert_total(:,:) * 10000 + ndeposition(:,:)
      N_fert_total(:,1) = 0.0
      DO j=2,nvm
        IF ((management_intensity(j) .EQ. 2).AND. (.NOT. is_c4(j))) THEN
          N_fert_total(:,mcut_C3)=N_fert_total(:,j)
          N_fert_total(:,mgraze_C3)=N_fert_total(:,j)
        ENDIF
        IF ((management_intensity(j) .EQ. 2).AND. (is_c4(j))) THEN
          N_fert_total(:,mcut_C4)=N_fert_total(:,j)
          N_fert_total(:,mgraze_C4)=N_fert_total(:,j)
        ENDIF

      ENDDO
      !JCADD new fertilization effect
      ! linear
      !N_limfert(:,:) = 1.0 + (1.60-1.0)/320 * N_fert_total(:,:) 
      ! index
      N_limfert(:,:) = 1. + N_effect - N_effect * (0.75 ** (N_fert_total(:,:)/30))

      WHERE (N_limfert(:,:) .LT. 1.0) 
        N_limfert(:,:) = 1.0
      ELSEWHERE (N_limfert(:,:) .GT. 2.0)
        N_limfert(:,:) = 1.+N_effect
      ENDWHERE

    END SUBROUTINE calc_N_limfert

END MODULE Grassland_Management
