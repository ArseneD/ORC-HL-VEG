! =================================================================================================================================
! MODULE       : lpj_establish
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Establish pft's
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: 
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!        plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!        global vegetation model, Global Change Biology, 9, 161-185.\n
!! - Haxeltine, A. and I. C. Prentice (1996), BIOME3: An equilibrium
!!        terrestrial biosphere model based on ecophysiological constraints, 
!!        resource availability, and competition among plant functional types,
!!        Global Biogeochemical Cycles, 10(4), 693-709.\n
!! - Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!        dynamics in the modelling of terrestrial ecosystems: comparing two
!!        contrasting approaches within European climate space,
!!        Global Ecology and Biogeography, 10, 621-637.\n
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_stomate/lpj_establish.f90 $
!! $Date: 2014-02-11 17:16:59 +0100 (Tue, 11 Feb 2014) $
!! $Revision: 1892 $
!! \n
!_ ================================================================================================================================

MODULE lpj_establish

  ! modules used:
  USE xios_orchidee
  USE ioipsl_para
  USE stomate_data
  USE constantes

  IMPLICIT NONE

  ! private & public routines
  PRIVATE
  PUBLIC establish,establish_clear

  LOGICAL, SAVE                          :: firstcall = .TRUE.           !! first call
!$OMP THREADPRIVATE(firstcall)
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : fire_clear 
!!
!>\BRIEF       Set the firstcall flag to .TRUE. and activate initialization
!_ ================================================================================================================================

  SUBROUTINE establish_clear
    firstcall = .TRUE.
  END SUBROUTINE establish_clear


! =================================================================================================================================
! SUBROUTINE   : establish 
!
!>\BRIEF       Calculate sstablishment of new woody PFT and herbaceous PFTs
!!
!! DESCRIPTION : Establishments of new woody and herbaceous PFT are simulated. 
!! Maximum establishment rate (0.12) declines due to competition for light (space).
!! There are two establishment estimates: one for the for DGVM and one for the 
!! static cases.\n
!! In the case of DGVM, competitive process of establishment for the area of 
!! available space is represented using more detailed description compared with static 
!! one. Biomass and distribution of plant age are updated on the basis of changes 
!! in number of individuals. Finally excess sapwood of is converted to heartwood.
!!
!! \latexonly
!! \input{equation_lpj_establish.tex}
!! \endlatexonly
!! \n
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)    :
!! Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!    dynamics in the modelling of terrestrial ecosystems: comparing two
!!    contrasting approaches within European climate space,
!!    Global Ecology and Biogeography, 10, 621-637.
!!
!! FLOWCHART       : 
!! \latexonly
!! \includegraphics[scale = 0.7]{establish.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================
 
  SUBROUTINE establish (npts, dt, PFTpresent, regenerate, &
       neighbours, resolution, need_adjacent, herbivores, &
       precip_annual, gdd0, lm_lastyearmax, &
       cn_ind, lai, avail_tree, avail_grass,  npp_longterm, &
       leaf_age, leaf_frac, &
       ind, biomass, age, everywhere, co2_to_bm,veget_max, woodmass_ind, &
!JCADD
       sla_calc,height,dia_cut) !! Arsene 26-08-2015 - Add height (in) & !! Arsene 01-09-2015 - Add dia_cut 
!ENDJCADD 
    !! 0. Variable and parameter declaration
    
    !! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                                :: npts            !! Domain size - number of pixels (dimensionless)    
    REAL(r_std), INTENT(in)                                   :: dt              !! Time step of vegetation dynamics for stomate 
                                                                                 !! (days)            
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                  :: PFTpresent      !! Is pft there (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: regenerate      !! Winter sufficiently cold (unitless)   
    INTEGER(i_std), DIMENSION(npts,8), INTENT(in)             :: neighbours      !! indices of the 8 neighbours of each grid point 
                                                                                 !! (unitless);  
                                                                                 !! [1=N, 2=NE, 3=E, 4=SE, 5=S, 6=SW, 7=W, 8=NW]  
    REAL(r_std), DIMENSION(npts,2), INTENT(in)                :: resolution      !! resolution at each grid point (m); 1=E-W, 2=N-S     
    LOGICAL, DIMENSION(npts,nvm), INTENT(in)                  :: need_adjacent   !! in order for this PFT to be introduced, does it
                                                                                 !! have to be present in an adjacent grid box?  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: herbivores      !! time constant of probability of a leaf to 
                                                                                 !! be eaten by a herbivore (days)     
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: precip_annual   !! annual precipitation (mm year^{-1}) 
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: gdd0            !! growing degree days (degree C)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: lm_lastyearmax  !! last year's maximum leaf mass for each PFT 
                                                                                 !! (gC m^{-2 })
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: cn_ind          !! crown area of individuals (m^2)        
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: lai             !! leaf area index OF an individual plant 
                                                                                 !! (m^2 m^{-2})           
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: avail_tree      !! space availability for trees (unitless)     
    REAL(r_std), DIMENSION(npts), INTENT(in)                  :: avail_grass     !! space availability for grasses (unitless)  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: npp_longterm    !! longterm NPP, for each PFT (gC m^{-2})   
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: veget_max       !! "maximal" coverage fraction of a PFT 
                                                                                 !! (LAI -> infinity) on ground (unitless)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)                :: height          !! Arsene 26-08-2015 - Add for shrub if not_like_trees!! Height of vegetation (m)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: dia_cut         !! Arsene 01-09-2015 - Add dia_cut 


    !! 0.2 Output variables
    
    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_age        !! leaf age (days)     
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac       !! fraction of leaves in leaf age class (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind             !! Number of individuals (individuals m^{-2})          
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass   !! biomass (gC m^{-2 })     
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age             !! mean age (years)       
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere      !! is the PFT everywhere in the grid box or very 
                                                                                 !! localized (unitless)    
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm       !! biomass up take for establishment i.e. 
                                                                                 !! pseudo-photosynthesis (gC m^{-2} day^{-1}) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: woodmass_ind    !! woodmass of the individual, needed to calculate
                                                                                 !! crownarea in lpj_crownarea (gC m^{-2 })
!JCADD
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)            :: sla_calc
!ENDJCADD
    !! 0.4 Local variables

    REAL(r_std)                                               :: tau_eatup       !! time during which a sapling can be entirely 
                                                                                 !! eaten by herbivores (days) 
    REAL(r_std), DIMENSION(npts,nvm)                          :: fpc_nat         !! new fpc, foliage projective cover: fractional
                                                                                 !! coverage (unitless)       
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_climate_tree  !! maximum tree establishment rate, 
                                                                                              !! based on climate only (unitless)  
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_climate_grass !! maximum grass establishment rate,
                                                                                              !! based on climate only (unitless)
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_tree          !! maximum tree establishment rate, 
                                                                                              !! based on climate and fpc 
                                                                                              !! (unitless) 
    REAL(r_std), DIMENSION(npts)                              :: estab_rate_max_grass         !! maximum grass establishment rate,
                                                                                              !! based on climate and fpc 
                                                                                              !! (unitless) 
    REAL(r_std), DIMENSION(npts)                              :: sumfpc          !! total natural fpc (unitless)
    REAL(r_std), DIMENSION(npts)                              :: fracnat         !! total fraction occupied by natural 
                                                                                 !! vegetation (unitless)  
    REAL(r_std), DIMENSION(npts)                              :: sumfpc_wood     !! total woody fpc (unitless)    
    REAL(r_std), DIMENSION(npts)                              :: spacefight_tree !! for trees, measures the total concurrence for
                                                                                 !! available space (unitless)      
    REAL(r_std), DIMENSION(npts)                              :: spacefight_grass!! for grasses, measures the total concurrence 
                                                                                 !! for available space (unitless)
    REAL(r_std), DIMENSION(npts,nvm)                          :: d_ind           !! change in number of individuals per time step 
                                                                                 !! (individuals m^{-2} day{-1})         
    REAL(r_std), DIMENSION(npts)                              :: bm_new          !! biomass increase (gC m^{-2 })        
    REAL(r_std), DIMENSION(npts)                              :: dia             !! stem diameter (m)    
    REAL(r_std), DIMENSION(npts)                              :: b1              !! temporary variable           
    REAL(r_std), DIMENSION(npts)                              :: woodmass        !! woodmass of an individual (gC m^{-2})  
    REAL(r_std), DIMENSION(npts)                              :: leaf_mass_young !! carbon mass in youngest leaf age class 
                                                                                 !! (gC m^{-2})
    REAL(r_std), DIMENSION(npts)                              :: factor          !! reduction factor for establishment if many 
                                                                                 !! trees or grasses are present (unitless)  
    REAL(r_std), DIMENSION(npts)                              :: total_bm_c      !! Total carbon mass for all pools (gC m^{-2})    
    REAL(r_std), DIMENSION(npts,nelements)                    :: total_bm_sapl   !! Total sappling biomass for all pools 
                                                                                 !! (gC m^{-2})  
    INTEGER(i_std)                                            :: nfrontx         !! from how many sides is the grid box invaded
                                                                                 !! (unitless?)   
    INTEGER(i_std)                                            :: nfronty         !! from how many sides is the grid box invaded
                                                                                 !! (unitless?)   
   !LOGICAL, DIMENSION(npts)                                  :: many_new        !! daily establishment rate is large compared to 
                                                                                 !! present number of individuals
    REAL(r_std), DIMENSION(npts)                              :: vn              !! flow due to new individuals veget_max after 
                                                                                 !! establishment, to get a proper estimate of 
                                                                                 !! carbon and nitrogen 
    REAL(r_std), DIMENSION(npts)                              :: lai_ind         !! lai on each PFT surface (m^2 m^{-2})   
    INTEGER(i_std)                                            :: i,j,k,m         !! indices (unitless)       
!
!    REAL(r_std)                                               :: ind_max         !! 2. * fraction of real individu by ind (special for shrub) !! 22-05-2015 Arsene

    REAL(r_std)                                       :: volume, signe, fact, num_it, volume1 !! Arsene 11-08-2015 - Add for iteration -signe_presc, 
    LOGICAL                                           :: dia_ok, init_ok !! Arsene 11-08-2015 - Add for iteration

!
!_ ================================================================================================================================

    IF (bavard.GE.3) WRITE(numout,*) 'Entering establish'

  !! 1. messages and initialization

    ! Assumption: time durasion that sapling is completely eaten by hervioures is a half year?   
    ! No reference
    tau_eatup = one_year/2.

    !! 1.1 First call only
    IF ( firstcall ) THEN

       WRITE(numout,*) 'establish:'

       WRITE(numout,*) '   > time during which a sapling can be entirely eaten by herbivores (d): ', &
            tau_eatup

       firstcall = .FALSE.

    ENDIF

  !! 2. recalculate fpc

    IF (control%ok_dgvm) THEN
       fracnat(:) = un

       !! 2.1 Only natural part of the grid cell
       do j = 2,nvm ! Loop over # PFTs
          
          IF ( .NOT. natural(j) ) THEN
             fracnat(:) = fracnat(:) - veget_max(:,j)
          ENDIF
       ENDDO ! Loop over # PFTs
       
       sumfpc(:) = zero

       !! 2.2 Total natural fpc on grid
       !      The overall fractional coverage of a PFT in a grid is calculated here.
       !      FPC is related to mean individual leaf area index by the Lambert-Beer law.
       !      See Eq. (1) in tex file.\n
       DO j = 2,nvm ! Loop over # PFTs
          IF ( natural(j) ) THEN
             WHERE(fracnat(:).GT.min_stomate)
                WHERE (lai(:,j) == val_exp) 
                   fpc_nat(:,j) = cn_ind(:,j) * ind(:,j) / fracnat(:)
                ELSEWHERE
                   fpc_nat(:,j) = cn_ind(:,j) * ind(:,j) / fracnat(:) & 
!JCMODIF
!                        * ( un - exp( - lm_lastyearmax(:,j) * sla(j) * ext_coeff(j) ) )
                        * ( un - exp( - lm_lastyearmax(:,j) * sla_calc(:,j) * ext_coeff(j) ) )
!ENDJCMODIF
                ENDWHERE
             ENDWHERE

             WHERE ( PFTpresent(:,j) )
                sumfpc(:) = sumfpc(:) + fpc_nat(:,j)
             ENDWHERE
          ELSE

             fpc_nat(:,j) = zero

          ENDIF

       ENDDO ! Loop over # PFTs
       
       !! 2.3 Total woody fpc on grid and number of regenerative tree pfts
       !      Total woody FPC increases by adding new FPC.
       !      Under the condition that temperature in last winter is higher than a threshold, 
       !      woody plants is exposed in higher competitive environment.
       sumfpc_wood(:) = zero
       spacefight_tree(:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          IF ( ( is_tree(j) .OR. is_shrub(j) ) .AND. natural(j) ) THEN       !! Arsene 31-07-2014 modifications

             ! total woody fpc
             WHERE ( PFTpresent(:,j) )
                sumfpc_wood(:) = sumfpc_wood(:) + fpc_nat(:,j)
             ENDWHERE

             ! how many trees are competing? Count a PFT fully only if it is present
             ! on the whole grid box.
!             WHERE ( PFTpresent(:,j) .AND. ( regenerate(:,j) .GT. regenerate_crit ) )
!                spacefight_tree(:) = spacefight_tree(:) + everywhere(:,j)
!             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

       !! 2.4 Total number of natural grasses on grid\n
       !     Grass increment equals 'everywhere'\n
       spacefight_grass(:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          IF ( .NOT. is_tree(j) .AND. .NOT. is_shrub(j) .AND. natural(j) ) THEN    !! Arsene 31-07-2014 modifications

             ! Count a PFT fully only if it is present on a grid.
             WHERE ( PFTpresent(:,j) )
                spacefight_grass(:) = spacefight_grass(:) + everywhere(:,j)
             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

       !! 2.5 Maximum establishment rate, based on climate only\n
       WHERE ( ( precip_annual(:) .GE. precip_crit ) .AND. ( gdd0(:) .GE. gdd_crit_estab ) )

          estab_rate_max_climate_tree(:) = estab_max_tree ! 'estab_max_*'; see 'stomate_constants.f90'
          estab_rate_max_climate_grass(:) = estab_max_grass

       ELSEWHERE

          estab_rate_max_climate_tree(:) = zero
          estab_rate_max_climate_grass(:) = zero

       ENDWHERE

       !! 2.6 Reduce maximum tree establishment rate if many trees are present.
       !      In the original DGVM, this is done using a step function which yields a
       !      reduction by factor 4 if sumfpc_wood(i) .GT.  fpc_crit - 0.05.
       !      This can lead to small oscillations (without consequences however).
       !      Here, a steady linear transition is used between fpc_crit-0.075 and
       !      fpc_crit-0.025.
       !      factor(:) = 1. - 15. * ( sumfpc_wood(:) - (fpc_crit-.075))
       !      factor(:) = MAX( 0.25_r_std, MIN( 1._r_std, factor(:)))
       !      SZ modified according to Smith et al. 2001
       !      See Eq. (2) in header
       factor(:)=(un - exp(- establish_scal_fact * (un - sumfpc_wood(:))))*(un - sumfpc_wood(:))
       estab_rate_max_tree(:) = estab_rate_max_climate_tree(:) * factor(:)

       !! 2.7 Modulate grass establishment rate.
       !      If canopy is not closed (fpc < fpc_crit-0.05), normal establishment.
       !      If canopy is closed, establishment is reduced by a factor 4.
       !      Factor is linear between these two bounds.
       !      This is different from the original DGVM where a step function is
       !      used at fpc_crit-0.05 (This can lead to small oscillations,
       !      without consequences however).
       !      factor(:) = 1. - 15. * ( sumfpc(:) - (fpc_crit-.05))
       !      factor(:) = MAX( 0.25_r_std, MIN( 1._r_std, factor(:)))
       !      estab_rate_max_grass(:) = estab_rate_max_climate_grass(:) * factor(:)
       !      SZ modified to true LPJ formulation, grasses are only allowed in the
       !      fpc fraction not occupied by trees..., 080806
       !      estab_rate_max_grass(:)=MAX(0.98-sumfpc(:),zero)
       !      See Eq. (3) in header
       estab_rate_max_grass(:) = MAX(MIN(estab_rate_max_climate_grass(:), max_tree_coverage - sumfpc(:)),zero)

       !! 2.8 Longterm grass NPP for competition between C4 and C3 grasses
       !      to avoid equal veget_max, the idea is that more reestablishment
       !      is possible for the more productive PFT
       factor(:) = min_stomate
       DO j = 2,nvm ! Loop over # PFTs
          IF ( natural(j) .AND. .NOT.is_tree(j) .AND. .NOT.is_shrub(j) ) &           !! Arsene 31-07-2014 modifications
               factor(:) = factor(:) + npp_longterm(:,j) * &
!JCMODIF
!               lm_lastyearmax(:,j) * sla(j)
               lm_lastyearmax(:,j) * sla_calc(:,j)
!ENDJCMODIF
       ENDDO ! Loop over # PFTs

       !! 2.9 Establish natural PFTs
       d_ind(:,:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          IF ( natural(j) ) THEN ! only for natural PFTs

             !! 2.9.1 PFT expansion across the grid box. Not to be confused with areal coverage.
             IF ( treat_expansion ) THEN

                ! only treat plants that are regenerative and present and still can expand
                DO i = 1, npts ! Loop over # pixels - domain size

                   IF ( PFTpresent(i,j) .AND. &
                        ( everywhere(i,j) .LT. un ) .AND. &
                        ( regenerate(i,j) .GT. regenerate_crit ) ) THEN

                      ! from how many sides is the grid box invaded (separate x and y directions
                      ! because resolution may be strongly anisotropic)
                      ! For the moment we only look into 4 direction but that can be expanded (JP) 
                      nfrontx = 0
                      IF ( neighbours(i,3) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,3),j) .GT. 1.-min_stomate ) nfrontx = nfrontx+1
                      ENDIF
                      IF ( neighbours(i,7) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,7),j) .GT. 1.-min_stomate ) nfrontx = nfrontx+1
                      ENDIF

                      nfronty = 0
                      IF ( neighbours(i,1) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,1),j) .GT. 1.-min_stomate ) nfronty = nfronty+1
                      ENDIF
                      IF ( neighbours(i,5) .GT. 0 ) THEN
                         IF ( everywhere(neighbours(i,5),j) .GT. 1.-min_stomate ) nfronty = nfronty+1
                      ENDIF
                      
                      everywhere(i,j) = &
                           everywhere(i,j) + migrate(j) * dt/one_year * &
                           ( nfrontx / resolution(i,1) + nfronty / resolution(i,2) )
                      
                      IF ( .NOT. need_adjacent(i,j) ) THEN
                         
                         ! in that case, we also assume that the PFT expands from places within
                         ! the grid box (e.g., oasis).
                         ! What is this equation? No reference.
                         everywhere(i,j) = &
                              everywhere(i,j) + migrate(j) * dt/one_year * &
                              2. * SQRT( pi*everywhere(i,j)/(resolution(i,1)*resolution(i,2)) )

                      ENDIF

                      everywhere(i,j) = MIN( everywhere(i,j), un )

                   ENDIF

                ENDDO ! Loop over # pixels - domain size

             ENDIF ! treat expansion?

             !! 2.9.2 Establishment rate
             !      - Is lower if the PFT is only present in a small part of the grid box
             !        (after its introduction), therefore multiplied by "everywhere".
             !      - Is divided by the number of PFTs that compete ("spacefight").
             !      - Is modulated by space availability (avail_tree, avail_grass).

             !! 2.9.2.1 present and regenerative trees
             IF ( is_tree(j) .OR. is_shrub(j) ) THEN       !! Arsene 31-07-2014 modifications A vérif pour d_ind   utilise avail_tree et estab_rate_max_tree et spacefight_tree

                WHERE ( PFTpresent(:,j) .AND. ( regenerate(:,j) .GT. regenerate_crit ) )
                   
                   
                   !d_ind(:,j) = estab_rate_max_tree(:)*everywhere(:,j)/spacefight_tree(:) * &
                   d_ind(:,j) = estab_rate_max_tree(:)*everywhere(:,j) * &
                        avail_tree(:) * dt/one_year

                ENDWHERE

             !! 2.9.2.2 present and regenerative grasses
             ELSE

                WHERE ( PFTpresent(:,j) .AND. ( regenerate(:,j) .GT. regenerate_crit )  & 
                     .AND.factor(:).GT.min_stomate .AND. spacefight_grass(:).GT. min_stomate) 
                   
                   d_ind(:,j) = estab_rate_max_grass(:)*everywhere(:,j)/spacefight_grass(:) * &
!JCMODIF
!                        MAX(min_stomate,npp_longterm(:,j)*lm_lastyearmax(:,j)*sla(j)/factor(:)) * fracnat(:) * dt/one_year
                   MAX(min_stomate,npp_longterm(:,j)*lm_lastyearmax(:,j) * &
                   sla_calc(:,j)/factor(:)) * fracnat(:) * dt/one_year
!ENDJCMODIF                   
                ENDWHERE
                
             ENDIF  ! tree/grass
             
          ENDIF ! if natural
       ENDDO  ! Loop over # PFTs
       
  !! 3. Lpj establishment in static case 

    !     Lpj establishment in static case, SZ 080806, account for real LPJ dynamics in
    !     prescribed vegetation, i.e. population dynamics within a given area of the grid cell.
    ELSE 

       d_ind(:,:) = zero

       DO j = 2,nvm ! Loop over # PFTs

          WHERE(ind(:,j)*cn_ind(:,j).GT.min_stomate)
!JCMODIF
!             lai_ind(:) = sla(j) * lm_lastyearmax(:,j)/(ind(:,j)*cn_ind(:,j))
             lai_ind(:) = sla_calc(:,j) * lm_lastyearmax(:,j)/(ind(:,j)*cn_ind(:,j))
!ENDJCMODIF
          ELSEWHERE
             lai_ind(:) = zero
          ENDWHERE

          !! 3.1 For natural woody PFTs
          IF ( natural(j) .AND. ( is_tree(j) .OR. is_shrub(j) ) ) THEN           !! Arsene 31-07-2014 modifications ok

             ! See Eq. (4) in tex file.            
             fpc_nat(:,j) =  MIN(un, cn_ind(:,j) * ind(:,j) * & 
                  MAX( ( un - exp( - ext_coeff(j) * lai_ind(:) ) ), min_cover ) )

!             IF ( is_tree(j) ) THEN              !! Arsene 22-05-2015 Add
!                 ind_max = 2.                    !! Arsene 22-05-2015 Add - like orignal = 2.
!             ELSE                                !! Arsene 22-05-2015 Add
!                 ind_max = 2. * shrub_ind_frac   !! Arsene 22-05-2015 Add - original * [ratio branches (ind) & real ind]
!             ENDIF                               !! Arsene 22-05-2015 Add

             WHERE (veget_max(:,j).GT.min_stomate.AND.ind(:,j).LE.2.)    !! Arsene 22-05-2015. Original. Note that this where and the next one are in "competition"
             !! Arsene 20-07-2015 On pourrait envisager que le nombre d'individus soit >> à 2.

!             WHERE (veget_max(:,j).GT.min_stomate.AND.ind(:,j).LE.ind_max)   !! Arsene 22-05-2015 Rempalce 2. by ind_max

                !! 3.1.1 Only establish into growing stands 
                !        Only establish into growing stands, ind can become very
                !        large in the static mode because LAI is very low in poor 
                !        growing conditions, favouring continuous establishment. 
                !        To avoid this a maximum IND is set. BLARPP: This should be
                !        replaced by a better stand density criteria.
                factor(:)=(un - exp(-establish_scal_fact * (un - fpc_nat(:,j))))*(un - fpc_nat(:,j))

                estab_rate_max_tree(:) = estab_max_tree * factor(:) !! Arsene 31-07-2014 modifications a vérif les variables associées a tree

                !! 3.1.2 do establishment for natural PFTs\n
                d_ind(:,j) = MAX( zero, estab_rate_max_tree(:) * dt/one_year)      !! Arsene 22-05-2015 Maybe have to create news parameters for shrubs !!!

             ENDWHERE

             !SZ: quickfix: to simulate even aged stand, uncomment the following lines...
             !where (ind(:,j) .LE. min_stomate)
             !d_ind(:,j) = 0.1 !MAX( 0.0, estab_rate_max_tree(:) * dt/one_year)
!             IF ( .NOT. is_shrub(j) ) THEN                                                      !! Arsene 22-05-2015. Change the number of start indivual number with branche ratio (for shrubs)
                WHERE (veget_max(:,j).GT.min_stomate .AND. ind(:,j).EQ.zero)
                   d_ind(:,j) = ind_0_estab
                ENDWHERE
!             ELSE                                                                               !! Arsene 22-05-2015. Change the number of start indivual number with branche ratio (for shrubs)
!                WHERE (veget_max(:,j).GT.min_stomate .AND. ind(:,j).EQ.zero)                    !! Arsene 22-05-2015. Change the number of start indivual number with branche ratio (for shrubs)
!                   d_ind(:,j) = ind_0_estab * shrub_ind_frac                                    !! Arsene 22-05-2015. Change the number of start indivual number with branche ratio (for shrubs)
!                ENDWHERE                                                                        !! Arsene 22-05-2015. Change the number of start indivual number with branche ratio (for shrubs)



!             ENDIF                                                                              !! Arsene 22-05-2015. Change the number of start indivual number with branche ratio (for shrubs)

          !! 3.2 For natural grass PFTs
          ELSEIF ( natural(j) .AND. .NOT.is_tree(j) .AND. .NOT.is_shrub(j) ) THEN         !! Arsene 31-07-2014 modifications ok

             WHERE (veget_max(:,j).GT.min_stomate)

                fpc_nat(:,j) =  cn_ind(:,j) * ind(:,j) * & 
                     MAX( ( un - exp( - ext_coeff(j) * lai_ind(:) ) ), min_cover )

                d_ind(:,j) = MAX(zero , (un - fpc_nat(:,j)) * dt/one_year )

             ENDWHERE

             WHERE (veget_max(:,j).GT.min_stomate .AND. ind(:,j).EQ. zero)
                d_ind(:,j) = ind_0_estab 
             ENDWHERE

          ENDIF

       ENDDO ! Loop over # PFTs

    ENDIF ! DGVM OR NOT

  !! 4. Biomass calculation

    DO j = 2,nvm ! Loop over # PFTs

       IF ( natural(j) ) THEN ! only for natural PFTs

          !! 4.1 Herbivores reduce establishment rate
          !      We suppose that saplings are vulnerable during a given time after establishment.
          !      This is taken into account by preventively reducing the establishment rate.
          IF ( ok_herbivores ) THEN

             d_ind(:,j) = d_ind(:,j) * EXP( - tau_eatup/herbivores(:,j) )

          ENDIF

          !! 4.2 Total biomass.
          !      Add biomass only if d_ind, over one year, is of the order of ind.
          !      save old leaf mass to calculate leaf age
          leaf_mass_young(:) = leaf_frac(:,j,1) * biomass(:,j,ileaf,icarbon)

          ! total biomass of existing PFT to limit biomass added from establishment
          total_bm_c(:) = zero

          DO k = 1, nparts
             total_bm_c(:) = total_bm_c(:) + biomass(:,j,k,icarbon)
          ENDDO
          IF(control%ok_dgvm) THEN
             vn(:) = veget_max(:,j)
          ELSE
             vn(:) = un
          ENDIF

          total_bm_sapl(:,:)=zero
          DO k = 1, nparts
             WHERE(d_ind(:,j).GT.min_stomate.AND.vn(:).GT.min_stomate)

                total_bm_sapl(:,icarbon) = total_bm_sapl(:,icarbon) + & 
                     bm_sapl(j,k,icarbon) * d_ind(:,j) / vn(:)
             ENDWHERE
          ENDDO


          !! 4.3 Woodmass calculation

          !! 4.3.1 with DGVM
          IF(control%ok_dgvm) THEN

             ! SZ calculate new woodmass_ind and veget_max after establishment (needed for correct scaling!)
             ! essential correction for MERGE!
             IF(is_tree(j) .OR. is_shrub(j))THEN                      !! Arsene 31-07-2014 modifications Vérif woodmass
                DO i=1,npts ! Loop over # pixels - domain size
                   IF((d_ind(i,j)+ind(i,j)).GT.min_stomate) THEN

                      IF((total_bm_c(i).LE.min_stomate) .OR. (veget_max(i,j) .LE. min_stomate)) THEN

                         ! new wood mass of PFT
                         woodmass_ind(i,j) = &
                              (((biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) &
                              + biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon))*veget_max(i,j)) &
                              + (bm_sapl(j,isapabove,icarbon) + bm_sapl(j,isapbelow,icarbon) &
                              + bm_sapl(j,iheartabove,icarbon) + bm_sapl(j,iheartbelow,icarbon))*d_ind(i,j))/(ind(i,j) + d_ind(i,j))

                      ELSE
 
                         ! new biomass is added to the labile pool, hence there is no change 
                         ! in CA associated with establishment
                         woodmass_ind(i,j) = &
                              & (biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) &
                              & +biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon))*veget_max(i,j) &
                              & /(ind(i,j) + d_ind(i,j))

                      ENDIF

                      ! new diameter of PFT

                      IF (is_tree(j)) THEN !! Arsene 03-08-2015 - Add allometry for shrubs
                          dia(i) = (woodmass_ind(i,j)/(pipe_density*pi/4.*pipe_tune2)) &
                              &       **(1./(2.+pipe_tune3))
                          vn(i) = (ind(i,j) + d_ind(i,j))*pipe_tune1* &
                                  & MIN(MAX(dia(i),dia_cut(i,j)),maxdia(j))**pipe_tune_exp_coeff           !! Arsene 01-09-2015 - Add dia_cut

                      ELSEIF (is_shrub(j) .AND. shrubs_like_trees) THEN                 !! Arsene 03-08-2015 - Add allometry for shrubs
                          dia(i) = (woodmass_ind(i,j)/(pipe_density_shrub*pi/4.*pipe_tune2_for_shrub)) &
                               &            **(1./(2.+pipe_tune3_for_shrub))
                          vn(i) = (ind(i,j) + d_ind(i,j))*pipe_tune1_for_shrub* &
                                  & MIN(MAX(dia(i),dia_cut(i,j)),maxdia(j))**pipe_tune_exp_coeff_for_shrub !! Arsene 01-09-2015 - Add dia_cut

                      ELSE !! Arsene 12-08-2015 - If shrub and New Allometry
                          !! To calculate diameter and cover... it's more complicate...

                          !! 1. Calculate the exact volume
                           volume = woodmass_ind(i,j) / pipe_density_shrub

                          !! 2. Stem diameter from Aiba & Kohyama
                 !          Stem diameter is calculated by allometory... but no analytique solution for this equation:
                 !! volume(i) = pi/4 * height_presc(j) * pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia(i)**(pipe_tune_shrub3+2) &
                 !!               & / ( height_presc(j) + pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia(i)**pipe_tune_shrub3 )

!                 !! On part de la hauteur corresonpodant au minheight... ==> min dia !! Peut être que ça marche "bien" du premier coup (on peut toujours réver...)
!                 dia(i) = mindia(j)
                 !! On part de la hauteur initiale, et on essaye de trouver le bon diamètre...
                 IF ( height(i,j) .LT. height_presc(j) ) THEN  !! Arsene 11-08-2015 Garde fou
                     dia(i) = (height(i,j)*height_presc(j) / (pipe_tune_shrub2*(height_presc(j)-height(i,j))) ) &
                                 & **(1/pipe_tune_shrub3) /100
                 ELSE
                     dia(i) = maxdia(j)
                 ENDIF

                 num_it = 0
                 signe = 0
                 init_ok = .false.
                 dia_ok = .false.
                 DO WHILE ( .NOT.dia_ok )
                    IF ( num_it .EQ. 1) init_ok=.TRUE.
                    volume1 = pi/4 * height_presc(j) * pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia(i)**(pipe_tune_shrub3+2.) &
                                & / ( height_presc(j) + pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia(i)**pipe_tune_shrub3 )

                    IF ( ABS(volume1-volume)/volume.GT.1. .AND. (num_it.EQ.0)) THEN
                       fact = 0.3
                    ELSEIF ( ABS(volume1-volume)/volume.GT.0.5 .AND. (num_it.EQ.0)) THEN
                       fact = dia(i)
                    ELSEIF ( ABS(volume1-volume)/volume.GT.0.1 .AND. (num_it.EQ.0)) THEN
                       fact = dia(i)*0.1
                    ELSEIF ( num_it.EQ.0 ) THEN !IF ( ABS(volume2-volume1).GT.0.00001 ) THEN
                       fact = dia(i)*0.01
                    ENDIF

                    IF  ( (volume1-volume) .GT. (accept_sigma_it*volume) ) THEN  !! So Dia to important or min_stomate
                       IF ( signe.EQ.1. .AND. .NOT.init_ok) THEN
                          fact = fact / factor_div_it 
                           signe = -1.
                          IF ( fact .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                             dia_ok = .true.
                          ENDIF
                       ENDIF
                       DO WHILE ( (dia(i)-fact).LT.min_stomate .AND. ((signe.EQ.-1.).OR.(.NOT.init_ok)))     !! On ne peut pas utiliser ce facteur..., car trop important ou changement de sens...
                          fact = fact / factor_div_it
                          IF ( fact .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                             fact = 2*dia(i)
                             dia_ok = .true.
                          ENDIF
                       ENDDO

                       IF ( .NOT.dia_ok ) THEN
                          IF ( signe.EQ.1. .AND. init_ok  ) THEN !! Arsene 02-09-2015 - Alors on reviens au début et on diminue the factor
                             dia(i) = dia(i) - fact !! on revient à l'état préscedent ;)
                             fact = fact / factor_div_it
                             IF ( fact .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                                dia_ok = .true.
                             ENDIF
                          ELSE
                             IF (init_ok .AND. num_it.GE.1) init_ok=.FALSE.
                             signe = -1
                             dia(i) = dia(i) - fact
                          ENDIF
                       ENDIF
                    ELSEIF ( (volume1-volume) .LT. (accept_sigma_it*volume) ) THEN !! dia too low  min_stomate
                       IF ( dia(i).GT.maxdia(j) ) THEN
                          dia(i) = maxdia(j)
                          dia_ok = .true.
                       ELSEIF ( signe.EQ.-1. .AND. num_it.NE.1) THEN
                          fact = fact / factor_div_it
                          IF ( fact .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                             dia_ok = .true.
                          ENDIF
                       ENDIF
                       IF ( .NOT.dia_ok ) THEN
                          IF ( signe.EQ.-1. .AND. init_ok ) THEN !! Arsene 02-09-2015 - Alors on reviens au début et on diminue the factor
                             dia(i) = dia(i) + fact !! on revient à l'état préscedent ;)
                             fact = fact / factor_div_it
                             IF ( fact .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                                dia_ok = .true.
                             ENDIF
                          ELSE
                             IF (init_ok .AND. num_it.GE.1.) init_ok=.FALSE.
                             signe = 1.
                             dia(i) = dia(i) + fact
                          ENDIF
                       ENDIF
                    ELSE !! Good dia
                       dia_ok = .true.
                    ENDIF

!write(*,*) "num_it", num_it, "dia", dia(1)

                    num_it = num_it+1
                    IF ((num_it .GE. 100 ) .AND. (.NOT.dia_ok) ) THEN ! Si trop de boucle... limit à 100 ?1
                       dia_ok = .true.
                       write(*,*) "The iteration in lpj_establish.f90 need probably to be check (Arsene)"
                       !! Arsene 11-08-2015 - By default: Dia = last dia...
                    ENDIF

                 ENDDO  !! While loops!! Arsene 11-08-2015 Attention à tester l'itération, notemment pour ajuster the "accept_sigma"


                           vn(i) = pipe_tune_shrub1 * ((ind(i,j)+d_ind(i,j))*pi/4* &
                                   & MAX(MIN(dia(i),maxdia(j)),dia_cut(i,j))**2)**pipe_tune_shrub_exp_coeff   !! Arsene 01-09-2015 - Add dia_cut

                      ENDIF                !! Arsene 03-08-2015 - Add allometry for shrubs

                   ENDIF
                ENDDO ! Loop over # pixels - domain size
             ELSE ! for grasses, cnd=1, so the above calculation cancels
                vn(:) = ind(:,j) + d_ind(:,j)
             ENDIF

          !! 4.3.2 without DGVM (static)\n
          ELSE 
             DO i=1,npts ! Loop over # pixels - domain size
                IF( (is_tree(j) .OR. is_shrub(j)) .AND. (d_ind(i,j)+ind(i,j)).GT.min_stomate) THEN     !! Arsene 31-07-2014 modifications woodmass...

                   IF(total_bm_c(i).LE.min_stomate) THEN

                      ! new wood mass of PFT
                      woodmass_ind(i,j) = &
                           & (((biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) &
                           & + biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon))) &
                           & + (bm_sapl(j,isapabove,icarbon) + bm_sapl(j,isapbelow,icarbon) &
                           & + bm_sapl(j,iheartabove,icarbon) + bm_sapl(j,iheartbelow,icarbon))*d_ind(i,j))/(ind(i,j)+d_ind(i,j))

                   ELSE
 
                      ! new biomass is added to the labile pool, hence there is no change 
                      ! in CA associated with establishment
                      woodmass_ind(i,j) = &
                           & (biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) &
                           & + biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon)) &
                           & /(ind(i,j) + d_ind(i,j))

                   ENDIF
                ENDIF
             ENDDO ! Loop over # pixels - domain size

             vn(:) = un ! cannot change in static!, and veget_max implicit in d_ind

          ENDIF

          !! 4.4 total biomass of PFT added by establishment defined over veget_max ...
          total_bm_sapl(:,:) = zero
          DO k = 1, nparts ! Loop over # litter tissues (nparts=8); see 'stomate_constants.f90'
             WHERE(d_ind(:,j).GT.min_stomate.AND.total_bm_c(:).GT.min_stomate.AND.vn(:).GT.min_stomate)

                total_bm_sapl(:,icarbon) = total_bm_sapl(:,icarbon) + & 
                     bm_sapl(j,k,icarbon) * d_ind(:,j) / vn(:)
             ENDWHERE
          ENDDO ! Loop over # litter tissues

          !! 4.5 Update biomass at each component
          DO k = 1, nparts ! Loop over # litter tissues

             bm_new(:) = zero

             ! first ever establishment, C flows
             WHERE( d_ind(:,j).GT.min_stomate .AND. &
                  total_bm_c(:).LE.min_stomate .AND. &
                  vn(:).GT.min_stomate)
                ! WHERE ( many_new(:) )

                !bm_new(:) = d_ind(:,j) * bm_sapl(j,k) / veget_max (:,j)
                bm_new(:) = d_ind(:,j) * bm_sapl(j,k,icarbon) / vn(:)

                biomass(:,j,k,icarbon) = biomass(:,j,k,icarbon) + bm_new(:)

                co2_to_bm(:,j) = co2_to_bm(:,j) + bm_new(:) / dt

             ENDWHERE

             ! establishment into existing population, C flows
             WHERE(d_ind(:,j).GT.min_stomate.AND.total_bm_c(:).GT.min_stomate)

                bm_new(:) = total_bm_sapl(:,icarbon) * biomass(:,j,k,icarbon) / total_bm_c(:)

                biomass(:,j,k,icarbon) = biomass(:,j,k,icarbon) + bm_new(:)
                co2_to_bm(:,j) = co2_to_bm(:,j) + bm_new(:) / dt

             ENDWHERE
          ENDDO ! Loop over # litter tissues


          !! 4.6 Decrease leaf age in youngest class if new leaf biomass is higher than old one.
          WHERE ( d_ind(:,j) * bm_sapl(j,ileaf,icarbon) .GT. min_stomate )
 
             ! reset leaf ages. Should do a real calculation like in the npp routine, 
             ! but this case is rare and not worth messing around.
             ! SZ 080806, added real calculation now, because otherwise leaf_age/leaf_frac
             ! are not initialised for the calculation of vmax, and hence no growth at all.
             ! logic follows that of stomate_npp.f90, just that it's been adjusted for the code here
             leaf_age(:,j,1) = leaf_age(:,j,1) * leaf_mass_young(:) / &
                  ( leaf_mass_young(:) + d_ind(:,j) * bm_sapl(j,ileaf,icarbon) )

          ENDWHERE

          leaf_mass_young(:) = leaf_mass_young(:) + d_ind(:,j) * bm_sapl(j,ileaf,icarbon)   

          !! 4.7 Youngest class: new mass in youngest class divided by total new mass
          WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )
             ! new age class fractions (fraction in youngest class increases)
             leaf_frac(:,j,1) = leaf_mass_young(:) / biomass(:,j,ileaf,icarbon)

          ENDWHERE

          !! 4.8 Other classes: old mass in leaf age class divided by new mass
          DO m = 2, nleafages

             WHERE ( biomass(:,j,ileaf,icarbon) .GT. min_stomate )

                leaf_frac(:,j,m) = leaf_frac(:,j,m) * & 
                     ( biomass(:,j,ileaf,icarbon) + d_ind(:,j) * bm_sapl(j,ileaf,icarbon) ) /  biomass(:,j,ileaf,icarbon)

             ENDWHERE

          ENDDO

          !! 4.9 Update age and number of individuals

          WHERE ( d_ind(:,j) .GT. min_stomate )

             age(:,j) = age(:,j) * ind(:,j) / ( ind(:,j) + d_ind(:,j) )

             ind(:,j) = ind(:,j) + d_ind(:,j)
             !! Arsene 01-09-2015 - If dia_cut, it's ok to increase ind number if height it's not obtain ?
             !!                       ==> chose: Yes, because we can find land-scape with low shrubs...
          ENDWHERE

          !! 4.10 Convert excess sapwood to heartwood
          !!      No longer done : supressed by S. Zaehle given that the LPJ logic of carbon allocation was 
          !!      contradictory to SLAVE allocation. See CVS tag 1_5 for initial formulation. 


       ENDIF ! natural

    ENDDO ! Loop over # PFTs


  !! 5. history

    d_ind = d_ind / dt

    CALL xios_orchidee_send_field("IND_ESTAB",d_ind)
    CALL xios_orchidee_send_field("ESTABTREE",estab_rate_max_tree)
    CALL xios_orchidee_send_field("ESTABGRASS",estab_rate_max_grass)

    CALL histwrite_p (hist_id_stomate, 'IND_ESTAB', itime, d_ind, npts*nvm, horipft_index)
    CALL histwrite_p (hist_id_stomate, 'ESTABTREE', itime, estab_rate_max_tree, npts, hori_index)
    CALL histwrite_p (hist_id_stomate, 'ESTABGRASS', itime, estab_rate_max_grass, npts, hori_index)
!!DZADD
!    CALL histwrite_p (hist_id_stomate, 'EST_LEAF_FRAC1', itime, leaf_frac(:,:,1), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'EST_LEAF_FRAC2', itime, leaf_frac(:,:,2), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'EST_LEAF_FRAC3', itime, leaf_frac(:,:,3), npts*nvm, horipft_index)
!    CALL histwrite_p (hist_id_stomate, 'EST_LEAF_FRAC4', itime, leaf_frac(:,:,4), npts*nvm, horipft_index)
!!ENDDZADD

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving establish'

  END SUBROUTINE establish

END MODULE lpj_establish
