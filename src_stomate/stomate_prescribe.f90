! =================================================================================================================================
! MODULE       : stomate_prescribe
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Initialize and update density, crown area.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	:
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_prescribe.f90 $
!! $Date: 2013-02-06 14:46:33 +0100 (Wed, 06 Feb 2013) $
!! $Revision: 1170 $
!! \n
!_ ================================================================================================================================

MODULE stomate_prescribe

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC prescribe,prescribe_clear

    ! first call
    LOGICAL, SAVE                                              :: firstcall = .TRUE.
!$OMP THREADPRIVATE(firstcall)

CONTAINS

! =================================================================================================================================
!! SUBROUTINE   : prescribe_clear
!!
!>\BRIEF        : Set the firstcall flag back to .TRUE. to prepare for the next simulation.
!_=================================================================================================================================

  SUBROUTINE prescribe_clear
    firstcall=.TRUE.
  END SUBROUTINE prescribe_clear

!! ================================================================================================================================
!! SUBROUTINE   : prescribe
!!
!>\BRIEF         Works only with static vegetation and agricultural PFT. Initialize biomass,
!!               density, crown area in the first call and update them in the following.
!!
!! DESCRIPTION (functional, design, flags): \n
!! This module works only with static vegetation and agricultural PFT.
!! In the first call, initialize density of individuals, biomass, crown area,
!! and leaf age distribution to some reasonable value. In the following calls,
!! these variables are updated.
!!
!! To fulfill these purposes, pipe model are used:
!! \latexonly 
!!     \input{prescribe1.tex}
!!     \input{prescribe2.tex}
!! \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLES(S): ::ind, ::cn_ind, ::leaf_frac
!!
!! REFERENCES   :
!! - Krinner, G., N. Viovy, et al. (2005). "A dynamic global vegetation model 
!!   for studies of the coupled atmosphere-biosphere system." Global 
!!   Biogeochemical Cycles 19: GB1015, doi:1010.1029/2003GB002199.
!! - Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!   plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!   global vegetation model, Global Change Biology, 9, 161-185.
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE prescribe (npts, &
                        veget_max, dt, PFTpresent, everywhere, when_growthinit, &
                        biomass, leaf_frac, ind, cn_ind, co2_to_bm, height, dia_cut)   !! Arsene 04-09-2014 - add height, !! 27-08-2015 Arsene - Add dia_cut

!! 0. Parameters and variables declaration

   !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                :: npts            !! Domain size (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_max       !! "maximal" coverage fraction of a PFT (LAI -> infinity) on ground (unitless;0-1)
    REAL(r_std), INTENT(in)                                   :: dt              !! time step (dt_days)

    !! 0.2 Output variables 

    !! 0.3 Modified variables

    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent      !! PFT present (0 or 1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere      !! is the PFT everywhere in the grid box or very localized (after its introduction) (?)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit !! how many days ago was the beginning of the growing season (days)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass   !! biomass (gC/(m^2 of ground))
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac       !! fraction of leaves in leaf age class (unitless;0-1)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind             !! density of individuals (1/(m^2 of ground))
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind          !! crown area per individual (m^2)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm       !! co2 taken up by carbohydrate 
                                                                                 !! reserve at the beginning of the
                                                                                 !! growing season @tex ($gC m^{-2} 
                                                                                 !! of total ground/day$) @endtex 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: height          !! Height of vegetation (m) 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: dia_cut         !! 27-08-2015 Arsene - Add dia_cut 

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts)                              :: dia             !! stem diameter (m)
    REAL(r_std), DIMENSION(npts)                              :: woodmass        !! woodmass (gC/(m^2 of ground))
!    REAL(r_std), DIMENSION(npts)                              :: woodmass_ind    !! woodmass of an individual (gC) !! Arsene 01-09-2015 - Remove (not use)
    INTEGER(i_std)                                            :: i,j             !! index (unitless)
    REAL(r_std)                                               :: height2, woodmass2!, ratio         !! local height (m) !! Arsene 04-09-2014
!    REAL(r_std), DIMENSION(nvm)                               :: maxdia2, mindia !! Arsene 29-10-2014 - !! Arsene 11-08-2015 - Remove...

    REAL(r_std)                                   :: signe, factor, num_it, num_it2, volume1, volume2, dia2 !! Arsene 03-08-2015 - Add for iteration
    LOGICAL                                       :: ind_ok !! Arsene 03-08-2015 - Add for iteration
    REAL(r_std)                                   :: pt1, pt2, pt3, ptcoeff, pdensity  !! Arsene 03-08-2015 - Change pipe_tune for shrubs

!_ ================================================================================================================================

    DO j = 2,nvm

      ! only when the DGVM is not activated or agricultural PFT.



      IF ( ( .NOT. control%ok_dgvm .AND. lpj_gap_const_mort ) .OR. ( .NOT. natural(j) ) ) THEN

        !
        !! 1.Update crown area
        !

        cn_ind(:,j) = zero

        IF ( is_tree(j) .OR. is_shrub(j) ) THEN        !! Arsene 31-07-2014 modifications 


          IF ( is_tree(j) ) THEN                        !! Arsene 03-08-2015 - Change pipe_tune for shrubs
             pt1 = pipe_tune1
             pt2 = pipe_tune2
             pt3 = pipe_tune3
             ptcoeff = pipe_tune_exp_coeff
             pdensity = pipe_density
          ELSEIF  ( shrubs_like_trees ) THEN!! Arsene 25-08-2015 - We could use pt1 and other if shrubs like trees or not...
             pt1 = pipe_tune1_for_shrub
             pt2 = pipe_tune2_for_shrub
             pt3 = pipe_tune3_for_shrub
             ptcoeff = pipe_tune_exp_coeff_for_shrub
             pdensity = pipe_density_shrub
          ELSE !! Shrub not like tree
             pt1 = pipe_tune_shrub1
             pt2 = pipe_tune_shrub2
             pt3 = pipe_tune_shrub3
             ptcoeff = pipe_tune_shrub_exp_coeff
             pdensity = pipe_density_shrub
          ENDIF                                         !! Arsene 03-08-2015 - Change pipe_tune for shrubs

          !
          !! 1.1 treat for trees
          !

          dia(:) = zero

          DO i = 1, npts ! loop over grid points

            IF ( veget_max(i,j) .GT. zero ) THEN

              !! 1.1.1 calculate wood mass on an area basis, which include sapwood and heartwood aboveground and belowground.

              woodmass(i) = (biomass(i,j,isapabove,icarbon) + biomass(i,j,isapbelow,icarbon) + &
                   biomass(i,j,iheartabove,icarbon) + biomass(i,j,iheartbelow,icarbon)) * veget_max(i,j)  

              IF ( woodmass(i) .GT. min_stomate ) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Arsene 30-10-2014 -START- New Version... Because the version before was tooooooooo full of bugs

!! Arsene 27-08-2015 - Old version of save vegetation under snow - START
!                ratio=2.
!                IF ( is_shrub(j) .AND. ind(i,j).GT.min_stomate ) THEN

                !
                !! 1. We have to compare the height save and the height from number of ind ==> height2 
                !

!                  IF (  shrubs_like_trees ) THEN
!                      height2 = pt2 * (veget_max(i,j) / (ind(i,j) * pt1)) **(pt3/ptcoeff)

                      !! Intermediate Calcul
                      !cn_ind(i,j) = veget_max(i,j) / ind(i,j)
                      !pipe_tune1 * dia(i) ** pipe_tune_exp_coeff = veget_max(i,j) / ind(i,j)
                      !dia(i) = (veget_max(i,j) / (ind(i,j) * pipe_tune1))**(1/pipe_tune_exp_coeff)

                      !With relation:
                      ! (1) ind(i,j) = veget_max(i,j) / cn_ind(i,j) WITHOUT DGVM veget_max=cst
                      ! (2) cn_ind(i,j) = pipe_tune1 * dia(i) ** pipe_tune_exp_coeff
                      ! (3) height2 = pipe_tune2 * dia(i)**(pipe_tune3)

!                  ELSE
!                      dia2 = ( veget_max(i,j) / pt1 )**(1./(2.*ptcoeff)) / (pi/4. * ind(i,j))**0.5
!                      height2 = height_presc(j) * pt2 * 100.**pt3 * dia2**pt3 &
!                             & / ( height_presc(j) + pt2 * 100.**pt3 * dia2**pt3 )
                      !! Cf to allometry of shrub
!                  ENDIF


                !
                !! 2. If the relation is not true, that mean the shrub was cut
                !
                !!      Care ! Can also be if height > height_presc(j) or if height < height min (from mindia)
                !!           => for the second point, we need height min (from mindia) < height min (cut on snow) -> stomate_lpj

!                  IF ((height2-height(i,j)) .GT. 0.001 .AND. height(i,j).LT.height_presc(j)  &
!                              .AND. height(i,j).GT.(height_presc(j)/fact_min_height) ) THEN !! .AND. height(i,j).GT.minheight) ??

                  !! 3. Calculation of cn_ind
                  !! In this case the shrub maintain the basis of its cylinder, but the lost in only on height
                  !!      That mean: ind=cst & cn_ind=cst [& dia=cst]. We have ind ==> we can recalculate the both other
                  !
!                     cn_ind(i,j) = veget_max(i,j) / ind(i,j)
                  !
                  !! 4. With a homogeneous vertical density on shrub, the lost of biomass is a linear fonction with height
                  !!      We can re-calcul the height with simple comparison between theorical woodmass2 (from the cylinder equation) and the woodmass
                  !!      woodmass2 = ind * volume * desity
                  !

!                     IF (  shrubs_like_trees ) THEN
!!                        dia(i) = (cn_ind(i,j)/pt1)**(1/ptcoeff) !! Arsene 25-08-2015 - A priori non utilisé
!                        woodmass2 = ind(i,j) * pdensity * pi/4. * pt2 * &
!                            & (height2/pt2)**((2.+pt3)/ pt3)
!                     ELSE
!                        woodmass2 = ind(i,j) * pdensity * pi/4. * dia2**2 * height2
!                     ENDIF

                  !!      After we deduce ce ratio between bot and we calculate the new height
!                    ratio=woodmass(i)/woodmass2
!                    height(i,j)=ratio*height2
!write(*,*) "ratio", ratio
!                  ENDIF
!                ENDIF   ! IF is_shrub
!! Arsene 27-08-2015 - Old version of save vegetation under snow - END

                IF ( is_shrub(j) .AND. .NOT.shrubs_like_trees) THEN
                   !
                   !! 5. Calculate the provisionnal vegetation area
                   !

                   !! 5.1 Calculate stem diameter per individual shrub (dia)        !! Arsene 22-07-2015
                   ! Equation from Aiba & Kohyama (1996), Constant value from Martinez et Lopez (2003)
                   ! No need "max diameter", because is already take into account in the equation
                   ! dia(i) = ( 1 / (pipe_tune_shrub2 * ( 1/height(i,j) - 1/height_presc(j) )))**(1/pipe_tune_shrub3) /100

                   signe = 0
                   ind_ok = .false.
                   num_it = 0
                   DO WHILE ( .NOT.ind_ok ) ! ( ind_ok .EQV. .false. )

                      ! two way to calculate volume:
                      ! 1. by the woodmass, density and ind number
                      volume1 = woodmass(i) / ( pdensity * ind(i,j) )

                      ! 2. by the diameter... (calcul by ind number & crown area= veget_max)
                      dia(i) = ( veget_max(i,j) / pt1 )**(1/(2*ptcoeff)) / (pi/4 * ind(i,j))**0.5 ! pipe_tune_shrub1 = beta ( et non log beta)==> beta = 10**(log(beta))
                      !! Arsene 27-08-2015 - Note that if shrub is cut, initial "ind" can be "far aways" the real one...
                      volume2 = pi/4 * height_presc(j) * pt2 * 100**pt3 * dia(i)**(pt3+2) &
                                & / ( height_presc(j) + pt2 * 100**pt3 * dia(i)**pt3 )

                      !! Suivant la valeur de la différence "vol2-vol1" on détermine les facteurs initiaux (methode prescriptive)
                      IF (num_it.EQ.0) THEN
                         IF ( (ABS(volume2-volume1)/(volume2+volume1)).GT.0.5 ) THEN
                            factor = 0.3
                         ELSEIF ( (ABS(volume2-volume1)/(volume2+volume1)).GT.0.1 ) THEN
                            factor = ind(i,j)
                         ELSEIF (( ABS(volume2-volume1)/(volume2+volume1)).GT.0.01 ) THEN
                            factor = ind(i,j)*0.1
                         ELSE !IF ( ABS(volume2-volume1).GT.0.00001 ) THEN
                            factor = ind(i,j)*0.01
                         ENDIF
                      ENDIF

                    !! Suivant la valeur de "vol2-vol1" on détermine si le diamètre proposé est trop petit/grand ou pas mal
                    !! On distingue 3 cas:
                    !! Cas 1. Vol2>vol1 ==> Alors ind(i,j) utilisee < ind(i,j) réel
                      IF  ( (volume2-volume1)/(volume2+volume1) .GT. accept_sigma_it ) THEN ! Si les 2 volume sont 'positivement' different, par rapport à une valeur définie...

                          !! a. Si l'on doit augmenter ind, alors on doit diminuer dia...  On vérifie si on ne dépasse pas min_dia
                         IF ( (dia(i).LE.mindia(j) )) THEN
                            ind_ok = .true.
                            ind(i,j) = (veget_max(i,j)/pt1)**(1/ptcoeff) / (pi/4*mindia(j)**2)

                         !! a. Si on avait prescedement le cas "vol2<vol1", on revient à l'état prescendent, on diminue le facteur, et on retest ce facteur
                         ELSEIF ( signe.EQ.-1. ) THEN 
                            ind(i,j) = ind(i,j) + factor     !! On revient à l'état préscedent (avant les nouveaux calculs de volume)
                            factor = factor / factor_div_it  !! On diminue le facteur préscedement utilisé
                            IF ( factor .LT. (0.1 * accept_sigma_it * ind(i,j)) ) THEN  !! Si le facteur devient trop petit, c'est qu'on a trouvé le bon diamétre !
                               ind_ok = .true.
                            ELSE
                               ind(i,j) = ind(i,j) - factor  !! On va donc pouvoir tester le nouveau facteur
                            ENDIF

                         !! b. le nombre d'individus est trop faible, on l'augmente donc...
                         ELSE
                            IF ( signe .NE. 1. ) signe = 1.
                            ind(i,j) = ind(i,j) + factor
                         ENDIF

                      !! Cas 2. Vol2<Vol2 ==> Alors ind(i,j) utilisee > ind(i,j) réel
                      ELSEIF ( (volume2-volume1)/(volume2+volume1) .LT. accept_sigma_it ) THEN
                          
                         !! a. Si l'on doit diminuer ind, alors on doit augmenter dia.. On vérifie si on ne dépasse pas max_dia
                         IF ( dia(i).GE.maxdia(j) ) THEN
                            ind_ok = .true.  !! Dia max est atteint
                            ind(i,j) = (veget_max(i,j)/pt1)**(1/ptcoeff) / (pi/4*maxdia(j)**2)

                         !! a. Si on avait prescedement le cas "vol2>vol1", on revient à l'état prescendent, on diminue le facteur, et on retest ce facteur
                         ELSEIF ( signe.EQ.1. ) THEN !! Arsene 02-09-2015 - Alors on reviens au début et on diminue the factor
                            ind(i,j) = ind(i,j) - factor     !! On revient à l'état préscedent (avant les nouveaux calculs de volume)
                            factor = factor / factor_div_it  !! On diminue le facteur préscedement utilisé
                            IF ( factor .LT. (0.1 * accept_sigma_it * ind(i,j)) ) THEN  !! Si le facteur devient trop petit, c'est qu'on a trouvé le bon diamétre !
                               ind_ok = .true.
                            ELSE
                               ind(i,j) = ind(i,j) + factor   !! On va donc pouvoir tester le nouveau facteur
                            ENDIF

                         !! c. le nombre d'individus est trop important, on le diminue donc...
                         ELSE

                            !! b.1 On vérifie que si l'on diminue ind par le facteur, le résultat n'est pas négatif ou nul
                            num_it2 = 0
                            DO WHILE ( (ind(i,j)-factor).LT.min_stomate .AND. (num_it2.LT.5.) )  !! On ne peut pas utiliser ce facteur..., car trop important ou changement de sens...
                               factor = factor / factor_div_it
                               IF ( factor .LT. (0.1 * accept_sigma_it * ind(i,j)) ) THEN
                                  num_it2 = 10.
                                  ind_ok = .true.
                               ENDIF
                               num_it2 = num_it2 + 1.
                               IF ( num_it2.EQ.5. ) write(*,*) "small iteration in stomate_prescribe.f90 need to be check (Arsene)"
                            ENDDO

                            !! b.2 Après vérification du facteur, la valeur du diamètre est trop grande... on la diminue donc...
                            IF ( .NOT.ind_ok ) THEN
                               IF ( signe .NE. -1. ) signe = -1.
                               ind(i,j) = ind(i,j) - factor
                            ENDIF
                         ENDIF

                      !! Cas 3. Vol2~=Vol1 ==> On s'arrete là !
                      ELSE !! Good ind number
                         ind_ok = .true.
                      ENDIF

                      num_it = num_it+1
                      IF ((num_it .GE. 100 ) .AND. (.NOT.ind_ok) ) THEN ! Si trop de boucle... limit à 100 ?
                         ind_ok = .true.
!                         dia(i) = (maxdia(j)+mindia(j))/2 !! Au pif... Ou alors on prend le dernier dia calculer... le plus proche possible du bon résultat...
                         write(*,*) "The iteration in stomate_prescribe.f90 need probably to be check (Arsene)"
                      ENDIF
                   ENDDO !! FIN DE BOUCLE WHILE
!write(*,*) "stomate_presc", num_it

                   !! Une fois la boucle achevée, on a notre "ind"
                   !! On recalcule donc tous les termes
                   dia(i) = MIN(MAX(mindia(j),( veget_max(i,j) / pt1 )**(1/(2*ptcoeff)) / (pi/4 * ind(i,j))**0.5),maxdia(j))

                   height(i,j) = height_presc(j) * pt2 * 100**pt3 * dia(i)**pt3 &
                         & / ( height_presc(j) + pt2 * 100**pt3 * dia(i)**pt3 )

                   !! On vérifie que le buisson n'a pas été coupé... Et s'il a été coupé on calcul la hauteur...
                   IF ( dia(i).LE.dia_cut(i,j) .AND. dia_cut(i,j).GT.min_stomate ) THEN 
                        !! On fix le diametre (identique à celui qu'il y avait lors de la coupe
                        dia(i) = dia_cut(i,j)
                        !! On en déduit me nombre d'individue (car cov = veget_max)
                        ind(i,j) =  (veget_max(i,j)/pt1)**(1/ptcoeff) / (pi/4*dia(i)**2)
                        !! Calcul de la hauteur théorique (pour le diametre dia_cut)
                        height2 = height_presc(j) * pt2 * 100**pt3 * dia(i)**pt3 &
                            & / ( height_presc(j) + pt2 * 100**pt3 * dia(i)**pt3 )  
                        !! On Déduit la woodmasse théorique d'un tel cylindre (fonction de densité, height, dia)
                        woodmass2 = ind(i,j) * pdensity * pi/4. * dia(i)**2 * height2
                        !! On peut faire alors le ratio du woodmass théorique et réel
                        height(i,j)=MIN(woodmass(i)/woodmass2*height2,height_presc(j)) !! 21-09-2015 ==>0 pb, mais au cas où...
                    ELSEIF ( dia_cut(i,j).GT.min_stomate ) THEN
                        dia_cut(i,j) = zero
                    ENDIF

                    cn_ind(i,j) = pt1 * (ind(i,j)*pi/4*dia(i)**2)**ptcoeff / ind(i,j)

                ELSE!IF ( ratio .GE. 1.) THEN     !! if no biomass & height (no ind) loss due to frost... (shrubs only) => (.NOT.is_shrub(j) .OR. shrubs_like_trees)
                   !
                   !! 5. Calculate the provisionnal vegetation area
                   !

                   !! 5.1 Calculate stem diameter per individual tree (dia)        !! Arsene 03-09-2014
                   ! With min and max dia (to don't have no realistique value when DGVM is desactivate)

                   dia(i) = MAX(MIN( maxdia(j), ( woodmass(i) / (veget_max(i,j) * pdensity*pi/4. &
                         & *pt2/pt1) )**(1/(2.+pt3-ptcoeff)) ),mindia(j))

                   !! Because 1. ind(i,j) = woodmass(i) / ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) ) 
                   !!         2. cn_ind(i,j) * ind(i,j) =  veget_max(i,j) because without DGVM
                   !!         3. cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff avec MIN( maxdia(j), dia(i)) = dia(i)
                   !!         4. simplification
                   !!         5. dia(i) = MAX( (MIN( maxdia(j), dia(i))), mindia )

                   !! 5.2 So, 1. give: !! We can probably use 3.
                   ind(i,j) = woodmass(i) / &
                           ( pdensity*pi/4.*pt2 * dia(i)**(2.+pt3) )

                   !! 5.4 Calculate provisional tree crown area for per individual tree

                   ! equation: CrownArea=pipe_tune1 * Diameter^{1.6}
                   cn_ind(i,j) = pt1 * dia(i)** ptcoeff

                   !! 5.5 Calculate height from dia
                   height(i,j)=pt2 * dia(i)**(pt3)

                   !! On vérifie que le buisson n'a pas été coupé... Et s'il a été coupé on calcul la hauteur...
                   IF ( dia(i).LE.dia_cut(i,j) .AND. is_shrub(j)) THEN
                       !! On fix le diametre (identique à celui qu'il y avait lors de la coupe
                       dia(i) = dia_cut(i,j)
                       !! On peut donc déduire le nombre d'individus:
                       cn_ind(i,j) = pt1 * dia(i)** ptcoeff
                       ind(i,j) = veget_max(i,j) / cn_ind(i,j)
                       woodmass2 = ind(i,j) * pdensity * pi/4. * pt2 * dia(i)**(2+pt3)
                       !! On peut faire alors le ratio du woodmass théorique et réel
                       height(i,j)=woodmass(i)/woodmass2* pt2 * dia(i)**pt3
                   ELSEIF ( dia_cut(i,j).GT.min_stomate ) THEN
                       dia_cut(i,j) = zero
                   ENDIF                        
                        
                   !! 5.3 Individual biomass corresponding to this critical density of individuals

                   !So : woodmass_ind(i) =  ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) )

                   !
                   !! 6. If total tree crown area for this tree PFT exceeds its veget_max, tree density is recalculated.
                   !! Arsene 27-08-2015: If all good... only if maxdia or mindia use
  
                   IF ( cn_ind(i,j) * ind(i,j) .GT. 1.002* veget_max(i,j) ) THEN
                       ind(i,j) = veget_max(i,j) / cn_ind(i,j)
                   ELSEIF ( cn_ind(i,j) * ind(i,j) .LT. 1.002* veget_max(i,j) ) THEN
                       cn_ind(i,j) = veget_max(i,j) / ind(i,j)
                   ENDIF

                 ENDIF  !! Arsene 30-10-2014

!! Arsene 30-10-2014 -END- New Version... Because the version before was tooooooooo full of bugs
!!!!!!!!!!!!!!!!!!!!!!!!!!

              ELSE !woodmas=0  => impose some value

                IF ( is_tree(j) .OR. shrubs_like_trees ) THEN  !! 05-08-2015 FOr shrub Add
                    cn_ind(i,j) = pt1 * maxdia(j)** ptcoeff      !! Arsene 20-05-2015 Use mindia or maxdia2 ???
                ELSE                    !! 05-08-2015 FOr shrub Add
                    cn_ind(i,j) =  pt1 * (pi/4*maxdia(j)**2)**ptcoeff  !! Arsene 12-08-2015 With ind=1
                ENDIF                   !! 05-08-2015 FOr shrub Add
              ENDIF ! IF ( woodmass(i) .GT. min_stomate )

            ENDIF    ! veget_max .GT. 0.

          ENDDO      ! loop over grid points

        ELSE !grass

          !
          !! 1.2 grasses: crown area always set to 1m**2
          !

          WHERE ( veget_max(:,j) .GT. zero )
            cn_ind(:,j) = un
          ENDWHERE

        ENDIF   !IF ( is_tree(j) )

        !
        !! 2 density of individuals
        !

        WHERE ( veget_max(:,j) .GT. zero )

          ind(:,j) = veget_max(:,j) / cn_ind(:,j)  

        ELSEWHERE

          ind(:,j) = zero

        ENDWHERE

      ENDIF     ! IF ( ( .NOT. control%ok_dgvm .AND. lpj_gap_const_mort ) .OR. ( .NOT. natural(j) ) )

    ENDDO       ! loop over PFTs


    !
    !!? it's better to move the code for first call at the beginning of the module.
    !! 2 If it's the first call for this module, 
    !

    IF (( firstcall ) .AND. (TRIM(stom_restname_in) == 'NONE')) THEN

      WRITE(numout,*) 'prescribe:'

      ! impose some biomass if zero and PFT prescribed

      WRITE(numout,*) '   > Imposing initial biomass for prescribed trees, '// &
                      'initial reserve mass for prescribed grasses.'
      WRITE(numout,*) '   > Declaring prescribed PFTs present.'

      DO j = 2,nvm ! loop over PFTs

        DO i = 1, npts ! loop over grid points

          ! is vegetation static or PFT agricultural?
          ! Static vegetation or agricultural PFT
          IF ( ( .NOT. control%ok_dgvm ) .OR. &
               ( ( .NOT. natural(j) ) .AND. ( veget_max(i,j) .GT. min_stomate ) ) ) THEN

            !
            !! 2.1 if tree biomass is extremely small, prescribe the biomass by assuming they have sapling biomass, which is a constant in the model.
            !!     then set all the leaf age as 1.
            !
            ! if tree PFT and biomass too small, prescribe the biomass to a value.
            IF ( ( is_tree(j) .OR. is_shrub(j)) .AND. &           !! Arsene 31-07-2014 modifications - A verif
                 ( veget_max(i,j) .GT. min_stomate ) .AND. &
                 ( SUM( biomass(i,j,:,icarbon) ) .LE. min_stomate ) ) THEN
               !!? here the code is redundant, as "veget_max(i,j) .GT. min_stomate" is already met in the above if condition.

               IF (veget_max(i,j) .GT. min_stomate) THEN
                  biomass(i,j,:,:) = (bm_sapl_rescale * bm_sapl(j,:,:) * ind(i,j)) / veget_max(i,j)
               ELSE
                  biomass(i,j,:,:) = zero
               ENDIF

              ! set leaf age classes
              leaf_frac(i,j,:) = zero
              leaf_frac(i,j,1) = un

              ! set time since last beginning of growing season
              when_growthinit(i,j) = large_value

              ! seasonal trees: no leaves at beginning

              IF ( pheno_model(j) .NE. 'none' ) THEN

                biomass(i,j,ileaf,icarbon) = zero
                leaf_frac(i,j,1) = zero

              ENDIF

              co2_to_bm(i,j) = co2_to_bm(i,j) + ( SUM(biomass(i,j,:,icarbon))  / dt )
            ENDIF

            !
            !! 2.2 for grasses, set only the carbon reserve pool to "sapling" carbon reserve pool.
            !!     and set all leaf age to 1.

            IF ( (( .NOT. is_tree(j) ) .AND. (.NOT. is_shrub(j)) ) .AND. &          !! Arsene 31-07-2014 modifications
                 ( veget_max(i,j) .GT. min_stomate ) .AND. &
                 ( SUM( biomass(i,j,:,icarbon) ) .LE. min_stomate ) ) THEN

              biomass(i,j,icarbres,:) = bm_sapl(j,icarbres,:) * ind(i,j) / veget_max(i,j)

              IF ( .NOT. vascular(j)) THEN  !! Arsene 18-04-2014
                 biomass(i,j,ileaf,:) = bm_sapl(j,ileaf,:) * ind(i,j) / veget_max(i,j)
              ENDIF

              ! set leaf age classes
              leaf_frac(i,j,:) = zero
              leaf_frac(i,j,1) = un

              ! set time since last beginning of growing season
              when_growthinit(i,j) = large_value

              co2_to_bm(i,j) = co2_to_bm(i,j) + ( biomass(i,j,icarbres,icarbon)  / dt )
            ENDIF

            !
            !! 2.3 declare all PFTs with positive veget_max as present everywhere in that grid box
            !

            IF ( veget_max(i,j) .GT. min_stomate ) THEN
              PFTpresent(i,j) = .TRUE.
              everywhere(i,j) = un
            ENDIF

          ENDIF   ! not control%ok_dgvm  or agricultural

        ENDDO ! loop over grid points
      ENDDO ! loop over PFTs

      firstcall = .FALSE.

    ENDIF

  END SUBROUTINE prescribe

END MODULE stomate_prescribe
