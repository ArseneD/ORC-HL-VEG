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
                        biomass, leaf_frac, ind, cn_ind, co2_to_bm, height)   !! Arsene 04-09-2014 - add height

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

    !! 0.4 Local variables

    REAL(r_std), DIMENSION(npts)                              :: dia             !! stem diameter (m)
    REAL(r_std), DIMENSION(npts)                              :: woodmass        !! woodmass (gC/(m^2 of ground))
    REAL(r_std), DIMENSION(npts)                              :: woodmass_ind    !! woodmass of an individual (gC)
    INTEGER(i_std)                                            :: i,j             !! index (unitless)
    REAL(r_std)                                               :: height2, woodmass2, ratio         !! local height (m) !! Arsene 04-09-2014
!    REAL(r_std), DIMENSION(nvm)                               :: maxdia2, mindia !! Arsene 29-10-2014 - !! Arsene 11-08-2015 - Remove...

    REAL(r_std)                                   :: signe, signe_presc, factor, num_it, accept_sigma, volume1, volume2, dia2 !! Arsene 03-08-2015 - Add for iteration
    LOGICAL                                                   :: ind_ok !! Arsene 03-08-2015 - Add for iteration
    REAL(r_std)                                       :: pt1, pt2, pt3, ptcoeff, pdensity  !! Arsene 03-08-2015 - Change pipe_tune for shrubs

!_ ================================================================================================================================

    DO j = 2,nvm

      ! only when the DGVM is not activated or agricultural PFT.

                IF ( is_tree(j) ) THEN                        !! Arsene 03-08-2015 - Change pipe_tune for shrubs
                   pt1 = pipe_tune1
                   pt2 = pipe_tune2
                   pt3 = pipe_tune3
                   ptcoeff = pipe_tune_exp_coeff
                   pdensity = pipe_density
                ELSE
                   pt1 = pipe_tune1_for_shrub
                   pt2 = pipe_tune2_for_shrub
                   pt3 = pipe_tune3_for_shrub
                   ptcoeff = pipe_tune_exp_coeff_for_shrub
                   pdensity = pipe_density_shrub
                ENDIF                                         !! Arsene 03-08-2015 - Change pipe_tune for shrubs

      IF ( ( .NOT. control%ok_dgvm .AND. lpj_gap_const_mort ) .OR. ( .NOT. natural(j) ) ) THEN

        !
        !! 1.Update crown area
        !

        cn_ind(:,j) = zero

        IF ( is_tree(j) .OR. is_shrub(j) ) THEN        !! Arsene 31-07-2014 modifications 

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
!write(*,*) "PFT n° ", j, "woodmass avant", woodmass(i), "ind", ind(i,j)

!                IF ( is_tree(j) ) THEN                        !! Arsene 03-08-2015 - Change pipe_tune for shrubs
!                   pt1 = pipe_tune1
!                   pt2 = pipe_tune2
!                   pt3 = pipe_tune3
!                   ptcoeff = pipe_tune_exp_coeff
!                   pdensity = pipe_density
!                ELSE
!                   pt1 = pipe_tune1_for_shrub
!                   pt2 = pipe_tune2_for_shrub
!                   pt3 = pipe_tune3_for_shrub
!                   ptcoeff = pipe_tune_exp_coeff_for_shrub
!                   pdensity = pipe_density_shrub
!                ENDIF                                         !! Arsene 03-08-2015 - Change pipe_tune for shrubs

                ratio=2.
                IF ( is_shrub(j) .AND. ind(i,j).GT.min_stomate .AND. shrubs_like_trees ) THEN        !! 03-08-2015 Rajout de "shrubs_like_trees"

                !
                !! 1. We have to compare the height save and the height from number of ind ==> height2 
                !

                  height2 = pt2 * (veget_max(i,j) / (ind(i,j) * pt1)) **(pt3/ptcoeff)

                  !! Intermediate Calcul
                  !cn_ind(i,j) = veget_max(i,j) / ind(i,j)
                  !pipe_tune1 * dia(i) ** pipe_tune_exp_coeff = veget_max(i,j) / ind(i,j)
                  !dia(i) = (veget_max(i,j) / (ind(i,j) * pipe_tune1))**(1/pipe_tune_exp_coeff)

                  !With relation:
                  ! (1) ind(i,j) = veget_max(i,j) / cn_ind(i,j) WITHOUT DGVM veget_max=cst
                  ! (2) cn_ind(i,j) = pipe_tune1 * dia(i) ** pipe_tune_exp_coeff
                  ! (3) height2 = pipe_tune2 * dia(i)**(pipe_tune3)

                !
                !! 2. If the relation is not true, that mean the shrub was cut
                !
                !!      Care ! Can also be if height > height_presc(j) or if height < height min (from mindia)
                !!           => for the second point, we need height min (from mindia) < height min (cut on snow) -> stomate_lpj

                  IF ((height2-height(i,j)) .GT. 0.001 .AND. height(i,j).LT.height_presc(j)  &
                              .AND. height(i,j).GT.(height_presc(j)/fact_min_height) ) THEN !! .AND. height(i,j).GT.minheight) ??
!write(*,*) "On modifie la hauteur du buisson"
                  !! 3. In this case the shrub maintain the basis of its cylinder, but the lost in only on height
                  !!      That mean: ind=cst & cn_ind=cst & dia=cst. We have ind ==> we can recalculate the both other
                  !

                    cn_ind(i,j) = veget_max(i,j) / ind(i,j)
                    dia = (cn_ind(i,j)/pt1)**(1/ptcoeff)

                  !
                  !! 4. With a homogeneous vertical density on shrub, the lost of biomass is a linear fonction with height
                  !!      We can calcul the height with simple comparison between theorical woodmass (with the dia) and the woodmass
                  !
                    woodmass2 = ind(i,j) * pdensity * pi/4. * pt2 * &
                       & (height2/pt2)**((2.+pt3)/ pt3)

                  !!      After we deduce ce ratio between bot and we calculate the new height
                    ratio=woodmass(i)/woodmass2
                    height(i,j)=ratio*height2
!write(*,*) "ratio", ratio
                  ENDIF
                ENDIF   ! IF is_shrub

! Nouveau buissons !
                IF ( is_shrub(j) .AND. ratio .GE. 1. .AND. .NOT.shrubs_like_trees) THEN

                   !
                   !! 5. Calculate the provisionnal vegetation area
                   !

                   !! 5.1 Calculate stem diameter per individual shrub (dia)        !! Arsene 22-07-2015
                   ! Equation from Aiba & Kohyama (1996), Constant value from Martinez et Lopez (2003)
                   ! No need "max diameter", because is already take into account in the equation
                   ! dia(i) = ( 1 / (pipe_tune_shrub2 * ( 1/height(i,j) - 1/height_presc(j) )))**(1/pipe_tune_shrub3) /100

!                   mindia(j)= ( 1 / (pipe_tune_shrub2 * &
!                      & ( 1/(height_presc(j)/fact_min_height) - 1/height_presc(j) )))**(1/pipe_tune_shrub3) / 100
!min dia ok ==> via min height
!write(*,*) "ind initial = ", ind(i,j), "woodmass(i)", woodmass(i), "veget_max(i,j)", veget_max(i,j)

              !! On définit maxdia à 95% de max height (car fonction infinie, ce qui peut donner des chiffres impossibles...
!              maxdia2(j) = ( height_presc(j) * 0.93 / (pipe_tune_shrub2*0.07) )**(1/pipe_tune_shrub3) / 100 
!              mindia(j) = ( height_presc(j) / (pipe_tune_shrub2 * (fact_min_height-1)) )**(1/pipe_tune_shrub3) / 100
!write(*,*) "maxdia2", maxdia2(j), "height", height_presc(j)*pipe_tune_shrub2*100**pipe_tune_shrub3*maxdia2(j)**pipe_tune_shrub3 &
!                         & / ( height_presc(j) + pipe_tune_shrub2 * 100**pipe_tune_shrub3 * maxdia2(j)**pipe_tune_shrub3 ), &
!          & "mindia(j)", mindia(j), height_presc(j) * pipe_tune_shrub2 * 100**pipe_tune_shrub3 * mindia(j)**pipe_tune_shrub3 &
!                         & / ( height_presc(j) + pipe_tune_shrub2 * 100**pipe_tune_shrub3 * mindia(j)**pipe_tune_shrub3 )
               
              signe = 0
              signe_presc = 0
              factor = ind(i,j)*0.1 

              ind_ok = .false.
              num_it = 0

!              accept_sigma = woodmass(i) / ( pipe_density_shrub * 0.000001 )
              accept_sigma = 0.002
!write(*,*) "accept_sigma", accept_sigma
              DO WHILE ( .NOT.ind_ok ) ! ( ind_ok .EQV. .false. )

!                   signe_presc = signe

! two way to calculate volume:
! 1. by the woodmass, density and ind number
                   volume1 = woodmass(i) / ( pipe_density_shrub * ind(i,j) )

! 2. by the diameter... (calcul by ind number & crown area= veget_max)
                   dia2 = ( veget_max(i,j) / pipe_tune_shrub1 )**(1/(2*pipe_tune_shrub_exp_coeff)) / (pi/4 * ind(i,j))**0.5 ! pipe_tune_shrub1 = beta ( et non log beta)==> beta = 10**(log(beta)) 
                   volume2 = pi/4 * height_presc(j) * pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia2**(pipe_tune_shrub3+2) &
                                & / ( height_presc(j) + pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia2**pipe_tune_shrub3 )

               !    volume(i) = woodmass(i) * 10**(pipe_tune_shrub1) * dia(i)**pipe_tune_shrub_exp_coeff / (veget_max(i,j) * pipe_density)  ==> Pour les arbres !
!write(*,*) "pipe_tune_shrub1", pipe_tune_shrub1, "pipe_tune_shrub_exp_coeff", pipe_tune_shrub_exp_coeff, "ind(i,j)", ind(i,j)

!write(*,*) "première partie = ", ( veget_max(i,j) / pipe_tune_shrub1 )**(1/(2*pipe_tune_shrub_exp_coeff)) 
!write(*,*) "deuxième  partie = ", (pi/4 * ind(i,j))**(1/2), "pi/4 =", pi/4, "ind =", ind(i,j), "esn:", (pi/4 * ind(i,j))
!write(*,*) "3**(1/2)", 3**(1/2),  (pi/4 * ind(i,j)), (pi/4 * ind(i,j))**(1/2), (1/2), 0.5, (1./2.) 

!write(*,*) "volume1", volume1, "dia2", dia2,  "volume2", volume2,"vol1-vol2= ", volume2-volume1
!read(*,*)

!! Arsene ==> on pourrait faire dépendre le facteur "factor", dépendement de la première différence de volume... mais attention ^^
!! TEST :
                   IF ( ABS(volume2-volume1).GT.0.001 .AND. (num_it.EQ.0)) THEN
                       factor = ind(i,j)
                   ELSEIF ( ABS(volume2-volume1).GT.0.0001 .AND. (num_it.EQ.0)) THEN
                       factor = ind(i,j)*0.1
                   ELSEIF ( num_it.EQ.0 ) THEN !IF ( ABS(volume2-volume1).GT.0.00001 ) THEN
                       factor = ind(i,j)*0.01
                   ENDIF

                   IF  ( (volume2-volume1) .GT. min_stomate ) THEN ! Si les 2 volume sont 'positivement' different, par rapport à une valeur définie...
                        ! Alors cela signifie que ind(i,j) utilisee < ind(i,j) réel
                        ind(i,j) = ind(i,j) + factor
                        signe = 1
!write(*,*) "(volume2-volume1) .GT. accept_sigma"
                   ELSEIF ( (volume2-volume1).LT.-min_stomate .AND. (ind(i,j) - factor).LT.min_stomate ) THEN
 !                      write(*,*) "on est ici..."
                       signe = -signe_presc 
                   ELSEIF ( (volume2-volume1) .LT. -min_stomate ) THEN
                        ind(i,j) = ind(i,j) - factor
                        signe = -1
!write(*,*) "(volume2-volume1) .LT. accept_sigma )"
                   ELSE !! Alors on est < accept sigma = > on à la bonne valeur de "ind"
                        ind_ok = .true.
!write(*,*) "on a trouvé le bon nombre d'individus"
                   ENDIF

                   IF ( ((signe + signe_presc) .EQ. 0) .AND. .NOT.ind_ok ) THEN
                        factor = factor / 10
!write(*,*) "Dimin ution of factor..."
                        IF ( factor.LT.(accept_sigma*ind(i,j)) ) THEN ! Precision de ind...
                            ind_ok = .true.
!write(*,*) "########################################################################################################"
                        ENDIF
                   ELSEIF ( signe.EQ.-1 .AND. dia2.GE.maxdia(j) ) THEN
                        ind_ok = .true.  !! Dia max est atteint
!                        dia(i) = maxdia2(j)@
                        ind(i,j) = (veget_max(i,j)/pipe_tune_shrub1)**(1/pipe_tune_shrub_exp_coeff) / &
                                          & (pi/4*maxdia(j)**2)
!write(*,*) "max_dia atteind", "ind =", ind(i,j)
!read(*,*)
                   ELSEIF ( signe.EQ.1 .AND. dia2.LE.mindia(j) ) THEN
                        ind_ok = .true.
!                        dia(i) = mindia(j)
                        ind(i,j) = (veget_max(i,j)/pipe_tune_shrub1)**(1/pipe_tune_shrub_exp_coeff) / &
                                          & (pi/4*mindia(j)**2)
!write(*,*) "min_dia atteind", "ind =", ind(i,j) 
!read(*,*)          
                   ELSEIF (.NOT.ind_ok .AND. signe_presc.EQ.0) THEN
                        signe_presc = signe
                   !! Pas besoin de redéfinir le signe prescedent car on revint tj au dernier "bon"
!                   ELSEIF (.NOT.ind_ok .AND. dia2.GT.maxdia(j) THEN
                        
                   ENDIF
                        

                   
                   num_it = num_it+1
!                   write(*,*) "number of iteration = " , num_it, "ind = ", ind(i,j)
                   IF ((num_it .GE. 100 ) .AND. (.NOT.ind_ok) ) THEN ! Si trop de boucle... limit à 100 ?
                        ind_ok = .true.
                        dia(i) = (maxdia(j)+mindia(j))/2 !! Au pif...
!                        write(*,*) "iteration shrubs not efficient"
!read(*,*)
                   ENDIF

              ENDDO !! FIN DE BOUCLE WHILE

!read(*,*)
                    !! Une fois la boucle achevée, on a notre "ind"
                    !! On recalcule donc tous les termes
                    dia(i) = ( veget_max(i,j) / pipe_tune_shrub1 )**(1/(2*pipe_tune_shrub_exp_coeff)) / (pi/4 * ind(i,j))**0.5
                    !! ==> Penser à rajouter les min /max

                    woodmass_ind(i) = woodmass(i) / ind(i,j)

                    cn_ind(i,j) = pipe_tune_shrub1 * (ind(i,j)*pi/4*dia(i)**2)**pipe_tune_shrub_exp_coeff / ind(i,j)
!                   cn_ind(i,j) = veget_max(i,j) / ind

                    height(i,j) = height_presc(j) * pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia(i)**pipe_tune_shrub3 &
                         & / ( height_presc(j) + pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia(i)**pipe_tune_shrub3 )

!write(*,*) "couronne", cn_ind(i,j) * ind(i,j), cn_ind(i,j), "dia", dia(i), "height", height(i,j)

!read(*,*)



!                   dia(i) = MAX( ( woodmass(i) / (veget_max(i,j) * pipe_density*pi/4. &
!                         & *pipe_tune2/pipe_tune1) )**(1/(2.+pipe_tune3-pipe_tune_exp_coeff)), mindia(j))

!write(*,*) "PFT n°",j,"  dia0", dia(i), "maxdia", maxdia2(j), "mindia", mindia(j)
                   !! Because 1. ind(i,j) = woodmass(i) / ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) ) 
                   !!         2. cn_ind(i,j) * ind(i,j) =  veget_max(i,j) because without DGVM
                   !!         3. cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff avec MIN( maxdia(j), dia(i)) = dia(i)
                   !!         4. simplification
                   !!         5. dia(i) = MAX( (MIN( maxdia(j), dia(i))), mindia )

                   !! 5.2 So, 1. give:
!                   ind(i,j) = woodmass(i) / &
!                           ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) )

                   !! 5.3 Individual biomass corresponding to this critical density of individuals

!                   woodmass_ind(i) = woodmass(i) / ind(i,j)
                   !So : woodmass_ind(i) =  ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) )


                   !! 5.4 Calculate provisional tree crown area for per individual tree

                   ! equation: CrownArea=pipe_tune1 * Diameter^{1.6}
!                   cn_ind(i,j) = pipe_tune1 * dia(i)** pipe_tune_exp_coeff

                   !! 5.5 Calculate height from dia
!                   height(i,j)=pipe_tune2 * dia(i)**(pipe_tune3)
!write(*,*) "ind", ind(i,j),"height", height(i,j)



                   !
                   !! 6. If total tree crown area for this tree PFT exceeds its veget_max, tree density is recalculated.
                   !

!                   IF ( cn_ind(i,j) * ind(i,j) .GT. 1.002* veget_max(i,j) ) THEN
!                       ind(i,j) = veget_max(i,j) / cn_ind(i,j)
!                       woodmass_ind(i) = woodmass(i) / ind(i,j)
!                   ELSEIF ( cn_ind(i,j) * ind(i,j) .LT. 1.002* veget_max(i,j) ) THEN
!                       cn_ind(i,j) = veget_max(i,j) / ind(i,j)
!                   ENDIF



                ELSEIF ( .NOT.is_shrub(j) .OR. shrubs_like_trees ) THEN     !! if no biomass & height (no ind) loss due to frost... (shrubs only)
                   !
                   !! 5. Calculate the provisionnal vegetation area
                   !

                   !! 5.1 Calculate stem diameter per individual tree (dia)        !! Arsene 03-09-2014
                   ! With min and max dia (to don't have no realistique value when DGVM is desactivate)
                   !dia(i) = ( height(i,j) / pipe_tune2 ) **(1/ pipe_tune3 ) !! Arsene 29-10-2014
!                   maxdia2(j)= (height_presc(j)/(pt2))**(1/ pt3)
!                   mindia(j)= (height_presc(j)/(fact_min_height*pt2))**(1/ pt3)       !! Arsene 02-04-2015 Pour le moment, height_max=height_presc & height_min = height_presc/fact_min_height


!#Ajout pour test
!                   dia(i) =  ( woodmass(i) / (veget_max(i,j) * pdensity*pi/4. &
!                         & *pt2/pt1) )**(1/(2.+pt3-ptcoeff))
!                   ind(i,j) = woodmass(i) / &
!                           ( pdensity*pi/4.*pt2 * dia(i)**(2.+pt3) )
!height(i,j)=pt2 * dia(i)**(pt3)
!write(*,*) "PFT n°", j, "dia1", dia(i), "ind", ind(i,j), "height", height(i,j) 
!#Ajout pour test

                   dia(i) = MAX(MIN( maxdia(j), ( woodmass(i) / (veget_max(i,j) * pdensity*pi/4. &
                         & *pt2/pt1) )**(1/(2.+pt3-ptcoeff)) ),mindia(j))

!write(*,*) "PFT n°",j,"  dia2", dia(i), "maxdia", maxdia(j), "mindia", mindia(j)
                   !! Because 1. ind(i,j) = woodmass(i) / ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) ) 
                   !!         2. cn_ind(i,j) * ind(i,j) =  veget_max(i,j) because without DGVM
                   !!         3. cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff avec MIN( maxdia(j), dia(i)) = dia(i)
                   !!         4. simplification
                   !!         5. dia(i) = MAX( (MIN( maxdia(j), dia(i))), mindia )

                   !! 5.2 So, 1. give:
                   ind(i,j) = woodmass(i) / &
                           ( pdensity*pi/4.*pt2 * dia(i)**(2.+pt3) )

                   !! 5.3 Individual biomass corresponding to this critical density of individuals

                   woodmass_ind(i) = woodmass(i) / ind(i,j)
                   !So : woodmass_ind(i) =  ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) )


                   !! 5.4 Calculate provisional tree crown area for per individual tree

                   ! equation: CrownArea=pipe_tune1 * Diameter^{1.6}
                   cn_ind(i,j) = pt1 * dia(i)** ptcoeff

                   !! 5.5 Calculate height from dia
                   height(i,j)=pt2 * dia(i)**(pt3)
!write(*,*) "ind", ind(i,j),"height", height(i,j)



                   !
                   !! 6. If total tree crown area for this tree PFT exceeds its veget_max, tree density is recalculated.
                   !
  
                   IF ( cn_ind(i,j) * ind(i,j) .GT. 1.002* veget_max(i,j) ) THEN
                       ind(i,j) = veget_max(i,j) / cn_ind(i,j)
                       woodmass_ind(i) = woodmass(i) / ind(i,j)
                   ELSEIF ( cn_ind(i,j) * ind(i,j) .LT. 1.002* veget_max(i,j) ) THEN
                       cn_ind(i,j) = veget_max(i,j) / ind(i,j)
                   ENDIF

                 ENDIF  !! Arsene 30-10-2014

!! Arsene 30-10-2014 -END- New Version... Because the version before was tooooooooo full of bugs
!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*) "PFT n° ", j ," woodmass(i)", woodmass(i)
              ELSE !woodmas=0  => impose some value

!!                dia(:) = maxdia(j)                                            !! Arsene 20-05-2015 Remove
!! mindia(j)= (height_presc(j)/(fact_min_height*pipe_tune2))**(1/ pipe_tune3)   !! Arsene 20-05-2015 A other solution ? 

                IF ( is_tree(j) .OR. shrubs_like_trees ) THEN  !! 05-08-2015 FOr shrub Add
                    cn_ind(i,j) = pt1 * maxdia(j)** ptcoeff      !! Arsene 20-05-2015 Use mindia or maxdia2 ???
                ELSE                    !! 05-08-2015 FOr shrub Add
                    cn_ind(i,j) =  pipe_tune_shrub1 * (pi/4*maxdia(j)**2)**pipe_tune_shrub_exp_coeff  !! Arsene 12-08-2015 With ind=1
!cn_ind(i,j) = 67. ! 5.     !! 05-08-2015 FOr shrub Add ==> à vérifier comment caler ça... Donnees de Martinez et lopezz ==> Maxdia ~= 0.18 ==> Crown area = 67. Marche avec 5..
                ENDIF                   !! 05-08-2015 FOr shrub Add
!write(*,*) "pft", j, "cn_ind(i,j)", cn_ind(i,j)
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
!write(*,*) "pft=", j, "ind 00002 = ", ind(1,j)        
        WHERE ( veget_max(:,j) .GT. zero )

          ind(:,j) = veget_max(:,j) / cn_ind(:,j)  

        ELSEWHERE

          ind(:,j) = zero

        ENDWHERE
!write(*,*) "pft=", j, "ind 00003 = ", ind(1,j)
      ENDIF     ! IF ( ( .NOT. control%ok_dgvm .AND. lpj_gap_const_mort ) .OR. ( .NOT. natural(j) ) )

    ENDDO       ! loop over PFTs
!read(*,*) !! Arsene 20-07-2015

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

                IF ( is_tree(j) ) THEN                        !! Arsene 03-08-2015 - Change pipe_tune for shrubs
                   pt1 = pipe_tune1
                   pt2 = pipe_tune2
                   pt3 = pipe_tune3
                   ptcoeff = pipe_tune_exp_coeff
                   pdensity = pipe_density
                ELSE
                   pt1 = pipe_tune1_for_shrub
                   pt2 = pipe_tune2_for_shrub
                   pt3 = pipe_tune3_for_shrub
                   ptcoeff = pipe_tune_exp_coeff_for_shrub
                   pdensity = pipe_density_shrub
                ENDIF                                         !! Arsene 03-08-2015 - Change pipe_tune for shrubs

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
