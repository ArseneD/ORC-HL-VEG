! =================================================================================================================================
! MODULE       : lpj_crown
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Calculate individual crown area from stem mass
!!
!! \n DESCRIPTION : Calculating crown area of individual tree by diameter and tree height
!!
!! REFERENCE(S) :
!! - Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!  dynamics in the modelling of terrestrial ecosystems: comparing two
!!  contrasting approaches within European climate space,
!!  Global Ecology and Biogeography, 10, 621-637.
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_stomate/lpj_crown.f90 $
!! $Date: 2013-11-13 16:02:26 +0100 (Wed, 13 Nov 2013) $
!! $Revision: 1615 $
!! \n
!_ ================================================================================================================================

MODULE lpj_crown

  USE ioipsl_para
  USE stomate_data
  USE constantes
  USE pft_parameters
  
  IMPLICIT NONE
  
  ! private & public routines

  PRIVATE
  PUBLIC crown
  
CONTAINS
  
  
!! ================================================================================================================================
!! SUBROUTINE    : lpj_crown
!!
!>\BRIEF         Calculate individual crown area from stem mass
!!
!! DESCRIPTION   : Calculating crown area of individual tree by diameter and tree height
!! which are also calculated internally within this program from stem mass and allometory.
!! Calculations for diameter, height and crown area originate from eqns 1, 2, and 3 in 
!! Appendix B, Smith et al. (2001) following Huang et al. 1992 and Zeide 1993.
!! \latexonly
!!  \input{lpj_crown1.tex}
!!  \input{lpj_crown2.tex}
!!  \input{lpj_crown3.tex}
!! \endlatexonly
!! \n
!! where \f$k_{allom1}(=100.)\f$, \f$k_{allom2}(=40.)\f$, \f$k_{allom3}(=0.85)\f$ and \f$k_{rp}(=1.6)\f$ are 
!! constants, \f$WD\f$ is wood density (\f$=2 \times 10^5\f$ gC m\f$^3\f$) and \f$CA_{max}\f$ is maximum 
!! crown area (\f$=27.3\f$ m\f$^2\f$).
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::cn_ind (crown area per individual, @tex $m^2 $ @endtex) and ::height (m)
!!
!! REFERENCE(S)   :
!! - Huang, S., Titus, S.J. and Wiens, D.P. (1992) Comparison of nonlinear height–diameter functions for major 
!! Alberta tree species. Canadian Journal of Forest Research, 22, 1297–1304.\n
!! - Zeide, B. (1993) Primary unit of the tree crown. Ecology, 74, 1598–1602.\n
!! - Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation dynamics in the modelling of 
!! terrestrial ecosystems: comparing two contrasting approaches within European climate space,
!! Global Ecology and Biogeography, 10, 621-637.\n
!! 
!! FLOWCHART : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE crown &
       &  (npts, PFTpresent, ind, biomass, woodmass_ind, veget_max, cn_ind, height, dia_cut) !! Arsene 27-08-2015 - Add dia_cut

  !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std),INTENT(in)                         :: npts              !! Domain size (unitless) 
    LOGICAL,DIMENSION(npts,nvm),INTENT(in)            :: PFTpresent        !! Is pft there (unitless)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)        :: ind               !! [DISPENSABLE] Density of individuals 
                                                                           !! @tex $(m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm,nparts,nelements),INTENT(in) :: biomass !! [DISPENSABLE] Biomass @tex $(gC.m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,nvm),INTENT(in)        :: woodmass_ind      !! Woodmass of the individual, needed to calculate 
                                                                           !! crownarea in lpj_crown (gC)

    !! 0.2 Output variables
  
    REAL(r_std),DIMENSION(npts,nvm),INTENT(out)       :: cn_ind            !! Crown area per individual @tex $(m^{2})$ @endtex    

    !! 0.3 Modified variables

    REAL(r_std),DIMENSION(npts,nvm),INTENT(inout)     :: veget_max        !! [DISPENSABLE] "Maximal" coverage fraction of a PFT 
                                                                          !! infinity) on ground (unitless)
    REAL(r_std),DIMENSION(npts,nvm),INTENT(inout)     :: height           !! Height of vegetation (m)           
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)   :: dia_cut          !! Arsene 27-08-2015 - Add dia_cut

    !! 0.4 Local variables

!   REAL(r_std),DIMENSION(npts)                       :: woodmass        !! Wood mass of an individual (gC)
    INTEGER(i_std)                                    :: j               !! Index
    REAL(r_std),DIMENSION(npts)                       :: dia             !! Stem diameter (m)
    REAL(r_std),DIMENSION(nvm)                        :: height_presc_12 !! [DISPENSABLE] Prescribed height of each pfts (m)

    REAL(r_std)                                       :: pt1, pt2, pt3, ptcoeff, pdensity  !! Arsene 03-08-2015 - Change pipe_tune for shrubs
    REAL(r_std),DIMENSION(npts)                       :: volume, woodmass_ind2     !! Arsene 11-08-2015 - New shrub allometry (from Aiba & Kohyama, 1996)
    REAL(r_std)                                       :: signe, factor, num_it, num_it2, volume1 !! Arsene 11-08-2015 - Add for iteration -, signe_presc
    LOGICAL                                           :: dia_ok, init_ok !! Arsene 11-08-2015 - Add for iteration
    INTEGER(i_std)                                    :: i               !! Arsene 11-08-2015 - Add for iteration - index (unitless)

!_ ================================================================================================================================
    
  !! 1. Initializations
    
    !! 1.1 Check if DGVM is activated
    IF (.NOT.control%ok_dgvm .AND. lpj_gap_const_mort) THEN
       STOP 'crown: not to be called with static vegetation.'
    ENDIF
    
    !! 1.2 Initialize output to zero
    cn_ind(:,:) = zero

    !! 1.3 Copy prescribed height to new variable**3 !![DISPENSABLE]
    height_presc_12(1:nvm) = height_presc(1:nvm)     !![DISPENSABLE]
    
  !! 2. Calculate (or prescribe) crown area

    DO j = 2,nvm ! loop over PFTs
       IF (is_tree(j) .OR. is_shrub(j)) THEN     !! Arsene 31-07-2014 modifications    A VERIFFFFF ==> def cn_ind
          
          !! 2.1 Trees
          IF (natural(j)) THEN

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

             !! 2.1.1 Natural trees
             !WHERE (PFTpresent(:,j) .AND.ind(:,j).GT.min_stomate)
             IF ( is_tree(j) .OR. shrubs_like_trees ) THEN     !! Arsene 11-08-2015 - Change for shrub allometry
              WHERE (PFTpresent(:,j) .AND.woodmass_ind(:,j).GT.min_stomate)

                !! 2.1.1.1 Calculate individual wood mass**2

                !! SZ note that woodmass_ind needs to be defined on the individual, hence
                !! biomass*veget_max/ind, not as stated here, correction MERGE 
                !!         woodmass(:) = &
                !! &         (biomass(:,j,isapabove,icarbon) + biomass(:,j,isapbelow,icarbon) &
                !! &         +biomass(:,j,iheartabove,icarbon) + biomass(:,j,iheartbelow,icarbon))/ind(:,j)
          
                !! 2.1.1.2 Stem diameter from pipe model
                !          Stem diameter (pipe model) is calculated by allometory (eqn 1, Appdx B, Smith et al. (2001))
                !!!$          dia(:) = (woodmass(:)/(pipe_density*pi/4.*pipe_tune2)) &
                dia(:) = (woodmass_ind(:,j)/(pdensity*pi/4.*pt2))**(1./(2.+pt3))          !! Arsene 03-08-2015 - Change pipe_tune for shrubs

                !! 2.1.1.3 Individual tree height from pipe model
                !          Individual tree height (eqn 2, Appdx B, Smith et al. (2001))
                height(:,j) = pt2*(dia(:)**pt3)                                           !! Arsene 03-08-2015 - Change pipe_tune for shrubs

                !!!$SZ: The constraint on height has nothing to do with LPJ (for that purpose there's dia_max
                !!!$ cannot see why this is necessary - it also blurrs the output, hence I leave it commented
                !!!$ WHERE (height(:,j) > height_presc_12(j))
                !!!$    dia(:) = (height_presc_12(j)/pipe_tune2)**(1./pipe_tune3)
                !!!$    height(:,j) = height_presc_12(j)
                !!!$ ENDWHERE

                !! 2.1.1.3bis Individual tree height if shrub under snow...               !! Arsene 01-09-2015 - Add
                !! On vérifie que le buisson n'a pas été coupé... Et s'il a été coupé on REcalcul la hauteur... !! Arsene 01-09-2015
                WHERE ( dia(:).LE.dia_cut(:,j) .AND. is_shrub(j))
                   !! On fix le diametre (identique à celui qu'il y avait lors de la coupe
                   dia(:) = dia_cut(:,j)
                   !! On a déjà le nombre d'individus fixé... ATTENTION QUE CELA SOIT BIEN COMPATIBLE ==> PAS D'augmentation de ind tant que dia<dia_cut
                   !! On calcul la wodmass_ind "comme si on avait height correspondant à dia_cut
                   woodmass_ind2(:) = pdensity * pi/4. * pt2 * dia(:)**(2+pt3)
                   !! On peut donx comparer la woodmass réel par rapport à la woodmass attendu, et en déduire height !
                   height(:,j)=woodmass_ind(:,j)/woodmass_ind2(:)* pt2 * dia(:)**pt3
                ELSEWHERE ( dia_cut(:,j).GT.min_stomate )
                   dia_cut(:,j) = zero
                ENDWHERE

                !! 2.1.1.4 Crown area of individual tree  
                !          Calculate crown area, truncate crown area for trunks with large diameters 
                ! crown area cannot exceed a certain value, prescribed through maxdia 
                ! (eqn 3, Appdx B, Smith et al. (2001))
                cn_ind(:,j) = pt1*MIN(dia(:),maxdia(j))**ptcoeff                           !! Arsene 03-08-2015 - Change pipe_tune for shrubs

              ENDWHERE

             ELSE                                             !! Arsene 11-08-2015 - Change for shrub allometry
             !! First, we start by first estimation ==> Ind number is good (normalement c'est bon !)
               WHERE (PFTpresent(:,j) .AND.woodmass_ind(:,j).GT.min_stomate)
               
                 !! We start diffr-erently because is not easy to calculate Dia.
                 !! Input: woodmass_ind, ind (so woodmass), and biomass

                 !! 2.1.1.1.bis. Calculate ind volume
                 volume(:) = woodmass_ind(:,j) / pdensity
                 
                 !! After we need to calculate the good diameter... by iteration...
               ENDWHERE

               DO i = 1, npts ! loop over grid points
                IF (PFTpresent(i,j) .AND. woodmass_ind(i,j).GT.min_stomate) THEN
                 !! 2.1.1.2.bis. Stem diameter from Aiba & Kohyama
                 !          Stem diameter is calculated by allometory... but no analytique solution for:
                 !! volume(i) = pi/4 * height_presc(j) * pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia(i)**(pipe_tune_shrub3+2) &
                 !!               & / ( height_presc(j) + pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia(i)**pipe_tune_shrub3 )
                 
                 !! On part de la hauteur initiale, et on essaye de trouver le bon diamètre...
                 IF ( height(i,j) .LT. height_presc(j) ) THEN  !! Arsene 11-08-2015 Garde fou
                     dia(i) = (height(i,j)*height_presc(j) / (pt2*(height_presc(j)-height(i,j))) ) &
                                 & **(1/pt3) /100
                 ELSE
                     dia(i) = maxdia(j)
                 ENDIF

                 init_ok = .FALSE.
                 num_it = 0
                 signe = 0
                 dia_ok = .false.
                 DO WHILE ( .NOT.dia_ok )
                    IF ( num_it .EQ. 1) init_ok=.TRUE.
                    volume1 = pi/4 * height_presc(j) * pt2 * 100**pt3 * dia(i)**(pt3+2.) &
                                & / ( height_presc(j) + pt2 * 100**pt3 * dia(i)**pt3 )

                    IF ( ABS(volume1-volume(i))/volume(i).GT.1. .AND. (num_it.EQ.0)) THEN  !! We know the good volume volume(i)
                       factor = 0.3
                    ELSEIF ( ABS(volume1-volume(i))/volume(i).GT.0.5 .AND. (num_it.EQ.0)) THEN
                       factor = dia(i)
                    ELSEIF ( ABS(volume1-volume(i))/volume(i).GT.0.1 .AND. (num_it.EQ.0)) THEN
                       factor = dia(i)*0.1
                    ELSEIF ( num_it.EQ.0 ) THEN !IF ( ABS(volume2-volume1).GT.0.00001 ) THEN
                       factor = dia(i)*0.01
                    ENDIF

                    IF  ( (volume1-volume(i)) .GT. (accept_sigma_it*volume(i)) ) THEN  !! So Dia to important or min_stomate
                       IF ( signe.EQ.1. .AND. .NOT.init_ok) THEN
                          factor = factor / factor_div_it
                          signe = -1.
                          IF ( factor .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                             dia_ok = .true.
                          ENDIF
                       ENDIF
                       num_it2 = 0
                       DO WHILE ( (dia(i)-factor).LT.min_stomate .AND. (num_it2.LT.5.) .AND. ((signe.EQ.-1.).OR.(.NOT.init_ok)) )    !! On ne peut pas utiliser ce facteur..., car trop important ou changement de sens...
                          factor = factor / factor_div_it
                          IF ( factor .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                             factor = 2*dia(i) 
                             dia_ok = .true.
                          ENDIF
                          num_it2 = num_it2 + 1.
                          IF ( num_it2.GE.5. ) write(*,*) "small iteration in lpj_crown.f90 need to be check (Arsene)"
                       ENDDO

                       IF ( .NOT.dia_ok ) THEN
                          IF ( signe.EQ.1. .AND. init_ok  ) THEN !! Arsene 02-09-2015 - Alors on reviens au début et on diminue the factor
                             dia(i) = dia(i) - factor !! on revient à l'état préscedent ;)
                             factor = factor / factor_div_it
                             IF ( factor .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                                dia_ok = .true.
                             ENDIF
                          ELSE
                             IF (init_ok .AND. num_it.GE.1) init_ok=.FALSE.
                             signe = -1
                             dia(i) = dia(i) - factor
                          ENDIF
                       ENDIF
                    ELSEIF ( (volume1-volume(i)) .LT. (accept_sigma_it*volume(i)) ) THEN !! dia too low  min_stomate
                       IF ( dia(i).GT.maxdia(j) ) THEN
                          dia(i) = maxdia(j)
                          dia_ok = .true.
                       ELSEIF ( signe.EQ.-1. .AND. num_it.NE.1 ) THEN
                          factor = factor / factor_div_it
                          IF ( factor .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                             dia_ok = .true.
                          ENDIF
                       ENDIF
                       IF ( .NOT.dia_ok ) THEN
                          IF ( signe.EQ.-1. .AND. init_ok ) THEN !! Arsene 02-09-2015 - Alors on reviens au début et on diminue the factor
                             dia(i) = dia(i) + factor !! on revient à l'état préscedent ;)
                             factor = factor / factor_div_it
                             IF ( factor .LT. (0.1 * accept_sigma_it * dia(i)) ) THEN
                                dia_ok = .true.
                             ENDIF
                          ELSE
                             IF (init_ok .AND. num_it.GE.1.) init_ok=.FALSE.
                             signe = 1.
                             dia(i) = dia(i) + factor
                          ENDIF
                       ENDIF
                    ELSE !! Good dia
                       dia_ok = .true.
                    ENDIF

                    num_it = num_it+1
                    IF ((num_it .GE. 100 ) .AND. (.NOT.dia_ok) ) THEN ! Si trop de boucle... limit à 100 ?1
                       dia_ok = .true.
                       write(*,*) "The iteration in lpj_crown.f90 need probably to be check (Arsene)"
                       !! Arsene 11-08-2015 - By default: Dia = last dia...
                    ENDIF

                 ENDDO  !! While loops!! Arsene 11-08-2015 Attention à tester l'itération, notemment pour ajuster the "accept_sigma"
!write(*,*) "lpj_crov it=", num_it, "dia=", dia(i), "dia_cut", dia_cut(:,j), "height before", height(:,j), "vol", volume(i)

                 dia(i)=MIN(dia(i),maxdia(j))

                ENDIF   !! if PFT present...
               ENDDO    !! Map loops

               WHERE (PFTpresent(:,j) .AND.woodmass_ind(:,j).GT.min_stomate)
                  !! 2.1.1.3.bis. Individual tree height from Aiba & Kohyama
                  WHERE ( dia(:).GT.dia_cut(:,j) ) !! 2.1.1.3.bis. Individual tree height from Aiba & Kohyama
                     height(:,j) = height_presc(j) * pt2 * 100**pt3 * dia(:)**pt3 &
                         & / ( height_presc(j) + pt2 * 100**pt3 * dia(:)**pt3 )

                     !! If dia_cut is below dia and no null ==> desactivate (= null)
                     WHERE ( dia_cut(:,j).GT.min_stomate )
                        dia_cut(:,j) = zero
                     ENDWHERE

                  ELSEWHERE ( dia(:).LT.dia_cut(:,j) .AND. dia_cut(:,j).GT.min_stomate )
                     !! On fix le diametre (identique à celui qu'il y avait lors de la coupe)
                     dia(:) = dia_cut(:,j)
                     !! On a déjà le nombre d'individus fixé... ATTENTION QUE CELA SOIT BIEN COMPATIBLE ==> PAS D'augmentation de ind tant que dia<dia_cut
                     !! On calcul la wodmass_ind "comme si on avait height correspondant à dia_cut
                     height(:,j) = height_presc(j) *  pt2* 100**pt3 * dia(:)**pt3 &
                            & / ( height_presc(j) + pt2 * 100**pt3 * dia(:)**pt3 )
                     woodmass_ind2(:) = pdensity * pi/4. * dia(:)**2 * height(:,j)
                     !! On peut donx comparer la woodmass réel par rapport à la woodmass attendu, et en déduire height !
                     height(:,j)=woodmass_ind(:,j)/woodmass_ind2(:) * height(:,j)
                  ELSEWHERE 
                     height(:,j) = zero      
                  ENDWHERE

                  !! 2.1.1.4 Crown area of individual tree
                  !          Calculate crown area, truncate crown area for trunks with large diameters
                  ! crown area cannot exceed a certain value, prescribed through maxdia
                  cn_ind(:,j) = pt1 * (ind(:,j)*pi/4*dia(:)**2)**ptcoeff / ind(:,j)
               ENDWHERE

             ENDIF                                            !! Arsene 11-08-2015 - Change for shrub allometry

          ELSE

             !! 2.1.2 Agricultural tree
             !        To be developped if needed
             STOP 'crown: cannot treat agricultural trees.'
          ENDIF
       ELSE
          
       !! 2.2 Grasses
          
          WHERE (PFTpresent(:,j))

             !! 2.2.1 Crown area of grass
             !        An "individual" is 1 m^2 of grass
             cn_ind(:,j) = un
          ENDWHERE
       ENDIF
       
       !! 2.3 Recalculate vegetation cover 
       
       !!!$SZ: since now all state variables are defined on veget_max it is very
       !!!$ dangerous to change this several times in stomate_lpj, as then GPP, turnover and allocated 
       !!!$ biomass are not defined on the same space! Hence, veget_max is now kept constant
       !!!$ and updated at the end of stomate_lpj in lpj_cover.f90
       !!!$ Eventually, this routine should only be called once at the beginning and the end of stomate_lpj
       !!!$ or prefereably cn_ind made a saved state variable!
       !!!$ IF (natural(j).AND.control%ok_dgvm) THEN
       !!!$   veget_max(:,j) = ind(:,j) * cn_ind(:,j)
       !!!$ ENDIF

    ENDDO ! loop over PFTs

  END SUBROUTINE crown

END MODULE lpj_crown
