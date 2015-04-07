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
    REAL(r_std), DIMENSION(nvm)                               :: mindia !! Arsene 29-10-2014
!_ ================================================================================================================================

    DO j = 2,nvm

      ! only when the DGVM is not activated or agricultural PFT.

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

                ratio=2.
                IF ( is_shrub(j) .AND. ind(i,j).GT.min_stomate) THEN

                !
                !! 1. We have to compare the height save and the height from number of ind ==> height2 
                !

                  height2 = pipe_tune2 * (veget_max(i,j) / (ind(i,j) * pipe_tune1))**(pipe_tune3/pipe_tune_exp_coeff)

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
                              .AND. height(i,j).GT.(height_presc(j)/5) ) THEN !! .AND. height(i,j).GT.minheight) ??
!write(*,*) "on redimensionne"
!read(*,*)
                  !! 3. In this case the shrub maintain the basis of its cylinder, but the lost in only on height
                  !!      That mean: ind=cst & cn_ind=cst & dia=cst. We have ind ==> we can recalculate the both other
                  !

                    cn_ind(i,j) = veget_max(i,j) / ind(i,j)
                    dia = (cn_ind(i,j)/pipe_tune1)**(1/pipe_tune_exp_coeff)

                  !
                  !! 4. With a homogeneous vertical density on shrub, the lost of biomass is a linear fonction with height
                  !!      We can calcul the height with simple comparison between theorical woodmass (with the dia) and the woodmass
                  !
                    woodmass2 = ind(i,j) * pipe_density*pi/4.*pipe_tune2 * (height2/(pipe_tune2))**((2.+pipe_tune3)/ pipe_tune3)

                  !!      After we deduce ce ratio between bot and we calculate the new height
                    ratio=woodmass(i)/woodmass2
                    height(i,j)=ratio*height2

                  ENDIF
                ENDIF

                IF ( ratio .GE. 1. ) THEN
                   !
                   !! 5. Calculate the provisionnal vegetation area
                   !

                   !! 5.1 Calculate stem diameter per individual tree (dia)        !! Arsene 03-09-2014
                   ! With min and max dia (to don't have no realistique value when DGVM is desactivate)
                   !dia(i) = ( height(i,j) / pipe_tune2 ) **(1/ pipe_tune3 ) !! Arsene 29-10-2014
                   maxdia(j)= (height_presc(j)/(pipe_tune2))**(1/ pipe_tune3)
                   mindia(j)= (height_presc(j)/(10*pipe_tune2))**(1/ pipe_tune3)       !! Arsene 02-04-2015 Pour le moment, height_max=height_presc & height_min = height_presc/10 (see stomate_lpj)

                   dia(i) = MAX(MIN( maxdia(j), ( woodmass(i) / (veget_max(i,j) * pipe_density*pi/4. &
                         & *pipe_tune2/pipe_tune1) )**(1/(2.+pipe_tune3-pipe_tune_exp_coeff)) ),mindia(j))

                   !! Because 1. ind(i,j) = woodmass(i) / ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) ) 
                   !!         2. cn_ind(i,j) * ind(i,j) =  veget_max(i,j) because without DGVM
                   !!         3. cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff avec MIN( maxdia(j), dia(i)) = dia(i)
                   !!         4. simplification
                   !!         5. dia(i) = MAX( (MIN( maxdia(j), dia(i))), mindia )

                   !! 5.2 So, 1. give:
                   ind(i,j) = woodmass(i) / &
                           ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) )

                   !! 5.3 Individual biomass corresponding to this critical density of individuals

                   woodmass_ind(i) = woodmass(i) / ind(i,j)
                   !So : woodmass_ind(i) =  ( pipe_density*pi/4.*pipe_tune2 * dia(i)**(2.+pipe_tune3) )


                   !! 5.4 Calculate provisional tree crown area for per individual tree

                   ! equation: CrownArea=pipe_tune1 * Diameter^{1.6}
                   cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

                   !
                   !! 6. If total tree crown area for this tree PFT exceeds its veget_max, tree density is recalculated.
                   !
  
                   IF ( cn_ind(i,j) * ind(i,j) .GT. 1.002* veget_max(i,j) ) THEN
                       ind(i,j) = veget_max(i,j) / cn_ind(i,j)
                   ELSEIF ( cn_ind(i,j) * ind(i,j) .LT. 1.002* veget_max(i,j) ) THEN
                       cn_ind(i,j) = veget_max(i,j) / ind(i,j)
                   ENDIF

                   !
                   !! 7. In this case, we recalculate the height
                   !

                   woodmass_ind(i) = woodmass(i) / ind(i,j)

                   dia(i) = ( woodmass_ind(i) / ( pipe_density * pi/4. * pipe_tune2 ) ) ** &
                              ( un / ( 2. + pipe_tune3 ) )

                   height(i,j)=pipe_tune2 * dia(i)**(pipe_tune3)

                   IF ( height(i,j) .GT. height_presc(j) ) THEN
                       height(i,j)=height_presc(j)
                   ENDIF

                 ENDIF  !! Arsene 30-10-2014

!! Arsene 30-10-2014 -END- New Version... Because the version before was tooooooooo full of bugs
!!!!!!!!!!!!!!!!!!!!!!!!!!

              ELSE !woodmas=0  => impose some value

                dia(:) = maxdia(j)

                cn_ind(i,j) = pipe_tune1 * MIN( maxdia(j), dia(i) ) ** pipe_tune_exp_coeff

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
