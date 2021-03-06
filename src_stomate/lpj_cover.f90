! =================================================================================================================================
! MODULE       : lpj_cover
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
!                This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : None
!!
!! RECENT CHANGE(S) : Including permafrost carbon
!!
!! REFERENCE(S) : 
!!        Sitch, S., B. Smith, et al. (2003), Evaluation of ecosystem dynamics,
!!        plant geography and terrestrial carbon cycling in the LPJ dynamic 
!!        global vegetation model, Global Change Biology, 9, 161-185.\n
!!        Smith, B., I. C. Prentice, et al. (2001), Representation of vegetation
!!        dynamics in the modelling of terrestrial ecosystems: comparing two
!!        contrasting approaches within European climate space,
!!        Global Ecology and Biogeography, 10, 621-637.\n
!!
!! SVN :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE/ORCHIDEE/src_stomate/lpj_cover.f90 $
!! $Date: 2013-10-14 15:38:24 +0200 (Mon, 14 Oct 2013) $
!! $Revision: 1536 $
!! \n
!_ ================================================================================================================================

MODULE lpj_cover

  ! modules used:

  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes_soil_var

  IMPLICIT NONE

  ! private & public routines

  PRIVATE
  PUBLIC cover

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE     : lpj_cover
!!
!>\BRIEF          Recalculate vegetation cover and LAI
!!
!!\n DESCRIPTION : Veget_max is first renewed here according to newly calculated foliage biomass in this calculation step 
!! Then, litter, soil carbon, and biomass are also recalcuted with taking into account the changes in Veget_max (i.e. delta_veg)
!! Grid-scale fpc (foliage projected coverage) is calculated to obtain the shadede ground area by leaf's light capture
!! Finally, grid-scale fpc is adjusted not to exceed 1.0
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::lai (leaf area index, @tex $(m^2 m^{-2})$ @endtex), 
!! :: veget (fractional vegetation cover, unitless)
!!
!! REFERENCE(S)   : None
!! 
!! FLOWCHART :
!! \latexonly 
!!     \includegraphics[scale=0.5]{lpj_cover_flowchart.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE cover (npts, cn_ind, ind, biomass, &
       veget_max, veget_max_old, lai, litter, litter_avail, litter_not_avail, &
       carbon, turnover_daily, bm_to_litter,&
       deepC_a, deepC_s, deepC_p)

!! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                                  :: npts             !! Domain size (unitless)  
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: cn_ind           !! Crown area 
                                                                                    !! @tex $(m^2)$ @endtex per individual
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: ind              !! Number of individuals 
                                                                                    !! @tex $(m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)                :: veget_max_old    !! "Maximal" coverage fraction of a PFT (LAI-> 
                                                                                    !! infinity) on ground at beginning of time 

    !! 0.2 Output variables

    !! 0.3 Modified variables

    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)             :: lai                 !! Leaf area index OF AN INDIVIDUAL PLANT 
                                                                                       !! @tex $(m^2 m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout) :: litter    !! Metabolic and structural litter, above and 
                                                                                       !! below ground @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nlitt,nvm), INTENT(inout):: litter_avail
    REAL(r_std), DIMENSION(npts,nlitt,nvm) , INTENT(inout):: litter_not_avail
    REAL(r_std), DIMENSION(npts,ncarb,nvm), INTENT(inout)             :: carbon        !! Carbon pool: active, slow, or passive 
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: biomass        !! Biomass @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)                  :: veget_max      !! "Maximal" coverage fraction of a PFT (LAI->
                                                                                       !! infinity) on ground (unitless)
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: turnover_daily !! Turnover rates (gC m^{-2} day^{-1})
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter   !! Conversion of biomass to litter 
                                                                                       !! @tex $(gC m^{-2} day^{-1})$ @endtex
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_a           !! Permafrost soil carbon (g/m**3) active
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_s           !! Permafrost soil carbon (g/m**3) slow
    REAL(r_std), DIMENSION(npts,ndeep,nvm), INTENT(inout)         :: deepC_p           !! Permafrost soil carbon (g/m**3) passive

    !! 0.4 Local variables

    INTEGER(i_std)                                              :: i,j,k,m          !! Index (unitless)
    REAL(r_std), DIMENSION(npts,nlitt,nlevs,nelements)          :: dilu_lit         !! Litter dilution @tex $(gC m^{-2})$ @endtex
    REAL(r_std), DIMENSION(npts,ncarb)                          :: dilu_soil_carbon !! Soil Carbondilution 
                                                                                    !! @tex $(gC m^{-2})$ @endtex
    REAL(r_std),DIMENSION(npts,ndeep,ncarb)                     :: dilu_soil_carbon_vertres !!vertically-resolved Soil Carbondilution (gC/m²)

    REAL(r_std), DIMENSION(nvm)                                 :: delta_veg        !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(nvm)                                 :: reduct           !! Conversion factors (unitless)
    REAL(r_std)                                                 :: delta_veg_sum    !! Conversion factors (unitless)
    REAL(r_std)                                                 :: diff             !! Conversion factors (unitless)
    REAL(r_std)                                                 :: sr               !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: frac_nat         !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: sum_vegettree    !! Conversion factors (unitless)
    REAL(r_std), DIMENSION(npts)                                :: sum_vegetgrass   !! Conversion factors (unitless) 
    REAL(r_std), DIMENSION(npts)                                :: sum_veget_natveg !! Conversion factors (unitless)

!_ ================================================================================================================================

 !! 1. If the vegetation is dynamic, calculate new maximum vegetation cover for natural plants
  
    IF ( control%ok_dgvm ) THEN

       !! 1.1  Calculate initial values of vegetation cover
       frac_nat(:) = un
       sum_veget_natveg(:) = zero
       veget_max(:,ibare_sechiba) = un

       DO j = 2,nvm ! loop over PFTs

          IF ( natural(j) ) THEN
	     
             ! Summation of individual tree crown area to get total foliar projected coverage
             veget_max(:,j) = ind(:,j) * cn_ind(:,j)
             sum_veget_natveg(:) = sum_veget_natveg(:) + veget_max(:,j)

          ELSE
             
             !fraction occupied by agriculture needs to be substracted for the DGVM
             !this is used below to constrain veget for natural vegetation, see below
             frac_nat(:) = frac_nat(:) - veget_max(:,j)

          ENDIF

       ENDDO ! loop over PFTs

       DO i = 1, npts ! loop over grid points
	  
          ! Recalculation of vegetation projected coverage when ::frac_nat was below ::sum_veget_natveg
          ! It means that non-natural vegetation will recover ::veget_max as natural vegetation
          IF (sum_veget_natveg(i) .GT. frac_nat(i) .AND. frac_nat(i) .GT. min_stomate) THEN

             DO j = 2,nvm ! loop over PFTs
                IF( natural(j) ) THEN
                   veget_max(i,j) =  veget_max(i,j) * frac_nat(i) / sum_veget_natveg(i)
                ENDIF
             ENDDO ! loop over PFTs

          ENDIF
       ENDDO ! loop over grid points
	
       ! Renew veget_max of bare soil as 0 to difference of veget_max (ibare_sechiba) 
       ! to current veget_max
       DO j = 2,nvm ! loop over PFTs
          veget_max(:,ibare_sechiba) = veget_max(:,ibare_sechiba) - veget_max(:,j)
       ENDDO ! loop over PFTs
       veget_max(:,ibare_sechiba) = MAX( veget_max(:,ibare_sechiba), zero )

       !! 1.2 Calculate carbon fluxes between PFTs to maintain mass balance
       !      Recalculate the litter and soil carbon with taking into accout the change in 
       !      veget_max (delta_veg)
       DO i = 1, npts ! loop over grid points
          
          ! calculate the change in veget_max between previous time step and current time step
          delta_veg(:) = veget_max(i,:)-veget_max_old(i,:)
          delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.zero)

          dilu_lit(i,:,:,:) = zero
          dilu_soil_carbon(i,:) = zero
          dilu_soil_carbon_vertres(i,:,:) = zero
          DO j=1, nvm ! loop over PFTs
             IF ( delta_veg(j) < -min_stomate ) THEN
                dilu_lit(i,:,:,:) =  dilu_lit(i,:,:,:) + delta_veg(j) * litter(i,:,j,:,:) / delta_veg_sum
                IF ( ok_pc ) THEN
                    dilu_soil_carbon_vertres(i,:,iactive)= dilu_soil_carbon_vertres(i,:,iactive) - &
                         delta_veg(j) * deepC_a(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,islow)= dilu_soil_carbon_vertres(i,:,islow) - &
                         delta_veg(j) * deepC_s(i,:,j) / delta_veg_sum
                    dilu_soil_carbon_vertres(i,:,ipassive)= dilu_soil_carbon_vertres(i,:,ipassive) - &
                         delta_veg(j) * deepC_p(i,:,j) / delta_veg_sum
                ELSE
                    dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
                ENDIF
             ENDIF
          ENDDO ! loop over PFTs

          DO j=1, nvm ! loop over PFTs
             IF ( delta_veg(j) > min_stomate) THEN

                ! Dilution of reservoirs
                ! Recalculate the litter and soil carbon with taking into accout the change in 
                ! veget_max (delta_veg)
                ! Litter
                litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_max_old(i,j) + dilu_lit(i,:,:,:) * delta_veg(j)) / veget_max(i,j)
                !JCADD available and not available litter for grazing
                ! only not available litter change, available litter will not
                ! change, because tree litter can not be eaten
               IF (is_grassland_manag(j) .AND. is_grassland_grazed(j)) THEN
                 litter_avail(i,:,j) = litter_avail(i,:,j) * veget_max_old(i,j) / veget_max(i,j)
                 litter_not_avail(i,:,j) = litter(i,:,j,iabove,icarbon) - litter_avail(i,:,j)
               ENDIF
                !ENDJCADD    
                IF ( ok_pc ) THEN
                   deepC_a(i,:,j)=(deepC_a(i,:,j) * veget_max_old(i,j) + &
                        dilu_soil_carbon_vertres(i,:,iactive) * delta_veg(j)) / veget_max(i,j)
                   deepC_s(i,:,j)=(deepC_s(i,:,j) * veget_max_old(i,j) + &
                        dilu_soil_carbon_vertres(i,:,islow) * delta_veg(j)) / veget_max(i,j)
                   deepC_p(i,:,j)=(deepC_p(i,:,j) * veget_max_old(i,j) + &
                        dilu_soil_carbon_vertres(i,:,ipassive) * delta_veg(j)) / veget_max(i,j)
                ELSE
                   ! Soil carbon
                   carbon(i,:,j)=(carbon(i,:,j) * veget_max_old(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_max(i,j)
                ENDIF

             ENDIF

             IF((j.GE.2).AND.(veget_max_old(i,j).GT.min_stomate).AND.(veget_max(i,j).GT.min_stomate)) THEN

                ! Correct biomass densities (i.e. also litter fall) to conserve mass 
                ! since it's defined on veget_max
                biomass(i,j,:,:) = biomass(i,j,:,:) * veget_max_old(i,j) / veget_max(i,j)
                turnover_daily(i,j,:,:) = turnover_daily(i,j,:,:) * veget_max_old(i,j) / veget_max(i,j)
                bm_to_litter(i,j,:,:) = bm_to_litter(i,j,:,:) * veget_max_old(i,j) / veget_max(i,j)

             ENDIF
          ENDDO ! loop over PFTs
       ENDDO ! loop over grid points
    ENDIF

  END SUBROUTINE cover

END MODULE lpj_cover
