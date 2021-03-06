! =================================================================================================================================
! MODULE 	: stomate_data
!
! CONTACT      : orchidee-help _at_ ipsl.jussieu.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "stomate_data" module defines the values about the PFT parameters. It will print
!! the values of the parameters for STOMATE in the standard outputs. 
!!
!!\n DESCRIPTION: None 
!!
!! RECENT CHANGE(S): Sonke Zaehle: Reich et al, 1992 find no statistically significant differences 
!!                  between broadleaved and coniferous forests, specifically, the assumption that grasses grow 
!!                  needles is not justified. Replacing the function with the one based on Reich et al. 1997. 
!!                  Given that sla=100cm2/gDW at 9 months, sla is:
!!                  sla=exp(5.615-0.46*ln(leaflon in months))
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/trunk/ORCHIDEE/src_stomate/stomate_data.f90 $
!! $Date: 2014-09-04 14:46:14 +0200 (Thu, 04 Sep 2014) $
!! $Revision: 2282 $
!! \n
!_ ================================================================================================================================

MODULE stomate_data

  ! modules used:

  USE constantes
  USE pft_parameters
  USE defprec
  

  IMPLICIT NONE

  INTEGER(i_std),SAVE :: bavard=1                 !! Level of online diagnostics in STOMATE (0-4, unitless)
!$OMP THREADPRIVATE(bavard)

  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:) :: hori_index     !! Move to Horizontal indices
!$OMP THREADPRIVATE(hori_index)

  INTEGER(i_std),ALLOCATABLE,SAVE,DIMENSION(:) :: horipft_index  !! Horizontal + PFT indices
!$OMP THREADPRIVATE(horipft_index)

  ! Land cover change

  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip10_index   !! Horizontal + P10 indices
!$OMP THREADPRIVATE(horip10_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip100_index  !! Horizontal + P100 indice
!$OMP THREADPRIVATE(horip100_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip11_index   !! Horizontal + P11 indices
!$OMP THREADPRIVATE(horip11_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horip101_index  !! Horizontal + P101 indices
!$OMP THREADPRIVATE(horip101_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horideep_index
!$OMP THREADPRIVATE(horideep_index)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION(:) :: horisnow_index
!$OMP THREADPRIVATE(horisnow_index)
  INTEGER(i_std),SAVE :: itime                 !! time step
!$OMP THREADPRIVATE(itime)
  INTEGER(i_std),SAVE :: hist_id_stomate       !! STOMATE history file ID
!$OMP THREADPRIVATE(hist_id_stomate)
  INTEGER(i_std),SAVE :: hist_id_stomate_IPCC  !! STOMATE history file ID for IPCC output
!$OMP THREADPRIVATE(hist_id_stomate_IPCC)
  INTEGER(i_std),SAVE :: rest_id_stomate       !! STOMATE restart file ID
!$OMP THREADPRIVATE(rest_id_stomate)

  REAL(r_std),PARAMETER :: adapted_crit = 1. - ( 1. / euler ) !! critical value for being adapted (1-1/e) (unitless)
  REAL(r_std),PARAMETER :: regenerate_crit = 1. / euler       !! critical value for being regenerative (1/e) (unitless)

  ! private & public routines

  PUBLIC data

CONTAINS

!! ================================================================================================================================
!! SUBROUTINE 	: data
!!
!>\BRIEF         This routine defines the values of the PFT parameters. It will print the values of the parameters for STOMATE
!!               in the standard outputs of ORCHIDEE. 
!!
!! DESCRIPTION : This routine defines PFT parameters. It initializes the pheno_crit structure by tabulated parameters.\n
!!               Some initializations are done for parameters. The SLA is calculated according *to* Reich et al (1992).\n
!!               Another formulation by Reich et al(1997) could be used for the computation of the SLA.
!!               The geographical coordinates might be used for defining some additional parameters
!!               (e.g. frequency of anthropogenic fires, irrigation of agricultural surfaces, etc.). \n
!!               For the moment, this possibility is not used. \n
!!               The specifc leaf area (SLA) is calculated according Reich et al, 1992 by :
!!               \latexonly
!!               \input{stomate_data_SLA.tex}
!!               \endlatexonly
!!               The sapling (young) biomass for trees and for each compartment of biomass is calculated by :
!!               \latexonly
!!               \input{stomate_data_sapl_tree.tex}
!!               \endlatexonly
!!               The sapling biomass for grasses and for each compartment of biomass is calculated by :
!!               \latexonly
!!               \input{stomate_data_sapl_grass.tex}
!!               \endlatexonly
!!               The critical stem diameter is given by the following formula :
!!               \latexonly
!!               \input{stomate_data_stem_diameter.tex}
!!               \endlatexonly
!!
!! RECENT CHANGE(S): Sonke Zaehle: Reich et al, 1992 find no statistically significant differences 
!!                  between broadleaved and coniferous forests, specifically, the assumption that grasses grow 
!!                  needles is not justified. Replacing the function with the one based on Reich et al. 1997. 
!!                  Given that sla=100cm2/gDW at 9 months, sla is:
!!                  sla=exp(5.615-0.46*ln(leaflon in months)) 
!!                   \latexonly
!!                   \input{stomate_data_SLA_Reich_97.tex}
!!                   \endlatexonly
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) :
!! - Reich PB, Walters MB, Ellsworth DS, (1992), Leaf life-span in relation to leaf, plant and 
!! stand characteristics among diverse ecosystems. Ecological Monographs, Vol 62, pp 365-392.
!! - Reich PB, Walters MB, Ellsworth DS (1997) From tropics to tundra: global convergence in plant 
!!  functioning. Proc Natl Acad Sci USA, 94:13730 13734
!!
!! FLOWCHART    :
!! \n
!_ ================================================================================================================================

  SUBROUTINE data (npts, lalo)


    !! 0. Variables and parameter declaration


    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                   :: npts    !! [DISPENSABLE] Domain size (unitless)
    REAL(r_std),DIMENSION (npts,2), INTENT (in)  :: lalo    !! [DISPENSABLE] Geographical coordinates (latitude,longitude)

    !! 0.4 Local variables

    INTEGER(i_std)                               :: j       !! Index (unitless)
    REAL(r_std)                                  :: alpha   !! alpha's : (unitless)
    REAL(r_std)                                  :: dia     !! stem diameter (m)
    REAL(r_std)                                  :: csa_sap !! Crown specific area sapling @tex $(m^2.ind^{-1})$ @endtex

!! Arsene 16-10-2015 Add for shrub allometry
    INTEGER(i_std)                               :: ii, shrub_height_i
    REAL(r_std)                                  :: shrub_min_h, shrub_max_h, shrub_h_x, shrub_h_cst
    REAL(r_std)                                  :: vol_min, vol_max, shrub_vol_cst
    REAL(r_std)                                  :: vol_last, vol_next, height_last, height_next
!! Arsene 29-12-2015 - ADD - LUT for new litter moist dependence
    REAL(r_std)                                  :: moist_func, moist_max

!    REAL(r_std)                                  :: ind_frac!! fraction of real individu by ind (special for shrub) !! 22-05-2015 Arsene

!_ ================================================================================================================================

    IF ( bavard .GE. 1 ) WRITE(numout,*) 'data: PFT characteristics'

    !- pheno_gdd_crit
    pheno_gdd_crit(:,1) = pheno_gdd_crit_c(:)
    pheno_gdd_crit(:,2) = pheno_gdd_crit_b(:)         
    pheno_gdd_crit(:,3) = pheno_gdd_crit_a(:) 
    !
    !- senescence_temp
    senescence_temp(:,1) = senescence_temp_c(:)
    senescence_temp(:,2) = senescence_temp_b(:)
    senescence_temp(:,3) = senescence_temp_a(:)
    !
    !- maint_resp_slope
    maint_resp_slope(:,1) = maint_resp_slope_c(:)              
    maint_resp_slope(:,2) = maint_resp_slope_b(:)
    maint_resp_slope(:,3) = maint_resp_slope_a(:)
    !
    !-coeff_maint_zero
    coeff_maint_zero(:,ileaf) = cm_zero_leaf(:)
    coeff_maint_zero(:,isapabove) = cm_zero_sapabove(:)
    coeff_maint_zero(:,isapbelow) = cm_zero_sapbelow(:)
    coeff_maint_zero(:,iheartabove) = cm_zero_heartabove(:)
    coeff_maint_zero(:,iheartbelow) = cm_zero_heartbelow(:)
    coeff_maint_zero(:,iroot) = cm_zero_root(:)
    coeff_maint_zero(:,ifruit) = cm_zero_fruit(:)
    coeff_maint_zero(:,icarbres) = cm_zero_carbres(:)


    IF ( bavard .GE. 1 ) WRITE(numout,*) 'data: PFT characteristics'

    DO j = 2,nvm ! Loop over # PFTS 

       IF ( bavard .GE. 1 ) WRITE(numout,'(a,i3,a,a)') '    > PFT#',j,': ', PFT_name(j)

       !
       ! 1 tree? (true/false)
       !
       IF ( bavard .GE. 1 ) WRITE(numout,*) '       tree: (::is_tree) ', is_tree(j)
       !
       ! 1.2 shrub? (true/false)                                                             !! Arsene 31-07-2014 modifications
       !                                                                                     !! Arsene 31-07-2014 modifications
       IF ( bavard .GE. 1 ) WRITE(numout,*) '       shrub: (::is_shrub) ', is_shrub(j)       !! Arsene 31-07-2014 modifications

       !
       ! 2 flamability (0-1, unitless)
       !

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       litter flamability (::flam) :', flam(j)

       !
       ! 3 fire resistance (unitless)
       !

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       fire resistance (::resist) :', resist(j)

       !
       ! 4 specific leaf area per mass carbon = 2 * sla / dry mass (m^2.gC^{-1})
       !

       ! SZ: Reich et al, 1992 find no statistically significant differences between broadleaved and coniferous
       ! forests, specifically, the assumption that grasses grow needles is not justified. Replacing the function
       ! with the one based on Reich et al. 1997. Given that sla=100cm2/gDW at 9 months, sla is:
       ! sla=exp(5.615-0.46*ln(leaflon in months))

       ! Oct 2010 : sla values are prescribed by values given by N.Viovy 

       ! includes conversion from 
       !!       sla(j) = 2. * 1e-4 * EXP(5.615 - 0.46 * log(12./leaflife_tab(j)))
       !!\latexonly
       !!\input{stomate_data_SLA.tex}
       !!\endlatexonly
!       IF ( leaf_tab(j) .EQ. 2 ) THEN
!
!          ! needle leaved tree
!          sla(j) = 2. * ( 10. ** ( 2.29 - 0.4 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!
!       ELSE
!
!          ! broad leaved tree or grass (Reich et al 1992)
!          sla(j) = 2. * ( 10. ** ( 2.41 - 0.38 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!
!       ENDIF

!!!$      IF ( leaf_tab(j) .EQ. 1 ) THEN
!!!$
!!!$        ! broad leaved tree
!!!$
!!!$        sla(j) = 2. * ( 10. ** ( 2.41 - 0.38 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!!!$
!!!$      ELSE
!!!$
!!!$        ! needle leaved or grass (Reich et al 1992)
!!!$
!!!$        sla(j) = 2. * ( 10. ** ( 2.29 - 0.4 * LOG10(12./leaflife_tab(j)) ) ) *1e-4
!!!$
!!!$      ENDIF
!!!$
!!!$      IF ( ( leaf_tab(j) .EQ. 2 ) .AND. ( pheno_type_tab(j) .EQ. 2 ) ) THEN
!!!$
!!!$        ! summergreen needle leaf
!!!$
!!!$        sla(j) = 1.25 * sla(j)
!!!$
!!!$      ENDIF

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       specific leaf area (m**2/gC) (::sla):', sla(j), 12./leaflife_tab(j)

       !
       ! 7 critical stem diameter: beyond this diameter, the crown area no longer
       !     increases (m)
       !

       IF ( is_tree(j) ) THEN       !! Arsene 31-07-2014 modifications ATTENTION. A vérif

          !!\latexonly
          !!\input{stomate_data_stem_diameter.tex}
          !!\endlatexonly

!! Arsene 03-09-2014 - It's possible to define maxdia with - TreeHeight=pipe_tune2 * Diameter^{pipe_tune3} ==> Definition in stomate_establish.f90 or lpj_crown.f90
!           maxdia(j) = ( height_presc(j) / pipe_tune2 ) **(1/ pipe_tune3) 

!!          maxdia(j) = ( ( pipe_tune4 / ((pipe_tune2*pipe_tune3)/(maxdia_coeff(1)**pipe_tune3)) ) &     !! Arsene 11-08-2015 BEFORE CHANGES - NOTE: Constant, not PFT dependent...
!!               ** ( un / ( pipe_tune3 - un ) ) ) * maxdia_coeff(2)                                     !! Arsene 11-08-2015 BEFORE CHANGES - Note: maxdia_coeff & pipe_tune4 not use now
!! Arsene 30-03-2015 ==> I think better to change, but not yet !!!!!!!!!!!!!!!!!!! ==> !! Arsene 11-08-2015 - YET !

          maxdia(j) = ( height_presc(j) / pipe_tune2 ) **(1/ pipe_tune3)         !! Arsene 11-08-2015 - New ! - Note:  height_presc(j) = max height
          mindia(j)= (height_presc(j)/(fact_min_height*pipe_tune2))**(1/ pipe_tune3)

          cn_sapl(j) = cn_sapl_init !crown of individual tree, first year

       ELSEIF ( is_shrub(j) .AND. shrubs_like_trees ) THEN  !! Arsene 03-08-2015 - Change pipe_tune for shrubs

!          maxdia(j) = ( ( pipe_tune4 / ((pipe_tune2_for_shrub*pipe_tune3_for_shrub)/ &    !! Arsene 04-08-2015 - No pipe_tune_4_shrub
!              & (maxdia_coeff(1)**pipe_tune3_for_shrub)) ) ** ( un / ( pipe_tune3_for_shrub - un ) ) ) * maxdia_coeff(2)

          maxdia(j) = ( height_presc(j) / pipe_tune2_for_shrub ) **(1./ pipe_tune3_for_shrub)
          mindia(j)= (height_presc(j)/(fact_min_height*pipe_tune2_for_shrub))**(1./ pipe_tune3_for_shrub)

          cn_sapl(j) = cn_sapl_init !crown of individual tree, first year (like trees)

       ELSEIF ( is_shrub(j) .AND. .NOT.shrubs_like_trees ) THEN  !! Arsene 03-08-2015 - Change pipe_tune for shrubs

          maxdia(j) = ( height_presc(j) * shrub_lim_maxdia / (pipe_tune_shrub2*(1.-shrub_lim_maxdia)) )**(1./pipe_tune_shrub3) / 100. !! Arsene 11-08-2015 - 0.93 correspond au observatio [90 - 96]
          mindia(j) = ( height_presc(j) / (pipe_tune_shrub2 * (fact_min_height-1.)) )**(1./pipe_tune_shrub3) / 100.

          cn_sapl(j) = cn_sapl_init !crown of individual tree, first year (like trees)

       ELSE
          maxdia(j) = undef
          cn_sapl(j)=1
       ENDIF !( is_tree(j) )

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       critical stem diameter (m): (::maxdia(j))', maxdia(j)
       IF ( bavard .GE. 1 ) WRITE(numout,*) '       critical stem diameter (m): (::mindia(j))', mindia(j)


!!! Arsene 16-10-2015 - ADD - START Look Up Table (LUT)
       !
       ! 4.bis - For shrubs (.NOT.shrubsliketrees): Calcul of array
       ! This part is use only in lpj_crown.f90 and lpj_establish.f90 (not in stomate_prescribe) ==> only if DGVM or not lpj_cst
       !
       IF ( is_shrub(j) .AND. .NOT.shrubs_like_trees .AND. (control%ok_dgvm .OR. .NOT.lpj_gap_const_mort) ) THEN

           !! 4.bis.1 Compute the "virtual array" volume=fn(height)
           shrub_min_h = 0.01                                !! in m
           shrub_max_h = height_presc(j) * shrub_lim_maxdia  !! real max height (height_presc(j) = theorical max)
           shrub_h_x = 1.                                    !! fraction of "line" for "array" volume=fn(height), compute but not save...
           shrub_h_cst = ( shrub_max_h / shrub_min_h )**( 1./ (shrub_h_x * shrub_allom_lig -1.))
           !! ==> We have for each "ligne" i: Height = shrub_h_cst**(i-1.) * shrub_min_h
           !!                          puis   Vol(i) = (pi/4) * Height(i) * dia(i)**2

           !! 4.bis.2 Invert array dia=fn(vol): ligne discretisation
           vol_min = (pi/4) * shrub_min_h * ((1./(pipe_tune_shrub2 * (1./shrub_min_h-1./height_presc(j)))) &
                         **(1./pipe_tune_shrub3) /100. )**2
           vol_max = (pi/4) * shrub_max_h * ((1./(pipe_tune_shrub2 * (1./shrub_max_h-1./height_presc(j)))) &
                         **(1./pipe_tune_shrub3) /100. )**2
           shrub_vol_cst = ( vol_max / vol_min )**( 1./ (shrub_allom_lig -1.))
           !! ==> We have for each "ligne" ii: Vol = shrub_vol_cst**(ii-1) * vol_min

           !! 4.bis.3 Fist and last value of array
           shrub_allom_array(j,1,1) = vol_min
           shrub_allom_array(j,1,2) = shrub_min_h
           shrub_allom_array(j,shrub_allom_lig,1) = vol_max
           shrub_allom_array(j,shrub_allom_lig,2) = shrub_max_h
           shrub_allom_array(j,shrub_allom_lig+1,1) = shrub_vol_cst
           shrub_allom_array(j,shrub_allom_lig+1,2) = shrub_h_cst

           !! 4.bis.4 Initialise last and next value, to compute array
           height_last = 0.                ! not use ?
           height_next = shrub_min_h

           vol_last = 0.                   ! not use ?
           vol_next = vol_min

           !! 4.bis.5 Start to compute heigt value for each 
           shrub_height_i = 1
           DO ii = 2,shrub_allom_lig-1

                !! 4.bis.6 Compute next value of volume (for dia=fn(vol))
                shrub_allom_array(j,ii,1) = shrub_vol_cst**(ii-1) * vol_min

                !! 4.bis.7 Check if shrub_allom_array(j,ii,1) is between vol_last and vol_next (from vol=fn(dia)). If not, find goof ones.
                DO WHILE ((shrub_allom_array(j,ii,1) .GT. vol_next) .AND. ( shrub_height_i.LE.(shrub_h_x * shrub_allom_lig) ))
                     shrub_height_i = shrub_height_i + 1
                     vol_last = vol_next
                     height_last = height_next
                     height_next = shrub_h_cst**(shrub_height_i-1.) * shrub_min_h
                     vol_next = (pi/4) * height_next * ((1./(pipe_tune_shrub2 * (1./height_next-1./height_presc(j)))) &
                         **(1./pipe_tune_shrub3) /100. )**2
                ENDDO

                !! 4.bis.8 Compute and save good dia for dia=fn(vol) with linear interpolation
                !! It is possible to go outside the array... ?
                !! First we se "where" is the volume:       vol_fract = ( shrub_allom_array(j,ii,1) - vol_last ) / (vol_next - vol_last)
                !! Add the linear interpolation for height: shrub_allom_array(j,ii) = height_last + vol_fract * (height_next - height_last)
                shrub_allom_array(j,ii,2) = height_last + (( shrub_allom_array(j,ii,1) - vol_last ) / (vol_next - vol_last)) &
                        & * (height_next - height_last)
           ENDDO

       ENDIF

!!! Arsene 16-10-2015 - ADD - END       

!!! Arsene 29-12-2015 - ADD - START LUT for new litter and soil-carbon moist decomposition dependence

       IF ( new_moist_func ) THEN

           !! We calculate each value, to 0. until 1. with moist_inerval, from Moyano et al., 2012.
           moist_func = 0.           

           litter_moist_array(1) = moist_coeff_new(4)
           moist_max = litter_moist_array(1)
           DO ii = 2, int(1./moist_interval+1.)
               moist_func = moist_func + moist_interval
               litter_moist_array(ii) = (moist_coeff_new(1)*moist_func**3 + moist_coeff_new(2)*moist_func**2 + &
                                  &     moist_coeff_new(3)*moist_func + moist_coeff_new(4)) * litter_moist_array(ii-1)
               moist_max = MAX(moist_max,litter_moist_array(ii))
           ENDDO

           litter_moist_array(:) = litter_moist_array(:) / moist_max 

       ENDIF

!!! Arsene 29-12-2015 - ADD - END LUT


       !
       ! 5 sapling characteristics
       !

       IF ( is_tree(j) .OR. is_shrub(j) ) THEN   !! Arsene 31-07-2014 modifications  A REVOIR l'adaptation pour les shrubs

          !> 5.1 trees

          !!\latexonly
          !!\input{stomate_data_sapl_tree.tex}
          !!\endlatexonly
   
!! Arsene 11-08-2015 - Change all, beginning with dia = mindia, and reverse some équations...

!!*          IF ( is_tree(j) ) THEN          !! Arsene 31-07-2014 modifications
!!*              alpha = alpha_tree
!!*
!!*          bm_sapl(j,ileaf,icarbon) = &
!!*               &     (((bm_sapl_leaf(1)*pipe_tune1*(mass_ratio_heart_sap *bm_sapl_leaf(2)*sla(j)/(pi*pipe_k1)) &
!!*               &     **bm_sapl_leaf(3))/sla(j))**bm_sapl_leaf(4))
!!*          ELSE                            !! Arsene 31-07-2014 modifications
!!*              alpha = alpha_shrub         !! Arsene 31-07-2014 modifications
!!*
!!*          bm_sapl(j,ileaf,icarbon) = &     !! Arsene 03-08-2015 - Change pipe_tune for shrubs
!!*               &     (((bm_sapl_leaf(1)*pipe_tune1_for_shrub*(mass_ratio_heart_sap *bm_sapl_leaf(2)*sla(j)/(pi*pipe_k1_shrub)) &   !! Arsene 03-08-2015 - Change pipe_tune for shrubs
!!*               &     **bm_sapl_leaf(3))/sla(j))**bm_sapl_leaf(4))
!!*
!!*          ENDIF                           !! Arsene 31-07-2014 modifications
!!*
!!*          IF ( pheno_type(j) .NE. 1 ) THEN
!!*             ! not evergreen
!!*             bm_sapl(j,icarbres,icarbon) = bm_sapl_carbres * bm_sapl(j,ileaf,icarbon)
!!*          ELSE
!!*             bm_sapl(j,icarbres,icarbon) = zero
!!*          ENDIF ! (pheno_type_tab(j) .NE. 1 )
!!*
!!*          csa_sap = bm_sapl(j,ileaf,icarbon) / ( pipe_k1 / sla(j) )
!!*
!!*          dia = (mass_ratio_heart_sap * csa_sap * dia_coeff(1) / pi ) ** dia_coeff(2)


          !! 5.1.1 - Compute sapl diameter  = min dia:
          dia= mindia(j)*5.     !! Like that bm_sapl_leaf and dia_coeff can be remove !

          !! 5.2.2 - Compute cas_sap: Reverse the original equation:
          !!    dia = (mass_ratio_heart_sap * csa_sap * dia_coeff(1) / pi ) ** dia_coeff(2)
          csa_sap = pi * dia**(1/dia_coeff(2)) / ( mass_ratio_heart_sap * dia_coeff(1) )

          !! 5.2.3 - Compute bm_sapl for leaf: reverse the original equation
          !!    csa_sap = bm_sapl(j,ileaf,icarbon) / ( pipe_k1 / sla(j)
          !!    Note: before - bm_sapl(j,ileaf,icarbon) = &
          !!         &     (((bm_sapl_leaf(1)*pipe_tune1*(mass_ratio_heart_sap *bm_sapl_leaf(2)*sla(j)/(pi*pipe_k1)) &
          !!         &     **bm_sapl_leaf(3))/sla(j))**bm_sapl_leaf(4))
          IF ( is_tree(j) ) THEN
              bm_sapl(j,ileaf,icarbon) = csa_sap * pipe_k1 / sla(j)
              alpha = alpha_tree
          ELSE
              bm_sapl(j,ileaf,icarbon) = csa_sap * pipe_k1_shrub / sla(j)
              alpha = alpha_shrub
          ENDIF

          !! 5.2.4 - Compute bm_sapl for icarbres (like original)
          IF ( pheno_type(j) .NE. 1 ) THEN
             ! not evergreen
             bm_sapl(j,icarbres,icarbon) = bm_sapl_carbres * bm_sapl(j,ileaf,icarbon)
          ELSE
             bm_sapl(j,icarbres,icarbon) = zero
          ENDIF ! (pheno_type_tab(j) .NE. 1 )

          !! 5.2.4 - Compute bm_sapl for isapabove (quite like original.. exept for shrub because of "pipe_tune2 * dia ** pipe_tune3"= height)
          IF ( is_tree(j) ) THEN  !! Arsene 03-08-2015
             bm_sapl(j,isapabove,icarbon) = &
                & bm_sapl_sapabove * pipe_density * csa_sap * pipe_tune2 * dia ** pipe_tune3

          ELSEIF ( shrubs_like_trees) THEN  !! Arsene 03-08-2015 - Change pipe_tune for shrubs
             bm_sapl(j,isapabove,icarbon) = &  !! Arsene 03-08-2015 - Change pipe_tune for shrubs
                & bm_sapl_sapabove * pipe_density_shrub * csa_sap * pipe_tune2_for_shrub * dia ** pipe_tune3_for_shrub  !! Arsene 03-08-2015 - Change pipe_tune for shrubs


          ELSE
             bm_sapl(j,isapabove,icarbon) = &
                & bm_sapl_sapabove * pipe_density_shrub * csa_sap * &
                & height_presc(j) * pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia**(pipe_tune_shrub3) &
                                & / ( height_presc(j) + pipe_tune_shrub2 * 100**pipe_tune_shrub3 * dia**pipe_tune_shrub3 )


          ENDIF
!! Arsene 11-08-2015 - Change all, beginning with dia = mindia, and reverse some équations...

          !! 5.2.6 - Compute bm_sapl for other "i" (original)
          bm_sapl(j,isapbelow,icarbon) = bm_sapl(j,isapabove,icarbon)

          bm_sapl(j,iheartabove,icarbon) = bm_sapl_heartabove * bm_sapl(j,isapabove,icarbon)
          bm_sapl(j,iheartbelow,icarbon) = bm_sapl_heartbelow * bm_sapl(j,isapbelow,icarbon)

       ELSE

          !> 5.2 grasses

          !!\latexonly
          !!\input{stomate_data_sapl_grass.tex}
          !!\endlatexonly

          alpha = alpha_grass

          IF ( natural(j) .AND. vascular(j)) THEN       !! Arsene 03-03-2015
             bm_sapl(j,ileaf,icarbon) = init_sapl_mass_leaf_nat / sla(j)
          ELSEIF ( .NOT. vascular(j) ) THEN             !! Arsene 05-03-14
             bm_sapl(j,ileaf,icarbon) = init_sapl_mass_leaf_novasc / sla(j)     !! Arsene 18-08-15 - = .05 ?, cf src_parameters/constantes_var.f90
          ELSE
             bm_sapl(j,ileaf,icarbon) = init_sapl_mass_leaf_agri / sla(j)
          ENDIF

          IF ( vascular(j) ) THEN                   !! Arsene 05-03-14
             bm_sapl(j,icarbres,icarbon) = init_sapl_mass_carbres *bm_sapl(j,ileaf,icarbon)
          ELSE                                      !! Arsene 05-03-14
             bm_sapl(j,icarbres,icarbon) = zero     !! Arsene 05-03-14 Pas de reserves
          ENDIF                                     !! Arsene 05-03-14

          bm_sapl(j,isapabove,icarbon) = zero
          bm_sapl(j,isapbelow,icarbon) = zero

          bm_sapl(j,iheartabove,icarbon) = zero
          bm_sapl(j,iheartbelow,icarbon) = zero

       ENDIF !( is_tree(j) )

       IF ( vascular(j) ) THEN                   !! Arsene 05-03-14
          bm_sapl(j,iroot,icarbon) = init_sapl_mass_root * (1./alpha) * bm_sapl(j,ileaf,icarbon)
          bm_sapl(j,ifruit,icarbon) = init_sapl_mass_fruit  * bm_sapl(j,ileaf,icarbon)   !! init_sapl_mass_fruit = 0.3
       ELSE                                      !! Arsene 05-03-14
          bm_sapl(j,iroot,icarbon) = zero        !! Arsene 05-03-14
          bm_sapl(j,ifruit,icarbon) = 0.05  * bm_sapl(j,ileaf,icarbon)  !! Arsene 05-03-14 On prend init_sapl_mass_fruit=0.05 à vérif, cf src_parameters/constantes_var.f90
       ENDIF                                     !! Arsene 05-03-14

       IF ( bavard .GE. 1 ) THEN
          WRITE(numout,*) '       sapling biomass (gC):'
          WRITE(numout,*) '         leaves: (::bm_sapl(j,ileaf,icarbon))',bm_sapl(j,ileaf,icarbon)
          WRITE(numout,*) '         sap above ground: (::bm_sapl(j,ispabove,icarbon)):',bm_sapl(j,isapabove,icarbon)
          WRITE(numout,*) '         sap below ground: (::bm_sapl(j,isapbelow,icarbon))',bm_sapl(j,isapbelow,icarbon)
          WRITE(numout,*) '         heartwood above ground: (::bm_sapl(j,iheartabove,icarbon))',bm_sapl(j,iheartabove,icarbon)
          WRITE(numout,*) '         heartwood below ground: (::bm_sapl(j,iheartbelow,icarbon))',bm_sapl(j,iheartbelow,icarbon)
          WRITE(numout,*) '         roots: (::bm_sapl(j,iroot,icarbon))',bm_sapl(j,iroot,icarbon)
          WRITE(numout,*) '         fruits: (::bm_sapl(j,ifruit,icarbon))',bm_sapl(j,ifruit,icarbon)
          WRITE(numout,*) '         carbohydrate reserve: (::bm_sapl(j,icarbres,icarbon))',bm_sapl(j,icarbres,icarbon)
       ENDIF !( bavard .GE. 1 ) 

       !
       ! 6 migration speed (m/year)
       !

       IF ( is_tree(j) ) THEN

          migrate(j) = migrate_tree

       ELSEIF ( is_shrub(j) ) THEN               !! Arsene 31-07-2014 modifications
          migrate(j) = migrate_shrub             !! Arsene 31-07-2014 modifications

       ELSE

          ! can be any value as grasses are, per *definition*, everywhere (big leaf).
          migrate(j) = migrate_grass

       ENDIF !( is_tree(j) )

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       migration speed (m/year): (::migrate(j))', migrate(j)

       !
!*       ! 7 critical stem diameter: beyond this diameter, the crown area no longer
       !     increases (m)
       !

!*       IF ( is_tree(j) ) THEN       !! Arsene 31-07-2014 modifications ATTENTION. A vérif

          !!\latexonly
          !!\input{stomate_data_stem_diameter.tex}
          !!\endlatexonly

!! Arsene 03-09-2014 - It's possible to define maxdia with - TreeHeight=pipe_tune2 * Diameter^{pipe_tune3} ==> Definition in stomate_establish.f90 or lpj_crown.f90
!           maxdia(j) = ( height_presc(j) / pipe_tune2 ) **(1/ pipe_tune3)

!!          maxdia(j) = ( ( pipe_tune4 / ((pipe_tune2*pipe_tune3)/(maxdia_coeff(1)**pipe_tune3)) ) &     !! Arsene 11-08-2015 BEFORE CHANGES - NOTE: Constant, not PFT dependent...
!!               ** ( un / ( pipe_tune3 - un ) ) ) * maxdia_coeff(2)                                     !! Arsene 11-08-2015 BEFORE CHANGES - Note: maxdia_coeff & pipe_tune4 not use now
!! Arsene 30-03-2015 ==> I think better to change, but not yet !!!!!!!!!!!!!!!!!!! ==> !! Arsene 11-08-2015 - YET !

!*          maxdia(j) = ( height_presc(j) / pipe_tune2 ) **(1/ pipe_tune3)         !! Arsene 11-08-2015 - New ! - Note:  height_presc(j) = max height
!*          mindia(j)= (height_presc(j)/(fact_min_height*pipe_tune2))**(1/ pipe_tune3)

!          maxdia(j) = ( ( pipe_tune4 / ((pipe_tune2*pipe_tune3)/(maxdia_coeff(1)**pipe_tune3)) ) &     !! Arsene 11-08-2015 BEFORE CHANGES - NOTE: Constant, not PFT dependent...
!               ** ( un / ( pipe_tune3 - un ) ) ) * maxdia_coeff(2)

!*          cn_sapl(j) = cn_sapl_init !crown of individual tree, first year

!*       ELSEIF ( is_shrub(j) .AND. shrubs_like_trees ) THEN  !! Arsene 03-08-2015 - Change pipe_tune for shrubs
!          maxdia(j) = ( ( pipe_tune4 / ((pipe_tune2_for_shrub*pipe_tune3_for_shrub)/ &    !! Arsene 04-08-2015 - No pipe_tune_4_shrub
!              & (maxdia_coeff(1)**pipe_tune3_for_shrub)) ) ** ( un / ( pipe_tune3_for_shrub - un ) ) ) * maxdia_coeff(2)

!*          maxdia(j) = ( height_presc(j) / pipe_tune2_for_shrub ) **(1/ pipe_tune3_for_shrub)
!*          mindia(j)= (height_presc(j)/(fact_min_height*pipe_tune2_for_shrub))**(1/ pipe_tune3_for_shrub)

!*          cn_sapl(j) = cn_sapl_init !crown of individual tree, first year (like trees)

!*       ELSEIF ( is_shrub(j) .AND. .NOT.shrubs_like_trees ) THEN  !! Arsene 03-08-2015 - Change pipe_tune for shrubs

!*          maxdia(j) = ( height_presc(j) * 0.93 / (pipe_tune_shrub2*0.07) )**(1/pipe_tune_shrub3) / 100 !! Arsene 11-08-2015 - 0.93 correspond au observatio [90 - 96]
!*          mindia(j) = ( height_presc(j) / (pipe_tune_shrub2 * (fact_min_height-1)) )**(1/pipe_tune_shrub3) / 100

!*       ELSE
!*          maxdia(j) = undef
!*          cn_sapl(j)=1

!*       ENDIF !( is_tree(j) )
!*       IF ( bavard .GE. 1 ) WRITE(numout,*) '       critical stem diameter (m): (::maxdia(j))', maxdia(j)

       !
       ! 8 Coldest tolerable temperature (K)
       !

       IF ( ABS( tmin_crit(j) - undef ) .GT. min_stomate ) THEN
          tmin_crit(j) = tmin_crit(j) + ZeroCelsius
       ELSE
          tmin_crit(j) = undef
       ENDIF 

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       coldest tolerable temperature (K): (::tmin_crit(j))', tmin_crit(j)

       !
       ! 9 Maximum temperature of the coldest month: need to be below this temperature
       !      for a certain time to regrow leaves next spring *(vernalization)* (K)
       !

       IF ( ABS ( tcm_crit(j) - undef ) .GT. min_stomate ) THEN
          tcm_crit(j) = tcm_crit(j) + ZeroCelsius
       ELSE
          tcm_crit(j) = undef
       ENDIF

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       vernalization temperature (K): (::tcm_crit(j))', tcm_crit(j)

       !
       ! 10 critical values for phenology
       !

       ! 10.1 model used

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       phenology model used: (::pheno_model(j)) ',pheno_model(j)

       ! 10.2 growing degree days. What kind of gdd is meant (i.e. threshold 0 or -5 deg C
       !        or whatever), depends on how this is used in stomate_phenology.


       IF ( ( bavard .GE. 1 ) .AND. ( ALL(pheno_gdd_crit(j,:) .NE. undef) ) ) THEN
          WRITE(numout,*) '         critical GDD is a function of long term T (C): (::gdd)'
          WRITE(numout,*) '          ',pheno_gdd_crit(j,1), &
               ' + T *',pheno_gdd_crit(j,2), &
               ' + T^2 *',pheno_gdd_crit(j,3)
       ENDIF

       ! consistency check

       IF ( ( ( pheno_model(j) .EQ. 'moigdd' ) .OR. &
            ( pheno_model(j) .EQ. 'humgdd' )       ) .AND. &
            ( ANY(pheno_gdd_crit(j,:) .EQ. undef) )                      ) THEN
          STOP 'problem with phenology parameters, critical GDD. (::pheno_model)'
       ENDIF

       ! 10.3 number of growing days

       IF ( ( bavard .GE. 1 ) .AND. ( ngd_crit(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical NGD: (::ngd_crit(j))', ngd_crit(j)

       ! 10.4 critical temperature for ncd vs. gdd function in phenology (C)

       IF ( ( bavard .GE. 1 ) .AND. ( ncdgdd_temp(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical temperature for NCD vs. GDD (C): (::ncdgdd_temp(j))', &
            ncdgdd_temp(j)

       ! 10.5 humidity fractions (0-1, unitless)

       IF ( ( bavard .GE. 1 ) .AND. ( hum_frac(j) .NE. undef ) ) &
            WRITE(numout,*) '         critical humidity fraction: (::hum_frac(j))', &
            &  hum_frac(j)


       ! 10.6 minimum time elapsed since moisture minimum (days)

       IF ( ( bavard .GE. 1 ) .AND. ( hum_min_time(j) .NE. undef ) ) &
            WRITE(numout,*) '         time to wait after moisture min (d): (::hum_min_time(j))', &
        &    hum_min_time(j)

       !
       ! 11 critical values for senescence
       !

       ! 11.1 type of senescence

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       type of senescence: (::senescence_type(j))',senescence_type(j)

       ! 11.2 critical temperature for senescence (C)

       IF ( ( bavard .GE. 1 ) .AND. ( ALL(senescence_temp(j,:) .NE. undef) ) ) THEN
          WRITE(numout,*) '         critical temperature for senescence (C) is'
          WRITE(numout,*) '          a function of long term T (C): (::senescence_temp)'
          WRITE(numout,*) '          ',senescence_temp(j,1), &
               ' + T *',senescence_temp(j,2), &
               ' + T^2 *',senescence_temp(j,3)
       ENDIF

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'cold' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )      ) .AND. &
            ( ANY(senescence_temp(j,:) .EQ. undef ) )           ) THEN
          STOP 'problem with senescence parameters, temperature. (::senescence_type)'
       ENDIF

       ! 11.3 critical relative moisture availability for senescence

       IF ( ( bavard .GE. 1 ) .AND. ( senescence_hum(j) .NE. undef ) ) &
            WRITE(numout,*)  ' max. critical relative moisture availability for' 
            WRITE(numout,*)  ' senescence: (::senescence_hum(j))',  &
            & senescence_hum(j)

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'dry' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )     ) .AND. &
            ( senescence_hum(j) .EQ. undef )                   ) THEN
          STOP 'problem with senescence parameters, humidity.(::senescence_type)'
       ENDIF

       ! 14.3 relative moisture availability above which there is no moisture-related
       !      senescence (0-1, unitless)

       IF ( ( bavard .GE. 1 ) .AND. ( nosenescence_hum(j) .NE. undef ) ) &
            WRITE(numout,*) '         relative moisture availability above which there is' 
            WRITE(numout,*) '             no moisture-related senescence: (::nosenescence_hum(j))', &
            &  nosenescence_hum(j)

       ! consistency check

       IF ( ( ( senescence_type(j) .EQ. 'dry' ) .OR. &
            ( senescence_type(j) .EQ. 'mixed' )     ) .AND. &
            ( nosenescence_hum(j) .EQ. undef )                   ) THEN
          STOP 'problem with senescence parameters, humidity. (::senescence_type)'
       ENDIF

       !
       ! 12 sapwood -> heartwood conversion time (days)
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       sapwood -> heartwood conversion time (d): (::tau_sap(j))', tau_sap(j)

       !
       ! 13 fruit lifetime (days)
       !

       IF ( bavard .GE. 1 ) WRITE(numout,*) '       fruit lifetime (d): (::tau_fruit(j))', tau_fruit(j)

       !
       ! 14 length of leaf death (days)
       !      For evergreen trees, this variable determines the lifetime of the leaves.
       !      Note that it is different from the value given in leaflife_tab.
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       length of leaf death (d): (::leaffall(j))', leaffall(j)

       !
       ! 15 maximum lifetime of leaves (days)
       !

       IF ( ( bavard .GE. 1 ) .AND. ( leafagecrit(j) .NE. undef ) ) &
            WRITE(numout,*) '       critical leaf age (d): (::leafagecrit(j))', leafagecrit(j)

       !
       ! 16 time constant for leaf age discretisation (days)
       !

       leaf_timecst(j) = leafagecrit(j) / REAL( nleafages,r_std )

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       time constant for leaf age discretisation (d): (::leaf_timecst(j))', &
            leaf_timecst(j)

       !
       ! 17 minimum lai, initial (m^2.m^{-2})
       !

       IF ( is_tree(j) ) THEN
          lai_initmin(j) = lai_initmin_tree
       ELSEIF ( is_shrub(j) ) THEN                     !! Arsene 31-07-2014 modifications
          lai_initmin(j) = lai_initmin_shrub           !! Arsene 31-07-2014 modifications
       ELSE
          lai_initmin(j) = lai_initmin_grass
       ENDIF !( is_tree(j) )

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       initial LAI: (::lai_initmin(j))', lai_initmin(j)

       !
       ! 19 maximum LAI (m^2.m^{-2})
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       critical LAI above which no leaf allocation: (::lai_max(j))', lai_max(j)

       !
       ! 20 fraction of primary leaf and root allocation put into reserve (0-1, unitless)
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       reserve allocation factor: (::ecureuil(j))', ecureuil(j)

       !
       ! 21 maintenance respiration coefficient (g/g/day) at 0 deg C
       !

       IF ( bavard .GE. 1 ) THEN

          WRITE(numout,*) '       maintenance respiration coefficient (g/g/day) at 0 deg C:'
          WRITE(numout,*) '         . leaves: (::coeff_maint_zero(j,ileaf))',coeff_maint_zero(j,ileaf)
          WRITE(numout,*) '         . sapwood above ground: (::coeff_maint_zero(j,isapabove)) ',&
                        & coeff_maint_zero(j,isapabove)
          WRITE(numout,*) '         . sapwood below ground: (::coeff_maint_zero(j,isapbelow))  ',&
                       & coeff_maint_zero(j,isapbelow)
          WRITE(numout,*) '         . heartwood above ground: (::coeff_maint_zero(j,iheartabove)) ',&
                       & coeff_maint_zero(j,iheartabove)
          WRITE(numout,*) '         . heartwood below ground: (::coeff_maint_zero(j,iheartbelow)) ',&
                       & coeff_maint_zero(j,iheartbelow)
          WRITE(numout,*) '         . roots: (::coeff_maint_zero(j,iroot))',coeff_maint_zero(j,iroot)
          WRITE(numout,*) '         . fruits: (::coeff_maint_zero(j,ifruit)) ',coeff_maint_zero(j,ifruit)
          WRITE(numout,*) '         . carbohydrate reserve: (::coeff_maint_zero(j,icarbres)) ',&
                       & coeff_maint_zero(j,icarbres)

       ENDIF !( bavard .GE. 1 )

       !
       ! 22 parameter for temperature sensitivity of maintenance respiration
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       temperature sensitivity of maintenance respiration (1/K) is'
       WRITE(numout,*) '          a function of long term T (C): (::maint_resp_slope)'
       WRITE(numout,*) '          ',maint_resp_slope(j,1),' + T *',maint_resp_slope(j,2), &
            ' + T^2 *',maint_resp_slope(j,3)

       !
       ! 23 natural ?
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       Natural: (::natural(j))', natural(j)

       !
       ! 24 Vcmax et Vjmax (umol.m^{-2}.s^{-1}) 
       !

       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       Maximum rate of carboxylation: (::Vcmax_25(j))', vcmax25(j)
       !
       ! 25 constants for photosynthesis temperatures
       !

       IF ( bavard .GE. 1 ) THEN


          !
          ! 26 Properties
          !

          WRITE(numout,*) '       C4 photosynthesis: (::is_c4(j))', is_c4(j)
          WRITE(numout,*) '       Depth constant for root profile (m): (::1./humcste(j))', 1./humcste(j)

       ENDIF !( bavard .GE. 1 ) 

       !
       ! 27 extinction coefficient of the Monsi and Saeki (1953) relationship 
       !
       IF ( bavard .GE. 1 ) THEN
          WRITE(numout,*) '       extinction coefficient: (::ext_coeff(j))', ext_coeff(j)
       ENDIF !( bavard .GE. 1 )

       !
       ! 30 fraction of allocatable biomass which is lost as growth respiration (0-1, unitless)
       !
       IF ( bavard .GE. 1 ) &
            WRITE(numout,*) '       growth respiration fraction: (::frac_growthresp(j))', frac_growthresp(j)

    ENDDO ! Loop over # PFTS 

    !
    ! 29 time scales for phenology and other processes (in days)
    !

    tau_longterm = coeff_tau_longterm * one_year

    IF ( bavard .GE. 1 ) THEN

       WRITE(numout,*) '   > time scale for ''monthly'' moisture availability (d): (::tau_hum_month)', &
            tau_hum_month
       WRITE(numout,*) '   > time scale for ''weekly'' moisture availability (d): (::tau_hum_week)', &
           tau_hum_week
       WRITE(numout,*) '   > time scale for ''monthly'' 2 meter temperature (d): (::tau_t2m_month)', &
            tau_t2m_month
       WRITE(numout,*) '   > time scale for ''weekly'' 2 meter temperature (d): (::tau_t2m_week)', &
            tau_t2m_week
       WRITE(numout,*) '   > time scale for ''weekly'' GPP (d): (::tau_gpp_week)', &
            tau_gpp_week
       WRITE(numout,*) '   > time scale for ''monthly'' soil temperature (d): (::tau_tsoil_month)', &
            tau_tsoil_month
       WRITE(numout,*) '   > time scale for ''monthly'' soil humidity (d): (::tau_soilhum_month)', &
            tau_soilhum_month
       WRITE(numout,*) '   > time scale for vigour calculations (y): (::tau_longterm / one_year)', &
            tau_longterm / one_year

    ENDIF !( bavard .GE. 1 ) 

    IF (bavard.GE.4) WRITE(numout,*) 'Leaving data'

  END SUBROUTINE data

END MODULE stomate_data
