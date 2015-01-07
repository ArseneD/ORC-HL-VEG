!
MODULE extrac_cti
!
!       This routine compute the topographic index in each grid cell of the
!       studied domain using the file from the HYDRO1K dataset (CTI). Two variables 
!       must be change when the domain change: NCELL and ZRES
!
        USE ioipsl
        USE constantes
        USE mod_orchidee_para
!!        USE reqdprec
!
        IMPLICIT NONE
        INCLUDE "netcdf.inc"
        INTEGER, PARAMETER              :: ISIMPLE=KIND(1.0)   ! Simple precision
        INTEGER, PARAMETER              :: IDOUBLE=KIND(1.D0)  ! Double precision
!       RADIUS OF INFLUENCE
!
        REAL(r_std), PARAMETER   :: ZRADIUS=6370.997
        REAL(r_std), PARAMETER   :: ZRADIUS_AF=6378.137
!        
!       MISSING VALUES 
!        
        INTEGER, PARAMETER              :: IMISSING=-9999
        REAL(r_std), PARAMETER   :: ZMISSING=-99.99
!
CONTAINS
        SUBROUTINE extrac_cti_main(kjpindex, lalo, resolution, ZMIN, ZMAX, ZMEAN, ZSTDT, ZSKEW)
        INTEGER(i_std), INTENT(in)                         :: kjpindex         !! Domain size
        REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo       !! Geogr. coordinates (latitude,longitude) (degrees)
        REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution !! size in x an y of the grid (m)
        REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: ZMIN
        REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: ZMAX
        REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: ZMEAN
        REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: ZSTDT
        REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: ZSKEW
!        REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: NB_PIXE

        CHARACTER(LEN=100) :: AF_CTI_file, AS_CTI_file, AU_CTI_file, EU_CTI_file, NA_CTI_file, SA_CTI_file
        CHARACTER(LEN=100) :: TOPM_DATA_File

!       AFRIQUE
!
        INTEGER, PARAMETER              :: NCOL_AF=8736
        INTEGER, PARAMETER              :: NROW_AF=9194 
        INTEGER, PARAMETER              :: NPIXEL_AF=NCOL_AF*NROW_AF
!
        INTEGER, PARAMETER              :: IX_FIRST_AF=-4368
        INTEGER, PARAMETER              :: IY_FIRST_AF=4149
!
        REAL(r_std), PARAMETER   :: ZLON_AF=20.
        REAL(r_std), PARAMETER   :: ZLAT_AF=5.
!
!       ASIE
!        
        INTEGER, PARAMETER              :: NCOL_AS=9341
        INTEGER, PARAMETER              :: NROW_AS=11882 
        INTEGER, PARAMETER              :: NPIXEL_AS=NCOL_AS*NROW_AS
!
        INTEGER, PARAMETER              :: IX_FIRST_AS=-4355
        INTEGER, PARAMETER              :: IY_FIRST_AS=6443
!
        REAL(r_std), PARAMETER   :: ZLON_AS=100.
        REAL(r_std), PARAMETER   :: ZLAT_AS=45.
!
!       AUSTRALIE
!        
        INTEGER, PARAMETER              :: NCOL_AU=15142
        INTEGER, PARAMETER              :: NROW_AU=11396 
        INTEGER, PARAMETER              :: NPIXEL_AU=NCOL_AU*NROW_AU
!
        INTEGER, PARAMETER              :: IX_FIRST_AU=-5219
        INTEGER, PARAMETER              :: IY_FIRST_AU=5307
!
        REAL(r_std), PARAMETER   :: ZLON_AU=135.
        REAL(r_std), PARAMETER   :: ZLAT_AU=-15.
!
!       EUROPE
!        
        INTEGER, PARAMETER              :: NCOL_EU=8319
        INTEGER, PARAMETER              :: NROW_EU=7638 
        INTEGER, PARAMETER              :: NPIXEL_EU=NCOL_EU*NROW_EU
!
        INTEGER, PARAMETER              :: IX_FIRST_EU=-4091
        INTEGER, PARAMETER              :: IY_FIRST_EU=3293
!
        REAL(r_std), PARAMETER   :: ZLON_EU=20.
        REAL(r_std), PARAMETER   :: ZLAT_EU=55.
!
!       AMERIQUE DU NORD
!        
        INTEGER, PARAMETER              :: NCOL_NA=9102
        INTEGER, PARAMETER              :: NROW_NA=8384 
        INTEGER, PARAMETER              :: NPIXEL_NA=NCOL_NA*NROW_NA
!
        INTEGER, PARAMETER              :: IX_FIRST_NA=-4462
        INTEGER, PARAMETER              :: IY_FIRST_NA=4384
!
        REAL(r_std), PARAMETER   :: ZLON_NA=-100.
        REAL(r_std), PARAMETER   :: ZLAT_NA=45.
!
!       AMERIQUE DU SUD
!        
        INTEGER, PARAMETER              :: NCOL_SA=7736
        INTEGER, PARAMETER              :: NROW_SA=9094 
        INTEGER, PARAMETER              :: NPIXEL_SA=NCOL_SA*NROW_SA
!
        INTEGER, PARAMETER              :: IX_FIRST_SA=-3776
        INTEGER, PARAMETER              :: IY_FIRST_SA=3835
!
        REAL(r_std), PARAMETER   :: ZLON_SA=-60.
        REAL(r_std), PARAMETER   :: ZLAT_SA=-15.
!
!        
!       TABLEAU DE LECTURE DES INDICES (1 PAR CONTINENT)
!
        INTEGER, DIMENSION(NPIXEL_AF) :: ICTI_READ_AF
        INTEGER, DIMENSION(NPIXEL_AS) :: ICTI_READ_AS
        INTEGER, DIMENSION(NPIXEL_AU) :: ICTI_READ_AU
        INTEGER, DIMENSION(NPIXEL_EU) :: ICTI_READ_EU
        INTEGER, DIMENSION(NPIXEL_NA) :: ICTI_READ_NA
        INTEGER, DIMENSION(NPIXEL_SA) :: ICTI_READ_SA
!
        INTEGER, DIMENSION(NPIXEL_AF)  :: ICTI_AF
        INTEGER, DIMENSION(NPIXEL_AS)  :: ICTI_AS
        INTEGER, DIMENSION(NPIXEL_AU)  :: ICTI_AU
        INTEGER, DIMENSION(NPIXEL_EU)  :: ICTI_EU
        INTEGER, DIMENSION(NPIXEL_NA)  :: ICTI_NA
        INTEGER, DIMENSION(NPIXEL_SA)  :: ICTI_SA
!
        REAL(r_std)                               :: XDIST, YDIST, xx_lon, yy_lat
        REAL(r_std), DIMENSION(kjpindex)          :: ZLON_G,ZLAT_G,ZMEAN2
!        REAL(r_std), DIMENSION(kjpindex,80000)    :: ZCTI
        REAL(r_std), DIMENSION(kjpindex,320000)    :: ZCTI
!        REAL(r_std), DIMENSION(kjpindex,3000000)    :: ZCTI
        INTEGER,    DIMENSION(kjpindex)          :: INBR, DISTANCE
        INTEGER                                  :: II, III,IV(1),IV_lev,NCELL
        INTEGER :: err, du_na, du_sa, du_af, du_as, du_au, du_eu
        INTEGER :: afstanid, asstanid, austanid, eustanid
        INTEGER :: nastanid, sastanid
        INTEGER :: afctiid, asctiid, auctiid, euctiid
        INTEGER :: nactiid, sactiid
!
PRINT*,'INITIALISATION !'
write(*,*) 'LECTURE DE EXTRAC_CTI!!!!!'
!
ZCTI(:,:) = ZMISSING
!
INBR(:)=1
NCELL=kjpindex
!
!xx_lon=7.5
xx_lon=1.0
CALL getin_p('ZONAL_RES',xx_lon)
!yy_lat=5.0
yy_lat=1.0
CALL getin_p('MERID_RES',yy_lat)
write(*,*) "lecture resol extrac_cti", xx_lon, yy_lat
YDIST = yy_lat/2.0
XDIST = xx_lon/2.0
!YDIST = resolution(1,1)/2.0 
!XDIST = resolution(1,2)/2.0
write(*,*) 'xdist', XDIST, 'ydist', YDIST
!
!LECTURE DES FICHIERS BINAIRE
!
PRINT*,''
PRINT*,'LECTURE DES FICHIERS BINAIRE !'
!

!pss:+
!IF (is_root_prc) THEN
AF_CTI_file = 'af_cti.nc'
CALL getin_p('AF_CTI', AF_CTI_file)
du_af = NF_OPEN(AF_CTI_file, nf_nowrite, afstanid)
!pss:-

write(*,*)NF_strerror(du_af)
du_af = nf_inq_varid(afstanid, "CTI", afctiid)
write(*,*)NF_strerror(du_af)
du_af = NF_GET_VAR_INT(afstanid, afctiid, ICTI_READ_AF)
write(*,*)NF_strerror(du_af)
du_af = NF_ClOSE(afstanid)
write(*,*)NF_strerror(du_af)
ICTI_AF(:)=ICTI_READ_AF(:)
!

!pss:+
AS_CTI_file = 'as_cti.nc'
CALL getin_p('AS_CTI', AS_CTI_file)
du_as = NF_OPEN(AS_CTI_file, nf_nowrite, afstanid)
!pss:-

!du_as = NF_OPEN('../../as_cti.nc', nf_nowrite, asstanid)
du_as = nf_inq_varid(asstanid, "CTI", asctiid)
du_as = NF_GET_VAR_INT(asstanid, asctiid, ICTI_READ_AS)
du_as = NF_ClOSE(asstanid)
ICTI_AS(:)=ICTI_READ_AS(:)
!

!pss:+
AU_CTI_file = 'au_cti.nc'
CALL getin_p('AU_CTI', AU_CTI_file)
du_au = NF_OPEN(AU_CTI_file, nf_nowrite, afstanid)
!pss:-
!du_au = NF_OPEN('../../au_cti.nc', nf_nowrite, austanid)
du_au = nf_inq_varid(austanid, "CTI", auctiid)
du_au = NF_GET_VAR_INT(austanid, auctiid, ICTI_READ_AU)
du_au = NF_ClOSE(austanid)
ICTI_AU(:)=ICTI_READ_AU(:)
!

!pss:+
EU_CTI_file = 'eu_cti.nc'
CALL getin_p('EU_CTI', EU_CTI_file)
du_eu = NF_OPEN(EU_CTI_file, nf_nowrite, afstanid)
!pss:-
!du_eu = NF_OPEN('../../eu_cti.nc', nf_nowrite, eustanid)
du_eu = nf_inq_varid(eustanid, "CTI", euctiid)
du_eu = NF_GET_VAR_INT(eustanid, euctiid, ICTI_READ_EU)
du_eu = NF_ClOSE(eustanid)
ICTI_EU(:)=ICTI_READ_EU(:)
!

!pss:+
NA_CTI_file = 'na_cti.nc'
CALL getin_p('NA_CTI', NA_CTI_file)
du_na = NF_OPEN(NA_CTI_file, nf_nowrite, afstanid)
!pss:-
!du_na = NF_OPEN('../../na_cti.nc', nf_nowrite, nastanid)
du_na = nf_inq_varid(nastanid, "CTI", nactiid)
du_na = NF_GET_VAR_INT(nastanid, nactiid, ICTI_READ_NA)
du_na = NF_ClOSE(nastanid)
ICTI_NA(:)=ICTI_READ_NA(:)
!

!pss:+
SA_CTI_file = 'sa_cti.nc'
CALL getin_p('SA_CTI', SA_CTI_file)
du_sa = NF_OPEN(SA_CTI_file, nf_nowrite, afstanid)
!pss:-
!du_sa = NF_OPEN('../../sa_cti.nc', nf_nowrite, sastanid)
du_sa = nf_inq_varid(sastanid, "CTI", sactiid)
du_sa = NF_GET_VAR_INT(sastanid, sactiid, ICTI_READ_SA)
du_sa = NF_ClOSE(sastanid)
ICTI_SA(:)=ICTI_READ_SA(:)
!
PRINT*,MAXVAL(ICTI_SA(:))
PRINT*,MINVAL(ICTI_SA(:))

!
!LECTURE DU MASK
!
PRINT*,''
PRINT*,'!!! READ LATITUDE LONGITUDE OF THE DOMAIN !!!'
PRINT*,''
!
ZLAT_G(:)=lalo(:,1)
ZLON_G(:)=lalo(:,2)
!
!LECTURE DES INDICES ET STATS
!
PRINT*,''
PRINT*,'!!! CALCUL DES COORDONÉES DE CHAQUE MAILLE,   !!!'
PRINT*,'!!! EXTRACTION DES PIXELS CONTENU DANS CHAQUE !!!'
PRINT*,'!!! MAILLE ET CALCUL DES STATISTIQUES         !!!'
PRINT*,''
!
DO II=1,NCELL
!
   PRINT*,''
   PRINT*,'MAILLE :',II
   PRINT*,''
!
!  TROUVER LES PIXELS CONTENU DANS CHAQUE MAILLE
!
   IF(ZLAT_G(II)<=45.0.AND.ZLAT_G(II)>=-45.0.AND.ZLON_G(II)<=60.0.AND.ZLON_G(II)>=-25.0)THEN
!
   CALL BOUCLE(ICTI_AF,ZCTI(II,:),ZLAT_AF,ZLON_AF,ZLAT_G(II),ZLON_G(II),&
               NROW_AF,NCOL_AF,IY_FIRST_AF,IX_FIRST_AF,XDIST,YDIST,ZRADIUS_AF,INBR(II))
!
   ENDIF
!
   IF(ZLAT_G(II)<=85.0.AND.ZLAT_G(II)>=-10.0.AND.ZLON_G(II)<=180.0.AND.ZLON_G(II)>=50.0)THEN
!
   CALL BOUCLE(ICTI_AS,ZCTI(II,:),ZLAT_AS,ZLON_AS,ZLAT_G(II),ZLON_G(II),&
               NROW_AS,NCOL_AS,IY_FIRST_AS,IX_FIRST_AS,XDIST,YDIST,ZRADIUS,INBR(II))
!
   ENDIF
!   
   IF(ZLAT_G(II)<=85.0.AND.ZLAT_G(II)>=60.0.AND.ZLON_G(II)<=-160.0.AND.ZLON_G(II)>=-180.0)THEN
!
   CALL BOUCLE(ICTI_AS,ZCTI(II,:),ZLAT_AS,ZLON_AS,ZLAT_G(II),ZLON_G(II),&
               NROW_AS,NCOL_AS,IY_FIRST_AS,IX_FIRST_AS,XDIST,YDIST,ZRADIUS,INBR(II))
!
   ENDIF
!
   IF(ZLAT_G(II)<=30.0.AND.ZLAT_G(II)>=-57.0.AND.ZLON_G(II)<=180.0.AND.ZLON_G(II)>=90.0)THEN
!
   CALL BOUCLE(ICTI_AU,ZCTI(II,:),ZLAT_AU,ZLON_AU,ZLAT_G(II),ZLON_G(II),&
               NROW_AU,NCOL_AU,IY_FIRST_AU,IX_FIRST_AU,XDIST,YDIST,ZRADIUS,INBR(II))
!
   ENDIF
!
   IF(ZLAT_G(II)<=85.0.AND.ZLAT_G(II)>=10.0.AND.ZLON_G(II)<=80.0.AND.ZLON_G(II)>=-30.0)THEN
!
   CALL BOUCLE(ICTI_EU,ZCTI(II,:),ZLAT_EU,ZLON_EU,ZLAT_G(II),ZLON_G(II),&
               NROW_EU,NCOL_EU,IY_FIRST_EU,IX_FIRST_EU,XDIST,YDIST,ZRADIUS,INBR(II))
!
   ENDIF
!
   IF(ZLAT_G(II)<=85.0.AND.ZLAT_G(II)>=0.0.AND.ZLON_G(II)<=-45.0.AND.ZLON_G(II)>=-180.0)THEN
!
   CALL BOUCLE(ICTI_NA,ZCTI(II,:),ZLAT_NA,ZLON_NA,ZLAT_G(II),ZLON_G(II),&
               NROW_NA,NCOL_NA,IY_FIRST_NA,IX_FIRST_NA,XDIST,YDIST,ZRADIUS,INBR(II))
!
   ENDIF
!
   IF(ZLAT_G(II)<=15.0.AND.ZLAT_G(II)>=-60.0.AND.ZLON_G(II)<=-30.0.AND.ZLON_G(II)>=-90.0)THEN
!
   CALL BOUCLE(ICTI_SA,ZCTI(II,:),ZLAT_SA,ZLON_SA,ZLAT_G(II),ZLON_G(II),&
               NROW_SA,NCOL_SA,IY_FIRST_SA,IX_FIRST_SA,XDIST,YDIST,ZRADIUS,INBR(II))
!
   ENDIF
!
  PRINT*,'nbr :',(INBR(II)-1)
!
! CALCUL DES STATISTIQUES
!
  CALL STATS(ZCTI(II,:),ZMISSING,ZMEAN(II),ZSTDT(II),ZSKEW(II),ZMIN(II),ZMAX(II))
!
  PRINT*,'mean :',ZMEAN (II), ZCTI(II,1), ZCTI(II,2)
!
ENDDO

!DO II=1,NCELL
!  NB_PIXE(II) = COUNT(ZCTI(II,:) > (ZMISSING+1))
!ENDDO

ZMEAN2(:)=ZMEAN(:)!ZMEAN2 reste constant au cours du temps alors que ZMEAn evolue
!cela evite d avoir des propagations de nouvelles valeurs dans les zones ou plusieurs grid cells sont missing (Asutralie) 

!je cree des vecteurs ou j enleve les grid-cells pour lesquels j ai des valeurs missing
DO II=1,NCELL
   if(ZMEAN2(II).EQ.ZMISSING)then
      !il faut trouver le grid cell le plus proche pour lequel je n ai pas de valeurs missing
      do III=1,NCELL
         DISTANCE(III) =  SQRT((lalo(III,1)-lalo(II,1))**2 + (lalo(III,2)-lalo(II,2))**2)
      enddo
      ! ~vecteur de longeur kjpindex et qui donne la distance avec le point II 

      IV=INT( MINLOC(DISTANCE(:)) )!rang du grid-cell dont la distance avec le grid-cell est la plus proche 
      IV_lev=IV(1)
      do while(ZMEAN2(IV_lev).EQ.ZMISSING)
         DISTANCE(IV_lev)=999999
         IV=INT(MINLOC(DISTANCE(:)))
         IV_lev=IV(1)
      enddo
      ZMEAN(II)=ZMEAN2(IV_lev)
      ZSTDT(II) = ZSTDT(IV_lev)
      ZSKEW(II) = ZSKEW(IV_lev)
      ZMIN(II) = ZMIN(IV_lev)
      ZMAX(II) = ZMAX(IV_lev)     
   endif
ENDDO

!ENDIF

!pss:modify
!OPEN(UNIT=11, FILE='TOPMODEL_G.DAT',  FORM='FORMATTED')
!WRITE(11,'(10F13.5)') ZMIN  (:)
!WRITE(11,'(10F13.5)') ZMAX  (:)
!WRITE(11,'(10F13.5)') ZMEAN (:)
!WRITE(11,'(10F13.5)') ZSTDT (:)
!WRITE(11,'(10F13.5)') ZSKEW (:)
!CLOSE(UNIT=11)
!pss:-
!
!PRINT*,'This is the end !'
END SUBROUTINE extrac_cti_main
!
!---------------------------------------------------------------!
!
!
SUBROUTINE PROJECTION(PLAT0,PLON0,PXDIST,PYDIST,PLAT,PLON,PRADIUS,PX_UR,PY_UR,PX_LL,PY_LL)
!
IMPLICIT NONE
!
!        
        REAL(r_std), PARAMETER   :: XEPSLN=1.0e-10
        REAL(r_std), PARAMETER   :: XPI=3.141592653589793238
!
        REAL(r_std), INTENT(IN)  :: PLAT0,PLON0,PLAT,PLON,PXDIST,PYDIST,PRADIUS
        REAL(r_std), INTENT(OUT) :: PX_UR,PY_UR,PX_LL,PY_LL
!
        REAL(r_std)              :: ZLON,ZLAT,ZLAT0,ZLON0
!        
        REAL(r_std)              :: ZXUR,ZYUR,ZXLL,ZYLL,ZXUL,ZYUL,ZXLR,ZYLR
!        
        REAL(r_std)              :: ZG,ZDLON,ZKSP
!
ZLAT0=PLAT0*XPI/180.
ZLON0=PLON0*XPI/180.
!
! UPPER RIGHT
!
ZLAT=PLAT+PYDIST
ZLON=PLON+PXDIST
!
ZLON=ZLON*XPI/180.
ZLAT=ZLAT*XPI/180.
!
ZDLON=ZLON-PLON0*XPI/180.
!
ZG=SIN(ZLAT0)*SIN(ZLAT)+COS(ZLAT0)*COS(ZLAT)*COS(ZDLON)
!
IF(ZG==-1)PRINT*,'!!! Point projects to a circle of radius = 2*RADIUS !!!'
!
ZKSP=PRADIUS*(2./(1.+ZG))**0.5
!
ZXUR=ZKSP*COS(ZLAT)*SIN(ZDLON)
ZYUR=ZKSP*(COS(ZLAT0)*SIN(ZLAT)-SIN(ZLAT0)*COS(ZLAT)*COS(ZDLON))
!
! LOWER LEFT
!
ZLAT=PLAT-PYDIST
ZLON=PLON-PXDIST
!
ZLON=ZLON*XPI/180.
ZLAT=ZLAT*XPI/180.
!
ZDLON=ZLON-PLON0*XPI/180.
!
ZG=SIN(ZLAT0)*SIN(ZLAT)+COS(ZLAT0)*COS(ZLAT)*COS(ZDLON)
!
IF(ZG==-1)PRINT*,'!!! Point projects to a circle of radius = 2*RADIUS !!!'
!
ZKSP=PRADIUS*(2./(1.+ZG))**0.5
!
ZXLL=ZKSP*COS(ZLAT)*SIN(ZDLON)
ZYLL=ZKSP*(COS(ZLAT0)*SIN(ZLAT)-SIN(ZLAT0)*COS(ZLAT)*COS(ZDLON))
!
! LOWER RIGHT
!
ZLAT=PLAT-PYDIST
ZLON=PLON+PXDIST
!
ZLON=ZLON*XPI/180.
ZLAT=ZLAT*XPI/180.
!
ZDLON=ZLON-PLON0*XPI/180.
!
ZG=SIN(ZLAT0)*SIN(ZLAT)+COS(ZLAT0)*COS(ZLAT)*COS(ZDLON)
!
IF(ZG==-1)PRINT*,'!!! Point projects to a circle of radius = 2*RADIUS !!!'
!
ZKSP=PRADIUS*(2./(1.+ZG))**0.5
!
ZXLR=ZKSP*COS(ZLAT)*SIN(ZDLON)
ZYLR=ZKSP*(COS(ZLAT0)*SIN(ZLAT)-SIN(ZLAT0)*COS(ZLAT)*COS(ZDLON))
!
! UPPER LEFT
!
ZLAT=PLAT+PYDIST
ZLON=PLON-PXDIST
!
ZLON=ZLON*XPI/180.
ZLAT=ZLAT*XPI/180.
!
ZDLON=ZLON-PLON0*XPI/180.
!
ZG=SIN(ZLAT0)*SIN(ZLAT)+COS(ZLAT0)*COS(ZLAT)*COS(ZDLON)
!
IF(ZG==-1)PRINT*,'!!! Point projects to a circle of radius = 2*RADIUS !!!'
!
ZKSP=PRADIUS*(2./(1.+ZG))**0.5
!
ZXUL=ZKSP*COS(ZLAT)*SIN(ZDLON)
ZYUL=ZKSP*(COS(ZLAT0)*SIN(ZLAT)-SIN(ZLAT0)*COS(ZLAT)*COS(ZDLON))
!
! BILAN
!
PX_UR=MAX(ZXUR,ZXLR)
PY_UR=MAX(ZYUR,ZYUL)
PX_LL=MIN(ZXUL,ZXLL)
PY_LL=MIN(ZYLL,ZYLR)
!
END SUBROUTINE PROJECTION
!
!-----------------------------------------------------------!
!
SUBROUTINE BOUCLE(KCTI,PCTI,PLAT,PLON,PLAT_G,PLON_G,NROW,   &
                  NCOL,KY_FIRST,KX_FIRST,PXDIST,PYDIST,PRADIUS,KNBR )
!
IMPLICIT NONE
!
!
        INTEGER, PARAMETER                :: IMISSING=-9999
        INTEGER, DIMENSION(:), INTENT(IN) :: KCTI
!        
        REAL(r_std), DIMENSION(:), INTENT(OUT):: PCTI
!        
        REAL(r_std), INTENT(IN)    :: PLAT,PLON
!        
        REAL(r_std), INTENT(IN)    :: PLAT_G,PLON_G,PXDIST,PYDIST,PRADIUS
!
        INTEGER,            INTENT(INOUT) :: KNBR
!        
        INTEGER,            INTENT(IN)    :: NROW,NCOL,KY_FIRST,KX_FIRST
!
        REAL(r_std)                :: ZX_UR,ZY_UR,ZX_LL,ZY_LL
!
        INTEGER                           :: IX,IY,IK,IR,IC,INBR
!
!CALCUL DES COORDONÉES (X,Y) DE CHAQUE MAILLE
!
CALL PROJECTION(PLAT,PLON,PXDIST,PYDIST,PLAT_G,PLON_G,PRADIUS,ZX_UR,ZY_UR,ZX_LL,ZY_LL)
!
!write(*,*) ZX_UR, ZY_UR, ZX_LL, ZY_LL   
IK=1
!
IY=KY_FIRST
!
DO IR=NROW,1,-1 !BOUCLE SUR LES LIGNES (Y)
!
   IX=KX_FIRST   
!
!   write(*,*) IY, IX, KCTI(IK), KNBR   
   DO IC=1,NCOL  !BOUCLE SUR LES COLONNES (X)
!
      IF(KCTI(IK).GT.-200)THEN
!
      IF(IX<=ZX_UR.AND.IY<=ZY_UR &
      .AND.IX>=ZX_LL.AND.IY>=ZY_LL)THEN
!
!         write(*,*) 'ok', IY, IX, KNBR, KCTI(IK)
         PCTI(KNBR)=(REAL(KCTI(IK))/100.)
!        
         KNBR=KNBR+1
!            
      ENDIF                 
!
      ENDIF                 
!
      IK=IK+1
      IX=IX+1
!
   ENDDO
!
   IY=IY-1
!
ENDDO
!write(*,*) 'ok', KNBR
!
END SUBROUTINE BOUCLE
!
!-----------------------------------------------------------!
!
SUBROUTINE STATS(PVALS,PMISSING,PMEAN,PSTDT,PSKEW,PMIN,PMAX)
!
IMPLICIT NONE
!
!
REAL(r_std), DIMENSION(:),INTENT(INOUT):: PVALS
!
REAL(r_std), INTENT(IN)                :: PMISSING
!
REAL(r_std), INTENT(OUT)               :: PMEAN,PSTDT,PSKEW,PMIN,PMAX
!
REAL(r_std)            :: ZVAR
!
INTEGER         :: INBR,JVAR,II
!
!INITIALISATION
!
JVAR=SIZE(PVALS(:))
!
INBR=COUNT(PVALS(:) > (PMISSING+1))
!
IF(INBR>50)THEN
!
!MIN ET MAX
!
PMIN=MINVAL(PVALS(:),PVALS(:)> (PMISSING+1))
PMAX=MAXVAL(PVALS(:),PVALS(:)> (PMISSING+1))
!
!MEAN VALUE
!
PMEAN=SUM(PVALS(:),PVALS(:)> (PMISSING+1))/INBR
!
! STANDARD DEVIATION
!
ZVAR=SUM((PVALS(:)-PMEAN)**2,PVALS(:)> (PMISSING+1))
ZVAR=ZVAR/(INBR-1)
!
PSTDT=ZVAR**0.5
!
! SKEWNESS
!
PSKEW=SUM((((PVALS(:)-PMEAN)/PSTDT)**3),PVALS(:)> (PMISSING+1))/(INBR-1)
!
ELSE
!
  PMIN =PMISSING
  PMAX =PMISSING
  PMEAN=PMISSING
  PSTDT=PMISSING
  PSKEW=PMISSING
!
ENDIF
!
END SUBROUTINE STATS
!
!-----------------------------------------------------------!
!
END MODULE extrac_cti 
