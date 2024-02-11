!#ifdef _ACCEL
!#  include "module_mp_wsm5_accel.F"
!#else
!#if ( RWORDSIZE == 4 )
!#  define VREC vsrec
!#  define VSQRT vssqrt
!#else
!#  define VREC vrec
!#  define VSQRT vsqrt
!#endif

!Including inline expansion statistical function 
MODULE module_mp_wsm5

!$acc routine(wsm52D)
!$acc routine(refl10cm_wsm5)
!$acc routine(effectRad_wsm5)
!$acc routine(wsm5x)
!$acc routine(slope_wsm5)
!$acc routine(nislfv_rain_plm)
!$acc routine(slope_rain)
!$acc routine(slope_snow)
!$acc routine(vsrec)
!$acc routine(vssqrt)

   USE module_mp_radar
   !-srf
   ! USE module_model_constants, only : RE_QC_BG, RE_QI_BG, RE_QS_BG
   REAL    , PARAMETER :: RE_QC_BG     = 2.51E-6     ! effective radius of cloud for background (m)
   REAL    , PARAMETER :: RE_QI_BG     = 5.01E-6     ! effective radius of ice for background (m)
   REAL    , PARAMETER :: RE_QS_BG     = 10.01E-6    ! effective radius of snow for background (m)
   !-srf
!
   REAL, PARAMETER, PRIVATE :: dtcldcr     = 120. ! maximum time step for minor loops
   REAL, PARAMETER, PRIVATE :: n0r = 8.e6         ! intercept parameter rain
   REAL, PARAMETER, PRIVATE :: avtr = 841.9       ! a constant for terminal velocity of rain
   REAL, PARAMETER, PRIVATE :: bvtr = 0.8         ! a constant for terminal velocity of rain
   REAL, PARAMETER, PRIVATE :: r0 = .8e-5         ! 8 microm  in contrast to 10 micro m
   REAL, PARAMETER, PRIVATE :: peaut = .55        ! collection efficiency
   REAL, PARAMETER, PRIVATE :: xncr = 3.e8        ! maritime cloud in contrast to 3.e8 in tc80
   REAL, PARAMETER, PRIVATE :: xmyu = 1.718e-5    ! the dynamic viscosity kgm-1s-1
   REAL, PARAMETER, PRIVATE :: avts = 11.72       ! a constant for terminal velocity of snow
   REAL, PARAMETER, PRIVATE :: bvts = .41         ! a constant for terminal velocity of snow
   REAL, PARAMETER, PRIVATE :: n0smax =  1.e11    ! maximum n0s (t=-90C unlimited)
   REAL, PARAMETER, PRIVATE :: lamdarmax = 8.e4   ! limited maximum value for slope parameter of rain
   REAL, PARAMETER, PRIVATE :: lamdasmax = 1.e5   ! limited maximum value for slope parameter of snow
   REAL, PARAMETER, PRIVATE :: lamdagmax = 6.e4   ! limited maximum value for slope parameter of graupel
   REAL, PARAMETER, PRIVATE :: dicon = 11.9       ! constant for the cloud-ice diamter
   REAL, PARAMETER, PRIVATE :: dimax = 500.e-6    ! limited maximum value for the cloud-ice diamter
   REAL, PARAMETER, PRIVATE :: n0s = 2.e6         ! temperature dependent intercept parameter snow
   REAL, PARAMETER, PRIVATE :: alpha = .12        ! .122 exponen factor for n0s
   REAL, PARAMETER, PRIVATE :: pfrz1 = 100.       ! constant in Biggs freezing
   REAL, PARAMETER, PRIVATE :: pfrz2 = 0.66       ! constant in Biggs freezing
   REAL, PARAMETER, PRIVATE :: qcrmin = 1.e-9     ! minimun values for qr, qs, and qg
   REAL, PARAMETER, PRIVATE :: eacrc = 1.0        ! Snow/cloud-water collection efficiency
   ! REAL, SAVE ::  !30/11/2023  
	REAL ::     &
             qc0, qck1, pidnc,                        &
             bvtr1,bvtr2,bvtr3,bvtr4,g1pbr,           &
             g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,    &
             precr1,precr2,xmmax,roqimax,bvts1,       &
             bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,     &
             g5pbso2,pvts,pacrs,precs1,precs2,pidn0r, &
             pidn0s,xlv1,pacrc,pi,                    &
             rslopermax,rslopesmax,rslopegmax,        &
             rsloperbmax,rslopesbmax,rslopegbmax,     &
             rsloper2max,rslopes2max,rslopeg2max,     &
             rsloper3max,rslopes3max,rslopeg3max

!newCode begin
!$acc declare create(qc0, qck1, pidnc)
!$acc declare create(bvtr1,bvtr2,bvtr3,bvtr4,g1pbr)
!$acc declare create(g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr)
!$acc declare create(precr1,precr2,xmmax,roqimax,bvts1)
!$acc declare create(bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs)
!$acc declare create(g5pbso2,pvts,pacrs,precs1,precs2,pidn0r)
!$acc declare create(pidn0s,xlv1,pacrc,pi)
!$acc declare create(rslopermax,rslopesmax,rslopegmax)
!$acc declare create(rsloperbmax,rslopesbmax,rslopegbmax)
!$acc declare create(rsloper2max,rslopes2max,rslopeg2max)
!$acc declare create(rsloper3max,rslopes3max,rslopeg3max)
!newCode end

! Specifies code-inlining of fpvs function in WSM52D below. JM 20040507
!
CONTAINS

!===================================================================
!
  SUBROUTINE wsm5(th, q, qc, qr, qi, qs                            & !(private - todas)
                 ,den, pii, p, delz                                & !(private - todas)
                 ,delt, g, cpd, cpv, rd, rv, t0c                   & !(private - somente delt)
                 ,ep1, ep2, qmin                                   &
                 ,XLS, XLV0, XLF0, den0, denr                      &
                 ,cliq,cice,psat                                   &
                 ,rain, rainncv                                    & !(private - todas)
                 ,snow, snowncv                                    & !(private - todas)
                 ,sr                                               & !(private)
                 ,refl_10cm, diagflag, do_radar_ref                & !(private - somente refl_10cm)
                 ,has_reqc, has_reqi, has_reqs                     &  ! for radiation
                 ,re_cloud, re_ice, re_snow                        &  ! for radiation (private - todas)      
                 ,ids,ide, jds,jde, kds,kde                        &
                 ,ims,ime, jms,jme, kms,kme                        &
                 ,its,ite, jts,jte, kts,kte                        &
                                                                   )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
   !-srf
   ! USE module_model_constants, only : RE_QC_BG, RE_QI_BG, RE_QS_BG
   REAL    , PARAMETER :: RE_QC_BG     = 2.51E-6     ! effective radius of cloud for background (m)
   REAL    , PARAMETER :: RE_QI_BG     = 5.01E-6     ! effective radius of ice for background (m)
   REAL    , PARAMETER :: RE_QS_BG     = 10.01E-6    ! effective radius of snow for background (m)
   !-srf
   
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                          &
                                                             th,  &
                                                              q,  &
                                                              qc, &
                                                              qi, &
                                                              qr, &
                                                              qs
  ! REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        ! INTENT(INOUT) ::                                          &
                                                             ! th,  &
                                                              ! q,  &
                                                              ! qc, &
                                                              ! qi, &
                                                              ! qr, &
                                                              ! qs
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                             pii, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                              rd, &
                                                              rv, &
                                                             t0c, &
                                                            den0, &
                                                             cpd, &
                                                             cpv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
! for radiation connecting
  INTEGER, INTENT(IN)::                                           &
                                                        has_reqc, &
                                                        has_reqi, &
                                                        has_reqs
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme),                     &
        INTENT(INOUT)::                                           &
                                                        re_cloud, &
                                                          re_ice, &
                                                         re_snow
!+---+-----------------------------------------------------------------+
  REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::     &  ! GT
                                                       refl_10cm
!+---+-----------------------------------------------------------------+

  REAL, DIMENSION( ims:ime , jms:jme ), OPTIONAL,                 &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
! LOCAL VAR
!#ifdef XEON_OPTIMIZED_WSM5
!#  include "mic-wsm5-3-5-locvar.h"
!#else
  REAL, DIMENSION( its:ite , kts:kte ) ::   t
  REAL, DIMENSION( its:ite , kts:kte, 2 ) ::   qci, qrs
  CHARACTER*256 :: emess
  INTEGER ::               i,j,k
!#endif

!+---+-----------------------------------------------------------------+
      REAL, DIMENSION(kts:kte):: qv1d, t1d, p1d, qr1d, qs1d, dBZ
      LOGICAL, OPTIONAL, INTENT(IN) :: diagflag
      INTEGER, OPTIONAL, INTENT(IN) :: do_radar_ref
!+---+-----------------------------------------------------------------+

! to calculate effective radius for radiation
  REAL, DIMENSION( kts:kte ) :: den1d
  REAL, DIMENSION( kts:kte ) :: qc1d
  REAL, DIMENSION( kts:kte ) :: qi1d
  REAL, DIMENSION( kts:kte ) :: re_qc, re_qi, re_qs

!+---+-----------------------------------------------------------------+
!$acc routine seq
! return
!TEMPO COM TUDO COMENTADO DAQUI PRA BAIXO 0.5547960

! #ifndef XEON_OPTIMIZED_WSM5
! print *,"ifndef XEON_OPTIMIZED_WSM5 OK"
      !! DO j=jts,jte
	  j=jts ! BLOCK 1 - 0.7280159 !! DESTA FORMA GPU JA ESTA PIOR QUE CPU :-(
         ! DO k=kts,kte
         ! DO i=its,ite
            ! t(i,k)=th(i,k,j)*pii(i,k,j)
            ! qci(i,k,1) = qc(i,k,j)
            ! qci(i,k,2) = qi(i,k,j)
            ! qrs(i,k,1) = qr(i,k,j)
            ! qrs(i,k,2) = qs(i,k,j)
         ! ENDDO
         ! ENDDO
          ! Sending array starting locations of optional variables may cause
          ! troubles, so we explicitly change the call.
		  ! th(1,kts:kte,j), pii(1,kts:kte,j)
         CALL wsm52D(&!t, q(ims,kms,j), qci, qrs &
					th(1,kts:kte,j), pii(1,kts:kte,j), q(ims,kms,j),& !(private - todas)
					qc(1,kts:kte,j),qi(1,kts:kte,j),& !(private - todas)
					qr(1,kts:kte,j), qs(1,kts:kte,j) & !(private - todas)
                    ,den(ims,kms,j)                                & !(private - todas)
                    ,p(ims,kms,j), delz(ims,kms,j),                & !(private - todas)
                    delt,g, cpd, cpv, rd, rv, t0c                  & !(private - somente delt)
                    ,ep1, ep2, qmin                                &
                    ,XLS, XLV0, XLF0, den0, denr                   &
                    ,cliq,cice,psat                                &
                    ,j                                             &
                    ,rain(ims,j),rainncv(ims,j)                    & !(private - todas)
                    ,sr(ims,j)                                     & !(private - todas)
                    ,ids,ide, jds,jde, kds,kde                     &
                    ,ims,ime, jms,jme, kms,kme                     &
                    ,its,ite, jts,jte, kts,kte                     &
                    ,snow,snowncv                                  & !(private - todas)
                                                                  )
														   
         ! DO K=kts,kte
         ! DO I=its,ite
            ! th(i,k,j)=t(i,k)/pii(i,k,j)
            ! qc(i,k,j) = qci(i,k,1)
            ! qi(i,k,j) = qci(i,k,2)
            ! qr(i,k,j) = qrs(i,k,1)
            ! qs(i,k,j) = qrs(i,k,2)
         ! ENDDO
         ! ENDDO
		! BLOCK 1 - 0.7280159
! +---+-----------------------------------------------------------------+

! +---+-----------------------------------------------------------------+
  END SUBROUTINE wsm5

!#ifndef XEON_OPTIMIZED_WSM5
!===================================================================
! t(i,k)=th(i,k,j)*pii(i,k,j)
! th(i,k,j)=t(i,k)/pii(i,k,j)
!th, pii
  SUBROUTINE wsm52D(& !t, q, qci, qrs, &
					th, pii, q,                                   & 
                   qc, qi, & 
				   qr, qs, &
				   den, p, delz,                        &
                   delt,g, cpd, cpv, rd, rv, t0c                 &
                   ,ep1, ep2, qmin                                &
                   ,XLS, XLV0, XLF0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,rain,rainncv                                  &
                   ,sr                                            &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte                     &
                   ,snow,snowncv                                  &
                                                                  )
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
!  This code is a 5-class mixed ice microphyiscs scheme (WSM5) of the 
!  Single-Moment MicroPhyiscs (WSMMP). The WSMMP assumes that ice nuclei
!  number concentration is a function of temperature, and seperate assumption
!  is developed, in which ice crystal number concentration is a function
!  of ice amount. A theoretical background of the ice-microphysics and related
!  processes in the WSMMPs are described in Hong et al. (2004).
!  Production terms in the WSM6 scheme are described in Hong and Lim (2006).
!  All units are in m.k.s. and source/sink terms in kgkg-1s-1.
!
!  WSM5 cloud scheme
!
!  Coded by Song-You Hong (Yonsei Univ.)
!             Jimy Dudhia (NCAR) and Shu-Hua Chen (UC Davis)
!             Summer 2002
!
!  Implemented by Song-You Hong (Yonsei Univ.) and Jimy Dudhia (NCAR)
!             Summer 2003
!
!  further modifications :
!        semi-lagrangian sedimentation (JH,2010),hong, aug 2009
!        ==> higher accuracy and efficient at lower resolutions
!        reflectivity computation from greg thompson, lim, jun 2011
!        ==> only diagnostic, but with removal of too large drops
!        effective radius of hydrometeors, bae from kiaps, jan 2015
!        ==> consistency in solar insolation of rrtmg radiation
!        bug fix in melting terms, bae from kiaps, nov 2015
!        ==> density of air is divided, which has not been
!
!  Reference) Hong, Dudhia, Chen (HDC, 2004) Mon. Wea. Rev.
!             Rutledge, Hobbs (RH83, 1983) J. Atmos. Sci.
!             Hong and Lim (HL, 2006) J. Korean Meteor. Soc.
!
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte,  &
                                      lat
  REAL, DIMENSION( its:ite , kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                              th, pii
		! REAL, DIMENSION( its:ite , kts:kte ),                           &
        ! INTENT(INOUT) ::                                          &
                                                              ! t
		! REAL, DIMENSION( its:ite , kts:kte, 2 ),                        &
        ! INTENT(INOUT) ::                                          &
                                                             ! qci, &
                                                              ! qrs
  REAL, DIMENSION( its:ite , kts:kte),                        &
        INTENT(INOUT) ::                                          &
															  qc, qi, &
															  qr, qs
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL, DIMENSION( ims:ime , kms:kme ),                           &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                               p, &
                                                            delz
  REAL, INTENT(IN   ) ::                                    delt, &
                                                               g, &
                                                             cpd, &
                                                             cpv, &
                                                             t0c, &
                                                            den0, &
                                                              rd, &
                                                              rv, &
                                                             ep1, &
                                                             ep2, &
                                                            qmin, &
                                                             XLS, &
                                                            XLV0, &
                                                            XLF0, &
                                                            cliq, &
                                                            cice, &
                                                            psat, &
                                                            denr
  REAL, DIMENSION( ims:ime ),                                     &
        INTENT(INOUT) ::                                    rain, &
                                                         rainncv, &
                                                              sr
  REAL, DIMENSION( ims:ime, jms:jme ),     OPTIONAL,              &
        INTENT(INOUT) ::                                    snow, &
                                                         snowncv
! LOCAL VAR
  REAL, DIMENSION( its:ite , kts:kte , 2) ::                      &
                                                              rh, &
                                                              qss, &
                                                          rslope, &
                                                         rslope2, &
                                                         rslope3, &
                                                         rslopeb, &
                                                         qrs_tmp, &
                                                            falk, &
                                                            fall, &
                                                           work1
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                           falkc, &
                                                           fallc, &
                                                              xl, &
                                                             cpm, &
                                                          denfac, &
                                                             xni, &
                                                         denqrs1, &
                                                         denqrs2, &
                                                          denqci, &
                                                          n0sfac, &
                                                           work2, &
                                                           workr, &
                                                           works, &
                                                          work1c, &
                                                          work2c
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                         den_tmp, &
                                                        delz_tmp
  REAL, DIMENSION( its:ite ) ::                                   &
                                                         delqrs1, &
                                                         delqrs2, &
                                                           delqi
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                           pigen, &
                                                           pidep, &
                                                           psdep, &
                                                           praut, &
                                                           psaut, &
                                                           prevp, &
                                                           psevp, &
                                                           pracw, &
                                                           psacw, &
                                                           psaci, &
                                                           pcond, &
                                                           psmlt
  INTEGER, DIMENSION( its:ite ) ::                                &
                                                           mstep, &
                                                           numdt
  REAL, DIMENSION(its:ite) ::                             tstepsnow
  REAL, DIMENSION(its:ite) ::                             rmstep
  REAL dtcldden, rdelz, rdtcld
  LOGICAL, DIMENSION( its:ite ) ::                        flgcld
!#define WSM_NO_CONDITIONAL_IN_VECTOR
!#ifdef WSM_NO_CONDITIONAL_IN_VECTOR
!  REAL, DIMENSION(its:ite) :: xal, xbl
!#endif
  REAL  ::                                                        &
            cpmcal, xlcal, diffus,                                &
            viscos, xka, venfac, conden, diffac,                  &
            x, y, z, a, b, c, d, e,                               &
            qdt, holdrr, holdrs, supcol, supcolt, pvt,            &
            coeres, supsat, dtcld, xmi, eacrs, satdt,             &
            vt2i,vt2s,acrfac,                                     &
            qimax, diameter, xni0, roqi0,                         &
            fallsum, fallsum_qsi, xlwork2, factor, source,        &
            value, xlf, pfrzdtc, pfrzdtr, supice,  holdc, holdci
! variables for optimization
  REAL, DIMENSION( its:ite )           ::                    tvec1
  REAL ::                                                    temp 
  INTEGER :: i, j, k, mstepmax,                                   &
            iprt, latd, lond, loop, loops, ifsat, n, idim, kdim
! Temporaries used for inlining fpvs function
  REAL  :: dldti, xb, xai, tr, xbi, xa, hvap, cvap, hsub, dldt, ttp
  REAL  :: logtr
!$acc routine seq
return

  END SUBROUTINE wsm52d

! ...................................................................
!--------------------------------------------------------------------------
      REAL FUNCTION fpvs(t,ice,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c)
!--------------------------------------------------------------------------
      IMPLICIT NONE
!--------------------------------------------------------------------------
      REAL t,rd,rv,cvap,cliq,cice,hvap,hsub,psat,t0c,dldt,xa,xb,dldti,         &
           xai,xbi,ttp,tr
      INTEGER ice
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)
      tr=ttp/t
      if(t.lt.ttp.and.ice.eq.1) then
        fpvs=psat*exp(log(tr)*(xai))*exp(xbi*(1.-tr))
      else
        fpvs=psat*exp(log(tr)*(xa))*exp(xb*(1.-tr))
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END FUNCTION fpvs
!------------------------------------------------------------------------------
      subroutine slope_wsm5(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  REAL, DIMENSION( its:ite , kts:kte,2) ::                                     &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
  REAL, DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  REAL       ::  lamdar, lamdas,  x, y, z, supcol
  integer :: i, j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
!
      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k,1).le.qcrmin)then
            rslope(i,k,1) = rslopermax
            rslopeb(i,k,1) = rsloperbmax
            rslope2(i,k,1) = rsloper2max
            rslope3(i,k,1) = rsloper3max
          else
            rslope(i,k,1) = 1./lamdar(qrs(i,k,1),den(i,k))
            rslopeb(i,k,1) = exp(log(rslope(i,k,1))*(bvtr))
            rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
            rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
          endif
          if(qrs(i,k,2).le.qcrmin)then
            rslope(i,k,2) = rslopesmax
            rslopeb(i,k,2) = rslopesbmax
            rslope2(i,k,2) = rslopes2max
            rslope3(i,k,2) = rslopes3max
          else
            rslope(i,k,2) = 1./lamdas(qrs(i,k,2),den(i,k),n0sfac(i,k))
            rslopeb(i,k,2) = exp(log(rslope(i,k,2))*(bvts))
            rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
            rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
          endif
          vt(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)
          if(qrs(i,k,1).le.0.0) vt(i,k,1) = 0.0
          if(qrs(i,k,2).le.0.0) vt(i,k,2) = 0.0
        enddo
      enddo
  END subroutine slope_wsm5
!-----------------------------------------------------------------------------
      subroutine slope_rain(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  REAL, DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  REAL       ::  lamdar, x, y, z, supcol
  integer :: i, j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdar(x,y)=   sqrt(sqrt(pidn0r/(x*y)))      ! (pidn0r/(x*y))**.25
!
      do k = kts, kte
        do i = its, ite
          if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopermax
            rslopeb(i,k) = rsloperbmax
            rslope2(i,k) = rsloper2max
            rslope3(i,k) = rsloper3max
          else
            rslope(i,k) = 1./lamdar(qrs(i,k),den(i,k))
            rslopeb(i,k) = rslope(i,k)**bvtr
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
          endif
          vt(i,k) = pvtr*rslopeb(i,k)*denfac(i,k)
          if(qrs(i,k).le.0.0) vt(i,k) = 0.0
        enddo
      enddo
  END subroutine slope_rain
!------------------------------------------------------------------------------
      subroutine slope_snow(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,   &
                            vt,its,ite,kts,kte)
  IMPLICIT NONE
  INTEGER       ::               its,ite, jts,jte, kts,kte
  REAL, DIMENSION( its:ite , kts:kte) ::                                       &
                                                                          qrs, &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  REAL, PARAMETER  :: t0c = 273.15
  REAL, DIMENSION( its:ite , kts:kte ) ::                                      &
                                                                       n0sfac
  REAL       ::  lamdas, x, y, z, supcol
  integer :: i, j, k
!----------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
      lamdas(x,y,z)= sqrt(sqrt(pidn0s*z/(x*y)))    ! (pidn0s*z/(x*y))**.25
!
      do k = kts, kte
        do i = its, ite
          supcol = t0c-t(i,k)
!---------------------------------------------------------------
! n0s: Intercept parameter for snow [m-4] [HDC 6]
!---------------------------------------------------------------
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          if(qrs(i,k).le.qcrmin)then
            rslope(i,k) = rslopesmax
            rslopeb(i,k) = rslopesbmax
            rslope2(i,k) = rslopes2max
            rslope3(i,k) = rslopes3max
          else
            rslope(i,k) = 1./lamdas(qrs(i,k),den(i,k),n0sfac(i,k))
            rslopeb(i,k) = rslope(i,k)**bvts
            rslope2(i,k) = rslope(i,k)*rslope(i,k)
            rslope3(i,k) = rslope2(i,k)*rslope(i,k)
          endif
          vt(i,k) = pvts*rslopeb(i,k)*denfac(i,k)
          if(qrs(i,k).le.0.0) vt(i,k) = 0.0
        enddo
      enddo
  END subroutine slope_snow
!-------------------------------------------------------------------
      SUBROUTINE nislfv_rain_plm(im,km,denl,denfacl,tkl,dzl,wwl,rql,precip,dt,id,iter)
!-------------------------------------------------------------------
!
! for non-iteration semi-Lagrangain forward advection for cloud
! with mass conservation and positive definite advection
! 2nd order interpolation with monotonic piecewise linear method
! this routine is under assumption of decfl < 1 for semi_Lagrangian
!
! dzl    depth of model layer in meter
! wwl    terminal velocity at model layer m/s
! rql    cloud density*mixing ration
! precip precipitation
! dt     time step
! id     kind of precip: 0 test case; 1 raindrop  2: snow
! iter   how many time to guess mean terminal velocity: 0 pure forward.
!        0 : use departure wind for advection
!        1 : use mean wind for advection
!        > 1 : use mean wind after iter-1 iterations
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         implemented by song-you hong
!
      implicit none
      integer  im,km,id
      real  dt
      real  dzl(im,km),wwl(im,km),rql(im,km),precip(im)
      real  denl(im,km),denfacl(im,km),tkl(im,km)
!
      integer  i,k,n,m,kk,kb,kt,iter
      real  tl,tl2,qql,dql,qqd
      real  th,th2,qqh,dqh
      real  zsum,qsum,dim,dip,c1,con1,fa1,fa2
      real  allold, allnew, zz, dzamin, cflmax, decfl
      real  dz(km), ww(km), qq(km), wd(km), wa(km), was(km)
      real  den(km), denfac(km), tk(km)
      real  wi(km+1), zi(km+1), za(km+1)
      real  qn(km), qr(km),tmp(km),tmp1(km),tmp2(km),tmp3(km)
      real  dza(km+1), qa(km+1), qmi(km+1), qpi(km+1)
!
      precip(:) = 0.0
!
      i_loop : do i=1,im
! -----------------------------------
      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)
! skip for no precipitation for all layers
      allold = 0.0
      do k=1,km
        allold = allold + qq(k)
      enddo
      if(allold.le.0.0) then
        cycle i_loop
      endif
!
! compute interface values
      zi(1)=0.0
      do k=1,km
        zi(k+1) = zi(k)+dz(k)
      enddo
!
! save departure wind
      wd(:) = ww(:)
      n=1
 100  continue
! plm is 2nd order, we can use 2nd order wi or 3rd order wi
! 2nd order interpolation to get wi
      wi(1) = ww(1)
      wi(km+1) = ww(km)
      do k=2,km
        wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
      enddo
! 3rd order interpolation to get wi
      fa1 = 9./16.
      fa2 = 1./16.
      wi(1) = ww(1)
      wi(2) = 0.5*(ww(2)+ww(1))
      do k=3,km-1
        wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
      enddo
      wi(km) = 0.5*(ww(km)+ww(km-1))
      wi(km+1) = ww(km)
!
! terminate of top of raingroup
      do k=2,km
        if( ww(k).eq.0.0 ) wi(k)=ww(k-1)
      enddo
!
! diffusivity of wi
      con1 = 0.05
      do k=km,1,-1
        decfl = (wi(k+1)-wi(k))*dt/dz(k)
        if( decfl .gt. con1 ) then
          wi(k) = wi(k+1) - con1*dz(k)/dt
        endif
      enddo
! compute arrival point
      do k=1,km+1
        za(k) = zi(k) - wi(k)*dt
      enddo
!
      do k=1,km
        dza(k) = za(k+1)-za(k)
      enddo
      dza(km+1) = zi(km+1) - za(km+1)
!
! computer deformation at arrival point
      do k=1,km
        qa(k) = qq(k)*dz(k)/dza(k)
        qr(k) = qa(k)/den(k)
      enddo
      qa(km+1) = 0.0
!     call maxmin(km,1,qa,' arrival points ')
!
! compute arrival terminal velocity, and estimate mean terminal velocity
! then back to use mean terminal velocity
      if( n.le.iter ) then
        if (id.eq.1) then
          call slope_rain(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,1,1,1,km)
        else
          call slope_snow(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa,1,1,1,km)
        endif 
        if( n.ge.2 ) wa(1:km)=0.5*(wa(1:km)+was(1:km))
        do k=1,km
!#ifdef DEBUG
!        print*,' slope_wsm3 ',qr(k)*1000.,den(k),denfac(k),tk(k),tmp(k),tmp1(k),tmp2(k),ww(k),wa(k)
!#endif
! mean wind is average of departure and new arrival winds
          ww(k) = 0.5* ( wd(k)+wa(k) )
        enddo
        was(:) = wa(:)
        n=n+1
        go to 100
      endif
!
! estimate values at arrival cell interface with monotone
      do k=2,km
        dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
        dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
        if( dip*dim.le.0.0 ) then
          qmi(k)=qa(k)
          qpi(k)=qa(k)
        else
          qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
          qmi(k)=2.0*qa(k)-qpi(k)
          if( qpi(k).lt.0.0 .or. qmi(k).lt.0.0 ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          endif
        endif
      enddo
      qpi(1)=qa(1)
      qmi(1)=qa(1)
      qmi(km+1)=qa(km+1)
      qpi(km+1)=qa(km+1)
!
! interpolation to regular point
      qn = 0.0
      kb=1
      kt=1
      intp : do k=1,km
             kb=max(kb-1,1)
             kt=max(kt-1,1)
! find kb and kt
             if( zi(k).ge.za(km+1) ) then
               exit intp
             else
               find_kb : do kk=kb,km
                         if( zi(k).le.za(kk+1) ) then
                           kb = kk
                           exit find_kb
                         else
                           cycle find_kb
                         endif
               enddo find_kb
               find_kt : do kk=kt,km
                         if( zi(k+1).le.za(kk) ) then
                           kt = kk
                           exit find_kt
                         else
                           cycle find_kt
                         endif
               enddo find_kt
               kt = kt - 1
! compute q with piecewise constant method
               if( kt.eq.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 th=(zi(k+1)-za(kb))/dza(kb)
                 tl2=tl*tl
                 th2=th*th
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qqh=qqd*th2+qmi(kb)*th
                 qql=qqd*tl2+qmi(kb)*tl
                 qn(k) = (qqh-qql)/(th-tl)
               else if( kt.gt.kb ) then
                 tl=(zi(k)-za(kb))/dza(kb)
                 tl2=tl*tl
                 qqd=0.5*(qpi(kb)-qmi(kb))
                 qql=qqd*tl2+qmi(kb)*tl
                 dql = qa(kb)-qql
                 zsum  = (1.-tl)*dza(kb)
                 qsum  = dql*dza(kb)
                 if( kt-kb.gt.1 ) then
                 do m=kb+1,kt-1
                   zsum = zsum + dza(m)
                   qsum = qsum + qa(m) * dza(m)
                 enddo
                 endif
                 th=(zi(k+1)-za(kt))/dza(kt)
                 th2=th*th
                 qqd=0.5*(qpi(kt)-qmi(kt))
                 dqh=qqd*th2+qmi(kt)*th
                 zsum  = zsum + th*dza(kt)
                 qsum  = qsum + dqh*dza(kt)
                 qn(k) = qsum/zsum
               endif
               cycle intp
             endif
!
       enddo intp
!
! rain out
      sum_precip: do k=1,km
                    if( za(k).lt.0.0 .and. za(k+1).lt.0.0 ) then
                      precip(i) = precip(i) + qa(k)*dza(k)
                      cycle sum_precip
                    else if ( za(k).lt.0.0 .and. za(k+1).ge.0.0 ) then
                      precip(i) = precip(i) + qa(k)*(0.0-za(k))
                      exit sum_precip
                    endif
                    exit sum_precip
      enddo sum_precip
!
! replace the new values
      rql(i,:) = qn(:)
!
! ----------------------------------
      enddo i_loop
!
  END SUBROUTINE nislfv_rain_plm

!# else
!#  include "mic-wsm5-3-5-code.h"
!# endif
! end of #ifndef XEON_OPTIMIZED_WSM5

! remainder of routines are common to original and MIC version
      subroutine refl10cm_wsm5_x (qv1d, qr1d, qs1d,                       &
                       th, pii, p1d, dBZ, kts, kte, ii, jj)
      !newCode begin
      !$acc routine(rayleigh_soak_wetgraupel)
      !newCode end

      IMPLICIT NONE
!$acc routine seq
!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                      qv1d, qr1d, qs1d, th, pii, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ

!..Local variables
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
      REAL, DIMENSION(kts:kte):: rr, rs
      REAL:: temp_C

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilams
      DOUBLE PRECISION, DIMENSION(kts:kte):: N0_r, N0_s
      DOUBLE PRECISION:: lamr, lams
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow
      DOUBLE PRECISION:: fmelt_s

      INTEGER:: i, k, k_0, kbot, n
      LOGICAL:: melti

      DOUBLE PRECISION:: cback, x, eta, f_d
      REAL, PARAMETER:: R=287.

!+---+

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = th(k)*pii(k) !t1d(k) 04/12/2023
         temp_C = min(-0.001, temp(K)-273.15)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))

         if (qr1d(k) .gt. 1.E-9) then
            rr(k) = qr1d(k)*rho(k)
            N0_r(k) = n0r
            lamr = (xam_r*xcrg(3)*N0_r(k)/rr(k))**(1./xcre(1))
            ilamr(k) = 1./lamr
            L_qr(k) = .true.
         else
            rr(k) = 1.E-12
            L_qr(k) = .false.
         endif

         if (qs1d(k) .gt. 1.E-9) then
            rs(k) = qs1d(k)*rho(k)
            N0_s(k) = min(n0smax, n0s*exp(-alpha*temp_C))
            lams = (xam_s*xcsg(3)*N0_s(k)/rs(k))**(1./xcse(1))
            ilams(k) = 1./lams
            L_qs(k) = .true.
         else
            rs(k) = 1.E-12
            L_qs(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt.273.15) .and. L_qr(k) .and. L_qs(k+1) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue

!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
      enddo


!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.kts+1) then
       do k = k_0-1, kts, -1

!..Reflectivity contributed by melting snow
          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           lams = 1./ilams(k)
           do n = 1, nrbins
              x = xam_s * xxDs(n)**xbm_s
              ! call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                    ! fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    ! CBACK, mixingrulestring_s, matrixstring_s,          &
                    ! inclusionstring_s, hoststring_s,                    &
                    ! hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif
       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k))*1.d18)
      enddo

      end subroutine refl10cm_wsm5_x
!+---+-----------------------------------------------------------------+
      subroutine refl10cm_wsm5 (qv1d, qr1d, qs1d,                       &
                       t1d, p1d, dBZ, kts, kte, ii, jj)
      !newCode begin
      !$acc routine(rayleigh_soak_wetgraupel)
      !newCode end

      IMPLICIT NONE

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                      qv1d, qr1d, qs1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ

!..Local variables
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
      REAL, DIMENSION(kts:kte):: rr, rs
      REAL:: temp_C

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilams
      DOUBLE PRECISION, DIMENSION(kts:kte):: N0_r, N0_s
      DOUBLE PRECISION:: lamr, lams
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow
      DOUBLE PRECISION:: fmelt_s

      INTEGER:: i, k, k_0, kbot, n
      LOGICAL:: melti

      DOUBLE PRECISION:: cback, x, eta, f_d
      REAL, PARAMETER:: R=287.

!+---+

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         temp_C = min(-0.001, temp(K)-273.15)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))

         if (qr1d(k) .gt. 1.E-9) then
            rr(k) = qr1d(k)*rho(k)
            N0_r(k) = n0r
            lamr = (xam_r*xcrg(3)*N0_r(k)/rr(k))**(1./xcre(1))
            ilamr(k) = 1./lamr
            L_qr(k) = .true.
         else
            rr(k) = 1.E-12
            L_qr(k) = .false.
         endif

         if (qs1d(k) .gt. 1.E-9) then
            rs(k) = qs1d(k)*rho(k)
            N0_s(k) = min(n0smax, n0s*exp(-alpha*temp_C))
            lams = (xam_s*xcsg(3)*N0_s(k)/rs(k))**(1./xcse(1))
            ilams(k) = 1./lams
            L_qs(k) = .true.
         else
            rs(k) = 1.E-12
            L_qs(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt.273.15) .and. L_qr(k) .and. L_qs(k+1) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue

!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
      enddo


!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.kts+1) then
       do k = k_0-1, kts, -1

!..Reflectivity contributed by melting snow
          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           lams = 1./ilams(k)
           do n = 1, nrbins
              x = xam_s * xxDs(n)**xbm_s
              ! call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                    ! fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    ! CBACK, mixingrulestring_s, matrixstring_s,          &
                    ! inclusionstring_s, hoststring_s,                    &
                    ! hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif
       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k))*1.d18)
      enddo

      end subroutine refl10cm_wsm5
!+---+-----------------------------------------------------------------+
!-------------------------------------------------------------------
  SUBROUTINE wsm5init(den0,denr,dens,cl,cpv,allowed_to_read)
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!.... constants which may not be tunable
   REAL, INTENT(IN) :: den0,denr,dens,cl,cpv
   LOGICAL, INTENT(IN) :: allowed_to_read
!
   pi = 4.*atan(1.)
   xlv1 = cl-cpv
!
   qc0  = 4./3.*pi*denr*r0**3*xncr/den0  ! 0.419e-3 -- .61e-3
   qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu*den0**(4./3.) ! 7.03
   pidnc = pi*denr/6.        ! syb
!
   bvtr1 = 1.+bvtr
   bvtr2 = 2.5+.5*bvtr
   bvtr3 = 3.+bvtr
   bvtr4 = 4.+bvtr
   g1pbr = rgmma(bvtr1)
   g3pbr = rgmma(bvtr3)
   g4pbr = rgmma(bvtr4)            ! 17.837825
   g5pbro2 = rgmma(bvtr2)          ! 1.8273
   pvtr = avtr*g4pbr/6.
   eacrr = 1.0
   pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
   precr1 = 2.*pi*n0r*.78
   precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2
   xmmax = (dimax/dicon)**2
   roqimax = 2.08e22*dimax**8
!
   bvts1 = 1.+bvts
   bvts2 = 2.5+.5*bvts
   bvts3 = 3.+bvts
   bvts4 = 4.+bvts
   g1pbs = rgmma(bvts1)    !.8875
   g3pbs = rgmma(bvts3)
   g4pbs = rgmma(bvts4)    ! 12.0786
   g5pbso2 = rgmma(bvts2)
   pvts = avts*g4pbs/6.
   pacrs = pi*n0s*avts*g3pbs*.25
   precs1 = 4.*n0s*.65
   precs2 = 4.*n0s*.44*avts**.5*g5pbso2
   pidn0r =  pi*denr*n0r
   pidn0s =  pi*dens*n0s
   pacrc = pi*n0s*avts*g3pbs*.25*eacrc
!
   rslopermax = 1./lamdarmax
   rslopesmax = 1./lamdasmax
   rsloperbmax = rslopermax ** bvtr
   rslopesbmax = rslopesmax ** bvts
   rsloper2max = rslopermax * rslopermax
   rslopes2max = rslopesmax * rslopesmax
   rsloper3max = rsloper2max * rslopermax
   rslopes3max = rslopes2max * rslopesmax
!
!+---+-----------------------------------------------------------------+
!..Set these variables needed for computing radar reflectivity.  These
!.. get used within radar_init to create other variables used in the
!.. radar module.
   xam_r = PI*denr/6.
   xbm_r = 3.
   xmu_r = 0.
   xam_s = PI*dens/6.
   xbm_s = 3.
   xmu_s = 0.
   xam_g = PI*dens/6.      !  These 3 variables for graupel are set but unused.
   xbm_g = 3.
   xmu_g = 0.

   call radar_init
!+---+-----------------------------------------------------------------+

  END SUBROUTINE wsm5init
! ...................................................................
      REAL FUNCTION rgmma(x)
!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!     rgmma function:  use infinite product form
      REAL :: euler
      PARAMETER (euler=0.577215664901532)
      REAL :: x, y
      INTEGER :: i
      if(x.eq.1.)then
        rgmma=0.
          else
        rgmma=x*exp(euler*x)
        do i=1,10000
          y=float(i)
          rgmma=rgmma*(1.000+x/y)*exp(-x/y)
        enddo
        rgmma=1./rgmma
      endif
      END FUNCTION rgmma

     subroutine effectRad_wsm5_x (th, pii, qc, qi, qs, rho, qmin, t0c,        &
                                re_qc, re_qi, re_qs, kts, kte, ii, jj)

!-----------------------------------------------------------------------
!  Compute radiation effective radii of cloud water, ice, and snow for 
!  single-moment microphysics.
!  These are entirely consistent with microphysics assumptions, not
!  constant or otherwise ad hoc as is internal to most radiation
!  schemes.  
!  Coded and implemented by Soo ya Bae, KIAPS, January 2015.
!-----------------------------------------------------------------------

      implicit none
!$acc routine seq
!..Sub arguments
      integer, intent(in) :: kts, kte, ii, jj
      real, intent(in) :: qmin
      real, intent(in) :: t0c
      ! real, dimension( kts:kte ), intent(in)::  t !30/11/2023
	  real, dimension( kts:kte ), intent(in)::  th, pii
      real, dimension( kts:kte ), intent(in)::  qc
      real, dimension( kts:kte ), intent(in)::  qi
      real, dimension( kts:kte ), intent(in)::  qs
      real, dimension( kts:kte ), intent(in)::  rho
      real, dimension( kts:kte ), intent(inout):: re_qc
      real, dimension( kts:kte ), intent(inout):: re_qi
      real, dimension( kts:kte ), intent(inout):: re_qs
!..Local variables
      integer:: i,k
      integer :: inu_c
      real, dimension( kts:kte ):: ni
      real, dimension( kts:kte ):: rqc
      real, dimension( kts:kte ):: rqi
      real, dimension( kts:kte ):: rni
      real, dimension( kts:kte ):: rqs
      real :: temp
      real :: lamdac
      real :: supcol, n0sfac, lamdas
      real :: diai      ! diameter of ice in m
      logical :: has_qc, has_qi, has_qs
!..Minimum microphys values
      real, parameter :: R1 = 1.E-12
      real, parameter :: R2 = 1.E-6
!..Mass power law relations:  mass = am*D**bm
      real, parameter :: bm_r = 3.0
      real, parameter :: obmr = 1.0/bm_r
      real, parameter :: nc0  = 3.E8
!-----------------------------------------------------------------------
      has_qc = .false.
      has_qi = .false.
      has_qs = .false.

      do k = kts, kte
        ! for cloud
        rqc(k) = max(R1, qc(k)*rho(k))
        if (rqc(k).gt.R1) has_qc = .true.
        ! for ice
        rqi(k) = max(R1, qi(k)*rho(k))
        temp = (rho(k)*max(qi(k),qmin))
        temp = sqrt(sqrt(temp*temp*temp))
        ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
        rni(k)= max(R2, ni(k)*rho(k))
        if (rqi(k).gt.R1 .and. rni(k).gt.R2) has_qi = .true.
        ! for snow
        rqs(k) = max(R1, qs(k)*rho(k))
        if (rqs(k).gt.R1) has_qs = .true.
      enddo

      if (has_qc) then
        do k=kts,kte
          if (rqc(k).le.R1) CYCLE
          lamdac   = (pidnc*nc0/rqc(k))**obmr
          re_qc(k) =  max(2.51E-6,min(1.5*(1.0/lamdac),50.E-6))
        enddo
      endif

     if (has_qi) then
        do k=kts,kte
          if (rqi(k).le.R1 .or. rni(k).le.R2) CYCLE
          diai = 11.9*sqrt(rqi(k)/ni(k))
          re_qi(k) = max(10.01E-6,min(0.75*0.163*diai,125.E-6))
        enddo
      endif

      if (has_qs) then
        do k=kts,kte
          if (rqs(k).le.R1) CYCLE
          ! supcol = t0c-t(k) th(i,k,j)*pii(i,k,j) !30/11/2023
		  supcol = t0c-(th(k)*pii(k))
          n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs(k)))
          re_qs(k) = max(25.E-6,min(0.5*(1./lamdas), 999.E-6))
        enddo
      endif

      end subroutine effectRad_wsm5_x
!-----------------------------------------------------------------------
     subroutine effectRad_wsm5 (t, qc, qi, qs, rho, qmin, t0c,        &
                                re_qc, re_qi, re_qs, kts, kte, ii, jj)

!-----------------------------------------------------------------------
!  Compute radiation effective radii of cloud water, ice, and snow for 
!  single-moment microphysics.
!  These are entirely consistent with microphysics assumptions, not
!  constant or otherwise ad hoc as is internal to most radiation
!  schemes.  
!  Coded and implemented by Soo ya Bae, KIAPS, January 2015.
!-----------------------------------------------------------------------

      implicit none

!..Sub arguments
      integer, intent(in) :: kts, kte, ii, jj
      real, intent(in) :: qmin
      real, intent(in) :: t0c
      real, dimension( kts:kte ), intent(in)::  t
      real, dimension( kts:kte ), intent(in)::  qc
      real, dimension( kts:kte ), intent(in)::  qi
      real, dimension( kts:kte ), intent(in)::  qs
      real, dimension( kts:kte ), intent(in)::  rho
      real, dimension( kts:kte ), intent(inout):: re_qc
      real, dimension( kts:kte ), intent(inout):: re_qi
      real, dimension( kts:kte ), intent(inout):: re_qs
!..Local variables
      integer:: i,k
      integer :: inu_c
      real, dimension( kts:kte ):: ni
      real, dimension( kts:kte ):: rqc
      real, dimension( kts:kte ):: rqi
      real, dimension( kts:kte ):: rni
      real, dimension( kts:kte ):: rqs
      real :: temp
      real :: lamdac
      real :: supcol, n0sfac, lamdas
      real :: diai      ! diameter of ice in m
      logical :: has_qc, has_qi, has_qs
!..Minimum microphys values
      real, parameter :: R1 = 1.E-12
      real, parameter :: R2 = 1.E-6
!..Mass power law relations:  mass = am*D**bm
      real, parameter :: bm_r = 3.0
      real, parameter :: obmr = 1.0/bm_r
      real, parameter :: nc0  = 3.E8
!-----------------------------------------------------------------------
      has_qc = .false.
      has_qi = .false.
      has_qs = .false.

      do k = kts, kte
        ! for cloud
        rqc(k) = max(R1, qc(k)*rho(k))
        if (rqc(k).gt.R1) has_qc = .true.
        ! for ice
        rqi(k) = max(R1, qi(k)*rho(k))
        temp = (rho(k)*max(qi(k),qmin))
        temp = sqrt(sqrt(temp*temp*temp))
        ni(k) = min(max(5.38e7*temp,1.e3),1.e6)
        rni(k)= max(R2, ni(k)*rho(k))
        if (rqi(k).gt.R1 .and. rni(k).gt.R2) has_qi = .true.
        ! for snow
        rqs(k) = max(R1, qs(k)*rho(k))
        if (rqs(k).gt.R1) has_qs = .true.
      enddo

      if (has_qc) then
        do k=kts,kte
          if (rqc(k).le.R1) CYCLE
          lamdac   = (pidnc*nc0/rqc(k))**obmr
          re_qc(k) =  max(2.51E-6,min(1.5*(1.0/lamdac),50.E-6))
        enddo
      endif

     if (has_qi) then
        do k=kts,kte
          if (rqi(k).le.R1 .or. rni(k).le.R2) CYCLE
          diai = 11.9*sqrt(rqi(k)/ni(k))
          re_qi(k) = max(10.01E-6,min(0.75*0.163*diai,125.E-6))
        enddo
      endif

      if (has_qs) then
        do k=kts,kte
          if (rqs(k).le.R1) CYCLE
          supcol = t0c-t(k)
          n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs(k)))
          re_qs(k) = max(25.E-6,min(0.5*(1./lamdas), 999.E-6))
        enddo
      endif

      end subroutine effectRad_wsm5
!-----------------------------------------------------------------------
!-srf 
 subroutine vssqrt(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=sqrt(x(j))
   10 continue
      return
      end
 subroutine vsrec(y,x,n)
      real*4 x(*),y(*)
      do 10 j=1,n
      y(j)=1.e0/x(j)
   10 continue
      return
      end
!-srf
!-----------------------------------------------------------------------
                                                         
END MODULE module_mp_wsm5
!#endif
