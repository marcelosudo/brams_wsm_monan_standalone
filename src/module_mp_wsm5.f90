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
!$acc routine(cpmcal)
!$acc routine(xlcall)
!!$acc routine(provaconc3d)

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
             rsloperbmax,rslopesbmax,     &
             rsloper2max,rslopes2max,     &
             rsloper3max,rslopes3max

!newCode begin
!$acc declare create(qc0, qck1, pidnc)
!$acc declare create(bvtr1,bvtr2,bvtr3,bvtr4,g1pbr)
!$acc declare create(g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr)
!$acc declare create(precr1,precr2,xmmax,roqimax,bvts1)
!$acc declare create(bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs)
!$acc declare create(g5pbso2,pvts,pacrs,precs1,precs2,pidn0r)
!$acc declare create(pidn0s,xlv1,pacrc,pi)
!$acc declare create(rslopermax,rslopesmax,rslopegmax)
!$acc declare create(rsloperbmax,rslopesbmax)
!$acc declare create(rsloper2max,rslopes2max)
!$acc declare create(rsloper3max,rslopes3max)
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

  REAL, DIMENSION(kms:kme),                 &
        INTENT(INOUT) ::                                          &
                                                             th, pii,  &
                                                              q,  &
                                                              qc, &
                                                              qi, &
                                                              qr, &
                                                              qs
  REAL, DIMENSION(kms:kme),                 &
        INTENT(IN   ) ::                                          &
                                                             den, &
                                                             ! pii, &
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
  ! REAL, DIMENSION(kms:kme) ::   t
  ! REAL, DIMENSION( its:ite , kts:kte, 2 ) ::   qci, qrs
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

! #ifndef XEON_OPTIMIZED_WSM5
! print *,"ifndef XEON_OPTIMIZED_WSM5 OK"
      !! DO j=jts,jte
	  j=jts ! BLOCK 1
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
					th(kts:kte), pii(kts:kte), q(kms),& !(private - todas)
					qc(kts:kte),qi(kts:kte),& !(private - todas)
					qr(kts:kte), qs(kts:kte) & !(private - todas)
                    ,den(kms)                                & !(private - todas)
                    ,p(kms), delz(kms),                & !(private - todas)
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
		! BLOCK 1
! +---+-----------------------------------------------------------------+
		! BLOCK 2
         IF ( PRESENT (diagflag) ) THEN
         if (diagflag .and. do_radar_ref == 1) then
      ! WRITE(emess,*)'calling refl10cm_wsm5 ',its, jts
      ! CALL wrf_debug ( 0, emess )
            ! DO I=its,ite
               ! DO K=kts,kte
                  ! t1d(k)=th(i,k,j)*pii(i,k,j)
                  ! p1d(k)=p(i,k,j)
                  ! qv1d(k)=q(i,k,j)
                  ! qr1d(k)=qr(i,k,j)
                  ! qs1d(k)=qs(i,k,j)
               ! ENDDO
			   ! print *,"call refl10cm_wsm5 OK"
               call refl10cm_wsm5_x (q(kts:kte), qr(kts:kte), qs(kts:kte),  &
                       th(kts:kte), pii(kts:kte), p(kts:kte), dBZ, kts, kte, i, j)
			  ! call refl10cm_wsm5 (qv1d, qr1d, qs1d,                    &
                      ! t1d, p1d, dBZ, kts, kte, i, j)
               do k = kts, kte
                  refl_10cm(i,k,j) = MAX(-35., dBZ(k))
               enddo
            ! ENDDO
         endif
         ENDIF
		 ! BLOCK 2

		 ! BLOCK 3
         if (has_reqc.ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) then
           ! do i=its,ite
             do k=kts,kte
               re_qc(k) = RE_QC_BG
               re_qi(k) = RE_QI_BG 
               re_qs(k) = RE_QS_BG

               ! t1d(k)  = th(i,k,j)*pii(i,k,j)
               ! den1d(k)= den(i,k,j)
               ! qc1d(k) = qc(i,k,j)
               ! qi1d(k) = qi(i,k,j)
               ! qs1d(k) = qs(i,k,j)
             enddo
             call effectRad_wsm5_x(th(kts:kte), pii(kts:kte), qc(kts:kte), qi(kts:kte), qs(kts:kte), den(kts:kte),                 &
                                 qmin, t0c, re_qc, re_qi, re_qs,               &
                                 kts, kte, i, j)
								 
		     ! call effectRad_wsm5(t1d, qc1d, qi1d, qs1d, den1d,                 &
                                ! qmin, t0c, re_qc, re_qi, re_qs,               &
                                ! kts, kte, i, j)
             do k=kts,kte
               re_cloud(i,k,j) = MAX(RE_QC_BG, MIN(re_qc(k),  50.E-6))
               re_ice(i,k,j)   = MAX(RE_QI_BG, MIN(re_qi(k), 125.E-6))
               re_snow(i,k,j)  = MAX(RE_QS_BG, MIN(re_qs(k), 999.E-6))
             enddo
           ! enddo
         endif     ! has_reqc, etc...
		 ! BLOCK 3
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
  REAL, DIMENSION(kts:kte ),                           &
        INTENT(INOUT) ::                                          &
                                                              th, pii, qc, qi, qr, qs
  REAL, DIMENSION(kts:kte ) ::                         t, &
		                                                     xl, &
                                                             cpm
  REAL, DIMENSION(kms:kme ),                           &
        INTENT(INOUT) ::                                          &
                                                               q
  REAL, DIMENSION(kms:kme ),                           &
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
  REAL, DIMENSION(kts:kte , 2) ::                      &
                                                              rh, &
                                                              qss!, &
															! falk, &
                                                            ! fall
  REAL, DIMENSION( its:ite , kts:kte , 2) ::                      &
                                                          rslope, &
                                                         rslope2, &
                                                         rslope3, &
                                                         rslopeb, &
                                                         qrs_tmp, &
                                                           work1
  ! REAL, DIMENSION(kts:kte ) ::    &
                                                           ! falkc, &
                                                           ! fallc, &
                                                           ! xni
  REAL, DIMENSION( its:ite , kts:kte ) ::                         &
                                                          denfac, &
                                                         denqrs1, &
                                                         denqrs2, &
                                                          denqci, &
                                                          n0sfac, &
                                                           work2, &
                                                           workr, &
                                                           works, &
                                                          work1c, &
                                                          work2c
  REAL, DIMENSION( its:ite ) ::                                   &
                                                         delqrs1, &
                                                         delqrs2, &
                                                           delqi
  ! REAL, DIMENSION(kts:kte ) ::                         &
                                                           ! pigen, &
                                                           ! pidep, &
                                                           ! psdep, &
                                                           ! praut, &
                                                           ! psaut, &
                                                           ! prevp, &
                                                           ! psevp, &
                                                           ! pracw, &
                                                           ! psacw, &
                                                           ! psaci, &
                                                           ! pcond, &
                                                           ! psmlt
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
! return
!=================================================================
!   compute internal functions
!
      cpmcal(x) = cpd*(1.-max(x,qmin))+max(x,qmin)*cpv
      xlcal(x) = xlv0-xlv1*(x-t0c)
!----------------------------------------------------------------
!     diffus: diffusion coefficient of the water vapor
!     viscos: kinematic viscosity(m2s-1)
!     diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y        ! 8.794e-5*x**1.81/y
!     viscos(x,y) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/y  ! 1.496e-6*x**1.5/(x+120.)/y
!     xka(x,y) = 1.414e3*viscos(x,y)*y
!     diffac(a,b,c,d,e) = d*a*a/(xka(c,d)*rv*c*c)+1./(e*diffus(c,b))
!     venfac(a,b,c) = exp(log((viscos(b,c)/diffus(b,a)))*((.3333333)))         &
!                    /sqrt(viscos(b,c))*sqrt(sqrt(den0/c))
!     conden(a,b,c,d,e) = (max(b,qmin)-c)/(1.+d*d/(rv*e)*c/(a*a))
!
!----------------------------------------------------------------
      idim = ite-its+1
      kdim = kte-kts+1
!
!----------------------------------------------------------------
!     paddint 0 for negative values generated by dynamics
! BLOCK 1
      do k = kts, kte
        ! do i = its, ite
          ! qc(k) = max(qc(k),0.0) !15/02/2024 - desnecessário
          ! qr(k) = max(qr(k),0.0)
          ! qi(k) = max(qi(k),0.0)
          ! qs(k) = max(qs(k),0.0)
		  t(k)=th(k)*pii(k) !10/12/2023
        ! enddo
      enddo
! /BLOCK 1
!----------------------------------------------------------------
!     latent heat for phase changes and heat capacity. neglect the
!     changes during microphysical process calculation
!     emanuel(1994)
! BLOCK 2
      do k = kts, kte
        ! do i = its, ite
          cpm(k) = cpmcal(q(k))
          xl(k) = xlcal(t(k))		  
        ! enddo
      enddo
! /BLOCK 2
! BLOCK 3
      ! do k = kts, kte
        ! do i = its, ite
          ! delz_tmp(i,k) = delz(k) !! 19/02 USAR variáveis diretamente, e não tmp 
          ! den_tmp(i,k) = den(k)
        ! enddo
      ! enddo
! / BLOCK 3
!----------------------------------------------------------------
!    initialize the surface rain, snow
! BLOCK 4
      ! do i = its, ite
		i = 1
        rainncv(i) = 0.
        if(PRESENT (snowncv) .AND. PRESENT (snow)) snowncv(i,lat) = 0.
        sr(i) = 0.
! new local array to catch step snow
        tstepsnow(i) = 0.
      ! enddo
! /BLOCK 4
!----------------------------------------------------------------
!     compute the minor time steps.
!
      loops = max(nint(delt/dtcldcr),1)
      dtcld = delt/loops
      if(delt.le.dtcldcr) dtcld = delt
	! print*,'loops',loops - sempre 1
     ! do loop = 1,loops DESCOMENTAR
!
!----------------------------------------------------------------
!     initialize the large scale variables
! BLOCK 5
      ! do i = its, ite
        mstep(i) = 1
        flgcld(i) = .true.
      ! enddo
!/BLOCK 5
!BLOCK 6 !0.8575728
      do k = kts, kte
        CALL VSREC( tvec1(its), den(k), ite-its+1)
        do i = its, ite
          tvec1(i) = tvec1(i)*den0
        enddo
        CALL VSSQRT( denfac(its,k), tvec1(its), ite-its+1)
      enddo
!/BLOCK 6
!! Inline expansion for fpvs
!!         qss(i,k,1) = fpvs(t(i,k),0,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
!!         qss(i,k,2) = fpvs(t(i,k),1,rd,rv,cpv,cliq,cice,xlv0,xls,psat,t0c)
      hsub = xls
      hvap = xlv0
      cvap = cpv
      ttp=t0c+0.01
      dldt=cvap-cliq
      xa=-dldt/rv
      xb=xa+hvap/(rv*ttp)
      dldti=cvap-cice
      xai=-dldti/rv
      xbi=xai+hsub/(rv*ttp)

!BLOCK 7 !0.8629489 
      do k = kts, kte
        ! do i = its, ite
          tr=ttp/t(k)
          logtr=log(tr)
          qss(k,1)=psat*exp(logtr*(xa)+xb*(1.-tr))
          qss(k,1) = min(qss(k,1),0.99*p(k))
          qss(k,1) = ep2 * qss(k,1) / (p(k) - qss(k,1))
          qss(k,1) = max(qss(k,1),qmin)
          rh(k,1) = max(q(k) / qss(k,1),qmin)
          if(t(k).lt.ttp) then
            qss(k,2)=psat*exp(logtr*(xai)+xbi*(1.-tr))
          else
            qss(k,2)=psat*exp(logtr*(xa)+xb*(1.-tr))
          endif
          qss(k,2) = min(qss(k,2),0.99*p(k))
          qss(k,2) = ep2 * qss(k,2) / (p(k) - qss(k,2))
          qss(k,2) = max(qss(k,2),qmin)
          rh(k,2) = max(q(k) / qss(k,2),qmin)
        ! enddo
      enddo
!/BLOCK 7
!----------------------------------------------------------------
!     initialize the variables for microphysical physics
!
!BLOCK 8 !0.8294189
      !!! do k = kts, kte    !02/03/2024 variáveis zeradas antes do loop, na inicialização
        !!! do i = its, ite
          ! prevp(kts:kte) = 0.
          ! psdep(kts:kte) = 0.
          ! praut(kts:kte) = 0.
          ! psaut(kts:kte) = 0.
          ! pracw(kts:kte) = 0.
          ! psaci(kts:kte) = 0.
          ! psacw(kts:kte) = 0.
          ! pigen(kts:kte) = 0.
          ! pidep(kts:kte) = 0.
          ! pcond(kts:kte) = 0.
          ! psmlt(kts:kte) = 0.
          ! psevp(kts:kte) = 0.
          ! falk(kts:kte,1) = 0. !BAIXOU p/ 6.804927 sem essas 4 variáveis
          ! falk(kts:kte,2) = 0. ! estava em 7,773113
          ! fall(kts:kte,1) = 0. 
          ! fall(kts:kte,2) = 0.
          ! fallc(kts:kte) = 0.
          ! falkc(kts:kte) = 0.
          ! xni(kts:kte) = 1.e3
        !!! enddo
      !!! enddo
!/BLOCK 8
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
!BLOCK 9 !0.8618581
      ! do k = kts, kte
        ! do i = its, ite
          ! temp = (den(k)*max(qci(i,k,2),qmin))
          ! temp = sqrt(sqrt(temp*temp*temp))
          ! xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
        ! enddo
      ! enddo
!/BLOCK 9
!----------------------------------------------------------------
!     compute the fallout term:
!     first, vertical terminal velosity for minor loops
!----------------------------------------------------------------
!BLOCK 10 !0.8670120
      ! do k = kts, kte
        ! do i = its, ite
          ! qrs_tmp(i,k,1) = qrs(i,k,1)
          ! qrs_tmp(i,k,2) = qrs(i,k,2)
        ! enddo
      ! enddo
!/BLOCK 10
      ! call slope_wsm5(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     ! work1,its,ite,kts,kte) DESCOMENTAR
!BLOCK 11 !0.8136392
      ! do k = kte, kts, -1
        ! do i = its, ite
          ! workr(i,k) = work1(i,k,1) 
          ! works(i,k) = work1(i,k,2) 
          ! denqrs1(i,k) = den(k)*qrs(i,k,1)
          ! denqrs2(i,k) = den(k)*qrs(i,k,2)
          ! if(qrs(i,k,1).le.0.0) workr(i,k) = 0.0
          ! if(qrs(i,k,2).le.0.0) works(i,k) = 0.0
        ! enddo
      ! enddo
!/BLOCK 11
	! t = 1.
      ! call nislfv_rain_plm(idim,kdim,den_tmp,denfac,t,delz_tmp,workr,denqrs1,  &
                           ! delqrs1,dtcld,1,1) DESCOMENTAR
      ! call nislfv_rain_plm(idim,kdim,den_tmp,denfac,t,delz_tmp,works,denqrs2,  &
                           ! delqrs2,dtcld,2,1) DESCOMENTAR
!BLOCK 12 !0.8296371 BUG rrp, rtp, thp
      ! do k = kts, kte
        ! do i = its, ite
          ! qrs(i,k,1) = max(denqrs1(i,k)/den(k),0.)
          ! qrs(i,k,2) = max(denqrs2(i,k)/den(k),0.)
          ! fall(i,k,1) = denqrs1(i,k)*workr(i,k)/delz(k)
          ! fall(i,k,2) = denqrs2(i,k)*works(i,k)/delz(k)
        ! enddo
      ! enddo
!/BLOCK 12
!BLOCK 13 !0.8716211 - VALIDADO
      ! do i = its, ite
        ! fall(i,1,1) = delqrs1(i)/delz(i,1)/dtcld
        ! fall(i,1,2) = delqrs2(i)/delz(i,1)/dtcld
      ! enddo
!/BLOCK 13
!BLOCK 14 !0.7921262 - VALIDADO
      ! do k = kts, kte
        ! do i = its, ite
          ! qrs_tmp(i,k,1) = qrs(i,k,1)
          ! qrs_tmp(i,k,2) = qrs(i,k,2)
        ! enddo
      ! enddo
!/BLOCK 14
      ! call slope_wsm5(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     ! work1,its,ite,kts,kte) ! VALIDADO
!BLOCK 15 !0.8217640 !! BUG rrp
      ! do k = kte, kts, -1
        ! do i = its, ite
          ! supcol = t0c-t(i,k)
          ! n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          ! if(t(i,k).gt.t0c.and.qrs(i,k,2).gt.0.) then
! ----------------------------------------------------------------
! psmlt: melting of snow [HL A33] [RH83 A25]
      ! (T>T0: S->R)
! ----------------------------------------------------------------
            ! xlf = xlf0
            ! work2(i,k)= (exp(log(((1.496e-6*((t(i,k))*sqrt(t(i,k)))        &
                        ! /((t(i,k))+120.)/(den(k)))/(8.794e-5             &
                        ! *exp(log(t(i,k))*(1.81))/p(k))))                 &
                        ! *((.3333333)))/sqrt((1.496e-6*((t(i,k))            &
                        ! *sqrt(t(i,k)))/((t(i,k))+120.)/(den(k))))        &
                        ! *sqrt(sqrt(den0/(den(k)))))
            ! coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
            ! psmlt(i,k) = (1.414e3*(1.496e-6*((t(i,k))*sqrt(t(i,k)))        &
                        ! /((t(i,k))+120.)/(den(k)) )*(den(k)))          &
                        ! /xlf*(t0c-t(i,k))*pi/2.                            &
                        ! *n0sfac(i,k)*(precs1*rslope2(i,k,2)+precs2         &
                        ! *work2(i,k)*coeres)/den(k)
            ! psmlt(i,k) = min(max(psmlt(i,k)*dtcld/mstep(i),                &
                        ! -qrs(i,k,2)/mstep(i)),0.)
            ! qrs(i,k,2) = qrs(i,k,2) + psmlt(i,k)
            ! qrs(i,k,1) = qrs(i,k,1) - psmlt(i,k)
            ! t(i,k) = t(i,k) + xlf/cpm(i,k)*psmlt(i,k)
          ! endif
        ! enddo
      ! enddo
!/BLOCK 15
!---------------------------------------------------------------
! Vice [ms-1] : fallout of ice crystal [HDC 5a]
!---------------------------------------------------------------
!BLOCK 16 !0.8442681 
      ! do k = kte, kts, -1
        ! do i = its, ite
          ! if(qci(i,k,2).le.0.) then
            ! work1c(i,k) = 0.
          ! else
            ! xmi = den(k)*qci(i,k,2)/xni(i,k)
            ! diameter  = max(min(dicon * sqrt(xmi),dimax), 1.e-25)
            ! work1c(i,k) = 1.49e4*exp(log(diameter)*(1.31))
          ! endif
        ! enddo
      ! enddo
!/BLOCK 16
!  forward semi-laglangian scheme (JH), PCM (piecewise constant),  (linear)
!BLOCK 17 !0.8326032
      ! do k = kte, kts, -1
        ! do i = its, ite
          ! denqci(i,k) = den(k)*qci(i,k,2)
        ! enddo
      ! enddo
!/BLOCK 17
      ! call nislfv_rain_plm(idim,kdim,den_tmp,denfac,t,delz_tmp,work1c,denqci,  &
                           ! delqi,dtcld,1,0)
!BLOCK 18 !0.8687339
      ! do k = kts, kte
        ! do i = its, ite
          ! qci(i,k,2) = max(denqci(i,k)/den(k),0.)
        ! enddo
      ! enddo
!/BLOCK 18
!BLOCK 19 !0.8241761
      ! do i = its, ite
        ! fallc(i,1) = delqi(i)/delz(i,1)/dtcld
      ! enddo
!/BLOCK 19
!----------------------------------------------------------------
!      rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
!BLOCK 20 !0.8583159
      ! do i = its, ite
        ! fallsum = fall(i,1,1)+fall(i,1,2)+fallc(i,1)
        ! fallsum_qsi = fall(i,1,2)+fallc(i,1)
        ! if(fallsum.gt.0.) then
          ! rainncv(i) = fallsum*delz(i,1)/denr*dtcld*1000. + rainncv(i)
          ! rain(i) = fallsum*delz(i,1)/denr*dtcld*1000. + rain(i)
        ! endif
        ! if(fallsum_qsi.gt.0.) then
          ! tstepsnow(i)   = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            &
                           ! +tstepsnow(i) 
        ! IF ( PRESENT (snowncv) .AND. PRESENT (snow)) THEN
          ! snowncv(i,lat) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000.            & 
                           ! +snowncv(i,lat)    
          ! snow(i,lat) = fallsum_qsi*delz(i,kts)/denr*dtcld*1000. + snow(i,lat)
        ! ENDIF
        ! endif
        ! IF ( PRESENT (snowncv) ) THEN
          ! if(fallsum.gt.0.)sr(i)=snowncv(i,lat)/(rainncv(i)+1.e-12)
        ! ELSE
          ! if(fallsum.gt.0.)sr(i)=tstepsnow(i)/(rainncv(i)+1.e-12)
        ! ENDIF
      ! enddo
!/BLOCK 20
!---------------------------------------------------------------
! pimlt: instantaneous melting of cloud ice [HL A47] [RH83 A28]
!       (T>T0: I->C)
!---------------------------------------------------------------
!BLOCK 21 !0.8219631
      ! do k = kts, kte
        ! do i = its, ite
          ! supcol = t0c-t(i,k)
          ! xlf = xls-xl(i,k)
          ! if(supcol.lt.0.) xlf = xlf0
          ! if(supcol.lt.0.and.qci(i,k,2).gt.0.) then
            ! qci(i,k,1) = qci(i,k,1) + qci(i,k,2)
            ! t(i,k) = t(i,k) - xlf/cpm(i,k)*qci(i,k,2)
            ! qci(i,k,2) = 0.
          ! endif
!---------------------------------------------------------------
! pihmf: homogeneous freezing of cloud water below -40c [HL A45]
!        (T<-40C: C->I)
!---------------------------------------------------------------
          ! if(supcol.gt.40..and.qci(i,k,1).gt.0.) then
            ! qci(i,k,2) = qci(i,k,2) + qci(i,k,1)
            ! t(i,k) = t(i,k) + xlf/cpm(i,k)*qci(i,k,1)
            ! qci(i,k,1) = 0.
          ! endif
!---------------------------------------------------------------
! pihtf: heterogeneous freezing of cloud water [HL A44]
!        (T0>T>-40C: C->I)
!---------------------------------------------------------------
          ! if(supcol.gt.0..and.qci(i,k,1).gt.0.) then
            ! supcolt=min(supcol,50.)
!           pfrzdtc = min(pfrz1*(exp(pfrz2*supcol)-1.)                         &
!              *den(k)/denr/xncr*qci(i,k,1)**2*dtcld,qci(i,k,1))
            ! pfrzdtc = min(pfrz1*(exp(pfrz2*supcolt)-1.)                        &
            ! *den(k)/denr/xncr*qci(i,k,1)*qci(i,k,1)*dtcld,qci(i,k,1))
            ! qci(i,k,2) = qci(i,k,2) + pfrzdtc
            ! t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc
            ! qci(i,k,1) = qci(i,k,1)-pfrzdtc
          ! endif
!---------------------------------------------------------------
! psfrz: freezing of rain water [HL A20] [LFO 45]
!        (T<T0, R->S)
!---------------------------------------------------------------
          ! if(supcol.gt.0..and.qrs(i,k,1).gt.0.) then
            ! supcolt=min(supcol,50.)
            ! temp = rslope(i,k,1)
            ! temp = temp*temp*temp*temp*temp*temp*temp
            ! pfrzdtr = min(20.*(pi*pi)*pfrz1*n0r*denr/den(k)                  &
                  ! *(exp(pfrz2*supcolt)-1.)*temp*dtcld,                         &
                  ! qrs(i,k,1))
            ! qrs(i,k,2) = qrs(i,k,2) + pfrzdtr
            ! t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            ! qrs(i,k,1) = qrs(i,k,1)-pfrzdtr
          ! endif
        ! enddo
      ! enddo
!/BLOCK 21
!----------------------------------------------------------------
!     update the slope parameters for microphysics computation
!BLOCK 22 !0.8470581
      ! do k = kts, kte
        ! do i = its, ite
          ! qrs_tmp(i,k,1) = qrs(i,k,1)
          ! qrs_tmp(i,k,2) = qrs(i,k,2)
        ! enddo
      ! enddo
!/BLOCK 22
      ! call slope_wsm5(qrs_tmp,den_tmp,denfac,t,rslope,rslopeb,rslope2,rslope3, &
                     ! work1,its,ite,kts,kte)
!------------------------------------------------------------------
!     work1:  the thermodynamic term in the denominator associated with
!             heat conduction and vapor diffusion
!             (ry88, y93, h85)
!     work2: parameter associated with the ventilation effects(y93)
!BLOCK 23 !0.8421991
      ! do k = kts, kte
        ! do i = its, ite
          ! work1(i,k,1) = ((((den(i,k))*(xl(i,k))*(xl(i,k)))*((t(i,k))+120.)    &
                       ! *(den(i,k)))/(1.414e3*(1.496e-6*((t(i,k))*sqrt(t(i,k))))&
                       ! *(den(i,k))*(rv*(t(i,k))*(t(i,k)))))                    &
                      ! +  p(k)/((qss(k,1))*(8.794e-5*exp(log(t(i,k))*(1.81))))
          ! work1(i,k,2) = ((((den(i,k))*(xls)*(xls))*((t(i,k))+120.)*(den(i,k)))&
                       ! /(1.414e3*(1.496e-6*((t(i,k))*sqrt(t(i,k))))*(den(i,k)) &
                       ! *(rv*(t(i,k))*(t(i,k))))                                &
                      ! + p(k)/(qss(k,2)*(8.794e-5*exp(log(t(i,k))*(1.81)))))
          ! work2(i,k) = (exp(.3333333*log(((1.496e-6 * ((t(i,k))*sqrt(t(i,k)))) &
                     ! *p(k))/(((t(i,k))+120.)*den(i,k)*(8.794e-5              &
                     ! *exp(log(t(i,k))*(1.81))))))*sqrt(sqrt(den0/(den(i,k))))) &
                      ! /sqrt((1.496e-6*((t(i,k))*sqrt(t(i,k))))                 &
                      ! /(((t(i,k))+120.)*den(i,k)))
        ! enddo 
      ! enddo 
!/BLOCK 23
!===============================================================
!
! warm rain processes
!
! - follows the processes in RH83 and LFO except for autoconcersion
!
!===============================================================
!BLOCK 24 !0.8380830
      ! do k = kts, kte
        ! do i = its, ite
          ! supsat = max(q(i,k),qmin)-qss(k,1)
          ! satdt = supsat/dtcld
!---------------------------------------------------------------
! praut: auto conversion rate from cloud to rain [HDC 16]
!        (C->R)
!---------------------------------------------------------------
          ! if(qci(i,k,1).gt.qc0) then
            ! praut(i,k) = qck1*exp(log(qci(i,k,1))*((7./3.)))
            ! praut(i,k) = min(praut(i,k),qci(i,k,1)/dtcld)
          ! endif
!---------------------------------------------------------------
! pracw: accretion of cloud water by rain [HL A40] [LFO 51]
!        (C->R)
!---------------------------------------------------------------
          ! if(qrs(i,k,1).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            ! pracw(i,k) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1)               & 
                         ! *qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          ! endif
!---------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [HDC 14]
!        (V->R or R->V)
!---------------------------------------------------------------
          ! if(qrs(i,k,1).gt.0.) then
            ! coeres = rslope2(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            ! prevp(i,k) = (rh(k,1)-1.)*(precr1*rslope2(i,k,1)                 &
                         ! +precr2*work2(i,k)*coeres)/work1(i,k,1)
            ! if(prevp(i,k).lt.0.) then
              ! prevp(i,k) = max(prevp(i,k),-qrs(i,k,1)/dtcld)
              ! prevp(i,k) = max(prevp(i,k),satdt/2)
            ! else
              ! prevp(i,k) = min(prevp(i,k),satdt/2)
            ! endif
          ! endif
        ! enddo
      ! enddo
!/BLOCK 24
!===============================================================
!
! cold rain processes
!
! - follows the revised ice microphysics processes in HDC
! - the processes same as in RH83 and RH84  and LFO behave
!   following ice crystal hapits defined in HDC, inclduing
!   intercept parameter for snow (n0s), ice crystal number
!   concentration (ni), ice nuclei number concentration
!   (n0i), ice diameter (d)
!
!===============================================================
!BLOCK 25 !0.8388951
      ! rdtcld = 1./dtcld
      ! do k = kts, kte
        ! do i = its, ite
          ! supcol = t0c-t(i,k)
          ! n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),1.)
          ! supsat = max(q(i,k),qmin)-qss(k,2)
          ! satdt = supsat/dtcld
          ! ifsat = 0
!-------------------------------------------------------------
! Ni: ice crystal number concentraiton   [HDC 5c]
!-------------------------------------------------------------
          ! temp = (den(i,k)*max(qci(i,k,2),qmin))
          ! temp = sqrt(sqrt(temp*temp*temp))
          ! xni(i,k) = min(max(5.38e7*temp,1.e3),1.e6)
          ! eacrs = exp(0.07*(-supcol))

          ! if(supcol.gt.0) then
            ! if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,2).gt.qmin) then
              ! xmi = den(i,k)*qci(i,k,2)/xni(i,k)
              ! diameter  = min(dicon * sqrt(xmi),dimax)
              ! vt2i = 1.49e4*diameter**1.31
              ! vt2s = pvts*rslopeb(i,k,2)*denfac(i,k)
!-------------------------------------------------------------
! psaci: Accretion of cloud ice by rain [HDC 10]
!        (T<T0: I->S)
!-------------------------------------------------------------
              ! acrfac = 2.*rslope3(i,k,2)+2.*diameter*rslope2(i,k,2)            &
                      ! +diameter**2*rslope(i,k,2)
              ! psaci(i,k) = pi*qci(i,k,2)*eacrs*n0s*n0sfac(i,k)                 &
                           ! *abs(vt2s-vt2i)*acrfac/4.
            ! endif
          ! endif
!-------------------------------------------------------------
! psacw: Accretion of cloud water by snow  [HL A7] [LFO 24]
!        (T<T0: C->S, and T>=T0: C->R)
!-------------------------------------------------------------
          ! if(qrs(i,k,2).gt.qcrmin.and.qci(i,k,1).gt.qmin) then
            ! psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)                  &
                           ! *rslopeb(i,k,2)*qci(i,k,1)*denfac(i,k)              &
!                          ,qci(i,k,1)/dtcld)
                           ! ,qci(i,k,1)*rdtcld)
          ! endif
          ! if(supcol .gt. 0) then
!-------------------------------------------------------------
! pidep: Deposition/Sublimation rate of ice [HDC 9]
!       (T<T0: V->I or I->V)
!-------------------------------------------------------------
            ! if(qci(i,k,2).gt.0.and.ifsat.ne.1) then
              ! xmi = den(i,k)*qci(i,k,2)/xni(i,k)
              ! diameter = dicon * sqrt(xmi)
              ! pidep(i,k) = 4.*diameter*xni(i,k)*(rh(k,2)-1.)/work1(i,k,2)
              ! supice = satdt-prevp(i,k)
              ! if(pidep(i,k).lt.0.) then
                ! pidep(i,k) = max(max(pidep(i,k),satdt*.5),supice)
                ! pidep(i,k) = max(pidep(i,k),-qci(i,k,2)*rdtcld)
              ! else
                ! pidep(i,k) = min(min(pidep(i,k),satdt*.5),supice)
              ! endif
              ! if(abs(prevp(i,k)+pidep(i,k)).ge.abs(satdt)) ifsat = 1
            ! endif
!-------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [HDC 14]
!        (V->S or S->V)
!-------------------------------------------------------------
            ! if(qrs(i,k,2).gt.0..and.ifsat.ne.1) then
              ! coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              ! psdep(i,k) = (rh(k,2)-1.)*n0sfac(i,k)                          &
                           ! *(precs1*rslope2(i,k,2)+precs2                      &
                           ! *work2(i,k)*coeres)/work1(i,k,2)
              ! supice = satdt-prevp(i,k)-pidep(i,k)
              ! if(psdep(i,k).lt.0.) then
                ! psdep(i,k) = max(psdep(i,k),-qrs(i,k,2)*rdtcld)
                ! psdep(i,k) = max(max(psdep(i,k),satdt*.5),supice)
              ! else
                ! psdep(i,k) = min(min(psdep(i,k),satdt*.5),supice)
              ! endif
              ! if(abs(prevp(i,k)+pidep(i,k)+psdep(i,k)).ge.abs(satdt))          &
                ! ifsat = 1
            ! endif
!-------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [HL A50] [HDC 7-8]
!       (T<T0: V->I)
!-------------------------------------------------------------
            ! if(supsat.gt.0.and.ifsat.ne.1) then
              ! supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
              ! xni0 = 1.e3*exp(0.1*supcol)
              ! roqi0 = 4.92e-11*exp(log(xni0)*(1.33))
              ! pigen(i,k) = max(0.,(roqi0/den(i,k)-max(qci(i,k,2),0.))          &
                         ! *rdtcld)
              ! pigen(i,k) = min(min(pigen(i,k),satdt),supice)
            ! endif
!
!-------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [HDC 12]
!       (T<T0: I->S)
!-------------------------------------------------------------
            ! if(qci(i,k,2).gt.0.) then
              ! qimax = roqimax/den(k)
              ! psaut(i,k) = max(0.,(qci(i,k,2)-qimax)*rdtcld)
            ! endif
          ! endif
!-------------------------------------------------------------
! psevp: Evaporation of melting snow [HL A35] [RH83 A27]
!       (T>T0: S->V)
!-------------------------------------------------------------
          ! if(supcol.lt.0.) then
            ! if(qrs(i,k,2).gt.0..and.rh(k,1).lt.1.)                           &
              ! psevp(i,k) = psdep(i,k)*work1(i,k,2)/work1(i,k,1)
              ! psevp(i,k) = min(max(psevp(i,k),-qrs(i,k,2)*rdtcld),0.)
          ! endif
        ! enddo
      ! enddo
!/BLOCK 25
!
!----------------------------------------------------------------
!     check mass conservation of generation terms and feedback to the
!     large scale
!BLOCK 26 !0.8072529
      ! do k = kts, kte
        ! do i = its, ite
          ! if(t(i,k).le.t0c) then

!    cloud water

            ! value = max(qmin,qci(i,k,1))
            ! source = (praut(i,k)+pracw(i,k)+psacw(i,k))*dtcld
            ! if (source.gt.value) then
              ! factor = value/source
              ! praut(i,k) = praut(i,k)*factor
              ! pracw(i,k) = pracw(i,k)*factor
              ! psacw(i,k) = psacw(i,k)*factor
            ! endif

!    cloud ice

            ! value = max(qmin,qci(i,k,2))
            ! source = (psaut(i,k)+psaci(i,k)-pigen(i,k)-pidep(i,k))*dtcld
            ! if (source.gt.value) then
              ! factor = value/source
              ! psaut(i,k) = psaut(i,k)*factor
              ! psaci(i,k) = psaci(i,k)*factor
              ! pigen(i,k) = pigen(i,k)*factor
              ! pidep(i,k) = pidep(i,k)*factor
            ! endif

!    rain


            ! value = max(qmin,qrs(i,k,1))
            ! source = (-praut(i,k)-pracw(i,k)-prevp(i,k))*dtcld
            ! if (source.gt.value) then
              ! factor = value/source
              ! praut(i,k) = praut(i,k)*factor
              ! pracw(i,k) = pracw(i,k)*factor
              ! prevp(i,k) = prevp(i,k)*factor
            ! endif

!   snow

            ! value = max(qmin,qrs(i,k,2))
            ! source = (-psdep(i,k)-psaut(i,k)-psaci(i,k)-psacw(i,k))*dtcld  
            ! if (source.gt.value) then
              ! factor = value/source
              ! psdep(i,k) = psdep(i,k)*factor
              ! psaut(i,k) = psaut(i,k)*factor
              ! psaci(i,k) = psaci(i,k)*factor
              ! psacw(i,k) = psacw(i,k)*factor
            ! endif

            ! work2(i,k)=-(prevp(i,k)+psdep(i,k)+pigen(i,k)+pidep(i,k))
!    update
            ! q(i,k) = q(i,k)+work2(i,k)*dtcld
            ! qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                        ! +psacw(i,k))*dtcld,0.)
            ! qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                        ! +prevp(i,k))*dtcld,0.)
            ! qci(i,k,2) = max(qci(i,k,2)-(psaut(i,k)+psaci(i,k)                 &
                        ! -pigen(i,k)-pidep(i,k))*dtcld,0.)
            ! qrs(i,k,2) = max(qrs(i,k,2)+(psdep(i,k)+psaut(i,k)                 &
                        ! +psaci(i,k)+psacw(i,k))*dtcld,0.)
            ! xlf = xls-xl(i,k)
            ! xlwork2 = -xls*(psdep(i,k)+pidep(i,k)+pigen(i,k))                  &
                      ! -xl(i,k)*prevp(i,k)-xlf*psacw(i,k)
            ! t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          ! else

!    cloud water

            ! value = max(qmin,qci(i,k,1))
            ! source=(praut(i,k)+pracw(i,k)+psacw(i,k))*dtcld
            ! if (source.gt.value) then
              ! factor = value/source
              ! praut(i,k) = praut(i,k)*factor
              ! pracw(i,k) = pracw(i,k)*factor
              ! psacw(i,k) = psacw(i,k)*factor
            ! endif

!    rain

            ! value = max(qmin,qrs(i,k,1))
            ! source = (-praut(i,k)-pracw(i,k)-prevp(i,k)-psacw(i,k))*dtcld
            ! if (source.gt.value) then
              ! factor = value/source
              ! praut(i,k) = praut(i,k)*factor
              ! pracw(i,k) = pracw(i,k)*factor
              ! prevp(i,k) = prevp(i,k)*factor
              ! psacw(i,k) = psacw(i,k)*factor
            ! endif  

!    snow

            ! value = max(qcrmin,qrs(i,k,2))
            ! source=(-psevp(i,k))*dtcld
            ! if (source.gt.value) then
              ! factor = value/source
              ! psevp(i,k) = psevp(i,k)*factor
            ! endif
            ! work2(i,k)=-(prevp(i,k)+psevp(i,k))
!    update
            ! q(i,k) = q(i,k)+work2(i,k)*dtcld
            ! qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k)                 &
                        ! +psacw(i,k))*dtcld,0.)
            ! qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)                 &
                        ! +prevp(i,k) +psacw(i,k))*dtcld,0.)
            ! qrs(i,k,2) = max(qrs(i,k,2)+psevp(i,k)*dtcld,0.)
            ! xlf = xls-xl(i,k)
            ! xlwork2 = -xl(i,k)*(prevp(i,k)+psevp(i,k))
            ! t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          ! endif
        ! enddo
      ! enddo
!/BLOCK 26
! Inline expansion for fpvs
      ! hsub = xls DESCOMENTAR
      ! hvap = xlv0
      ! cvap = cpv
      ! ttp=t0c+0.01
      ! dldt=cvap-cliq
      ! xa=-dldt/rv
      ! xb=xa+hvap/(rv*ttp)
      ! dldti=cvap-cice
      ! xai=-dldti/rv
      ! xbi=xai+hsub/(rv*ttp)
!BLOCK 27 !0.8587410
      ! do k = kts, kte
      ! do i = its, ite
          ! tr=ttp/t(i,k)
          ! logtr = log(tr)
          ! qss(k,1)=psat*exp(logtr*(xa)+xb*(1.-tr))
          ! qss(k,1) = min(qss(k,1),0.99*p(k))
          ! qss(k,1) = ep2 * qss(k,1) / (p(k) - qss(k,1))
          ! qss(k,1) = max(qss(k,1),qmin)
        ! enddo
      ! enddo
!/BLOCK 27
!----------------------------------------------------------------
!  pcond: condensational/evaporational rate of cloud water [HL A46] [RH83 A6]
!     if there exists additional water vapor condensated/if
!     evaporation of cloud water is not enough to remove subsaturation
!BLOCK 28 !0.8564219
      ! do k = kts, kte
        ! do i = its, ite
          ! work1(i,k,1) = ((max(q(i,k),qmin)-(qss(k,1)))/(1.+(xl(i,k))         &   
                        ! *(xl(i,k))/(rv*(cpm(i,k)))*(qss(k,1))                 &
                        ! /((t(i,k))*(t(i,k)))))
          ! work2(i,k) = qci(i,k,1)+work1(i,k,1)
          ! pcond(i,k) = min(max(work1(i,k,1)/dtcld,0.),max(q(i,k),0.)/dtcld)
          ! if(qci(i,k,1).gt.0..and.work1(i,k,1).lt.0.)                          &
            ! pcond(i,k) = max(work1(i,k,1),-qci(i,k,1))/dtcld
          ! q(i,k) = q(i,k)-pcond(i,k)*dtcld
          ! qci(i,k,1) = max(qci(i,k,1)+pcond(i,k)*dtcld,0.)
          ! t(i,k) = t(i,k)+pcond(i,k)*xl(i,k)/cpm(i,k)*dtcld
        ! enddo
      ! enddo
!/BLOCK 28
!
!----------------------------------------------------------------
!     padding for small values
!BLOCK 29
      ! do k = kts, kte
        ! do i = its, ite
          ! if(qci(i,k,1).le.qmin) qci(i,k,1) = 0.0
          ! if(qci(i,k,2).le.qmin) qci(i,k,2) = 0.0
		  !! th(i,k)=t(i,k)/pii(i,k) !10/12/2023
        ! enddo
      ! enddo
!/BLOCK 29
   ! enddo                  ! big loops
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
! subroutine provaconc(kts, kte)
  ! INTEGER,      INTENT(IN   )    ::   kts,kte
  ! REAL, DIMENSION(kts:kte ) ::                         &
                                                           ! pigen, &
                                                           ! pidep, &
                                                           ! psdep, &
                                                           ! praut, &
                                                           ! psaut, &
                                                           ! prevp, &
                                                           ! psevp, &
                                                           ! pracw, &
                                                           ! psacw, &
                                                           ! psaci, &
                                                           ! pcond, &
                                                           ! psmlt
  ! REAL, DIMENSION(kts:kte , 2) ::                      &
															! falk, &
                                                            ! fall
  ! REAL, DIMENSION(kts:kte ) ::    &
                                                           ! falkc, &
                                                           ! fallc, &
                                                           ! xni
  !!$acc routine seq
  ! prevp(kts:kte) = 0.
  ! psdep(kts:kte) = 0.
  ! praut(kts:kte) = 0.
  ! psaut(kts:kte) = 0.
  ! pracw(kts:kte) = 0.
  ! psaci(kts:kte) = 0.
  ! psacw(kts:kte) = 0.
  ! pigen(kts:kte) = 0.
  ! pidep(kts:kte) = 0.
  ! pcond(kts:kte) = 0.
  ! psmlt(kts:kte) = 0.
  ! psevp(kts:kte) = 0.
  ! falk(kts:kte,1) = 0. !BAIXOU p/ 6.804927 sem essas 4 variáveis
  ! falk(kts:kte,2) = 0. ! estava em 7,773113
  ! fall(kts:kte,1) = 0. 
  ! fall(kts:kte,2) = 0.
  ! fallc(kts:kte) = 0.
  ! falkc(kts:kte) = 0.
  ! xni(kts:kte) = 1.e3
! end

subroutine provaconc3d(prevp,psdep,praut,psaut,pracw,psaci,psacw,pigen,pidep,pcond,psmlt, &
          psevp,falk1,falk2,fall1,fall2,fallc,falkc,xni, kts, kte, i, j,m1,m2,m3)
! subroutine provaconc3d(prevp, kts, kte, i, j,m1,m2,m3)
		  !$acc routine seq
  INTEGER,      INTENT(IN   )    ::   i,j, kts,kte, m1,m2,m3
    REAL, DIMENSION(:,:,:), INTENT(INOUT) ::                         &
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
  REAL, DIMENSION(:,:,:), INTENT(INOUT) ::                      &
															falk1,falk2, &
                                                            fall1,fall2
  REAL, DIMENSION(:,:,:), INTENT(INOUT) ::    &
                                                           falkc, &
                                                           fallc, &
                                                           xni
  INTEGER :: k
      do k = kts, kte
          prevp(k,i,j) = 0.
          psdep(k,i,j) = 0.
          praut(k,i,j) = 0.
          psaut(k,i,j) = 0.
          pracw(k,i,j) = 0.
          psaci(k,i,j) = 0.
          psacw(k,i,j) = 0.
          pigen(k,i,j) = 0.
          pidep(k,i,j) = 0.
          pcond(k,i,j) = 0.
          psmlt(k,i,j) = 0.
          psevp(k,i,j) = 0.
          falk1(k,i,j) = 0.
          falk2(k,i,j) = 0.
          fall1(k,i,j) = 0.
          fall2(k,i,j) = 0.
          fallc(k,i,j) = 0.
          falkc(k,i,j) = 0.
          xni(k,i,j) = 1.e3
      enddo
	  ! prevp(kts:kte,i,j) = 0.
end
		  
END MODULE module_mp_wsm5
!#endif
