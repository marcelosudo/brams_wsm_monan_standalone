program wsm_test
!  !$acc routine(wsm5, wsm52D)
!!!$acc routine(wsm5, wsm5x)
  !!$acc routine(refl10cm_wsm5)
!! !$acc routine(effectRad_wsm5_x)

   USE module_mp_radar
!!  use module_mp_wsm3, only: wsm3init, wsm3
  use module_mp_wsm5, only: wsm5init, wsm5
!!  use module_mp_wsm6, only: wsm6init, wsm6
!!  use module_mp_wsm7, only: wsm7init, wsm7
  use rconstants, only: p00,cpor,alvl,alvi,cpi,cpi4,cp253i

  implicit none     

!
   REAL, PARAMETER :: dtcldcr     = 120. ! maximum time step for minor loops
   REAL, PARAMETER :: n0r = 8.e6         ! intercept parameter rain
   REAL, PARAMETER :: avtr = 841.9       ! a constant for terminal velocity of rain
   REAL, PARAMETER :: bvtr = 0.8         ! a constant for terminal velocity of rain
   REAL, PARAMETER :: r0 = .8e-5         ! 8 microm  in contrast to 10 micro m
   REAL, PARAMETER :: peaut = .55        ! collection efficiency
   REAL, PARAMETER :: xncr = 3.e8        ! maritime cloud in contrast to 3.e8 in tc80
   REAL, PARAMETER :: xmyu = 1.718e-5    ! the dynamic viscosity kgm-1s-1
   REAL, PARAMETER :: avts = 11.72       ! a constant for terminal velocity of snow
   REAL, PARAMETER :: bvts = .41         ! a constant for terminal velocity of snow
   REAL, PARAMETER :: n0smax =  1.e11    ! maximum n0s (t=-90C unlimited)
   REAL, PARAMETER :: lamdarmax = 8.e4   ! limited maximum value for slope parameter of rain
   REAL, PARAMETER :: lamdasmax = 1.e5   ! limited maximum value for slope parameter of snow
   REAL, PARAMETER :: lamdagmax = 6.e4   ! limited maximum value for slope parameter of graupel
   REAL, PARAMETER :: dicon = 11.9       ! constant for the cloud-ice diamter
   REAL, PARAMETER :: dimax = 500.e-6    ! limited maximum value for the cloud-ice diamter
   REAL, PARAMETER :: n0s = 2.e6         ! temperature dependent intercept parameter snow
   REAL, PARAMETER :: alpha = .12        ! .122 exponen factor for n0s
   REAL, PARAMETER :: pfrz1 = 100.       ! constant in Biggs freezing
   REAL, PARAMETER :: pfrz2 = 0.66       ! constant in Biggs freezing
   REAL, PARAMETER :: qcrmin = 1.e-9     ! minimun values for qr, qs, and qg
   REAL, PARAMETER :: eacrc = 1.0        ! Snow/cloud-water collection efficiency
   REAL, SAVE ::                                      &
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
			 
  integer :: m1,m2,m3, mynum,ia,iz,ja,jz
  integer :: ids,ide, jds,jde, kds,kde 
  integer :: ims,ime, jms,jme, kms,kme 
  integer :: its,ite, jts,jte, kts,kte 
  integer :: mcphys_type , ilwrtyp, iswrtyp
  real ::  time, dtlt
  real ::   g      
  real ::  cp      
  real ::  cpv     
  real ::  r_d     
  real ::  r_v     
  real ::  svpt0   
  real ::  ep_1    
  real ::  ep_2    
  real ::  epsilon 
  real ::  xls     
  real ::  xlv     
  real ::  xlf     
  real ::  rhoair0 
  real ::  rhowater, rhosnow  
  real ::  cliq    
  real ::  cice    
  real ::  psat     
  integer ,parameter :: nzpmax=100
  integer :: l_unit = 4, i,j,k,kr
  character(len=5) :: ctime,cmynum
  character(len=2) :: cmcphys_type
  character(len=180) :: filein, fileout
  integer :: rec_size, irec, nz,nlev,nz1,nz2,int_byte_size 
  logical :: first_read = .true.

  real :: tempk,rliq,rice,til,qhydm,tairstr,real_byte_size
  integer :: do_radar_ref = 0 
  integer :: has_reqc, has_reqi, has_reqs 
  real :: dx, dy              ! grid spacing (m)
  real :: dt                  ! model timestep (s)
  integer  :: ke_diag            
  logical  :: wetscav_on = .false.     
  logical :: start_of_simulation =.true. 
  logical :: diagflag=.false. 
  character(len=180) :: prefix
  
  !- this for the namelist wsm.inp
  namelist /run/ prefix, time, mcphys_type, mynum

  real, allocatable, dimension(:) :: dzt, zt

  real, allocatable, dimension(:,:) :: &
       rtgt   &
       ,accpr  &
       ,pcprr  &    
       ,accph  &
       ,pcprh  &
       ,accpg  &
       ,pcprg  &
       ,accps  &
       ,pcprs  &
       ,glon   &
       ,glat

  real, allocatable, dimension(:,:,:) :: &
       thp    &
       ,theta  &
       ,pp     & 
       ,rtp    &
       ,rv     &
       ,wp     &
       ,dn0    &
       ,pi0    &
       ,rcp    &
       ,rrp    &
       ,rhp    & 
       ,rgp    &
       ,rpp    &
       ,rsp    !&
       ! ,rei    &
       ! ,rel
  !--- local vars
  real :: &
       rtgt_1d     !&
       ! ,accpr_1d    &
       ! ,pcprr_1d    &    
       ! ,accph_1d    &
       ! ,pcprh_1d    &
       ! ,accpg_1d    &
       ! ,pcprg_1d    &
       ! ,accps_1d    &
       ! ,pcprs_1d  

!   real, dimension(nzpmax) :: & !!Comentado ALF 26/10
      !  theta_1d&
      !  ,pp_1d   & 
      !  ,wp_1d   &
      !  ,dn0_1d  &
      !  ,pi0_1d  &
      !  ,thp_1d  &
      !  ,rtp_1d  &
      !  ,rv_1d   &
      !  ,rcp_1d  &
      !  ,rrp_1d  &
      !  ,rhp_1d  & 
      !  ,rgp_1d  &
      !  ,rpp_1d  &
      !  ,rsp_1d  &
      !  ,rel_1d  &
      !  ,rei_1d             

  real, allocatable, dimension(:,:,:) :: &
       th       &
       ,dz8w     &
       ,pi_phy   &
       ,p        &
       ,air_dens &
       ! ,w 		  &
       ,qv_curr,qc_curr,qr_curr,qi_curr,qs_curr,qg_curr &
       ,qh_curr,re_cloud, re_ice, re_snow, orho	      &
       ,hgt, refl_10cm , rainprod,evapprod

  real, allocatable, dimension(:,:) :: &
       rainnc	    &
       ,rainncv	 &
       ,snownc	    &
       ,snowncv	 &
       ! ,graupelnc  &
       ! ,graupelncv &
       ,hail	    & 
       ! ,hailncv	 &
       ,sr

  !- 
  integer :: nvar_out
  real :: t1, t2
  real :: wp_sum, thp_sum, theta_sum, rtp_sum, rv_sum, rrp_sum, rcp_sum
  real :: accpr_sum, pcprr_sum
  
  !! $acc routine (wsm3)
!!$acc routine (wsm5)
  !- read namelist  
  open(15,file='wsm.inp',status='old',form='formatted')    
  read(15,nml=run)
  close(15)

  !time=18000
  !mynum=365
  !mcphys_type=30
  !------------------------ reading input data
  write(ctime ,fmt="(I5.5)") int(time)
  write(cmynum,fmt="(I5.5)") mynum
  write(cmcphys_type,fmt="(I2.2)") mcphys_type

  filein= "ref/"//"WSM"//trim(cmcphys_type)//"_dataIn-"//ctime//"-"//trim(cmynum)//".bin"

  print*,"opening file for reading: ",trim(filein)

  open(newunit = l_unit,file =trim(filein),ACCESS = "stream",action="read", status="old")

  read(l_unit) m1,m2,m3, mynum
  print*,'dim1=', m1,m2,m3, mynum
  read(l_unit) ia,iz,ja,jz

  read(l_unit) ids,ide, jds,jde, kds,kde 
  read(l_unit) ims,ime, jms,jme, kms,kme
  print *, "DEBUG-ALF: ims, ime, jms, jme=", ims, ime, jms, jme
  read(l_unit) its,ite, jts,jte, kts,kte 

  if(first_read) then
     first_read = .false.
     allocate(dzt   (m1))       ;dzt    = 0.0
     allocate(zt    (m1))       ;zt     = 0.0
     allocate(rtgt  (m2,m3))    ;rtgt   = 0.0 
     allocate(accpr (m2,m3))    ;accpr  = 0.0
     allocate(pcprr (m2,m3))    ;pcprr  = 0.0    
     allocate(accph (m2,m3))    ;accph  = 0.0 
     allocate(pcprh (m2,m3))    ;pcprh  = 0.0 
     allocate(accpg (m2,m3))    ;accpg  = 0.0 
     allocate(pcprg (m2,m3))    ;pcprg  = 0.0 
     allocate(accps (m2,m3))    ;accps  = 0.0 
     allocate(pcprs (m2,m3))    ;pcprs  = 0.0
     allocate(glon  (m2,m3))    ;glon   = 0.0
     allocate(glat  (m2,m3))    ;glat   = 0.0


     allocate(thp  (m1,m2,m3))   ;thp   = 0.0 
     allocate(theta(m1,m2,m3))   ;theta = 0.0
     allocate(pp   (m1,m2,m3))   ;pp    = 0.0 
     allocate(rtp  (m1,m2,m3))   ;rtp   = 0.0 
     allocate(rv   (m1,m2,m3))   ;rv    = 0.0 
     allocate(wp   (m1,m2,m3))   ;wp    = 0.0 
     allocate(dn0  (m1,m2,m3))   ;dn0   = 0.0 
     allocate(pi0  (m1,m2,m3))   ;pi0   = 0.0 
     allocate(rcp  (m1,m2,m3))   ;rcp   = 0.0 
     allocate(rrp  (m1,m2,m3))   ;rrp   = 0.0 
     allocate(rhp  (m1,m2,m3))   ;rhp   = 0.0 
     allocate(rgp  (m1,m2,m3))   ;rgp   = 0.0 
     allocate(rpp  (m1,m2,m3))   ;rpp   = 0.0 
     allocate(rsp  (m1,m2,m3))   ;rsp   = 0.0 
     ! allocate(rei  (m1,m2,m3))   ;rei   = 0.0 
     ! allocate(rel  (m1,m2,m3))   ;rel   = 0.0 

     !---------  Local vars
     allocate(th        ( ims:ime, kms:kme, jms:jme )) ;th       = 0.0 
     allocate(dz8w      ( ims:ime, kms:kme, jms:jme )) ;dz8w     = 0.0 
     allocate(pi_phy    ( ims:ime, kms:kme, jms:jme )) ;pi_phy   = 0.0 
     allocate(p         ( ims:ime, kms:kme, jms:jme )) ;p        = 0.0 
     allocate(air_dens  ( ims:ime, kms:kme, jms:jme )) ;air_dens = 0.0 
     ! allocate(w         ( ims:ime, kms:kme, jms:jme )) ;w        = 0.0 
     allocate(qv_curr   ( ims:ime, kms:kme, jms:jme )) ;qv_curr  = 0.0 
     allocate(qc_curr   ( ims:ime, kms:kme, jms:jme )) ;qc_curr  = 0.0 
     allocate(qr_curr   ( ims:ime, kms:kme, jms:jme )) ;qr_curr  = 0.0
     allocate(qi_curr   ( ims:ime, kms:kme, jms:jme )) ;qi_curr  = 0.0
     allocate(qs_curr   ( ims:ime, kms:kme, jms:jme )) ;qs_curr  = 0.0
     allocate(qg_curr   ( ims:ime, kms:kme, jms:jme )) ;qg_curr  = 0.0
     allocate(qh_curr   ( ims:ime, kms:kme, jms:jme )) ;qh_curr  = 0.0
     allocate(re_cloud  ( ims:ime, kms:kme, jms:jme )) ;re_cloud = 0.0
     allocate(re_ice    ( ims:ime, kms:kme, jms:jme )) ;re_ice   = 0.0
     allocate(re_snow   ( ims:ime, kms:kme, jms:jme )) ;re_snow  = 0.0
     allocate(orho      ( ims:ime, kms:kme, jms:jme )) ;orho     = 0.0   
     allocate(hgt       ( ims:ime, kms:kme, jms:jme )) ;hgt      = 0.0   
     allocate(refl_10cm ( ims:ime, kms:kme, jms:jme )) ;refl_10cm= 0.0   
     allocate(rainprod  ( ims:ime, kms:kme, jms:jme )) ;rainprod = 0.0   
     allocate(evapprod  ( ims:ime, kms:kme, jms:jme )) ;evapprod = 0.0   

     allocate(rainnc    ( ims:ime,jms:jme )) ;rainnc       = 0.0
     allocate(rainncv   ( ims:ime,jms:jme )) ;rainncv      = 0.0
     allocate(snownc    ( ims:ime,jms:jme )) ;snownc       = 0.0
     allocate(snowncv   ( ims:ime,jms:jme )) ;snowncv      = 0.0
     ! allocate(graupelnc ( ims:ime,jms:jme )) ;graupelnc    = 0.0
     ! allocate(graupelncv( ims:ime,jms:jme )) ;graupelncv   = 0.0
     ! allocate(hail      ( ims:ime,jms:jme )) ;hail         = 0.0
     ! allocate(hailncv   ( ims:ime,jms:jme )) ;hailncv      = 0.0
     allocate(sr        ( ims:ime,jms:jme )) ;sr           = 0.0

  endif

  read(l_unit) mcphys_type , ilwrtyp, iswrtyp
  print*,'mcyphys=', mcphys_type
  read(l_unit) time, dtlt
  read(l_unit)  g      
  read(l_unit) cp      
  read(l_unit) cpv     
  read(l_unit) r_d     
  read(l_unit) r_v     
  read(l_unit) svpt0   
  read(l_unit) ep_1    
  read(l_unit) ep_2    
  read(l_unit) epsilon 
  read(l_unit) xls     
  read(l_unit) xlv     
  read(l_unit) xlf     
  read(l_unit) rhoair0 
  read(l_unit) rhowater  
  read(l_unit) rhosnow      
  read(l_unit) cliq    
  read(l_unit) cice    
  read(l_unit) psat
  read(l_unit) dzt
  read(l_unit) zt
  read(l_unit) glon
  read(l_unit) glat

  !-- 1st section all microphysics
  read(l_unit) thp
  read(l_unit) theta
  read(l_unit) pp   
  read(l_unit) rtp  
  read(l_unit) rv   
  read(l_unit) wp   
  read(l_unit) dn0  
  read(l_unit) pi0  
  read(l_unit) rcp    
  read(l_unit) rrp    


  read(l_unit) rtgt 
  read(l_unit) accpr
  read(l_unit) pcprr

  if(mcphys_type == 7 )  then 
     read(l_unit)  rhp  
     read(l_unit)  accph
     read(l_unit)  pcprh
  endif

  if(mcphys_type == 6 .or. mcphys_type == 7) then
     read(l_unit)  rgp  
     read(l_unit)  accpg
     read(l_unit)  pcprg
  endif

  if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then
     read(l_unit)  rpp  
     read(l_unit)  rsp  
     read(l_unit)  accps
     read(l_unit)  pcprs
  endif
  close(l_unit)


  !------------------------------ process cloud microphysics -----------------------
  IF(start_of_simulation) THEN !.or.restart.)   
     ! IF(mcphys_type == 30 ) &   
          ! CALL wsm3init(rhoair0,rhowater,rhosnow,cliq,cpv,   .false. )
     IF(mcphys_type ==  5 ) &   
          CALL wsm5init(rhoair0,rhowater,rhosnow,cliq,cpv,   .false. )
     ! IF(mcphys_type ==  6 ) &   
          ! CALL wsm6init(rhoair0,rhowater,rhosnow,cliq,cpv, 0,.false. )
     ! IF(mcphys_type ==  7 ) &   
          ! CALL wsm7init(rhoair0,rhowater,rhosnow,cliq,cpv,   .false. )

     start_of_simulation =.false.
  ENDIF

  print*,'ja, jz',ja,jz
  print*,'ia, iz',ia,iz
  print*,'its,ite, jts,jte, kts,kte',its,ite, jts,jte, kts,kte

  !--- coupling with the WSM MPs
  ! flags to calculate effec radius
  !   IF( (ilwrtyp==6 .or. iswrtyp==6)) then           
  !      has_reqc= 1 ; has_reqi= 1 ; has_reqs= 1 
  !   ELSE
  has_reqc= 0 ; has_reqi= 0 ; has_reqs= 0 
  !   ENDIF
  call cpu_time(t1)
  !$acc parallel loop private(k, kr, qv_curr, qc_curr, qr_curr, &
  !$acc qi_curr, qs_curr, qg_curr, qh_curr, pi_phy, P, dz8w, &
  !$acc dt, SR, refl_10cm, tempk, til, rliq, rice, qhydm, tairstr, TH, air_dens, &
  !$acc RAINNC, RAINNCV, SNOWNC, SNOWNCV, &
  !$acc re_cloud, re_ice, re_snow) !&
  !!!$acc collapse(2)
  do j = ja,jz
     do i = ia,iz

      !   ! **Local copies per each column** COMENTADA 26/10
      !   !- column quantities
      !   thp_1d   (1:m1)= thp  (1:m1,i,j)
      !   theta_1d (1:m1)= theta(1:m1,i,j)
      !   pp_1d    (1:m1)= pp   (1:m1,i,j)
      !   rtp_1d   (1:m1)= rtp  (1:m1,i,j)
      !   rv_1d    (1:m1)= rv   (1:m1,i,j)
      !   wp_1d    (1:m1)= wp   (1:m1,i,j)
      !   dn0_1d   (1:m1)= dn0  (1:m1,i,j)
      !   pi0_1d   (1:m1)= pi0  (1:m1,i,j)

      !   !--- mass mixing ratio
      !   rcp_1d   (1:m1)= rcp    (1:m1,i,j)
      !   rrp_1d   (1:m1)= rrp    (1:m1,i,j)

        !- surface quantities
        !rtgt_1d  = rtgt (i,j)
        ! accpr_1d = accpr(i,j) !21/11/2023
        ! pcprr_1d = pcprr(i,j) !21/11/2023

        !if(mcphys_type == 7 )  then 
        !   rhp_1d  (1:m1)= rhp  (1:m1,i,j)
           !accph_1d      = accph(i,j)
           !pcprh_1d      = pcprh(i,j)
        !else
           !rhp_1d  (1:m1)= 0.
           rhp(1:m1,i,j) = 0
           ! accph_1d      = 0. !Descomentar 21/11/2023
           ! pcprh_1d      = 0. !Descomentar 21/11/2023
        !endif

        !if(mcphys_type == 6 .or. mcphys_type == 7) then
        !   rgp_1d  (1:m1)= rgp  (1:m1,i,j)
        !   accpg_1d      = accpg(i,j)
        !   pcprg_1d      = pcprg(i,j)
        !else
           !rgp_1d  (1:m1)= 0.
           rgp(1:m1,i,j) = 0.
           ! accpg_1d      = 0. !21/11/2023
           ! pcprg_1d      = 0. !21/11/2023
        !endif

        ! if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then !21/11/2023
         !   rpp_1d  (1:m1)= rpp  (1:m1,i,j) !Comentado 26/10
         !   rsp_1d  (1:m1)= rsp  (1:m1,i,j) !Comentado 26/10
           ! accps_1d      = accps(i,j) 21/11/2023
           ! pcprs_1d      = pcprs(i,j) 21/11/2023     
      !   else
      !      rpp_1d  (1:m1)= 0.
      !      rsp_1d  (1:m1)= 0.
      !      accps_1d      = 0.
      !      pcprs_1d      = 0.  
        ! endif !21/11/2023

      !MSUDO - 23/11
      !   !--- coupling with the WSM MPs
      !   ! flags to calculate effec radius
      ! !   IF( (ilwrtyp==6 .or. iswrtyp==6)) then           
      ! !      has_reqc= 1 ; has_reqi= 1 ; has_reqs= 1 
      ! !   ELSE
      !      has_reqc= 0 ; has_reqi= 0 ; has_reqs= 0 
      ! !   ENDIF


        ! ** initial definitions **
        dt= dtlt        ! time step            (s)
        !rainprod  =0.0  ! for scaveging   aerosols/gases !MSUDO - nao usado em WSM5
        !evapprod  =0.0  ! for evaporation aerosols/gases !MSUDO - nao usado em WSM5
        SR        =0.0  ! fraction of snow of the total water
        ! ( for land surface models)
        refl_10cm =0.0  ! 
        !ke_diag   = kte !MSUDO - nao usado em WSM5

        ! **apparently another copies** ??
        !- surface precipitation (total accumulated)
        !!!!$acc serial
        RAINNC    (1,1)=  accpr(i,j) !accpr_1d 21/11/2023 !- rain+ice+snow+graupel+hail
        SNOWNC    (1,1)=  accps(i,j) !accps_1d 21/11/2023!- ice+snow
        ! GRAUPELNC (1,1)=  accpg_1d !- graupel 21/11/2023
        ! HAIL      (1,1)=  accph_1d !- hail 21/11/2023
        !!!!$acc end serial
        ! **First Loop over a vertical column**
        do k = 1,kme-1

           kr = k + 1
           qv_curr (1,k,1)= max(1.e-12,rtp(kr,i,j) - &    ! QV
                (rcp(kr,i,j)+rrp(kr,i,j)+rpp(kr,i,j)+rsp(kr,i,j)+rgp(kr,i,j)+rhp(kr,i,j)))   
           qc_curr (1,k,1)= max(0.0,rcp(kr,i,j))       ! QC     
           qr_curr (1,k,1)= max(0.0,rrp(kr,i,j))       ! QR   
           qi_curr (1,k,1)= max(0.0,rpp(kr,i,j))       ! QI   
           qs_curr (1,k,1)= max(0.0,rsp(kr,i,j))       ! QS   
           qg_curr (1,k,1)= max(0.0,rgp(kr,i,j))       ! QG

           qh_curr (1,k,1)= max(0.0,rhp(kr,i,j))       ! QH   

           pi_phy  (1,k,1)= (pp(kr,i,j)+pi0(kr,i,j))*cpi ! Exner function/cp (dimensionless)

           P   (1,k,1)= ( (pp(kr,i,j)+pi0(kr,i,j))*cpi )** cpor * p00      ! pressure(Pa)
           ! W   (1,k,1)= wp(kr,i,j)    ! vertical velocity (m/s) ! must be at center or face? ASK !21/11/2023

           dz8w   (1,k,1)= rtgt_1d/dzt(kr) ! layer thickness (m) 
           !print*,'dz8',k,dz8w   (1,k,1)
        enddo 

        ! **Second Loop over a vertical column**
        !- get potential temperature (theta) from theta_il (thp) and condensates
        DO k=1,kme -1 
           kr = k + 1
           tempK = theta(kr,i,j)* (pp(kr,i,j)+pi0(kr,i,j))*cpi 
           til   = thp(kr,i,j)* (pp(kr,i,j)+pi0(kr,i,j))*cpi 

           rliq  =  qc_curr(1,k,1) + qr_curr(1,k,1)                
           rice  =  qi_curr(1,k,1) + qs_curr(1,k,1) + qg_curr(1,k,1) + qh_curr(1,k,1)
           qhydm =  alvl * rliq + alvi * rice

           if (tempK > 253.) then
              tairstr = 0.5 * (til + sqrt(til * (til + cpi4 * qhydm)))
           else
              tairstr = til * (1. + qhydm * cp253i)
           endif
           !- updated potential temperature TH in Kelvin (adv+dif+rad+conv+)
           TH (1,k,1) = tairstr / pi_phy(1,k,1)

           !- air density
           air_dens(1,k,1) = P(1,k,1)/(287.04*tempK*(1.+0.608*qv_curr(1,k,1)))
           !air_dens(1,k,1)= dn0(kr) 

        ENDDO !PAREI 26/10

      !   ! ** WSM5 main routine **
        IF(mcphys_type == 5)                    &   
             
             CALL wsm5(                         &
             TH,                        &! potential temperature    (K)
             qv_curr,                   &! QV=qv_curr,     
             qc_curr,                   &! QC=qc_curr,     
             qr_curr,                   &! QR=qr_curr,     
             qi_curr,                   &! QI=qi_curr,     
             qs_curr,                   &! QS=qs_curr,     
                                
             air_dens,                  &          
             pi_phy,                    &! exner function (dimensionless)
             P,                         &! pressure(Pa)
             dz8w,                      &! deltaz
                                
             dt,      &                  ! time step              (s)
             g,       &
             cp,      &
             cpv,     &		! COMENTAR - marcar o q é privado
             r_d,     &
             r_v,     &     
             svpt0,   &
             ep_1,    &
             ep_2,    &
             epsilon, &
             xls,     &
             xlv,     &
             xlf,     &
             rhoair0, &
             rhowater,&  
             cliq,    &
             cice,    &
             psat,    &     
             RAINNC,                    &
             RAINNCV,                   &
             SNOWNC,                    &
             SNOWNCV,                   &
             SR,                        &
                                
             refl_10cm,                 &
             diagflag,                  &
             do_radar_ref,              &
                                
             has_reqc,                  & 
             has_reqi,                  &  
             has_reqs,                  & 
                                
             re_cloud,                  & 
             re_ice,                    &
             re_snow,                   &
             IDS,IDE, JDS,JDE, KDS,KDE, &
             IMS,IME, JMS,JME, KMS,KME, &
             ITS,ITE, JTS,JTE, KTS,KTE  &
             )

      !   ! ** Third loop over a colunm K, returning values **
      !   !- updated variables after microphysics processes (from WSM to BRAMS)
        DO k=1,kme-1
           kr=k+1
           rtp(kr,i,j)= qv_curr(1,k,1) + &
                qc_curr(1,k,1) + &     
                qr_curr(1,k,1) + &    
                qi_curr(1,k,1) + &    
                qs_curr(1,k,1) + &    
                qg_curr(1,k,1) + &    
                qh_curr(1,k,1)  

           rcp(kr,i,j)= qc_curr(1,k,1)
           rrp(kr,i,j)= qr_curr(1,k,1)
           rpp(kr,i,j)= qi_curr(1,k,1)
           rsp(kr,i,j)= qs_curr(1,k,1)
           rgp(kr,i,j)= qg_curr(1,k,1)
           rhp(kr,i,j)= qh_curr(1,k,1)  

           rv(kr,i,j)= max(1.0e-12, rtp(kr,i,j) -(rcp(kr,i,j)+rrp(kr,i,j)+rpp(kr,i,j)+rsp(kr,i,j)+rgp(kr,i,j)+rhp(kr,i,j)))

           theta(kr,i,j) =  TH(1,k,1)
           tempK        =  TH(1,k,1)*pi_phy(1,k,1)

           rliq     =  qc_curr(1,k,1) + qr_curr(1,k,1)       
           rice     =  qi_curr(1,k,1) + qs_curr(1,k,1) + qg_curr(1,k,1) + qh_curr(1,k,1)

           !- update liq-ice potential temperature THP in Kelvin including microphysics processes
           thp(kr,i,j)  =   TH(1,k,1)*(1. + alvl * rliq/(cp * max(tempK,253.))  &
                + alvi * rice/(cp * max(tempK,253.)) ) **(-1.0)      
        ENDDO
        !- definition for k=1
        rtp(1,i,j)  = rtp(2,i,j)
        rcp(1,i,j)  = rcp(2,i,j)
        rrp(1,i,j)  = rrp(2,i,j)
        rpp(1,i,j)  = rpp(2,i,j)
        rsp(1,i,j)  = rsp(2,i,j)
        rgp(1,i,j)  = rgp(2,i,j)
        rhp(1,i,j)  = rhp(2,i,j)
        rv(1,i,j)  = rv(2,i,j)
        thp(1,i,j)  = thp(2,i,j)
        theta(1,i,j) = theta(2,i,j)

        IF( (ilwrtyp==6 .or. iswrtyp==6)) then           
           DO k=1,kme-1
              kr=k+1
              ! rel (kr,i,j) = re_cloud (1,k,1) * 1.e+6 ! RRTM requires in micrometer !22/11/2023
              ! rei (kr,i,j) = re_ice   (1,k,1) * 1.e+6 ! RRTM requires in micrometer !22/11/2023
           ENDDO
           ! rel (1,i,j) =rel (2,i,j) !;  rel (kme) =rel (kme-1) !22/11/2023
           ! rei (1,i,j) =rei (2,i,j) !;  rei (kme) =rei (kme-1) !22/11/2023
        ENDIF
        !- surface precipitation (units are kg/m^2 = mm)
        !- RAINNC and RAINNCV constains all precipitation hidrometeors (rain, graupel, snow, ...)
        ! accpr_1d = RAINNC    (1,1) ! a = accum !21/11/2023
		accpr(i,j) = RAINNC    (1,1) ! a = accum
        ! pcprr_1d = RAINNCV   (1,1) ! p = for each dt  (or per time step)
		pcprr(i,j) = RAINNCV   (1,1) ! p = for each dt  (or per time step)
        ! accps_1d = SNOWNC    (1,1) !21/11/2023
		accps(i,j) = SNOWNC    (1,1)
        pcprs(i,j) = SNOWNCV   (1,1) 
        ! accpg_1d = GRAUPELNC (1,1) 21/11/2023
        ! pcprg_1d = GRAUPELNCV(1,1) 21/11/2023
        ! accph_1d = HAIL      (1,1) 21/11/2023 
        ! pcprh_1d = HAILNCV   (1,1) 

        !!- column quantities - NOT NECESSARY (Marcelo 30/10)
        ! thp  (1:m1,i,j) =thp_1d  (1:m1) 
        ! theta(1:m1,i,j) =theta_1d(1:m1)
        ! rtp  (1:m1,i,j) =rtp_1d  (1:m1)
        ! rv   (1:m1,i,j) =rv_1d   (1:m1)  
        ! rcp     (1:m1,i,j) =rcp_1d  (1:m1)   
        ! rrp     (1:m1,i,j) =rrp_1d  (1:m1)
        ! if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then   
           ! rpp    (1:m1,i,j) =rpp_1d (1:m1)   
           ! rsp    (1:m1,i,j) =rsp_1d (1:m1)   
        ! endif
        ! if(mcphys_type == 6 .or. mcphys_type == 7) rgp(1:m1,i,j) =rgp_1d  (1:m1)  
        ! if(mcphys_type == 7           )            rhp(1:m1,i,j) =rhp_1d  (1:m1)
        ! if( ilwrtyp==6 .or. iswrtyp==6 ) then
           ! rei   (1:m1,i,j) =rei_1d  (1:m1)
           ! rel   (1:m1,i,j) =rel_1d  (1:m1)
        ! endif
		!! END - NOT NECESSARY (Marcelo 30/10)

        !- surface quantities
        ! accpr(i,j) = accpr_1d !constains all precipitation hidrometeors (rain, graupel, snow, ...) !21/11/2023
        ! pcprr(i,j) = pcprr_1d !constains all precipitation hidrometeors (rain, graupel, snow, ...) !21/11/2023
        ! if(mcphys_type == 5 .or. mcphys_type == 6 .or. mcphys_type == 7 ) then !21/11/2023
           ! accps(i,j) = accps_1d !21/11/2023
           ! pcprs(i,j) = pcprs_1d !21/11/2023
        ! endif !21/11/2023
        ! if(mcphys_type == 6 .or. mcphys_type == 7) then
           ! accpg(i,j) = accpg_1d
           ! pcprg(i,j) = pcprg_1d
        ! endif
        ! if(mcphys_type == 7) then
           ! accph(i,j) = accph_1d
           ! pcprh(i,j) = pcprh_1d
        ! endif
     enddo; enddo ! loop i,j
     !$acc end parallel loop
!!!!     !$acc end data
	 call cpu_time(t2)
	 print *,"times DO-DO",t2,t1,t2-t1
     !=======================================================================
     inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list
     rec_size = m2*m3*int_byte_size    

     nz1 = 2
     nz2 = m1-1
	 
	!MÉTRICA - medição de acurácia
	wp_sum = 0.0
	thp_sum = 0.0
	theta_sum = 0.0
	rtp_sum = 0.0
	rv_sum = 0.0
	rrp_sum = 0.0
	rcp_sum = 0.0
	do nz=nz1,nz2
	  do i=1,m2
	    do j=1,m3
			wp_sum = wp_sum + wp(nz,i,j)
			thp_sum = thp_sum + thp(nz,i,j)
			theta_sum = theta_sum + theta(nz,i,j)
			rtp_sum = rtp_sum + 1000.*rtp(nz,i,j)
			rv_sum = rv_sum + 1000.*rv(nz,i,j)
			rrp_sum = rrp_sum + 1000.*rrp(nz,i,j)
			rcp_sum = rcp_sum + 1000.*rcp(nz,i,j)
		enddo
	  enddo
	enddo
	print*,"wp sum avg", wp_sum, wp_sum/((nz2-nz1+1)*m2*m3)
	print*,"thp sum avg", thp_sum, thp_sum/((nz2-nz1+1)*m2*m3)
	print*,"theta sum avg", theta_sum, theta_sum/((nz2-nz1+1)*m2*m3)
	print*,"rtp sum avg", rtp_sum, rtp_sum/((nz2-nz1+1)*m2*m3)
	print*,"rv sum avg", rv_sum, rv_sum/((nz2-nz1+1)*m2*m3)
	print*,"rrp sum avg", rrp_sum, rrp_sum/((nz2-nz1+1)*m2*m3)
	print*,"rcp sum avg", rcp_sum, rcp_sum/((nz2-nz1+1)*m2*m3)

   accpr_sum = 0.0
   pcprr_sum = 0.0
   do i=1,m2
      do j=1,m3
         accpr_sum = accpr_sum + 3600*accpr(i,j)
         pcprr_sum = pcprr_sum + 3600*pcprr(i,J)
      enddo
   enddo
   print*,"accpr sum avg", accpr_sum, accpr_sum/(m2*m3)
   print*,"pcprr sum avg", pcprr_sum, pcprr_sum/(m2*m3)
	!FIM MÉTRICA
	
     fileout=trim(prefix)//"_WSM"//trim(cmcphys_type)//"_dataOut-"//ctime//"-"//trim(cmynum)
     print*," ===> writing output in file: ",trim(fileout)


     open(newunit = l_unit, file=trim(fileout)//".gra", &
          form='unformatted', access='direct', status='replace', recl=rec_size)
     irec=1
     nlev= m1-1
     nvar_out = 0
     !-- 3d updated fields
     do nz=nz1,nz2;  write(l_unit,rec=irec) wp         (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
     do nz=nz1,nz2;  write(l_unit,rec=irec) thp        (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
     do nz=nz1,nz2;  write(l_unit,rec=irec) theta      (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
     do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*rtp  (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
     do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*rv   (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
     do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*rrp  (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1
     do nz=nz1,nz2;  write(l_unit,rec=irec) 1000.*rcp  (nz,1:m2,1:m3); irec=irec+1 ;enddo; nvar_out = nvar_out + 1

     !-- 2d 
     write(l_unit,rec=irec) 3600*accpr(1:m2,1:m3); irec=irec+1 ; nvar_out = nvar_out + 1
     write(l_unit,rec=irec) 3600*pcprr(1:m2,1:m3); irec=irec+1 ; nvar_out = nvar_out + 1

     close(l_unit)  
     !== number of levels 2: m1 -1 
                          nlev = nz2-nz1+1

                          open(newunit = l_unit, file=trim(fileout)//".ctl", action='write', status='replace')

                          fileout=trim(fileout)//".gra"
                          write(l_unit,2001) trim(fileout) 
                          !writing others infos to ctl
                          write(l_unit,*) 'undef -0.9990000E+34'
                          write(l_unit,*) 'title WSM '
                          write(l_unit,*) 'xdef ',m2,' linear ',glon(1,1),glon(2,1)-glon(1,1)
                          write(l_unit,*) 'ydef ',m3,' linear ',glat(1,1),glat(1,2)-glat(1,1)
                          write(l_unit,2005) nlev, zt(nz1:nz2)
                          write(l_unit,*) 'tdef 1 linear 00:00Z01JAN200 1mo'
                          write(l_unit,*) 'vars ' , nvar_out
                          write(l_unit,*) 'wp'    ,nlev,'99 ','K'
                          write(l_unit,*) 'thp'   ,nlev,'99 ','K'
                          write(l_unit,*) 'theta' ,nlev,'99 ','K'
                          write(l_unit,*) 'rtp'   ,nlev,'99 ','g/kg'
                          write(l_unit,*) 'rv'    ,nlev,'99 ','g/kg'
                          write(l_unit,*) 'rrp'   ,nlev,'99 ','g/kg'
                          write(l_unit,*) 'rcp'   ,nlev,'99 ','g/kg'
                          write(l_unit,*) 'apr',' 0 ','99 ','K'
                          write(l_unit,*) 'ppr',' 0 ','99 ','K'

                          write(l_unit,*) 'endvars'

                          close(l_unit)

2005                      format('zdef ',i4,' levels ',100f8.1)
2001                      format('dset ',A100)
                        end program wsm_test