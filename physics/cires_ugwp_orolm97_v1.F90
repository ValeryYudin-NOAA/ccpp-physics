module cires_ugwp_orolm97_v1
!
! rename  cires_ugwp_oro_v1
!
contains

   subroutine orogw_v1 (im, km,  imx,   me,  master, dtp, kdt, do_tofd,   &
                      grav, con_omega, rd, cpd, rv, pi, arad, fv,         &
		      xlatd, sinlat, coslat, sparea,                      &
                      cdmbgwd, hprime, oc, oa4, clx4, theta, sigmad,      &
		      gammad, elvmaxd, sgh30,  kpbl,                      &					    	     
                      u1 ,v1, t1, q1, prsi,del,prsl,prslk, zmeti, zmet,   &		      	      
                      pdvdt, pdudt, pdtdt, pkdis, dusfc, dvsfc,rdxzb  ,   &		      
                      zobl, zlwb, zogw, tau_ogw, dudt_ogw, dvdt_ogw,      &
                      dudt_obl, dvdt_obl,dudt_ofd, dvdt_ofd,              &
                      du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol,         &
                      du_ofdcol, dv_ofdcol)
!64
!
!     subroutine orogw_v1(im,  km,  imx, do_tofd,                          &   
!         pdvdt, pdudt, pdtdt, pkdis, u1,v1,t1,q1,kpbl,                    &  
!         prsi,del,prsl,prslk, zmeti, zmet, dtp, kdt, hprime,              &
!         oc, oa4, clx4, theta, sigmad, gammad, elvmaxd,                   &
!         grav, con_omega, rd, cpd, rv, pi, arad, fv, sgh30,               &
!         dusfc, dvsfc,  xlatd, sinlat, coslat, sparea,                    &
!         cdmbgwd, me, master, rdxzb,  zobl, zlwb, zogw, tau_ogw,          &	 
!         dudt_ogw, dvdt_ogw,dudt_obl,dvdt_obl,dudt_ofd, dvdt_ofd,         &
!        du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol, du_ofdcol dv_ofdcol)
	 
!---------------------------------------------------------------------------
! ugwp_v1: orogw_v1 following recent updates of Lott & Miller 1997
!          eventually will be replaced with more "advanced" LLWB
!          and multi-wave solver that produce competitive FV3GFS-skills  
! 
!          computation of kref for ogw + coorde diagnostics
!          all constants/parameters inside cires_ugwp_initialize.f90 
!
! 10/2020   main updates 
!         (a)   introduce extra diagnostics of x-y obl-ofd-ogw as in the GSL-drag
!                for intercomparisons
!
!         (b)  quit with cdmbgwd(1:2)
!              cdmbgwd(1) = 1 for all resolutions, number of hills control SA-effects
!              cdmbgwd(2) = 1 ...............number of hills control SA-effects
!
!          (c) cleff = pi2/(nk*dx)    lheff = nk*dx  (nk = 6,4,2, 1)  
!                alternative       lheff = min( dogw=hprime/sigma*gamma, dx)  
!                we still not use  the "broad spectral solver"    
!
!          (d)  hefff = (nsig * hprime -znlk)/nsig, orchestrating MB and OGW
!
!          (e)  for linsat-solver "eddy" damping Ked = Ked * Nhills, scale-aware
!               amplification of the momentum deposition for low-res simulations
!----------------------------------------    

      use machine ,      only : kind_phys
      use ugwp_common,   only : dw2min, velmin
 
      use ugwp_oro_init, only : rimin,  ric,     efmin,     efmax   ,  &
                                hpmax,  hpmin,   sigfaci => sigfac  ,  &
                                dpmin,  minwnd,  hminmt,    hncrit  ,  &
                                rlolev, gmax,    veleps,    factop  ,  &
                                frc,    ce,      ceofrc,    frmax, cg, &
                                fdir,   mdir,    nwdir,                &
                                cdmb,   cleff,   fcrit_gfs, fcrit_mtb, &
                                n_tofd, ze_tofd, ztop_tofd

      use cires_ugwp_module, only : kxw,  max_kdis, max_axyz

      use cires_orowam2017, only     : oro_wam_2017

!      use cires_vert_orodis_v1, only : ugwp_tofd1d   ! now inside the module "cires_ugwp_orolm97_v1"
      
!----------------------------------------
      implicit none
        
      real(kind=kind_phys), parameter :: sigfac = 3     ! N*hprime height of Subgrid Hill over which SSO-flow is specified
      real(kind=kind_phys), parameter :: sigfacs = 0.25 ! M*hprime height is the low boundary of the hill 
          
      character(len=8)                :: strsolver='pss-1986'  ! current operational Ri-solver or  'spect_2020'
      

      real(kind=kind_phys)            :: gammin = 0.00999999       ! a/b = gammma_min =1% <====>      
      real(kind=kind_phys), parameter :: nhilmax = 25.             ! max number of SSO-hills in grid-box
      real(kind=kind_phys), parameter :: sso_min = 1500.           ! min-lenghth of the hill
      
      logical, parameter              :: do_adjoro = .false.       ! 
!----------------------------------------      
      
      integer, intent(in) :: im, km, imx, kdt
      integer, intent(in) :: me, master
      logical, intent(in) :: do_tofd
      
      integer, intent(in)              :: kpbl(im)    ! index for the pbl top layer!
      real(kind=kind_phys), intent(in) :: dtp         !  time step
      real(kind=kind_phys), intent(in) :: cdmbgwd(2)
      
      real(kind=kind_phys), intent(in) :: hprime(im), oc(im), oa4(im,4),       &
                                          clx4(im,4), theta(im),               &
                                          sigmad(im), gammad(im), elvmaxd(im)

      real(kind=kind_phys), intent(in) :: grav, con_omega, rd, cpd, rv,        &
                                          pi, arad, fv
					  
      real(kind=kind_phys), intent(in) :: sgh30(im)       
      real(kind=kind_phys), intent(in), dimension(im,km) ::                    &
                                  u1,  v1,   t1, q1,del, prsl, prslk, zmet
     
      real(kind=kind_phys), intent(in),dimension(im,km+1):: prsi, zmeti
      
      real(kind=kind_phys), intent(in) :: xlatd(im),sinlat(im), coslat(im)
      real(kind=kind_phys), intent(in) :: sparea(im)

!    
!output -phys-tend
      real(kind=kind_phys),dimension(im,km),intent(out) ::                     &
                           pdvdt,    pdudt,    pkdis, pdtdt
! output - diag-coorde
      real(kind=kind_phys),dimension(im,km),intent(out) ::                      &
                 dudt_ogw,dvdt_ogw,  dudt_obl,dvdt_obl, dudt_ofd,dvdt_ofd 
			   
      real(kind=kind_phys),dimension(im),intent(out) ::  dusfc,   dvsfc,        & 
               du_ogwcol,dv_ogwcol,  du_oblcol,dv_oblcol, du_ofdcol,dv_ofdcol                     
!                     
      real(kind=kind_phys),dimension(im),intent(out) :: rdxzb
      real(kind=kind_phys),dimension(im),intent(out) :: zobl, zogw, zlwb, tau_ogw
		    

!
!---------------------------------------------------------------------
! # of permissible sub-grid orography hills for "any" resolution  < 25
!    correction for "elliptical" hills based on shilmin-area =sgrid/25 
!     4.*gamma*b_ell*b_ell  >=  shilmin
!     give us limits on [b_ell & gamma *b_ell] > 5 km =sso_min
!     gamma_min = 1/4*shilmin/sso_min/sso_min
!23.01.2019:  cdmb = 4.*192/768_c192=1 x 0.5
!     192: cdmbgwd        = 0.5, 2.5
!     cleff = 2.5*0.5e-5 * sqrt(192./768.) => lh_eff = 1004. km
!      6*dx = 240 km 8*dx = 320. ~ 3-5 more effective OGW-lin
!---------------------------------------------------------------------
! 
! locals vars for SSO
!

      real(kind=kind_phys), dimension(im)    :: oa,  clx
      real(kind=kind_phys), dimension(im)    :: sigma, gamma, elvmax  ! corrected sigmaD, gammaD, elvmaxD
            
      real(kind=kind_phys)            :: shilmin, sgrmax, sgrmin
      real(kind=kind_phys)            :: belpmin, dsmin,  dsmax
      
      real(kind=kind_phys)            :: arhills(im), mkd05_hills(im)   ! number of hills in the grid                             
! 
! locals        mean flow  ...etc
!
      real(kind=kind_phys), dimension(im,km) :: ri_n, bnv2, ro
      real(kind=kind_phys), dimension(im,km) :: vtk, vtj, velco
!==================      
!mtb  
!==================
      real(kind=kind_phys)            :: ztoph,zlowh,ph_blk, dz_blk  
      real(kind=kind_phys), dimension(im)    :: wk, pe, ek, up
      
      real(kind=kind_phys), dimension(im,km) :: db, ang, uds

      real(kind=kind_phys) :: zlen, dbtmp, r, phiang, dbim, zr
      real(kind=kind_phys) :: eng0, eng1, cosang2, sinang2
      real(kind=kind_phys) :: bgam, cgam, gam2, rnom, rdem    
         
!==================
! tofd
!     some constants now in "use ugwp_oro_init" +   "use ugwp_common"
!
!==================

      real(kind=kind_phys)   :: unew, vnew,  zpbl,  sigflt, zsurf
      real(kind=kind_phys), dimension(km)    :: utofd1, vtofd1
      real(kind=kind_phys), dimension(km)    :: epstofd1, krf_tofd1
      real(kind=kind_phys), dimension(km)    :: up1, vp1, zpm
!==================
! ogw
!==================
      real(kind=kind_phys)            :: xlingfs
      logical                         :: icrilv(im)
!
      real(kind=kind_phys), dimension(im) :: xn, yn, ubar, vbar, ulow, &
                      roll,  bnv2bar, scor, dtfac, xlinv, delks, delks1
!
      real(kind=kind_phys) :: taup(im,km+1), taud(im,km)
      real(kind=kind_phys) :: taub(im), taulin(im), ahdxres(im)
      real(kind=kind_phys) :: heff, hsat, hdis
      integer, dimension(im) :: kref, idxzb, ipt, kreflm, iwklm, iwk, izlow
!   
! local scalars     
!check what we need
!
      real(kind=kind_phys) ::   bnv,  fr, ri_gw, brvf, fr2
      real(kind=kind_phys) ::   tem,   tem1,  tem2, temc, temv
      real(kind=kind_phys) ::   ti,     rdz,   dw2,   shr2, bvf2
      real(kind=kind_phys) ::   rdelks, efact, coefm, gfobnv
      real(kind=kind_phys) ::   scork,  rscor, hd,    fro,  sira
      real(kind=kind_phys) ::   dtaux,  dtauy, zmetp, zmetk
      
      real(kind=kind_phys) ::   grav2, rcpdt, windik, wdir
      real(kind=kind_phys) ::   sigmin, dxres,sigres,hdxres, cdmb4, mtbridge

      real(kind=kind_phys) ::   kxridge, inv_b2eff, zw1, zw2
      real(kind=kind_phys) ::   belps, aelps, nhills, selps

      real(kind=kind_phys) ::   rgrav, rcpd, rcpd2, rad_to_deg, deg_to_rad
      real(kind=kind_phys) ::   pi2, pi2h, rdi, gor, grcp, gocp, gr2, bnv2min
      
      real(kind=kind_phys) ::   cleff_max  ! resolution-aware max-wn  
      real(kind=kind_phys) ::   nonh_fact  ! non-hydroststic factor 1.-(kx/kz_hh)**2        
!       
!      
! local integers
!     
      integer ::   kmm1, kmm2, lcap, lcapp1
      integer ::   npt,   kbps, kbpsp1,kbpsm1
      integer ::   kmps,  idir, nwd,  klcap, kp1, kmpbl, kmll
      integer ::   k_mtb, k_zlow, ktrial, klevm1
      integer ::   i, j, k
      
      
!===========================
! First step Check do we have sub-grid hills
!      
!
! out-arrays are zreoed in unified_ugwp.F90
!
      do i=1,im
        rdxzb(i)    = 0.0	
        dusfc(i)    = 0.0
        dvsfc(i)    = 0.0
        ipt(i) = 0 
      enddo     
 
! ----  for lm and gwd calculation points
!cires_ugwp_initialize.F90:      real, parameter :: hpmax=2400.0, hpmin=25.0  
!cires_ugwp_initialize.F90:      real,      parameter :: hminmt=50.     ! min mtn height (*j*) 
!----  for lm and gwd calculation points  


      npt = 0
      
      do i = 1,im
        if ( elvmaxd(i) >= hminmt .and. hprime(i)  >= hpmin ) then          
          npt      = npt + 1
          ipt(npt) = i
	endif  
      enddo
      
      if (npt == 0) then
!         print *,  'oro-npt = 0 elvmax ', maxval(elvmaxd), hminmt
!         print *,  'oro-npt = 0 hprime ', maxval(hprime), hpmin	     	    
        return      ! no ogw/mbl calculation done
      endif   
!===========================
! scalars from phys-contants
!===========================     
      rcpdt = 1.0 / (cpd*dtp)
      grav2 = grav + grav
!
      rgrav = 1.0/grav
      rcpd = 1.0/cpd
      rcpd2 = 0.5/cpd
      rad_to_deg=180.0/pi
      deg_to_rad=pi/180.0
      pi2 = 2.*pi
      pi2h = 0.5*pi
      rdi = 1.0/rd
      gor = grav/rd
      grcp = grav*rcpd
      gocp = grcp
      gr2  = grav*gor
      bnv2min = (pi2/1800.)*(pi2/1800.)          ! tau_BV_max = 30 min !      
!===========================
! Start         
! 
! initialize gamma and sigma
!
      gamma(:) = gammad(:)
      sigma(:) = sigmad(:)
!

      
!=======================================================================       
! mtb-blocking  sigma_min and dxres => cires_initialize (best way ....)
!  
      sgrmax = maxval(sparea) ; sgrmin = minval(sparea)
      dsmax  = sqrt(sgrmax)   ; dsmin  = sqrt(sgrmin)
      
      cleff_max = pi2/(dsmin/5.)          ! maxval for kx = 6.28/(dx_min/5. ~2.5 km)
      
      dxres   = pi2*arad/float(imx)
      hdxres  = 0.5*dxres
!     shilmin = sgrmin/nhilmax            ! 
      gammin = min(sso_min/dxres, 1.)     ! 


      sigmin = 2.*hpmin/dxres             ! min-slope Hmin= 2*hpmin, dxres=Lmax
      kxridge = float(imx)/arad           !     Lh=~6.28*dx* cdmbgwd(2)

      if (me == master .and. kdt == 1) then
        print *, ' orogw_v1 kxridge ', kxridge
        print *, ' orogw_v1 scale2 ', cdmbgwd(2)
        print *, ' orogw_v1 imx ', imx
        print *, ' orogw_v1 gam_min ', gammin
        print *, ' orogw_v1 sso_min ', sso_min
      endif
      
!============================================================
! Purpose to adjust oro-specification on the fly
!     needs to be done 1-time during init-n for each block
!  hprime sigma gamma and grid-length must be "related"
!  width_mount_a = hprime/sigma *2 < dxres
!  width_mount_b = width_mount_a * gamma
!
!   Sellipse=  4 a*b = (width_mount_a)^2 *gamma < Sarea
!   Limiters on "elongated" hills gamma= b/a    < gam_min
!   Limiters on "longest"   hills (b,  a)       <  sqrt(area)
!   sigma = hprime/a_elliptical_hill
!           0.01 <  gamma=a_hill/b_hill < 1    
!     2*hpmin/dx <  sigma               < 1. 
!  Nhills = (dx*dy=Sarea)/(pi*  a_hill *b_hill)
!=============================================================  

           arhills(:) =0.
       mkd05_hills(:) =0.
               	  
        do j = 1,npt	  
          i = ipt(j)
          dxres = sqrt(sparea(i))
	  hdxres=0.5*dxres
	  ahdxres(i) = hdxres
!	  
! min-adjustment: 1) abs(gamma(i)) ; 2) sigres = max(sigmin, sigma(i))
!	  
	  sigres = max(sigmin, sigma(i))
	  sigma(i) =sigres
	  aelps = min( hprime(i)/sigres, hdxres)
	  belps = min(aelps/abs(gamma(i)), hdxres)
	  
        if (do_adjoro ) then   
!	
! more adjustments "lengths", gamma and sigma
!	      
          if (hprime(i) > hdxres*sigres) sigres= hprime(i)/hdxres
          aelps = min( hprime(i)/sigres, hdxres)
	  sigma(i) = sigres
          if (gamma(i) > 0.0 ) belps = min(aelps/gamma(i), hdxres)
!
! small-scale "turbulent" oro-hills < sso_min,     sso_min_dx = 1.5 km
! will be treated as "circular" elevations
!
          if( aelps < sso_min ) then             
! 
! a, b > sso_min upscale ellipse  a/b > 0.1 a>sso_min & h/b=>new_sigm
!
            aelps = sso_min 
            if (belps < sso_min ) then
              gamma(i) = 1.0
              belps = aelps*gamma(i)
            else
              gamma(i) = min(aelps/belps, 1.0)
            endif
            
	    sigma(i) = hprime(i)/aelps
            gamma(i) = min(aelps/belps, 1.0)
	    
          endif      !aelps < sso_min
         endif                                    !  ============== (do_adjoro ) 
	 
          selps      = belps*belps*gamma(i)*pi    ! area of the elliptical hill 
	  
          nhills     = min(nhilmax, sparea(i)/selps)
          arhills(i) = max(nhills, 1.0)

!333   format( ' nhil: ', i6, 4(2x, f9.3), 2(2x, e9.3))	    
!      if (kdt==1 )
!     & write(6,333) nint(nhills)+1,xlatd(i), hprime(i),aelps*1.e-3,
!     &   belps*1.e-3, sigma(i),gamma(i)
        
       enddo

!=======================================================================       
! mtb-blocking : LM-1997; Zadra et al. 2004 ;metoffice dec 2010 H Wells
!=======================================================================          

      do i=1,npt
        iwklm(i)  = 2
        idxzb(i)  = 0 
        kreflm(i) = 0
      enddo
 
      do k=1,km
        do i=1,im
          db(i,k)  = 0.0
          ang(i,k) = 0.0
          uds(i,k) = 0.0 
        enddo
      enddo

      kmm1 = km - 1 ;  kmm2   = km - 2 ; kmll   = kmm1
      lcap = km     ;  lcapp1 = lcap + 1 
      
      cdmb4 = 0.25*cdmb 
       
      do i = 1, npt
        j = ipt(i)
        elvmax(j) = min (elvmaxd(j)*0. + sigfac * hprime(j), hncrit)
        izlow(i)  = 1                                                  ! surface-level
      enddo
      
      
!===================================================================
! below iwklm-level  H= 3*hp, and izlow = 0.5*Hp or the "first" layer
!       are used tp estimate "Mean" Flow that interact with SG-HILL
! if sig*HP < Hpbl => GWs-> above PBL
! WRF:    ( 1 to max(2*Hp or H_pbl)
! GFS-15/16: OGWs (1 to max(Kpbl+1, or K_dPs=(Ps-Pk=50hPa) ~ 950 mb)
!          excitation above Kref
!    BLOCKING:  ZDOMAIN (1  - Kaver => ELVMAX(J) + sigfac * hp)
!===================================================================


      do k = 1, kmm1
        do i = 1, npt
          j = ipt(i)
	  
          ztoph   = sigfac * hprime(j)
          zlowh   = sigfacs* hprime(j) 
          zmetp   =  zmet(j,k+1) 
          zmetk   =  zmet(j,k)
!	    
! GFSv15/16: izlow=1
! elvmax(j)=elvmaxd(J) + sig*hp: if (( elvmax(j) <= zmetp) .and. (elvmax(j).ge.zmetk) ) iwklm(i)  =  max(iwklm(i), k+1 ) 
!   

          if (( ztoph <= zmetp) .and. (ztoph >= zmetk) ) iwklm(i)  =  max(iwklm(i), k+1 )
          if (zlowh <= zmetp .and. zlowh >= zmetk)       izlow(i)  =  max(izlow(i),k)
    
        enddo
      enddo
!
      do k = 1,km
        do i =1,npt
          j         = ipt(i)
          vtj(i,k)  = t1(j,k)  * (1.+fv*q1(j,k))
          vtk(i,k)  = vtj(i,k) / prslk(j,k)
          ro(i,k)   = rdi * prsl(j,k) / vtj(i,k)       ! density mid-levels
          taup(i,k) = 0.0
        enddo
      enddo
!
! perform ri_n or ri_mf computation for both OGW and OBL
!
      do k = 1,kmm1
        do i =1,npt
          j         = ipt(i)
          rdz       = 1.   / (zmet(j,k+1) - zmet(j,k))
          tem1      = u1(j,k) - u1(j,k+1)
          tem2      = v1(j,k) - v1(j,k+1)
          dw2       = tem1*tem1 + tem2*tem2
          shr2      = max(dw2,dw2min) * rdz * rdz
!          ti        = 2.0 / (t1(j,k)+t1(j,k+1))
!          bvf2      = grav*(gocp+rdz*(vtj(i,k+1)-vtj(i,k)))* ti
!          ri_n(i,k) = max(bvf2/shr2,rimin)   ! richardson number
!
          bvf2 = grav2 * rdz * (vtk(i,k+1)-vtk(i,k))/ (vtk(i,k+1)+vtk(i,k))
     
          bnv2(i,k+1) = max( bvf2, bnv2min )
          ri_n(i,k+1) = bnv2(i,k)/shr2        ! richardson number consistent with bnv2	
!
! add here computation for "ktur" and ogw-dissipation for the spectral ORO-scheme
!	    
        enddo
      enddo
      k = 1
      do i = 1, npt
        bnv2(i,k) = bnv2(i,k+1)
      enddo
!		
! level iwklm => zmet(j,k) < sigfac * hprime(j) < zmet(j,k+1) 
!
      do i = 1, npt
        j   = ipt(i)
        k_zlow = izlow(i)
        if (k_zlow == iwklm(i)) k_zlow = 1
        delks(i)   = 1.0 / (prsi(j,k_zlow) - prsi(j,iwklm(i)))
!       delks1(i)  = 1.0 /(prsl(j,k_zlow) - prsl(j,iwklm(i)))
        ubar (i)   = 0.0
        vbar (i)   = 0.0
        roll (i)   = 0.0
        pe   (i)   = 0.0
        ek   (i)   = 0.0
        bnv2bar(i) = 0.0  
		
        do k = k_zlow, iwklm(i)-1                      !   kreflm(i)= iwklm(i)-1   
          rdelks  = del(j,k) * delks(i)
          ubar(i) = ubar(i)  + rdelks * u1(j,k)          	  
          vbar(i) = vbar(i)  + rdelks * v1(j,k)          	  
          roll(i) = roll(i)  + rdelks * ro(i,k)          	   
          bnv2bar(i) = bnv2bar(i) + .5*(bnv2(i,k)+bnv2(i,k+1))* rdelks
        enddo		 
      enddo      
!
      do i = 1, npt
        j = ipt(i)
!
! integrate from ztoph = sigfac*hprime  down to zblk if exists
! find ph_blk, dz_blk as introduced in  LM-97 and ifs
!	
        ph_blk =0.  
        do k = iwklm(i), 1, -1
	
          phiang   =  atan2(v1(j,k),u1(j,k))
          phiang = theta(j)*rad_to_deg - phiang
	  
          if ( phiang >  pi2h ) phiang = phiang - pi
          if ( phiang < -pi2h ) phiang = phiang + pi
          ang(i,k) = phiang
          uds(i,k) = max(sqrt(u1(j,k)*u1(j,k) + v1(j,k)*v1(j,k)), velmin)
!
          if (idxzb(i) == 0 ) then
            dz_blk = zmeti(j,k+1) - zmeti(j,k)
            pe(i)  =  pe(i) + bnv2(i,k) *( elvmax(j) - zmet(j,k) ) * dz_blk

            up(i)  =  max(uds(i,k) * cos(ang(i,k)), velmin)  
            ek(i)  = 0.5 *  up(i) * up(i) 

            ph_blk = ph_blk + dz_blk*sqrt(bnv2(i,k))/up(i)

! --- dividing stream lime  is found when pe =exceeds ek. oper-l gfs
!           if ( pe(i) >=  ek(i) ) then
! --- LM97
            if ( ph_blk >=  fcrit_gfs ) then
               idxzb(i) = k
               zobl (j) = zmet(j, k)
               rdxzb(j) = real(k, kind=kind_phys)
            endif

          endif
        enddo
!
!       fcrit_gfs/fr_flow	
! 
        goto 788
!
! alternative expression for blocking:
!         zobl = max(heff*(1. -fcrit_gfs/fr_Flow), 0)
! 
!

          bnv     = sqrt( bnv2bar(i) )
          heff    = 2.*min(hprime(j),hpmax)
          zw2     = ubar(i)*ubar(i)+vbar(i)*vbar(i)
          ulow(i) = sqrt(max(zw2,dw2min))
          fr      = heff*bnv/ulow(i)
          zw1     = max(heff*(1. -fcrit_gfs/fr), 0.0)
          zw2     = zmet(j,2)
	
         if (fr > fcrit_gfs .and. zw1 > zw2 ) then 
           do k=2, kmm1
            zmetp =  zmet(j,k+1) 
            zmetk   =  zmet(j,k)   
            if (zw1 <= zmetp .and. zw1 >= zmetk)  exit
           enddo
            idxzb(i) = k
            zobl (j) = zmet(j, k)
         endif	
788     continue
!
! --- the drag for the blocked flow
!
        if ( idxzb(i) > 0 ) then
!	
! (4.16)-ifs description
!	  
          gam2 = gamma(j)*gamma(j)
          bgam = 1.0 - 0.18*gamma(j) - 0.04*gam2
          cgam =       0.48*gamma(j) + 0.30*gam2
	  
          do k = idxzb(i)-1, 1, -1
!
! empirical height dep-nt "blocking" length from LM-1997
!	  
            zlen = sqrt( (zobl(j)-zmet(j,k) )/(zmet(j,k ) + hprime(j)) )
!	    
!	    
            tem     = cos(ang(i,k))
            cosang2 = tem * tem
            sinang2 = 1.0 - cosang2 
!	      
!  cos =1 sin =0 =>   1/r= gam     zr = 2.-gam 
!  cos =0 sin =1 =>   1/r= 1/gam   zr = 2.- 1/gam
!
            rdem = cosang2      +  gam2 * sinang2
            rnom = cosang2*gam2 +         sinang2
!	       
! metoffice dec 2010
! correction of H. Wells & A. Zadra for the
! aspect ratio  of the hill seen by mean flow
! (1/r , r-inverse below: 2-r)

            rdem = max(rdem, 1.e-6)       
            r    = sqrt(rnom/rdem)
            zr   =  max( 2. - r, 0. )
            sigres = max(sigmin, sigma(j))
	    	    
            mtbridge = zr * sigres*zlen / hprime(j)
! (4.15)-ifs 	   
!           dbtmp = cdmb4 * mtbridge *                               &
!     &           max(cos(ang(i,k)), gamma(j)*sin(ang(i,k)))
! (4.16)-ifs
            dbtmp  = cdmb4*mtbridge*(bgam* cosang2 +cgam* sinang2)
!
! linear damping due to OBL [1/sec]=[U/L_block_orthogonal]
! more accurate along 2-axes of ellipse, here zr-factor is based on Phillips' analytics
!	    
            db(i,k)= dbtmp * uds(i,k)
          enddo
!                  
        endif
      enddo
!.............................
!.............................
!   end  mtn blocking section
!.............................
!.............................
!
!--- OGW section
!     
!  scale cleff between im=384*2 and 192*2 for t126/t170 and t62
!  inside "cires_ugwp_initialize.f90" now
!
      kmpbl  = km / 2 
      iwk(1:npt) = 2
!
! meto/UK-scheme: 
! k_mtb = max(k_zmtb, k_n*hprime/2] to reduce diurnal variations taub_ogw 
!     
      do k=3,kmpbl
        do i=1,npt
          j   = ipt(i)
          tem = (prsi(j,1) - prsi(j,k))
          if (tem < dpmin) iwk(i) = k           ! dpmin=50 mb	
        enddo	  
      enddo
!
! iwk - adhoc gfs-parameter to select ogw-launch level between
!      level ~0.4-0.5 km from surface or/and  HPBL-top
!
! in all UGWP-schemes: zogw > zobl
!      in ugwp-v1: options to modify as  htop ~ (2-3)*hprime > zmtb
! 
!

      kbps  = 1
      kmps  = km
      k_mtb = 1
      
      do i=1,npt
        j         = ipt(i)
        k_mtb     = max(1, idxzb(i))
                                                        !  WRF/GSL: kogw .ge. kbl
        kref(i)   = max(iwk(i),  kpbl(j)+1 )            ! reference level pbl or smt-else ???? Zogw= sigfac*Hprime
        kref(i)   = max(kref(i), iwklm(i) )             ! iwklm => sigfac*hprime
!	
! zogw > zobl
!
        if (kref(i) <= k_mtb)  kref(i) = k_mtb + 1      ! layer above blocking
        kbps      = max(kbps,  kref(i))
        kmps      = min(kmps,  kref(i))
!
        delks(i)  = 1.0 / (prsi(j,k_mtb) - prsi(j,kref(i)))
        ubar (i)  = 0.0
        vbar (i)  = 0.0
        roll (i)  = 0.0
        bnv2bar(i)= 0.0
      enddo
!
      kbpsp1 = kbps + 1
      kbpsm1 = kbps - 1
      k_mtb  = 1
!
!====================== we estimate MF-parameters from k= k_mtb to [kref~kpbl] > k_mtb
!
      do i = 1,npt
        k_mtb = max(1, idxzb(i))
        do k = k_mtb,kbps                              !kbps = max(kref) ;kmps= min(kref)
          if (k < kref(i)) then
            j          = ipt(i)
            rdelks     = del(j,k) * delks(i)
            ubar(i)    = ubar(i)  + rdelks * u1(j,k)   ! mean u below kref
            vbar(i)    = vbar(i)  + rdelks * v1(j,k)   ! mean v below kref
            roll(i)    = roll(i)  + rdelks * ro(i,k)   ! mean ro below kref
            bnv2bar(i) = bnv2bar(i) + .5*(bnv2(i,k)+bnv2(i,k+1))* rdelks
          endif
        enddo
      enddo
!
! orographic asymmetry parameters (oa), and (clx)  [Kim & Arakawa Kim & Doyle]
!
      do i = 1,npt
        j      = ipt(i)
        wdir   = atan2(ubar(i),vbar(i))    + pi   ! not sure about "+pi" due to "nwdir"-Kim OA/CLX-processing
        idir   = mod(nint(fdir*wdir),mdir) + 1
        nwd    = nwdir(idir)
        oa(i)  = (1-2*int( (nwd-1)/4 )) * oa4(j,mod(nwd-1,4)+1)
        clx(i) = clx4(j,mod(nwd-1,4)+1)   ! number of "eefective" hills in the grid-box KA-95/KD-05
!
!      
       dtfac(i)  = 1.0
       icrilv(i) = .false.                      ! initialize critical level control Logic
       
       ulow(i) = max(sqrt(ubar(i)*ubar(i)+vbar(i)*vbar(i)),velmin)
       xn(i)  = ubar(i) / ulow(i)
       yn(i)  = vbar(i) / ulow(i)
      enddo
!
      do  k = 1, kmm1
        do  i = 1,npt
          j            = ipt(i)
          velco(i,k)   = 0.5 * ((u1(j,k)+u1(j,k+1))*xn(i)+  (v1(j,k)+v1(j,k+1))*yn(i))
        enddo
      enddo
      
       do  i = 1,npt
	 velco(i,km) =  velco(i,kmm1)
       enddo      
!
!------------------------------------------------------------------------
! v0: incorporates latest modifications for kxridge and heff/hsat
!             and taulin for fr <=fcrit_gfs 
!             and concept of "clipped" hill if zmtb > 0. to make
! the integrated "tau_sso = tau_ogw +tau_mtb" close to reanalysis data
!      it is still used the "single-orowave"-approach along ulow-upwind
!
! in contrast to the 2-orthogonal wave (2otw) schemes of ifs/meto/e-canada  
! 2-wave scheme requires "aver angle" and wind projections on axes of ell-hill
!     with 2-stresses:  taub_a & taub_b as of Phillips  (1984)
!------------------------------------------------------------------------
      
      taub(:)  = 0. ; taulin(:)= 0.
      
      do i = 1,npt
        j    = ipt(i)
        bnv  = sqrt( bnv2bar(i) )
        heff = min(hprime(j),hpmax)

        if( zobl(j) > 0.) heff = max(sigfac*heff-zobl(j), 0.)/sigfac  
	
        if (heff <= 0) cycle
        zw1 = ulow(i)/bnv
        hsat = fcrit_gfs *zw1
        heff = min(heff, hsat)

        fr   = min(heff/zw1, frmax)
	fr2 = fr*fr
!
! [Kim & Doyle, 2005]
!
        efact    = (oa(i) + 2.) ** (ceofrc*fr)     ! enhnancement factor due to the resonance ampification/downstream
        efact    = min( max(efact,efmin), efmax )
	
        mkd05_hills(i)    = (1. + clx(i)) ** (oa(i)+1.)   ! number of mountains with some anizoropy of KD-2005
	
!
! cleff_max(C768 = 6.28/2.5 km) .....
!	
!       xlinv(i) = min(coefm * cleff, cleff_max)  
!
        xlinv(i) = min(cleff, cleff_max)  		 
        gfobnv   = efact* gmax * oc(j)/(fr2*oc(j) + cg) 
!	
! tem = fr2*oc(j) ;	 gfobnv   = gmax * tem / ((tem + cg)*bnv(i))	
! kx =or max(kxridge, inv_b2eff)  ! 6.28/lx ..0.5*sigma(j)/heff = 1./lridge  
!
        sigres = sigma(j)
	
        inv_b2eff =  pi*sigres/heff              ! pi2/(2b)
	kxridge   =  pi /ahdxres(i)              ! pi2/(2*dx)	
		
        xlingfs  =  max(inv_b2eff, kxridge)
!		
	xlinv(i) = max(xlingfs,  xlinv(i) )
	
	nonh_fact = 1. - xlinv(i)*zw1 * xlinv(i)*zw1 
	
	if ( nonh_fact <= 0.) cycle                        ! non-hydrostatic trapping kx =kz = N/U
!	
        taulin(i) = xlinv(i)*roll(i)*bnv*ulow(i)*heff*heff * nonh_fact              	
!	
!        zw2  = xlinv(i) * roll(i) *ulow(i)*ulow(i)* zw1 * gfobnv
!
! fr2 = (bnv*heff/Ulow)**2 + non-hydrostatic trapping effects  Fr2_nh = Fr2 - kx2*heff^2
! 	
        if ( fr > fcrit_gfs ) then
	   gfobnv = max(gfobnv, 2.0)                        ! can be 2-times higher than taulin
           taub(i)  =  taulin(i) * gfobnv                   ! nonlinear flux tau0...xlinv(i)		             
        else
           taub(i)  = taulin(i)                             !  linear flux for fr <= fcrit_gfs	 
        endif		 
!
        k       = max(1, kref(i)-1)
        tem     = max(velco(i,k)*velco(i,k), dw2min)
        scor(i) = bnv2(i,k) / tem           ! scorer parameter below kref level Bn2/U2
!
! diagnostics for zogw, tau_ogw
!
        zogw(j)    = zmeti(j, kref(i) )
	tau_ogw(j) = taub(i)
      enddo
!                                                                       
!----set up bottom values of stress
!
        do i = 1,npt
          taup(i, 1:kref(i) ) = taub(i)
        enddo
      
      if (strsolver == 'pss-1986') then     
       
!======================================================
!   v0-gfs orogw-solver of palmer et al 1986 -"pss-1986"
!     modified by KD05 with the expression (11):below k=kref ???
!     tau(k+1) = tau(k)*Scorer(K+1)/Scorer(K) 
!
!   in v1-orogw  linsatdis of                 "wam-2017"
!     with llwb-mechanism for
!     rotational/non-hydrostat ogws important for 
!     highres-fv3gfs with dx < 10 km
!======================================================
      
        do k = kmps, kmm1                   ! vertical level loop from min(kref)
          kp1 = k + 1
	  
          do i = 1, npt
             if (k >= kref(i)) then
              icrilv(i) = icrilv(i) .or. ( ri_n(i,k) < ric).or. (velco(i,k) <= 0. )
             endif
          enddo
!
          do i = 1,npt
            if (k >= kref(i))   then
             if (.not.icrilv(i) .and. taup(i,k) > 0.0 ) then
                temv = 1.0 / max(velco(i,k), velmin)
!===============
! Condition for levels below kref(i): k+1 < kref(i))  ??? see KD05 expression (11) for LLWB ??? only OA >0
!
                if (oa(i) > 0. .and. kp1 < kref(i)) then				
                  scork   = bnv2(i,k) * temv * temv
                  rscor   = min(1.0, scork / scor(i))
                  scor(i) = scork
                else 
                  rscor   = 1.
                endif
!===============
                brvf = sqrt(bnv2(i,k))        ! brent-vaisala frequency interface
!               tem1 = xlinv(i)*(ro(i,kp1)+ro(i,k))*brvf*velco(i,k)*0.5

                tem1 = xlinv(i)*(ro(i,kp1)+ro(i,k))*brvf*0.5* max(velco(i,k), velmin)
                hd   = sqrt(taup(i,k) / tem1)
                fro  = brvf * hd * temv
!
!    rim is the  "wave"-richardson number byPalmer,Shutts & Swinbank 1986 , PSS-1986
!
                tem2   = sqrt(ri_n(i,k))
                tem    = 1. + tem2 * fro
                ri_gw  = ri_n(i,k) * (1.0-fro) / (tem * tem)
!
!    check Ri-stability to employ the 'dynamical saturation hypothesis' PSS-1986
!    assuming co-existence of Dyn-Ins and Conv-Ins
!                                       
                if (ri_gw <= ric .and.(oa(i) <= 0. .or.  kp1 >= kref(i) )) then
                   temc = 2.0 + 1.0 / tem2
                   hd   = velco(i,k) * (2.*sqrt(temc)-temc) / brvf
                   taup(i,kp1) = tem1 * hd * hd
                else 
		
                   taup(i,kp1) = taup(i,k) * rscor
                endif
!	      
                taup(i,kp1) = min(taup(i,kp1), taup(i,k))
              endif         !  if (.not.icrilv(i) .and. taup(i,k) > 0.0 )
            endif           !  k >= kref(i))
          enddo             !  oro-points
        enddo               ! do k = kmps, kmm1  vertical level loop 
!     
!  zero momentum deposition at the top model layer: taup(k+1) = taup(k)
!      
        taup(1:npt,km+1) = taup(1:npt,km)      
!
!     calculate wave acc-n: - (grav)*d(tau)/d(p) = taud
!
        do k = 1,km
          do i = 1,npt
            taud(i,k) = grav*(taup(i,k+1) - taup(i,k))/del(ipt(i),k)
!======================================================================================
! we estimated "impact" of the single sub-grid hill, we have "arhills" in the grid-box
! 2-estimations of "nhills":  1) geometry-arhills and 2) KDO5 mkd05_hills
!           for OBL we used:  1) nhills=Grid_Area/Hill_area 
!            nhills = max(mkd05_hills(i), arhills(i))
! Trapped "Lee" downslope wave regimes are not properly modelled: vertical shear +NH/Nonlin
!               tau(z) = const =>  tau(z)/m2(z) = const (empirical mesoscale)
!====================================================================================== 	    
	    taud(i,k) = taud(i,k)*arhills(i)    ! simple scale-awareness nhills=Grid_Area/Hill_area 
          enddo
        enddo
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!------if the gravity wave drag would force a critical line in the
!------layers below sigma=rlolev during the next deltim timestep,
!------then only apply drag until that critical line is reached.
! empirical implementation of the llwb-mechanism: lower level wave breaking
! by limiting "ax = dtfac*ax" due to possible llwb around kref and 500 mb
! critical line [v - ax*dtp = 0.] is smt like "llwb" for stationary ogws
!2019:  this option limits sensitivity of taux/tauy to variations in "taub"
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        do k = 1,kmm1
          do i = 1,npt
           if (k  >= kref(i) .and. prsi(ipt(i),k) >= rlolev) then

               if(taud(i,k) /= 0.) then
                 tem = dtp * taud(i,k)                          ! tem = du/dt-oro*dt => U/dU vs 1
                 dtfac(i) = min(dtfac(i),abs(velco(i,k)/tem))   ! reduce Ax= Ax*(1, or U/dU <=1)
!	         dtfac(i) = 1.0
               endif
            endif
          enddo
        enddo
!
!--------- orogw-solver of gfs PSS-1986 is performed 
  
      else 

!-----------Unified orogw-solver of wam2017 out :   taup, taud, pkdis
   
        dtfac(:) =  1.0  
	     
        call oro_wam_2017(im, km, npt, ipt, kref, kdt, me, master,             &
           dtp, dxres, taub, u1, v1, t1, xn, yn, bnv2, ro, prsi,prsl,          &
           grav, con_omega, rd,                                                &
           del, sigma, hprime, gamma, theta, sinlat, xlatd, taup, taud, pkdis)  
	      
      endif            !  oro_linsat - linsatdis-solver for stationary OGWs
!      
!---- above orogw-solver of wam2017------------
! 
! tofd as in Beljaars-2004 IFS sep-scale ~5km
!                   CESM ~ 6km (TMS + OGW/OBL)
!  sgh30 = varss of GSL (?)
! ----------------------------------------------
    
      if( do_tofd ) then
	
        if ( kdt == 1 .and. me == 0) then
          print *, 'ugwp-v1 do_tofd  from surface to ', ztop_tofd 
        endif
	
        do i = 1,npt          
          j = ipt(i)
          zpbl  = zmet( j, kpbl(j) )
 
          sigflt = min(sgh30(j), 0.3*hprime(j))    ! cannot exceed 30% of ls-sso ! GSL-limit 250 m ; var_maxfd =150m
          zsurf = zmeti(j,1)
          do k=1,km
            zpm(k) = zmet(j,k)
            up1(k) = u1(j,k)
            vp1(k) = v1(j,k)
          enddo
 
          call ugwp_tofd1d(km, cpd, dtp, sigflt, zsurf, zpbl,  &
               up1, vp1, zpm,  utofd1, vtofd1, epstofd1, krf_tofd1)
     
          do k=1,km
            dudt_ofd(j,k) = utofd1(k)
            dvdt_ofd(j,k) = vtofd1(k)
!	 
! add tofd to gw-tendencies
!	 
            pdvdt(j,k)  = pdvdt(j,k) + utofd1(k)
            pdudt(j,k)  = pdudt(j,k) + vtofd1(k)
	    pdtdt(j,k)  = pdtdt(j,k) + epstofd1(k)
          enddo
!2018-diag
           du_ofdcol(j) = sum( utofd1(1:km)* del(j,1:km))
           dv_ofdcol(j) = sum( vtofd1(1:km)* del(j,1:km))
	   
	    dusfc(j)   = dusfc(j) + du_ofdcol(j)
            dvsfc(j)   = dvsfc(j) + dv_ofdcol(j)	  
        enddo
      endif             ! do_tofd 

!--------------------------------------------
! combine oro-drag effects MB +TOFD + OGWs +  diag-3d
!--------------------------------------------  
! 
 
      do k = 1,km
        do i = 1,npt
          j    = ipt(i)
!
          eng0 = 0.5*(u1(j,k)*u1(j,k)+v1(j,k)*v1(j,k))
!
          if ( k < idxzb(i) .and. idxzb(i) /= 0 ) then
!
! if blocking layers -- no ogws
!	  
            dbim       = db(i,k) / (1.+db(i,k)*dtp)
	    
            dudt_obl(j,k) = -dbim * u1(j,k)
            dvdt_obl(j,k) = -dbim * v1(j,k)
	    	    
            pdvdt(j,k) = dudt_obl(j,k) +pdvdt(j,k)
            pdudt(j,k) = dvdt_obl(j,k) +pdudt(j,k)	         	    
!2018-diag 	    	    
            du_oblcol(j)    = du_oblcol(j) + dudt_obl(j,k)* del(j,k)
            dv_oblcol(j)    = dv_oblcol(j) + dvdt_obl(j,k)* del(j,k)
	    
            dusfc(j)   = dusfc(j) + du_oblcol(j)
            dvsfc(j)   = dvsfc(j) + dv_oblcol(j)
	    	    
          else
!
! ogw-s above blocking height
!	  
            taud(i,k)  = taud(i,k) * dtfac(i)
            dtaux      = taud(i,k) * xn(i) 
            dtauy      = taud(i,k) * yn(i) 
!	    
            dudt_ogw(j,k) = dtaux
            dvdt_ogw(j,k) = dtauy 
!	       
            pdvdt(j,k)   = dtauy  +pdvdt(j,k)
            pdudt(j,k)   = dtaux  +pdudt(j,k)
 
!	    
	    du_ogwcol(j) = du_ogwcol(j) + dtaux * del(j,k)
	    dv_ogwcol(j) = dv_ogwcol(j) + dtauy * del(j,k)	    
!
            dusfc(j) = dusfc(j)  + du_ogwcol(j)
            dvsfc(j) = dvsfc(j)  + dv_ogwcol(j)
          endif
!============ 	  
! local energy deposition sso-heat due to loss of kinetic energy
!============ 
            unew   = u1(j,k) + pdudt(j,k)*dtp             !   pdudt(j,k)*dtp
            vnew   = v1(j,k) + pdvdt(j,k)*dtp             !   pdvdt(j,k)*dtp
            eng1   = 0.5*(unew*unew + vnew*vnew)	
            pdtdt(j,k) = max(eng0-eng1,0.)*rcpdt + pdtdt(j,k) 
	    
        enddo
      enddo
! dusfc w/o tofd  sign as in the era-i, merra  and cfsr
      do i = 1,npt
         j           = ipt(i)
         dusfc(j)    = -rgrav * dusfc(j)
         dvsfc(j)    = -rgrav * dvsfc(j)
         du_oblcol(j)  = -rgrav *du_oblcol (j)
         tau_ogw(j)  = -rgrav * tau_ogw(j)
         du_ofdcol(j) = -rgrav * du_ofdcol(j)
       enddo

       return


!============ debug ------------------------------------------------
       if (kdt <= 2 .and. me == 0) then
        print *, 'vgw-oro done gwdps_v0 in ugwp-v0 step-proc ', kdt, me
!
        print *, maxval(pdudt)*86400.,  minval(pdudt)*86400, 'vgw_axoro'
        print *, maxval(pdvdt)*86400.,  minval(pdvdt)*86400, 'vgw_ayoro'
!       print *, maxval(kdis),  minval(kdis),  'vgw_kdispro m2/sec'
        print *, maxval(pdtdt)*86400.,  minval(pdtdt)*86400,'vgw_epsoro'
!        print *, maxval(zobl), ' z_mtb ',  maxval(tau_mtb), ' tau_mtb '
        print *, maxval(zogw), ' z_ogw ',  maxval(tau_ogw), ' tau_ogw '
!       print *, maxval(tau_tofd),  ' tau_tofd '
!       print *, maxval(axtms)*86400.,  minval(axtms)*86400, 'vgw_axtms'
!       print *,maxval(dudt_mtb)*86400.,minval(dudt_mtb)*86400,'vgw_axmtb'
        if (maxval(abs(pdudt))*86400. > 100.) then

          print *, maxval(u1),  minval(u1),  ' u1 gwdps-v1 '
          print *, maxval(v1),  minval(v1),  ' v1 gwdps-v1 '
          print *, maxval(t1),  minval(t1),  ' t1 gwdps-v1 '
          print *, maxval(q1),  minval(q1),  ' q1 gwdps-v1 '
          print *, maxval(del), minval(del), ' del gwdps-v1 '
          print *, maxval(zmet),minval(zmet), 'zmet'
          print *, maxval(zmeti),minval(zmeti), 'zmeti'
          print *, maxval(prsi), minval(prsi), ' prsi '
          print *, maxval(prsl), minval(prsl), ' prsl '
          print *, maxval(ro), minval(ro), ' ro-dens '
          print *, maxval(bnv2(1:npt,:)), minval(bnv2(1:npt,:)),' bnv2 '
          print *, maxval(kpbl), minval(kpbl), ' kpbl '
          print *, maxval(sgh30), maxval(hprime), maxval(elvmax),'oro-d'
          print *
          do i =1, npt
            j= ipt(i)
            print *,zogw(j)/hprime(j), zobl(j)/hprime(j), &
                   zmet(j,1)*1.e-3, nint(hprime(j)/sigma(j))
!    
          enddo
          print *
          stop
        endif
       endif
       
      return
      end subroutine orogw_v1 
!
!      
     subroutine ugwp_tofd1d(levs, con_cp, dtp, sigflt, zsurf, zpbl,  u, v, &
                            zmid, utofd, vtofd, epstofd, krf_tofd)
			    
       use machine ,      only : kind_phys 
       use ugwp_oro_init, only : n_tofd, const_tofd, ze_tofd, a12_tofd, ztop_tofd
!
! adding the implicit tendency estimate
!       
     implicit none
      integer,         intent(in) ::  levs
      real(kind_phys), intent(in) :: con_cp
      real(kind_phys), intent(in) :: dtp  
          
      real(kind_phys), intent(in), dimension(levs)  ::   u, v, zmid
      real(kind_phys), intent(in)                   ::  sigflt, zpbl, zsurf
      
      real(kind_phys), intent(out), dimension(levs) ::  utofd, vtofd, epstofd, krf_tofd

         
!
! locals
!
     integer :: i, k
     real    :: rcpd2
     real(kind_phys)    ::  unew, vnew, eknew          
     real, parameter    :: sghmax = 5.        ! dz(1)/5= 25/5 m dz-of the first layer
     real, parameter    :: tend_imp = 0.          
     
     real    :: sgh2, ekin, zdec, rzdec, umag, zmet, zarg, ztexp, krf
!     
       utofd =0.0 ; vtofd = 0.0 ;  epstofd =0.0 ; krf_tofd =0.0    
       rcpd2 = 0.5/con_cp
!         
       zdec = max(n_tofd*sigflt, zpbl)          ! ntimes*sgh_turb or Zpbl
       zdec = min(ze_tofd, zdec)                ! cannot exceed ~15 km
       rzdec = 1.0/zdec
       sgh2 = max(sigflt*sigflt, sghmax*sghmax) ! 25 meters dz-of the first layer
       
      do k=1, levs  
         zmet = zmid(k)-zsurf   
      if (zmet > ztop_tofd) cycle
         ekin = u(k)*u(k) + v(k)*v(k)
         umag = sqrt(ekin)
         zarg = zmet*rzdec
         ztexp = exp(-zarg*sqrt(zarg))
         krf   = const_tofd* a12_tofd *sgh2* zmet ** (-1.2) *ztexp
	 
         if (tend_imp == 1.) then
	    krf = krf/(1.+krf*dtp)
	 endif
	 
         utofd(k)    = -krf*u(k)
         vtofd(k)    = -krf*v(k)
         if (tend_imp == 1.)      then
	   unew =u(k)+ utofd(k)*dtp ; vnew =v(k)+ vtofd(k)*dtp 
	   eknew =unew*unew + vnew*vnew
	   epstofd(k)  = rcpd2*(ekin-eknew)
	  else
	   epstofd(k)  = rcpd2*krf*ekin
	 endif
	                               ! more accurate heat/mom form using "implicit tend-solver"
                                       ! to update momentum and temp-re; epstofd(k) can be skipped
         krf_tofd(k) = krf             ! can be used as addition to the mesoscale blocking
     enddo
!                
     end subroutine ugwp_tofd1d    

end module cires_ugwp_orolm97_v1
