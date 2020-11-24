!>  \file unified_ugwp.F90
!! This file combines three gravity wave drag schemes under one ("unified_ugwp") suite:
!!      1) The "V0 CIRES UGWP" scheme (cires_ugwp.F90) as implemented in the FV3GFSv16 atmosphere model, which includes:
!!            a) the "traditional" EMC orograhic gravity wave drag and flow blocking scheme of gwdps.f
!!            b) the v0 cires ugwp non-stationary GWD scheme
!!      2) The GSL orographic drag suite (drag_suite.F90), as implmeneted in the RAP/HRRR, which includes:
!!            a) large-scale gravity wave drag and low-level flow blocking -- active at horizontal scales
!!               down to ~5km (Kim and Arakawa, 1995 \cite kim_and_arakawa_1995; Kim and Doyle, 2005 \cite kim_and_doyle_2005)
!!            b) small-scale gravity wave drag scheme -- active typically in stable PBL at horizontal grid resolutions down to ~1km
!!               (Steeneveld et al, 2008 \cite steeneveld_et_al_2008; Tsiringakis et al, 2017 \cite tsiringakis_et_al_2017)
!!            c) turbulent orographic form drag -- active at horizontal grid ersolutions down to ~1km
!!               (Beljaars et al, 2004 \cite beljaars_et_al_2004)
!!      3) The "V1 CIRES UGWP" scheme developed by Valery Yudin (University of Colorado, CIRES)
!! See Valery Yudin's presentation at 2017 NGGPS PI meeting:
!! Gravity waves (GWs): Mesoscale GWs transport momentum, energy (heat) , and create eddy mixing in the whole atmosphere domain; Breaking and dissipating GWs deposit: (a) momentum; (b) heat (energy); and create (c) turbulent mixing of momentum, heat, and tracers
!! To properly incorporate GW effects (a-c) unresolved by DYCOREs we need GW physics
!! "Unified": a) all GW effects due to both dissipation/breaking; b) identical GW solvers for all GW sources; c) ability to replace solvers.
!! Unified Formalism:
!! 1. GW Sources: Stochastic and physics based mechanisms for GW-excitations in the lower atmosphere, calibrated by the high-res analyses/forecasts, and observations (3 types of GW sources: orography, convection, fronts/jets).
!! 2. GW Propagation: Unified solver for "propagation, dissipation and breaking" excited from all type of GW sources.
!! 3. GW Effects: Unified representation of GW impacts on the "resolved" flow for all sources (energy-balanced schemes for momentum, heat and mixing).
!! https://www.weather.gov/media/sti/nggps/Presentations%202017/02%20NGGPS_VYUDIN_2017_.pdf
!!
!! The unified_ugwp scheme is activated by gwd_opt = 2 in the namelist.
!! The choice of schemes is activated at runtime by the following namelist options (boolean):
!!       do_ugwp_v0           -- activates V0 CIRES UGWP scheme - both orographic and non-stationary GWD
!!       do_ugwp_v0_orog_only -- activates V0 CIRES UGWP scheme - orographic GWD only
!!       do_gsl_drag_ls_bl    -- activates RAP/HRRR (GSL) large-scale GWD and blocking
!!       do_gsl_drag_ss       -- activates RAP/HRRR (GSL) small-scale GWD
!!       do_gsl_drag_tofd     -- activates RAP/HRRR (GSL) turbulent orographic drag
!!       do_ugwp_v1           -- activates V1 CIRES UGWP scheme - both orographic and non-stationary GWD
!!       do_ugwp_v1_orog_only -- activates V1 CIRES UGWP scheme - orographic GWD only
!! Note that only one "large-scale" scheme can be activated at a time.
!!

module unified_ugwp

    use machine, only: kind_phys

    use cires_ugwp_module, only: knob_ugwp_version, cires_ugwp_dealloc
    use cires_ugwp_module, only: cires_ugwp_init_modv1, cires_ugwp_init_modv0 
    use cires_ugwp_module, only: calendar_ugwp, ngwflux_update

    use cires_ugwp_module, only: tamp_mpa, tau_min

    use gwdps, only: gwdps_run

    use drag_suite, only: drag_suite_run

    use cires_ugwp_orolm97_v1, only: orogw_v1
    use cires_ugwp_triggers, only:  slat_geos5_2020, slat_geos5_tamp_v1
    use cires_ugwp_triggers, only:  emc_modulation, slat_geos5_tamp_v0
    use cires_ugwp_triggers, only:  get_spectra_tau_convgw,  get_spectra_tau_okw, get_spectra_tau_nstgw  
    
    ! use cires_ugwp_ngw_utils, only: tau_limb_advance
    
    use cires_ugwp_solv2_v1_mod, only: cires_ugwp_solv2_v1

    implicit none

    private

    public unified_ugwp_init, unified_ugwp_run, unified_ugwp_finalize

    logical :: is_initialized = .False.
!    use cires_ugwp_module, only: pa_rf_in, tau_rf_in
contains

! ------------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme v0
! ------------------------------------------------------------------------
!>@brief The subroutine initializes the unified UGWP
!> \section arg_table_unified_ugwp_init Argument Table
!! \htmlinclude unified_ugwp_init.html
!!
! -----------------------------------------------------------------------
!
    subroutine unified_ugwp_init (me, master, nlunit, input_nml_file, logunit, &
                fn_nml2, jdat, lonr, latr, levs, ak, bk, dtp, cdmbgwd, cgwf,   &
                con_pi, con_rerth, con_p0,                                     &
                do_ugwp,do_ugwp_v0, do_ugwp_v0_orog_only, do_gsl_drag_ls_bl,   &
                do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1,                  &
                do_ugwp_v1_orog_only, errmsg, errflg)

!----  initialization of unified_ugwp
    implicit none

    integer,              intent (in) :: me
    integer,              intent (in) :: master
    integer,              intent (in) :: nlunit
    character(len=*),     intent (in) :: input_nml_file(:)
    integer,              intent (in) :: logunit
    integer,              intent(in)  :: jdat(8)
    integer,              intent (in) :: lonr
    integer,              intent (in) :: levs
    integer,              intent (in) :: latr
    real(kind=kind_phys), intent (in) :: ak(levs+1), bk(levs+1)
    real(kind=kind_phys), intent (in) :: dtp
    real(kind=kind_phys), intent (in) :: cdmbgwd(4), cgwf(2) ! "scaling" controls for "old" GFS-GW schemes

!    real(kind=kind_phys), intent (in) :: pa_rf_in, tau_rf_in

    real(kind=kind_phys), intent (in) :: con_p0, con_pi, con_rerth
    logical,              intent (in) :: do_ugwp
    logical,              intent (in) :: do_ugwp_v0, do_ugwp_v0_orog_only,  &
                                         do_gsl_drag_ls_bl, do_gsl_drag_ss, &
                                         do_gsl_drag_tofd, do_ugwp_v1,      &
                                         do_ugwp_v1_orog_only

    character(len=*), intent (in) :: fn_nml2
    !character(len=*), parameter   :: fn_nml='input.nml'

    integer :: ios
    logical :: exists
    real    :: dxsg
    integer :: k

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0


    ! Test to make sure that at most only one large-scale/blocking
    ! orographic drag scheme is chosen
    if ( (do_ugwp_v0.and.(do_ugwp_v0_orog_only.or.do_gsl_drag_ls_bl.or.    &
                          do_ugwp_v1.or.do_ugwp_v1_orog_only))        .or. &
         (do_ugwp_v0_orog_only.and.(do_gsl_drag_ls_bl.or.do_ugwp_v1.or.    &
                                    do_ugwp_v1_orog_only))            .or. &
         (do_gsl_drag_ls_bl.and.(do_ugwp_v1.or.do_ugwp_v1_orog_only)) .or. &
         (do_ugwp_v1.and.do_ugwp_v1_orog_only) ) then

       write(errmsg,'(*(a))') "Logic error: Only one large-scale&
          &/blocking scheme (do_ugwp_v0,do_ugwp_v0_orog_only,&
          &do_gsl_drag_ls_bl,do_ugwp_v1 or &
          &do_ugwp_v1_orog_only) can be chosen"
       errflg = 1
       return

    end if


    if (is_initialized) return


    if ( do_ugwp_v0 ) then
       ! if (do_ugwp .or. cdmbgwd(3) > 0.0) then (deactivate effect of do_ugwp)
       if (cdmbgwd(3) > 0.0) then
         call cires_ugwp_init_modv0 (me, master, nlunit, input_nml_file, logunit, &
                                fn_nml2, lonr, latr, levs, ak, bk, con_p0, dtp, &
                                cdmbgwd(1:2), cgwf)
       else
         write(errmsg,'(*(a))') "Logic error: cires_ugwp_mod_init called but &
               &do_ugwp_v0 is true and cdmbgwd(3) <= 0"
         errflg = 1
         return
       end if
    end if


    if ( do_ugwp_v1 ) then
       call cires_ugwp_init_modv1 (me, master, nlunit, logunit, jdat, con_pi,      &
                                con_rerth, fn_nml2, lonr, latr, levs, ak, bk,   &
                                con_p0, dtp, cdmbgwd(1:2), cgwf, errmsg, errflg)  

    end if

    is_initialized = .true.

    end subroutine unified_ugwp_init


! -----------------------------------------------------------------------
! finalize of unified_ugwp   (_finalize)
! -----------------------------------------------------------------------

!>@brief The subroutine finalizes the CIRES UGWP

!> \section arg_table_unified_ugwp_finalize Argument Table
!! \htmlinclude unified_ugwp_finalize.html
!!

    subroutine unified_ugwp_finalize(do_ugwp_v0,do_ugwp_v1,errmsg, errflg)

    implicit none
!
    logical,          intent (in) :: do_ugwp_v0, do_ugwp_v1
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

    call cires_ugwp_dealloc
    
    is_initialized = .false.

    end subroutine unified_ugwp_finalize


! -----------------------------------------------------------------------
!    originally from ugwp_driver_v0.f
!    driver of cires_ugwp   (_driver)
! -----------------------------------------------------------------------
!   driver is called after pbl & before chem-parameterizations
! -----------------------------------------------------------------------
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
! -----------------------------------------------------------------------
!>@brief These subroutines and modules execute the CIRES UGWP Version 0
!>\defgroup unified_ugwp_run Unified Gravity Wave Physics General Algorithm
!> @{
!! The physics of NGWs in the UGWP framework (Yudin et al. 2018 \cite yudin_et_al_2018) is represented by four GW-solvers, which is introduced in Lindzen (1981) \cite lindzen_1981, Hines (1997) \cite hines_1997, Alexander and Dunkerton (1999) \cite alexander_and_dunkerton_1999, and Scinocca (2003) \cite scinocca_2003. The major modification of these GW solvers is represented by the addition of the background dissipation of temperature and winds to the saturation criteria for wave breaking. This feature is important in the mesosphere and thermosphere for WAM applications and it considers appropriate scale-dependent dissipation of waves near the model top lid providing the momentum and energy conservation in the vertical column physics (Shaw and Shepherd 2009 \cite shaw_and_shepherd_2009). In the UGWP-v0, the modification of Scinocca (2003) \cite scinocca_2003 scheme for NGWs with non-hydrostatic and rotational effects for GW propagations and background dissipation is represented by the subroutine \ref fv3_ugwp_solv2_v0. In the next release of UGWP, additional GW-solvers will be implemented along with physics-based triggering of waves and stochastic approaches for selection of GW modes characterized by horizontal phase velocities, azimuthal directions and magnitude of the vertical momentum flux (VMF).
!!
!! In UGWP-v0, the specification for the VMF function is adopted from the GEOS-5 global atmosphere model of GMAO NASA/GSFC, as described in Molod et al. (2015) \cite molod_et_al_2015 and employed in the MERRRA-2 reanalysis (Gelaro et al., 2017 \cite gelaro_et_al_2017). The Fortran subroutine \ref slat_geos5_tamp describes the latitudinal shape of VMF-function as displayed in Figure 3 of Molod et al. (2015) \cite molod_et_al_2015. It shows that the enhanced values of VMF in the equatorial region gives opportunity to simulate the QBO-like oscillations in the equatorial zonal winds and lead to more realistic simulations of the equatorial dynamics in GEOS-5 operational and MERRA-2 reanalysis products. For the first vertically extended version of FV3GFS in the stratosphere and mesosphere, this simplified function of VMF allows us to tune the model climate and to evaluate multi-year simulations of FV3GFS with the MERRA-2 and ERA-5 reanalysis products, along with temperature, ozone, and water vapor observations of current satellite missions. After delivery of the UGWP-code, the EMC group developed and tested approach to modulate the zonal mean NGW forcing by 3D-distributions of the total precipitation as a proxy for the excitation of NGWs by convection and the vertically-integrated  (surface - tropopause) Turbulent Kinetic Energy (TKE). The verification scores with updated NGW forcing, as reported elsewhere by EMC researchers, display noticeable improvements in the forecast scores produced by FV3GFS configuration extended into the mesosphere. However, the EMC-updates are not employed in GFSv16
!!
!> \section arg_table_unified_ugwp_run Argument Table
!! \htmlinclude unified_ugwp_run.html
!!
!> \section gen_unified_ugwp CIRES UGWP Scheme General Algorithm
!! @{
!     subroutine unified_ugwp_run(me,  master, im,  levs, ntrac, dtp, fhzero, kdt,      &
!         lonr, oro, oro_uf, hprime, nmtvr, oc, theta, sigma, gamma, elvmax, clx, oa4,  &
!         varss,oc1ss,oa4ss,ol4ss,dx,                                                   &
!         dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,dusfc_ss,                                 &
!         dvsfc_ss,dusfc_fd,dvsfc_fd,dtaux2d_ls,dtauy2d_ls,dtaux2d_bl,dtauy2d_bl,       &
!         dtaux2d_ss,dtauy2d_ss,dtaux2d_fd,dtauy2d_fd,br1,hpbl,slmsk,                   &
!         do_tofd, ldiag_ugwp, cdmbgwd, jdat, xlat, xlat_d, sinlat, coslat, area,       &
!         ugrs, vgrs, tgrs, q1, prsi, prsl, prslk, phii, phil,                          &
!         del, kpbl, dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,                &
!         tau_tofd, tau_mtb, tau_ogw, tau_ngw, zmtb, zlwb, zogw,                        &
!         dudt_mtb,dudt_ogw, dudt_tms, du3dt_mtb, du3dt_ogw, du3dt_tms,                 &
!         dudt, dvdt, dtdt, rdxzb, con_g, con_omega, con_pi, con_cp, con_rd, con_rv,    &
!         con_rerth, con_fvirt, rain, ntke, q_tke, dqdt_tke, lprnt, ipr,                &
!         ldu3dt_ogw, ldv3dt_ogw, ldt3dt_ogw, ldu3dt_cgw, ldv3dt_cgw, ldt3dt_cgw,       &
!         ldiag3d, lssav, flag_for_gwd_generic_tend, do_ugwp_v0, do_ugwp_v0_orog_only,  &
!         do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1,              &
!         do_ugwp_v1_orog_only, gwd_opt, errmsg, errflg)

     subroutine unified_ugwp_run(me,  master, im,  levs, ntrac, lonr, dtp, fhzero, kdt, &  
          ldiag3d, lssav, flag_for_gwd_generic_tend, do_ugwp_v0, do_ugwp_v0_orog_only,  &
          do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd, do_ugwp_v1,              &
          do_ugwp_v1_orog_only, gwd_opt, do_tofd, ldiag_ugwp, cdmbgwd, jdat,            &
          con_g, con_omega, con_pi, con_cp, con_rd, con_rv, con_rerth, con_fvirt,       &
          nmtvr, oro, oro_uf, hprime, oc, theta, sigma, gamma, elvmax, clx, oa4,        &
          varss,oc1ss,oa4ss,ol4ss, dx,  xlat, xlat_d, sinlat, coslat, area,             &
          rain, br1, hpbl, kpbl, slmsk,                                                 &
          ugrs, vgrs, tgrs, q1, prsi, prsl, prslk, phii, phil,  del,                    &
	  jindx1_tau, jindx2_tau, ddy_j1tau, ddy_j2tau,                                 &
          dudt_ogw, dvdt_ogw, dtdt_sso,  du_ogwcol, dv_ogwcol,                          &
          dudt_obl, dvdt_obl, du_oblcol, dv_oblcol,                                     &
          dudt_oss, dvdt_oss, du_osscol, dv_osscol,                                     &
          dudt_ofd, dvdt_ofd, du_ofdcol, dv_ofdcol,                                     &
          dudt_ngw, dvdt_ngw, dtdt_ngw,   kdis_ngw,                                     &
          dudt_gw,  dvdt_gw,  dtdt_gw,    kdis_gw,   tau_ogw, tau_ngw,  tau_oss,        &
          zogw,  zlwb,  zobl,  zngw,   dusfcg, dvsfcg,  dudt, dvdt, dtdt, rdxzb,        &
          ldu3dt_ogw, ldv3dt_ogw, ldt3dt_ogw, ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw,       &
          lprnt, ipr,   ntke, q_tke,   dqdt_tke,  errmsg, errflg)
	  
!	  
! cap:               dt3dt(i,k) = dt3dt(i,k) - dtdt(i,k)*dtf
!
!
!########################################################################
!  Attention New Arrays and Names must be ADDED inside
!
!  a) /FV3/gfsphysics/GFS_layer/GFS_typedefs.meta
!  b) /FV3/gfsphysics/GFS_layer/GFS_typedefs.F90
!  c) /FV3/gfsphysics/GFS_layer/GFS_diagnostics.F90
!########################################################################
![ccpp-table-properties]
!  name = GFS_interstitial_type
!  type = ddt
!########################################################################
!
!
    implicit none

!   Preference use    (im,levs) rather than (:,:) to avoid memory-leaks
! order description control-logical 
! other in-variables
!       out-variables
!       local-variables
!       unified diagnostics inside CCPP and GFS_typedefs.F90/GFS_diagnostics.F90
!
!
! interface variables
    logical,                 intent(in) :: ldiag3d, lssav
    logical,                 intent(in) :: do_tofd, ldiag_ugwp
    logical,                 intent(in) :: flag_for_gwd_generic_tend  
   ! flags for choosing combination of GW drag schemes to run
    logical,                intent (in) :: do_ugwp_v0, do_ugwp_v0_orog_only,  &
                                           do_gsl_drag_ls_bl, do_gsl_drag_ss, &
                                           do_gsl_drag_tofd, do_ugwp_v1,      &
                                           do_ugwp_v1_orog_only

    logical,                 intent(in) :: lprnt
    integer,                 intent(in) :: ipr
    integer,                 intent(in) :: gwd_opt


    integer,                 intent(in) :: me, master, im, levs, ntrac, kdt, lonr
    real(kind=kind_phys),    intent(in) :: dtp, fhzero
    integer,                 intent(in) :: jdat(8)

! SSO parameters and variables

    integer,                 intent(in) :: nmtvr
    real(kind=kind_phys),    intent(in) :: cdmbgwd(4)
    
    real(kind=kind_phys),    intent(in), dimension(im)       :: oro, oro_uf
    real(kind=kind_phys),    intent(in), dimension(im)       :: hprime, oc, theta, sigma, gamma
!=====
! "semi-bug" in gwdps.f
! only in  must edit gwdps.f with passing "elvmax_in" and modifying
!                   ELVMAX --- local
! line: ELVMAX(J) = min (ELVMAX_IN(J) + sigfac * hprime(j), hncrit)
!=====
    real(kind=kind_phys),    intent(inout), dimension(im)       :: elvmax
    real(kind=kind_phys),    intent(in), dimension(im, 4)    :: clx, oa4

    real(kind=kind_phys),    intent(in), dimension(im)       :: varss,oc1ss,oa4ss,ol4ss,dx

!
!ccpp-style passing constants
!
    real(kind=kind_phys),    intent(in) :: con_g, con_omega, con_pi, con_cp, con_rd, &
                                           con_rv, con_rerth, con_fvirt
! grids 

    real(kind=kind_phys),    intent(in), dimension(im)       :: xlat, xlat_d, sinlat, coslat, area

! State vars + PBL/slmsk +rain 

    real(kind=kind_phys),    intent(in), dimension(im, levs)   :: del, ugrs, vgrs, tgrs, prsl, prslk, phil
    real(kind=kind_phys),    intent(in), dimension(im, levs+1) :: prsi, phii
    real(kind=kind_phys),    intent(in), dimension(im, levs)   :: q1
    integer,                 intent(in), dimension(im)         :: kpbl
    
    real(kind=kind_phys),    intent(in), dimension(im) :: rain
    real(kind=kind_phys),    intent(in), dimension(im) :: br1, hpbl,  slmsk   
    real(kind=kind_phys),    intent(in), dimension(im) :: ddy_j1tau, ddy_j2tau
    integer,                 intent(in), dimension(im) :: jindx1_tau, jindx2_tau 
!
! Arrays that proposed EMC without the proof of concept: GW excitation in PBL by "tke-turbulence"???
!  q_tke, dqdt_tke, ntke
!
    integer,                 intent(in) :: ntke
    real(kind=kind_phys),    intent(in), dimension(im, levs) :: q_tke, dqdt_tke
!Output (optional):

    real(kind=kind_phys), intent(out), dimension(im)  ::                  &
                            du_ogwcol,  dv_ogwcol,  du_oblcol, dv_oblcol, &       
                            du_osscol,  dv_osscol,  du_ofdcol, dv_ofdcol  
!
! we may add later but due to launch in the upper layes ~ mPa comparing to ORO Pa*(0.1)			    
!                            du_ngwcol, dv_ngwcol

    real(kind=kind_phys), intent(out), dimension(im)  :: dusfcg, dvsfcg
    real(kind=kind_phys), intent(out), dimension(im)  :: tau_ogw, tau_ngw, tau_oss

    real(kind=kind_phys), intent(out) , dimension(im, levs) ::    &
                          dudt_ogw, dvdt_ogw, dudt_obl, dvdt_obl, &
                          dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd

    real(kind=kind_phys), intent(out) , dimension(im, levs) ::              &
                          dudt_ngw, dvdt_ngw, kdis_ngw,                     &
                          dudt_gw,  dvdt_gw, kdis_gw

    real(kind=kind_phys), intent(out) , dimension(im) ::  zogw,  zlwb,  zobl,  zngw
!  GFS-style of tendencies (accumulated or chain) = dudt = dudt_prev + dudt_cur
!
    real(kind=kind_phys), intent(out) , dimension(im, levs) :: dtdt_sso, dtdt_ngw, dtdt_gw
    real(kind=kind_phys), intent(inout), dimension(im, levs):: dudt, dvdt, dtdt

!
! These arrays are only allocated if ldiag=.true.
!
! Version of COORDE updated by CCPP-dev for time-aver
!
    real(kind=kind_phys),    intent(inout), dimension(:,:)   :: ldu3dt_ogw, ldv3dt_ogw, ldt3dt_ogw
    real(kind=kind_phys),    intent(inout), dimension(:,:)   :: ldu3dt_ngw, ldv3dt_ngw, ldt3dt_ngw
 


    real(kind=kind_phys),    intent(out), dimension(im)      :: rdxzb     ! for stoch phys. mtb-level

    character(len=*),        intent(out) :: errmsg
    integer,                 intent(out) :: errflg

! local variables
    integer :: i, k
    real(kind=kind_phys), dimension(im)       :: sgh30
    real(kind=kind_phys), dimension(im, levs) :: Pdvdt, Pdudt
    real(kind=kind_phys), dimension(im, levs) :: Pdtdt, Pkdis
!------------
!
! from ugwp_driver_v0.f -> cires_ugwp_initialize.F90 -> module ugwp_wmsdis_init
!  now in the namelist of cires_ugwp "knob_ugwp_tauamp" controls tamp_mpa 
!
!        tamp_mpa =knob_ugwp_tauamp           !amplitude for GEOS-5/MERRA-2
!------------
    real(kind=kind_phys), parameter :: tamp_mpa_v0=30.e-3  ! large flux to help "GFS-ensembles" in July 2019

! switches that activate impact of OGWs and NGWs
    
    real(kind=kind_phys), parameter :: fw1_tau=1.0

    integer :: nmtvr_temp

    real(kind=kind_phys) :: inv_g
    
    real(kind=kind_phys), dimension(im, levs)   :: zmet  ! geopotential height at model Layer centers
    real(kind=kind_phys), dimension(im, levs+1) :: zmeti ! geopotential height at model layer interfaces


! ugwp_v1 local variables
    integer :: y4, month, day,  ddd_ugwp, curdate, curday
    integer :: hour
    real(kind=kind_phys) :: hcurdate, hcurday, fhour, fhrday
    integer :: kdtrest
    integer :: curday_ugwp
    integer :: curday_save=20150101
    logical :: first_qbo=.true.
    real    ::  hcurday_save =20150101.00

    save curday_save, hcurday_save


!  ugwp_v1 temporary (local) diagnostic variables from cires_ugwp_solv2_v1
!  diagnostics for wind and temp rms to compare with space-borne data and metrics
!   in the Middle atmosphere: 20-110 km ( not active in CCPP-style, oct 2020)
!    real(kind=kind_phys) :: tauabs(im,levs), wrms(im,levs), trms(im,levs)


    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

! 1) ORO stationary GWs
!    ------------------
!    
! for all oro-suites can uze geo-meters having "hpbl"
!    
       inv_g = 1./con_g
!
! All GW-schemes operate with Zmet =phil*inv_g, passing Zmet/Zmeti is more robust
!       
       zmeti  = phii*inv_g
       zmet   = phil*inv_g
       
       
!===============================================================
! ORO-diag
		 
      dudt_ogw  = 0. ; dvdt_ogw=0. ; dudt_obl=0. ; dvdt_obl=0.            
      dudt_oss  = 0. ; dvdt_oss=0. ; dudt_ofd=0. ; dvdt_ofd=0.  
              
      dusfcg   = 0.  ;  dvsfcg =0.    
                                     
      du_ogwcol=0. ; dv_ogwcol=0. ; du_oblcol=0. ; dv_oblcol=0. 
      du_osscol=0. ; dv_osscol=0. ;du_ofdcol=0.  ; dv_ofdcol=0. 
      
! ngw-diag
      
      dudt_ngw=0. ; dvdt_ngw=0. ; dtdt_ngw=0. ; kdis_ngw=0.  
      
! ngw+ogw - diag            
      
      dudt_gw=0. ;  dvdt_gw=0. ;  dtdt_gw=0. ;  kdis_gw=0. 
	  
! source fluxes
	  	    
      tau_ogw(:)=0. ; tau_ngw(:)=0. ;  tau_oss(:)=0.  
      
! launch layers
              
      zlwb(:)= 0.  ; zogw(:)=0. ;  zobl(:)=0. ;  zngw(:)=0.
!===============================================================
!  Accumulated tendencies due to 3-SSO schemes (all ORO-physics)
!  ogw + obl +oss +ofd ..... no explicit "lee wave trapping"
!===============================================================
     do k=1,levs
        do i=1,im
          Pdvdt(i,k) = 0.0
          Pdudt(i,k) = 0.0
          Pdtdt(i,k) = 0.0
          Pkdis(i,k) = 0.0
        enddo
      enddo
      
!   ------------------
!
!  Also zero all ORO diag-c arrays to avoid "special ifs and zeros"
!       like old GFS-ORO gwdps_run has limited diagnostics
!
!   ------------------

    ! Run the appropriate large-scale (large-scale GWD + blocking) scheme
    ! Note:  In case of GSL drag_suite, this includes ss and tofd

    if ( do_gsl_drag_ls_bl.or.do_gsl_drag_ss.or.do_gsl_drag_tofd ) then
!
! the zero value assigned inside "drag_suite_run":
! dudt_ogw, dvdt_ogw, dudt_obl, dvdt_obl,dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd
!du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol, du_osscol, dv_osscol, du_ofdcol dv_ofdcol
! dusfcg,  dvsfcg
!                    gsd_diss_ht_opt =0 => Pdtdt = bl+ls +(Pdtdt=0)
!
       call drag_suite_run(im,levs, Pdvdt, Pdudt, Pdtdt,             &
                 ugrs,vgrs,tgrs,q1,                                  &
                 kpbl,prsi,del,prsl,prslk,phii,phil,dtp,             &
                 kdt,hprime,oc,oa4,clx,varss,oc1ss,oa4ss,            &
                 ol4ss,theta,sigma,gamma,elvmax,                     &
                  dudt_ogw, dvdt_ogw, dudt_obl, dvdt_obl,            &
                  dudt_oss, dvdt_oss, dudt_ofd, dvdt_ofd,            &
                  dusfcg,  dvsfcg,                                   &
                  du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol,        &
                  du_osscol, dv_osscol, du_ofdcol, dv_ofdcol,         &
                 slmsk,br1,hpbl,con_g,con_cp,con_rd,con_rv,          &
                 con_fvirt,con_pi,lonr,                              &
                 cdmbgwd(1:2),me,master,lprnt,ipr,rdxzb,dx,gwd_opt,  &
                 do_gsl_drag_ls_bl,do_gsl_drag_ss,do_gsl_drag_tofd,  &
                 errmsg,errflg)
! 		 
! dusfcg = du_ogwcol + du_oblcol + du_osscol + du_ofdcol
    end if

    if ( do_ugwp_v1.or.do_ugwp_v1_orog_only ) then

       ! Valery's TOFD
       ! topo paras
                                        ! w/ orographic effects
       if(nmtvr == 14)then
                                        ! calculate sgh30 for TOFD
!         sgh30 = abs(oro - oro_uf)
	 sgh30 = varss                  ! small-scale "turbulent" oro
                                        ! w/o orographic effects
       else
         sgh30   = varss
       endif
!
! only sum of integrated ORO+GW effects (dusfcg and dvsfcg) = sum(ogw + obl + oss*0 + ofd + ngw)
!
! OROGW_V1 suggests "orchestration" between OGW-effects and Mountain Blocking
!      it provides options for the Scale-Aware (SA)formulation of SSO-effects
!
       call orogw_v1 (im, levs,  lonr,  me, master,dtp, kdt, do_tofd,     &
                      con_g, con_omega, con_rd, con_cp, con_rv,con_pi,    &
                      con_rerth, con_fvirt,xlat_d, sinlat, coslat, area,  &
                      cdmbgwd(1:2), hprime, oc, oa4, clx, theta,          &
                      sigma, gamma, elvmax,  sgh30,  kpbl,                &		     
                ugrs ,vgrs, tgrs, q1, prsi,del,prsl,prslk, zmeti, zmet,   &		      	      
                      Pdvdt, Pdudt, Pdtdt, Pkdis, DUSFCg, DVSFCg,rdxzb,   &		      
                      zobl, zlwb, zogw, tau_ogw, dudt_ogw, dvdt_ogw,      &
                      dudt_obl, dvdt_obl,dudt_ofd, dvdt_ofd,              &
                      du_ogwcol, dv_ogwcol, du_oblcol, dv_oblcol,         &
                      du_ofdcol, dv_ofdcol)
!
! dusfcg = du_ogwcol + du_oblcol  + du_ofdcol                                  only 3 terms
!
! add complete diag like WRF: dv_ofdcol, dvdt_obl, dvdt_ogw, dvdt_ofd
!

    end if

    if ( do_ugwp_v0.or.do_ugwp_v0_orog_only ) then


      if (cdmbgwd(1) > 0.0 .or. cdmbgwd(2) > 0.0) then

! Override nmtvr with nmtvr_temp = 14 for passing into gwdps_run if necessary
        if ( nmtvr == 24 ) then           ! gwd_opt = 2, 22, 3, or 33
           nmtvr_temp = 14
        else
           nmtvr_temp = nmtvr
        end if
!============================================================================
! we may add full diagnostics of GFSv16 operational scheme:  and Fix elvmax
!
! only MB + OGW:   Pdvdt = dvdt_ogw + dvdt_obl .....
!                  can be passed as separate terms
!        adding    zogw + zobl,  tau_ogw, 
!============================================================================
        call gwdps_run(im, levs, Pdvdt, Pdudt, Pdtdt,                  &
                   ugrs, vgrs, tgrs, q1,                               &
                   kpbl, prsi, del, prsl, prslk, phii, phil, dtp, kdt, &
                   hprime, oc, oa4, clx, theta, sigma, gamma,          &
                   elvmax, dusfcg, dvsfcg,                             &
                   con_g,  con_cp, con_rd, con_rv, lonr,               &
                   nmtvr_temp, cdmbgwd, me, lprnt, ipr, rdxzb,         &
                   errmsg, errflg)
        if (errflg/=0) return
      endif
  
    end if
!
! Assign all SSO-tendenices   Pdvdt, Pdudt, Pdtdt and Pkdis
!        
!
!  GFS-style diag dt3dt(:.:, 1:14)
!
     if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
        do k=1,levs
          do i=1,im
             ldu3dt_ogw(i,k) = ldu3dt_ogw(i,k) + Pdudt(i,k)*dtp
             ldv3dt_ogw(i,k) = ldv3dt_ogw(i,k) + Pdvdt(i,k)*dtp
             ldt3dt_ogw(i,k) = ldt3dt_ogw(i,k) + Pdtdt(i,k)*dtp
          enddo
        enddo
      endif


      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin non-stationary GW schemes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! ugwp_v0 non-stationary GW drag
    !
    if (do_ugwp_v0) then

      if (cdmbgwd(3) > 0.0) then

! 2) non-stationary GW-scheme with GMAO/MERRA GW-forcing

        call slat_geos5_tamp_v0(im, tamp_mpa_v0, xlat_d, tau_ngw)
!
! below the EMC-impact without rigorous tests reported elsewhere       
!
      if (cdmbgwd(3) .lt. 1.0 ) &
         call emc_modulation(im , levs,  ntke, tau_ngw, cdmbgwd(3), cdmbgwd(4), dtp,  &
                            q_tke, dqdt_tke, del, rain)

        call fv3_ugwp_solv2_v0(im, levs, dtp, tgrs, ugrs, vgrs, q1,                            &
             prsl, prsi, phil, xlat_d, sinlat, coslat, dudt_ngw, dvdt_ngw, dtdt_ngw, kdis_ngw, &
             tau_ngw, me, master, kdt)


      else  ! .not.(cdmbgwd(3) > 0.0)
	  
            dtdt_ngw(:,:) = 0.
            dudt_ngw(:,:) = 0.
            dvdt_ngw(:,:) = 0.
            kdis_ngw(:,:) =  0.

      endif  ! cdmbgwd(3) > 0.0


    end if  ! do_ugwp_v0 


!
! ugwp_v1 non-stationary GW drag
!
    if (do_ugwp_v1) then

       call slat_geos5_tamp_v1(im, tamp_mpa, xlat_d, tau_ngw)       
       call    slat_geos5_2020(im, tamp_mpa, xlat_d, tau_ngw)
       
       y4 = jdat(1); month = jdat(2); day = jdat(3) ; hour = jdat(5)
       
! fhour = float(hour)+float(jdat(6))/60. + float(jdat(7))/3600.
!       fhour = (kdt-1)*dtp/3600.
!       fhrday  = fhour/24.  - nint(fhour/24.)
     
                
       call calendar_ugwp(y4, month, day, ddd_ugwp)       
       curdate = y4*1000 + ddd_ugwp
       
       call ngwflux_update(me, master, im, levs, kdt, ddd_ugwp,curdate, &
            jindx1_tau, jindx2_tau, ddy_j1tau, ddy_j2tau,               &
	    xlat_d, sinlat,coslat, rain, tau_ngw)
                      
       call cires_ugwp_solv2_v1(me, master, im,   levs,  kdt, dtp,      &
                      tau_ngw, tgrs, ugrs,  vgrs,   q1, prsl, prsi,     &
                      zmet, zmeti,prslk,   xlat_d, sinlat, coslat,      &
        con_g, con_cp, con_rd, con_rv, con_omega,  con_pi, con_fvirt,   &
                      dudt_ngw, dvdt_ngw, dtdt_ngw, kdis_ngw, zngw)
		      

       if (me == master .and. kdt < 2) then
         print *
         write(6,*)'FV3GFS finished fv3_ugwp_solv2_v1  in ugwp_driver_v0 '
         write(6,*) ' non-stationary GWs with GMAO/MERRA GW-forcing '
         print *
       endif

     
    end if   ! do_ugwp_v1
    
!
!  GFS-style diag dt3dt(:.:, 1:14)  time-averaged
!   
      if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
        do k=1,levs
          do i=1,im
             ldu3dt_ngw(i,k) = ldu3dt_ngw(i,k) + dudt_ngw(i,k)*dtp
             ldv3dt_ngw(i,k) = ldv3dt_ngw(i,k) + dvdt_ngw(i,k)*dtp
             ldt3dt_ngw(i,k) = ldt3dt_ngw(i,k) + dtdt_ngw(i,k)*dtp
          enddo
        enddo
      endif
         
!
! get total sso-OGW + NGW
!
     dudt_gw =  Pdudt +dudt_ngw
     dvdt_gw =  Pdvdt +dvdt_ngw
     dtdt_gw =  Pdtdt +dtdt_ngw  
     kdis_gw =  Pkdis +kdis_ngw               
! 
! add to previous phys-tendencies
! Weird accumulation of GFS ( pbl+gw+rf )  

     dudt  = dudt  + dudt_gw
     dvdt  = dvdt  + dvdt_gw
     dtdt  = dtdt  + dvdt_gw
!=======================================================     
! GFS_GWD_generic_post_run   in   GFS_GWD_generic.F90
! ( five diag-averaged arrays)
!        dugwd(:) = dugwd(:) + dusfcg(:)*dtf       time-averaged  where ? no entries in ugwp_post.F90
!        dvgwd(:) = dvgwd(:) + dvsfcg(:)*dtf
!
!        if (ldiag3d .and. flag_for_gwd_generic_tend) then
!          du3dt(:,:) = du3dt(:,:) + dudt(:,:) * dtf
!          dv3dt(:,:) = dv3dt(:,:) + dvdt(:,:) * dtf
!          dt3dt(:,:) = dt3dt(:,:) + dtdt(:,:) * dtf
!        endif
!===========================================   
    end subroutine unified_ugwp_run
!! @}
!>@}
end module unified_ugwp

!==============================================================
!
! cap auto_code sequence:     1) call GFS_GWD_generic_pre_run(im=
!                             2) call unified_ugwp_run     
!                             3) call unified_ugwp_post_run(
!                             4) call GFS_GWD_generic_post_run(
!                             5) call rayleigh_damp_run
!                             6) call GFS_suite_stateout_update_run
!
!      gt0(:,:)   = tgrs(:,:)   + dtdt(:,:)   * dtp
!      gu0(:,:)   = ugrs(:,:)   + dudt(:,:)   * dtp
!      gv0(:,:)   = vgrs(:,:)   + dvdt(:,:)   * dtp
!      gq0(:,:,:) = qgrs(:,:,:) + dqdt(:,:,:) * dtp
!==============================================================
!
!
!
! extra operations
! cap:               dt3dt(i,k) = dt3dt(i,k) - dtdt(i,k)*dtf
!types: GFS_Control%lssav
!       GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,2)
!       GFS_Interstitial(cdata%thrd_no)%dtdt
!
!                  lssav=GFS_Control%lssav,ldiag3d=GFS_Control%ldiag3d,
!                  dudt=GFS_Interstitial(cdata%thrd_no)%dudt, &
!                  dvdt=GFS_Interstitial(cdata%thrd_no)%dvdt,
!                  dtdt=GFS_Interstitial(cdata%thrd_no)%dtdt, &
!
!                  du3dt=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,2),
!                  dv3dt=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,2), &
!                  dt3dt=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,7)
!
!ldu3dt_ogw=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,2), &
!ldv3dt_ogw=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,2),
!ldt3dt_ogw=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,7), &
!ldu3dt_cgw=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,4),
!ldv3dt_cgw=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,4), &
!ldt3dt_cgw=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,9)
!ldiag3d=GFS_Control%ldiag3d, &
!
! diagnostics instant vs time-aver
! dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,dusfc_ss,  dvsfc_ss,  usfc_fd,dvsfc_fd 
! dtaux2d_ls,dtauy2d_ls, dtaux2d_bl,dtauy2d_bl,  dtaux2d_ss,dtauy2d_ss,   dtaux2d_fd,dtauy2d_fd
!      
!                  dusfc_ls=GFS_Data(cdata%blk_no)%Intdiag%dusfc_ls,dvsfc_ls=GFS_Data(cdata%blk_no)%Intdiag%dvsfc_ls, &
!                  dusfc_bl=GFS_Data(cdata%blk_no)%Intdiag%dusfc_bl,dvsfc_bl=GFS_Data(cdata%blk_no)%Intdiag%dvsfc_bl, &
!                  dusfc_ss=GFS_Data(cdata%blk_no)%Intdiag%dusfc_ss,dvsfc_ss=GFS_Data(cdata%blk_no)%Intdiag%dvsfc_ss, &
!                  dusfc_fd=GFS_Data(cdata%blk_no)%Intdiag%dusfc_fd,dvsfc_fd=GFS_Data(cdata%blk_no)%Intdiag%dvsfc_fd, &
!                  !dtaux2d_ls=GFS_Data(cdata%blk_no)%Intdiag%dtaux2d_ls,dtauy2d_ls=GFS_Data(cdata%blk_no)%Intdiag%dtauy2d_ls, &
!                  !dtaux2d_bl=GFS_Data(cdata%blk_no)%Intdiag%dtaux2d_bl,dtauy2d_bl=GFS_Data(cdata%blk_no)%Intdiag%dtauy2d_bl, &
!                  !dtaux2d_ss=GFS_Data(cdata%blk_no)%Intdiag%dtaux2d_ss,dtauy2d_ss=GFS_Data(cdata%blk_no)%Intdiag%dtauy2d_ss, &
!                  !dtaux2d_fd=GFS_Data(cdata%blk_no)%Intdiag%dtaux2d_fd,dtauy2d_fd=GFS_Data(cdata%blk_no)%Intdiag%dtauy2d_fd, &
!
! dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,                
! tau_tofd, tau_mtb, tau_ogw, tau_ngw, zmtb, zlwb, zogw,                        
! dudt_mtb,dudt_ogw, dudt_tms, du3dt_mtb, du3dt_ogw, du3dt_tms,     
!
!           
!                  dusfcg=GFS_Interstitial(cdata%thrd_no)%dusfcg,dvsfcg=GFS_Interstitial(cdata%thrd_no)%dvsfcg, &
!                  gw_dudt=GFS_Interstitial(cdata%thrd_no)%gw_dudt,gw_dvdt=GFS_Interstitial(cdata%thrd_no)%gw_dvdt, &
!                  gw_dtdt=GFS_Interstitial(cdata%thrd_no)%gw_dtdt,gw_kdis=GFS_Interstitial(cdata%thrd_no)%gw_kdis, &
!                  tau_tofd=GFS_Interstitial(cdata%thrd_no)%tau_tofd,tau_mtb=GFS_Interstitial(cdata%thrd_no)%tau_mtb, &
!                  tau_ogw=GFS_Interstitial(cdata%thrd_no)%tau_ogw,tau_ngw=GFS_Interstitial(cdata%thrd_no)%tau_ngw, &
!                  zmtb=GFS_Interstitial(cdata%thrd_no)%zmtb,zlwb=GFS_Interstitial(cdata%thrd_no)%zlwb, &
!                  zogw=GFS_Interstitial(cdata%thrd_no)%zogw,dudt_mtb=GFS_Interstitial(cdata%thrd_no)%dudt_mtb, &
!                  dudt_ogw=GFS_Interstitial(cdata%thrd_no)%dudt_ogw,dudt_tms=GFS_Interstitial(cdata%thrd_no)%dudt_tms, &
!
!                  du3dt_mtb=GFS_Data(cdata%blk_no)%Intdiag%du3dt_mtb,
!                  du3dt_ogw=GFS_Data(cdata%blk_no)%Intdiag%du3dt_ogw, 
!                  du3dt_tms=GFS_Data(cdata%blk_no)%Intdiag%du3dt_tms,
! dudt, dvdt, dtdt
!                  dudt=GFS_Interstitial(cdata%thrd_no)%dudt, 
!                  dvdt=GFS_Interstitial(cdata%thrd_no)%dvdt,
!                  dtdt=GFS_Interstitial(cdata%thrd_no)%dtdt, 
!
!ldu3dt_ogw, ldv3dt_ogw, ldt3dt_ogw, ldu3dt_cgw, ldv3dt_cgw, ldt3dt_cgw =============> OLD-style-GFS/GSM diagnostics
!                  ldu3dt_ogw=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,2), &
!                  ldv3dt_ogw=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,2),
!                  ldt3dt_ogw=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,7), &
!                  ldu3dt_cgw=GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,4),
!                  ldv3dt_cgw=GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,4), &
!                  ldt3dt_cgw=GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,9)
!=============================================
