!>  \file unified_ugwp_post.F90
!! This file contains
module unified_ugwp_post

contains

!>\defgroup unified_ugwp_post unified_UGWP Scheme Post
!! @{
!> \section arg_table_unified_ugwp_post_init Argument Table
!!
    subroutine unified_ugwp_post_init ()
    end subroutine unified_ugwp_post_init

!>@brief The subroutine initializes the unified UGWP

!> \section arg_table_unified_ugwp_post_run Argument Table
!! \htmlinclude unified_ugwp_post_run.html
!!



     subroutine unified_ugwp_post_run (ldiag_ugwp, dtf, im, levs,     &
         dtdt_gw, dudt_gw, dvdt_gw, du_ofdcol, du_oblcol, tau_ogw,    &
         tau_ngw, zobl, zlwb, zogw, dudt_obl, dudt_ogw, dudt_ofd,     &
         tot_zmtb, tot_zlwb, tot_zogw,                                &
         tot_tofd, tot_mtb, tot_ogw, tot_ngw,                         &
         du3dt_mtb,du3dt_ogw, du3dt_tms, du3dt_ngw, dv3dt_ngw,        &
         dtdt, dudt, dvdt, errmsg, errflg)

        use machine,                only: kind_phys

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: ldiag_ugwp      !< flag for CIRES UGWP Diagnostics

        real(kind=kind_phys), intent(in),    dimension(:)   :: zobl, zlwb, zogw
        real(kind=kind_phys), intent(in),    dimension(:)   :: du_ofdcol, tau_ogw, du_oblcol, tau_ngw
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_mtb, tot_ogw, tot_tofd, tot_ngw
        real(kind=kind_phys), intent(inout), dimension(:)   :: tot_zmtb, tot_zlwb, tot_zogw
        real(kind=kind_phys), intent(in),    dimension(:,:) :: dtdt_gw, dudt_gw, dvdt_gw
        real(kind=kind_phys), intent(in),    dimension(:,:) :: dudt_obl, dudt_ogw, dudt_ofd
        real(kind=kind_phys), intent(inout), dimension(:,:) :: du3dt_mtb, du3dt_ogw, du3dt_tms, du3dt_ngw, dv3dt_ngw
	
        real(kind=kind_phys), intent(inout), dimension(:,:) :: dtdt, dudt, dvdt

        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0
!
! post creates the "time-averaged" diagnostics"
!

        if (ldiag_ugwp) then
          tot_zmtb =  tot_zmtb + dtf *zobl
          tot_zlwb =  tot_zlwb + dtf *zlwb
          tot_zogw =  tot_zogw + dtf *zogw
    
          tot_tofd  = tot_tofd + dtf *du_ofdcol
          tot_mtb   = tot_mtb +  dtf *du_oblcol
          tot_ogw   = tot_ogw +  dtf *tau_ogw
          tot_ngw   = tot_ngw +  dtf *tau_ngw
    
          du3dt_mtb = du3dt_mtb + dtf *dudt_obl
          du3dt_tms = du3dt_tms + dtf *dudt_ofd
          du3dt_ogw = du3dt_ogw + dtf *dudt_ogw
          du3dt_ngw = du3dt_ngw + dtf *dudt_gw
          dv3dt_ngw = dv3dt_ngw + dtf *dvdt_gw
        endif
	
!=====================================================================
! Updates inside the unified_ugwp.F90
!
!        dtdt = dtdt + gw_dtdt
!        dudt = dudt + gw_dudt
!        dvdt = dvdt + gw_dvdt
!
!       "post" may  also create the "time-averaged" diagnostics"
!            
!     if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
!        do k=1,levs
!          do i=1,im
!             ldu3dt_ngw(i,k) = ldu3dt_ngw(i,k) + dudt_ngw(i,k)*dtf
!             ldv3dt_ngw(i,k) = ldv3dt_ngw(i,k) + dvdt_ngw(i,k)*dtf
!             ldt3dt_ngw(i,k) = ldt3dt_ngw(i,k) + dtdt_ngw(i,k)*dtf
!	  
!             ldu3dt_ogw(i,k) = ldu3dt_ogw(i,k) + dudt_ogw(i,k)*dtf
!             ldv3dt_ogw(i,k) = ldv3dt_ogw(i,k) + dvdt_ogw(i,k)*dtf
!             ldt3dt_ogw(i,k) = ldt3dt_ogw(i,k) + dtdt_ogw(i,k)*dtf
!          enddo
!        enddo
!      endif
! 
!=====================================================================
      end subroutine unified_ugwp_post_run

!> \section arg_table_unified_ugwp_post_finalize Argument Table
!!
      subroutine unified_ugwp_post_finalize ()
      end subroutine unified_ugwp_post_finalize

!! @}
end module unified_ugwp_post
