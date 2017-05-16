subroutine cal_crstm(sig1_ev1,sig1_ev2,sig2_ev1,sig2_ev2,npts,f12,f21,npt)
implicit none
integer :: npts,npt,nptsmax
parameter(nptsmax=1000000)
real(4) :: sig1_ev1(nptsmax),sig2_ev2(nptsmax)
real(4) :: sig1_ev2(nptsmax),sig2_ev1(nptsmax)
real(4) :: f12(nptsmax),f21(nptsmax)
!write(*,*)'npts=',npts,'npts=',npts
call do_conv(sig1_ev1,sig2_ev2,f12,npts,npts,npt) ! fft
call do_conv(sig1_ev2,sig2_ev1,f21,npts,npts,npt) ! fft
return
end subroutine
