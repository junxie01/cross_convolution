! npts is the length of sig1 and sig2
subroutine do_conv(sig1,sig2,sig_out,npts1,npts2,shft,id)
implicit none
integer nptsmax
parameter (nptsmax=1000000)
integer shft,nk,id
real sig1(nptsmax),sig2(nptsmax),dt
real sig_out(nptsmax)
real(8) ::df
complex(4),dimension(:),allocatable :: ss1,ss2
integer nft,i,npts1,npow,nft2,npts2
nft=2;npow=1
shft=npts1+npts2-1
do while(nft.lt.shft)
   nft=nft*2;npow=npow+1
enddo
dt=1
!df=dble(1.0/nft/dt)
nk=nft/2+1
allocate(ss1(nft),ss2(nft))
ss1=cmplx(0,0);ss2=cmplx(0,0)
ss1(1:npts1)=cmplx(sig1(1:npts1),0)
ss2(1:npts2)=cmplx(sig2(1:npts2),0)
! do fft
call clogc(npow,ss1,1,dt)
call clogc(npow,ss2,1,dt)
! do convolution 
ss1=ss1*ss2
! do ifft
call clogc(npow,ss1,-1,dt)
sig_out(1:shft)=real(ss1(1:shft))
deallocate(ss1,ss2)
return
end subroutine
