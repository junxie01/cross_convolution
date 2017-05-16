subroutine conv(sig1,sig2,npts,dt,sig)
implicit none
integer:: nptsmax
parameter(nptsmax=4000000)
integer:: i,npts,npow,nn
real(4):: dt
real(4):: sig1(nptsmax),sig2(nptsmax),sig(nptsmax)
complex(4) :: sig1c(nptsmax),sig2c(nptsmax)
nn=2;npow=1
do while(nn.le.npts)
   nn=nn*2;npow=npow+1
enddo
nk=nn/2+1
sig1c=cmplx(0.0,0.0);sig2c=cmplx(0.0,0.0)
do i=1,npts
   sig1c(npts-i+1)=cmplx(sig1(i),0.0)
enddo
sig2c(1:npts)=cmplx(sig2(1:npts),0.0)
call clogc(npow,sig1c,1,dt)
call clogc(npow,sig2c,1,dt)
call cor(sig1c,sig2c,nft,nk,npow,dt)
end subroutine
