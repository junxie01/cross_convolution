! use cross concolutio to measure the time shift
program cross_conv_measure
use sacio
implicit none
type(sac_head) :: sachead(4),sacheado
integer:: ncross
integer:: nn,nk,npt,iflag,shft,id
integer i,j,npts,nerr,ncros,nb,ne,n1,n2
character(180):: file(4),outpt,par
integer,parameter:: nptsmax=1000000
real(4):: f12(nptsmax),f21(nptsmax)
real(4):: sig(nptsmax,4),dt,beg,inta,intb
real(4):: sig1_ev1(nptsmax),sig1_ev2(nptsmax)
real(4):: sig2_ev1(nptsmax),sig2_ev2(nptsmax)
real(4):: cc(nptsmax),amax,aa(nptsmax),bb(nptsmax)
if(iargc().ne.4)stop 'Usage: cross_conv_measure sta1_ev1 sta1_ev2 sta2_ev1 sta2_ev2'
do i=1,4
   call getarg(i,file(i))
enddo
!call getarg(5,par)
!read(par,'(bn,i20)')id
do i=1,4
   call read_sachead(file(i),sachead(i),nerr)
   if(nerr.eq.-1)stop 'Read file error'
   call read_sac(file(i),sig(:,i),sachead(i),nerr)
   if(nerr.eq.-1)stop 'Read file error'
enddo
npts=sachead(1)%npts
dt=sachead(1)%delta
call initial_sachead(sacheado)
sacheado%delta=dt
call cal_crstm(sig(:,1),sig(:,2),sig(:,3),sig(:,4),sachead(1)%npts,f12,f21,npt)
sacheado%npts=npt
sacheado%b=0
shft=npt-1
sacheado%t1=0
!call write_sac('f12.sac',f12,sacheado,nerr)
!call write_sac('f21.sac',f21,sacheado,nerr)
!call do_ncc(npt,f12,f21,cc,shft,dt)
ncros=npt-sachead(1)%npts-2
nb=npt/2-npts/2
ne=nb+npts-1
aa(1:npts)=f21(nb:ne)
cc=0
do i=1,ncros
   n1=i
   n2=npts+i-1
   bb(1:npts)=f12(n1:n2)        
   call cross(aa,bb,npts,cc(i))
enddo
!call cross_td(f12,f21,cc,npt,shft)
!beg=-shft*sacheado%delta
beg=-nb*sachead(1)%delta
sacheado%b=beg
sacheado%npts=2*shft+1
!call write_sac('test.sac',cc,sacheado,nerr)
write(*,*)(maxloc(cc)-1)*sacheado%delta+beg,maxval(cc)
end program
