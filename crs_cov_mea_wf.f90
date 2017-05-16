! use cross concolutio to measure the time shift
program cross_conv_measure
use sacio
implicit none
type(sac_head) :: sachead(4),sacheado
integer:: ncross
integer:: nn,nk,npt,iflag,shft,nid(1)
integer i,j,npts,nerr,ncros,nb,ne,n1,n2
character(180):: file(4),outpt,para
integer,parameter:: nptsmax=1000000
real(4):: tb,te,f1,f2,fre(4)
real(4):: f12(nptsmax),f21(nptsmax)
real(4):: sig(nptsmax,4),dt,beg,inta,intb
real(4):: sigo(nptsmax,4),sigf(nptsmax,4)
real(4):: sig1_ev1(nptsmax),sig1_ev2(nptsmax)
real(4):: sig2_ev1(nptsmax),sig2_ev2(nptsmax)
real(4):: cc(nptsmax),amax,aa(nptsmax),bb(nptsmax)
logical:: ext
if(iargc().ne.1)then
   write(*,*) 'Usage: cross_conv_measure param.dat '
   write(*,*) 'param.dat:'
   write(*,*) 'sta1_ev1 sta1_ev2 sta2_ev1 sta2_ev2'
!   write(*,*) 'f1 f2'
   write(*,*) 'tb te'
   stop
endif
call getarg(1,para)
inquire(file=trim(para),exist=ext)
if(.not.ext)stop 'param.dat file does not exist'
open(11,file=para)
read(11,*)(file(i),i=1,4)
!read(11,*)fre(2),fre(3)
!fre(1)=0.90*fre(2)
!fre(4)=1.10*fre(3)
!if(fre(2).ge.fre(3))stop 'Please check the corner frequency'
read(11,*)tb,te
if(tb.ge.te)stop 'Please check the time window'
close(11)
call initial_sachead(sacheado)
call read_sachead(file(1),sachead(1),nerr)
if(nerr.eq.-1)stop 'Read file error'
dt=sachead(1)%delta
nb=int((sachead(1)%t1+tb-sachead(1)%b)/dt)
ne=int((sachead(1)%t1+te-sachead(1)%b)/dt)
npts=int((te-tb)/dt)+1
sigo=0
sacheado%b=tb
sacheado%delta=dt
sacheado%npts=npts
call read_sac(file(1),sig(:,1),sachead(1),nerr)
if(nerr.eq.-1)stop 'Read file error'
!call filter(sig(:,1),sigf(:,1),fre,dt,sachead(1)%npts)
!   write(*,*)nb,ne
sigo(1:npts,1)=sig(nb:ne,1)
call write_sac('1.sac',sigo(:,1),sacheado,nerr)
do i=2,4
   write(outpt,'(1i1,".sac")')i
   call read_sachead(file(i),sachead(i),nerr)
   if(nerr.eq.-1)stop 'Read file error'
   nb=int((sachead(i)%t1+tb-sachead(i)%b)/dt)
   ne=int((sachead(i)%t1+te-sachead(i)%b)/dt)
   if(sachead(i)%delta.ne.dt)stop 'Please check the sampling rate'
   call read_sac(file(i),sig(:,i),sachead(i),nerr)
   if(nerr.eq.-1)stop 'Read file error'
   !call filter(sig(:,i),sigf(:,i),fre,dt,sachead(i)%npts)
!   write(*,*)nb,ne
   sigo(1:npts,i)=sig(nb:ne,i)
   call write_sac(outpt,sigo(:,i),sacheado,nerr)
enddo
sacheado%delta=dt
call do_conv(sigo(:,1),sigo(:,4),f12,npts,npts,npt)
call do_conv(sigo(:,2),sigo(:,3),f21,npts,npts,npt)
sacheado%npts=npt
sacheado%b=0
ncros=npt-npts-2
ncros=2*int(5/dt)
nb=npt/2-npts/2-1
ne=nb+npts-1
sacheado%t1=nb*dt
sacheado%t2=ne*dt
call write_sac('f21.sac',f21,sacheado,nerr)
aa(1:npts)=f21(nb:ne) ! take f21 as reference
sacheado%npts=npts
do i=1,ncros+1
   bb=0
   cc(i)=0
   n1=nb-ncros/2+i-1
   n2=n1+npts-1
   bb(1:npts)=f12(n1:n2)        
   call cross(aa,bb,npts,cc(i))
!   write(*,*)cc(i)
enddo
sacheado%npts=ncros
beg=-ncros/2*dt
sacheado%b=beg
call write_sac('test.sac',cc,sacheado,nerr)
write(*,*)(maxloc(cc)-1)*dt+beg,maxval(cc)
nid=maxloc(cc)-1
sacheado%npts=npt
sacheado%b=nid(1)*dt+beg
call write_sac('f12.sac',f12,sacheado,nerr)
!call write_sac('f12.sac',f12,sacheado,nerr)
!call write_sac('f12.sac',f12,sacheado,nerr)
end program
