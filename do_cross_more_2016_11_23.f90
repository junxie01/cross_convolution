program do_cross
use sacio
implicit none
integer,parameter :: nmax=1000000
type(sac_head) :: sachead1,sachead2,sacheado
integer :: ia,ib,i,ncros,n1,n2,npts
integer :: nerr,nft,id,nid(1),nn,na1,nb1,na2,nb2
real :: st,tb,te,dt,aa(nmax),bb(nmax),cc(nmax)
real :: sig1(nmax),sig2(nmax),sig3(nmax),tw,tw1
character(180) :: file1,file2,par,output
if(iargc().ne.6)stop 'Usage: do_cross file1 file2 tb te tw tw1'
call getarg(1,file1)
call getarg(2,file2)
call getarg(3,par)
read(par,'(bn,f20.0)')tb
call getarg(4,par)
read(par,'(bn,f20.0)')te
call getarg(5,par)
read(par,'(bn,f20.0)')tw
call getarg(6,par)
read(par,'(bn,f20.0)')tw1
if(tw1.gt.tw)stop 'the subwindow is bigger than main window'
call read_sachead(file1,sachead1,nerr)
call read_sac(file1,sig1,sachead1,nerr)
dt=sachead1%delta
call read_sachead(file2,sachead2,nerr)
if(sachead2%delta.ne.dt)stop 'Sampling rate not match'
call read_sac(file2,sig2,sachead2,nerr)

ia=int((sachead1%t1-sachead1%b+tb)/dt)
ib=int((sachead1%t1-sachead1%b+te)/dt)
npts=int((te-tb)/dt)+1
aa(1:npts)=sig1(ia:ib)
!ncros=2*npts
ncros=2*int(tw/dt)
cc=0
do i=1,ncros+1
   ia=int((sachead2%t1+tb-sachead2%b)/dt)
   n1=ia-ncros/2+i-1
   n2=n1+npts-1
   bb(1:npts)=sig2(n1:n2)
   call cross(aa,bb,npts,cc(i))
enddo
call initial_sachead(sacheado)
sacheado%delta=sachead1%delta
sacheado%b=-ncros*sacheado%delta/2.0
sacheado%npts=ncros+1
output=trim(file2)//".cc"
n1=int((tw-tw1)/dt)
n2=int((tw+tw1)/dt)
call write_sac(output,cc,sacheado,nerr)
nid=maxloc(cc(n1:n2))
nn=int(27/dt)
na1=nid(1)-nn+n1
if(na1.lt.1)na1=1
nb1=nid(1)-200+n1
na2=nid(1)+200+n1
nb2=nid(1)+nn+n1
if(nb2.gt.ncros)nb2=ncros+1
!write(*,*)n1,n2,na1,nb1,na2,nb2,nid(1),nn
write(*,*)maxval(cc(n1:n2)),(nid(1)-1)*dt-ncros/2*dt+(n1-1)*dt,maxval(cc(na1:nb1)),&
(maxloc(cc(na1:nb1))-1)*dt-ncros/2*dt+(na1-1)*dt,maxval(cc(na2:nb2)),&
(maxloc(cc(na2:nb2))-1)*dt-ncros/2*dt+(na2-1)*dt
nid=maxloc(cc(1:ncros))
!write(*,*)maxval(cc),sacheado%b+(nid(1)-1)*sacheado%delta
!sachead1%b=sachead1%b+sacheado%b+(nid(1)-1)*sacheado%delta
sachead2%t1=sachead2%t1-ncros/2*dt+(nid(1)-1)*dt
output=trim(file2)//'_shift'
call write_sac(output,sig2,sachead2,nerr)
end program
