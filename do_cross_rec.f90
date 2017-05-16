program do_cross
use sacio
implicit none
integer,parameter :: nmax=1000000
type(sac_head) :: sachead1,sachead2,sacheado
integer :: nerr,nft,id,nid(1)
integer :: ia,ib,i,ncros,n1,n2,npts
real :: sig1(nmax),sig2(nmax),sig3(nmax),tw
real :: st,tb,te,dt,aa(nmax),bb(nmax),cc(nmax),rr
character(180) :: file1,file2,par,output
if(iargc().ne.4)stop 'Usage: do_cross file1 file2 tb te'
call getarg(1,file1)
call getarg(2,file2)
call getarg(3,par)
read(par,'(bn,f20.0)')tb
call getarg(4,par)
read(par,'(bn,f20.0)')te
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
bb(1:npts)=sig2(ia:ib)
call cross(aa,bb,npts,rr)
write(*,*)trim(file2),rr
end program
