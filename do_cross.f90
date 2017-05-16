program do_cross
use sacio
implicit none
integer,parameter :: nmax=1000000
type(sac_head) :: sachead1,sachead2,sacheado
integer :: nerr,nft,id,nid(1)
integer :: ia,ib,i,ncros,n1,n2,npts
real :: sig1(nmax),sig2(nmax),sig3(nmax),tw
real :: st,tb,te,dt,aa(nmax),bb(nmax),cc(nmax)
character(180) :: file1,file2,par,output
if(iargc().ne.5)stop 'Usage: do_cross file1 file2 tb te tw'
call getarg(1,file1)
call getarg(2,file2)
call getarg(3,par)
read(par,'(bn,f20.0)')tb
call getarg(4,par)
read(par,'(bn,f20.0)')te
call getarg(5,par)
read(par,'(bn,f20.0)')tw
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
call write_sac(output,cc,sacheado,nerr)
write(*,*)maxval(cc),(maxloc(cc)-1)*dt-ncros/2*dt
nid=maxloc(cc(1:ncros))
!write(*,*)maxval(cc),sacheado%b+(nid(1)-1)*sacheado%delta
!sachead1%b=sachead1%b+sacheado%b+(nid(1)-1)*sacheado%delta
sachead2%t1=sachead2%t1-ncros/2*dt+(nid(1)-1)*dt
output=trim(file2)//'_shift'
call write_sac(output,sig2,sachead2,nerr)
end program
