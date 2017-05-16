program calculate_convolution
use sacio
implicit none
integer nptsmax
type(sac_head):: sachead1,sachead2,sacheado
parameter (nptsmax=1000000)
integer i,npts,nerr,npt,id
character(80):: file1,file2,outpt,param
real(4) :: sig1(nptsmax),sig2(nptsmax)
real(4) :: sigo(nptsmax)
if(iargc().ne.4)stop 'Usage: do_convolve file1 file2 output id'
call getarg(1,file1)
call getarg(2,file2)
call getarg(3,outpt)
call getarg(4,param)
read(param,'(bn,i20)')id
call read_sachead(file1,sachead1,nerr)
if(nerr.eq.-1)stop "Something wrong with file 1"
call read_sac(file1,sig1,sachead1,nerr)
write(*,*)(maxloc(sig1)-1)*sachead1%delta+sachead1%b
write(*,*)(maxloc(sig1)-1)*sachead1%delta
if(nerr.eq.-1)stop "Something wrong with file 1"

call read_sachead(file2,sachead2,nerr)
if(nerr.eq.-1)stop "Something wrong with file 2"
call read_sac(file2,sig2,sachead2,nerr)
if(nerr.eq.-1)stop "Something wrong with file 2"
if(sachead1%delta.ne.sachead2%delta)stop "the sampling rate is different dude."
if(id.eq.1)call do_conv(sig1,sig2,sigo,sachead1%npts,sachead2%npts,npt)
if(id.eq.2)call conv1(sig1,sig2,sigo,sachead1%npts,sachead2%npts,npt)
call initial_sachead(sacheado)
sacheado%delta=sachead2%delta
sacheado%npts=npt
sacheado%b=0
call write_sac(outpt,sigo,sacheado,nerr)
write(*,*)(maxloc(sigo)-1)*sacheado%delta
end program
