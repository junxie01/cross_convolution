program add_noise
use sacio
implicit none
type(sac_head) :: sachead
integer i,npts,n1,nerr
character(80):: input1,input2,input3,input4,par,output
integer,parameter:: nsmax=4000000
real,dimension(nsmax) :: sig
real:: rad,per,depmax
if(iargc().ne.5)stop "Usage: add_noise input1 input2 input3 input4 percentage"
call getarg(1,input1) 
call getarg(2,input2) 
call getarg(3,input3) 
call getarg(4,input4) 
call getarg(5,par)
read(par,'(bn,f20.0)')per
call read_sachead(input1,sachead,nerr)
call read_sac(input1,sig,sachead,nerr)
call init_random_seed()
depmax=sachead%depmax
do i=1,sachead%npts
   call random_number(rad)
   sig(i)=sig(i)+(2*rad-1)*depmax*per/100.0
enddo
output=trim(input1)//"_no"
call write_sac(output,sig,sachead,nerr)

call read_sachead(input2,sachead,nerr)
call read_sac(input2,sig,sachead,nerr)
do i=1,sachead%npts
   call random_number(rad)
   sig(i)=sig(i)+(2*rad-1)*depmax*per/100.0
enddo
output=trim(input2)//"_no"
call write_sac(output,sig,sachead,nerr)

call read_sachead(input3,sachead,nerr)
call read_sac(input3,sig,sachead,nerr)
do i=1,sachead%npts
   call random_number(rad)
   sig(i)=sig(i)+(2*rad-1)*depmax*per/100.0
enddo
output=trim(input3)//"_no"
call write_sac(output,sig,sachead,nerr)

call read_sachead(input4,sachead,nerr)
call read_sac(input4,sig,sachead,nerr)
do i=1,sachead%npts
   call random_number(rad)
   sig(i)=sig(i)+(2*rad-1)*depmax*per/100.0
enddo
output=trim(input4)//"_no"
call write_sac(output,sig,sachead,nerr)
end program

subroutine init_random_seed()
integer :: i, n, clock
integer, dimension(:), allocatable :: seed
call random_seed(size = n)
allocate(seed(n))
call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
call random_seed(put = seed)
deallocate(seed)
end subroutine init_random_seed
