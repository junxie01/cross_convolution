program add_noise
use sacio
implicit none
type(sac_head) :: sachead
integer i,npts,n1,nerr
character(80):: input, ouput,par
integer,parameter:: nsmax=4000000
real,dimension(nsmax) :: sig
real:: rad,per
if(iargc().ne.3)stop "Usage: add_noise input output percentage"
call getarg(1,input)
call getarg(2,ouput)
call getarg(3,par)
read(par,'(bn,f20.0)')per
call read_sachead(input,sachead,nerr)
call read_sac(input,sig,sachead,nerr)
call init_random_seed()
do i=1,sachead%npts
   call random_number(rad)
   sig(i)=sig(i)+(2*rad-1)*sachead%depmax*per/100.0
enddo
call write_sac(ouput,sig,sachead,nerr)
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
