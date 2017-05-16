program generate_triagular_function
use sacio
implicit none
type(sac_head) :: sachead
integer nsmax,i,npts,n1,nerr
parameter (nsmax=4000000)
character(80) file1,par
real(4) sig(nsmax),t1
if(iargc().ne.1)stop 'Usage: generate_triangular_function half_dualration'
call getarg(1,par)
read(par,'(bn,f20.0)')t1
call initial_sachead(sachead)
sachead%npts=1000
sachead%delta=0.1
sachead%b=0
sachead%e=(sachead%npts-1)*sachead%delta
n1=int(t1/sachead%delta)
sig(1:1000)=0
!write(*,*)sachead
do i=1,n1
   sig(i+500)=i
   sig(2*n1+1+500-i)=i
enddo
write(*,*)sum(sig*sig)
!do i=1,sachead%npts
!   write(*,*)0.1*(i-1),sig(i)
!enddo
sachead%t1=49.9
call write_sac('sta1_ev1.sac',sig,sachead,nerr)
sachead%t1=49.7
call write_sac('sta2_ev1.sac',sig,sachead,nerr)
sachead%t1=49.8
call write_sac('sta1_ev2.sac',sig,sachead,nerr)
sachead%t1=49.6
call write_sac('sta2_ev2.sac',sig,sachead,nerr)
end program
