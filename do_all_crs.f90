program main
use sacio
implicit none
integer,parameter :: nn=4000000,nstmax=1000,nfmax=1000,neqmax=1000
type(sac_head) :: sachead(4),sacheado
integer :: nb,ne,n1,n2,ncros,nnpts
integer :: ia,ib,npts,npt,shft,nerr
integer ieq1,ieq2,i,nsta,nevent,ieq,ist1,ist2
integer :: year(neqmax),month(neqmax),day(neqmax),hour(neqmax),mmin(neqmax),sec(neqmax),msec(neqmax)
character(8) :: sta(nstmax)
character(2) :: net(nstmax)
character(180) :: sac(4),nm1,nm2
character(80):: para,sta_list,eq_list,output,outpu
real(4) :: f12(nn),f21(nn)
real(4) :: t1(4),f1,f2,sig(nn,4),sigo(nn,4),sigt(nn,4)
real(4) :: fre(4),tb,te,az
real(4) :: cc(nn),amax,dt,aa(nn),bb(nn)
real(4) :: cc1,cc2,cc3,sf1(1),sf2(1),sf3(1)
real(4) :: evla(neqmax),evlo(neqmax),stlo1,stla1,stlo2,stla2,evdp(neqmax),mag(neqmax)
real(4) :: dist_ev,dist_st
logical :: ext(4)
if(iargc().ne.1)then
   write(*,'("Usage: do_all_crs param.dat")')
   write(*,'("param.dat:")')
   write(*,'("station list")')
   write(*,'("earthquake list")')
   write(*,'("t_before t_after")')
   write(*,'("output")')
   stop
endif
call getarg(1,para)
open(9,file=para)
read(9,'(a80)')sta_list
read(9,'(a80)')eq_list
read(9,*)tb,te
read(9,*)output
close(9)
open(10,file=sta_list)
do i=1,nstmax
   read(10,*,end=12,err=12)sta(i),net(i)
enddo
12  close(10)
nsta=i-1
open(12,file=eq_list)
do ieq=1,neqmax                                  ! loop over year
   read(12,*,end=14,err=14)year(ieq),month(ieq),day(ieq),hour(ieq),mmin(ieq),&
   sec(ieq),msec(ieq),evla(ieq),evlo(ieq),evdp(ieq),mag(ieq)
enddo
14 close(12)
nevent=ieq-1
write(*,*)'There are ',nevent,' events.'

open(17,file=output)
do ieq1=1,nevent-1
   do ieq2=ieq1+1,nevent
      call cal_az_dist(evlo(ieq1),evla(ieq1),evlo(ieq2),evla(ieq2),dist_ev)
      do ist1=1,nsta-1
         write(sac(1),'(i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i3.3,"_",1a,"_",1a,"_BHZ.SAC")')&
         year(ieq1),month(ieq1),day(ieq1),hour(ieq1),mmin(ieq1),sec(ieq1),msec(ieq1),trim(net(ist1)),trim(sta(ist1))
         write(sac(2),'(i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i3.3,"_",1a,"_",1a,"_BHZ.SAC")')&
         year(ieq2),month(ieq2),day(ieq2),hour(ieq2),mmin(ieq2),sec(ieq2),msec(ieq2),trim(net(ist1)),trim(sta(ist1))
!         write(*,*)'sac1=',trim(sac(1))
!         write(*,*)'sac2=',trim(sac(2))
         do ist2=ist1+1,nsta
            write(sac(3),'(i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i3.3,"_",1a,"_",1a,"_BHZ.SAC")')&
            year(ieq1),month(ieq1),day(ieq1),hour(ieq1),mmin(ieq1),sec(ieq1),msec(ieq1),trim(net(ist2)),trim(sta(ist2))
            write(sac(4),'(i4.4,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i2.2,"_",i3.3,"_",1a,"_",1a,"_BHZ.SAC")')&
            year(ieq2),month(ieq2),day(ieq2),hour(ieq2),mmin(ieq2),sec(ieq2),msec(ieq2),trim(net(ist2)),trim(sta(ist2))
!            write(*,*)'sac3=',trim(sac(3))
!            write(*,*)'sac4=',trim(sac(4))
            do i=1,4
               inquire(file=sac(i),exist=ext(i)) 
            enddo
            if(.not.ext(1))cycle
            if(.not.ext(2))cycle
            if(.not.ext(3))cycle
            if(.not.ext(4))cycle
            call initial_sachead(sacheado)
            do i=1,4
               call read_sachead(sac(i),sachead(i),nerr)
               if(nerr.eq.-1)stop 'Read file error'
               call read_sac(sac(i),sig(:,i),sachead(i),nerr)
               if(nerr.eq.-1)stop 'Read file error'
               ia=int((sachead(i)%t1-sachead(i)%b+tb)/sachead(i)%delta)
               ib=int((sachead(i)%t1-sachead(i)%b+te)/sachead(i)%delta)
               npts=int((te-tb)/sachead(i)%delta)+1
               sigo(1:npts,i)=sig(ia:ib,i)
            enddo
            dt=sachead(1)%delta
            call cal_az_dist(sachead(1)%stlo,sachead(1)%stla,sachead(4)%stlo,sachead(4)%stla,dist_st)
            call do_conv(sigo(:,1),sigo(:,4),f12,npts,npts,npt)
            call do_conv(sigo(:,2),sigo(:,3),f21,npts,npts,npt)
            sacheado%npts=npt
            sacheado%b=0
            cc=0
! use f12 f21 to calculate the time shift
!            call cross_td(f12,f21,cc,npt,shft)
            ncros=2*int(5/dt)
            nb=npt/2-npts/2-1
            ne=nb+npts-1
            aa(1:npts)=f21(nb:ne)
            do i=1,ncros+1
               bb=0
               n1=nb-ncros/2+i-1
               n2=npts+n1-1
               bb(1:npts)=f12(n1:n2)        
               call cross(aa,bb,npts,cc(i))
            enddo
            sf1=(maxloc(cc)-1)*dt-ncros/2*dt
            cc1=maxval(cc)


! use direct cross correlation to calculate the time shift
            cc=0
            ia=int((sachead(1)%t1-sachead(1)%b+tb)/dt)
            ib=int((sachead(1)%t1-sachead(1)%b+te)/dt)
            aa(1:npts)=sig(ia:ib,1)
            do i=1,ncros+1
               ia=int((sachead(3)%t1+tb-sachead(3)%b)/dt)
               n1=ia-ncros/2+i-1
               n2=n1+npts-1
               bb(1:npts)=sig(n1:n2,3)
               call cross(aa,bb,npts,cc(i))
            enddo
            !call cross_td(sigo(:,3),sigo(:,1),cc,npts,shft) 
            sf2=-ncros/2*dt+(maxloc(cc)-1)*dt
            cc2=maxval(cc) 



            cc=0
            ia=int((sachead(2)%t1-sachead(2)%b+tb)/dt)
            ib=int((sachead(2)%t1-sachead(2)%b+te)/dt)
            aa(1:npts)=sig(ia:ib,2)
            do i=1,ncros+1
               n1=int((sachead(4)%t1+tb-sachead(4)%b)/dt)-ncros/2+i-1
               n2=n1+npts-1
               bb(1:npts)=sig(n1:n2,4)
               call cross(aa,bb,npts,cc(i))
            enddo
            sf3=-ncros/2*dt+(maxloc(cc)-1)*dt
            cc3=maxval(cc) 



            !write(17,'(8f10.4)')dist_ev,dist_st,cc1,sf1(1),cc2,sf2(1),cc3,sf3(1)
            write(nm1,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,".",i3.3)')&
            year(ieq1),month(ieq1),day(ieq1),hour(ieq1),mmin(ieq1),sec(ieq1),msec(ieq1)
            write(nm2,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,".",i3.3)')&
            year(ieq2),month(ieq2),day(ieq2),hour(ieq2),mmin(ieq2),sec(ieq2),msec(ieq2)
            !write(17,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,".",i3.3,1x,i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",&
            !i2.2,".",i3.3,1x,2a,1x,2a,1x,8f10.4)')year(ieq1),month(ieq1),day(ieq1),hour(ieq1),mmin(ieq1),sec(ieq1),msec(ieq1)&
            !,year(ieq2),month(ieq2),day(ieq2),hour(ieq2),mmin(ieq2),sec(ieq2),msec(ieq2),net(ist1),trim(sta(ist1)), &
            !net(ist2),trim(sta(ist2)),dist_ev,dist_st,cc1,sf1(1),cc2,sf2(1),cc3,sf3(1)
            write(17,'(1a,1x,1a,1x,2a,1x,2a,1x,8f10.4)')&
            trim(nm1),trim(nm2),net(ist1),trim(sta(ist1)),net(ist2),trim(sta(ist2)),dist_ev,dist_st,cc1,sf1(1),cc2,sf2(1),cc3,sf3(1)
            write(*,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,".",i3.3,1x,i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",&
            i2.2,".",i3.3,1x,2a,1x,2a,1x,8f10.4)')year(ieq1),month(ieq1),day(ieq1),hour(ieq1),mmin(ieq1),sec(ieq1),msec(ieq1)&
            ,year(ieq2),month(ieq2),day(ieq2),hour(ieq2),mmin(ieq2),sec(ieq2),msec(ieq2),net(ist1),trim(sta(ist1)), &
            net(ist2),trim(sta(ist2)),dist_ev,dist_st,cc1,sf1(1),cc2,sf2(1),cc3,sf3(1)
            !write(*, '(8f10.4)')dist_ev,dist_st,cc1,sf1(1),cc2,sf2(1),cc3,sf3(1)
         enddo ! loop over station1
      enddo ! loop over stataion2 
   enddo ! loop over event2
enddo ! loop over event1
close(17)
end program
