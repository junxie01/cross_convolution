subroutine cross_td(x,y,zz,n,n21)
implicit none
integer nsmax
parameter(nsmax=1000000)
integer i,j,n,ii,nn,jj,n21
real x(nsmax),y(nsmax),zz(nsmax)
n21=int((n-1)/2)
nn=2*n21+1;zz(1:nn)=0
do ii=1,nn
   do i=1,n
      j=i+ii-n21-1
      if(j.le.n.and.j.gt.0)zz(ii)=zz(ii)+x(i)*y(j)
   enddo
enddo
zz(1:nn)=zz(1:nn)/sqrt(sum(x(1:n)**2))/sqrt(sum(y(1:n)**2))
end subroutine
