subroutine smooth(sig, npts, halfn)
integer*4 namx
parameter(nmax=4000000)
integer*4 ii, jj, k,num/0/,halfn,n1
real sig(nmax),temp(nmax),sig_out(nmax)
real max_val/0/,summ
do ii=1,npts
   k=0;summ=0
   do jj=ii-halfn,ii+halfn
      if (jj.gt.0.and. jj.lt.npts+1) then
         summ=summ+abs(sig(jj))
         k=k+1
      endif
   enddo
   temp(ii)=summ/k
   if(max_val.lt.temp(ii))max_val=temp(ii)
enddo
sig=temp+0.00001*max_val
end subroutine
