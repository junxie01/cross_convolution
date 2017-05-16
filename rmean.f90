subroutine rmean(sig,npts)
integer npts,i
real sig(npts),avg/0/
avg=sum(sig(1:npts))/npts
sig(1:npts)=sig(1:npts)-avg
return
end subroutine
