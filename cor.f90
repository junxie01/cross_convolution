subroutine cor(src,data,nft,nk,npow,dt)
complex src(nft),data(nft)
integer nft,npow,i,nk
real aa,dt
aa=-1
do i=1,nft
   src(i)=aa*src(i)*conjg(data(i))
   aa=-aa
enddo
call clogc(npow,src,-1,dt)
return
end subroutine cor
