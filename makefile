#makefile
FC=gfortran
FFLAG=-fbounds-check
all:sacio.mod do_conv crs_cov_mea do_cross add_noise generate_trai_func add_noise1 do_all_crs crs_cov_mea_wf do_cross_more
objects1=do_conv.o do_convolution.o sacio.o conv1.o zfour.o clogc.o
objects2=crs_cov_mea.o do_conv.o sacio.o clogc.o cal_crst.o rmean.o conv1.o cor.o zfour.o clogc.o cross_td.o cross.o
objects3=do_cross.o cross_td.o sacio.o cross_fd.o cor.o clogc.o cross.o
objects4=add_noise.o sacio.o
objects5=generate_trai_func.o sacio.o
objects6=add_noise1.o sacio.o
objects7=do_all_crs.o sacio.o cross_td.o clogc.o cor.o do_conv.o cal_crst.o rmean.o  cal_dist_and_az.o filter.o taperf.o conv1.o cross.o
objects8=crs_cov_mea_wf.o do_conv.o sacio.o clogc.o cal_crst.o rmean.o conv1.o cor.o zfour.o clogc.o cross_td.o cross.o taperf.o filter.o
objects9=do_cross_more.o cross_td.o sacio.o cross_fd.o cor.o clogc.o cross.o
do_conv: $(objects1)
	$(FC) $^ -o $@ 
crs_cov_mea: $(objects2)
	$(FC) $^ -o $@ 
crs_cov_mea_wf: $(objects8)
	$(FC) $^ -o $@ 
do_cross: $(objects3)
	$(FC) $^ -o $@ 
do_cross_more: $(objects9)
	$(FC) $^ -o $@ 
add_noise:$(objects4)
	$(FC) $^ -o $@ 
generate_trai_func:$(objects5)
	$(FC) $^ -o $@ 
add_noise1:$(objects6)
	$(FC) $^ -o $@ 
do_all_crs:$(objects7)
	$(FC) $^ -o $@ 
%.o: %.f90
	$(FC) -c $^ $(FFLAG)
sacio.mod:sacio.f90
	$(FC) -c $^ 
install:
	cp do_conv crs_cov_mea crs_cov_mea_wf do_cross add_noise add_noise1 do_all_crs do_cross_more ../bin
clean:
	-rm *.o *.mod
