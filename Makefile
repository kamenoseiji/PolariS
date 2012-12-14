#----------------- LIBRARY -------------------
PGDIR = /usr/local/pgplot
CPGDIR = $(PGDIR)
PGPLOT= $(PGDIR)/libpgplot.a
CPGPLOT= $(CPGDIR)/libcpgplot.a
BINDIR = /usr/local/custom/bin
#----------------- LINK OPTIONS -------------------
CFLAGS= -I/usr/local/pgplot -I/usr/X11R6/include
CCOMPL=gcc $(CFLAGS)
NVCC=nvcc -I/usr/local/cuda/include -I/usr/local/cuda/common/inc
FCOMPL=gfortran 
#------- Followings are PASS or DIRECTORY -------
PROGS=	polaris_start shm_alloc shm_init shm_param_view k5sample_store cuda_fft_xspec shm_segdata
GRLIBS= -L/usr/X11R6/lib -lX11
MATH=	-lm
FFTLIB= -lcufft
CUDALIB= -lcutil
#----------------- MAPPING ------------------------
OBJ_start= polaris_start.o shm_access.o
OBJ_shm_alloc= shm_alloc.o shm_init_create.o shm_access.o erase_shm.o
OBJ_shm_init = shm_init.o shm_access.o
OBJ_shm_view = shm_param_view.o shm_access.o timesystem.o
OBJ_k5sample_store = k5sample_store.o shm_access.o
OBJ_shm_segdata = shm_segdata.o
OBJ_cuda_fft = cuda_fft_xspec.o shm_access.o
#----------------- Compile and link ------------------------
polaris_start : $(OBJ_start)
	$(CCOMPL) -o $@ $(OBJ_start)

k5sample_store : $(OBJ_k5sample_store)
	$(CCOMPL) -o $@ $(OBJ_k5sample_store)

shm_alloc : $(OBJ_shm_alloc)
	$(CCOMPL) -o $@ $(OBJ_shm_alloc)

shm_init : $(OBJ_shm_init)
	$(CCOMPL) -o $@ $(OBJ_shm_init)

shm_param_view : $(OBJ_shm_view)
	$(CCOMPL) -o $@ $(OBJ_shm_view)

shm_segdata : $(OBJ_shm_segdata)
	$(CCOMPL) -o $@ $(OBJ_shm_segdata)

cuda_fft_xspec : $(OBJ_cuda_fft)
	$(NVCC) -o $@ $(OBJ_cuda_fft) $(FFTLIB) $(CUDALIB)

clean :
	\rm $(PROGS) *.o a.out core *.trace

all :	$(PROGS)

install:
	@mv $(PROGS) $(BINDIR)

#----------------- Objects ------------------------
#.cu.o:
#	$(NVCC) -c $*.cu
cuda_fft_xspec.o:	cuda_fft_xspec.cu	shm_k5data.inc
	$(NVCC) -c cuda_fft_xspec.cu

.c.o:
	$(CCOMPL) -c $*.c

k5sample_store.o:	k5sample_store.c	shm_k5data.inc
polaris_start.o:	polaris_start.c		shm_k5data.inc
shm_alloc.o:		shm_alloc.c			shm_k5data.inc
shm_init.o:			shm_init.c			shm_k5data.inc
erase_shm.o:		erase_shm.c			shm_k5data.inc
shm_param_view.o:	shm_param_view.c	shm_k5data.inc
shm_segdata.o:		shm_segdata.c		shm_k5data.inc
shm_access.o:		shm_access.c		shm_k5data.inc
shm_init_create.o:	shm_init_create.c	shm_k5data.inc
timesystem.o:		timesystem.c

#----------------- End of File --------------------
