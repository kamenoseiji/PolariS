//	cuda_fft_segmeng.c : 1-segment FFT using CuFFT
//
//	Author : Seiji Kameno
//	Created: 2012/12/6
//
#include "shm_k5data.inc"
#include "k5dict.inc"
#include <stdio.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/types.h>
#include <sys/shm.h>
#include <cuda.h>
#include <cufft.h>
#define	NFFT		262144
#define	NSEG		128

// int	shm_access(key_t shm_key, size_t shm_size, int *shrd_id, int *shm_ptr);

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	unsigned char		*k5head_ptr;	// Pointer to the shared K5 header
	unsigned char		*k5data_ptr;	// Pointer to the shared K5 data
	float	*host_segdata;
	int		index, index_seg, offset;

	cufftHandle		cufft_plan;
	cufftComplex	*seg_data;
	cufftReal		*seg_real_data;
/*
-------------------------------------------- Prepare for CuFFT
*/
	host_segdata = (float *)malloc(NSEG* NFFT* sizeof(float));

	cudaMalloc( (void **)&seg_data, NSEG* sizeof(cufftComplex)* (NFFT/2+1) );
	cudaMalloc( (void **)&seg_data, NSEG* sizeof(cufftComplex)* (NFFT/2+1) );
	cudaMalloc( (void **)&seg_real_data, NSEG* sizeof(cufftReal)* NFFT );
	if(cudaGetLastError() != cudaSuccess){
		fprintf(stderr, "Cuda Error : Failed to allocate memory.\n"); return(-1);
	}

	if(cufftPlan1d(&cufft_plan, NFFT, CUFFT_R2C, NSEG) != CUFFT_SUCCESS){
		fprintf(stderr, "Cuda Error : Failed to create plan.\n"); return(-1);
	}
	
/*
-------------------------------------------- Access to the SHARED MEMORY
*/
	shrd_param_id = shmget( SHM_PARAM_KEY, sizeof(struct SHM_PARAM), 0444);
	param_ptr = (struct SHM_PARAM *)shmat(shrd_param_id, NULL, 0);
	k5head_ptr = (unsigned char *)shmat(param_ptr->shrd_k5head_id, NULL, SHM_RDONLY);
	k5data_ptr = (unsigned char *)shmat(param_ptr->shrd_k5data_id, NULL, SHM_RDONLY);
/*
-------------------------------------------- K5 Header and Data
*/
	setvbuf(stdout, (char *)NULL, _IONBF, 0);   // Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }
		for(index=0; index < NFFT; index ++){
			for(index_seg=0; index_seg<NSEG; index_seg += 4){
				offset = index_seg* NFFT + index;
				host_segdata[          offset] = wt4[k5data_ptr[          offset] & 0x0f];
				host_segdata[   NFFT + offset] = wt4[k5data_ptr[   NFFT + offset] & 0x0f];
				host_segdata[2* NFFT + offset] = wt4[k5data_ptr[2* NFFT + offset] & 0x0f];
				host_segdata[3* NFFT + offset] = wt4[k5data_ptr[3* NFFT + offset] & 0x0f];
			}
		}
		printf("p\n");
		cudaMemcpy( seg_real_data, host_segdata, NSEG* NFFT* sizeof(float), cudaMemcpyHostToDevice);
		cufftExecR2C(cufft_plan, seg_real_data, seg_data);
		cudaMemcpy( host_segdata, (cufftReal *)seg_data, NSEG* NFFT* sizeof(float), cudaMemcpyDeviceToHost);

		printf("%e %e %e %e\n", host_segdata[0], host_segdata[1], host_segdata[2], host_segdata[3]);
		usleep(200);
	}
/*
-------------------------------------------- RELEASE the SHM
*/
	free(host_segdata);
	cufftDestroy(cufft_plan);
	cudaFree(seg_data);
    return(0);
}
