//	cuda_fft_xspec.c : FFT using CuFFT
//
//	Author : Seiji Kameno
//	Created: 2012/12/6
//
#include "shm_k5data.inc"
#include "k5dict.inc"
#include <cuda.h>
#include <cutil.h>
#include <cufft.h>
#define	NFFT		262144
#define	NSEG		256

__device__ float2 complexMult(float2 a, float2 b)		// output a x b
{
	return make_float2( a.x* b.x - a.y* b.y, a.x* b.y + a.y* b.x );
}

__device__ float2 complexMultConj(float2 a, float2 b)		// output a x b*
{
	return make_float2( a.x* b.x + a.y* b.y,  a.y* b.x - a.x* b.y );
}

__device__ float complexMod( float2 a )				// output |a|^2
{
	return  a.x* a.x + a.y* a.y 
}

__global__ void complexMultConjVec(		// calculate a x b*
	float2	*vec_in_a,			// Input vector
	float2	*vec_in_b,			// Input vector
	float2	*vec_out_c,			// Output vector
	int		length)				// Vector length
{
	int tid = blockIdx.x* blockDim.x + threadIdx.x;
	if((tid >= 0) && (tid < length)){
		vec_out_c[tid] = complexMultConj(vec_in_a[tid], vec_in_b[tid]);
	}
}

__global__ void complexPowerVec(		// calculate a x a*
	float2	*vec_in,		// Input vector
	float	*vec_out,
	int		length)
{
	int tid = blockIdx.x* blockDim.x + threadIdx.x;
	if((tid >= 0) && (tid < length)){
		vec_out[tid] = complexMod(vec_in_a[tid]);
	}
}


main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data access
	float	*segdata_ptr;				// Pointer to sampled data

	cufftHandle		cufft_plan;
	cufftComplex	*cuSpec_data;
	cufftReal		*cuReal_data;
/*
-------------------------------------------- Prepare for CuFFT
*/
	cudaMalloc( (void **)&cuReal_data, SEGDATA_SIZE/2 );
	cudaMalloc( (void **)&cuSpec_data, SEGDATA_SIZE/2 );
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
	segdata_ptr = (float *)shmat(param_ptr->shrd_seg_id, NULL, SHM_RDONLY);
/*
-------------------------------------------- K5 Header and Data
*/
	setvbuf(stdout, (char *)NULL, _IONBF, 0);   // Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		printf("FFT status EW\r");
		sops.sem_num = (ushort)8; sops.sem_op = (short)-4; sops.sem_flg = (short)0;
        semop( param_ptr->sem_data_id, &sops, 1);

		printf("FFT status ES\r");
		cudaMemcpy( cuReal_data, segdata_ptr, SEGDATA_SIZE/2, cudaMemcpyHostToDevice);
		printf("FFT status EF\r");
		cufftExecR2C(cufft_plan, cuReal_data, cuSpec_data);

		printf("FFT status OW\r");
		segdata_ptr += MAX_IF* MAX_seg_len* MAX_seg_sec / 2;
		sops.sem_num = (ushort)9; sops.sem_op = (short)-4; sops.sem_flg = (short)0;
        semop( param_ptr->sem_data_id, &sops, 1);
		printf("FFT status OS\r");
		cudaMemcpy( cuReal_data, segdata_ptr, SEGDATA_SIZE/2, cudaMemcpyHostToDevice);
		segdata_ptr -= MAX_IF* MAX_seg_len* MAX_seg_sec / 2;
		printf("FFT status OF\r");

		cufftExecR2C(cufft_plan, cuReal_data, cuSpec_data);
	}
/*
-------------------------------------------- RELEASE the SHM
*/
	cufftDestroy(cufft_plan);
	cudaFree(&cuReal_data);
	cudaFree(&cuSpec_data);
    return(0);
}
