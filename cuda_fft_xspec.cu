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
	return  a.x* a.x + a.y* a.y;
}

__global__ void complexMultConjVec(		// calculate a x b*
	float2	*vec_in_a,			// Input vector
	float2	*vec_in_b,			// Input vector
	float2	*vec_out,			// Output vector
	int		length)				// Vector length
{
	int tid = blockIdx.x* blockDim.x + threadIdx.x;
	if((tid >= 0) && (tid < length)){
		vec_out[tid] = complexMultConj(vec_in_a[tid], vec_in_b[tid]);
	}
}

__global__ void complexPowerVec(		// calculate a x a*
	float2	*vec_in,		// Input vector
	float	*vec_out,		// Output vector
	int		length)			// Number of elements
{
	int tid = blockIdx.x* blockDim.x + threadIdx.x;
	if((tid >= 0) && (tid < length)){
		vec_out[tid] = complexMod(vec_in[tid]);
	}
}

__global__ void accumReal(	// a <- a + b
	float	*vec_in_a,		// Accumuration Results
	float	*vec_in_b,		// to be accumulated
	int		length)
{
    int tid = blockIdx.x* blockDim.x + threadIdx.x;
    if((tid >= 0) && (tid < length)){
        vec_in_a[tid] += vec_in_b[tid];
    }
}

__global__ void accumComplex(	// a <- a + b
	float2	*vec_in_a,		// Accumuration Results
	float2	*vec_in_b,		// to be accumulated
	int		length)
{
    int tid = blockIdx.x* blockDim.x + threadIdx.x;
    if((tid >= 0) && (tid < length)){
        vec_in_a[tid].x += vec_in_b[tid].x;
        vec_in_a[tid].y += vec_in_b[tid].y;
    }
}

__global__ void accumPowerSpec(
	float2	*vec_in,		// Input vector
	float	*vec_out,		// Output vector
	int	Nx, int	Ny)
{
    int ix = blockIdx.x* blockDim.x + threadIdx.x;
    int iy;

	for(iy=0; iy<Ny; iy++){ vec_out[ix] += complexMod(vec_in[Nx* iy + ix]); }
}

__global__ void accumCrossSpec(
	float2	*vec_in_a,		// Input vector
	float2	*vec_in_b,		// Input vector
	float2	*vec_out,		// Output vector
	int	Nx, int	Ny)
{
    int ix = blockIdx.x* blockDim.x + threadIdx.x;
    int iy;

	for(iy=0; iy<Ny; iy++){
		vec_out[ix].x += complexMultConj(vec_in_a[Nx* iy + ix], vec_in_b[Nx* iy + ix]).x;
		vec_out[ix].y += complexMultConj(vec_in_a[Nx* iy + ix], vec_in_b[Nx* iy + ix]).y;
	}
}

//-------- Convert SoD (Second of Day) into hour, min, and second
int sod2hms(
	int	sod,		// Second of Day
	int	*hour,		// Hour
	int	*min,		// Min
	int	*sec)		// Sec
{
	*hour = sod / 3600;
	*min  = (sod % 3600) / 60;
	*sec  = (sod % 60);
	return(*sec);
}

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	int		if_index;					// Index for IF
	int     sod=0, hour=0, min=0, sec=0;// Second of day, hour, min, sec
	unsigned char		*k5head_ptr;	// Pointer to the K5 header
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data access
	float	*init_segdata_ptr, *segdata_ptr;				// Pointer to sampled data
	float	*xspec_ptr;					// Pointer to 1-sec-integrated Power Spectrum
	FILE	*file_ptr;					// File Pointer to write

	dim3			Dg, Db(512,1, 1);	// Grid and Block size
	cufftHandle		cufft_plan;			// 1-D FFT Plan, to be used in cufft
	cufftComplex	*cuSpecData;		// FFTed spectrum, every IF, every segment
	cufftComplex	*cuCrossSpec;		// (crosscorrelation) Cross Spectrum
	cufftReal		*cuRealData;		// Time-beased data before FFT, every IF, every segment
	cufftReal		*cuPowerSpec;		// (autocorrelation) Power Spectrum

//------------------------------------------ Prepare for Output Data
	file_ptr = fopen("/DATA/PolariS.dat", "w");
//------------------------------------------ Prepare for CuFFT
	cudaMalloc( (void **)&cuRealData, Nif* NsegSec2* SegLen * sizeof(cufftReal) );
	cudaMalloc( (void **)&cuSpecData, Nif* NsegSec2* NFFTC* sizeof(cufftComplex) );
	cudaMalloc( (void **)&cuPowerSpec, Nif* NFFTC* sizeof(float));
	cudaMalloc( (void **)&cuCrossSpec, Nif/2* NFFTC* sizeof(float2));

	if(cudaGetLastError() != cudaSuccess){
		fprintf(stderr, "Cuda Error : Failed to allocate memory.\n"); return(-1);
	}

	if(cufftPlan1d(&cufft_plan, NFFT, CUFFT_R2C, NsegSec2* Nif ) != CUFFT_SUCCESS){
		fprintf(stderr, "Cuda Error : Failed to create plan.\n"); return(-1);
	}
//------------------------------------------ Access to the SHARED MEMORY
	shrd_param_id = shmget( SHM_PARAM_KEY, sizeof(struct SHM_PARAM), 0444);
	param_ptr = (struct SHM_PARAM *)shmat(shrd_param_id, NULL, 0);
	init_segdata_ptr = (float *)shmat(param_ptr->shrd_seg_id, NULL, SHM_RDONLY);
	xspec_ptr   = (float *)shmat(param_ptr->shrd_xspec_id, NULL, 0);
	k5head_ptr = (unsigned char *)shmat(param_ptr->shrd_k5head_id, NULL, SHM_RDONLY);
//------------------------------------------ K5 Header and Data
	setvbuf(stdout, (char *)NULL, _IONBF, 0);   // Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- Wait for the first half in the S-part
		sops.sem_num = (ushort)SEM_SEG_F; sops.sem_op = (short)-4; sops.sem_flg = (short)0;
        semop( param_ptr->sem_data_id, &sops, 1);
		segdata_ptr = init_segdata_ptr;								// Move to the head of segmant data

		//-------- FFT Time -> Freq for the first half
		cudaMemcpy(
			cuRealData,												// Segmant Data in GPU
			segdata_ptr,											// Segmant Data in Host
			Nif* NsegSec2* SegLen * sizeof(cufftReal),				// Number of segments
			cudaMemcpyHostToDevice);

		cufftExecR2C(cufft_plan, cuRealData, cuSpecData);			// FFT Time -> Freq

		//---- Auto Corr
		Dg.x=NFFT2/512; Dg.y=1; Dg.z=1;
		cudaMemset( cuPowerSpec, 0, Nif* NFFT2* sizeof(float));		// Clear Power Spec
		for(if_index=0; if_index<Nif; if_index++){
			accumPowerSpec<<<Dg, Db>>>(
				&cuSpecData[if_index* NsegSec2* NFFTC],
				&cuPowerSpec[if_index* NFFT2],
				NFFTC, NsegSec2);
		}

		//---- Cross Corr
		Dg.x=NFFT2/512; Dg.y=1; Dg.z=1;
		cudaMemset( cuCrossSpec, 0, Nif/2* NFFT2* sizeof(float2));	// Clear Power Spec

		//---- Cross Corr for IF0 * IF2
		accumCrossSpec<<<Dg, Db>>>(
			(float2 *)cuSpecData,								// FFTed Spectra (IF0)
			(float2 *)&cuSpecData[2* NsegSec2* NFFTC],			// FFTed Spectra (IF2)
			cuCrossSpec,										// Cross Spectra
			NFFTC, NsegSec2);

		//---- Cross Corr for IF1 * IF3
		accumCrossSpec<<<Dg, Db>>>(
			(float2 *)&cuSpecData[NsegSec2* NFFTC],				// FFTed Spectra (IF1)
			(float2 *)&cuSpecData[3* NsegSec2* NFFTC],			// FFTed Spectra (IF3)
			&cuCrossSpec[NFFT2],								// Cross Spectra
			NFFTC, NsegSec2);


		//-------- Wait for the Last half in the S-part
		segdata_ptr += Nif* NsegSec2* SegLen;
		sops.sem_num = (ushort)SEM_SEG_L; sops.sem_op = (short)-4; sops.sem_flg = (short)0;
        semop( param_ptr->sem_data_id, &sops, 1);

		//-------- FFT Time -> Freq for the last half
		cudaMemcpy(
			cuRealData,											// Segmant Data in GPU
			segdata_ptr,										// Segmant Data in Host
			Nif* NsegSec2* SegLen * sizeof(cufftReal),			// Number of segments
			cudaMemcpyHostToDevice);

		cufftExecR2C(cufft_plan, cuRealData, cuSpecData);		// FFT Time -> Freq

		//---- Auto Corr
		Dg.x=NFFT2/512; Dg.y=1; Dg.z=1;
		for(if_index=0; if_index<Nif; if_index++){
			accumPowerSpec<<<Dg, Db>>>(
				&cuSpecData[if_index* NsegSec2* NFFTC],
				&cuPowerSpec[if_index* NFFT2],
				NFFTC, NsegSec2);
		}

		//---- Cross Corr
		Dg.x=NFFT2/512; Dg.y=1; Dg.z=1;

		//---- Cross Corr for IF0 * IF2
		accumCrossSpec<<<Dg, Db>>>(
			(float2 *)cuSpecData,								// FFTed Spectra (IF0)
			(float2 *)&cuSpecData[2* NsegSec2* NFFTC],			// FFTed Spectra (IF2)
			cuCrossSpec,										// Cross Spectra
			NFFTC, NsegSec2);

		//---- Cross Corr for IF1 * IF3
		accumCrossSpec<<<Dg, Db>>>(
			(float2 *)&cuSpecData[NsegSec2* NFFTC],				// FFTed Spectra (IF1)
			(float2 *)&cuSpecData[3* NsegSec2* NFFTC],			// FFTed Spectra (IF3)
			&cuCrossSpec[NFFT2],								// Cross Spectra
			NFFTC, NsegSec2);


		//-------- Dump cross spectra to shared memory
		cudaMemcpy( xspec_ptr, cuPowerSpec, Nif* NFFT2* sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy( &xspec_ptr[Nif* NFFT2], cuCrossSpec, Nif/2* NFFT2* sizeof(float2), cudaMemcpyDeviceToHost);
		fwrite(xspec_ptr, sizeof(float), 2* Nif* NFFT2, file_ptr);

		sops.sem_num = (ushort)SEM_FX; sops.sem_op = (short)1; sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);
		memcpy(&sod, &k5head_ptr[4], 2);
		sod |= (k5head_ptr[6] & 0x01) << 17; sod2hms(sod, &hour, &min, &sec);
		printf("UT=%02d:%02d:%02d -- Succeeded\n", hour, min, sec);

	}
/*
-------------------------------------------- RELEASE the SHM
*/
	cufftDestroy(cufft_plan);
	cudaFree(cuRealData);
	cudaFree(cuSpecData);
	cudaFree(cuPowerSpec);
	cudaFree(cuCrossSpec);
	fclose(file_ptr);

    return(0);
}
