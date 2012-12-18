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
	float	*vec_out,
	int		length)
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


main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	int		if_index;					// Index for IF
	int		seg_index;					// Index for Segment
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data access
	float	*segdata_ptr;				// Pointer to sampled data
	float	*xspec_ptr;					// Pointer to 1-sec-integrated Power Spectrum

	cufftHandle		cufft_plan;
	cufftComplex	*cuSpec_data;
	cufftReal		*cuReal_data;
	cufftReal		*cuPowerSpecSeg;
	cufftComplex	*cuCrossSpecSeg;
	cufftReal		*cuPowerSpec;
	cufftComplex	*cuCrossSpec;
//------------------------------------------ Prepare for CuFFT
	cudaMalloc( (void **)&cuReal_data, MAX_IF* Nseg_halfsec* MAX_seg_len * sizeof(cufftReal) );
	cudaMalloc( (void **)&cuSpec_data, MAX_IF* Nseg_halfsec* NFFTC* sizeof(cufftComplex) );

	cudaMalloc( (void **)&cuPowerSpecSeg, MAX_IF* Nseg_halfsec* NFFTC* sizeof(float) );
	cudaMalloc( (void **)&cuCrossSpecSeg, MAX_IF/2* Nseg_halfsec* NFFTC* sizeof(float2) );

	cudaMalloc( (void **)&cuPowerSpec, MAX_IF* NFFTC* sizeof(float));
	cudaMalloc( (void **)&cuCrossSpec, MAX_IF/2* NFFTC* sizeof(float2));
	if(cudaGetLastError() != cudaSuccess){
		fprintf(stderr, "Cuda Error : Failed to allocate memory.\n"); return(-1);
	}

	if(cufftPlan1d(&cufft_plan, NFFT, CUFFT_R2C, Nseg_halfsec* MAX_IF ) != CUFFT_SUCCESS){
		fprintf(stderr, "Cuda Error : Failed to create plan.\n"); return(-1);
	}
//------------------------------------------ Access to the SHARED MEMORY
	shrd_param_id = shmget( SHM_PARAM_KEY, sizeof(struct SHM_PARAM), 0444);
	param_ptr = (struct SHM_PARAM *)shmat(shrd_param_id, NULL, 0);
	segdata_ptr = (float *)shmat(param_ptr->shrd_seg_id, NULL, SHM_RDONLY);
	xspec_ptr   = (float *)shmat(param_ptr->shrd_xspec_id, NULL, 0);
//------------------------------------------ K5 Header and Data
	setvbuf(stdout, (char *)NULL, _IONBF, 0);   // Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- Wait for the first half in the S-part
		printf("FFT status EW\n");
		sops.sem_num = (ushort)8; sops.sem_op = (short)-4; sops.sem_flg = (short)0;
        semop( param_ptr->sem_data_id, &sops, 1);

		//-------- FFT Time -> Freq for the first half
		printf("FFT status ES\n");
		cudaMemcpy(
			cuReal_data,											// Segmant Data in GPU
			segdata_ptr,											// Segmant Data in Host
			MAX_IF* Nseg_halfsec* MAX_seg_len * sizeof(cufftReal),	// Number of segments
			cudaMemcpyHostToDevice);

		cufftExecR2C(cufft_plan, cuReal_data, cuSpec_data);			// FFT Time -> Freq

		//---- Auto Corr
		complexPowerVec<<<NFFTC,  Nseg_halfsec* MAX_IF>>>(
			cuSpec_data,								// FFTed Spectra
			cuPowerSpecSeg,							// Power Spectra
			NFFTC* Nseg_halfsec* MAX_IF);

		//---- Cross Corr for IF0 * IF2
		complexMultConjVec<<<NFFTC, Nseg_halfsec>>>(
			(float2 *)cuSpec_data,								// FFTed Spectra (IF0)
			(float2 *)&cuSpec_data[2* NFFTC* Nseg_halfsec],		// FFTed Spectra (IF2)
			cuCrossSpecSeg,										// Cross Spectra
			NFFTC* Nseg_halfsec);

		//---- Cross Corr for IF1 * IF3
		complexMultConjVec<<<NFFTC, Nseg_halfsec>>>(
			(float2 *)&cuSpec_data[3* NFFTC* Nseg_halfsec],		// FFTed Spectra (IF1)
			(float2 *)&cuSpec_data[4* NFFTC* Nseg_halfsec],		// FFTed Spectra (IF3)
			&cuCrossSpecSeg[Nseg_halfsec* NFFTC],				// Cross Spectra
			NFFTC* Nseg_halfsec);

		printf("FFT status EF\n");

		//---- Accumulate power spectra
		cudaMemset( cuPowerSpec, 0, MAX_IF* NFFTC* sizeof(float));	// Clear Power Spec
		for(if_index=0; if_index<MAX_IF; if_index ++){
			for(seg_index=0; seg_index<Nseg_halfsec; seg_index++){
				accumReal<<<NFFTC, 1>>>(
					cuPowerSpec,								// Accumulator Register
					&cuPowerSpecSeg[(if_index* Nseg_halfsec +  seg_index)* NFFTC],
					NFFTC); 
			}
		}

		cudaMemcpy( xspec_ptr, cuPowerSpec, NFFTC* sizeof(float), cudaMemcpyDeviceToHost);
		printf("%e %e\n", xspec_ptr[0], xspec_ptr[1]);


#ifdef HIDOI
		//---- Accumulate Cross spectra
		cudaMemset( cuCrossSpec, 0, MAX_IF/2* NFFTC* sizeof(float2));	// Clear Cross Spec
		for(seg_index=0; seg_index<Nseg_halfsec; seg_index++){
			accumComplex<<<NFFTC, 1>>>(
				cuCrossSpec,
				&cuCrossSpecSeg[seg_index* NFFTC],
				NFFTC); 

			accumComplex<<<NFFTC, 1>>>(
				&cuCrossSpec[NFFTC],
				&cuCrossSpecSeg[(Nseg_halfsec + seg_index)* NFFTC],
				NFFTC); 
		}
#endif

		//-------- Latter Half --------
		segdata_ptr += MAX_IF* MAX_seg_len* Nseg_halfsec;
		//-------- Wait for the first half in the S-part
		printf("FFT status OW\n");
		sops.sem_num = (ushort)9; sops.sem_op = (short)-4; sops.sem_flg = (short)0;
        semop( param_ptr->sem_data_id, &sops, 1);

		printf("FFT status OS\n");
		cudaMemcpy(
			cuReal_data,
			segdata_ptr,
			MAX_IF* Nseg_halfsec* MAX_seg_len * sizeof(cufftReal),
			cudaMemcpyHostToDevice);
		segdata_ptr -= MAX_IF* MAX_seg_len* Nseg_halfsec;

		cufftExecR2C(cufft_plan, cuReal_data, cuSpec_data);			// FFT Time -> Freq

#ifdef HIDOI
		//---- Auto Corr
		complexPowerVec<<<NFFTC,  MAX_IF* Nseg_halfsec>>>(
			(float2 *)cuSpec_data,								// FFTed Spectra
			cuPowerSpecSeg,										// Power Spectra
			NFFTC* Nseg_halfsec* MAX_IF);

		printf("Length=%d\n", NFFTC* Nseg_halfsec* MAX_IF);

		//---- Cross Corr for IF0 * IF2
		complexMultConjVec<<<NFFTC, Nseg_halfsec>>>(
			(float2 *)cuSpec_data,								// FFTed Spectra (IF0)
			(float2 *)&cuSpec_data[2* NFFTC* Nseg_halfsec],		// FFTed Spectra (IF2)
			(float2 *)cuCrossSpecSeg,							// Cross Spectra
			NFFTC* Nseg_halfsec);

		//---- Cross Corr for IF1 * IF3
		complexMultConjVec<<<NFFTC, Nseg_halfsec>>>(
			(float2 *)&cuSpec_data[3* NFFTC* Nseg_halfsec],		// FFTed Spectra (IF1)
			(float2 *)&cuSpec_data[4* NFFTC* Nseg_halfsec],		// FFTed Spectra (IF3)
			(float2 *)&cuCrossSpecSeg[Nseg_halfsec* NFFTC],		// Cross Spectra
			NFFTC* Nseg_halfsec);

		printf("FFT status OF\n");

		//---- Accumulate power spectra
		for(if_index=0; if_index<MAX_IF; if_index ++){
			for(seg_index=0; seg_index<Nseg_halfsec; seg_index++){
				accumReal<<<NFFTC, 1>>>(
					cuPowerSpec,								// Accumulator Register
					&cuPowerSpecSeg[(if_index* Nseg_halfsec +  seg_index)* NFFTC],
					NFFTC); 
			}
		}
		//---- Accumulate Cross spectra
		for(seg_index=0; seg_index<Nseg_halfsec; seg_index++){
			accumComplex<<<NFFTC, 1>>>(
				cuCrossSpec,
				&cuCrossSpecSeg[seg_index* NFFTC],
				NFFTC); 

			accumComplex<<<NFFTC, 1>>>(
				&cuCrossSpec[NFFTC],
				&cuCrossSpecSeg[(Nseg_halfsec + seg_index)* NFFTC],
				NFFTC); 
		}
#endif

//		cudaMemcpy( xspec_ptr, cuPowerSpec, MAX_IF* NFFTC* sizeof(float), cudaMemcpyDeviceToHost);
//		cudaMemcpy( &xspec_ptr[MAX_IF* NFFTC], cuCrossSpec, MAX_IF/2* NFFTC* sizeof(float), cudaMemcpyDeviceToHost);

		sops.sem_num = (ushort)SEM_FX; sops.sem_op = (short)1; sops.sem_flg = (short)0;
        semop( param_ptr->sem_data_id, &sops, 1);
	}
/*
-------------------------------------------- RELEASE the SHM
*/
	cufftDestroy(cufft_plan);
	cudaFree(&cuReal_data);
	cudaFree(&cuSpec_data);
	cudaFree(&cuPowerSpecSeg);
	cudaFree(&cuCrossSpecSeg);
	cudaFree(&cuPowerSpec);
	cudaFree(&cuCrossSpec);

    return(0);
}
