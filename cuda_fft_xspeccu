//	cuda_fft_xspec.c : FFT using CuFFT
//
//	Author : Seiji Kameno
//	Created: 2012/12/6
//
#include <cuda.h>
#include <cufft.h>
#include </usr/local/cuda-5.0/samples/common/inc/timer.h>
#include "cuda_polaris.inc"

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	int		index;						// General Index
	int		part_index;					// First and Last Part
	int		seg_index;					// Index for Segment
	int		mean_offset, fraction, offset[128], overlap[128];
	int		sod = 0;						// Seconds of Day
	int		rec_index=0;
	unsigned char		*k5head_ptr;	// Pointer to the K5 header
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data access
	short	*k5data_ptr;				// Pointer to shared K5 data
	float	*xspec_ptr;					// Pointer to 1-sec-integrated Power Spectrum
	FILE	*file_ptr[6];				// File Pointer to write
	char	fname[24];					// File Name [YYYYDOYHHMMSSIF]
	char	fname_pre[16];

	dim3			Dg, Db(512,1, 1);	// Grid and Block size
	short			*cuk5data_ptr;		// Pointer to K5 data
	cufftHandle		cufft_plan;			// 1-D FFT Plan, to be used in cufft
	cufftReal		*cuRealData;		// Time-beased data before FFT, every IF, every segment
	cufftComplex	*cuSpecData;		// FFTed spectrum, every IF, every segment
	cufftReal		*cuPowerSpec;		// (autocorrelation) Power Spectrum
	float2			*cuXSpec;

//------------------------------------------ Prepare for CuFFT
	cudaMalloc( (void **)&cuk5data_ptr, MAX_SAMPLE_BUF);
	cudaMalloc( (void **)&cuRealData, Nif* NsegSec2* SegLen * sizeof(cufftReal) );
	cudaMalloc( (void **)&cuSpecData, Nif* NsegSec2* NFFTC* sizeof(cufftComplex) );
	cudaMalloc( (void **)&cuPowerSpec, Nif* NFFT2* sizeof(float));
	cudaMalloc( (void **)&cuXSpec, 2* NFFT2* sizeof(float2));

	if(cudaGetLastError() != cudaSuccess){
		fprintf(stderr, "Cuda Error : Failed to allocate memory.\n"); return(-1);
	}

	if(cufftPlan1d(&cufft_plan, NFFT, CUFFT_R2C, Nif* NsegSec2 ) != CUFFT_SUCCESS){
		fprintf(stderr, "Cuda Error : Failed to create plan.\n"); return(-1);
	}
//------------------------------------------ Access to the SHARED MEMORY
	shrd_param_id = shmget( SHM_PARAM_KEY, sizeof(struct SHM_PARAM), 0444);
	param_ptr  = (struct SHM_PARAM *)shmat(shrd_param_id, NULL, 0);
	k5data_ptr = (short *)shmat(param_ptr->shrd_k5data_id, NULL, SHM_RDONLY);
	xspec_ptr  = (float *)shmat(param_ptr->shrd_xspec_id, NULL, 0);
	k5head_ptr = (unsigned char *)shmat(param_ptr->shrd_k5head_id, NULL, SHM_RDONLY);
//------------------------------------------ Parameters for S-part format
	mean_offset = (param_ptr->fsample - param_ptr->segLen) / (param_ptr->segNum - 1);
	fraction = (param_ptr->fsample - param_ptr->segLen) % mean_offset;
	offset[0] = 0;  overlap[0] = 0;
	for(seg_index = 1; seg_index < param_ptr->segNum; seg_index ++){
		overlap[seg_index] = param_ptr->segLen - mean_offset;
		if( seg_index % (param_ptr->segNum / fraction) == 1 ){  overlap[seg_index] --;}
		offset[seg_index] = offset[seg_index-1] + param_ptr->segLen - overlap[seg_index];
    }
//------------------------------------------ K5 Header and Data
	setvbuf(stdout, (char *)NULL, _IONBF, 0);   // Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }
		cudaMemset( cuPowerSpec, 0, Nif* NFFT2* sizeof(float));		// Clear Power Spec
		cudaMemset( cuXSpec, 0, 2* NFFT2* sizeof(float2));		// Clear Power Spec

		//-------- UTC in the K5 header
<<<<<<< HEAD
		sod = 0;
		while(sod == 0){	// Wait until UTC to be valid
			usleep(100000);	// Wait 100 msec
			memcpy(&sod, &k5head_ptr[4], 2);
			sod |= ((k5head_ptr[6] & 0x01) << 16);
		}
=======
		sod = 0; memcpy(&sod, &k5head_ptr[4], 2);
		sod |= ((k5head_ptr[6] & 0x01) << 16);
>>>>>>> d461cc85d81d7650e34e7a394988ef9c5dc2d33a
		sod2hms(sod, &(param_ptr->hour), &(param_ptr->min), &(param_ptr->sec));
		param_ptr->doy  =  k5head_ptr[8] | ((k5head_ptr[9] & 0x01) << 8);
		param_ptr->year = 2000 + ((k5head_ptr[9] >> 1) & 0x3f);

		//-------- Open output files
		if(rec_index == 0){
			sprintf(fname_pre, "%04d%03d%02d%02d%02d", param_ptr->year, param_ptr->doy, param_ptr->hour, param_ptr->min, param_ptr->sec );
			for(index=0; index<Nif; index++){
				sprintf(fname, "%s.%s.%02d", fname_pre, "A", index);
				file_ptr[index] = fopen(fname, "w");
				fwrite( param_ptr, sizeof(SHM_PARAM), 1, file_ptr[index]);
			}
			sprintf(fname, "%s.%s.%02d", fname_pre, "C", 0);  file_ptr[Nif]   = fopen(fname, "w");
			sprintf(fname, "%s.%s.%02d", fname_pre, "C", 1);  file_ptr[Nif+1] = fopen(fname, "w");
			fwrite( param_ptr, sizeof(SHM_PARAM), 1, file_ptr[Nif]);
			fwrite( param_ptr, sizeof(SHM_PARAM), 1, file_ptr[Nif+1]);
		}

		for(part_index=0; part_index<2; part_index ++){
			//-------- Wait for the first half in the S-part
			sops.sem_num = (ushort)(4* part_index); sops.sem_op = (short)-1; sops.sem_flg = (short)0;
			semop( param_ptr->sem_data_id, &sops, 1);

			//-------- Segment data format
//			StartTimer();
			cudaMemcpy(
				&cuk5data_ptr[part_index* HALFSEC_OFFSET],
				&k5data_ptr[part_index* HALFSEC_OFFSET],
				MAX_SAMPLE_BUF/2,
				cudaMemcpyHostToDevice);
			cudaThreadSynchronize();

			//-------- FFT Real -> Complex spectrum
			Dg.x=SegLen/512; Dg.y=1; Dg.z=1;
			for(index=0; index < NsegSec2; index ++){
				seg_index = part_index* NsegSec2 + index;
				segform<<<Dg, Db>>>(
					&cuk5data_ptr[offset[seg_index]],
					&cuRealData[index* Nif* SegLen],
					SegLen);
			}
			cudaThreadSynchronize();
			cufftExecR2C(cufft_plan, cuRealData, cuSpecData);			// FFT Time -> Freq
			cudaThreadSynchronize();

			//---- Auto Corr
			Dg.x= NFFTC/512; Dg.y=1; Dg.z=1;
			for(index=0; index<Nif; index++){
				accumPowerSpec<<<Dg, Db>>>(
					&cuSpecData[index* NFFTC],
					&cuPowerSpec[index* NFFT2],  NFFT2);
			}
			//---- Cross Corr
			accumCrossSpec<<<Dg, Db>>>(cuSpecData, &cuSpecData[2*NFFTC], cuXSpec, NFFT2);
			accumCrossSpec<<<Dg, Db>>>(&cuSpecData[NFFTC], &cuSpecData[3*NFFTC], &cuXSpec[NFFT2], NFFT2);
//			printf("%lf [msec]\n", GetTimer());
			
		}

		//-------- Dump cross spectra to shared memory
		cudaMemcpy(xspec_ptr, cuPowerSpec, Nif* NFFT2* sizeof(float), cudaMemcpyDeviceToHost);
		for(index=0; index<Nif; index++){
			fwrite(&xspec_ptr[index* NFFT2], sizeof(float), NFFT2, file_ptr[index]);
		}
		cudaMemcpy(&xspec_ptr[4* NFFT2], cuXSpec, 2* NFFT2* sizeof(float2), cudaMemcpyDeviceToHost);
		fwrite(&xspec_ptr[4* NFFT2], sizeof(float2), NFFT2, file_ptr[4]);
		fwrite(&xspec_ptr[6* NFFT2], sizeof(float2), NFFT2, file_ptr[5]);

		//-------- Refresh output data file
		if(rec_index == MAX_FILE_REC - 1){
			for(index=0; index<Nif+2; index++){ fclose(file_ptr[index]); }
			rec_index = 0;
		} else { rec_index ++; }

		sops.sem_num = (ushort)SEM_FX; sops.sem_op = (short)1; sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);
		printf("%04d %03d SOD=%d UT=%02d:%02d:%02d -- Succeeded.\n",
			param_ptr->year, param_ptr->doy, sod,
			param_ptr->hour, param_ptr->min, param_ptr->sec);
	}
/*
-------------------------------------------- RELEASE the SHM
*/
	for(index=0; index<Nif+2; index++){ fclose(file_ptr[index]); }
	cufftDestroy(cufft_plan);
	cudaFree(cuk5data_ptr); cudaFree(cuRealData); cudaFree(cuSpecData); cudaFree(cuPowerSpec), cudaFree(cuXSpec);

    return(0);
}
