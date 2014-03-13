//	cuda_fft_xspec.c : FFT using CuFFT
//
//	Author : Seiji Kameno
//	Created: 2012/12/6
//
#include <cuda.h>
#include <cufft.h>
#include <string.h>
#include <math.h>
#include </usr/local/cuda-5.5/samples/common/inc/timer.h>
#include "cuda_polaris.inc"
#define	PARTNUM 2
#define SCALEFACT 2.0/(NFFT* NsegSec* PARTNUM)

extern int prob4bit( double *,	double *);
extern int initGauss4bit( double *,	double *);
extern int gaussBit(int, unsigned int *, double *, double *);
extern int sod2hms(int, int *, int *, int *);
extern int k5utc(unsigned char *,	struct SHM_PARAM *);
extern int fileRecOpen(struct SHM_PARAM *, int, int, char *, char *, FILE **);

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	int		index;						// General Index
	int		part_index;					// First and Last Part
	int		seg_index;					// Index for Segment
	int		offset[128];
	int		sod = 0;						// Seconds of Day
	int		sample_addr;
	int		bitmask = 0x000f;
	unsigned char		*k5head_ptr;	// Pointer to the K5 header
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data access
	short	*k5data_ptr;				// Pointer to shared K5 data
	float	*xspec_ptr;					// Pointer to 1-sec-integrated Power Spectrum
	FILE	*file_ptr[6];				// File Pointer to write
	FILE	*power_ptr[4];				// Power File Pointer to write
	char	fname_pre[16];
	unsigned int		bitDist[64];
	double	param[2], param_err[2];		// Gaussian parameters derived from bit distribution

	dim3			Dg, Db(512,1, 1);	// Grid and Block size
	short			*cuk5data_ptr;		// Pointer to K5 data
	cufftHandle		cufft_plan;			// 1-D FFT Plan, to be used in cufft
	cufftReal		*cuRealData;		// Time-beased data before FFT, every IF, every segment
	cufftComplex	*cuSpecData;		// FFTed spectrum, every IF, every segment
	float			*cuPowerSpec;		// (autocorrelation) Power Spectrum
	float2			*cuXSpec;

//------------------------------------------ Access to the SHARED MEMORY
	shrd_param_id = shmget( SHM_PARAM_KEY, sizeof(struct SHM_PARAM), 0444);
	param_ptr  = (struct SHM_PARAM *)shmat(shrd_param_id, NULL, 0);
	k5data_ptr = (short *)shmat(param_ptr->shrd_k5data_id, NULL, SHM_RDONLY);
	xspec_ptr  = (float *)shmat(param_ptr->shrd_xspec_id, NULL, 0);
	k5head_ptr = (unsigned char *)shmat(param_ptr->shrd_k5head_id, NULL, SHM_RDONLY);
//------------------------------------------ Prepare for CuFFT
	cudaMalloc( (void **)&cuk5data_ptr, MAX_SAMPLE_BUF);
	cudaMalloc( (void **)&cuRealData, Nif* NsegSec2* NFFT * sizeof(cufftReal) );
	cudaMalloc( (void **)&cuSpecData, Nif* NsegSec2* NFFTC* sizeof(cufftComplex) );
	cudaMalloc( (void **)&cuPowerSpec, Nif* NFFT2* sizeof(float));
	cudaMalloc( (void **)&cuXSpec, 2* NFFT2* sizeof(float2));

	if(cudaGetLastError() != cudaSuccess){
		fprintf(stderr, "Cuda Error : Failed to allocate memory.\n"); return(-1);
	}

	if(cufftPlan1d(&cufft_plan, NFFT, CUFFT_R2C, Nif* NsegSec2 ) != CUFFT_SUCCESS){
		fprintf(stderr, "Cuda Error : Failed to create plan.\n"); return(-1);
	}
//------------------------------------------ Parameters for S-part format
	for(seg_index = 0; seg_index < NsegSec; seg_index ++){
		offset[seg_index] = seg_index* (param_ptr->fsample - param_ptr->segLen) / (NsegSec - 1);
	}
//------------------------------------------ K5 Header and Data
	param_ptr->current_rec = 0;
	setvbuf(stdout, (char *)NULL, _IONBF, 0);   // Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }
		cudaMemset( cuPowerSpec, 0, Nif* NFFT2* sizeof(float));		// Clear Power Spec
		cudaMemset( cuXSpec, 0, 2* NFFT2* sizeof(float2));		// Clear Power Spec

		//-------- UTC in the K5 header
		while(k5utc(k5head_ptr, param_ptr) == 0){	usleep(1000000);}

		//-------- Open output files
		if(param_ptr->current_rec == 0){
			sprintf(fname_pre, "%04d%03d%02d%02d%02d", param_ptr->year, param_ptr->doy, param_ptr->hour, param_ptr->min, param_ptr->sec );
			for(index=0; index<Nif; index++){
				fileRecOpen(param_ptr, index, (A00_REC << index), fname_pre, "A", file_ptr);		// Autocorr
				fileRecOpen(param_ptr, index, (P00_REC << index), fname_pre, "P", power_ptr);		// Bitpower
			}
			for(index=0; index<Nif/2; index++){
				fileRecOpen(param_ptr, index, (C00_REC << index), fname_pre, "C", &file_ptr[Nif]);	// Crosscorr
			}
		}

		memset(bitDist, 0, sizeof(bitDist));
		for(part_index=0; part_index<PARTNUM; part_index ++){
			//-------- Wait for the first half in the S-part
			sops.sem_num = (ushort)(4* part_index); sops.sem_op = (short)-1; sops.sem_flg = (short)0;
			semop( param_ptr->sem_data_id, &sops, 1);

			//-------- Segment data format
			StartTimer();
			cudaMemcpy(
				&cuk5data_ptr[part_index* HALFSEC_OFFSET],
				&k5data_ptr[part_index* HALFSEC_OFFSET],
				MAX_SAMPLE_BUF/2,
				cudaMemcpyHostToDevice);

			//-------- FFT Real -> Complex spectrum
			cudaThreadSynchronize();
			Dg.x=NFFT/512; Dg.y=1; Dg.z=1;
			for(index=0; index < NsegSec2; index ++){
				seg_index = part_index* NsegSec2 + index;
				segform<<<Dg, Db>>>(
					&cuk5data_ptr[offset[seg_index]],
					&cuRealData[index* Nif* NFFT],
					NFFT);
			}

			//-------- Bit Distribution
			for(index=0; index<HALFSEC_OFFSET; index++){
				sample_addr = part_index* HALFSEC_OFFSET + index;
				bitDist[     ((k5data_ptr[sample_addr]      ) & bitmask)] ++;	// IF-0 bitdist
				bitDist[16 + ((k5data_ptr[sample_addr] >>  4) & bitmask)] ++;	// IF-1 bitdist
				bitDist[32 + ((k5data_ptr[sample_addr] >>  8) & bitmask)] ++;	// IF-2 bitdist
				bitDist[48 + ((k5data_ptr[sample_addr] >> 12) & bitmask)] ++;	// IF-3 bitdist
			}

			cudaThreadSynchronize();
			cufftExecR2C(cufft_plan, cuRealData, cuSpecData);			// FFT Time -> Freq
			cudaThreadSynchronize();

			//---- Auto Corr
			Dg.x= NFFTC/512; Dg.y=1; Dg.z=1;
			for(seg_index=0; seg_index<NsegSec2; seg_index++){
				for(index=0; index<Nif; index++){
					accumPowerSpec<<<Dg, Db>>>(
						&cuSpecData[(seg_index* Nif + index)* NFFTC],
						&cuPowerSpec[index* NFFT2],  NFFT2);
				}
			}
			//---- Cross Corr
			for(seg_index=0; seg_index<NsegSec2; seg_index++){
				accumCrossSpec<<<Dg, Db>>>(
					&cuSpecData[(seg_index* Nif)* NFFTC],
					&cuSpecData[(seg_index* Nif + 2)* NFFTC],
					cuXSpec, NFFT2);
				accumCrossSpec<<<Dg, Db>>>(
					&cuSpecData[(seg_index* Nif + 1)*NFFTC],
					&cuSpecData[(seg_index* Nif + 3)*NFFTC],
					&cuXSpec[NFFT2], NFFT2);
			}
			printf("%lf [msec]\n", GetTimer());
			
		}	// End of part loop
		Dg.x = Nif* NFFT2/512; Dg.y=1; Dg.z=1;
		scalePowerSpec<<<Dg, Db>>>(cuPowerSpec, SCALEFACT, Nif* NFFT2);
		scaleCrossSpec<<<Dg, Db>>>(cuXSpec, SCALEFACT, 2* NFFT2);

		//-------- Dump cross spectra to shared memory
		cudaMemcpy(xspec_ptr, cuPowerSpec, Nif* NFFT2* sizeof(float), cudaMemcpyDeviceToHost);
		for(index=0; index<Nif; index++){
			if(file_ptr[index] != NULL){fwrite(&xspec_ptr[index* NFFT2], sizeof(float), NFFT2, file_ptr[index]);}	// Save Pspec
			if(power_ptr[index] != NULL){fwrite(&bitDist[index* 16], sizeof(int), 16, power_ptr[index]);}			// Save Bitdist
			//-------- Total Power calculation
			gaussBit( 0x01<<(param_ptr->qbit), &bitDist[index*16], param, param_err );
			param_ptr->power[index] = 1.0/(param[0]* param[0]);
		}
		cudaMemcpy(&xspec_ptr[4* NFFT2], cuXSpec, 2* NFFT2* sizeof(float2), cudaMemcpyDeviceToHost);
		for(index=0; index<Nif/2; index++){
			if(file_ptr[Nif + index] != NULL){
				fwrite(&xspec_ptr[(Nif + index * 2)* NFFT2], sizeof(float2), NFFT2, file_ptr[Nif + index]);	// Save Xspec
			}
		}

		//-------- Refresh output data file
		if(param_ptr->current_rec == MAX_FILE_REC - 1){
			for(index=0; index<Nif+2; index++){ if( file_ptr[index] != NULL){	fclose(file_ptr[index]);} }
			for(index=0; index<Nif; index++){ if( power_ptr[index] != NULL){	fclose(power_ptr[index]);} }
			param_ptr->current_rec = 0;
		} else { param_ptr->current_rec ++; }

		sops.sem_num = (ushort)SEM_FX; sops.sem_op = (short)1; sops.sem_flg = (short)0; semop( param_ptr->sem_data_id, &sops, 1);
		sops.sem_num = (ushort)SEM_POWER; sops.sem_op = (short)1; sops.sem_flg = (short)0; semop( param_ptr->sem_data_id, &sops, 1);
		printf("%04d %03d SOD=%d UT=%02d:%02d:%02d Rec %d / %d -- Succeeded.\n",
			param_ptr->year, param_ptr->doy, sod, param_ptr->hour, param_ptr->min, param_ptr->sec, param_ptr->current_rec, param_ptr->integ_rec);
	}	// End of 1-sec loop
/*
-------------------------------------------- RELEASE the SHM
*/
	for(index=0; index<Nif+2; index++){ if( file_ptr[index] != NULL){	fclose(file_ptr[index]);} }
	for(index=0; index<Nif; index++){ if( power_ptr[index] != NULL){	fclose(power_ptr[index]);} }
	cufftDestroy(cufft_plan);
	cudaFree(cuk5data_ptr); cudaFree(cuRealData); cudaFree(cuSpecData); cudaFree(cuPowerSpec), cudaFree(cuXSpec);

    return(0);
}

