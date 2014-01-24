//	cuda_fft_xspec.c : FFT using CuFFT
//
//	Author : Seiji Kameno
//	Created: 2012/12/6
//
#include <cuda.h>
#include <cufft.h>
#include <string.h>
#include <math.h>
#include </usr/local/cuda-5.0/samples/common/inc/timer.h>
#include "cuda_polaris.inc"
#define	PARTNUM 2
#define SCALEFACT 1.0/(NFFT* NsegSec2* PARTNUM)
#define	MAX_LOOP	10		// Maximum number of iterations
#define	MAX(a,b)	a>b?a:b	// Larger Value

int prob4bit(
	double *param,	// IN: Gaussian mean and sigma
	double *prob)	// OUT:Probabilities in 16 levels
{
	int		index;	// General purpose index
	double	volt[] = {-7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

	for(index = 0; index < 15; index ++){ volt[index] *= param[0]; }		// Scale threshold
	//-------- Calculate probabilities
	prob[0] = 0.5* (erf(M_SQRT1_2*(volt[0] - param[1])) + 1.0);
	for(index = 1; index < 14; index ++){
		prob[index] = 0.5*(erf(M_SQRT1_2*(volt[index] - param[1])) - erf(M_SQRT1_2*(volt[index-1] - param[1])));
	}
	prob[15] = 0.5* (1.0 - erf(M_SQRT1_2*(volt[14] - param[1])));
	return(0);
}

int initGauss4bit(
	double	*prob,		// IN: Probabilities in 16 levels
	double	*param)		// OUT:Estimated parameters
{
	double	Vweight[] = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	double	Pweight[] = {56.25, 42.25, 30.25, 20.25, 12.25, 6.25, 2.25, 0.25, 0.25, 2.25, 6.25, 12.25, 20.25, 30.25, 42.25, 56.25};
	double	Average=0.0;
	double	Variance=0.0;
	int		index;			// General purpose index

	for(index=0; index<16; index++){
		Average += Vweight[index]* prob[index];
		Variance += Pweight[index]* prob[index];
	}
	param[0] = 1.0/sqrt(Variance); param[1] = Average* param[0];
	return(0);
}

int gauss4bit(
	unsigned int *nsample,	// IN : number of samples in 16 levels
	double	*param,			// OUT: Gaussian parameters 
	double	*param_err)		// OUT: Gaussian parameters 
{
	int		index;			// General index for loops
	int		loop_counter = 0;		// Loop Counter
	unsigned int	total_sample = 0;	// Total number of samples
	double	pwp[2][2];		// Weighted Partial Matrix
	double	prob[16];		// Probability in each state
	double	pred[16];		// Predicted probability in each state
	double	weight[16];		// Weight for Probability 
	double	resid[16];		// residual from trial Gaussian
	double	erfDeriv[16];	// Vector to produce partial matrix
	double	WerfDeriv[16];	// Vector to produce partial matrix
	double	wpr[2];			// WPr vector
	double	solution[2];	// correction vector for parameters
	double	expArg;			// 
	double	det;			// determinant of the partial matrix
	double	norm;			// Norm of the correction vector
	double	epsz;			// criterion for convergence

	//-------- Calculate probability in each state
	for(index=0; index<16; index++){ total_sample += nsample[index]; }	
	for(index=0; index<16; index++){ prob[index] = (double)nsample[index] / (double)total_sample; }	
	for(index=0; index<16; index++){ weight[index] = (double)nsample[index] / ((1.0 - prob[index])* (1.0 - prob[index]))  ; }	
	epsz = MAX(1.0e-6 / (total_sample* total_sample), 1.0e-29);		// Convergence

	initGauss4bit(prob, param);	// Initial parameter

	while(1){				// Loop for Least-Square Fit
		//-------- Calculate Residual Probability
		prob4bit(param, pred);
		for(index=0; index<16; index++){ resid[index] = prob[index] - pred[index]; }

		//-------- Calculate Elements of partial matrix
		erfDeriv[0] = 0.0; WerfDeriv[0] = 0.0;
		for(index=1; index<16; index++){
			expArg = ((double)index - 8.0)* param[0] - param[1];
			erfDeriv[index] = exp( -0.5* expArg* expArg);
			WerfDeriv[index] = ((double)index - 8.0)* erfDeriv[index];
		}
		for(index=0; index<15; index++){
			 erfDeriv[index] = 0.5* M_2_SQRTPI* M_SQRT1_2*( -erfDeriv[index + 1] +  erfDeriv[index]);
			WerfDeriv[index] = 0.5* M_2_SQRTPI* M_SQRT1_2*( WerfDeriv[index + 1] - WerfDeriv[index]);
		}
		erfDeriv[15] = 0.5* M_2_SQRTPI* M_SQRT1_2* erfDeriv[15];
		WerfDeriv[15] = -0.5* M_2_SQRTPI* M_SQRT1_2* WerfDeriv[15];

		//-------- Partial Matrix
		memset(pwp, 0, sizeof(pwp)); memset(wpr, 0, sizeof(wpr));
		for(index=0; index<16; index++){
			pwp[0][0] += (WerfDeriv[index]* WerfDeriv[index]* weight[index]);
			pwp[0][1] += (WerfDeriv[index]*  erfDeriv[index]* weight[index]);
			pwp[1][1] += ( erfDeriv[index]*  erfDeriv[index]* weight[index]);
			wpr[0] += (weight[index]* WerfDeriv[index]* resid[index]);
			wpr[1] += (weight[index]*  erfDeriv[index]* resid[index]);
		}
		pwp[1][0] = pwp[0][1];

		//-------- Solutions for correction vectors
		det = pwp[0][0]* pwp[1][1] - pwp[1][0]* pwp[0][1];
		if( fabs(det) < epsz ){	return(1);	}						// Too small determinant -> Error
		solution[0] = (pwp[1][1]* wpr[0] - pwp[0][1]* wpr[1])/ det;
		solution[1] =(-pwp[1][0]* wpr[0] + pwp[0][0]* wpr[1])/ det;

		//-------- Correction
		param[0] += solution[0];	param[1] += solution[1];	norm = solution[0]*solution[0] + solution[1]*solution[1];

		//-------- Converged?
		loop_counter ++;
		if( norm < epsz ){	break;	}
		if( loop_counter > MAX_LOOP ){	return(1);	}		// Doesn't converge
	}	// End of iteration loop

	//-------- Standard Error
	param_err[0] = sqrt(pwp[1][1] / det);
	param_err[1] = sqrt(pwp[0][0] / det);
	return(0);
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
	char	fname[24];					// File Name [YYYYDOYHHMMSSIF]
	char	fname_pre[16];
	unsigned int		bitDist[64];
	double	param[2], param_err[2];		// Gaussian parameters derived from bit distribution
	// float	bitPower;
	// float wt[] = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	// int		totalSample;
	// int		level_index;

	dim3			Dg, Db(512,1, 1);	// Grid and Block size
	short			*cuk5data_ptr;		// Pointer to K5 data
	cufftHandle		cufft_plan;			// 1-D FFT Plan, to be used in cufft
	cufftReal		*cuRealData;		// Time-beased data before FFT, every IF, every segment
	cufftComplex	*cuSpecData;		// FFTed spectrum, every IF, every segment
	float			*cuPowerSpec;		// (autocorrelation) Power Spectrum
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
	for(seg_index = 0; seg_index < param_ptr->segNum; seg_index ++){
		offset[seg_index] = seg_index* (param_ptr->fsample - param_ptr->segLen) / (param_ptr->segNum - 1);
	}
//------------------------------------------ K5 Header and Data
	param_ptr->current_rec = 0;
	setvbuf(stdout, (char *)NULL, _IONBF, 0);   // Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }
		cudaMemset( cuPowerSpec, 0, Nif* NFFT2* sizeof(float));		// Clear Power Spec
		cudaMemset( cuXSpec, 0, 2* NFFT2* sizeof(float2));		// Clear Power Spec

		//-------- UTC in the K5 header
		sod = 0;
		while(sod == 0){	// Wait until UTC to be valid
			usleep(100000);	// Wait 100 msec
			memcpy(&sod, &k5head_ptr[4], 2);
			sod |= ((k5head_ptr[6] & 0x01) << 16);
		}
		sod2hms(sod, &(param_ptr->hour), &(param_ptr->min), &(param_ptr->sec));
		param_ptr->doy  =  k5head_ptr[8] | ((k5head_ptr[9] & 0x01) << 8);
		param_ptr->year = 2000 + ((k5head_ptr[9] >> 1) & 0x3f);

		//-------- Open output files
		if(param_ptr->current_rec == 0){
			sprintf(fname_pre, "%04d%03d%02d%02d%02d", param_ptr->year, param_ptr->doy, param_ptr->hour, param_ptr->min, param_ptr->sec );

			for(index=0; index<Nif; index++){
				//-------- Autocorrelation File Record Switch
				if(param_ptr->validity & (A00_REC << index)){
					sprintf(fname, "%s.%s.%02d", fname_pre, "A", index);
					file_ptr[index] = fopen(fname, "w");
					fwrite( param_ptr, sizeof(SHM_PARAM), 1, file_ptr[index]);
				} else { file_ptr[index] = NULL;}

				//-------- Bit-distribution File Record Switch
				if(param_ptr->validity & (P00_REC << index)){
					sprintf(fname, "%s.%s.%02d", fname_pre, "P", index);
					power_ptr[index] = fopen(fname, "w");
					fwrite( param_ptr, sizeof(SHM_PARAM), 1, power_ptr[index]);
				} else { power_ptr[index] = NULL;}
			}
			for(index=0; index<Nif/2; index++){
				//-------- Crosscorrelation File Record Switch
				if(param_ptr->validity & (C00_REC << index)){
					sprintf(fname, "%s.%s.%02d", fname_pre, "C", index);  file_ptr[Nif + index]   = fopen(fname, "w");
					fwrite( param_ptr, sizeof(SHM_PARAM), 1, file_ptr[Nif + index]);
				} else { file_ptr[Nif + index] = NULL;}
			}
		}

		memset(bitDist, 0, sizeof(bitDist));
		for(part_index=0; part_index<PARTNUM; part_index ++){
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

			//-------- FFT Real -> Complex spectrum
			cudaThreadSynchronize();
			Dg.x=SegLen/512; Dg.y=1; Dg.z=1;
			for(index=0; index < NsegSec2; index ++){
				seg_index = part_index* NsegSec2 + index;
				segform<<<Dg, Db>>>(
					&cuk5data_ptr[offset[seg_index]],
					&cuRealData[index* Nif* SegLen],
					SegLen);
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
//			printf("%lf [msec]\n", GetTimer());
			
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
			// totalSample = 0;	bitPower = 0.0;
			// for(level_index=0; level_index<16; level_index++){
			// 	totalSample += bitDist[index* 16 + level_index];
			// 	bitPower	+= wt[level_index]* wt[level_index]* (float)bitDist[index* 16 + level_index];
			// }
			gauss4bit( &bitDist[index*16], param, param_err );
			// param_ptr->power[index] = bitPower / (float)totalSample;
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

