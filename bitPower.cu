//	bitPower.cu: Functions to estimate total power from bit distribution
//
//	Author : Seiji Kameno
//	Created: 2012/12/6
//
#include "shm_k5data.inc"
#include <stdio.h>
#include <string.h>
#include <math.h>
#define MAX_LEVEL	256     // Maximum number of digitized levels
#define MAX_LOOP    10      // Maximum number of iterations
#define MAX(a,b)    a>b?a:b // Larger Value

//-------- Expected probabilities in 4-bit 16-level
int probBit(
	int		nlevel,	// IN: Number of quantization levels
	double *param,	// IN: Gaussian mean and sigma
	double *prob)	// OUT:Probabilities in 16 levels
{
	int		index;	// General purpose index
	double	volt[MAX_LEVEL - 1];

	for(index = 0; index < (nlevel - 1); index ++){ volt[index] = param[0]* (double)(index - nlevel/2 + 1);}	// scaled thresh 
	//-------- Calculate probabilities
	prob[0] = 0.5* (erf(M_SQRT1_2*(volt[0] - param[1])) + 1.0);
	for(index = 1; index < (nlevel - 1); index ++){
		prob[index] = 0.5*(erf(M_SQRT1_2*(volt[index] - param[1])) - erf(M_SQRT1_2*(volt[index-1] - param[1])));
	}
	prob[nlevel-1] = 0.5* (1.0 - erf(M_SQRT1_2*(volt[nlevel-2] - param[1])));
	return(0);
}

//-------- Guess initial parameters of Gaussian distribution
int initGaussBit(
	int		nlevel,		// IN: Number of quantization levels
	double	*prob,		// IN: Probabilities in 16 levels
	double	*param)		// OUT:Estimated parameters
{
	double	Vweight;		// Weight for voltage level
	double	Average=0.0;
	double	Variance=0.0;
	int		index;			// General purpose index

	for(index=0; index<nlevel; index++){
		Vweight = (double)(index  - nlevel/2) + 0.5;
		Average += Vweight* prob[index];
		Variance += Vweight* Vweight* prob[index];
	}
	param[0] = 1.0/sqrt(Variance); param[1] = Average* param[0];
	return(0);
}

//-------- Estimate power and bias using bit distribution.
int gaussBit(
	int		nlevel,			// IN : Number of quantization levels
	unsigned int *nsample,	// IN : number of samples in each level
	double	*param,			// OUT: Gaussian parameters 
	double	*param_err)		// OUT: Gaussian parameters 
{
	int		index;					// General index for loops
	int		loop_counter = 0;		// Loop Counter
	unsigned int	total_sample = 0;	// Total number of samples
	double	pwp[2][2];				// Weighted Partial Matrix
	double	prob[MAX_LEVEL];		// Probability in each state
	double	pred[MAX_LEVEL];		// Predicted probability in each state
	double	weight[MAX_LEVEL];		// Weight for Probability 
	double	resid[MAX_LEVEL];		// residual from trial Gaussian
	double	erfDeriv[MAX_LEVEL];	// Vector to produce partial matrix
	double	WerfDeriv[MAX_LEVEL];	// Vector to produce partial matrix
	double	wpr[2];					// WPr vector
	double	solution[2];			// correction vector for parameters
	double	expArg;					// 
	double	det;					// determinant of the partial matrix
	double	norm;					// Norm of the correction vector
	double	epsz;					// criterion for convergence

	//-------- Calculate probability in each state
	for(index=0; index<nlevel; index++){ total_sample += nsample[index]; }	
	for(index=0; index<nlevel; index++){ prob[index] = (double)nsample[index] / (double)total_sample; }	
	for(index=0; index<nlevel; index++){ weight[index] = (double)nsample[index] / ((1.0 - prob[index])* (1.0 - prob[index]))  ; }	
	epsz = MAX(1.0e-6 / (total_sample* total_sample), 1.0e-29);		// Convergence

	initGaussBit(nlevel, prob, param);	// Initial parameter

	while(1){				// Loop for Least-Square Fit
		//-------- Calculate Residual Probability
		probBit(nlevel, param, pred);
		for(index=0; index<nlevel; index++){
			resid[index] = prob[index] - pred[index];
		}

		//-------- Calculate Elements of partial matrix
		erfDeriv[0] = 0.0; WerfDeriv[0] = 0.0;
		for(index=1; index<nlevel; index++){
			expArg = ((double)(index - nlevel/2))* param[0] - param[1];
			erfDeriv[index] = exp( -0.5* expArg* expArg);
			WerfDeriv[index] = ((double)(index - nlevel/2))* erfDeriv[index];
		}
		for(index=0; index<(nlevel-1); index++){
			 erfDeriv[index] = 0.5* M_2_SQRTPI* M_SQRT1_2*( -erfDeriv[index + 1] +  erfDeriv[index]);
			WerfDeriv[index] = 0.5* M_2_SQRTPI* M_SQRT1_2*( WerfDeriv[index + 1] - WerfDeriv[index]);
		}
		erfDeriv[nlevel-1] = 0.5* M_2_SQRTPI* M_SQRT1_2* erfDeriv[nlevel-1];
		WerfDeriv[nlevel-1] = -0.5* M_2_SQRTPI* M_SQRT1_2* WerfDeriv[nlevel-1];

		//-------- Partial Matrix
		memset(pwp, 0, sizeof(pwp)); memset(wpr, 0, sizeof(wpr));
		for(index=0; index<nlevel; index++){
			pwp[0][0] += (WerfDeriv[index]* WerfDeriv[index]* weight[index]);
			pwp[0][1] += (WerfDeriv[index]*  erfDeriv[index]* weight[index]);
			pwp[1][1] += ( erfDeriv[index]*  erfDeriv[index]* weight[index]);
			wpr[0] += (weight[index]* WerfDeriv[index]* resid[index]);
			wpr[1] += (weight[index]*  erfDeriv[index]* resid[index]);
		}
		pwp[1][0] = pwp[0][1];

		//-------- Solutions for correction vectors
		det = pwp[0][0]* pwp[1][1] - pwp[1][0]* pwp[0][1];
		if( fabs(det) < epsz ){	return(-1);	}						// Too small determinant -> Error
		solution[0] = (pwp[1][1]* wpr[0] - pwp[0][1]* wpr[1])/ det;
		solution[1] =(-pwp[1][0]* wpr[0] + pwp[0][0]* wpr[1])/ det;

		//-------- Correction
		param[0] += solution[0];	param[1] += solution[1];	norm = solution[0]*solution[0] + solution[1]*solution[1];

		//-------- Converged?
		loop_counter ++;
		if( norm < epsz ){	break;	}
		if( loop_counter > MAX_LOOP ){	return(-1);	}		// Doesn't converge
	}	// End of iteration loop

	//-------- Standard Error
	param_err[0] = sqrt(pwp[1][1] / det);
	param_err[1] = sqrt(pwp[0][0] / det);
	return(loop_counter);
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

//-------- UTC in the K5 header
int	k5utc(
	unsigned char		*k5head_ptr,	// IN: K5 header from K5/VSSP32
	struct SHM_PARAM	*param_ptr)		// OUT: UTC will be set in param_ptr
{
	int	sod = 0;		// Second of Day

	memcpy(&sod, &k5head_ptr[4], 2);
	if( (sod |= ((k5head_ptr[6] & 0x01) << 16)) == 0){	return(0);}

	sod2hms(sod, &(param_ptr->hour), &(param_ptr->min), &(param_ptr->sec));
	param_ptr->doy  =  k5head_ptr[8] | ((k5head_ptr[9] & 0x01) << 8);
	param_ptr->year = 2000 + ((k5head_ptr[9] >> 1) & 0x3f);
	return(sod);
}

//-------- Open Files to Record Data
int	fileRecOpen(
	struct SHM_PARAM	*param_ptr,		// IN: Shared Parameter
	int					file_index,		// IN: File index number
	int					file_flag,		// IN: File flag (A00_REC - C01_REC)
	char				*fname_pre,		// IN: File name prefix
	char				*fileType,		// IN: File type A/C/P
	FILE				**file_ptr)		//OUT: file pointer
{
	char	fname[24];
	if(param_ptr->validity & file_flag){
		sprintf(fname, "%s.%s.%02d", fname_pre, fileType, file_index);
		file_ptr[file_index] = fopen(fname, "w");
		fwrite( param_ptr, sizeof(SHM_PARAM), 1, file_ptr[file_index]);
	} else { file_ptr[file_index] = NULL;}
	return(0);
}
//-------- 4-Bit Distribution Counter
int bitDist4(
	int				nbytes,		// Number of bytes to examine
	unsigned char	*data_ptr,	// 4-bit quantized data stream (4 IF)
	unsigned int	*bitDist)	// Bit distribution counter	(4 IF x 16 levels)
{
	int	bitmask = 0x0f;			// 4-bit mask
	int	nlevel  = 16;			// Number of levels
	int index;					// Counter

	for(index=0; index<nbytes; index+=2){
		bitDist[             ((data_ptr[index  ]     ) & bitmask)] ++;	// IF-0 bitdist
		bitDist[    nlevel + ((data_ptr[index  ] >> 4) & bitmask)] ++;	// IF-1 bitdist
		bitDist[ 2* nlevel + ((data_ptr[index+1]     ) & bitmask)] ++;	// IF-2 bitdist
		bitDist[ 3* nlevel + ((data_ptr[index+1] >> 4) & bitmask)] ++;	// IF-3 bitdist
	}
	return(nbytes);
}

//-------- 8-Bit Distribution Counter
int bitDist8(
	int				nbytes,		// Number of bytes to examine
	unsigned char	*data_ptr,	// 4-bit quantized data stream (4 IF)
	unsigned int	*bitDist)	// Bit distribution counter	(4 IF x 256 levels)
{
//	int	bitmask = 0xff;			// 8-bit mask
	int	nlevel  = 256;			// Number of levels
	int index;					// Counter

	for(index=0; index<nbytes; index+=4){
		bitDist[             data_ptr[index  ] ] ++;	// IF-0 bitdist
		bitDist[    nlevel + data_ptr[index+1] ] ++;	// IF-1 bitdist
		bitDist[ 2* nlevel + data_ptr[index+2] ] ++;	// IF-2 bitdist
		bitDist[ 3* nlevel + data_ptr[index+3] ] ++;	// IF-3 bitdist
	}
//	for(index=0; index<nlevel; index++){
//		printf("%d, ", bitDist[index]);
//	}
//	printf("\n");
	return(nbytes);
}

//-------- Offset to the pointer of  segmant
int	segment_offset(
	struct SHM_PARAM	*param_ptr,	// Pointer to shared parameter
	int					*offset_ptr)
{
	int			seg_index;		// Index for segments
	long long	seg_offset;		// Offset to the segment head

	//-------- First Half
	for(seg_index = 0; seg_index < NsegSec2; seg_index ++){
		seg_offset = (long long)seg_index* ((long long)param_ptr->fsample/2 - (long long)param_ptr->segLen);
		seg_offset /= (long long)(NsegSec2 - 1);
		offset_ptr[seg_index] = (int)seg_offset;
	}

	//-------- Last Half
	for(seg_index = NsegSec2; seg_index < NsegSec; seg_index ++){
		seg_offset = (long long)(seg_index - 1)* ((long long)param_ptr->fsample/2 - (long long)param_ptr->segLen/2);
		seg_offset /= (long long)(NsegSec2 - 1);
		offset_ptr[seg_index] = (int)seg_offset;
	}

	return(NsegSec);
}

