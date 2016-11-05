//	shm_spec_view.c : Cross Power Spectrom Viewer
//
//	Author : Seiji Kameno
//	Created: 2012/12/18
//
#include "shm_k5data.inc"
// #include "k5dict.inc"
#include <stdlib.h>
#include <cpgplot.h>
#include <math.h>

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data area
	int		ch_index;					// Index for channel
	int		bunchNum;					// Number of bunching channels	
	float	*xspec_ptr;					// Pointer to Cross Spectra
	float	*spec_ptr;					// Pointer to spectrum to plot
	float	*freq_ptr;					// Pointer to Frequency
	float	freq_incr;					// Freqnecy increment
	char	pg_text[256];				// Text to plot
	char	xlabel[64];					// X-axis label
//------------------------------------------ Access to the SHARED MEMORY
    //------- SHARED PARAMETERS --------
    if(shm_access(
        SHM_PARAM_KEY,					// ACCESS KEY
        sizeof(struct SHM_PARAM),		// SIZE OF SHM
        &shrd_param_id,					// SHM ID
        &param_ptr) != -1){				// Pointer to the SHM
		printf("Succeeded to access the shared parameter [%d]!\n",  param_ptr->shrd_param_id);
	}
	xspec_ptr = shmat(param_ptr->shrd_xspec_id, NULL, SHM_RDONLY);
//------------------------------------------ K5 Header and Data
	sleep(1);									// Timing adjustment for PGPLOT device
	setvbuf(stdout, (char *)NULL, _IONBF, 0);	// Disable stdout cache
	cpgbeg(1, argv[1], 1, 1);
	sprintf(xlabel, "Frequency [MHz]"); cpg_setup(xlabel);
	freq_ptr = (float *)malloc(MAX_CH_VIEW* sizeof(float));
	spec_ptr = (float *)malloc(NST* MAX_CH_VIEW* sizeof(float));
	freq_incr = (double)(param_ptr->fsample) * 5.0e-7 / MAX_CH_VIEW;
	for(ch_index=0; ch_index<MAX_CH_VIEW; ch_index ++){
		freq_ptr[ch_index] = (float)((ch_index - 0.5)* freq_incr);
	}
	bunchNum = NFFT2 / MAX_CH_VIEW;
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- Wait for Semaphore
		sops.sem_num = (ushort)SEM_FX;	sops.sem_op = (short)-1;	sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);

		//-------- Plot spectra using PGPLOT
		bunchVec(NST*NFFT2, bunchNum, xspec_ptr, spec_ptr);
		cpg_spec(param_ptr, freq_ptr, spec_ptr);
	}
	cpgend();
//------------------------------------------ RELEASE the SHM
	free(freq_ptr);
	free(spec_ptr);
    return(0);
}

int bunchVec(
	int		vecLen,		// Length of the input vector
	int		bunchNum,	// Number of channels to bunch
	float	*vec_ptr,	// Pointer to the input vector
	float	*out_ptr)	// Pointer to the output vector
{
	int	outNum;
	int	index=0;
	int	out_index, avg_index;

	outNum = vecLen / bunchNum;
	memset(out_ptr, 0x00, outNum* sizeof(float));
	for(out_index=0; out_index<outNum; out_index++){
		for(avg_index=0; avg_index<bunchNum; avg_index++){
			out_ptr[out_index] += vec_ptr[index];
			index++;
		}
		out_ptr[out_index] /= bunchNum;
	}
	return(outNum);
}
