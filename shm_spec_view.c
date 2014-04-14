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
	int		index;
	float	*xspec_ptr;					// Pointer to Cross Spectra
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
	sprintf(xlabel, "Frequency [Hz]\0"); cpg_setup(xlabel);
	freq_ptr = (float *)malloc(NFFT2* sizeof(float));
	freq_incr = (double)(param_ptr->fsample) / 2 / NFFT2;
	for(index=0; index<NFFT2; index ++){
		freq_ptr[index] = (float)((index - 0.5)* freq_incr);
	}

	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- Wait for Semaphore
		sops.sem_num = (ushort)SEM_FX;	sops.sem_op = (short)-1;	sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);

		//-------- Plot spectra using PGPLOT
		cpg_spec(param_ptr, freq_ptr, xspec_ptr);
	}
	cpgend();
//------------------------------------------ RELEASE the SHM
	free(freq_ptr);
    return(0);
}
