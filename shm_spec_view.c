//	shm_spec_view.c : Cross Power Spectrom Viewer
//
//	Author : Seiji Kameno
//	Created: 2012/12/18
//
#include "shm_k5data.inc"
#include "k5dict.inc"
#include <stdlib.h>
#include <cpgplot.h>

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data area
	int		index;
	int		sod=0, hour=0, min=0, sec=0;// Second of day, hour, min, sec
	float	*xspec_ptr;					// Pointer to Cross Spectra
	float	*freq_ptr;					// Pointer to Frequency
	float	freq_incr;					// Freqnecy increment
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

	setvbuf(stdout, (char *)NULL, _IONBF, 0);	// Disable stdout cache
	cpgbeg(1, "/xserv", 1, 1);
	cpg_setup(param_ptr);
	freq_ptr = (float *)malloc(NFFT2* sizeof(float));
	freq_incr = (float)(param_ptr->fsample / 2 / NFFT2);
	for(index=0; index<NFFT2; index ++){
		freq_ptr[index] = -0.5* freq_incr + index* freq_incr;
	}

	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- Wait for Semaphore
		sops.sem_num = (ushort)SEM_FX;	sops.sem_op = (short)-1;	sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);

		printf("%7.2e %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e\n",
			xspec_ptr[0], xspec_ptr[1], xspec_ptr[2], xspec_ptr[3],
			xspec_ptr[4], xspec_ptr[5], xspec_ptr[6], xspec_ptr[7]);
		cpg_spec(param_ptr, freq_ptr, xspec_ptr);

	}
	cpgend();
//------------------------------------------ RELEASE the SHM
	free(freq_ptr);
    return(0);
}
