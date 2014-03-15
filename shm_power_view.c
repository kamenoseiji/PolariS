//	shm_power_view.c : Total Power Monitor
//
//	Author : Seiji Kameno
//	Created: 2013/12/23
//
#include "shm_k5data.inc"
// #include "k5dict.inc"
#include <stdlib.h>
#include <cpgplot.h>
#include <math.h>
#define timeNum 128

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data area
	int		IFindex;
	int		index;
	float	bitPower[Nif][timeNum];		// Power Monitor Data
	float	tempPower[timeNum];			// Buffer
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
	memset(bitPower, 0, Nif* timeNum* sizeof(float));
//------------------------------------------ K5 Header and Data
	setvbuf(stdout, (char *)NULL, _IONBF, 0);	// Disable stdout cache
	cpgbeg(1, argv[1], 1, 1);

	while(param_ptr->validity & ACTIVE){
		cpgbbuf();
		sprintf(xlabel, "Elapsed Time [sec]\0"); cpg_setup(xlabel);
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- Wait for Semaphore
		sops.sem_num = (ushort)SEM_POWER;	sops.sem_op = (short)-1;	sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);

		//-------- Plot Power Monitor
		for(IFindex=0; IFindex<Nif; IFindex++){
			memcpy(tempPower, bitPower[IFindex], timeNum*sizeof(float));
			bitPower[IFindex][0] = 10.0* log10(param_ptr->power[IFindex]);
			memcpy(&bitPower[IFindex][1], tempPower, (timeNum-1)*sizeof(float));
		}
		cpg_power(param_ptr, bitPower);
	}
	cpgend();
//------------------------------------------ RELEASE the SHM
    return(0);
}
