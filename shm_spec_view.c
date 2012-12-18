//	shm_spec_view.c : Cross Power Spectrom Viewer
//
//	Author : Seiji Kameno
//	Created: 2012/12/18
//
#include "shm_k5data.inc"
#include "k5dict.inc"

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
//------------------------------------------ Browse
#ifdef HIDOI
	printf("-------- Param IDs --------\n");
	printf(" PARAM   K5DATA  SEGDATA    XSPEC\n");
	printf("%06d   %06d   %06d   %06d\n \n", param_ptr->shrd_param_id, param_ptr->shrd_k5data_id, param_ptr->shrd_segdata_id, param_ptr->shrd_xspec_id);

#endif
//------------------------------------------ K5 Header and Data

	setvbuf(stdout, (char *)NULL, _IONBF, 0);	// Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- Wait for Semaphore
		printf("Wait for Semaphore\n");
		sops.sem_num = (ushort)SEM_FX;	sops.sem_op = (short)-1;	sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);

		printf("%7.2e %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e %7.2e\n",
			xspec_ptr[0], xspec_ptr[1], xspec_ptr[2], xspec_ptr[3],
			xspec_ptr[4], xspec_ptr[5], xspec_ptr[6], xspec_ptr[7]);

	}
//------------------------------------------ RELEASE the SHM
    return(0);
}
