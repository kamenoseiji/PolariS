//	shm_init.c : Initialize Shared Memory Area and finish all processes
//
//	Author : Seiji Kameno
//	Created: 2012/10/18
//
#include "shm_k5data.inc"

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
/*
-------------------------------------------- ALLOC SHARED MEMORY
*/
    /*-------- SHARED PARAMETERS --------*/
    shm_access(
        SHM_PARAM_KEY,					// ACCESS KEY
        sizeof(struct SHM_PARAM),		// SIZE OF SHM
        &shrd_param_id,					// SHM ID
        &param_ptr);					// Pointer to the SHM
/*
-------------------------------------------- WAIT UNTIL FINISH
*/
    /*-------- Immediate FINISH --------*/
	if( argc >= 2){
		param_ptr->validity ^= atoi(argv[1]);	// Toggle Flag

	} else {
		param_ptr->validity |= FINISH;
	}
	printf("SET VALIDITY %X\n", param_ptr->validity);
/*
-------------------------------------------- RELEASE the SHM
*/
    return(0);
}
