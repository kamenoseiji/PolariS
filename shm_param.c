//	shm_param.c : Open Shared Parameter Area
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
//------------------------------------------ ALLOC SHARED MEMORY
    //-------- SHARED PARAMETERS --------
    if(shm_init_create(
        SHM_PARAM_KEY,					// ACCESS KEY
        sizeof(struct SHM_PARAM),		// SIZE OF SHM
        &shrd_param_id,					// SHM ID
        &param_ptr) == -1){						// Pointer to the SHM
		perror("Can't Create Shared Parameter!!\n"); return(-1);
	}
    param_ptr->shrd_param_id = shrd_param_id;   // Store PARAM ID
	printf("Allocated %d bytes for Shared Param [%d]!\n",  sizeof(struct SHM_PARAM), param_ptr->shrd_param_id);
//------------------------------------------ WAIT UNTIL FINISH
    while((param_ptr->validity & ABSFIN) == 0 ){

		if( (param_ptr->integ_rec > 0) && (param_ptr->current_rec >= param_ptr->integ_rec - 1)){
			param_ptr->validity |= FINISH;
		}

        //---- SOFT FINISH (Release SHM after 5 seconds) ----
        if( (param_ptr->validity & FINISH) != 0){	// detect FINISH
            sleep(5);
            break;
        }
        sleep(1);
    }
//------------------------------------------ RELEASE the SHM
    erase_shm(param_ptr);
    return(0);
}
