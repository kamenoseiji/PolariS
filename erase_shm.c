//	erase_shm.c : Erase Shared Memory Area
//
//	Author : Seiji Kameno
//	Created: 2012/10/18
//
#include "shm_k5data.inc"

int erase_shm(
	struct SHM_PARAM	*param_ptr )	// Pointer to the Shared Param
{
//-------- RELEASE SHARED MEMORY
	int		index;						// General counter

	shmctl(param_ptr->shrd_k5head_id, IPC_RMID, 0);		// Release K5 header
	shmctl(param_ptr->shrd_k5data_id, IPC_RMID, 0);		// Release K5 data
//	shmctl(param_ptr->shrd_seg_id, IPC_RMID, 0);		// Release Segment data
	shmctl(param_ptr->shrd_xspec_id, IPC_RMID, 0);		// Release xspec data

	sleep(1);
	semctl(param_ptr->sem_data_id, 0, IPC_RMID, 0);		// Release Semaphore
	sleep(1);
	shmctl(param_ptr->shrd_param_id, IPC_RMID, 0);		// Release shared param 

    return(0);
}
