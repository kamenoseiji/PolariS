//	shm_init_create.c : Create a Shared Memory Area
//
//	Author : Seiji Kameno
//	Created: 2012/10/18
//
#include	<stdio.h>
#include	<unistd.h>
#include	<sys/types.h>		// InterProcess Communication 
#include	<sys/ipc.h>			// InterProcess Communication
#include	<sys/shm.h>			// Shared Memory
#include	<string.h>			// Shared Memory


int shm_init_create(
	key_t	shm_key,		// INPUT:   Keyword for Shared Memory
	size_t	shm_size,		// INPUT:   Size of the SHM [bytes]
	int		*shrd_id,		// OUTPUT:  Pointer to the Shared Memory ID 
	int		*shm_ptr)		// OUTPUT:  Pointer to the Top of SHM
{
/*
-------------------------------------------------- OPEN SHARED MEMORY
*/
	//-------- ACCESS TO CURRENT STATUS --------
	*shrd_id = shmget(shm_key, shm_size, IPC_CREAT|0666);
    if( *shrd_id < 0 ){
		printf("[shm_init_create]: Can't Access to the Shared Memory !! \n" );
		return(-1);
	}

	*shm_ptr = (int)shmat( *shrd_id, NULL, 0);
	memset( (void *)*shm_ptr, 0, shm_size );

	/*-------- ENDING --------*/
	return( *shrd_id);
}
