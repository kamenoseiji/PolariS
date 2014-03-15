//	shm_alloc.c : Open Shared Memory Area
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
	struct	sembuf	sops;				// Semaphore for data area
	unsigned char	*k5head_ptr;		// Pointer to the Shared Data 
	unsigned char	*k5data_ptr;		// Pointer to the Shared Data 
	float	*segdata_ptr;				// Pointer to the shared segment data
	float	*xspec_ptr;					// Pointer to the cross-power spectra
	int		index;
//------------------------------------------ Access to the Shared Param
	if(shm_access(SHM_PARAM_KEY, sizeof(struct SHM_PARAM), &shrd_param_id, &param_ptr) == -1){
		perror("  Error : shm_alloc() can't access to the shared memory!!");   return(-1);
	};
//------------------------------------------ ALLOC SHARED MEMORY
	//-------- Semaphore for data area
	param_ptr->sem_data_id = semget(SEM_DATA_KEY, SEM_NUM, IPC_CREAT | 0666);
	for(index=0; index<SEM_NUM; index++){
		sops.sem_num = (ushort)index;
		sops.sem_op = (short)0;
		sops.sem_flg = IPC_NOWAIT;
		semop( param_ptr->sem_data_id, &sops, 1);
	}

    //-------- SHARED K5-HEADER AREA --------
	if(shm_init_create(
		K5HEAD_KEY,						// ACCESS KEY
		K5HEAD_SIZE,					// Data Area Size
		&(param_ptr->shrd_k5head_id),	// SHM ID
		&k5head_ptr) == -1){			// Pointer to the shared data
		perror("Can't Create Shared K5 header!!\n"); return(-1);
	}
	memset(k5head_ptr, 0x00, K5HEAD_SIZE);
	printf("Allocated %d bytes for Shared K5 header [%d]!\n", K5HEAD_SIZE, param_ptr->shrd_k5head_id);

    //-------- SHARED K5-DATA AREA --------
	if(shm_init_create(
		K5DATA_KEY,						// ACCESS KEY
		MAX_SAMPLE_BUF,					// Data Area Size
		&(param_ptr->shrd_k5data_id),	// SHM ID
		&k5data_ptr) == -1){			// Pointer to the shared data
		perror("Can't Create Shared K5 data!!\n"); return(-1);
	}
	memset(k5data_ptr, 0x00, MAX_SAMPLE_BUF);
	printf("Allocated %d bytes for Shared K5 data [%d]!\n", MAX_SAMPLE_BUF, param_ptr->shrd_k5data_id);

    //-------- SHARED cross-power-spectra data area --------
	if(shm_init_create(
		XSPEC_KEY,						// ACCESS KEY
		XSPEC_SIZE,						// Data Area Size
		&(param_ptr->shrd_xspec_id),	// SHM ID
		&xspec_ptr) == -1){				// Pointer to the shared segment data
		perror("Can't Create Shared XSPEC data!!\n"); return(-1);
	}
	printf("Allocated %d bytes for Shared Xspec data [%d]!\n", XSPEC_SIZE, param_ptr->shrd_xspec_id);
//------------------------------------------ End of Process
    exit(0);
}
