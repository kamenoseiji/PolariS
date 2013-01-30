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

    //-------- SHARED Segment Data --------
#ifdef HIDOI
	if(shm_init_create(
		SEGDATA_KEY,					// ACCESS KEY
		SEGDATA_SIZE,					// Data Area Size
		&(param_ptr->shrd_seg_id),		// SHM ID
		&segdata_ptr) == -1){			// Pointer to the shared data
		perror("Can't Create Shared Segment data!!\n"); return(-1);
	}
	printf("Allocated %d bytes for Shared Segment data [%d]!\n", SEGDATA_SIZE, param_ptr->shrd_seg_id);
#endif

    //-------- SHARED cross-power-spectra data area --------
	if(shm_init_create(
		XSPEC_KEY,						// ACCESS KEY
		XSPEC_SIZE,						// Data Area Size
		&(param_ptr->shrd_xspec_id),	// SHM ID
		&xspec_ptr) == -1){				// Pointer to the shared segment data
		perror("Can't Create Shared XSPEC data!!\n"); return(-1);
	}
	printf("Allocated %d bytes for Shared Xspec data [%d]!\n", XSPEC_SIZE, param_ptr->shrd_xspec_id);
//------------------------------------------ WAIT UNTIL FINISH
    while((param_ptr->validity & ABSFIN) == 0 ){

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
