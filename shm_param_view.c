//	shm_param_view.c : Shared Parameter Browser
//
//	Author : Seiji Kameno
//	Created: 2012/11/8
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
	unsigned char		*k5head_ptr;	// Pointer to the shared K5 header
	unsigned char		*k5data_ptr;	// Pointer to the shared K5 data
	int		index;
	int		sod=0, hour=0, min=0, sec=0;// Second of day, hour, min, sec
//------------------------------------------ Access to the SHARED MEMORY
    //------- SHARED PARAMETERS --------
    if(shm_access(
        SHM_PARAM_KEY,					// ACCESS KEY
        sizeof(struct SHM_PARAM),		// SIZE OF SHM
        &shrd_param_id,					// SHM ID
        &param_ptr) != -1){				// Pointer to the SHM
		printf("Succeeded to access the shared parameter [%d]!\n",  param_ptr->shrd_param_id);
	}
//------------------------------------------ Browse
#ifdef HIDOI
	printf("-------- Param IDs --------\n");
	printf(" PARAM   K5DATA  SEGDATA    XSPEC\n");
	printf("%06d   %06d   %06d   %06d\n \n", param_ptr->shrd_param_id, param_ptr->shrd_k5data_id, param_ptr->shrd_segdata_id, param_ptr->shrd_xspec_id);

#endif
	printf("-------- Validitys --------\n");
	printf("validity = %X\n", param_ptr->validity);
//------------------------------------------ K5 Header and Data
	k5head_ptr = shmat(param_ptr->shrd_k5head_id, NULL, SHM_RDONLY);
//	for(index=0; index<32; index ++){ printf("%02X ", k5head_ptr[index]); }
	printf("\n");

	k5data_ptr = shmat(param_ptr->shrd_k5data_id, NULL, SHM_RDONLY);
//	for(index=0; index<32; index ++){ printf("%02X ", k5data_ptr[index]); }
	printf("\n");

	setvbuf(stdout, (char *)NULL, _IONBF, 0);	// Disable stdout cache
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- Wait for Semaphore
		sops.sem_num = (ushort)0;	sops.sem_op = (short)-1;	sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);

		memcpy(&sod, &k5head_ptr[4], 2);
		sod |= (k5head_ptr[6] & 0x01) << 17; sod2hms(sod, &hour, &min, &sec);
//		printf("UT=%02d:%02d:%02d : %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\r",
		printf("UT=%02d:%02d:%02d : %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n",
			hour, min, sec,
			wt4[ k5data_ptr[0] & 0x0f],
			wt4[(k5data_ptr[0] & 0xf0) >> 4],
			wt4[ k5data_ptr[1] & 0x0f],
			wt4[(k5data_ptr[1] & 0xf0) >> 4],
			wt4[ k5data_ptr[2] & 0x0f],
			wt4[(k5data_ptr[2] & 0xf0) >> 4],
			wt4[ k5data_ptr[3] & 0x0f],
			wt4[(k5data_ptr[3] & 0xf0) >> 4],
			wt4[ k5data_ptr[4] & 0x0f],
			wt4[(k5data_ptr[4] & 0xf0) >> 4],
			wt4[ k5data_ptr[5] & 0x0f],
			wt4[(k5data_ptr[5] & 0xf0) >> 4],
			wt4[ k5data_ptr[6] & 0x0f],
			wt4[(k5data_ptr[6] & 0xf0) >> 4],
			wt4[ k5data_ptr[7] & 0x0f],
			wt4[(k5data_ptr[7] & 0xf0) >> 4]);

		//-------- Wait for Semaphore
		sops.sem_num = (ushort)1;	sops.sem_op = (short)-1;	sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);
		printf("HIDOI\n");

	}
//------------------------------------------ RELEASE the SHM
    return(0);
}
