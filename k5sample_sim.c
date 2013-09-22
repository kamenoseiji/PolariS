//	k5sample_sim.c : Simulate K5-sampler
//
//	Author : Seiji Kameno
//	Created: 2012/11/8
//
#include "shm_k5data.inc"
#include "k5dict.inc"
#include <stdint.h>


main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	int		shrd_k5data_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore
	unsigned char	*k5head_ptr;		// Pointer to the shared K5 header
	unsigned char	*k5data_ptr;		// Pointer to the shared K5 data
	unsigned char	*shm_write_ptr;		// Writing Pointer

	int		year;
	int		index;
	int		num_read_cycle, read_fraction;	// K5/VSSP32 Reading cycle and fraction
	int		fsample, qbit, num_IF, filter;
	char	buf[K5FIFO_SIZE];
//------------------------------------------ Access to the SHARED MEMORY
    //-------- SHARED PARAMETERS --------
    if(shm_access(
        SHM_PARAM_KEY,					// ACCESS KEY
        sizeof(struct SHM_PARAM),		// SIZE OF SHM
        &shrd_param_id,					// SHM ID
        &param_ptr) == -1){				// Pointer to the SHM
		perror("  Error : Can't access to the shared memory!!");
		return(-1);
	}
	printf("Succeeded to access the shared parameter [%d]!\n",  param_ptr->shrd_param_id);

    //-------- SHARED K5 Header and Data --------
	k5head_ptr = shmat( param_ptr->shrd_k5head_id, NULL, 0 );
	k5data_ptr = shmat( param_ptr->shrd_k5data_id, NULL, 0 );
	param_ptr->validity |= ENABLE;		// Set Shared memory readiness bit to 1
//------------------------------------------ Start Sampling
	param_ptr->sd_len = 32e6;							// Size of 1-sec sampling data [bytes]
	param_ptr->fsample= K5HEAD_FS[fsample]*1000;		// Sampling frequency [Hz]
	param_ptr->num_st = K5HEAD_CH[num_IF];				// Number of IFs
	param_ptr->qbit	  = K5HEAD_QB[qbit];				// Quantization Bits
	param_ptr->segLen = SegLen;							// Segment length 
	param_ptr->segNum = NsegSec;						// Number of segments in 1 sec
	param_ptr->validity |= ACTIVE;		// Set Sampling Activity Bit to 1

 	while( param_ptr->validity & ACTIVE ){
		if( param_ptr->validity & (FINISH + ABSFIN) ){	break; }

		//-------- Read Data
		shm_write_ptr = k5data_ptr;		// Initialize pointer 
		//---- First half
		memset(k5data_ptr, 0x66, param_ptr->sd_len);
		usleep(500000);
	    //-------- Semaphore for Fiest Half--------
		for(index=0; index<4; index++){
			sops.sem_num = (ushort)index; sops.sem_op = (short)1; sops.sem_flg = (short)0;
			semop(param_ptr->sem_data_id, &sops, 1);
		}

		//---- Last half
		usleep(500000);
	    //-------- Semaphore for Fiest Half--------
		for(index=4; index<8; index++){
			sops.sem_num = (ushort)index; sops.sem_op = (short)1; sops.sem_flg = (short)0;
			semop(param_ptr->sem_data_id, &sops, 1);
		}
	}
//------------------------------------------ Stop Sampling
//	rv = ioctl(fd_in, TDSIO_SAMPLING_STOP);	close(fd_in);
	param_ptr->validity &= (~ACTIVE);		// Set Sampling Activity Bit to 0

    return(0);
}
