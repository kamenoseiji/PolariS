//	shm_segdata.c : Extend segment data from K5-sampled data
//
//	Author : Seiji Kameno
//	Created: 2012/12/12
//
#include "shm_k5data.inc"
#include "k5dict.inc"
#include <stdio.h>
#include <unistd.h>
#include <sys/ipc.h>
#include <sys/types.h>
#include <sys/shm.h>
#define	NFFT		262144
#define	NSEG		128
#define	IF_index	1		// Command line arguments
#define	ARGNUM		2		// Number of arguments

// int	shm_access(key_t shm_key, size_t shm_size, int *shrd_id, int *shm_ptr);

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore for data access
	unsigned char		*k5head_ptr;	// Pointer to the shared K5 header
	short 	*k5data_ptr;				// Pointer to the shared K5 data
	float	*segdata_ptr;
	int		index, index_seg, offset;
	int		if_index;
	short	bitmask, bitshift;			// Bit shift and mask
//-------- Argument Check
	if(argc < ARGNUM){
		perror(" IF number must be specified!!\n");
		return;
	}
//------------------------------------------ Access to the SHARED MEMORY
	shrd_param_id = shmget( SHM_PARAM_KEY, sizeof(struct SHM_PARAM), 0444);
	param_ptr = (struct SHM_PARAM *)shmat(shrd_param_id, NULL, 0);
	k5head_ptr = (unsigned char *)shmat(param_ptr->shrd_k5head_id, NULL, SHM_RDONLY);
	k5data_ptr = (short *)shmat(param_ptr->shrd_k5data_id, NULL, SHM_RDONLY);
	segdata_ptr = (float *)shmat(param_ptr->shrd_seg_id, NULL, 0);

	//-------- IF-related parameters
	if_index = atoi(argv[IF_index]);	// IF# (0, 1, 2, or 3)
	bitmask = k5bitmask[if_index];		// to pick 4-bit from K5 data
	bitshift= k5bitshift[if_index];		// bitshift for 4-bit data
	segdata_ptr += if_index* MAX_seg_len* MAX_seg_sec;	// Write Pointer Offset
//------------------------------------------ K5 Header and Data
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		//-------- First Segment --------
		for(smp_index=0; smp_index<MAX_seg_len; smp_index ++){
			*segdata_ptr = wt4[(*k5data_ptr & bitmask) >> bitshift];
			segdata_ptr ++;
			k5data_ptr ++;
		}
		//--------
		for( seg_index=1; seg_index<MAX_seg_sec/2; seg_index ++){
			memcpy( segdata_ptr, (segdata_ptr - MAX_seg_len + offset), (MAX_seg_len - offset)* sizeof(float));
			segdata_ptr,
		}
			

#ifdef HIDOI
		for(index=0; index<num_sample; index ++){
			*segdata_ptr = wt4[(k5data_ptr[index] & bitmask) >> bitshift ];
			segdata_ptr++;
		}
		segdata_ptr -= num_sample;
#endif
		// Wait for Semaphore
		sops.sem_num = (ushort)if_index; sops.sem_op = (short)-1; sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);
		printf(" SEM for IF %d !!\n", if_index);
	}
//------------------------------------------ RELEASE the SHM
    return(0);
}
