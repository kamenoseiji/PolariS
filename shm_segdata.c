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
	short 	*k5init_ptr, *k5data_ptr;	// Pointer to the shared K5 data
	float	*seginit_ptr, *segdata_ptr;	// Pointer to write the segment data
	int		index, index_seg;
	int		index_if, index_smp;
	int		overlap[128];				// Number of samples to overlap
	int		offset[128];				// Number of samples for segment offset
	int		mean_offset, fraction;
	short	bitmask, bitshift;			// Bit shift and mask
//-------- Argument Check
	if(argc < ARGNUM){
		perror(" IF number must be specified!!\n");
		return;
	}
	 setvbuf(stdout, (char *)NULL, _IONBF, 0);   // Disable stdout cache
//------------------------------------------ Access to the SHARED MEMORY
	shrd_param_id = shmget( SHM_PARAM_KEY, sizeof(struct SHM_PARAM), 0444);
	param_ptr = (struct SHM_PARAM *)shmat(shrd_param_id, NULL, 0);
	k5head_ptr = (unsigned char *)shmat(param_ptr->shrd_k5head_id, NULL, SHM_RDONLY);
	k5init_ptr = (short *)shmat(param_ptr->shrd_k5data_id, NULL, SHM_RDONLY);
	seginit_ptr = (float *)shmat(param_ptr->shrd_seg_id, NULL, 0);

	//-------- IF-related parameters
	index_if = atoi(argv[IF_index]);	// IF# (0, 1, 2, or 3)
	bitmask = k5bitmask[index_if];		// to pick 4-bit from K5 data
	bitshift= k5bitshift[index_if];		// bitshift for 4-bit data

	mean_offset = (param_ptr->fsample - param_ptr->seg_len) / (param_ptr->seg_num - 1);
	fraction = (param_ptr->fsample - param_ptr->seg_len) % mean_offset;
	offset[0] = 0;	overlap[0] = 0;
	printf("Mean offset = %d, fraction = %d\n", mean_offset, fraction);
	for( index_seg = 1; index_seg < param_ptr->seg_num; index_seg ++){
		offset[index_seg] = mean_offset;
		if( index_seg % (param_ptr->seg_num / fraction) == 1 ){	offset[index_seg] ++;}
		overlap[index_seg] = param_ptr->seg_len - offset[index_seg];

//		printf("segment %d : offset = %d  overlap = %d\n", index_seg, offset[index_seg], overlap[index_seg]);
	}

	seginit_ptr += (index_if* param_ptr->seg_num* param_ptr->seg_num);	// Write Pointer Offset
//------------------------------------------ K5 Header and Data
	while(param_ptr->validity & ACTIVE){
		if( param_ptr->validity & (FINISH + ABSFIN) ){  break; }

		segdata_ptr = seginit_ptr; k5data_ptr  = k5init_ptr;	// Return to the initial pointer

		// Wait for Semaphore
		sops.sem_num = (ushort)index_if; sops.sem_op = (short)-1; sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);

		//-------- First Segment --------
		for(index_smp=0; index_smp<param_ptr->seg_len; index_smp ++){
			*segdata_ptr = wt4[(*k5data_ptr & bitmask) >> bitshift]; segdata_ptr ++; k5data_ptr ++;
		}

		//-------- First Half
		for( index_seg=1; index_seg<param_ptr->seg_num/2; index_seg ++){
			memcpy(segdata_ptr, segdata_ptr - overlap[index_seg], overlap[index_seg]* sizeof(float));
//			printf("Seg %d : Copied duplicated overlap : addr=%X\n", index_seg, segdata_ptr);
			segdata_ptr += overlap[index_seg];
			for(index_smp=overlap[index_seg]; index_smp<param_ptr->seg_len; index_smp ++){
				*segdata_ptr = wt4[(*k5data_ptr & bitmask) >> bitshift];
				segdata_ptr ++; k5data_ptr ++;
			}
		}
			
		printf(" SEM for IF %d --- first half copied!!\n", index_if);

		// Wait for Semaphore
		sops.sem_num = (ushort)(4+index_if); sops.sem_op = (short)-1; sops.sem_flg = (short)0;
		semop( param_ptr->sem_data_id, &sops, 1);

		//-------- Last Half
		for( index_seg=param_ptr->seg_num/2; index_seg< param_ptr->seg_num; index_seg ++){
			memcpy(segdata_ptr, segdata_ptr - overlap[index_seg], overlap[index_seg]* sizeof(float));
//			printf("Seg %d : Copied duplicated overlap : addr=%X\n", index_seg, segdata_ptr);
			segdata_ptr += overlap[index_seg];
			for(index_smp=overlap[index_seg]; index_smp<param_ptr->seg_len; index_smp ++){
				*segdata_ptr = wt4[(*k5data_ptr & bitmask) >> bitshift];
				segdata_ptr ++; k5data_ptr ++;
			}
		}
		printf(" SEM for IF %d --- latter half copied!!\n", index_if);
	}
//------------------------------------------ RELEASE the SHM
    return(0);
}
