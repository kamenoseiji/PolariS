//	k5sample_store.c : Store K5-sampled data into shared memory
//
//	Author : Seiji Kameno
//	Created: 2012/11/8
//
#include "shm_k5data.inc"
// #include "k5dict.inc"
#include <stdint.h>
#include <sys/ioctl.h>
#include <sys/tdsio.h>
#define DEV_FILE "/dev/utds0"	// Device driver of VSSP32

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	int		shrd_k5data_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	struct	sembuf		sops;			// Semaphore
	char	dev[] = DEV_FILE;			// Device for K5/VSSP
	int		fd_in;						// ID of the K5/VSSP32 device
	int		rv;							// Return Value from K5/VSSP32
	unsigned char	*k5head_ptr;		// Pointer to the shared K5 header
	unsigned char	*k5data_ptr;		// Pointer to the shared K5 data
	unsigned char	*shm_write_ptr;		// Writing Pointer
//	FILE	*dumpfile_ptr;				// Dump File
	int		K5HEAD_CH[] = {1, 4};		// 0:1ch, 1:4ch

	int		year;
	int		index;
	int		num_read_cycle, read_fraction;	// K5/VSSP32 Reading cycle and fraction
	int		fsample, qbit, num_IF, filter;
	char	buf[K5FIFO_SIZE];
//------------------------------------------ Open K5/VSSP32 device
	if( (fd_in = open(dev, O_RDONLY)) == -1){		// Open VSSP device
		perror("  Error : Can't open K5/VSSP32 device");
		return(-1);
	}
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
//	dumpfile_ptr = fopen("K5.dump", "w");
//------------------------------------------ Start Sampling
	rv = ioctl(fd_in, TDSIO_INIT);	printf("TDSIO_INIT result in %d\n", rv);
	rv = ioctl(fd_in, TDSIO_BUFFER_CLEAR);  printf("TDSIO_BUFFER_CLEAR result in %d\n", rv);
	rv = ioctl(fd_in, TDSIO_GET_YEAR, &year);

	//---- Set sampling frequency
	index = param_ptr->fsample / 1000000; fsample=TDSIO_SAMPLING_1MHZ;
	while(index > 1){	fsample++; index = index >> 1;}
	rv = ioctl(fd_in, TDSIO_SET_FSAMPLE, &fsample);

	//---- Set Quantization Bits
	index = param_ptr->qbit;	qbit = TDSIO_SAMPLING_1BIT;
	while(index > 1){	qbit++; index = index >> 1;}
	rv = ioctl(fd_in, TDSIO_SET_RESOLUTIONBIT, &qbit);	

	//---- Set Antialias Filter
	if( (param_ptr->filter <= 0) || (param_ptr->filter > 16)){
		filter = TDSIO_SAMPLING_THRU;
	} else {
		index = 16 / param_ptr->filter;	filter = TDSIO_SAMPLING_16M;
		while(index > 1){	filter++; index = index >> 1;}
	}
	rv = ioctl(fd_in, TDSIO_SET_FILTER, &filter);

	//---- Set Number of IFs
	num_IF  = TDSIO_SAMPLING_4CH;   rv = ioctl(fd_in, TDSIO_SET_CHMODE, &num_IF);

	//---- Check Sampling Modes
	rv = ioctl(fd_in,TDSIO_GET_FSAMPLE, &fsample);  	// printf("Fsample=%d\n", fsample);
	rv = ioctl(fd_in,TDSIO_GET_RESOLUTIONBIT,  &qbit);	// printf("Qbit =%d\n", qbit);
	rv = ioctl(fd_in,TDSIO_GET_CHMODE,  &num_IF);		// printf("Num of IF =%d\n", num_IF);
	rv = ioctl(fd_in,TDSIO_GET_FILTER,  &filter);		// printf("Filter =%d\n", filter);
	rv = ioctl(fd_in, TDSIO_SAMPLING_START);			// printf("TDSIO_SAMPLING_START result in %d\n", rv);

	param_ptr->sd_len = MAX_SAMPLE_BUF;					// Size of 1-sec sampling data [bytes]
	param_ptr->num_st = K5HEAD_CH[num_IF];				// Number of IFs
	num_read_cycle = param_ptr->sd_len / K5FIFO_SIZE;	// Number of read cycles
	read_fraction  = param_ptr->sd_len % K5FIFO_SIZE;	// Fraction bytes
	param_ptr->validity |= ACTIVE;		// Set Sampling Activity Bit to 1

 	while( param_ptr->validity & ACTIVE ){
		if( param_ptr->validity & (FINISH + ABSFIN) ){	break; }
		//-------- Read Header
		rv = read(fd_in, k5head_ptr, K5HEAD_SIZE);
		// printf("READ result in %d : %X\n", rv, (int)k5head_ptr);

		//-------- Read Data
		shm_write_ptr = k5data_ptr;		// Initialize pointer 
		//---- First half
		for(index=0; index<num_read_cycle/2 + 2; index++){
			rv = read(fd_in, shm_write_ptr, K5FIFO_SIZE); shm_write_ptr += rv;
		}

	    //-------- Semaphore for Fiest Half--------
		for(index=0; index<Nif; index++){
			sops.sem_num = (ushort)index; sops.sem_op = (short)1; sops.sem_flg = (short)0;
			semop(param_ptr->sem_data_id, &sops, 1);
		}

		//---- Last half
		for(index=num_read_cycle/2 + 2; index<num_read_cycle; index++){
			rv = read(fd_in, shm_write_ptr, K5FIFO_SIZE); shm_write_ptr += rv;
		}
		rv = read(fd_in, shm_write_ptr, read_fraction); shm_write_ptr += rv;
	    //-------- Semaphore for Fiest Half--------
		for(index=Nif; index<2*Nif; index++){
			sops.sem_num = (ushort)index; sops.sem_op = (short)1; sops.sem_flg = (short)0;
			semop(param_ptr->sem_data_id, &sops, 1);
		}
//		fwrite(k5data_ptr, MAX_SAMPLE_BUF, 1, dumpfile_ptr);
	}
//------------------------------------------ Stop Sampling
//	fclose(dumpfile_ptr);
	rv = ioctl(fd_in, TDSIO_SAMPLING_STOP);	close(fd_in);
	param_ptr->validity &= (~ACTIVE);		// Set Sampling Activity Bit to 0

    return(0);
}
