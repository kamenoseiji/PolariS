//	polaris_start.c : Start Polaris-Related Processes
//
//	Author : Seiji Kameno
//	Created: 2012/10/18
//
#include "shm_k5data.inc"
#include <unistd.h>
#include <string.h>
#include <sys/wait.h>

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		ch_option;					// Option Charactors
	extern char	*optarg;				// Option Arguments
	extern int	optind, opterr;			// Option indicator and error
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	char	cmd[8][16];					// Command line arguments
	char	pgdev[16];					// PGPLOT Device
	int		pid;						// Process ID
	int		integPP=0;					// Maximum record number [sec]
	int		statusFlag = 0;				// Used for option parser
//------------------------------------------ Option Parser
	while(( ch_option = getopt(argc, argv, "a:c:hi:p:v:")) != -1){
		switch(ch_option){
			case 'a':	statusFlag |= (valid_bit(optarg) <<  8);	break;
			case 'p':	statusFlag |= (valid_bit(optarg) << 12);	break;
			case 'c':	statusFlag |= ((valid_bit(optarg) & 0x03) << 16);	break;
			case 'i':	integPP = atoi(optarg);	break;
			case 'v':	statusFlag |= PGPLOT; strcpy(pgdev, optarg);	break;
			case 'h':	usage();	return(0);
		}	
	}
//------------------------------------------ Start shm_alloc()
	if( fork() == 0){
		pid = getpid(); sprintf(cmd[0], "shm_alloc");
		printf(" Exec %s as Chiled Process [PID = %d]\n", cmd[0], pid);
		if( execl( SHM_ALLOC, cmd[0], (char *)NULL ) == -1){
			perror("Can't Create Chiled Proces!!\n"); return(-1);
		}
	}
	sleep(1);
//------------------------------------------ Access to the shared parameter
	if(shm_access(SHM_PARAM_KEY, sizeof(struct SHM_PARAM), &shrd_param_id, &param_ptr) == -1){
		 perror("  Error : Can't access to the shared memory!!");	return(-1);
	}
	param_ptr->pid_shm_alloc = pid;
//------------------------------------------ Set the validity bits
	param_ptr->validity |= statusFlag;
	param_ptr->integ_rec = integPP;
//------------------------------------------ Start K5 sampling
	if( fork() == 0){
		pid = getpid(); sprintf(cmd[0], "k5sample_store");
		printf(" Exec %s as Chiled Process [PID = %d]\n", cmd[0], pid);
		if( execl( K5_SAMPLE, cmd[0], (char *)NULL ) == -1){
			perror("Can't Create Chiled Proces!!\n"); return(-1);
		}
	}
//------------------------------------------ Start K5 simulator
#ifdef DEBUG
	if( fork() == 0){
		pid = getpid(); sprintf(cmd[0], "k5sim");
		printf(" Exec %s as Chiled Process [PID = %d]\n", cmd[0], pid);
		if( execl( K5_SIM, cmd[0], (char *)NULL ) == -1){
			perror("Can't Create Chiled Proces!!\n"); return(-1);
		}
	}
#endif
//------------------------------------------ Start Spectrum Viewer
	if( param_ptr->validity & PGPLOT ){
		strcpy(cmd[1], pgdev);
		if( fork() == 0){
			pid = getpid(); sprintf(cmd[0], "shm_spec_view");
			printf(" Exec %s as Chiled Process [PID = %d]\n", cmd[0], pid);
			if( execl( SPEC_VIEW, cmd[0], cmd[1], (char *)NULL ) == -1){
				perror("Can't Create Chiled Proces!!\n"); return(-1);
			}
		}
	}
//------------------------------------------ Start CUDA FFT
	sleep(1);		// Wait 1 sec
	if( fork() == 0){
		pid = getpid(); sprintf(cmd[0], "cuda_fft_xspec");
		printf(" Exec %s as Chiled Process [PID = %d]\n", cmd[0], pid);
		if( execl( CUDA_FFT, cmd[0], (char *)NULL ) == -1){
			perror("Can't Create Chiled Proces!!\n"); return(-1);
		}
	}
    return(0);
}

int usage(){
	fprintf(stderr, "USAGE: polaris_start [-chipv] \n");
	fprintf(stderr, "  -a : Specify autocorrelation files to save.\n");
	fprintf(stderr, "       0 -> CH0 is recorded, 12 -> CH1 and CH2 are recorded. Default: no autocorr, recorded. \n");
	fprintf(stderr, "  -c : Specify crosscorrelation files not saved.\n");
	fprintf(stderr, "       0 -> CH0xCH2 is recorded, 01 -> all Xcorrs are recorded. Default: no xcorr, recoreded. \n");
	fprintf(stderr, "  -h : Show help \n");
	fprintf(stderr, "  -i : Recording time [sec]. Unless specified, polaris keep recordeng until shm_init.\n"); 
	fprintf(stderr, "  -p : Specify bit-distibution files to save. Index is the same with -a option.\n");
	fprintf(stderr, "  -v : Specify PGPLOT window to display spectra. /xw -> X window, /gif -> GIF, /null -> no view.\n");
	return(0);
}

int valid_bit( char *option ){
	int		valid = 0;
	if(strstr( option, "0" ) != NULL)	valid |= 0x01;
	if(strstr( option, "1" ) != NULL)	valid |= 0x02;
	if(strstr( option, "2" ) != NULL)	valid |= 0x04;
	if(strstr( option, "3" ) != NULL)	valid |= 0x08;
	return(valid);
}
