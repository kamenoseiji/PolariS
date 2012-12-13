//	polaris_start.c : Start Polaris-Related Processes
//
//	Author : Seiji Kameno
//	Created: 2012/10/18
//
#include "shm_k5data.inc"
#include <unistd.h>
#include <sys/wait.h>

main(
	int		argc,			// Number of Arguments
	char	**argv )		// Pointer to Arguments
{
	int		shrd_param_id;				// Shared Memory ID
	struct	SHM_PARAM	*param_ptr;		// Pointer to the Shared Param
	char	cmd[8][16];					// Command line arguments
	int		pid;						// Process ID
	int		index;
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
//------------------------------------------ Start K5 sampling
	if( fork() == 0){
		pid = getpid(); sprintf(cmd[0], "k5sample_store");
		printf(" Exec %s as Chiled Process [PID = %d]\n", cmd[0], pid);
		if( execl( K5_SAMPLE, cmd[0], (char *)NULL ) == -1){
			perror("Can't Create Chiled Proces!!\n"); return(-1);
		}
	}


    return(0);
}
