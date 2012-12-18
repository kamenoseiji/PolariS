#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <error.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sem.h>
#include <sys/msg.h>

//-------- Shared Memory
#define	SEM_DATA_KEY	0x2000
#define	SHM_PARAM_KEY	0x1000
#define	K5HEAD_KEY 		0x0800
#define	K5DATA_KEY 		0x0400
#define	SEGDATA_KEY		0x0200
#define	XSPEC_KEY		0x0100

//-------- Shared Memory Size
#define SEM_NUM			16			// Number of Semaphores
#define	K5FIFO_SIZE		65536		// 64 kbytes
#define	K5HEAD_SIZE		32			// K5 Header size
#define	MAX_SAMPLE_BUF	33554432	// 256 Mbit = 32 MB
#define	MAX_seg_len		262144		// Segment size (2^18)
#define	NFFT			262144		// FFT points
#define	NFFTC			131073		// NFFT/2+1 = Number of spectral points
#define	MAX_seg_sec		128			// # of segments per sec (128)
#define	Nseg_halfsec	64			// # of segments per 0.5 sec (64)
#define	MAX_IF			4			// Number of IF
#define	MAX_CROSS		8			// <XX*>, <XY*>, <YX*>, <YY*>
#define	SEGDATA_SIZE	MAX_IF* MAX_seg_len* MAX_seg_sec* sizeof(float)
#define	XSPEC_SIZE		MAX_CROSS* NFFTC* sizeof(float)

//-------- Memory Map
#define	SEG_SAMPLE_OFFSET	MAX_SAMPLE_BUF

//-------- Shared Parameter
struct	SHM_PARAM{
	int		shrd_param_id;		// Shared Memory ID
	int		shrd_k5head_id;		// Shared Memory ID
	int		shrd_k5data_id;		// Shared Memory ID
	int		shrd_seg_id;		// Shared Memory ID
	int		shrd_xspec_id;		// Shared Memory ID
	int		sem_data_id;		// Semaphore
	int		dummy1;				// filler, reserved
	int		dummy2;				// filler, reserved

	//-------- Sampling Modes
	int		sd_len;				// Size of 1-sec sampling data [bytes]
	int		num_st;				// Total number of streams (IFs)
	int		qbit;				// Quantization bits (1, 2, 4, or 8)
	int		fsample;			// Sampling Frequency [Hz]
	int		seg_len;			// 1-segment length
	int		seg_num;			// Number of segments in 1 sec

	//-------- Process ID
	int		pid_shm_alloc;		// Process ID
	int		pid_k5sample;		// Process ID
	int		pid_segment;		// Process ID
	int		pid_fft;			// Process ID

	//-------- Status
	int		validity;			// Running / Ending signal
};

//-------- Dictionary of validity
#define	UNDEF		0			// 0000 0000 0000 0000 0000 0000 0000 0000
#define	DISABLE		1			// 0000 0000 0000 0000 0000 0000 0000 0001
#define	ENABLE		2			// 0000 0000 0000 0000 0000 0000 0000 0010
#define	ACTIVE		4			// 0000 0000 0000 0000 0000 0000 0000 0100
#define	RESERVED	8			// 0000 0000 0000 0000 0000 0000 0000 1000
#define	DONE		16			// 0000 0000 0000 0000 0000 0000 0001 0000
#define	CONFIRMED	32			// 0000 0000 0010 0000 0000 0000 0000 0000
#define	FINISH		0x40000000	// 0100 0000 0000 0000 0000 0000 0000 0000
#define	ABSFIN		0x80000000	// 1000 0000 0000 0000 0000 0000 0000 0000

//-------- Semaphore Map
#define	SEM_IF0_SMP_F	0			// A/D for IF0, first half of every sec. 
#define	SEM_IF1_SMP_F	1			// A/D for IF1, first half of every sec. 
#define	SEM_IF2_SMP_F	2			// A/D for IF2, first half of every sec. 
#define	SEM_IF3_SMP_F	3			// A/D for IF3, first half of every sec. 
#define	SEM_IF0_SMP_L	4			// A/D for IF0, last half of every sec. 
#define	SEM_IF1_SMP_L	5			// A/D for IF1, last half of every sec. 
#define	SEM_IF2_SMP_L	6			// A/D for IF2, last half of every sec. 
#define	SEM_IF3_SMP_L	7			// A/D for IF3, last half of every sec. 
#define	SEM_SEG_F		8			// S-part bit->weight for first half
#define	SEM_SEG_L		9			// S-part bit->weight for first half
#define	SEM_FX			10			// FFT and Cross Corr. for every sec.

//-------- Chiled Processes
#define SHM_ALLOC	"/usr/local/custom/bin/shm_alloc"
#define K5_SAMPLE	"/usr/local/custom/bin/k5sample_store"
#define SHM_SEGDATA	"/usr/local/custom/bin/shm_segdata"
#define CUDA_FFT	"/usr/local/custom/bin/cuda_fft_xspec"
