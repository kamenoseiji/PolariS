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
#define	HALFSEC_OFFSET	 8388608	// 8192k words per 1/2 sec
#define	SegLen			262144		// Segment size (2^18)
#define	NFFT			262144		// FFT points
#define	NFFT2			131072		// Half of FFT points
#define	NFFTC			131073		// NFFT/2+1 = Number of spectral points
#define	NsegSec			128			// # of segments per sec (128)
#define	NsegSec2		64			// # of segments per 1/2 sec (64)
#define	Nif				4			// Number of IF
#define	Ncross			4			// Number of correlations <XX*>, <XY*>, <YX*>, <YY*>
#define	SEGDATA_SIZE	Nif* SegLen* NsegSec* sizeof(float)
#define	XSPEC_SIZE		Nif* Ncross* NFFT2* sizeof(float)
#define	MAX_FILE_REC	1800		// Maximum number of records in output file

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
	int		segLen;				// 1-segment length
	int		segNum;				// Number of segments in 1 sec
	int		num_ch;				// Number of spectral channels
	int		dummy3;

	//-------- Process ID
	int		pid_shm_alloc;		// Process ID
	int		pid_k5sample;		// Process ID
	int		pid_segment;		// Process ID
	int		pid_fft;			// Process ID
	int		dummy4;
	int		dummy5;
	int		dummy6;
	int		dummy7;

	//-------- Status
	int		validity;			// Running / Ending signal
	int		year;				// Year
	int		doy;				// Day of Year
	int		hour;				// Hour (UT)
	int		min;				// Minute
	int		sec;				// Second
	int		current_rec;		// current record index
	int		integ_rec;			// Total integration numer (MAX = 1800) 
};

//-------- Dictionary of validity bits
#define	UNDEF		0x00000000	// 0000 0000 0000 0000 0000 0000 0000 0000
#define	DISABLE		0x00000001	// 0000 0000 0000 0000 0000 0000 0000 0001
#define	ENABLE		0x00000002	// 0000 0000 0000 0000 0000 0000 0000 0010
#define	ACTIVE		0x00000004	// 0000 0000 0000 0000 0000 0000 0000 0100
#define	RESERVED	0x00000008	// 0000 0000 0000 0000 0000 0000 0000 1000
#define	DONE		0x00000010	// 0000 0000 0000 0000 0000 0000 0001 0000
#define	CONFIRMED	0x00000020	// 0000 0000 0000 0000 0000 0000 0010 0000
#define A00_REC		0x00000100	// 0000 0000 0000 0000 0000 0001 0000 0000
#define A01_REC		0x00000200	// 0000 0000 0000 0000 0000 0010 0000 0000
#define A02_REC		0x00000400	// 0000 0000 0000 0000 0000 0100 0000 0000
#define A03_REC		0x00000800	// 0000 0000 0000 0000 0000 1000 0000 0000
#define P00_REC		0x00001000	// 0000 0000 0000 0000 0001 0000 0000 0000
#define P01_REC		0x00002000	// 0000 0000 0000 0000 0010 0000 0000 0000
#define P02_REC		0x00004000	// 0000 0000 0000 0000 0100 0000 0000 0000
#define P03_REC		0x00008000	// 0000 0000 0000 0000 1000 0000 0000 0000
#define C00_REC		0x00010000	// 0000 0000 0000 0001 0000 0000 0000 0000
#define C01_REC		0x00020000	// 0000 0000 0000 0010 0000 0000 0000 0000
#define PGPLOT		0x00040000	// 0000 0000 0000 0100 0000 0000 0000 0000
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
#define K5_SIM		"/usr/local/custom/bin/k5sim"
#define SHM_SEGDATA	"/usr/local/custom/bin/shm_segdata"
#define CUDA_FFT	"/usr/local/custom/bin/cuda_fft_xspec"
#define SPEC_VIEW	"/usr/local/custom/bin/shm_spec_view"
