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
#define	MAX_SAMPLE_BUF	32000000	// 256 Mbit = 32 MB
#define	HALFBUF			16000000	// 16Mbytes per 1/2 sec
#define PARTNUM			2			// Number of Parts in 1 sec
#define	NFFT		(param_ptr->segLen)		// FFT points
#define	NFFT2		(param_ptr->num_ch)		// Half of FFT points
#define	NFFTC		(param_ptr->num_ch + 1)	// NFFT/2+1 = Number of spectral points
#define	NsegSec		(param_ptr->segNum)		// # of segments per sec
#define	NsegPart	(param_ptr->segNum/ PARTNUM)	// # of segments per part
#define	Nif				4			// Number of IF
#define	Ncross			4			// Number of correlations <XX*>, <XY*>, <YX*>, <YY*>
#define	SEGDATA_SIZE	Nif* param_ptr->segLen* param_ptr->segNum* sizeof(float)
#define	XSPEC_SIZE		Nif* Ncross* param_ptr->num_ch* sizeof(float)
#define	MAX_FILE_REC	1800		// Maximum number of records in output file

//-------- Memory Map

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
	int		filter;				// Antialias Filter BW (MHz)
	int		fsample;			// Sampling Frequency [Hz]
	int		segLen;				// 1-segment length
	int		segNum;				// Number of segments in 1 sec
	int		num_ch;				// Number of spectral channels

	//-------- Process ID
	int		pid_shm_alloc;		// Process ID
	int		pid_k5sample;		// Process ID
	int		pid_segment;		// Process ID
	int		pid_fft;			// Process ID

	//-------- Total Power
	float	power[4];				// Power

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
#define	SEM_POWER		11			// Total Power Measurement

//-------- Chiled Processes
#define SHM_PARM	"/usr/local/custom/bin/shm_param"
#define SHM_ALLOC	"/usr/local/custom/bin/shm_alloc"
#define K5_SAMPLE	"/usr/local/custom/bin/k5sample_store"
#define K5_SIM		"/usr/local/custom/bin/k5sim"
#define SHM_SEGDATA	"/usr/local/custom/bin/shm_segdata"
#define CUDA_FFT	"/usr/local/custom/bin/cuda_fft_xspec"
#define SPEC_VIEW	"/usr/local/custom/bin/shm_spec_view"
#define POWER_VIEW	"/usr/local/custom/bin/shm_power_view"
