//-------- K5 Header Dictionary
int	K5HEAD_CH[] = {1, 4};	// 0:1ch, 1:4ch
int	K5HEAD_FS[] = {     40,    100,     200,     500,
					  1000,   2000,    4000,    8000,
					 16000,  32000,   64000,  128000,
					256000, 512000, 1024000, 2048000};	// sampling frequency [kHz]
int	K5HEAD_QB[] = {1, 2, 4, 8};	// 0:1-bit, 1:2-bit, 2:4-bit, 3:8-bit

int	  k5bitshift[] = {0, 4, 8, 12};
float wt1[] = {-1.0, 1.0};
float wt2[] = {-3.0, -1.0, 1.0, 3.0};
float wt4[] = {-7.5/NFFT, -6.5/NFFT, -5.5/NFFT, -4.5/NFFT, -3.5/NFFT, -2.5/NFFT, -1.5/NFFT, -0.5/NFFT, 0.5/NFFT, 1.5/NFFT, 2.5/NFFT, 3.5/NFFT, 4.5/NFFT, 5.5/NFFT, 6.5/NFFT, 7.5/NFFT}; // 4-bit quantization
