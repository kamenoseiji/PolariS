#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shm_k5data.inc"
#define	MAX(a,b)	a>b?a:b	// Larger Value
#define	FNAME	1
#define STARTPP	2
#define ENDPP	3
#define ARGNUM	4

int	hms2sod(
	int		hour,		// IN: Hour
	int		min,		// IN: Minute
	int		sec)		// IN: Sec
{
	return( sec + 60*(min + 60*hour));
}

int	sod2hms(
	int		sod,		// IN: Second of the day
	int		*hour,		// Hour
	int		*min,		// Minute
	int		*sec)		// Sec
{
	*hour = sod / 3600;
	*min  = (sod % 3600) / 60;
	*sec  = sod % 60;
	return(sod - hms2sod(*hour, *min, *sec));
}


int	fileType( char	*fname)
{
	return(
		strstr(fname, ".A.")?1:
		strstr(fname, ".C.")?2:
		strstr(fname, ".P.")?3:0);
}

int	fileExtract(
	char	*fname,			// IN: Input File Name
	int		StartPP,		// IN: Rec# to start extract
	int		EndPP )			// IN: Final Rec# to extract
{
	FILE	*file_ptr;				// File Pointer
	FILE	*file_out;				// File Pointer
	char	outName[24];			// Output File Name
	struct	SHM_PARAM	param;		// File Header
	float	*recdata;				// Power spectrum data record
//	int		recSize[] = {0, NFFT2* sizeof(float), 2* NFFT2* sizeof(float), 16* sizeof(int)};		// Record size [bytes]
	int		recSize;				// Record size [bytes]
	int		filetype;				// 0:undef, 1:Acorr, 2:Ccorr, 3:bitDist
	int		startSod, endSod;		// Start and End second of day
	int		startH, startM, startS, endH, endM, endS;
	int		index;

	StartPP = MAX(StartPP, 1);		// Avoid StartPP==0 (duplication of file name)

	//-------- Open File
	if((filetype = fileType(fname)) == 0){ return(-1);}
	if((file_ptr = fopen(fname, "r")) == NULL){ return(-1);}	// Open input file

	//-------- Read and Write Header
	fread(&param, sizeof(struct SHM_PARAM), 1, file_ptr);
	switch(filetype){
		case 0:		recSize = 0;									break;		// Not-available file
		case 1:		recSize = param.num_ch* sizeof(float);			break;		// Autocorr file
		case 2:		recSize = param.segLen* sizeof(float);			break;		// Crosscorr file
		case 3:		recSize = (0x01 << param.qbit)* sizeof(int);	break;		// Bit distribution
	}
	recdata = malloc(recSize);
	startSod = hms2sod(param.hour, param.min, param.sec);
	endSod = startSod + EndPP; startSod += StartPP;
	sod2hms(startSod, &startH, &startM, &startS);
	sod2hms(endSod,   &endH,   &endM,   &endS);
	printf("Extract %02d:%02d:%02d - %02d:%02d:%02d\n", startH, startM, startS, endH, endM, endS);

	//-------- Skip to the start position
	fseek(file_ptr, StartPP* recSize, SEEK_CUR);
	param.hour = startH;	param.min = startM;		param.sec = startS;
	sprintf(outName, "%04d%03d%02d%02d%02d.%c.%c%c", param.year, param.doy, param.hour, param.min, param.sec, fname[14], fname[16], fname[17]);
	printf("Output File Name = %s\n", outName);
	file_out = fopen(outName, "w");
	fwrite(&param, sizeof(struct SHM_PARAM), 1, file_out);

	//-------- Read and Write Records
	for(index=StartPP; index<=EndPP; index++){
		if(fread(recdata, recSize, 1, file_ptr) != 1){printf("File Read Error [%d]\n", index); break;}
		fwrite(recdata, recSize, 1, file_out);
	}
	
	//-------- Close File
	free(recdata);
	fclose(file_ptr);
	fclose(file_out);

	return(0);
}

int main(
	int		argc,		// Number of Arguments
	char	**argv)		// Pointer to Arguments
{
	if( argc < ARGNUM ){
		printf("USAGE: PolariSplit [file name] [Start PP] [End PP] !!\n");
		printf("  file name : input file name (e.g. 2014016154932.C.00)\n");
		printf("  Start PP  : Start record number to extract\n");
		printf("  End   PP  : Final record number to extract\n");
		exit(-1);
	}
	fileExtract(argv[FNAME], atoi(argv[STARTPP]), atoi(argv[ENDPP]));
	return(0);
}
