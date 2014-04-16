// PolariBunch: A tool to bunch spectrum
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shm_k5data.inc"
#define	MAX(a,b)	a>b?a:b	// Larger Value
#define	FNAME	1
#define	BUNCH	2
#define STARTPP	3
#define ENDPP	4
#define ARGNUM	5

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
		strstr(fname, ".C.")?2:0);
}

int	bunchReal(
	int		length,			// IN: Length of input original data
	int		bunchNum,		// IN: Number of bunch
	float	*inData,		// IN: Pointer to input original data
	float	*outData)		// IN: Pointer to output bunched data
{
	int		out_ch_index;	// Index of output bunched channel
	int		chIndex;		// Index within accumulation
	int		outLen;			// Length of output bunched data
	double	accumData;

	outLen = length / bunchNum;
	for(out_ch_index=0; out_ch_index<outLen; out_ch_index ++){
		accumData = 0.0; 
		for(chIndex=0; chIndex<bunchNum; chIndex++){
			accumData += (double)inData[bunchNum* out_ch_index + chIndex];
		}
		outData[out_ch_index] = (float)( accumData / (double)bunchNum );
	}
	if(length % bunchNum){
		accumData = 0.0; 
		for(chIndex=bunchNum* outLen; chIndex<length; chIndex++){
			accumData += (double)inData[chIndex];
		}
		outData[outLen] = (float)( accumData / (double)(length % bunchNum) );
		outLen ++;
	}
	return( outLen );
}

int	bunchComplex(
	int		length,			// IN: Length of input original data
	int		bunchNum,		// IN: Number of bunch
	float	*inData,		// IN: Pointer to input original data
	float	*outData)		// IN: Pointer to output bunched data
{
	int		out_ch_index;	// Index of output bunched channel
	int		chIndex;		// Index within accumulation
	int		outLen;			// Length of output bunched data
	double	accumReal, accumImag;

	outLen = length / bunchNum;
	for(out_ch_index=0; out_ch_index<outLen; out_ch_index ++){
		accumReal = 0.0;  accumImag = 0.0;
		for(chIndex=0; chIndex<bunchNum; chIndex++){
			accumReal += (double)inData[2* (bunchNum* out_ch_index + chIndex)];
			accumImag += (double)inData[2* (bunchNum* out_ch_index + chIndex)+1];
		}
		outData[2* out_ch_index]   = (float)( accumReal / (double)bunchNum );
		outData[2* out_ch_index+1] = (float)( accumImag / (double)bunchNum );
	}
	if(length % bunchNum){
		accumReal = 0.0;  accumImag = 0.0;
		for(chIndex=bunchNum* outLen; chIndex<length; chIndex++){
			accumReal += (double)inData[2* chIndex];
			accumImag += (double)inData[2* chIndex+1];
		}
		outData[2* outLen]   = (float)( accumReal / (double)(length % bunchNum) );
		outData[2* outLen+1] = (float)( accumImag / (double)(length % bunchNum) );
		outLen ++;
	}
	return( outLen );
}

int	fileExtractBunch(
	char	*fname,			// IN: Input File Name
	int		bunchNum,		// IN: Number of channels to bunch
	int		StartPP,		// IN: Rec# to start extract
	int		EndPP )			// IN: Final Rec# to extract
{
	FILE	*file_ptr;				// File Pointer
	FILE	*file_out;				// File Pointer
	char	outName[24];			// Output File Name
	struct	SHM_PARAM	param;		// File Header
	float	*recdata;				// Power spectrum data record
	float	*outdata;				// Power spectrum data record
	int		recSize;				// Original Record size [bytes]
	int		outSize;				// Output Record size [bytes]
	int		inCHnum, outCHnum;		// In/Output Channel Num
	int		filetype;				// 0:undef, 1:Acorr, 2:Ccorr, 3:bitDist
	int		startSod, endSod;		// Start and End second of day
	int		startH, startM, startS, endH, endM, endS;
	int		index;

	//-------- Open File
	if((filetype = fileType(fname)) == 0){ return(-1);}
	if((file_ptr = fopen(fname, "r")) == NULL){ return(-1);}	// Open input file

	//-------- Read and Write Header
	fread(&param, sizeof(struct SHM_PARAM), 1, file_ptr);
	inCHnum = param.num_ch; outCHnum = inCHnum / bunchNum; if( inCHnum % bunchNum){	outCHnum ++;}
	switch(filetype){
		case 0:		recSize = 0;									break;		// Not-available file
		case 1:
				recSize = inCHnum* sizeof(float);							// Autocorr file
				outSize = outCHnum* sizeof(float);						break;
		case 2:
				recSize = inCHnum* sizeof(double);							// Crosscorr file
				outSize = outCHnum* sizeof(double);					break;
		case 3:		recSize = 0;									break;		// Bit distribution
	}
	if( recSize == 0 ){	perror("Unavailable Data Format"); return(-1);}
	recdata = malloc(recSize);
	outdata = malloc(outSize);
	startSod = hms2sod(param.hour, param.min, param.sec);
	endSod = startSod + EndPP; startSod += StartPP;
	sod2hms(startSod, &startH, &startM, &startS);
	sod2hms(endSod,   &endH,   &endM,   &endS);
	printf("Extract %02d:%02d:%02d - %02d:%02d:%02d  CHnum %d -> %d \n", startH, startM, startS, endH, endM, endS, inCHnum, outCHnum);

	//-------- Skip to the start position
	fseek(file_ptr, StartPP* recSize, SEEK_CUR);
	param.hour = startH;	param.min = startM;		param.sec = startS;
	sprintf(outName, "%04d%03d%02d%02d%02d.%c.%c%cB", param.year, param.doy, param.hour, param.min, param.sec, fname[14], fname[16], fname[17]);
	param.num_ch = outCHnum;
	printf("Output File Name = %s\n", outName);
	file_out = fopen(outName, "w");
	fwrite(&param, sizeof(struct SHM_PARAM), 1, file_out);

	//-------- Read and Write Records
	for(index=StartPP; index<=EndPP; index++){
		if(fread(recdata, recSize, 1, file_ptr) != 1){printf("File Read Error [%d]\n", index); break;}
		if( filetype == 1 ){bunchReal(inCHnum, bunchNum, recdata, outdata);	}
		if( filetype == 2 ){bunchComplex(inCHnum, bunchNum, recdata, outdata);	}
		fwrite(outdata, outSize, 1, file_out);
	}
	
	//-------- Close File
	free(recdata);
	free(outdata);
	fclose(file_ptr);
	fclose(file_out);

	return(0);
}

int main(
	int		argc,		// Number of Arguments
	char	**argv)		// Pointer to Arguments
{
	if( argc < ARGNUM ){
		printf("USAGE: PolariBunch [file name] [BUNCH] [Start PP] [End PP] !!\n");
		printf("  file name : input file name (e.g. 2014016154932.C.00)\n");
		printf("  BUNCH     : Bunching number \n");
		printf("  Start PP  : Start record number to extract\n");
		printf("  End   PP  : Final record number to extract\n");
		exit(-1);
	}
	fileExtractBunch(argv[FNAME], atoi(argv[BUNCH]), atoi(argv[STARTPP]), atoi(argv[ENDPP]));
	return(0);
}
