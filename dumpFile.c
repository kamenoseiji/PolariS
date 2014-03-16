#include <stdio.h>
#include <stdlib.h>
#include "shm_k5data.inc"

#define BYTE	262144
main(
	int		argc,		// Number of arguments
	char	**argv)		// Pointer to arguments
{
	FILE	*file_ptr;	// Dump file
	struct	SHM_PARAM	param;
	unsigned char	k5data[BYTE];
	int		index, IF_index;

	IF_index = atoi(argv[2]);
	file_ptr = fopen(argv[1], "r");
	fread(&param, sizeof(struct SHM_PARAM), 1, file_ptr);		
	fread(k5data, BYTE, 1, file_ptr);

	for(index=0; index<BYTE; index += 4){
		printf("%03d\n", k5data[index + IF_index]);
		//printf("%02X\n", k5data[index + IF_index]);
	}

	fclose(file_ptr);
	return(0);
}
