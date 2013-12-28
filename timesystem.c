//	timesystem.c : functions for time system
//
//	Author : Seiji Kameno
//	Created: 2012/11/14
//
#include <stdio.h>
#include <stdlib.h>

//-------- Convert SoD (Second of Day) into hour, min, and second
int sod2hms(
	int	sod,		// Second of Day
	int	*hour,		// Hour
	int	*min,		// Min
	int	*sec)		// Sec
{
	*hour = sod / 3600;
	*min  = (sod % 3600) / 60;
	*sec  = (sod % 60);
	return(*sec);
}

