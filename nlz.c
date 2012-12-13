//	nlz.c : Number of Leading Zero
//
#include <stdio.h>

//-------- Number of Leading Zero
// e.g. If input is x = 0001010010110001, nlz(x) returns 3 because
//      there are three zeros at the head of the bit series.
//
int nlz( unsigned int x)
{
	int	n;		// Output = Number of Leading Zero
	union {
		unsigned long long	as_uint64;
		double				as_double;
	} data;

	data.as_double = (double)x + 0.5;
	n = 1054 - (data.as_uint64 >> 52);
	return(n);
}

//-------- Bit count for 32-bit integer
// See http://mono-comp.com/tech/programming/hispeed-bit-count2/ for the algorhithm
int bitcount32(long bits)
{
	bits = (bits & 0x55555555) + (bits >> 1 & 0x55555555);
	bits = (bits & 0x33333333) + (bits >> 2 & 0x33333333);
	bits = (bits & 0x0f0f0f0f) + (bits >> 4 & 0x0f0f0f0f);
	bits = (bits & 0x00ff00ff) + (bits >> 8 & 0x00ff00ff);
	return (bits & 0x0000ffff) + (bits >>16 & 0x0000ffff);
}

//-------- Bit count for 64-bit integer
// See http://mono-comp.com/tech/programming/hispeed-bit-count2/ for the algorhithm
int bitcount64(long long bits)
{
	bits = (bits & 0x5555555555555555L) + (bits >> 1 & 0x5555555555555555L);
	bits = (bits & 0x3333333333333333L) + (bits >> 2 & 0x3333333333333333L);
	bits = (bits & 0x0f0f0f0f0f0f0f0fL) + (bits >> 4 & 0x0f0f0f0f0f0f0f0fL);
	bits = (bits & 0x00ff00ff00ff00ffL) + (bits >> 8 & 0x00ff00ff00ff00ffL);
	bits = (bits & 0x0000ffff0000ffffL) + (bits >>16 & 0x0000ffff0000ffffL);
	return (bits & 0x00000000ffffffffL) + (bits >>32 & 0x00000000ffffffffL);
}

//-------- Round down to 2^n
int	pow2round(unsigned int x)
{
	return(0x0001 << (31 - nlz(x)));
}
