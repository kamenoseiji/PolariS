//	cuda_fft_xspec.c : FFT using CuFFT
//
//	Author : Seiji Kameno
//	Created: 2012/12/6
//
#include "shm_k5data.inc"
#include "k5dict.inc"
#include <cuda.h>
#include <cufft.h>

__device__ float2 complexMult(float2 a, float2 b)		// output a x b
{
	return make_float2( a.x* b.x - a.y* b.y, a.x* b.y + a.y* b.x );
}

__device__ float2 complexMultConj(float2 a, float2 b)		// output a x b*
{
	return make_float2( a.x* b.x + a.y* b.y,  a.y* b.x - a.x* b.y );
}

__device__ float complexMod( float2 a )				// output |a|^2
{
	return  a.x* a.x + a.y* a.y;
}

__global__ void complexMultConjVec(		// calculate a x b*
	float2	*vec_in_a,			// Input vector
	float2	*vec_in_b,			// Input vector
	float2	*vec_out,			// Output vector
	int		length)				// Vector length
{
	int tid = blockIdx.x* blockDim.x + threadIdx.x;
	if((tid >= 0) && (tid < length)){
		vec_out[tid] = complexMultConj(vec_in_a[tid], vec_in_b[tid]);
	}
}

__global__ void complexPowerVec(		// calculate a x a*
	float2	*vec_in,		// Input vector
	float	*vec_out,		// Output vector
	int		length)			// Number of elements
{
	int tid = blockIdx.x* blockDim.x + threadIdx.x;
	if((tid >= 0) && (tid < length)){
		vec_out[tid] = complexMod(vec_in[tid]);
	}
}

__global__ void accumReal(	// a <- a + b
	float	*vec_in_a,		// Accumuration Results
	float	*vec_in_b,		// to be accumulated
	int		length)
{
    int tid = blockIdx.x* blockDim.x + threadIdx.x;
    if((tid >= 0) && (tid < length)){
        vec_in_a[tid] += vec_in_b[tid];
    }
}

__global__ void accumComplex(	// a <- a + b
	float2	*vec_in_a,		// Accumuration Results
	float2	*vec_in_b,		// to be accumulated
	int		length)
{
    int tid = blockIdx.x* blockDim.x + threadIdx.x;
    if((tid >= 0) && (tid < length)){
        vec_in_a[tid].x += vec_in_b[tid].x;
        vec_in_a[tid].y += vec_in_b[tid].y;
    }
}

__global__ void accumPowerSpec(
	float2	*vec_in,		// Input vector to be accumulated
	float	*vec_out,		// In/Output vector to accumulate
	int		length)			// Length to accumulating vector
{
    int ix = blockIdx.x* blockDim.x + threadIdx.x;

	if( (ix >= 0) && (ix < length) ){
		vec_out[ix] += complexMod(vec_in[ix])/(length);
	}
}

__global__ void accumCrossSpec(
	float2	*vec_in_a,		// Input vector
	float2	*vec_in_b,		// Input vector
	float2	*vec_out,		// Output vector
	int		length)
{
    int ix = blockIdx.x* blockDim.x + threadIdx.x;

	if( (ix >= 0) && (ix < length) ){
		vec_out[ix].x += complexMultConj(vec_in_a[ix], vec_in_b[ix]).x / length;
		vec_out[ix].y += complexMultConj(vec_in_a[ix], vec_in_b[ix]).y / length;
	}
}

//-------- Format segment data using 4-bit quantized K5 data
__global__ void segform(
	short	*k5data_in,		// Input K5 Segmanet Data
	float	*segdata,		// Output Segment Data
	int		length)			// Length to copy
{
	int ix = blockIdx.x* blockDim.x + threadIdx.x;
	short bitshift[]={0, 4, 8, 12};
	short bitmask = 0x000f;
	float wt4[] = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	
	if( (ix >= 0) && (ix < length) ){
		segdata[           ix] = wt4[(k5data_in[ix] >> bitshift[0]) & bitmask];
		segdata[  length + ix] = wt4[(k5data_in[ix] >> bitshift[1]) & bitmask];
		segdata[2*length + ix] = wt4[(k5data_in[ix] >> bitshift[2]) & bitmask];
		segdata[3*length + ix] = wt4[(k5data_in[ix] >> bitshift[3]) & bitmask];
	}
}

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
