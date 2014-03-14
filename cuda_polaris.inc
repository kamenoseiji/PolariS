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
		vec_out[ix] += vec_in[ix].x *  vec_in[ix].x +  vec_in[ix].y *  vec_in[ix].y;
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
		vec_out[ix].x += vec_in_a[ix].x * vec_in_b[ix].x + vec_in_a[ix].y * vec_in_b[ix].y;
		vec_out[ix].y += vec_in_a[ix].y * vec_in_b[ix].x - vec_in_a[ix].x * vec_in_b[ix].y;
	}
}

__global__ void scalePowerSpec(
	float	*vec_in,		// Input vector to be accumulated
	float	scaleFact,		// Scaling Factor
	int		length)			// Length to accumulating vector
{
	int ix = blockIdx.x* blockDim.x + threadIdx.x;
	if( (ix >= 0) && (ix < length) ){
		 vec_in[ix] *= scaleFact;
	}
}

__global__ void scaleCrossSpec(
	float2	*vec_in,		// Input vector to be accumulated
	float	scaleFact,		// Scaling Factor
	int		length)			// Length to accumulating vector
{
	int ix = blockIdx.x* blockDim.x + threadIdx.x;
	if( (ix >= 0) && (ix < length) ){
		 vec_in[ix].x *= scaleFact;
		 vec_in[ix].y *= scaleFact;
	}
}

//-------- Format segment data using 4-bit quantized K5 data
__global__ void segform4bit(
	char	*k5data_in,		// Input K5 Segmanet Data
	float	*segdata,		// Output Segment Data
	int		length)			// Length to copy
{
	int ix = blockIdx.x* blockDim.x + threadIdx.x;
	unsigned char bitmask = 0x0f;				// 4-bit mask
	float wt4[] = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	
	if( (ix >= 0) && (ix < length) ){
		segdata[           ix] = wt4[(k5data_in[2*ix]       ) & bitmask];
		segdata[  length + ix] = wt4[(k5data_in[2*ix] >>  4 ) & bitmask];
		segdata[2*length + ix] = wt4[(k5data_in[2*ix+1]     ) & bitmask];
		segdata[3*length + ix] = wt4[(k5data_in[2*ix+1]>> 4 ) & bitmask];
	}
}

//-------- Format segment data using 8-bit quantized K5 data
__global__ void segform8bit(
	char	*k5data_in,		// Input K5 Segmanet Data
	float	*segdata,		// Output Segment Data
	int		length)			// Length to copy
{
	int ix = blockIdx.x* blockDim.x + threadIdx.x;
	unsigned char bitmask = 0xff;		// 8-bit mask
	float wt8[] = {
		-7.96875, -7.90625, -7.84375, -7.78125, -7.71875, -7.65625, -7.59375, -7.53125, -7.46875, -7.40625, -7.34375, -7.28125, -7.21875, -7.15625, -7.09375, -7.03125, 
		-6.96875, -6.90625, -6.84375, -6.78125, -6.71875, -6.65625, -6.59375, -6.53125, -6.46875, -6.40625, -6.34375, -6.28125, -6.21875, -6.15625, -6.09375, -6.03125, 
		-5.96875, -5.90625, -5.84375, -5.78125, -5.71875, -5.65625, -5.59375, -5.53125, -5.46875, -5.40625, -5.34375, -5.28125, -5.21875, -5.15625, -5.09375, -5.03125, 
		-4.96875, -4.90625, -4.84375, -4.78125, -4.71875, -4.65625, -4.59375, -4.53125, -4.46875, -4.40625, -4.34375, -4.28125, -4.21875, -4.15625, -4.09375, -4.03125, 
		-3.96875, -3.90625, -3.84375, -3.78125, -3.71875, -3.65625, -3.59375, -3.53125, -3.46875, -3.40625, -3.34375, -3.28125, -3.21875, -3.15625, -3.09375, -3.03125, 
		-2.96875, -2.90625, -2.84375, -2.78125, -2.71875, -2.65625, -2.59375, -2.53125, -2.46875, -2.40625, -2.34375, -2.28125, -2.21875, -2.15625, -2.09375, -2.03125, 
		-1.96875, -1.90625, -1.84375, -1.78125, -1.71875, -1.65625, -1.59375, -1.53125, -1.46875, -1.40625, -1.34375, -1.28125, -1.21875, -1.15625, -1.09375, -1.03125, 
		-0.96875, -0.90625, -0.84375, -0.78125, -0.71875, -0.65625, -0.59375, -0.53125, -0.46875, -0.40625, -0.34375, -0.28125, -0.21875, -0.15625, -0.09375, -0.03125, 
		 0.03125,  0.09375,  0.15625,  0.21875,  0.28125,  0.34375,  0.40625,  0.46875,  0.53125,  0.59375,  0.65625,  0.71875,  0.78125,  0.84375,  0.90625,  0.96875, 
		 1.03125,  1.09375,  1.15625,  1.21875,  1.28125,  1.34375,  1.40625,  1.46875,  1.53125,  1.59375,  1.65625,  1.71875,  1.78125,  1.84375,  1.90625,  1.96875, 
		 2.03125,  2.09375,  2.15625,  2.21875,  2.28125,  2.34375,  2.40625,  2.46875,  2.53125,  2.59375,  2.65625,  2.71875,  2.78125,  2.84375,  2.90625,  2.96875, 
		 3.03125,  3.09375,  3.15625,  3.21875,  3.28125,  3.34375,  3.40625,  3.46875,  3.53125,  3.59375,  3.65625,  3.71875,  3.78125,  3.84375,  3.90625,  3.96875, 
		 4.03125,  4.09375,  4.15625,  4.21875,  4.28125,  4.34375,  4.40625,  4.46875,  4.53125,  4.59375,  4.65625,  4.71875,  4.78125,  4.84375,  4.90625,  4.96875, 
		 5.03125,  5.09375,  5.15625,  5.21875,  5.28125,  5.34375,  5.40625,  5.46875,  5.53125,  5.59375,  5.65625,  5.71875,  5.78125,  5.84375,  5.90625,  5.96875, 
		 6.03125,  6.09375,  6.15625,  6.21875,  6.28125,  6.34375,  6.40625,  6.46875,  6.53125,  6.59375,  6.65625,  6.71875,  6.78125,  6.84375,  6.90625,  6.96875, 
		 7.03125,  7.09375,  7.15625,  7.21875,  7.28125,  7.34375,  7.40625,  7.46875,  7.53125,  7.59375,  7.65625,  7.71875,  7.78125,  7.84375,  7.90625,  7.96875}; 
	if( (ix >= 0) && (ix < length) ){
		segdata[           ix] = wt8[(k5data_in[4*ix]   ) & bitmask];
		segdata[  length + ix] = wt8[(k5data_in[4*ix+1] ) & bitmask];
		segdata[2*length + ix] = wt8[(k5data_in[4*ix+2] ) & bitmask];
		segdata[3*length + ix] = wt8[(k5data_in[4*ix+3] ) & bitmask];
	}
}
