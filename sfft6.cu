#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cufft.h>
#include <string.h>
#include <math.h>

#define M (int)pow(2,18)
#define MM (int)pow(2,17)
#define batch 1
#define REAL 10
#define	IMAGINARY 100

__device__
float2 complex_mult(float2 a,float2 b)
{
	return make_float2(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}

__device__
float2 complex_conjugate(float2 a)
{
	return make_float2(a.x,-a.y);
}

__global__
void complex_mult_conj(
	float2 *ia,	//input data1
	float2 *ib,	//input data2
	float2 *o,	//output data
	int size	//size of elements
	)
{
        int tid=blockIdx.x*blockDim.x+threadIdx.x;
	if((tid>=0)&&(tid<size)){
                o[tid]=complex_mult(ia[tid],complex_conjugate(ib[tid]));
	}
}

__global__
void complex_choose(
	float2 *a,	//input data
	float *b,	//output data
	int size,
	int type)
{
	int tid=blockIdx.x*blockDim.x+threadIdx.x;
	if((tid>=0)&&(tid<size)){
		if(type==REAL)	b[tid]=a[tid].x;
		else if(type==IMAGINARY)	b[tid]=a[tid].y;
	}
}

void f2wData(
	char *filename,
	float2 *data,
	int size)
{
        FILE *write;
        char str[128];

        write=fopen(filename,"w");
        if(write==NULL){
                printf("fopen Error\n");
                exit(0);
        }

        for(int i=0;i<size;i++){
                sprintf(str,"%f %f\n",data[i].x,data[i].y);
                fprintf(write,str);
        }
        fclose(write);
}

void fwData(
	char *filename,
	float *data,
	int size)
{
        FILE *write;
        char str[128];

        write=fopen(filename,"w");
        if(write==NULL){
                printf("fopen Error\n");
                exit(0);
        }

        for(int i=0;i<size;i++){
                sprintf(str,"%f\n",data[i]);
                fprintf(write,str);
        }
        fclose(write);
}

void makedata(
	float *data,
	int size,
	unsigned param)
{
	srand(param);
        for(int i=0;i<size;i++){
                data[i]=(float)(rand()%16);
        }
        return;
}

void SetUp(
        float *ia,
        float *ib,
        float *oa,
        float *ob,
        float *oabr,
	float *oabi,
        int size)
{
        cudaError_t err[6];

        err[0]=cudaMalloc((void **)&ia,2*size*sizeof(float));
        err[1]=cudaMalloc((void **)&ib,2*size*sizeof(float));
        err[2]=cudaMalloc((void **)&oa,size*sizeof(float));
        err[3]=cudaMalloc((void **)&ob,size*sizeof(float));
        err[4]=cudaMalloc((void **)&oabr,size*sizeof(float));
	err[5]=cudaMalloc((void **)&oabi,size*sizeof(float));

        for(int i=0;i<6;++i){
                if(err[i]!=cudaSuccess){
                        printf("error in setup err[%d] = code %d",i,err[i]);
                }
        }
}

void fftandcross(
	cufftHandle a,
	cufftHandle b,
	float *ia,
	float *ib,
	float *oa,
	float *ob,
	float *oabr,
	float *oabi,
	int size)
{
	float2 *ta,*tb,*tab;
	cudaMalloc((void **)&ta,size*sizeof(float2)*batch);
	cudaMalloc((void **)&tb,size*sizeof(float2)*batch);
	cudaMalloc((void **)&tab,size*sizeof(float2)*batch);

	cufftExecR2C(a,(cufftReal*)ia,ta);
        cufftExecR2C(b,(cufftReal*)ib,tb);
	complex_mult_conj<<<1024,512>>>(ta,ta,ta,size);
        complex_mult_conj<<<1024,512>>>(tb,tb,tb,size);
	complex_mult_conj<<<1024,512>>>(ta,tb,tab,size);
	cudaThreadSynchronize();

	complex_choose<<<1024,512>>>(ta,oa,size,REAL);
        complex_choose<<<1024,512>>>(tb,ob,size,REAL);
        complex_choose<<<1024,512>>>(tab,oabr,size,REAL);
	complex_choose<<<1024,512>>>(tab,oabi,size,IMAGINARY);
	cudaThreadSynchronize();

	cudaFree(ta);
	cudaFree(tb);
	cudaFree(tab);
}

void quit_fft(
	cufftHandle a,
	cufftHandle b,
	float *ia,
	float *ib,
	float *oa,
	float *ob,
	float *oabr,
	float *oabi)
{
	cufftDestroy(a);
        cufftDestroy(b);
        cudaFree(ia);
        cudaFree(ib);
        cudaFree(oa);
        cudaFree(ob);
        cudaFree(oabr);
        cudaFree(oabi);
}

void copy(float *a,float *da,int size)
{
	cudaMemcpy(a,da,size*sizeof(float),cudaMemcpyDeviceToHost);
}
int main(void){
        float iVdata[M],iHdata[M];
	float oVdata[MM],oHdata[MM],oHVdatar[MM],oHVdatai[MM];
        int data_size=M;
        int err,counter=1000;

	float elapsed_time_ms=0.0f;
	cudaEvent_t start, stop;

        cufftHandle planV,planH;
        float *gpuiV,*gpuiH,*odH,*odV,*odHVr,*odHVi;

//################# SET UP and MAKE DATA ###################
	cudaMalloc((void **)&gpuiV,data_size*sizeof(float));
        cudaMalloc((void **)&gpuiH,data_size*sizeof(float));
	cudaMalloc((void **)&odHVr,(int)data_size/2*sizeof(float)*batch);
	cudaMalloc((void **)&odHVi,(int)data_size/2*sizeof(float)*batch);
	cudaMalloc((void **)&odH,(int)data_size/2*sizeof(float)*batch);
	cudaMalloc((void **)&odV,(int)data_size/2*sizeof(float)*batch);
	cufftPlan1d(&planV,data_size,CUFFT_R2C,batch);
	cufftPlan1d(&planH,data_size,CUFFT_R2C,batch);

//	SetUp(gpuiH,gpuiV,odH,odV,odHVr,odHVi,(int)data_size/2);
        makedata(iVdata,data_size,10);
        makedata(iHdata,data_size,100);
	
	cudaEventCreate( &start );
	cudaEventCreate( &stop  );
	cudaEventRecord( start, 0 );

	for(int kk=0;kk<counter;kk++){
//################ COPY DATA from Host To Device ##########
        cudaMemcpy(gpuiV,iVdata,data_size*sizeof(float),cudaMemcpyHostToDevice);
        cudaMemcpy(gpuiH,iHdata,data_size*sizeof(float),cudaMemcpyHostToDevice);

//################ FFT and CROSS SPECTRUM  #####################
	fftandcross(planV,planH,gpuiV,gpuiH,odV,odH,odHVr,odHVi,(int)data_size/2);

//############### COPY DATA from Device To Host ########
	cudaMemcpy(oVdata,odV,data_size/2*batch*sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(oHdata,odH,data_size/2*batch*sizeof(float),cudaMemcpyDeviceToHost);
        cudaMemcpy(oHVdatar,odHVr,data_size/2*batch*sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemcpy(oHVdatai,odHVi,data_size/2*batch*sizeof(float),cudaMemcpyDeviceToHost);


	}
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );
	cudaEventElapsedTime( &elapsed_time_ms, start, stop );
//############### OUT PUT DATA ##########################

	printf("time of fft for 1000times = %f [msec]\n",elapsed_time_ms); 
	fwData("inH.dat",iHdata,data_size);
	fwData("inV.dat",iVdata,data_size);
	fwData("VVr.dat",oVdata,(int)data_size/2);
	fwData("HHr.dat",oHdata,(int)data_size/2);
	fwData("HVr.dat",oHVdatar,(int)data_size/2);
	fwData("HVi.dat",oHVdatai,(int)data_size/2);

//############## Destroy Memory #################
	quit_fft(planV,planH,gpuiV,gpuiH,odH,odV,odHVr,odHVi);
        cudaEventDestroy( start );
	cudaEventDestroy( stop );
	return 0;
}
