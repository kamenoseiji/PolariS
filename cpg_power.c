// cpg_power.c	: Total Power Monitor
//
// Author : Seiji Kameno
// Created: 2013/12/23
#include "shm_k5data.inc"
#include <cpgplot.h>
#include <math.h>
#define	timeNum		128
#define	MAX(a, b)	((a) > (b)?(a):(b))
#define	MIN(a, b)	((a) < (b)?(a):(b))

int	cpg_power(
	struct SHM_PARAM	*param_ptr,	// Shared Parameter
	float	*bitPower)				// Pointer to bitPower
{
	float	xmin, xmax;				// Plot Window Range
	float	ymin, ymax;				// Plot Window Range
	int		index;
	int		timeIndex;				// time index
	int		st_index;				// Index for Sub-Stream
	int		nxwin, nywin;			// Number of Panels in X and Y
	int		nx_index, ny_index;		// Index for Panels
	int		oldestIndex;
	float	plotX[timeNum], plotY[timeNum];
	float	xwin_incr,	ywin_incr;	// Position Increment of Panels
	float	x_text, y_text;			// Text Drawing Position
	char	text[256];				// Text to Write

	cpgsch(0.5);
	for(timeIndex=0; timeIndex<timeNum; timeIndex++){
		plotX[timeIndex] = (float)(-timeIndex);
	};
	nxwin   = (int)sqrt((double)param_ptr->num_st);
	nywin   = (param_ptr->num_st + nxwin - 1)/nxwin;
	xwin_incr = 0.9 / (float)nxwin;
	ywin_incr = 0.9 / (float)nywin;

	oldestIndex = MIN(timeNum-1, param_ptr->current_rec);
	for(st_index=0; st_index<param_ptr->num_st; st_index++){

		nx_index	= st_index % nxwin;
		ny_index	= st_index / nxwin;

		memcpy(plotY, &bitPower[st_index* timeNum], timeNum* sizeof(float));

		/*-------- PLOT WINDOW --------*/
		xmin = plotX[oldestIndex];	xmax = plotX[0];
		ymin = 30.0;	ymax = -30.0;
		for(index=0; index<oldestIndex; index++){
			ymin = MIN(ymin, plotY[index]);
			ymax = MAX(ymax, plotY[index]);
		}
		ymin -= 1.0; ymax += 1.0;

		cpgsvp(	0.067+xwin_incr*nx_index, 0.067+xwin_incr*(nx_index+0.9),
				0.067+ywin_incr*ny_index, 0.067+ywin_incr*(ny_index+0.9));
		cpgswin(xmin, xmax, ymin, ymax);

		cpgsci(2);	cpgrect(xmin, xmax, ymin, ymax);
		cpgsci(0);	cpgbox("G", 0.0, 0, "G", 10.0, 0);
		cpgsci(1);	cpgbox(	"BCNTS", 0.0, 0, "BCNTS", 10.0, 10);
		cpgsci(4);	cpgpt( timeNum, plotX, plotY, 16);

		//-------- IF number
		x_text = xmin*0.2 + xmax*0.8; y_text = ymin*0.1 + ymax*0.9;
		sprintf(text, "IF = %d\0", st_index); cpgsci(3);	cpgtext( x_text, y_text, text );
		y_text = ymin*0.15 + ymax*0.85;
		sprintf(text, "%7.3f [dB]\0", plotY[0]); cpgsci(1);	cpgtext( x_text, y_text, text );
	}
	//-------- UTC
	x_text = xmin*0.3 + xmax*0.7; y_text = 0.05*ymin + 0.95* ymax;
	sprintf(text, "%04d %03d %02d:%02d:%02d\0", param_ptr->year, param_ptr->doy, param_ptr->hour, param_ptr->min, param_ptr->sec); cpgsci(1);	cpgtext( x_text, y_text, text );

	cpgebuf();

	return(0);
}
