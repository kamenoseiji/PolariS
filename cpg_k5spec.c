/*****************************************************************
**	cpg_k5spec.c	: Plot Power Spaectrum of K5-recorded data	**
**																**
**	AUTHOR	: KAMENO Seiji										**
**	CREATED	: 2010/6/01											**
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cpgplot.h>
#include <math.h>
#include "vssphead.inc"

int	cpg_k5spec(
	char	*fname,					// File Name
	struct K5MODE	*k5mode,		// K-5 Recording Mode
	int		freq_num,				// Number of Freq. Channels
	float	*vis_max_ptr,			// Maximum Visibility
	float	**vis_ptr_ptr)			// Pointer of Visibility Data
{
	float	xmin, xmax;				/* Plot Window Range			*/
	float	ymin, ymax;				/* Plot Window Range			*/
	int		st_index;				/* Index for Sub-Stream			*/
	int		nxwin, nywin;			/* Number of Panels in X and Y	*/
	int		nx_index, ny_index;		/* Index for Panels				*/
	int		err_code;				/* Error Code					*/
	int		freq_index;				/* Frequency Points				*/
	float	*freq_ptr;				/* Pointer to Frequencies		*/
	float	*vis_amp_ptr;			/* Pointer to Spectral Amp.		*/
	float	freq_incr;				// Frequency Increment
	float	xwin_incr,	ywin_incr;	/* Position Increment of Panels	*/
	float	x_text, y_text;			/* Text Drawing Position		*/
	char	text[256];				/* Text to Write				*/
	char	pg_device[32];			/* PGPLOT Device Name			*/


	/*-------- Parameters --------*/
	freq_incr = (float)(k5mode->smp_clk / 2 / freq_num);

	err_code	= 32;
	cpgqinf("FILE", pg_device, &err_code);
	err_code = 0;

	/*-------- OPEN PGPLOT DEVICE --------*/
	if( strstr( pg_device, "cps")  ){
		cpgscrn(0, "White", &err_code);			/* COLOR DEFINISHON */
		cpgscrn(1, "Black", &err_code);			/* COLOR DEFINISHON */
		cpgscrn(2, "ivory", &err_code);			/* COLOR DEFINISHON */
		cpgscrn(3, "Blue", &err_code);			/* COLOR DEFINISHON */
		cpgscrn(4, "Pink", &err_code);			/* COLOR DEFINISHON */
	} else if( strstr( pg_device, "ps") == NULL ){
		cpgscrn(0, "DarkSlateGray", &err_code);	/* COLOR DEFINISHON */
		cpgscrn(1, "White", &err_code);			/* COLOR DEFINISHON */
		cpgscrn(2, "SlateGrey", &err_code);		/* COLOR DEFINISHON */
		cpgscrn(3, "Yellow", &err_code);		/* COLOR DEFINISHON */
		cpgscrn(4, "Cyan", &err_code);			/* COLOR DEFINISHON */
	} else {
		cpgscrn(0, "White", &err_code);			/* COLOR DEFINISHON */
		cpgscrn(1, "Black", &err_code);			/* COLOR DEFINISHON */
		cpgscrn(2, "LightGrey", &err_code);		/* COLOR DEFINISHON */
		cpgscrn(3, "Black", &err_code);			/* COLOR DEFINISHON */
		cpgscrn(4, "Black", &err_code);			/* COLOR DEFINISHON */
	}
	cpgeras();

	cpgbbuf();
	nxwin	= (int)sqrt((double)k5mode->stnum);
	nywin	= (k5mode->stnum + nxwin - 1)/nxwin;
	xwin_incr = 0.9 / (float)nxwin;
	ywin_incr = 0.9 / (float)nywin;

	cpgsvp( 0.0, 1.0, 0.0, 1.0 );
	cpgswin( 0.0, 1.0, 0.0, 1.0 );
	cpgsci(1);
	cpgsch(0.5);

	/*-------- OPEN PGPLOT DEVICE --------*/
	cpgsch(1.0);
	sprintf( text, "%s Bandpass", fname );
	cpgtext( 0.42, 0.975, text );
	cpgtext( 0.45, 0.025, "Frequency [Hz]" );

	cpgsch(0.5);
	for(st_index=0; st_index<k5mode->stnum; st_index++){

		nx_index	= st_index % nxwin;
		ny_index	= st_index / nxwin;

		/*-------- PLOT WINDOW --------*/
		xmin = - 0.5*freq_incr;
		xmax = xmin + ((double)freq_num - 0.5) * freq_incr;
		ymin = 0.0;		ymax = (float)(*vis_max_ptr)*1.2;
		cpgsvp(	0.067+xwin_incr*nx_index, 0.067+xwin_incr*(nx_index+0.9),
				0.067+ywin_incr*ny_index, 0.067+ywin_incr*(ny_index+0.9));
		cpgswin(xmin, xmax, ymin, ymax);

		/*-------- X-AXIS DATA --------*/
		vis_amp_ptr = vis_ptr_ptr[st_index];
		freq_ptr	= (float *)malloc( freq_num * sizeof(float) );
		for( freq_index=0; freq_index< freq_num; freq_index++){
			freq_ptr[freq_index] = xmin + freq_index * freq_incr;
		}

		cpgsci(2);	cpgrect(xmin, xmax, ymin, ymax);
		cpgsci(0);	cpgbox("G", 0.0, 0, "G", 0.0, 0);
		cpgsci(1);	cpgbox(	"BCNTS", 0.0, 0, "BCNTS", 0.0, 0 );
		cpgsci(3);	cpgline( freq_num, freq_ptr, vis_amp_ptr );

		x_text = xmin*0.2 + xmax*0.8;
		y_text = ymin*0.1 + ymax*0.9;
		sprintf(text, "IF = %d", st_index);
		cpgsci(3);	cpgtext( x_text, y_text, text );
		free(freq_ptr);

		vis_max_ptr++;
	}

	cpgebuf();

	return(0);
}
