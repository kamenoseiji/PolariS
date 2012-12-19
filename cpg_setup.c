// cpg_setup.c	: Set up for PGPLOT
//
// Author : Seiji Kameno
// Created: 2012/12/19
//
#include "shm_k5data.inc"
#include <cpgplot.h>
#include <math.h>

int	cpg_setup(struct SHM_PARAM	*param_ptr)
{
	float	xmin, xmax;				/* Plot Window Range			*/
	float	ymin, ymax;				/* Plot Window Range			*/
	int		nxwin, nywin;			/* Number of Panels in X and Y	*/
	int		nx_index, ny_index;		/* Index for Panels				*/
	int		err_code;				/* Error Code					*/
	float	xwin_incr,	ywin_incr;	/* Position Increment of Panels	*/

	/*-------- OPEN PGPLOT DEVICE --------*/
	cpgscrn(0, "DarkSlateGrey", &err_code);	/* COLOR DEFINISHON */
	cpgscrn(1, "White", &err_code);			/* COLOR DEFINISHON */
	cpgscrn(2, "SlateGrey", &err_code);		/* COLOR DEFINISHON */
	cpgscrn(3, "Yellow", &err_code);		/* COLOR DEFINISHON */
	cpgscrn(4, "Cyan", &err_code);			/* COLOR DEFINISHON */
	cpgeras();

	nxwin	= (int)sqrt((double)param_ptr->num_st);
	nywin	= (param_ptr->num_st + nxwin - 1)/nxwin;
	xwin_incr = 0.9 / (float)nxwin;
	ywin_incr = 0.9 / (float)nywin;

	cpgsvp( 0.0, 1.0, 0.0, 1.0 );
	cpgswin( 0.0, 1.0, 0.0, 1.0 );
	cpgsci(1);
	cpgsch(0.5);

	/*-------- OPEN PGPLOT DEVICE --------*/
	cpgsch(1.0);
	cpgtext( 0.42, 0.975, "PolariS" );
	cpgtext( 0.45, 0.025, "Frequency [Hz]" );

	return(0);
}
