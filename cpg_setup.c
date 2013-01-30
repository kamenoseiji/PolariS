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

	cpgsvp( 0.0, 1.0, 0.0, 1.0 );
	cpgswin( 0.0, 1.0, 0.0, 1.0 );
	cpgsci(1);

	/*-------- OPEN PGPLOT DEVICE --------*/
	cpgsch(1.0);
	cpgptxt( 0.5, 0.975, 0.0, 0.5, "PolariS" );
	cpgptxt( 0.5, 0.025, 0.0, 0.5, "Frequency [Hz]" );
	cpgptxt( 0.025, 0.5, 90.0, 0.5, "Relative Power [dB]" );

	return(0);
}
