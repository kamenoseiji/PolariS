#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define	MAX_LOOP	10		// Maximum number of iterations
#define	MAX(a,b)	a>b?a:b	// Larger Value

int prob4bit(
	double *param,	// IN: Gaussian mean and sigma
	double *prob)	// OUT:Probabilities in 16 levels
{
	int		index;	// General purpose index
	double	volt[] = {-7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

	for(index = 0; index < 15; index ++){ volt[index] *= param[0]; }		// Scale threshold
	//-------- Calculate probabilities
	prob[0] = 0.5* (erf(M_SQRT1_2*(volt[0] - param[1])) + 1.0);
	for(index = 1; index < 14; index ++){
		prob[index] = 0.5*(erf(M_SQRT1_2*(volt[index] - param[1])) - erf(M_SQRT1_2*(volt[index-1] - param[1])));
	}
	prob[15] = 0.5* (1.0 - erf(M_SQRT1_2*(volt[14] - param[1])));
	return(0);
}

int initGauss4bit(
	double	*prob,		// IN: Probabilities in 16 levels
	double	*param)		// OUT:Estimated parameters
{
	double	Vweight[] = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	double	Pweight[] = {56.25, 42.25, 30.25, 20.25, 12.25, 6.25, 2.25, 0.25, 0.25, 2.25, 6.25, 12.25, 20.25, 30.25, 42.25, 56.25};
	double	Average=0.0;
	double	Variance=0.0;
	int		index;			// General purpose index

	for(index=0; index<16; index++){
		Average += Vweight[index]* prob[index];
		Variance += Pweight[index]* prob[index];
	}
	param[0] = 1.0/sqrt(Variance); param[1] = Average* param[0];
	return(0);
}

int gauss4bit(
	unsigned int *nsample,	// IN : number of samples in 16 levels
	double	*param,			// OUT: Gaussian parameters 
	double	*param_err)		// OUT: Gaussian parameters 
{
	int		index;			// General index for loops
	int		loop_counter = 0;		// Loop Counter
	unsigned int	total_sample = 0;	// Total number of samples
	double	pwp[2][2];		// Weighted Partial Matrix
	double	prob[16];		// Probability in each state
	double	pred[16];		// Predicted probability in each state
	double	weight[16];		// Weight for Probability 
	double	resid[16];		// residual from trial Gaussian
	double	erfDeriv[16];	// Vector to produce partial matrix
	double	WerfDeriv[16];	// Vector to produce partial matrix
	double	wpr[2];			// WPr vector
	double	solution[2];	// correction vector for parameters
	double	expArg;			// 
	double	det;			// determinant of the partial matrix
	double	norm;			// Norm of the correction vector
	double	epsz;			// criterion for convergence

	//-------- Calculate probability in each state
	for(index=0; index<16; index++){ total_sample += nsample[index]; }	
	for(index=0; index<16; index++){ prob[index] = (double)nsample[index] / (double)total_sample; }	
	for(index=0; index<16; index++){ weight[index] = (double)nsample[index] / ((1.0 - prob[index])* (1.0 - prob[index]))  ; }	
	epsz = MAX(1.0e-6 / (total_sample* total_sample), 1.0e-29);		// Convergence

	initGauss4bit(prob, param);	// Initial parameter

	while(1){				// Loop for Least-Square Fit
		//-------- Calculate Residual Probability
		prob4bit(param, pred);
		for(index=0; index<16; index++){ resid[index] = prob[index] - pred[index]; }

		//-------- Calculate Elements of partial matrix
		erfDeriv[0] = 0.0; WerfDeriv[0] = 0.0;
		for(index=1; index<16; index++){
			expArg = ((double)index - 8.0)* param[0] - param[1];
			erfDeriv[index] = exp( -0.5* expArg* expArg);
			WerfDeriv[index] = ((double)index - 8.0)* erfDeriv[index];
		}
		for(index=0; index<15; index++){
			 erfDeriv[index] = 0.5* M_2_SQRTPI* M_SQRT1_2*( -erfDeriv[index + 1] +  erfDeriv[index]);
			WerfDeriv[index] = 0.5* M_2_SQRTPI* M_SQRT1_2*( WerfDeriv[index + 1] - WerfDeriv[index]);
		}
		erfDeriv[15] = 0.5* M_2_SQRTPI* M_SQRT1_2* erfDeriv[15];
		WerfDeriv[15] = -0.5* M_2_SQRTPI* M_SQRT1_2* WerfDeriv[15];

		//-------- Partial Matrix
		memset(pwp, 0, sizeof(pwp)); memset(wpr, 0, sizeof(wpr));
		for(index=0; index<16; index++){
			pwp[0][0] += (WerfDeriv[index]* WerfDeriv[index]* weight[index]);
			pwp[0][1] += (WerfDeriv[index]*  erfDeriv[index]* weight[index]);
			pwp[1][1] += ( erfDeriv[index]*  erfDeriv[index]* weight[index]);
			wpr[0] += (weight[index]* WerfDeriv[index]* resid[index]);
			wpr[1] += (weight[index]*  erfDeriv[index]* resid[index]);
		}
		pwp[1][0] = pwp[0][1];

		//-------- Solutions for correction vectors
		det = pwp[0][0]* pwp[1][1] - pwp[1][0]* pwp[0][1];
		if( fabs(det) < epsz ){	return(1);	}						// Too small determinant -> Error
		solution[0] = (pwp[1][1]* wpr[0] - pwp[0][1]* wpr[1])/ det;
		solution[1] =(-pwp[1][0]* wpr[0] + pwp[0][0]* wpr[1])/ det;

		//-------- Correction
		param[0] += solution[0];	param[1] += solution[1];	norm = solution[0]*solution[0] + solution[1]*solution[1];

		//-------- Converged?
		loop_counter ++;
		if( norm < epsz ){	break;	}
		if( loop_counter > MAX_LOOP ){	return(1);	}		// Doesn't converge
	}	// End of iteration loop

	//-------- Standard Error
	param_err[0] = sqrt(pwp[1][1] / det);
	param_err[1] = sqrt(pwp[0][0] / det);
	return(0);
}
