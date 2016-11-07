/*
 *	This file is part of qpDUNES.
 *
 *	qpDUNES -- A DUal NEwton Strategy for convex quadratic programming.
 *	Copyright (C) 2012 by Janick Frasch, Hans Joachim Ferreau et al. 
 *	All rights reserved.
 *
 *	qpDUNES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpDUNES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpDUNES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file examples/example1.c
 *	\author Janick Frasch, Hans Joachim Ferreau
 *	\version 1.0beta
 *	\date 2012
 *
 *	Very simple example for testing qpDUNES.
 */



#include <qpDUNES.h>

#define INFTY 1.0e12

int main( )
{	
	unsigned int i;
	boolean_t isLTI;
	return_t statusFlag;
	
	#ifdef __MEASURE_TIMINGS__
	real_t totalTime, start, ende;
	#endif

	qpOptions_t options = qpDUNES_setupDefaultOptions();
	options.cyclicReduction = QPDUNES_TRUE;
	options.nThreads = 4;
	
	/* set dimensions */
	unsigned int nI = 65535;		/* number of stages */
	unsigned int nX = 30;		/* number of states */
	unsigned int nU = 29;		/* number of controls */
	unsigned int* nD = 0;		/* number of affine constraints */
	
	
	/* specify problem data */
	double* Q = (double*)qpDUNES_calloc( nX*nX,sizeof(double) );
	for ( i = 0; i < nX*nX; i = i+nX+1 )
			Q[i] = 1.0;
	
	double* R = (double*)qpDUNES_calloc( nU*nU,sizeof(double) );
	for ( i = 0; i < nU*nU; i = i+nU+1 )
			R[i] = 1.0;
		
	double *S=0;
	
	double* P = Q;
	
	double* A = (double*)qpDUNES_calloc( nX*nX,sizeof(double) );
	for ( i = 0; i < nX*nX; i = i+nX+1 )
			A[i] = 1.0;
		
	double* B = (double*)qpDUNES_calloc( nX*nU,sizeof(double) );
	for ( i = 0; i < nU*nU; i = i+nU+1 )
			B[i] = 1.0;
		
	for ( i = nU*nU; i < nX*nU; ++i )
			B[i] = 1.0;
			
	double* c = (double*)qpDUNES_calloc( nX,sizeof(double) );
	for ( i = 0; i < nX; ++i )
			c[i] = 5.0;
	
	double* xLow = (double*)qpDUNES_calloc( nX,sizeof(double) );
	for ( i = 0; i < nX; ++i )
			xLow[i] = -INFTY;
	double* xUpp = (double*)qpDUNES_calloc( nX,sizeof(double) );
	for ( i = 0; i < nX; ++i )
			xUpp[i] = INFTY;
	double* uLow = (double*)qpDUNES_calloc( nU,sizeof(double) );
	for ( i = 0; i < nU; ++i )
			uLow[i] = -INFTY;
	double* uUpp = (double*)qpDUNES_calloc( nU,sizeof(double) );
	for ( i = 0; i < nU; ++i )
			uUpp[i] = INFTY;
	
	
	/* qpData struct */
	qpData_t qpData;

	/* memory allocation */
	
	/*** solve sequential (without CR)***/
	statusFlag = qpDUNES_setup( &qpData, nI, nX, nU, nD, 0 );	/* passing 0 in the last argument sets the default QP options */
	if (statusFlag != QPDUNES_OK)
	{
		printf("Setup of the QP solver failed\n");
		return (int)statusFlag;
	}
	
	/* manual setup of intervals */
	for( i=0; i<nI; ++i )
	{
		statusFlag = qpDUNES_setupSimpleBoundedInterval(  &qpData, qpData.intervals[i],Q,R,S, A,B,c, xLow,xUpp,uLow,uUpp );
		if (statusFlag != QPDUNES_OK)
		{
			printf("Setup of the QP solver failed\n");
			return (int)statusFlag;
		}
	}
	statusFlag = qpDUNES_setupSimpleBoundedInterval(  &qpData, qpData.intervals[nI], P,0,0, 0,0,0, xLow,xUpp,0,0 );
	if (statusFlag != QPDUNES_OK)
	{
		printf("Setup of the QP solver failed\n");
		return (int)statusFlag;
	}

	statusFlag = qpDUNES_setupAllLocalQPs( &qpData, isLTI=QPDUNES_TRUE );	/* determine local QP solvers and set up auxiliary data */
	if (statusFlag != QPDUNES_OK)
	{
		printf("Setup of the QP solver failed\n");
		return (int)statusFlag;
	}

	/* solve problem */
	printf("SOLVING QP (SEQUENTIAL)\n\n");
	
	#ifdef __MEASURE_TIMINGS__
	start = getTime();
	#endif
	
	statusFlag = qpDUNES_solve( &qpData );
	
	#ifdef __MEASURE_TIMINGS__
	ende = getTime();
	#endif
	
	if (statusFlag != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND)
	{
		printf("QP solver failed. The error code is: %d\n", statusFlag);
		return (int)statusFlag;
	}
	
	#ifdef __MEASURE_TIMINGS__
	totalTime = ende - start;
	printf("total time: %.10f s\n",totalTime);
	#endif
	
	/* write out solution */
	/*for( i=0; i<nI; ++i )
	{
		qpDUNES_printMatrixData( qpData.intervals[i]->z.data, 1, nX+nU, "z[%d]:", i );
	}
	qpDUNES_printMatrixData( qpData.intervals[nI]->z.data, 1, nX, "z[%d]:", i );*/
	
	qpDUNES_cleanup( &qpData );
	printf("\n\n");
	
	/*** solve parallel ***/
	statusFlag = qpDUNES_setup( &qpData, nI, nX, nU, nD, &options );	/* passing 0 in the last argument sets the default QP options */
	if (statusFlag != QPDUNES_OK)
	{
		printf("Setup of the QP solver failed\n");
		return (int)statusFlag;
	}
	
	/* manual setup of intervals */
	for( i=0; i<nI; ++i )
	{
		statusFlag = qpDUNES_setupSimpleBoundedInterval(  &qpData, qpData.intervals[i],Q,R,S, A,B,c, xLow,xUpp,uLow,uUpp );
		if (statusFlag != QPDUNES_OK)
		{
			printf("Setup of the QP solver failed\n");
			return (int)statusFlag;
		}
	}
	statusFlag = qpDUNES_setupSimpleBoundedInterval(  &qpData, qpData.intervals[nI], P,0,0, 0,0,0, xLow,xUpp,0,0 );
	if (statusFlag != QPDUNES_OK)
	{
		printf("Setup of the QP solver failed\n");
		return (int)statusFlag;
	}

	statusFlag = qpDUNES_setupAllLocalQPs( &qpData, isLTI=QPDUNES_TRUE );	/* determine local QP solvers and set up auxiliary data */
	if (statusFlag != QPDUNES_OK)
	{
		printf("Setup of the QP solver failed\n");
		return (int)statusFlag;
	}

	/* solve problem */
	printf("SOLVING QP (PARALLEL)\n\n");
	
	#ifdef __MEASURE_TIMINGS__
	start = getTime();
	#endif
	
	statusFlag = qpDUNES_solve( &qpData );
	
	#ifdef __MEASURE_TIMINGS__
	ende = getTime();
	#endif
	
	if (statusFlag != QPDUNES_SUCC_OPTIMAL_SOLUTION_FOUND)
	{
		printf("QP solver failed. The error code is: %d\n", statusFlag);
		return (int)statusFlag;
	}
	
	#ifdef __MEASURE_TIMINGS__
	totalTime = ende - start;
	printf("total time: %.10f s\n",totalTime);
	#endif
	
	/* write out solution */
	/*for( i=0; i<nI; ++i )
	{
		qpDUNES_printMatrixData( qpData.intervals[i]->z.data, 1, nX+nU, "z[%d]:", i );
	}
	qpDUNES_printMatrixData( qpData.intervals[nI]->z.data, 1, nX, "z[%d]:", i );*/
	
	qpDUNES_cleanup( &qpData );
	
	free (Q);
	free (R);
	free (B);
	free (A);
	free(c);
	free (xLow);
	free (xUpp);
	free (uLow);
	free (uUpp);
	
	printf( "example1 done.\n" );
	
	return 0;
}


/*
 *	end of file
 */
