# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>

int main ( int argc, char *argv[] )
{
  int M=500;
  int N=500;

  double diff;
  double err_tol = 0.001;
  int i;
  int itr;
  int itr_print;
  int j;
  double mean;
  double obtd_diff;
  double u[M][N];
  double w[M][N];
  double wtime;

  printf ( "\n" );
  printf ( "HEATED_PLATE_OPENMP\n" );
  printf ( "  Parallel version\n" );
  printf ( "A program to solve the steady state temperature distribution in a 2D rectangular plate\n" );
  
  printf ( "\n" );
  printf ( "  Spatial grid of %d by %d points.\n", M, N );
  printf ( "  The iteration will be repeated until the change is <= %e\n", err_tol ); 
  
  /*
   Set the Number of threads for parallelization
   */
  
  omp_set_num_threads(4);
   
#pragma omp parallel  
        #pragma omp master  
        {  
            printf("  Number of threads available =%d\n", omp_get_num_threads( ));  
        }  
/*
  Set the boundary values, which don't change. 
*/
  mean = 0.0;

#pragma omp parallel shared ( w ) private ( i, j )
  {
#pragma omp for
    for ( i = 1; i < M - 1; i++ )
    {
      w[i][0] = 100.0;
    }
#pragma omp for
    for ( i = 1; i < M - 1; i++ )
    {
      w[i][N-1] = 100.0;
    }
#pragma omp for
    for ( j = 0; j < N; j++ )
    {
      w[M-1][j] = 100.0;
    }
#pragma omp for
    for ( j = 0; j < N; j++ )
    {
      w[0][j] = 0.0;
    }
/*
  Average the boundary values, to come up with a reasonable
  initial value for the interior.
*/
#pragma omp for reduction ( + : mean )
    for ( i = 1; i < M - 1; i++ )
    {
      mean = mean + w[i][0] + w[i][N-1];
    }
#pragma omp for reduction ( + : mean )
    for ( j = 0; j < N; j++ )
    {
      mean = mean + w[M-1][j] + w[0][j];
    }
  }
/*
  OpenMP note:
  You cannot normalize MEAN inside the parallel region.  It
  only gets its correct value once you leave the parallel region.
  So we interrupt the parallel region, set MEAN, and go back in.
*/
  mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
  printf ( "\n" );
  printf ( "  Mean Temperature = %f\n", mean );
/* 
  Initialize the interior solution to the mean value.
*/
#pragma omp parallel shared ( mean, w ) private ( i, j )
  {
#pragma omp for
    for ( i = 1; i < M - 1; i++ )
    {
      for ( j = 1; j < N - 1; j++ )
      {
        w[i][j] = mean;
      }
    }
  }
/*
  iterate until the  new solution W differs from the old solution U
  by no more than err_tol.
*/
  itr = 0;
  itr_print = 1;
  printf ( "\n" );
  printf ( " Iteration  Change in temperature\n" );
  printf ( "\n" );
  wtime = omp_get_wtime ( );

  diff = err_tol;

  while ( err_tol <= diff )
  {
# pragma omp parallel shared ( u, w ) private ( i, j )
    {
/*
  Save the old solution in U.
*/
# pragma omp for
      for ( i = 0; i < M; i++ ) 
      {
        for ( j = 0; j < N; j++ )
        {
          u[i][j] = w[i][j];
        }
      }
/*
  Determine the new estimate of the solution at the interior points.
  The new solution W is the average of north, south, east and west neighbors.
*/
# pragma omp for
      for ( i = 1; i < M - 1; i++ )
      {
        for ( j = 1; j < N - 1; j++ )
        {
          w[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;
        }
      }
    }
/*
  C and C++ cannot compute a maximum as a reduction operation.

  Therefore, we define a private variable obtd_diff for each thread.
  Once they have all computed their values, we use a CRITICAL section
  to update DIFF.
*/
    diff = 0.0;
# pragma omp parallel shared ( diff, u, w ) private ( i, j, obtd_diff )
    {
      obtd_diff = 0.0;
# pragma omp for
      for ( i = 1; i < M - 1; i++ )
      {
        for ( j = 1; j < N - 1; j++ )
        {
          if ( obtd_diff < fabs ( w[i][j] - u[i][j] ) )
          {
            obtd_diff = fabs ( w[i][j] - u[i][j] );
          }
        }
      }
# pragma omp critical
      {
        if ( diff < obtd_diff )
        {
          diff = obtd_diff;
        }
      }
    }

    itr++;
    if ( itr == itr_print )
    {
      printf ( "  %8d  %f\n", itr, diff );
      itr_print = 2 * itr_print;
    }
  } 
  wtime = omp_get_wtime ( ) - wtime;

  printf ( "\n" );
  printf ( "  %8d  %f\n", itr, diff );
  printf ( "\n" );
  printf ( "  Error tolerance achieved.\n" );
  printf ( "  Wallclock time = %f\n", wtime );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HEATED_PLATE_OPENMP:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;    
}

