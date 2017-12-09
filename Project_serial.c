# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>

void central(int M,int N,double w[M][N],double u[M][N])
{  
    int i,j;
    for ( i = 1; i < M - 1; i++ )
    {
     for ( j = 1; j < N - 1; j++ )
     {
         w[i][j] = (u[i-1][j]+u[i+1][j]+u[i][j-1]+u[i][j+1])/4.0;
     }
    }
}

void update(int M,int N,double w[M][N],double u[M][N])
{
    int i,j;
    for ( i = 0; i < M; i++ ) 
    {
      for ( j = 0; j < N; j++ )
      {
        u[i][j] = w[i][j];
       }
    }
}

int main ( int argc, char *argv[] )
/******************************************************************************/
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
  double my_diff;
  double u[M][N];
  double w[M][N];
  double wtime;

  printf ( "\n" );
  printf ( "HEATED_PLATE_OPENMP:\n" );
  printf ( "  Serial version\n" );
  printf ( "  A program to solve for the steady state temperature distribution\n" );
  printf ( "  over a rectangular plate.\n" );
  printf ( "\n" );
  printf ( "  Spatial grid of %d by %d points.\n", M, N );
  printf ( "  The iteration will be repeated until the change is <= %e\n", err_tol ); 
/*
  Set the boundary values, which don't change. 
*/
  mean = 0.0;
    for ( i = 1; i < M - 1; i++ )
    {
      w[i][0] = 100.0;
    }
    for ( i = 1; i < M - 1; i++ )
    {
      w[i][N-1] =100.0;
    }
    for ( j = 0; j < N; j++ )
    {
      w[M-1][j] = 100.0;
    }
    for ( j = 0; j < N; j++ )
    {
      w[0][j] = 0.0;
    }
/*
  Average the boundary values, to come up with a reasonable
  initial value for the interior.
*/
    for ( i = 1; i < M - 1; i++ )
    {
      mean = mean + w[i][0] + w[i][N-1];
    }
    for ( j = 0; j < N; j++ )
    {
      mean = mean + w[M-1][j] + w[0][j];
    }
  mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
  printf ( "\n" );
  printf ( "  Mean temperature = %f\n", mean );
/* 
  Initialize the interior temperature to the mean value.
*/
    for ( i = 1; i < M - 1; i++ )
    {
      for ( j = 1; j < N - 1; j++ )
      {
        w[i][j] = mean;
      }
    }
/*
  iterate until the  new solution W differs from the old solution U
  by no more than err_tol.
*/
  itr = 0;
  itr_print = 1;
  printf ( "\n" );
  printf ( " Iteration  Change in temperature:\n" );
  printf ( "\n" );
  wtime = omp_get_wtime ( );

  diff = err_tol;
  while ( err_tol <= diff )
  {
/*
  Save the old solution in U.
*/
  update(M,N,w,u);
/*
  Determine the new estimate of the solution at the interior points.
  The new solution W is the average of north, south, east and west neighbours.
*/
   central(M,N,w,u);
  
/*Now take the difference of the previous temperature and the new one for every grid point
  and find the maximum from it*/   
   diff = 0.0;
      for ( i = 1; i < M - 1; i++ )
      {
        for ( j = 1; j < N - 1; j++ )
        {
          if ( diff < fabs ( w[i][j] - u[i][j] ) )
          {
            diff = fabs ( w[i][j] - u[i][j] );
          }
          
        }
      }
      
    itr++;
    
 /*print the difference for specific iterations not for every iteration*/
 
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
  printf ( "\nSteady State has been achieved.\n" );
  printf ( "HEATED_PLATE_OPENMP:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;

}