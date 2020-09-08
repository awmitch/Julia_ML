# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
//# include <gtk/gtk.h> 
# include <omp.h>
# include <sys/stat.h>
# include <sys/types.h>
# include "julia_gen.h"

#define BYTES_PER_PIXEL 3

/******************************************************************************/
/*
See https://people.sc.fsu.edu/~jburkardt/c_src/mandelbrot_openmp/mandelbrot_openmp.html
for original code to generate a mandelbrot set. Note the below should not be...
taken for grade as it only produces raw Julia set images
*/


int img_gen(int m, int n, int count_max, double scale, double rotation, 
            double c1, double c2, char *output_filename)
{  
  int r_start = 0;
  int r_end = 255;
  int g_start = 0;
  int g_end = 255;
  int b_start = 0;
  int b_end = 255;

  int *hist;
  int stride;
  int stride_adjust;

  int **r;
  int **g;
  int **b;
  double **Count;
  int c;
  int hist_adj;

  int i;
  int j;
  int jhi;
  int jlo;
  int k;
  int tid;

  unsigned char *rgb;
  FILE *output_unit;

  double wtime;
  double x_max =   1.75;
  double x_min = - 1.75;
  double x;
  double x1;
  double x2;
  double y_max =   1.75;
  double y_min = - 1.75;
  double y;
  double y1;
  double y2;

  //printf ( "\n" );
  //printf ( "MANDELBROT_OPENMP\n" );
  //printf ( "  C/OpenMP version\n" );
  //printf ( "\n" );
  //printf ( "  Create an ASCII PPM image of the Mandelbrot set.\n" );
  //printf ( "\n" );
  //printf ( "  For each point C = X + i*Y\n" );
  //printf ( "  with X range [%g,%g]\n", x_min, x_max );
  //printf ( "  and  Y range [%g,%g]\n", y_min, y_max );
  //printf ( "  carry out %d iterations of the map\n", count_max );
  //printf ( "  Z(n+1) = Z(n)^2 + C.\n" );
  //printf ( "  If the iterates stay bounded (norm less than 2)\n" );
  //printf ( "  then C is taken to be a member of the set.\n" );
  //printf ( "\n" );
  //printf ( "  An ASCII PPM image of the set is created using\n" );
  //printf ( "    M = %d pixels in the X direction and\n", m );
  //printf ( "    N = %d pixels in the Y direction.\n", n );
  wtime = omp_get_wtime ( );

  /*
    Carry out the iteration for each pixel, determining Count.
  */
    
  r = ( int ** ) malloc ( m * sizeof ( int * ) );
  g = ( int ** ) malloc ( m * sizeof ( int * ) );
  b = ( int ** ) malloc ( m * sizeof ( int * ) );
  Count = ( double ** ) malloc ( m * sizeof ( double * ) );
  for (i=0; i < m; i++) {
    Count[i] = ( double * ) malloc ( n * sizeof ( double ));
    r[i] = ( int * ) malloc ( n * sizeof ( int ));
    g[i] = ( int * ) malloc ( n * sizeof ( int ));
    b[i] = ( int * ) malloc ( n * sizeof ( int ));
  }

# pragma omp parallel \
  default(none) \
  shared ( Count, count_max, x_max, x_min, y_max, y_min, tid, n, m, r, g, b, scale, rotation, c1, c2) \
  private ( i, j, k, x, x1, x2, y, y1, y2, c)
  {
    int c;
    /*if (omp_get_thread_num () == 0) {
      tid = omp_get_num_threads (); 
      printf ( "  Number of threads is: %d\n",tid);
    }*/
# pragma omp for schedule(static,100)
    for ( i = 0; i < m; i++ ) {
      y = ( ( double ) (     i - 1 ) * y_max   
          + ( double ) ( m - i     ) * y_min ) 
          / ( double ) ( m     - 1 );

      for ( j = 0; j < n; j++ ) {
        x = ( ( double ) (     j - 1 ) * x_max   
            + ( double ) ( n - j     ) * x_min ) 
            / ( double ) ( n     - 1 );

        Count[i][j] = 0.0; 
        x1 = x;
        y1 = y;
        for ( k = 1.0; k <= count_max; k++ ) {
          x2 = x1 * x1 - y1 * y1 + c1;
          y2 = 2 * x1 * y1 + c2;
          if ( x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2 ) {
            Count[i][j] = k;
            break;
          }
          x1 = x2;
          y1 = y2;
        }
        if ( Count[i][j] == count_max ) {
          r[i][j] = 255;
          g[i][j] = 255;
          b[i][j] = 255;
        }
        else {
          c = ( int ) ( 255.0 * sqrt ( sqrt ( //sqrt ( 
              ( ( double ) ( Count[i][j] ) / ( double ) ( count_max ) ) ) ) );// );
          r[i][j] = c;
          g[i][j] = c;
          b[i][j] = c;  
        }
      }
    }
  }
  wtime = omp_get_wtime ( ) - wtime;
  //printf ( "\n" );
  printf ( "Wrote %s\n\t in %g seconds.\n", output_filename, wtime);

  FILE * fp;
  fp = fopen("test", "w+");
  if (fp == NULL){
    perror("Could not open file");
    fclose(fp);
  }
  else {
    for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
        fprintf(fp, "%d ", r[i][j]);
        //printf("it");
      }
      fprintf(fp,"\n");
    }
  }
  fclose(fp); 

  /*
    Write data to an ASCII PPM file.
  */

  output_unit = fopen ( output_filename, "wt" );
  fprintf ( output_unit, "P3\n" );
  fprintf ( output_unit, "%d  %d\n", n, m );
  fprintf ( output_unit, "%d\n", 255 );

  for ( i = 0; i < m; i++ )
  {
    for ( jlo = 0; jlo < n; jlo = jlo + 4 )
    {
      jhi = i4_min ( jlo + 4, n );
      for ( j = jlo; j < jhi; j++ )
      {
        fprintf ( output_unit, "  %d  %d  %d\n", r[i][j], g[i][j], b[i][j] );
      }
      //fprintf ( output_unit, "\n" );
    }
  }
  fclose(output_unit);
  //printf ( "\n" );
  //printf ( "  Graphics data written to \"%s\".\n", output_filename );
   
  for (i=0; i < m; i++) {
    free ( Count[i] );
    free ( r[i] );
    free ( g[i] );
    free ( b[i] );
  }
  free ( Count );
  free ( g );
  free ( r );
  free ( b );

  /*
    Terminate.
  */
  //printf ( "\n" );
  //printf ( "MANDELBROT_OPENMP\n" );
  //printf ( "  Normal end of execution.\n" );
  //printf ( "\n" );
  return 0;
}

/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}


