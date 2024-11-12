#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ctimer.h"

#define Cseq(i,j) Cseq[(i)+m*(j)]
#define C(i,j) C[(i)+m*(j)]
#define A(i,j) A[(i)+m*(j)]
#define B(i,j) B[(i)+k*(j)]

void print( const char *matriz, const int m, const int n, double *M ) {
  printf("%s = [\n",matriz);
  for( int i = 0; i<m; i++ ) {
    for( int j = 0; j<n; j++ ) {
      printf("%9.4f", M[i+m*j]);
    }
    printf("\n");
  }
  printf("]\n");
}

int main(int argc, char *argv[]){

  if( argc<2 ) {
    printf("Usage: <m> <n> <p> seq?\n");
    return 1;
  }
    
  int m, n, k, mc, kc, nc;
  sscanf(argv[1],"%d",&m);
  if( argc>2 ) {
    sscanf(argv[2],"%d",&n);
  } else { n = k = m; }
  if( argc>3 ) {
    sscanf(argv[3],"%d",&k);
  } else { k = m; }
  int seq = 0;
  if( argc>4 ) {
    sscanf(argv[4],"%d",&seq);
  } 
  
  printf("Pon tamaños de bloque de mc nc kc: ");
  scanf("%d %d %d", &mc, &nc, &kc);

  double *A = (double *) malloc (m*k*sizeof(double));
  double *B = (double *) malloc (k*n*sizeof(double));
  double *Cseq = (double *) malloc (m*n*sizeof(double));
  double *C = (double *) malloc (m*n*sizeof(double));

  //Inicializar matrices 
  for( int i = 0; i<m*k; i++ ) {
    A[i] = (double) rand()/RAND_MAX * 2.0 - 1.0;
  }
  for( int i = 0; i<k*n; i++ ) {
    B[i] = (double) rand()/RAND_MAX * 2.0 - 1.0;
  }
  for( int i = 0; i<k*n; i++ ) {
    C[i] = Cseq[i] = 0;
  }
  // print( "A", m, k, A );
  // print( "B", k, n, B );

  if( seq ) {
    double elapsed, ucpu, scpu;
    ctimer( &elapsed, &ucpu, &scpu );
    for( int j = 0; j < n; j++ ) {
       for( int p = 0; p < k; p++ ) {
          for( int i = 0; i < m; i++ ) {
             Cseq(i,j) = Cseq(i,j) + A(i,p) * B(p,j); 
          }
       }
    }
    ctimer( &elapsed, &ucpu, &scpu );
    // print( "Cseq", m, n, Cseq );
    printf("Tiempo = %.1f segundos de la multiplicación simple\n",elapsed);
  }

  /*************************************************************************/
  /* # Comienzo del código a optimizar #                                   */
  /*************************************************************************/
  double elapsed, ucpu, scpu;
  ctimer( &elapsed, &ucpu, &scpu );
  for(int jc = 0; jc < n; jc += nc) {	
  	for(int pc = 0; pc < k; pc += kc) {
  		for(int ic = 0; ic < m; ic += mc) {
  			for( int j = jc; j < jc + nc; j++ ) {
     				for( int p = pc; p < kc + pc; p++ ) {
        				for( int i = ic; i < ic + mc; i++ ) {
           					C(i, j) = C(i, j) + A(i, p) * B(p, j); 
           				}
           			}
           		}
        	}
     	}
  }
  
  
  ctimer( &elapsed, &ucpu, &scpu );
  // print( "C", m, n, C );
  printf("Tiempo = %.1f segundos de la multiplicación optimizada\n",elapsed);
  /*************************************************************************/
  /* # Fin del código a optimizar #                                        */
  /*************************************************************************/

  if( seq ) {
    /* Comprobación del resultado */
    /* Si la optimización es correcta el error debe ser pequeño o cero */
    double error = 0.0;
    for( int i = 0; i<m*n; i++ ) {
      error += fabs(C[i]-Cseq[i]);
    }
    printf("Error = %.2e\n",error);
  }
   
  free(A);
  free(B);
  free(Cseq);
  free(C);

  return 0;
}
