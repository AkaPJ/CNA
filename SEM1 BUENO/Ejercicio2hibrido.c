#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc,char *argv[]) { 
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  if (argc<2) {
    printf("Uso: %s radio [max_puntos]\n",argv[0]);
    MPI_Finalize();
    return 0;
  }
  double radio;
  sscanf(argv[1],"%lf",&radio);
  double max_puntos = 10E7;
  if( argc > 2 ) {
    sscanf(argv[2],"%lf",&max_puntos);
  }

  if(rank==0){
  	printf("Calculo del volumen de una esfera de radio %.2f (max_puntos = %.0e)\n",radio,max_puntos);
  }

  double total_aciertos = 0, aciertos = 0;
  #pragma omp parallel
  {
  	
  	double aciertos_hilo = 0;
  
  	#pragma omp parallel for
 	for( unsigned long int i = 0; i< (unsigned long int)(max_puntos/(size*2)); i++ ) {
     		double x = (((double) rand() ) / RAND_MAX ) * radio ;
     		double y = (((double) rand() ) / RAND_MAX ) * radio ;
     		double z = (((double) rand() ) / RAND_MAX ) * radio ;
    		if( ( x*x + y*y + z*z ) <= (radio*radio) ) aciertos_hilo+=1.0;
  	}
  	
  	#pragma omp atomic
 	aciertos += aciertos_hilo;
	  
  }
  MPI_Reduce(&aciertos, &total_aciertos, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0){
  	double volumen = radio * radio * radio * total_aciertos / max_puntos * 8;

  	printf("Volumen calculado = %f, volumen esfera = %f. \n", volumen,	 4.0*M_PI*radio*radio*radio/3.0);
  }
  
  MPI_Finalize();
  return 0;
}
