#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define MAX_NUMEROS 100000

#define SOLUCION_ENCONTRADA 100
#define BUSCAR_PRIMO 200
#define FIN 300
#define ESCLAVO_PREPARADO 400

int esPrimo(int n);
int main(int argc,char *argv[]) {
  int rank, size, continuar;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status status;
  int cantidad;
  
  
  if(rank == 0){
  	int i = 0;
  	int encontrado = 0;
  	int esclavos = size-1;
  
   	if (argc != 2) {
	    printf("Uso: %s <numero de primos a generar>\n", argv[0]);
	    MPI_Finalize();
	    return 1;
	}
	cantidad = atoi(argv[1]);

  	while (encontrado<cantidad){
  	
	
  		int recibido;
  		MPI_Recv(&recibido,1, MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
  		if(status.MPI_TAG == SOLUCION_ENCONTRADA ){
  			printf("El numero primo encontrado es: %d\n", recibido);
  			//esclavos--;
  			encontrado++;
  		}
  		
		if(encontrado<cantidad){
			int aleatorio = (int) rand() % 10000;
			MPI_Send(&aleatorio,1,MPI_INT,status.MPI_SOURCE,BUSCAR_PRIMO,MPI_COMM_WORLD);//etiqueta buscar primo
			i++;
		}else{
		
			for(int j=1; j<size;j++){
				int fin = 1;
				MPI_Send(&fin,1,MPI_INT,j,FIN,MPI_COMM_WORLD);//etiqueta fin	
                	}
		}
  		
  	
	  		
  	}
  }
  
  if (rank != 0) {
    MPI_Send(&rank,1,MPI_INT,0,ESCLAVO_PREPARADO,MPI_COMM_WORLD);//etiqueta ESCLAVO_PREPARADO
    continuar = 1;
    int numero;
    while(continuar!=0){
    	MPI_Recv(&numero,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    	if(status.MPI_TAG==FIN){
    		continuar = 0;
    	}
    	else{
    		if(esPrimo(numero)){
    			MPI_Send(&numero,1,MPI_INT,0,SOLUCION_ENCONTRADA,MPI_COMM_WORLD);
    		        //continuar=0;
    		}
    		else{
    			MPI_Send(&rank,1,MPI_INT,0,ESCLAVO_PREPARADO,MPI_COMM_WORLD);
    		}
    	}
    }
    
  }

  MPI_Finalize();
  return 0;
}

int esPrimo(int num){
	if(num<1){
		return 0;
	} else{
		for(int i =2; i<num; i++){
			if(num%i == 0) return 0;	
		}
		return 1;
	}
}
