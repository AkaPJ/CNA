#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>  // Asegúrate de incluir OpenMP

int main( int argc, char *argv[] ) {
  if( argc < 3 ) {
    fprintf(stderr,"Uso: %s piedras capacidad\n",argv[0]);
    return -1;
  }
  
  double valor = 0;
  int capacidad = 0;
  int n, c;
  
  sscanf(argv[1],"%d",&n);
  sscanf(argv[2],"%d",&c);
  printf("Problema de la mochila con %d piedras y capacidad = %d\n", n, c);

  int w[n];
  double p[n];
  int solucion[n];
  
  // Generando pesos y valores aleatoriamente
  #pragma omp parallel for
  for( int i = 0; i < n; i++ ) {
    w[i] = rand() % 10;
    p[i] = (double) rand() / RAND_MAX;
  }

#ifdef CHECK
  printf("Pesos = [");
  for( int i = 0; i < n; i++ ) {
    printf("%5d",w[i]);
  }
  printf("   ]\n");

  printf("Valores = [");
  for( int i = 0; i < n; i++ ) {
    printf("%6.2f",p[i]);
  }
  printf("   ]\n");
#endif

  memset(solucion, 0, n * sizeof(int));
  long max = (long) powl(2, n);
  
  if( max < 0 ) {
    printf("El problema no se puede calcular\n");
    return 0;
  }

  // Variables globales para la mejor solución
  double mejor_valor = 0;
  int mejor_capacidad = c;
  int mejor_solucion[n];
  memset(mejor_solucion, 0, n * sizeof(int));

  // Paralelización adecuada
  #pragma omp parallel
  {
    // Variables locales para cada hilo
    double valor_local = 0;
    int capacidad_local = 0;
    int solucion_local[n];
    int s[n];

    #pragma omp for
    for (long x = 0; x < max; x++) {
      int cap = c;
      double val = 0;
      int k = 0;
      memset(s, 0, n * sizeof(int));
      long y = x;
      
      // Explorar una posible combinación de la mochila
      while (y > 0 && cap >= 0) {
        int resto = y % 2;
        if (resto) {
          cap -= w[k];
          val += p[k];
        }
        s[k] = resto;
        y /= 2;
        k++;
      }
      
      if (cap >= 0) {
        if (val > valor_local) {
          valor_local = val;
          capacidad_local = cap;
          memcpy(solucion_local, s, n * sizeof(int));
        }
      }
    }

    // Sección crítica para actualizar las mejores soluciones
    #pragma omp critical
    {
      if (valor_local > mejor_valor) {
        mejor_valor = valor_local;
        mejor_capacidad = capacidad_local;
        memcpy(mejor_solucion, solucion_local, n * sizeof(int));
      }
    }
  }

#ifdef CHECK
  printf("Solucion = [");
  for (int i = 0; i < n; i++) {
    printf("%2d", mejor_solucion[i]);
  }
  printf(" ]\n");
#endif

  printf("Mejor valor = %f\n", mejor_valor);
  printf("Capacidad restante = %d\n", mejor_capacidad);

  return 0;
}
