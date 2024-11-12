
#define Bc(i,j) Bc[(i)+kc*(j)]
#define Ac(i,j) Ac[(i)+mc*(j)]
inline __attribute__((always_inline)) void packB( const int kc, const int nc, const int nr, double *B, const int k, double *Bc ) {
  int ind = 0;
  for( int j = 0; j < nc; j+=nr ) {
    for( int p = 0; p < kc; p++ ) {
      for( int jj = j; jj < j+nr; jj++ ) {
        Bc[ind++] = B(p,jj);
      }
    }
  }
}

inline __attribute__((always_inline)) void packA( const int mc, const int kc, const int mr, double *A, const int m, double *Ac ) {
  int ind = 0;
  for( int i = 0; i < mc; i+=mr ) {
    for( int p = 0; p < kc; p++ ) {
      for( int ii = i; ii < i+mr; ii++ ) {
        Ac[ind++] = A(ii,p);
      }
    }
  }
}


