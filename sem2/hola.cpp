#include <iostream>

using namespace std;

int main(int argc, char *argv[]){
  int edad;
  int anyos;
  int accidentes;
  int pago = 300; //300 es la base
  cout << "Introduce tu edad: ";
  cin >> edad;
  cout << "Introduce los años que llevas conduciendo: ";
  cin >> anyos;
  cout << "Introduce el número de accidentes que has tenido: ";
  cin >> accidentes;

  if(anyos < 3) {
    pago +=200;
  } else if (anyos >= 3 && anyos < 5) {
    pago += 100;
  } else if (anyos)
  return 0;
}