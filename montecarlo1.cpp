#include <iostream>
#include <cmath>
#include "random64.h"

const int L=8;
const int L2=L*L;

class SpinSystem{
private:
  int s[L][L],E,M;
public:
  void InicieTodosAbajo(void);
  void UnPasoDeMetropolis(double Beta, Crandom & Ran64);
  double GetE(void){return E;};
  double GetM(void){return std::fabs(M);};
};

void SpinSystem::InicieTodosAbajo(void){
  for(int i=0 ; i<L ; i++)
    for(int j=0 ; j<L; j++)
      s[i][j]=-1;
  M=-L2;E=-2*L2;
}

void SpinSystem::UnPasoDeMetropolis(double Beta, Crandom & Ran64){
  int n,i,j;
  double dE;
  n= (int) L2*Ran64.r(); i = n%L ; j=n/L ; //Escojo un spin al azar
  dE=2*s[i][j]*(s[i][(j+L-1)%L]+s[i][(j+1)%L]+s[(i+1)%L][j]+s[(i+L-1)%L][j]);  //calculo dE que se produciría si lo volteo;
  if(dE>=0){ //lo volteo
    s[i][j]*=-1;
    E += dE;
    M += 2*s[i][j]; 
  }
  else if(Ran64.r()<std::exp(-Beta*dE)){//lo volteo
    s[i][j]*=-1;
    E += dE;
    M += 2*s[i][j]; 
  }
}

int main(){
  SpinSystem Ising;
  int t;
  Ising.InicieTodosAbajo();
  std::cout<<"E="<<Ising.GetE()<<" , M="<<Ising.GetM()<<endl;
  return 0;
}
