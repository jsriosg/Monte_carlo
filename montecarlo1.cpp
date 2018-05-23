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

const int teq=(int)(600*std::pow(L/8.0,2.125));
const int tcorr=(int)(50*std::pow(L/8.0,2.125));
const int Nmuestras=10000;

int main(){
  SpinSystem Ising;
  Crandom Ran64(2605);
  double kT, Beta;
  int mcs,mcss;
  double E, M, sumM, sumM2, sumM4, sumE, sumE2;
  kT=2.2; //para cada T
  Beta=1.0/kT;
  Ising.InicieTodosAbajo(); //inicio sistema
  for(mcss=0;mcss<teq;mcss++){
    for(mcs=0;mcs<L2;mcs++){//1 mcs(monte carlo step)
      Ising.UnPasoDeMetropolis(Beta,Ran64);//llegar al equilibrio
    }
  }
  //tomar muestras
  sumM=sumM2=sumM4=sumE=sumE2=0; //inicio acumuladores en cero
  int cmuestras;
  for(cmuestras=0;cmuestras<Nmuestras;cmuestras++){
    E=Ising.GetE(); M=Ising.GetM();//tomo muestras
    sumM+=M; sumM2+=M*M; sumM4+=M*M*M*M; sumE+=E;sumE2+=E*E;
    //avanzar hasta la siguiente muestra
    for(mcss=0;mcss<tcorr;mcss++){
      for(mcs=0;mcs<L2;mcs++){//1 mcs(monte carlo step)
	Ising.UnPasoDeMetropolis(Beta,Ran64);//llegar al equilibrio
      }
    }
  }
  std::cout<<"E="<<Ising.GetE()<<" , M="<<Ising.GetM()<<endl;
  return 0;
}
