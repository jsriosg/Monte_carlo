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
  if(dE<=0){ //lo volteo
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
const double kB=1.0;

int main(){
  SpinSystem Ising;
  Crandom Ran64(2605);
  double kT, Beta;
  int mcs,mcss;
  double E, M, sumM, sumM2, sumM4, sumE, sumE2;
  double Mprom, Eprom, M2prom, E2prom, M4prom;
  double Cv, Xs, Ubinder;
  
  for(kT=0.8;kT<3;kT+=0.1){; //para cada T
    Beta=1.0/kT;
    
    Ising.InicieTodosAbajo(); //inicio sistema
    
    for(mcss=0;mcss<teq;mcss++){
      for(mcs=0;mcs<L2;mcs++){//1 mcs(monte carlo step)
	Ising.UnPasoDeMetropolis(Beta,Ran64);//llegar al equilibrio
      }
    }
    
    //tomar muestras
    sumM=sumM2=sumM4=sumE=sumE2=0.0; //inicio acumuladores en cero
    for(int cmuestras=0;cmuestras<Nmuestras;cmuestras++){

      E=Ising.GetE();
      M=Ising.GetM();//tomo muestras
      sumM+=M; sumM2+=M*M; sumM4+=M*M*M*M; sumE+=E;sumE2+=E*E;
      
      //avanzar hasta la siguiente muestra
      for(mcss=0;mcss<tcorr;mcss++){
	for(mcs=0;mcs<L2;mcs++){//1 mcs(monte carlo step)
	  Ising.UnPasoDeMetropolis(Beta,Ran64);//llegar al equilibrio
	}
      }
      
    }
    
    Mprom=sumM/Nmuestras; M2prom=sumM2/Nmuestras; M4prom=sumM4/Nmuestras;
    E2prom=sumE2/Nmuestras; Eprom=sumE/Nmuestras;
    
    Cv=kB/(kT*kT)*(E2prom - Eprom*Eprom);
    Xs=Beta*(M2prom - Mprom*Mprom);
    Ubinder=1 - 1.0/3*M4prom/(M2prom*M2prom);
  
    std::cout<<kT<<'\t'<<Eprom<<'\t'<<Mprom<<'\t'<<Cv<<'\t'<<Xs<<'\t'<<Ubinder<<std::endl;
  }
  return 0;
}
