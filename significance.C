#include <iostream>
#include <TMath.h>
using namespace std;

void significance()
{
	Double_t NLL1, NLL2;

	//nonr + K*892p
//	NLL1 = 264.89977;

	//nonr + K*892p + K*892zero
//	NLL1 = 89.81672947463;

	//nonr + K*892p + K*892zero + Kspi-swave
//	NLL2 = 13.18435471466;

	//nonr + K*892p + K*892zero + Kspi-swave + Kppi-swave
//	NLL1 = -7.608576473516;

	//nonr + K*892p + K*892zero + Kspi-swave + Kppi-swave + a980p
//	NLL2 = -13.79408945366;


	//nonr + K*892p + K*892zero + Kspi-swave + Kppi-swave + rho1450p
//	NLL2 = -25.99852101094;


	NLL1 = -9.299463732623;
	NLL2 = -25.99853869172;


	//calculate significance
	Double_t DeltaLL = 2*fabs(NLL1-NLL2);
	Double_t prob = TMath::Prob(DeltaLL, 2);
//	Double_t nSigma = sqrt(2)*TMath::ErfcInverse(prob);
	Double_t nSigma = RooStats::PValueToSignificance(prob*0.5);
	cout<<" delta 2lnL = " << DeltaLL <<"  ;  " << nSigma <<" =  sigma "<<endl;
}
