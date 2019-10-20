#include <goofit/fitting/FitManagerMinuit2.h>

#include <goofit/Color.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserCovariance.h>
#include "Minuit2/MnStrategy.h"
//#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PdfBase.h>

#include <CLI/Timer.hpp>

namespace GooFit {

FitManagerMinuit2::FitManagerMinuit2(PdfBase *dat)
    : upar_(*dat)
    , fcn_(upar_) {pdfPointer = dat;}

Minuit2::FunctionMinimum FitManagerMinuit2::fit() {
    auto val = Minuit2::MnPrint::Level();
    Minuit2::MnPrint::SetLevel(verbosity);

    // Setting global call number to 0
    host_callnumber = 0;

    CLI::Timer timer{"The minimization took"};

    Minuit2::MnMigrad migrad{fcn_, upar_};


    // Do the minimization
    if(verbosity > 0)
        std::cout << GooFit::gray << GooFit::bold;

    CLI::Timer avetimer{"Average time per call"};
    Minuit2::FunctionMinimum min = migrad(maxfcn_);

	matCov = migrad.Covariance();
//    Minuit2::MnUserCovariance tmp_Cov = migrad.Covariance();
//    std::cout << tmp_Cov << std::endl;
//    Minuit2::MnUserParameters tmp_Cov = migrad.Parameters();
//    std::cout << tmp_Cov << std::endl;

    // Print nice output
    if(verbosity > 0) {
        std::cout << GooFit::reset << (min.IsValid() ? GooFit::green : GooFit::red);
        std::cout << min << GooFit::reset;
        std::cout << GooFit::magenta << timer.to_string() << std::endl;
        std::cout << (avetimer / min.NFcn()).to_string() << GooFit::reset << std::endl;
    }

    if(min.IsValid() && min.HasCovariance() && !min.IsAboveMaxEdm() && !min.HasReachedCallLimit()) {
        retval_ = FitErrors::Valid;
    } else {
        if(verbosity > 0) {
            std::cout << GooFit::red;
            std::cout << "HesseFailed: " << min.HesseFailed() << std::endl;
            std::cout << "HasCovariance: " << min.HasCovariance() << std::endl;
            std::cout << "HasValidCovariance: " << min.HasValidCovariance() << std::endl;
            std::cout << "HasValidParameters: " << min.HasValidParameters() << std::endl;
            std::cout << "IsAboveMaxEdm: " << min.IsAboveMaxEdm() << std::endl;
            std::cout << "HasReachedCallLimit: " << min.HasReachedCallLimit() << std::endl;
            std::cout << "HasAccurateCovar: " << min.HasAccurateCovar() << std::endl;
            std::cout << "HasPosDefCovar : " << min.HasPosDefCovar() << std::endl;
            std::cout << "HasMadePosDefCovar : " << min.HasMadePosDefCovar() << std::endl;
            std::cout << GooFit::reset;
        }

        retval_ = FitErrors::InValid;
    }

    // Set the parameters in GooFit to the new values
    upar_.SetGooFitParams(min.UserState());

    Minuit2::MnPrint::SetLevel(val);
    return min;
}


// For fit errors calculation by xxh

void FitManagerMinuit2::printCovMat()
{
	std::cout << std::endl << matCov << std::endl;
}

//void FitManagerMinuit2::printParams()
//{
//	std::cout << std::endl << pdfPointer->getParameters() << std::endl;
//}


void FitManagerMinuit2::setRandMinuitValues (const int nSamples){
	rnd.SetSeed(nSamples+388);
	std::vector<double> floatVarVal;
	std::vector<double> floatVarErr;
	//std::vector<Variable> vec_vars = upar_.get_vars();
	std::vector<Variable> vec_vars = pdfPointer->getParameters();
//	fvarIndex.clear();
	for(Variable &var : vec_vars) {
		if (var.IsFixed()) continue;
		int counter = var.getFitterIndex();
//		fvarIndex.push_back(var.getFitterIndex());
		floatVarVal.push_back(var.getValue());
		floatVarErr.push_back(var.getError());
//		std::cout << "check for Index " << counter << ": value = " << var.getValue() << "  ;  error = " << var.getError() << std::endl;
	}
	const int nFPars = floatVarVal.size();
	VectorXd vmean(nFPars);
	for (int i = 0; i < nFPars; i++)
		vmean(i) = floatVarVal[i];

	MatrixXd A(nFPars,nFPars);
	const int n = matCov.Nrow();
	if(nFPars != n) std::cout <<"Error!!!!!!!!" << std::endl;
	for (int ii = 0; ii < n; ii++)
		for (int jj = 0; jj < n; jj++)
			A(ii,jj) = matCov(ii,jj);

	SelfAdjointEigenSolver<MatrixXd> es(A);
	sqrtCov = new MatrixXd(es.operatorSqrt());


	VectorXd vy(nFPars);
	samples.clear();
	for (int ii=0;ii<nSamples;ii++){
		for (int i=0;i<nFPars;i++) 
			vy(i) = rnd.Gaus(0,1);
		vy = vmean + (*sqrtCov) * vy;
//		std::cout<<"Mean and random: "<<std::endl<<vmean<<std::endl<<vy<<std::endl;
		samples.push_back(vy);
	}
}

void FitManagerMinuit2::loadSample (const int iSample){
//	std::vector<Variable> vec_vars = upar_.get_vars();
	std::vector<Variable> vec_vars = pdfPointer->getParameters();
	int counter = 0;
//	for(int i = 0; i < vec_vars.size(); ++i){
	for(Variable &var : vec_vars){
		if (var.IsFixed()) continue;
//		std::cout << "load value Index " << var.getFitterIndex() << " : " << samples[iSample][counter] << std::endl;
		pdfPointer->updateVariable(var,(fptype)samples[iSample][counter]);
		counter++;
	}
}


} // namespace GooFit