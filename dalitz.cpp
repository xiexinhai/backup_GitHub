//modified to analysis: D+ -> Ks K+ pi0
//xiexh 2019 04 13

// ROOT stuff
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TMath.h>
#include <TCut.h>

// System stuff
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>

// GooFit stuff
#include <goofit/Application.h>
#include <goofit/FitManager.h>
#include <goofit/PDFs/GooPdf.h>
#include <goofit/PDFs/basic/PolynomialPdf.h>
#include <goofit/PDFs/basic/TrigThresholdPdf.h>
#include <goofit/PDFs/combine/AddPdf.h>
#include <goofit/PDFs/combine/ProdPdf.h>
#include <goofit/PDFs/physics/Amp3Body.h>
#include <goofit/PDFs/physics/DalitzPlotPdf.h>
#include <goofit/PDFs/physics/DalitzPlotter.h>
#include <goofit/PDFs/physics/DalitzVetoPdf.h>
#include <goofit/PDFs/physics/ResonancePdf.h>
#include <goofit/PDFs/basic/SmoothHistogramPdf.h>
#include <goofit/UnbinnedDataSet.h>
#include <goofit/BinnedDataSet.h>
#include <goofit/Variable.h>
#include <goofit/utilities/Style.h>

using namespace std;
using namespace GooFit;


const fptype m_K892p = 0.89176;
//const fptype m_K892p = 0.89176+0.00025;
//const fptype m_K892p = 0.89176-0.00025;
const fptype w_K892p = 0.0503;
//const fptype w_K892p = 0.0503+0.0008;
//const fptype w_K892p = 0.0503-0.0008;

const fptype m_K892zero = 0.89555;
//const fptype m_K892zero = 0.89555+0.00020;
//const fptype m_K892zero = 0.89555-0.00020;
const fptype w_K892zero = 0.0473;
//const fptype w_K892zero = 0.0473+0.0005;
//const fptype w_K892zero = 0.0473-0.0005;

const fptype m_K1430 = 1.441;
//const fptype m_K1430 = 1.441+0.002;
//const fptype m_K1430 = 1.441-0.002;

const fptype w_K1430 = 0.193;
//const fptype w_K1430 = 0.193+0.004;
//const fptype w_K1430 = 0.193-0.004;

//old
//const fptype m_K1430 = 1.425;
//const fptype w_K1430 = 0.27;

const fptype m_a980p = 0.980;
//const fptype m_a980p = 0.980+0.020;
//const fptype m_a980p = 0.980-0.020;
const fptype c_g1 = 0.341;
//const fptype c_g1 = 0.341+0.004;
//const fptype c_g1 = 0.341-0.004;
const fptype c_g2og1 = 0.892;
//const fptype c_g2og1 = 0.892+0.022;
//const fptype c_g2og1 = 0.892-0.022;

const fptype _radius = 1.5;

const bool fitMasses = false;
char strbuffer[1000];
const fptype _mDp       = 1.86965;
const fptype _mDp2      = _mDp * _mDp;
const fptype _mDp2inv   = 1. / _mDp2;
const fptype KsMass = 0.497611;
const fptype KpMass = 0.493677;
const fptype pi0Mass = 0.134977;
const int BinNumsM12 = 3500;
const int BinNumsM13 = 3500;
const int drawBinM12 = 1000;
const int drawBinM13 = 1000;

const fptype m12_lower = (KpMass+pi0Mass)*(KpMass+pi0Mass);
const fptype m12_upper = (_mDp-KsMass)*(_mDp-KsMass);

const fptype m13_lower = (KsMass+pi0Mass)*(KsMass+pi0Mass);
const fptype m13_upper = (_mDp-KpMass)*(_mDp-KpMass);

const fptype m23_lower = (KsMass+KpMass)*(KsMass+KpMass);
const fptype m23_upper = (_mDp-pi0Mass)*(_mDp-pi0Mass);

bool m_draw_data = true;
bool m_effPoly = false;

bool m_float_init = false;
bool m_float_polyeff = false;

bool m_draw_polyeff = false;
bool m_draw_smootheff = false;

bool m_test_pdf = false;

bool m_fit_data = true;
bool m_fit_use_eff = true;

bool m_test_for_pull = true;
bool m_print_for_sys = true;

const fptype dE_min = -0.03;
const fptype dE_max = 0.02;
const fptype mrec_min = 1.8648;
const fptype mrec_max = 1.8772;
//const fptype dE_min = -0.03;
//const fptype dE_max = 0.02;
//const fptype mrec_min = 1.8648;
//const fptype mrec_max = 1.8772;


std::string eff_filename="weighted_acceptance_bin150.root";
//std::string eff_filename="weighted_acceptance_dE_plus.root";
//std::string eff_filename="weighted_acceptance_dE_mius.root";
//std::string eff_filename="weighted_acceptance_Mrec_plus.root";
//std::string eff_filename="weighted_acceptance_Mrec_mius.root";


std::string data_filename = "data_tagAll.root";
//std::string data_filename = "tagAll_DIY.root";

TRandom3 rnd;

void makeEffPlot(GooPdf* total, UnbinnedDataSet* data);
void saveToy(Amp3Body *signal, Observable &m12, Observable &m13, EventNumber &eventNumber,int,int);
// Constants used in more than one PDF component.
Variable motherM("motherM", _mDp);
Variable chargeM("chargeM", KpMass);
Variable neutrlMKs("neutrlMKs", KsMass);
Variable neutrlMpi0("neutrlMpi0", pi0Mass);
Variable massSum("massSum", _mDp * _mDp + KsMass * KsMass+ KpMass * KpMass + pi0Mass * pi0Mass); 
Variable constantOne("constantOne", 1);
Variable constantOne2("constantOne2", 1);
Variable constantZero("constantZero", 0);

fptype fRand(double fMin, double fMax)
{
     double f = rnd.Rndm();
     return fMin + f * (fMax - fMin);
}

fptype cpuGetM23(fptype massPZ, fptype massPM) {
    return (_mDp2 + KsMass * KsMass + KpMass * KpMass + pi0Mass * pi0Mass - massPZ - massPM);
}

bool cpuDalitz (fptype m12, fptype m13, fptype bigM = _mDp, fptype dm1 = pi0Mass, fptype dm2 = KpMass, fptype dm3 = KsMass) {
  if (m12 < TMath::Power(dm1 + dm2, 2)) return false; // This m12 cannot exist, it's less than the square of the (1,2) particle mass.
  if (m12 > TMath::Power(bigM - dm3, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter. 
  if (m13 < pow(dm1 + dm3, 2)) return false; // This m13 cannot exist, it's less than the square of the (1,2) particle mass.
  if (m13 > pow(bigM - dm2, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter. 

  // Calculate energies of 1 and 3 particles in m12 rest frame. 
  fptype e1star = 0.5 * (m12 - dm2*dm2 + dm1*dm1) / sqrt(m12);
  fptype e3star = 0.5 * (bigM*bigM - m12 - dm3*dm3) / sqrt(m12);

  // Bounds for m13 at this value of m12.
  fptype minimum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - dm1*dm1) + sqrt(e3star*e3star - dm3*dm3), 2);
  if (m13 < minimum) return false;
  fptype maximum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - dm1*dm1) - sqrt(e3star*e3star - dm3*dm3), 2);
  if (m13 > maximum) return false;
  // Calculate energies of 1 and 2 particles in m13 rest frame. 
  e1star = 0.5 * (m13 - dm3*dm3 + dm1*dm1) / sqrt(m13);
  fptype e2star = 0.5 * (bigM*bigM - m13 - dm2*dm2) / sqrt(m13);

  // Bounds for m12 at this value of m13.
  fptype m12minimum = pow(e1star + e2star, 2) - pow(sqrt(e1star*e1star - dm1*dm1) + sqrt(e2star*e2star - dm2*dm2), 2);
  if (m12 < m12minimum) return false;
  fptype m12maximum = pow(e1star + e2star, 2) - pow(sqrt(e1star*e1star - dm1*dm1) - sqrt(e2star*e2star - dm2*dm2), 2);
  if (m12 > m12maximum) return false;
//  if (m13< maximum - delta && m13> minimum+ delta 
//          && m12 < m12maximum - delta && m12> m12minimum+ delta ) return true;

//  return false; 
  return true;
}

struct BigBin {
	int xbin;
	int ybin;
	int width;
	int height; 
	double getContent (const TH2F* plot); 
};

double BigBin::getContent (const TH2F* plot) {
	double ret = 0;
	//std::cout << "getContent with " << width << " " << height << " " << xbin << " " << ybin <<std::endl;
	for (unsigned int i = 0; i < width; ++i) {
		for (unsigned int j = 0; j < height; ++j) {
			//std::cout << i << ", " << j << std::endl; 
			if (xbin+i > plot->GetNbinsX()) continue;
			if (ybin+j > plot->GetNbinsY()) continue;
			ret += plot->GetBinContent(xbin+i, ybin+j);
		}
	}
	//std::cout << "Total " << ret << std::endl; 
	return ret;
}

struct ChisqInfo {
	ChisqInfo (); 
	double chisq;
	int dof;
	TH2F* contribPlot; 
};

ChisqInfo::ChisqInfo () 
	: chisq(0), dof(0), contribPlot(0)
{}

ChisqInfo* getAdaptiveChisquare (const TH2F* datPlot, const TH2F* pdfPlot, bool useLH = false) {
	bool acceptable = false;
	int binSize = 1; 
	vector<BigBin> binlist;
	double limit = 10;
	while (!acceptable) {
		binlist.clear();
//		std::cout << "Attempting bin generation with size " << binSize << std::endl; 
	
		for (int xbin = 1; xbin <= datPlot->GetNbinsX(); xbin += binSize) {
			for (int ybin = 1; ybin <= datPlot->GetNbinsY(); ybin += binSize) {
				double lox = datPlot->GetXaxis()->GetBinLowEdge(xbin+0);
				double hix = datPlot->GetXaxis()->GetBinLowEdge(xbin+1+binSize);
				double loy = datPlot->GetYaxis()->GetBinLowEdge(ybin+0);
				double hiy = datPlot->GetYaxis()->GetBinLowEdge(ybin+1+binSize);
				bool corner = false; 
				if      (cpuDalitz(lox, loy)&&((!useLH)||(lox < loy))) corner = true;
				else if (cpuDalitz(lox, hiy)&&((!useLH)||(lox < hiy))) corner = true;
				else if (cpuDalitz(hix, loy)&&((!useLH)||(hix < loy))) corner = true;
				else if (cpuDalitz(hix, hiy)&&((!useLH)||(hix < hiy))) corner = true;
				if (!corner) continue; 
			
				BigBin curr; 
				curr.xbin = xbin;
				curr.ybin = ybin;
				curr.width = curr.height = binSize;
				binlist.push_back(curr); 
			}
		}
	    
		acceptable = true;
		for (vector<BigBin>::iterator bin = binlist.begin(); bin != binlist.end(); ++bin) {
			if ((*bin).getContent(datPlot) > limit) continue; 
			acceptable = false;
			binSize *= 2; 
//			std::cout << "Couldn't get good bins, retry.\n"; 
			break;
		}
	}

	std::cout << "Good bins at size " << binSize << ", beginning splits.\n";

	// Now attempt to split bins. 
	int numSplits = 1;
	while (0 < numSplits) {
	numSplits = 0; 
	vector<BigBin> newbins; 
	for (vector<BigBin>::iterator bin = binlist.begin(); bin != binlist.end(); ++bin) {
		if (1 == (*bin).width*(*bin).height) {
			newbins.push_back(*bin); 
			continue;
      	}
	      BigBin lolef;
	      BigBin lorig;
	      BigBin hilef;
	      BigBin hirig;
	      lolef.xbin = (*bin).xbin;
	      hilef.xbin = (*bin).xbin;
	      lorig.xbin = (*bin).xbin + (*bin).width/2;
	      hirig.xbin = (*bin).xbin + (*bin).width/2;
	      lolef.ybin = (*bin).ybin;
	      hilef.ybin = (*bin).ybin + (*bin).height/2;
	      lorig.ybin = (*bin).ybin;
	      hirig.ybin = (*bin).ybin + (*bin).height/2;
	
	      lolef.width =  (*bin).width/2;
	      lorig.width =  (*bin).width/2;
	      hilef.width =  (*bin).width/2;
	      hirig.width =  (*bin).width/2;
	      lolef.height =  (*bin).height/2;
	      lorig.height =  (*bin).height/2;
	      hilef.height =  (*bin).height/2;
	      hirig.height =  (*bin).height/2;
	      
	      int mask = 0;
	      if (limit < lolef.getContent(datPlot)) mask += 1;
	      if (limit < lorig.getContent(datPlot)) mask += 2;
	      if (limit < hilef.getContent(datPlot)) mask += 4;
	      if (limit < hirig.getContent(datPlot)) mask += 8;
	
	      if (mask != 15) {
		newbins.push_back(*bin);
	      }
	      else {
		newbins.push_back(lolef);
		newbins.push_back(lorig);
		newbins.push_back(hilef);
		newbins.push_back(hirig);
		numSplits++; 
	      }
    }
    binlist.clear(); 
    for (vector<BigBin>::iterator i = newbins.begin(); i != newbins.end(); ++i) binlist.push_back(*i); 
    std::cout << "Split " << numSplits << " bins.\n"; 
  }

  ChisqInfo* ret = new ChisqInfo(); 
  ret->dof = binlist.size(); 
  ret->contribPlot = new TH2F("contribPlot", "", 
			      datPlot->GetNbinsX(), 
			      datPlot->GetXaxis()->GetBinLowEdge(1),
			      datPlot->GetXaxis()->GetBinLowEdge(datPlot->GetNbinsX() + 1),
			      datPlot->GetNbinsY(), 
			      datPlot->GetYaxis()->GetBinLowEdge(1),
			      datPlot->GetYaxis()->GetBinLowEdge(datPlot->GetNbinsY() + 1));

  double totalDat = 0;
  double totalPdf = 0; 
  for (vector<BigBin>::iterator bin = binlist.begin(); bin != binlist.end(); ++bin) {
    double dat = (*bin).getContent(datPlot);
    double pdf = (*bin).getContent(pdfPlot);
//    double term = (dat - pdf) / sqrt(dat); 
    double term = (dat - pdf) / sqrt(pdf); 
    ret->chisq += term*term; 
    /*
    std::cout << "Bin (" << (*bin).xbin << ", " << (*bin).ybin << ") " 
	      << (*bin).width << " " << (*bin).height << " : " 
	      << dat << " " << pdf << " " 
	      << term << std::endl; 
    */
    for (int i = 0; i < (*bin).width; ++i) {
      for (int j = 0; j < (*bin).height; ++j) {
	ret->contribPlot->SetBinContent((*bin).xbin + i, (*bin).ybin + j, 0); 
	bool corner = false; 
	double lox = datPlot->GetXaxis()->GetBinLowEdge((*bin).xbin+i);
	double hix = datPlot->GetXaxis()->GetBinLowEdge((*bin).xbin+i+1); 
	double loy = datPlot->GetYaxis()->GetBinLowEdge((*bin).ybin+j);
	double hiy = datPlot->GetYaxis()->GetBinLowEdge((*bin).ybin+j+1);
	if      (cpuDalitz(lox, loy)&&((!useLH)||(lox < loy))) corner = true;
	else if (cpuDalitz(lox, hiy)&&((!useLH)||(lox < hiy))) corner = true;
	else if (cpuDalitz(hix, loy)&&((!useLH)||(hix < loy))) corner = true;
	else if (cpuDalitz(hix, hiy)&&((!useLH)||(hix < hiy))) corner = true;
	if (!corner) continue; 

	ret->contribPlot->SetBinContent((*bin).xbin + i, (*bin).ybin + j, term); 
      }
    }
    totalPdf += pdf;
    totalDat += dat;
  }

  return ret;
}


double calWeight(double xmin, double xmax, double ymin, double ymax, TRandom3 &r){
     double m12, m13;
     int isample = 10000;
     int ipass = 0;
     for(int i = 0; i < isample; ++i){
          m12 = xmin + r.Rndm() * (xmax - xmin);
          m13 = ymin + r.Rndm() * (ymax - ymin);
          if(cpuDalitz(m12,m13))
               ++ipass;
     }
     double ret;
     if (ipass != 0)
          ret = (double) isample / (double) ipass;
     else
          ret = 0;

     return ret;
}

GooPdf* getDPVeto(Observable m12, Observable m13) {
	Variable motherM("motherM", _mDp);
	Variable dau1M("dau1M", pi0Mass);
	Variable dau2M("dau2M", KpMass);
	Variable dau3M("dau3M", KsMass);

	GooPdf* dp_veto = new DalitzVetoPdf("dp_veto", m12, m13, motherM, dau1M, dau2M, dau3M, vector<VetoInfo>());
	return dp_veto;
}

fptype getM23LowerLimit(Observable m12, Observable m13){
	fptype ret = 999;
	fptype tmp = 999;
	for (int i = 0; i < m12.getNumBins(); ++i){
		fptype tmp_m12 = m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit())*(i) / (fptype)m12.getNumBins();
		for (int j = 0; j <= m13.getNumBins(); ++j){
			fptype tmp_m13 = m13.getLowerLimit() + (m13.getUpperLimit() - m13.getLowerLimit())*(j) / (fptype)m13.getNumBins();
			if(cpuDalitz(tmp_m12,tmp_m13))
				tmp = cpuGetM23(tmp_m12,tmp_m13);
			if(tmp<ret)
				ret = tmp;
		}
	}
	return (ret-0.1);
}

fptype getM23UpperLimit(Observable m12, Observable m13){
	fptype ret = -999;
	fptype tmp = -999;
	for (int i = 0; i <= m12.getNumBins(); ++i){
		fptype tmp_m12 = m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit())*(i) / (fptype)m12.getNumBins();
		for (int j = 0; j <= m13.getNumBins(); ++j){
			fptype tmp_m13 = m13.getLowerLimit() + (m13.getUpperLimit() - m13.getLowerLimit())*(j) / (fptype)m13.getNumBins();
			if(cpuDalitz(tmp_m12,tmp_m13))
				tmp = cpuGetM23(tmp_m12,tmp_m13);
			if(tmp>ret)
				ret = tmp;
		}
	}
	return (ret+0.1);
}

void getRealData(std::string toyFileName, GooFit::Application &app, DataSet &data) {
	toyFileName = app.get_filename(toyFileName, "examples/dalitz");
	
	auto obs               = data.getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

	//here we define: m12 = m^2(K+ pi0), m13 = m^2(Ks pi0)

	
	TFile *f = TFile::Open(toyFileName.c_str());
	TTree *t = (TTree *)f->Get("DTag");
	assert(t);
	std::cout<<"Entries: "<<t->GetEntries()<<std::endl;

	fptype m_Kppi0_sq, m_Kspi0_sq, m_recoil, dE_sig;
	t->SetBranchAddress("m_Kppi0_sq",&m_Kppi0_sq);
	t->SetBranchAddress("m_Kspi0_sq",&m_Kspi0_sq);
	t->SetBranchAddress("m_recoil",&m_recoil);
	t->SetBranchAddress("dE_sig",&dE_sig);

	for(int i = 0; i < t->GetEntries(); i++){
		t->GetEvent(i);
		if(m_recoil<mrec_min||m_recoil>mrec_max) continue;
		if(dE_sig<dE_min||dE_sig>dE_max) continue;

		if(!cpuDalitz(m_Kppi0_sq,m_Kspi0_sq)) {cout << "mKpi0 = " << m_Kppi0_sq << endl << "mKspi0 = " << m_Kspi0_sq << endl;}

		if(!cpuDalitz(m_Kppi0_sq,m_Kspi0_sq)) continue;
		m12.setValue(m_Kppi0_sq);
		m13.setValue(m_Kspi0_sq);
		eventNumber.setValue(data.getNumEvents());
		data.addEvent();
	}
	f->Close();

	GOOFIT_INFO("Read in {} events", data.getNumEvents());
}

void makeToyData(DalitzPlotter &dplotter, UnbinnedDataSet &data) {}

Amp3Body *makeSignalPdf(Observable m12, Observable m13, EventNumber eventNumber, vector <fptype> init_val, int i, GooPdf *eff = 0, bool fixAmps = true) {
    DecayInfo3 dtop0pp;
    dtop0pp.motherMass   = _mDp;
    dtop0pp.daug1Mass    = pi0Mass;
    dtop0pp.daug2Mass    = KpMass;
    dtop0pp.daug3Mass    = KsMass;
    dtop0pp.meson_radius = _radius;

//	bool fixAmps = false; // Takes ~400x longer

	//K*(892)+
	Variable K892pMass("K892p_mass", m_K892p, 0.00025, 0.85, 0.93);
	Variable K892pWidth("K892p_width", w_K892p, 0.0008, 0.04, 0.06);
	ResonancePdf *K892p = new Resonances::RBW(
		"K892p", Variable("K892p_amp_real", 1), Variable("K892p_amp_imag", 0), K892pMass, K892pWidth, 1, PAIR_12);


	//Kbar*(892)0
	Variable K892zeroMass("K892zero_mass", m_K892zero, 0.00020, 0.85, 0.93);
	Variable K892zeroWidth("K892zero_width", w_K892zero, 0.0005, 0.04, 0.06);

	sprintf(strbuffer, "K892zero_amp_real_%d", i);
	Variable K892zero_amp_real(strbuffer,init_val[0], 0.01, 0, 0);
	sprintf(strbuffer, "K892zero_amp_imag_%d", i);
	Variable K892zero_amp_imag(strbuffer,init_val[1], 0.01, 0, 0);
	ResonancePdf *K892zero = new Resonances::RBW(
		"K892zero",
		fixAmps ? Variable("K892zero_amp_real", -0.1567496650418) : K892zero_amp_real,
		fixAmps ? Variable("K892zero_amp_imag", -0.6644198096373) : K892zero_amp_imag,
		K892zeroMass,
		K892zeroWidth,
		1,
		PAIR_13);

	//K*(1410)+
	Variable K1410pMass("K1410p_mass", 1.414, 0.009, 1.2, 1.6);
	Variable K1410pWidth("K1410p_width", 0.232, 0.018, 0.1, 0.4);

	sprintf(strbuffer, "K1410p_amp_real_%d", i);
	Variable K1410p_amp_real(strbuffer,init_val[2], 0.01, 0, 0);
	sprintf(strbuffer, "K1410p_amp_imag_%d", i);
	Variable K1410p_amp_imag(strbuffer,init_val[3], 0.01, 0, 0);
	ResonancePdf *K1410p = new Resonances::RBW(
		"K1410p",
		fixAmps ? Variable("K1410p_amp_real", 0) : K1410p_amp_real,
		fixAmps ? Variable("K1410p_amp_imag", 0) : K1410p_amp_imag,
		K1410pMass,
		K1410pWidth,
		1,
		PAIR_12);

	//Kbar*(1410)0
	Variable K1410zeroMass("K1410zero_mass", 1.414, 0.009, 1.2, 1.6);
	Variable K1410zeroWidth("K1410zero_width", 0.232, 0.018, 0.1, 0.4);

	sprintf(strbuffer, "K1410zero_amp_real_%d", i);
	Variable K1410zero_amp_real(strbuffer,init_val[4], 0.01, 0, 0);
	sprintf(strbuffer, "K1410zero_amp_imag_%d", i);
	Variable K1410zero_amp_imag(strbuffer,init_val[5], 0.01, 0, 0);
	ResonancePdf *K1410zero = new Resonances::RBW(
		"K1410zero",
		fixAmps ? Variable("K1410zero_amp_real", 1.24683082311) : K1410zero_amp_real,
		fixAmps ? Variable("K1410zero_amp_imag", 0.4483709871808) : K1410zero_amp_imag,
		K1410zeroMass,
		K1410zeroWidth,
		1,
		PAIR_13);

	//a980+
	Variable a980pMass("a980pMass",m_a980p,0.020,0.9,1.06);

	sprintf(strbuffer, "a980p_amp_real_%d", i);
	Variable a980p_amp_real(strbuffer,init_val[6], 0.01, 0, 0);
	sprintf(strbuffer, "a980p_amp_imag_%d", i);
	Variable a980p_amp_imag(strbuffer,init_val[7], 0.01, 0, 0);
	ResonancePdf *a980p = new Resonances::FLATTE(
		"a980p",
		fixAmps ? Variable("a980p_amp_real", 80) : a980p_amp_real,
		fixAmps ? Variable("a980p_amp_imag", -25) : a980p_amp_imag,
		a980pMass,
		Variable("g1",c_g1),
		Variable("rg2og1",c_g2og1),
		PAIR_23,
		false);

	//SwaveKppi0
//	Variable SwaveKppi0Mass("SwaveKppi0_Mass", 1.425, 0.050, 1.2, 1.6);
//	Variable SwaveKppi0Width("SwaveKppi0_Width", 0.27, 0.08, 0.1, 0.3);
	Variable SwaveKppi0Mass("SwaveKppi0_Mass", m_K1430, 0.050, 1.2, 1.6);
	Variable SwaveKppi0Width("SwaveKppi0_Width", w_K1430, 0.08, 0.1, 0.3);

	sprintf(strbuffer, "SwaveKppi0_amp_real_%d", i);
	Variable SwaveKppi0_amp_real(strbuffer,init_val[8], 0.01, 0, 0);
	sprintf(strbuffer, "SwaveKppi0_amp_imag_%d", i);
	Variable SwaveKppi0_amp_imag(strbuffer,init_val[9], 0.01, 0, 0);
	ResonancePdf *SwaveKppi0 = new Resonances::LASS(
		"SwaveKppi0",
		fixAmps ? Variable("SwaveKppi0_amp_real", 1.869715774108) : SwaveKppi0_amp_real,
		fixAmps ? Variable("SwaveKppi0_amp_imag", 0.2625631972586) : SwaveKppi0_amp_imag,
		SwaveKppi0Mass,
		SwaveKppi0Width,
		0,
		PAIR_12);

	//SwaveKspi0
//	Variable SwaveKspi0Mass("SwaveKspi0_Mass", 1.425, 0.050, 1.2, 1.6);
//	Variable SwaveKspi0Width("SwaveKspi0_Width", 0.27, 0.08, 0.1, 0.3);
	Variable SwaveKspi0Mass("SwaveKspi0_Mass", m_K1430, 0.050, 1.2, 1.6);
	Variable SwaveKspi0Width("SwaveKspi0_Width", w_K1430, 0.08, 0.1, 0.3);

	sprintf(strbuffer, "SwaveKspi0_amp_real_%d", i);
	Variable SwaveKspi0_amp_real(strbuffer,init_val[10], 0.01, 0, 0);
	sprintf(strbuffer, "SwaveKspi0_amp_imag_%d", i);
	Variable SwaveKspi0_amp_imag(strbuffer,init_val[11], 0.01, 0, 0);
	ResonancePdf *SwaveKspi0 = new Resonances::LASS(
		"SwaveKspi0",
		fixAmps ? Variable("SwaveKspi0_amp_real", -2.61217416547) : SwaveKspi0_amp_real,
		fixAmps ? Variable("SwaveKspi0_amp_imag", -0.126650709178) : SwaveKspi0_amp_imag,
		SwaveKspi0Mass,
		SwaveKspi0Width,
		0,
		PAIR_13);

	//non resonance
	sprintf(strbuffer, "nonr_amp_real_%d", i);
	Variable nonr_amp_real(strbuffer,init_val[12], 0.01, 0, 0);
	sprintf(strbuffer, "nonr_amp_imag_%d", i);
	Variable nonr_amp_imag(strbuffer,init_val[13], 0.01, 0, 0);
	ResonancePdf *nonr = new Resonances::NonRes(
		"nonr",
		fixAmps ? Variable("nonr_amp_real", 16.82603957812) : nonr_amp_real,
		fixAmps ? Variable("nonr_amp_imag", 11.59309345464) : nonr_amp_imag);

	//K1430p
	Variable K1430pMass("K1430p_Mass", 1.425, 0.050, 1.2, 1.6);
	Variable K1430pWidth("K1430p_Width", 0.27, 0.08, 0.1, 0.3);

	sprintf(strbuffer, "K1430p_amp_real_%d", i);
	Variable K1430p_amp_real(strbuffer,init_val[14], 0.01, 0, 0);
	sprintf(strbuffer, "K1430p_amp_imag_%d", i);
	Variable K1430p_amp_imag(strbuffer,init_val[15], 0.01, 0, 0);
	ResonancePdf *K1430p = new Resonances::RBW(
		"K1430p",
		fixAmps ? Variable("K1430p_amp_real", 1.869715774108) : K1430p_amp_real,
		fixAmps ? Variable("K1430p_amp_imag", 0.2625631972586) : K1430p_amp_imag,
		K1430pMass,
		K1430pWidth,
		0,
		PAIR_12);

	//Kbar1430zero
	Variable Kbar1430zeroMass("Kbar1430zero_Mass", 1.425, 0.050, 1.2, 1.6);
	Variable Kbar1430zeroWidth("Kbar1430zero_Width", 0.27, 0.08, 0.1, 0.3);

	sprintf(strbuffer, "Kbar1430zero_amp_real_%d", i);
	Variable Kbar1430zero_amp_real(strbuffer,init_val[16], 0.01, 0, 0);
	sprintf(strbuffer, "Kbar1430zero_amp_imag_%d", i);
	Variable Kbar1430zero_amp_imag(strbuffer,init_val[17], 0.01, 0, 0);
	ResonancePdf *Kbar1430zero = new Resonances::RBW(
		"Kbar1430zero",
		fixAmps ? Variable("Kbar1430zero_amp_real", -2.61217416547) : Kbar1430zero_amp_real,
		fixAmps ? Variable("Kbar1430zero_amp_imag", -0.126650709178) : Kbar1430zero_amp_imag,
		Kbar1430zeroMass,
		Kbar1430zeroWidth,
		0,
		PAIR_13);

	//a0(1450)+
	Variable a1450pMass("a1450p_Mass", 1.474, 0.019, 1.2, 1.6);
	Variable a1450pWidth("a1450p_Width", 0.265, 0.013, 0.1, 0.3);

	sprintf(strbuffer, "a1450p_amp_real_%d", i);
	Variable a1450p_amp_real(strbuffer,init_val[18], 0.01, 0, 0);
	sprintf(strbuffer, "a1450p_amp_imag_%d", i);
	Variable a1450p_amp_imag(strbuffer,init_val[19], 0.01, 0, 0);
	ResonancePdf *a1450p = new Resonances::RBW(
		"a1450p",
		fixAmps ? Variable("a1450p_amp_real", -1.040242409497) : a1450p_amp_real,
		fixAmps ? Variable("a1450p_amp_imag", 1.452729026326) : a1450p_amp_imag,
		a1450pMass,
		a1450pWidth,
		0,
		PAIR_23);


	//Kbar*(1430)0_2
	Variable K1430zero2Mass("K1430zero2_mass", 1.4324, 0.0013, 1.2, 1.6);
	Variable K1430zero2Width("K1430zero2_width", 0.109, 0.005, 0.1, 0.4);
	sprintf(strbuffer, "K1430zero2_amp_real_%d", i);
	Variable K1430zero2_amp_real(strbuffer,init_val[20], 0.01, 0, 0);
	sprintf(strbuffer, "K1430zero2_amp_imag_%d", i);
	Variable K1430zero2_amp_imag(strbuffer,init_val[21], 0.01, 0, 0);
	ResonancePdf *K1430zero2 = new Resonances::RBW(
		"K1430zero2",
		fixAmps ? Variable("K1430zero2_amp_real", 1.714) : K1430zero2_amp_real,
		fixAmps ? Variable("K1430zero2_amp_imag", -0.125) : K1430zero2_amp_imag,
		K1430zero2Mass,
		K1430zero2Width,
		2,
		PAIR_13);

	//rho(1450)+
	Variable rho1450pMass("rho1450p_Mass", 1.465, 0.025, 1.3, 1.6);
	Variable rho1450pWidth("rho1450p_Width", 0.4, 0.06, 0.04, 0.06);

	sprintf(strbuffer, "rho1450p_amp_real_%d", i);
	Variable rho1450p_amp_real(strbuffer,init_val[22], 0.01, 0, 0);
	sprintf(strbuffer, "rho1450p_amp_imag_%d", i);
	Variable rho1450p_amp_imag(strbuffer,init_val[23], 0.01, 0, 0);
	ResonancePdf *rho1450p = new Resonances::GS(
		"rho1450p",
		fixAmps ? Variable("rho1450p_amp_real", -1.040242409497) : rho1450p_amp_real,
		fixAmps ? Variable("rho1450p_amp_imag", 1.452729026326) : rho1450p_amp_imag,
		rho1450pMass,
		rho1450pWidth,
		1,
		PAIR_23);


	//rho(1700)+
	Variable rho1700pMass("rho1700p_Mass", 1.541, 0.035, 1.4, 1.7);
	Variable rho1700pWidth("rho1700p_Width", 0.25, 0.1, 0.1, 0.4);

	sprintf(strbuffer, "rho1700p_amp_real_%d", i);
	Variable rho1700p_amp_real(strbuffer,init_val[24], 0.01, 0, 0);
	sprintf(strbuffer, "rho1450p_amp_imag_%d", i);
	Variable rho1700p_amp_imag(strbuffer,init_val[25], 0.01, 0, 0);
	ResonancePdf *rho1700p = new Resonances::GS(
		"rho1700p",
		fixAmps ? Variable("rho1700p_amp_real", -1.241) :rho1700p_amp_real,
		fixAmps ? Variable("rho1700p_amp_imag", 6.049) : rho1700p_amp_imag,
		rho1700pMass,
		rho1700pWidth,
		1,
		PAIR_23);


	dtop0pp.resonances.push_back(K892p);
	dtop0pp.resonances.push_back(K892zero);

//	dtop0pp.resonances.push_back(K1410p);
//	dtop0pp.resonances.push_back(K1410zero);
//	dtop0pp.resonances.push_back(a980p);

	dtop0pp.resonances.push_back(SwaveKppi0);
	dtop0pp.resonances.push_back(SwaveKspi0);

//	dtop0pp.resonances.push_back(K1430p);
//	dtop0pp.resonances.push_back(Kbar1430zero);
//	dtop0pp.resonances.push_back(nonr);

//	dtop0pp.resonances.push_back(K1430zero2);
//	dtop0pp.resonances.push_back(a1450p);
//	dtop0pp.resonances.push_back(rho1450p);
//	dtop0pp.resonances.push_back(rho1700p);


    if(!fitMasses) {
        for(vector<ResonancePdf *>::iterator res = dtop0pp.resonances.begin(); res != dtop0pp.resonances.end(); ++res) {
            (*res)->setParameterConstantness(true);
        }
    }

    if(!eff) {
        // By default create a constant efficiency.
        vector<Variable> offsets;
        vector<Observable> observables;
        vector<Variable> coefficients;

        observables.push_back(m12);
        observables.push_back(m13);
        offsets.push_back(constantZero);
        offsets.push_back(constantZero);
        coefficients.push_back(constantOne);
        eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
    }
	sprintf(strbuffer, "signalPDF_%d", i);
	return new Amp3Body(strbuffer, m12, m13, eventNumber, dtop0pp, eff);
}


GooPdf* makeHistogramPdf(UnbinnedDataSet &data, string filename = "", string histname = "", string pdfname = ""){
 	cout << "begin generate efficiency pdf" << endl;

	auto obs               = data.getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

	TFile *f = TFile::Open(filename.c_str()); //assert(f);
	TH2D* h2 = (TH2D*)f->Get(histname.c_str()); //assert(h2);
	
	int oldBins1 = m12.getNumBins();
	int oldBins2 = m13.getNumBins();
	m12.setNumBins(h2->GetNbinsX());
	m13.setNumBins(h2->GetNbinsY());

	BinnedDataSet *binEffData = new BinnedDataSet({m12, m13, eventNumber});

	int num = m12.getNumBins()*m13.getNumBins();

//	int tmp = 0;
	for (int i = 0;i < num; i++){
		if(i==Int_t(num*0.1)) cout<<"*******************************completed 10%"<<"**************************************"<<endl;
		if(i==Int_t(num*0.2)) cout<<"*******************************completed 20%"<<"**************************************"<<endl;
		if(i==Int_t(num*0.3)) cout<<"*******************************completed 30%"<<"**************************************"<<endl;
		if(i==Int_t(num*0.4)) cout<<"*******************************completed 40%"<<"**************************************"<<endl;
		if(i==Int_t(num*0.5)) cout<<"*******************************completed 50%"<<"**************************************"<<endl;
		if(i==Int_t(num*0.6)) cout<<"*******************************completed 60%"<<"**************************************"<<endl;
		if(i==Int_t(num*0.7)) cout<<"*******************************completed 70%"<<"**************************************"<<endl;
		if(i==Int_t(num*0.8)) cout<<"*******************************completed 80%"<<"**************************************"<<endl;
		if(i==Int_t(num*0.9)) cout<<"*******************************completed 90%"<<"**************************************"<<endl;
		if(i==Int_t(num*1-1)) cout<<"*******************************completed !!!"<<"**************************************"<<endl;
		int iy = i/m12.getNumBins();
		int ix = i%m12.getNumBins();
//		if (ix > iy){
//			int itt = iy; iy = ix; ix = itt;
//		}
		double content = h2->GetBinContent(ix+1,iy+1);
		binEffData->setBinContent(i, content);

//		if(h2->GetBinContent(i+1) != 0){
//			cout<<"Hist: "<<i<<"     "<<h2->GetBinContent(i+1)<<endl;
//			++tmp;
//		}
//		cout<<"Hist: "<<ix<<"     "<<iy<<"     "<<h2->GetBinContent(ix+1, iy+1)<<endl;

	}
//	cout << "total non-zero BinContent:" << tmp << endl;


//	Variable* effSmoothing = new Variable("effSmoothing", 0); 
	cout << "begin generate SmoothHistogramPdf" << endl;

//	SmoothHistogramPdf* ret = new SmoothHistogramPdf(pdfname.c_str(), binEffData, constantZero); 
	SmoothHistogramPdf* ret = new SmoothHistogramPdf(pdfname.c_str(), binEffData, Variable("smoothConst", 10));

	cout << "end generate SmoothHistogramPdf" << endl;

	m12.setNumBins(oldBins1);
	m13.setNumBins(oldBins2);

	f->Close();
	if(m_draw_smootheff)
		makeEffPlot(ret,&data);
 	cout << "end generate efficiency pdf" << endl;

	return ret; 
}

GooPdf* makeEfficiencyPoly (Observable m12, Observable m13, vector <fptype> init_val) {
	vector<Variable> offsets;
	offsets.clear();
	offsets.push_back(constantOne);
	offsets.push_back(constantOne2); 

	vector<Observable> observables;
	observables.clear();
	observables.push_back(m12);
	observables.push_back(m13); 
	
	/*
	coefficients.push_back(new Variable("x0y0",  1.0)); 
	Variable* x1y0 = new Variable("x1y0",  6.09353e-02, 0.001, -0.1, 0.5);
	coefficients.push_back(x1y0);
	Variable* x2y0 = new Variable("x2y0", -6.46116e-02, 0.001, -0.3, 0.5);
	coefficients.push_back(x2y0);
	Variable* x3y0 = new Variable("x3y0",  -4.66692e-03, 0.001, -1.5, 0.5);
	coefficients.push_back(x3y0);
	coefficients.push_back(x1y0);
	Variable* x1y1 = new Variable("x1y1", -1.58147e-01, 0.001, -0.5, 0.5);
	coefficients.push_back(x1y1);
	Variable* x2y1 = new Variable("x2y1",  8.99077e-02, 0.001, -1.5, 0.5);
	coefficients.push_back(x2y1);
	coefficients.push_back(x2y0); 
	coefficients.push_back(x2y1);
	coefficients.push_back(x3y0);
	*/

	/*
	coefficients.push_back(new Variable("x1y0",  0.07999, 0.01, -0.5, 0.5));
	coefficients.push_back(new Variable("x2y0", -0.23732, 0.01, -0.5, 0.5));
	coefficients.push_back(new Variable("x3y0",  0.10369, 0.01, -1.5, 0.5));
	coefficients.push_back(new Variable("x0y1",  0.10248, 0.01, -1.5, 0.5));
	coefficients.push_back(new Variable("x1y1", -0.28543, 0.01, -0.5, 0.5));
	coefficients.push_back(new Variable("x2y1",  0.15058, 0.01, -1.5, 0.5));
	coefficients.push_back(new Variable("x0y2", -0.20648, 0.01, -0.5, 0.5)); 
	coefficients.push_back(new Variable("x1y2",  0.14567, 0.01, -1.5, 0.5));
	coefficients.push_back(new Variable("x0y3",  0.06231, 0.01, -0.5, 0.5));
	*/
	vector<Variable> coefficients; 
	coefficients.clear();

	Variable x0y0("x0y0",1);
	coefficients.push_back(x0y0);
	Variable x1y0("x1y0", init_val[0], 0.01, -1, 1);
	coefficients.push_back(x1y0);
	Variable x2y0("x2y0", init_val[1], 0.01, -1, 1);
	coefficients.push_back(x2y0);
	Variable x3y0("x3y0", init_val[2], 0.01, -1, 1);
	coefficients.push_back(x3y0);
	Variable x0y1("x0y1", init_val[3], 0.01, -1, 1);
	coefficients.push_back(x0y1);
	Variable x1y1("x1y1", init_val[4], 0.01, -1, 1);
	coefficients.push_back(x1y1);
	Variable x2y1("x2y1", init_val[5], 0.01, -1, 1);
	coefficients.push_back(x2y1);
	Variable x0y2("x0y2", init_val[6], 0.01, -1, 1);
	coefficients.push_back(x0y2);
	Variable x1y2("x1y2", init_val[7], 0.01, -1, 1);
	coefficients.push_back(x1y2);
	Variable x0y3("x0y3", init_val[8], 0.01, -1, 1);
	coefficients.push_back(x0y3);

	PolynomialPdf* poly = new PolynomialPdf("efficiency", observables, coefficients, offsets, 3);
//	poly->setParameterConstantness(true); 

	Variable decXmax("decXmax", init_val[9], 0.001, 0, 5);
	Variable conXmax("conXmax", init_val[10], 0.001, 0, 1);

	Variable decYmax("decYmax", init_val[11], 0.001, 0, 5);
	Variable conYmax("conYmax", init_val[12], 0.001, 0, 1);

	Variable decZmax("decZmax", init_val[13], 0.001, 0, 5);
	Variable conZmax("conZmax", init_val[14], 0.001, 0, 1);
	
	Variable maxDalitzX("maxDalitzX", pow(_mDp - KsMass, 2));
	TrigThresholdPdf* hiX = new TrigThresholdPdf("hiX", m12, maxDalitzX, decXmax, conXmax, true); 
	
	Variable maxDalitzY("maxDalitzY", pow(_mDp - KpMass, 2));
	TrigThresholdPdf* hiY = new TrigThresholdPdf("hiY", m13, maxDalitzY, decYmax, conYmax, true); 

	Variable maxDalitzZ("maxDalitzZ", pow(_mDp - pi0Mass, 2));
	TrigThresholdPdf* hiZ = new TrigThresholdPdf("hiZ", m12, m13, maxDalitzZ, decZmax, conZmax, massSum, true); 

	Variable decXmin("decXmin", init_val[15], 0.001, 0, 50);
	Variable conXmin("conXmin", init_val[16], 0.001, 0, 1);

	Variable decYmin("decYmin", init_val[17], 0.001, 0, 50);
	Variable conYmin("conYmin", init_val[18], 0.001, 0, 1);

     Variable decZmin("decZmin", init_val[19], 0.001, 0, 50);
     Variable conZmin("conZmin", init_val[20], 0.001, 0, 1);


	Variable minDalitzX("minDalitzX", pow(KpMass + pi0Mass, 2));
	TrigThresholdPdf* loX = new TrigThresholdPdf("loX", m12, minDalitzX, decXmin, conXmin, false);

	Variable minDalitzY("minDalitzY", pow(KsMass + pi0Mass, 2));
	TrigThresholdPdf* loY = new TrigThresholdPdf("loY", m13, minDalitzY, decYmin, conYmin, false);

	Variable minDalitzZ("minDalitzZ", pow(KsMass + KpMass, 2));
	TrigThresholdPdf* loZ = new TrigThresholdPdf("loZ", m12, m13, minDalitzZ, decZmin, conZmin, massSum, false);

	vector<PdfBase*> comps;
	comps.clear();
	comps.push_back(poly);
	comps.push_back(hiX);
	comps.push_back(hiY);
	comps.push_back(hiZ);

	comps.push_back(loX);
	comps.push_back(loY);
	comps.push_back(loZ);
	comps.push_back(getDPVeto(m12,m13)); 
	ProdPdf* ret = new ProdPdf("efficiency_total", comps); 
	return ret; 
}

UnbinnedDataSet *loadEffData(UnbinnedDataSet &data,std::string toyFileName) {
	auto obs               = data.getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

//	EventNumber eventNumber("eventNumber");
//	BinnedDataSet *effdata = new BinnedDataSet({m12,m13,eventNumber});
	UnbinnedDataSet *effdata = new UnbinnedDataSet({m12,m13,eventNumber});
	std::cout<<"Reading file "<<toyFileName<<std::endl;

	TChain *t = new TChain("DTag");
	t->Add(toyFileName.c_str());
	std::cout<<"Entries: "<<t->GetEntries()<<std::endl;
	
	double m2_12, m2_13;
	t->SetBranchAddress("m_Kppi0_sq", &m2_12);
	t->SetBranchAddress("m_Kspi0_sq", &m2_13);
	for (int i=0;i<t->GetEntries();i++){
	    t->GetEntry(i);
	    if (!cpuDalitz(m2_12, m2_13)) continue;
	    m12.setValue(m2_12);
	    m13.setValue(m2_13);
	    eventNumber.setValue(effdata->getNumEvents()); 
	    effdata->addEvent(); 
	}
	delete t;
	std::cout<<"Passed Entries: "<<effdata->getNumEvents()<<std::endl;
	return effdata; 
}

void makeEffPlot(GooPdf* total, UnbinnedDataSet* data){
	auto obs               = data->getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

	int oldBins1 = m12.getNumBins();
	int oldBins2 = m13.getNumBins();
	m12.setNumBins(drawBinM12);
	m13.setNumBins(drawBinM13);


//	ProdPdf prodpdf{"effpdf", {total}};
//	DalitzPlotter plotter(&prodpdf, (Amp3Body *)total);
//	UnbinnedDataSet toyMC({m12, m13, eventNumber});
//	plotter.fillDataSetMC(toyMC, 5000000);

	//prepare drawing hist
	TH1::SetDefaultSumw2(1);
	TH2F *effpdfplot = new TH2F("effpdfplot",
                    "",
                    m12.getNumBins(),
                    m12.getLowerLimit(),
                    m12.getUpperLimit(),
                    m13.getNumBins(),
                    m13.getLowerLimit(),
                    m13.getUpperLimit());
	effpdfplot->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
	effpdfplot->GetYaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");

	TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	m12_pdf_hist.SetStats(false); 
	m12_pdf_hist.SetLineColor(kRed); 
	m12_pdf_hist.SetLineWidth(2); 

	TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m13_pdf_hist.SetStats(false); 
	m13_pdf_hist.SetLineColor(kRed); 
	m13_pdf_hist.SetLineWidth(2); 

	TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
//	TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), m23_lower, m23_upper);
	m23_pdf_hist.SetStats(false); 
	m23_pdf_hist.SetLineColor(kRed); 
	m23_pdf_hist.SetLineWidth(2); 

	TH2F *effplot = new TH2F("effplot",
                    "",
                    m12.getNumBins(),
                    m12.getLowerLimit(),
                    m12.getUpperLimit(),
                    m13.getNumBins(),
                    m13.getLowerLimit(),
                    m13.getUpperLimit());
	effplot->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
	effplot->GetYaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");


	TH1F *m12_data = new TH1F("m12_data","", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	m12_data->SetMarkerStyle(20);
	m12_data->SetMarkerSize(0.6);
	m12_data->SetLineWidth(2);
	m12_data->GetYaxis()->SetTitle("Efficiency");
	m12_data->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
	m12_data->SetLineColor(1);
	
	TH1F *m13_data = new TH1F("m13_data","", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m13_data->SetMarkerStyle(20);
	m13_data->SetMarkerSize(0.6);
	m13_data->SetLineWidth(2);
	m13_data->GetYaxis()->SetTitle("Efficiency");
	m13_data->GetXaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");
	m13_data->SetLineColor(1);
	
	TH1F *m23_data = new TH1F("m23_data","", m13.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
//	TH1F *m23_data = new TH1F("m23_data","", m12.getNumBins(), m23_lower, m23_upper);
	m23_data->SetMarkerStyle(20);
	m23_data->SetMarkerSize(0.6);
	m23_data->SetLineWidth(2);
	m23_data->GetYaxis()->SetTitle("Efficiency");
	m23_data->GetXaxis()->SetTitle("M^{2}_{K^{+}K_{S}^{0}} (GeV^{2}/c^{4})");
	m23_data->SetLineColor(1);


//	TChain *tr_data = new TChain("DTag");tr_data->Add(eff_filename.c_str());
//	tr_data->Project("m12_data","m_Kppi0_sq");
//	tr_data->Project("m13_data","m_Kspi0_sq");
//	tr_data->Project("m23_data","m_KpKs_sq");
//	tr_data->Draw("m_Kspi0_sq:m_Kppi0_sq>>effplot");

	TFile *f = TFile::Open(eff_filename.c_str());
	TTree *t = (TTree *)f->Get("DTag");
	assert(t);
	fptype m_Kppi0_sq, m_Kspi0_sq, m_KpKs_sq, m_recoil, dE_sig;
	t->SetBranchAddress("m_Kppi0_sq",&m_Kppi0_sq);
	t->SetBranchAddress("m_Kspi0_sq",&m_Kspi0_sq);
	t->SetBranchAddress("m_KpKs_sq",&m_KpKs_sq);
	t->SetBranchAddress("m_recoil",&m_recoil);
	t->SetBranchAddress("dE_sig",&dE_sig);

	for(int i = 0; i < t->GetEntries(); i++){
		t->GetEvent(i);
		if(m_recoil<mrec_min||m_recoil>mrec_max) continue;
		if(dE_sig<dE_min||dE_sig>dE_max) continue;
		if(!cpuDalitz(m_Kppi0_sq,m_Kspi0_sq)) continue;
		m12_data->Fill(m_Kppi0_sq);
		m13_data->Fill(m_Kspi0_sq);
		m23_data->Fill(m_KpKs_sq);
		effplot->Fill(m_Kppi0_sq,m_Kspi0_sq);
	}
	f->Close();


//     TRandom3 rr;
//     rr.SetSeed(233333);
//
//     int idx = effplot->GetXaxis()->GetNbins();
//     double xstep = (effplot->GetXaxis()->GetXmax()-effplot->GetXaxis()->GetXmin())/(double)idx;
//     int idy = effplot->GetYaxis()->GetNbins();
//     double ystep = (effplot->GetYaxis()->GetXmax()-effplot->GetYaxis()->GetXmin())/(double)idy;
//
//     double xleft, xright, ydown, yup;
//     double weight = 1;
//     for (int ix = 0; ix < idx; ++ix){
//          for (int iy = 0; iy < idy; ++iy){
//               weight = 1;
//               xleft = effplot->GetXaxis()->GetXmin() + ix*xstep;
//               xright = effplot->GetXaxis()->GetXmin() + (ix+1)*xstep;
//               ydown = effplot->GetYaxis()->GetXmin() + iy*ystep;
//               yup = effplot->GetYaxis()->GetXmin() + (iy+1)*ystep;
//               weight = calWeight(xleft,xright,ydown,yup,rr);
//
//               effplot->SetBinContent(ix+1,iy+1,weight*effplot->GetBinContent(ix+1,iy+1));
//          }
//     }

	double totalPdf = 0; 
	double totalDat = 0; 

	std::vector<Observable> vars;
	vars.push_back(m12);
	vars.push_back(m13);
	vars.push_back(eventNumber); 
	UnbinnedDataSet currData(vars); 

	int evtCounter = 0; 
	for (int i = 0; i < m12.getNumBins(); ++i) {
		m12.setValue(m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit())*(i + 0.5) / m12.getNumBins()); 
		for (int j = 0; j < m13.getNumBins(); ++j) {
			m13.setValue(m13.getLowerLimit()+ (m13.getUpperLimit()- m13.getLowerLimit())*(j + 0.5) / m13.getNumBins()); 
	
	          if (!cpuDalitz(m12.getValue(), m13.getValue())) continue;
	          eventNumber.setValue(evtCounter); 
	          evtCounter++;
	          currData.addEvent(); 
		}
	}

	total->setData(&currData);
	std::vector<std::vector<fptype>> pdfValues;
	pdfValues=total->getCompProbsAtDataPoints();
	for (unsigned int j = 0; j < pdfValues[0].size(); ++j) {
		double currm12 = currData.getValue(m12, j);
		double currm13 = currData.getValue(m13, j);
		double currm23 = cpuGetM23(currm12, currm13);

		effpdfplot->Fill(currm12, currm13, pdfValues[0][j]);
		m12_pdf_hist.Fill(currm12, pdfValues[0][j]);
		m13_pdf_hist.Fill(currm13, pdfValues[0][j]);
		m23_pdf_hist.Fill(currm23, pdfValues[0][j]); 

		totalPdf += pdfValues[0][j]; 
	}
	cout<<"ht integral: "<<totalPdf<<endl;

	//Drawing
	effpdfplot->Scale(effplot->Integral()/effpdfplot->Integral());
	TCanvas *foo = new TCanvas("c1","c1",800,700);
	foo->Divide(2,2);

	foo->cd(3);
	m12_data->Draw("P");
	m12_pdf_hist.Scale(m12_data->Integral()/m12_pdf_hist.Integral());
	for (unsigned int j = 0; j < pdfValues[0].size(); ++j)
		m12_pdf_hist.SetBinError(j,0);

	m12_pdf_hist.Draw("sameL");
	m12_data->Rebin(2);
	m12_pdf_hist.Rebin(2);

	foo->cd(4);
	m13_data->Draw("P");
	m13_pdf_hist.Scale(m13_data->Integral()/m13_pdf_hist.Integral());
	for (unsigned int j = 0; j < pdfValues[0].size(); ++j)
		m13_pdf_hist.SetBinError(j,0);

	m13_pdf_hist.Draw("sameL");
	m13_data->Rebin(2);
	m13_pdf_hist.Rebin(2);

//	foo->cd(4);
//	m23_data->Draw("P");
//	m23_pdf_hist.Scale(m23_data->Integral()/m23_pdf_hist.Integral());
//	for (unsigned int j = 0; j < pdfValues[0].size(); ++j)
//		m23_pdf_hist.SetBinError(j,0);
//
//	m23_pdf_hist.Draw("sameL");
//	m23_data->Rebin(10);
//	m23_pdf_hist.Rebin(10);


	TCanvas *foo2 = new TCanvas("c2","c2",700,600);
	foo2->cd();

	TH2F *pull = new TH2F("Pull","Pull Distribution",
		effplot->GetNbinsX(), effplot->GetXaxis()->GetXmin(), effplot->GetXaxis()->GetXmax(),
		effplot->GetNbinsY(), effplot->GetYaxis()->GetXmin(), effplot->GetYaxis()->GetXmax());
	pull->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
	pull->GetYaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");

	effplot->Rebin2D(2,2);
	effpdfplot->Rebin2D(2,2);
	pull->Rebin2D(2,2);

	int nonEmpty = 0;
	double Chi2 = 0;
	for (int i = 1; i <= effplot->GetNbinsX(); ++i){
		for (int j = 1; j <= effplot->GetNbinsY(); ++j){
			double dataTmp = effplot->GetBinContent(i,j);
			double fitTmp = effpdfplot->GetBinContent(i,j);
			if(dataTmp != 0){
				double val = (fitTmp-dataTmp)/sqrt(dataTmp);
				pull->SetBinContent(i,j,val);
				Chi2 += val*val;
				nonEmpty++;
			}
		}
	}
	cout << "Efficiency fit: chi and chi2 non-empty bin nums: " << nonEmpty << endl;
	cout << "chi2 = " << Chi2 << endl;
	pull->Draw("colz");
	foo2->Update();

	foo->cd(1);
	effplot->Draw("colz");

	foo->cd(2);
	effpdfplot->Draw("colz");
	if(m_draw_polyeff){
		foo2->SaveAs("plots/eff_pull.C");
		foo->SaveAs("plots/eff_fit_plot.C");
	}
	if(m_draw_smootheff){
		foo2->SaveAs("plots/smootheff_pull.C");
		foo->SaveAs("plots/smootheff_plot.C");
	}

	m12.setNumBins(oldBins1);
	m13.setNumBins(oldBins2);
}

GooPdf* fitEffPoly(UnbinnedDataSet &data){
	auto obs               = data.getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

//	Observable m12("m12",pow(pi0Mass+KpMass,2),pow(_mDp-KsMass,2));
//	Observable m13("m13",pow(pi0Mass+KsMass,2),pow(_mDp-KpMass,2));

	GooPdf* eff = NULL;
	UnbinnedDataSet *effdata = loadEffData(data,eff_filename);
	vector <fptype> init_val;

	if(m_float_polyeff){
		init_val.clear();
		for (int j = 0; j < 9; ++j)
			init_val.push_back(fRand(-1,1));
		init_val.push_back(fRand(0,50));
		init_val.push_back(fRand(0,1));
		init_val.push_back(fRand(0,50));
		init_val.push_back(fRand(0,1));
		init_val.push_back(fRand(0,50));
		init_val.push_back(fRand(0,1));

		init_val.push_back(fRand(0,50));
		init_val.push_back(fRand(0,1));
		init_val.push_back(fRand(0,50));
		init_val.push_back(fRand(0,1));
		init_val.push_back(fRand(0,50));
		init_val.push_back(fRand(0,1));
	}
	else{
		init_val.clear();
		//trial 1
          init_val.push_back(-0.3095447475032);
          init_val.push_back(0.1712820541676);
          init_val.push_back(-0.1219472651623);
          init_val.push_back(-0.168509897536);
          init_val.push_back(0.1357869374861);
          init_val.push_back(-0.04812401129686);
          init_val.push_back(-0.1413035963685);
          init_val.push_back(0.1079023080486);
          init_val.push_back(-0.2398233259839);
  
          init_val.push_back(3.503892613348);
          init_val.push_back(0.6510782063868);
          init_val.push_back(2.295115152977);
          init_val.push_back(0);
          init_val.push_back(0.03929480684439);
          init_val.push_back(0.280717494653);
	}

	eff = makeEfficiencyPoly(m12,m13,init_val);
	eff->setData(effdata);

	FitManager effpdf(eff);
	effpdf.fit();
	fptype nll = eff->calculateNLL();

	eff->setParameterConstantness(true);

	cout << "init values: " ;
	for (int k = 0; k < init_val.size(); ++k)
		cout << init_val[k] << "      ";
	cout << endl;
	cout << "NLL = " << endl << nll << endl;

	if(m_draw_polyeff)
		makeEffPlot(eff,&data);
	return eff;
}

double runDataFit(Amp3Body *signal, UnbinnedDataSet *data, bool m_err) {
	auto obs               = data->getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

	signal->setData(data);
	signal->setDataSize(data->getNumEvents());
	FitManager *datapdf= new FitManager(signal);
	datapdf->setMaxCalls(70000);
	datapdf->fit();

	cout << fixed << setprecision(8);

	std::vector <std::vector<double>> vec_phi = datapdf->printParams();

	//fit fractions
	vector <vector<fptype>> fracMat;
	fracMat = signal->fit_fractions();
	vector <fptype> fracList;
	fracList.clear();

	const int num_res = fracMat.size();
	cout << "Un-diagonal elements in fraction matrix may be meaningless!!!" << endl;
	for (int i = 0; i < num_res; ++i){
		for (int j = 0; j < num_res; ++j){
			cout << fracMat[i][j] << ",  ";
			if(i == j) fracList.push_back(fracMat[i][j]);
		}
		cout << endl;
	}

	double nll = signal->calculateNLL();
//	cout << "NLL = " << nll << endl;

	//Draw option
	if (m_draw_data){
		int oldBins1 = m12.getNumBins();
		int oldBins2 = m13.getNumBins();
		m12.setNumBins(drawBinM12);
		m13.setNumBins(drawBinM13);

//		double bin_m12[70] = {0.368, 0.444623188485189, 0.48379710157630235, 0.5341195654555669, 0.5724927539407559, 0.6149106282839634, 0.644619565455567, 0.6743449277585291, 0.6917884063052095, 0.7148795988820618, 0.7275951692938374, 0.7389613528898877, 0.7464611802603239, 0.7525502071867796, 0.7586392341132351, 0.7635826086629881, 0.7684538300247593, 0.7733250513865305, 0.7781352652505958, 0.7828711750822834, 0.7876070849139711, 0.7923429947456588, 0.7952943838247574, 0.7979583331050817, 0.800622282385406, 0.8032862316657303, 0.8059501809460546, 0.808614130226379, 0.8116693283791978, 0.8151487724286995, 0.8186282164782014, 0.8221076605277031, 0.8255871045772049, 0.8303753626681648, 0.83605845446619, 0.8417415462642153, 0.848466918961286, 0.8558796469780737, 0.864766044845469, 0.8769440986983802, 0.893571013571808, 0.9106202889658834, 0.9217934777249746, 0.936884056791478, 0.9566681151481101, 0.9749082117135399, 1.0168333324203267, 1.0424318826606742, 1.0731304334488252, 1.1073064173121259, 1.1316625250179482, 1.1476482208242176, 1.1730890261438787, 1.2023022767068436, 1.2301217381777325, 1.2990169074715696, 1.355765699795029, 1.409681158785155, 1.4437632844419475, 1.4657072460592444, 1.4925130429066396, 1.5334115936947907, 1.5625169074715695, 1.6143381640394963, 1.635362318613746, 1.6698057968473954, 1.7021304345444332, 1.7346028984236979, 1.7958357486864147, 1.983};
//
//		double bin_m13[70] = {0.385, 0.41336618357487925, 0.4240579710144927, 0.4340869565217391, 0.44519806763285025, 0.4568188405797102, 0.4689968944099379, 0.4830374396135266, 0.4964570791527313, 0.5112434782608696, 0.5276297760210804, 0.5442689210950081, 0.561159420289855, 0.5751081382385731, 0.5882229654403568, 0.6013377926421404, 0.6159894598155469, 0.6407615283267457, 0.6562608695652175, 0.673236231884058, 0.6888050065876152, 0.710463768115942, 0.7199355877616747, 0.7306666666666668, 0.7428447204968944, 0.7550227743271222, 0.769437417654809, 0.7825536231884058, 0.793657004830918, 0.8017757073844031, 0.8098944099378882, 0.8193486312399356, 0.8291845410628019, 0.8405507246376811, 0.8623961352657005, 0.8766038647342995, 0.9155478260869567, 0.9463079710144927, 0.9778178053830228, 1.043072463768116, 1.091177536231884, 1.179400966183575, 1.224898550724638, 1.2936376811594201, 1.3290851449275363, 1.3576956521739132, 1.391777777777778, 1.412596618357488, 1.4292065217391305, 1.4505181159420293, 1.4718297101449276, 1.4953754940711463, 1.5117359098228664, 1.5286243032329991, 1.544202898550725, 1.5777563405797104, 1.5884121376811595, 1.6081552795031058, 1.6291755233494365, 1.6468393719806766, 1.661570652173913, 1.6809214975845412, 1.6952318840579712, 1.7107312252964428, 1.7360893719806765, 1.7617184265010353, 1.7740372670807456, 1.7958357487922707, 1.8226438923395447, 1.983};
//
//		double bin_m23[70] = {0.9793836006789989, 1.0324554847369698, 1.0726679162989667, 1.112127252852912, 1.163643214205569, 1.2041931175872116, 1.2536949050268251, 1.2936231497933306, 1.3162248808722359, 1.3348722172797893, 1.3548603294574668, 1.3706631451924565, 1.3864659609274461, 1.4060139526458726, 1.4314246810479057, 1.448758057200738, 1.4625855209688539, 1.4751119420638623, 1.487914557200738, 1.5095829122732018, 1.5234103760413178, 1.5359669919833467, 1.5470289629978395, 1.55980671179011, 1.5745560064761004, 1.592880726282864, 1.6071250670046595, 1.6225764702442165, 1.6497305185533953, 1.664788866379482, 1.68322548473697, 1.7031492396776815, 1.7198551327700962, 1.743970412273202, 1.7674981724840055, 1.7993335523698202, 1.8287895137224772, 1.8746176793539473, 1.9401125861862454, 2.017642296331173, 2.1079604122732025, 2.2499117166210283, 2.351579542707985, 2.385697182459537, 2.4205090499543616, 2.4434429243504967, 2.4680250821604806, 2.490217487372016, 2.510330161943821, 2.5319671111459883, 2.5508094435470388, 2.5624536235622943, 2.578748132958314, 2.598860807530119, 2.618377658650014, 2.6374004122732018, 2.6578282078490982, 2.669472387864354, 2.6829324654132987, 2.7010427129978396, 2.7148701767659555, 2.728697640534072, 2.7425251043021874, 2.7534551948818975, 2.7639904053718904, 2.785568302772397, 2.813880412273202, 2.83231703063069, 2.8979641803891454, 3.075083600678999};


//		datapdf->printCovMat();
		ProdPdf prodpdf{"prodpdf", {signal}};
	
		DalitzPlotter plotter(&prodpdf, signal);
//		datapdf->printCovMat();
 
		TH2F *dalitzplot = plotter.make2D();
		dalitzplot->SetTitle("dalitzFit");
		dalitzplot->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
		dalitzplot->GetYaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");
	
		TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
//		TH1F m12_pdf_hist("m12_pdf_hist", "", 69, bin_m12);
		m12_pdf_hist.SetStats(false); 
		m12_pdf_hist.SetLineColor(kRed); 
		m12_pdf_hist.SetLineWidth(2); 
	
		TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
//		TH1F m13_pdf_hist("m13_pdf_hist", "", 69, bin_m13);
		m13_pdf_hist.SetStats(false); 
		m13_pdf_hist.SetLineColor(kRed); 
		m13_pdf_hist.SetLineWidth(2); 
	
		TH1F m23_pdf_hist("m23_pdf_hist", "", m13.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
//		TH1F m23_pdf_hist("m23_pdf_hist", "", 69, bin_m23);
		m23_pdf_hist.SetStats(false); 
		m23_pdf_hist.SetLineColor(kRed); 
		m23_pdf_hist.SetLineWidth(2); 
	
		fptype m12_tmp,m13_tmp,m23_tmp;
	
		UnbinnedDataSet toyMC({m12, m13, eventNumber});
		plotter.fillDataSetMC(toyMC, 5000000);

		cout << "toyMC size = " << eventNumber.getValue() << endl;
		for (int i = 0; i < eventNumber.getValue(); ++i){
			m12_tmp = toyMC.getValue(m12, i);
			m12_pdf_hist.Fill(m12_tmp);
	
			m13_tmp = toyMC.getValue(m13, i);
			m13_pdf_hist.Fill(m13_tmp);
	
			m23_tmp = cpuGetM23(m12_tmp, m13_tmp);
			m23_pdf_hist.Fill(m23_tmp);
		}

/*
		//trial for resonances distribution 
		vector<PdfBase*> comps;
		comps.clear();
		comps.push_back(signal);
		ProdPdf* sigtemp = new ProdPdf("signal_temp", comps);
//		GooPdf * sigtemp = (GooPdf *) signal;

		std::vector<Observable> vars;
		vars.push_back(m12);
		vars.push_back(m13);
		vars.push_back(eventNumber); 
		UnbinnedDataSet currData(vars); 
	
		int evtCounter = 0; 
		for (int i = 0; i < m12.getNumBins(); ++i) {
			m12.setValue(m12.getLowerLimit() + (m12.getUpperLimit() - m12.getLowerLimit())*(i + 0.5) / m12.getNumBins()); 
			for (int j = 0; j < m13.getNumBins(); ++j) {
				m13.setValue(m13.getLowerLimit()+ (m13.getUpperLimit()- m13.getLowerLimit())*(j + 0.5) / m13.getNumBins()); 
		
		          if (!cpuDalitz(m12.getValue(), m13.getValue())) continue;
		          eventNumber.setValue(evtCounter); 
		          evtCounter++;
		          currData.addEvent(); 
			}
		}
		sigtemp->setData(&currData);
		std::vector<std::vector<fptype>> pdfValues;
		pdfValues = sigtemp->getCompProbsAtDataPoints();

		const int nRes = pdfValues.size();
		cout << "num of nRes = " << nRes << endl;
		TH1F* m12_pdf_hist_res[nRes];
		TH1F* m13_pdf_hist_res[nRes];
		TH1F* m23_pdf_hist_res[nRes];
		for (int i=0;i<nRes;i++){
			sprintf(strbuffer, "%s_res%d", m12_pdf_hist.GetName(), i);
			m12_pdf_hist_res[i] = (TH1F*)m12_pdf_hist.Clone(strbuffer);
			m12_pdf_hist_res[i]->Reset();

			sprintf(strbuffer, "%s_res%d", m13_pdf_hist.GetName(), i);
			m13_pdf_hist_res[i] = (TH1F*)m13_pdf_hist.Clone(strbuffer);
			m13_pdf_hist_res[i]->Reset();

			sprintf(strbuffer, "%s_res%d", m23_pdf_hist.GetName(), i);
			m23_pdf_hist_res[i] = (TH1F*)m23_pdf_hist.Clone(strbuffer);
			m23_pdf_hist_res[i]->Reset();
		}
		 for (unsigned int j = 0; j < pdfValues[0].size(); ++j) {
			double currm12 = currData.getValue(m12, j);
			double currm13 = currData.getValue(m13, j);
			double currm23 = cpuGetM23(currm12, currm13);
			for (int i=0;i<nRes;i++){
				m12_pdf_hist_res[i]->Fill(currm12, pdfValues[i][j]);
				m13_pdf_hist_res[i]->Fill(currm13, pdfValues[i][j]);
				m23_pdf_hist_res[i]->Fill(currm23, pdfValues[i][j]);
//				cout<<"For Res #"<<i<<", "<<j<<": "<<currm12<<','<<currm13<<": "<<pdfValues[i][j]<<endl;
			}
		}
*/

		//for data comparison
		TString a("Events/"); TString c(" GeV^{2}/c^{4})");
		TH1F *m12_data = new TH1F("m12_data","", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
//		TH1F *m12_data = new TH1F("m12_data","", 69, bin_m12);
		char b_m12[20];  sprintf(b_m12, "(%.2f",(m12.getUpperLimit()-m12.getLowerLimit())/m12.getNumBins()*drawBinM12/50.0);//Rebin() later
		TString ytitle_m12 = a + b_m12 + c;
//		TString ytitle_m12 = a + c;

		m12_data->SetMarkerStyle(20);
		m12_data->SetMarkerSize(0.6);
		m12_data->SetLineWidth(2);
		m12_data->GetYaxis()->SetTitle(ytitle_m12);
		m12_data->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
		m12_data->SetLineColor(1);
	
		TH1F *m13_data = new TH1F("m13_data","", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
//		TH1F *m13_data = new TH1F("m13_data","", 69, bin_m13);
		char b_m13[20];  sprintf(b_m13, "(%.2f",(m13.getUpperLimit()-m13.getLowerLimit())/m13.getNumBins()*drawBinM13/50.0);//Rebin() later
		TString ytitle_m13 = a + b_m13 + c;
//		TString ytitle_m13 = a + c;
		m13_data->SetMarkerStyle(20);
		m13_data->SetMarkerSize(0.6);
		m13_data->SetLineWidth(2);
		m13_data->GetYaxis()->SetTitle(ytitle_m13);
		m13_data->GetXaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");
		m13_data->SetLineColor(1);
	
		TH1F *m23_data = new TH1F("m23_data","", m13.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
//		TH1F *m23_data = new TH1F("m23_data","", 69, bin_m23);
		char b_m23[20];  sprintf(b_m23, "(%.2f",(getM23UpperLimit(m12,m13)-getM23LowerLimit(m12,m13))/m13.getNumBins()*drawBinM13/50.0);//Rebin() later
		TString ytitle_m23 = a + b_m23 + c;
//		TString ytitle_m23 = a + c;
		m23_data->SetMarkerStyle(20);
		m23_data->SetMarkerSize(0.6);
		m23_data->SetLineWidth(2);
		m23_data->GetYaxis()->SetTitle(ytitle_m23);
		m23_data->GetXaxis()->SetTitle("M^{2}_{K^{+}K_{S}^{0}} (GeV^{2}/c^{4})");
		m23_data->SetLineColor(1);
	
		//Drawing 2D pull distribution
		TH2F *dalitzData = new TH2F ("dalitzData","dalitzData",
			dalitzplot->GetNbinsX(), dalitzplot->GetXaxis()->GetXmin(), dalitzplot->GetXaxis()->GetXmax(),
			dalitzplot->GetNbinsY(), dalitzplot->GetYaxis()->GetXmin(), dalitzplot->GetYaxis()->GetXmax());
	
		TFile *f = TFile::Open(data_filename.c_str());
		TTree *t = (TTree *)f->Get("DTag");
		assert(t);
		fptype m_Kppi0_sq, m_Kspi0_sq, m_KpKs_sq, m_recoil, dE_sig;
		t->SetBranchAddress("m_Kppi0_sq",&m_Kppi0_sq);
		t->SetBranchAddress("m_Kspi0_sq",&m_Kspi0_sq);
		t->SetBranchAddress("m_KpKs_sq",&m_KpKs_sq);
		t->SetBranchAddress("m_recoil",&m_recoil);
		t->SetBranchAddress("dE_sig",&dE_sig);

		for(int i = 0; i < t->GetEntries(); i++){
			t->GetEvent(i);
			if(m_recoil<mrec_min||m_recoil>mrec_max) continue;
			if(dE_sig<dE_min||dE_sig>dE_max) continue;
			if(!cpuDalitz(m_Kppi0_sq,m_Kspi0_sq)) continue;
			m12_data->Fill(m_Kppi0_sq);
			m13_data->Fill(m_Kspi0_sq);
			m23_data->Fill(m_KpKs_sq);
			dalitzData->Fill(m_Kppi0_sq,m_Kspi0_sq);
		}
		cout << dalitzData->GetEntries() << endl;
		f->Close();

		TFile *fnew = new TFile("plots/hist_results.root","recreate");
		dalitzData->Write();
		dalitzplot->Write();
		m12_data->Write();
		m12_pdf_hist.Write();
		m13_data->Write();
		m13_pdf_hist.Write();
		m23_data->Write();
		m23_pdf_hist.Write();
		fnew->Close();

		TCanvas tmp("c2","c2",800,700);
		tmp.Divide(2,2);
		tmp.cd(1);
		gStyle->SetPalette(52,0);
		dalitzData->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
		dalitzData->GetYaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");
		dalitzData->Draw("colz");
	
		tmp.cd(2);
		dalitzplot->Scale(dalitzData->Integral()/dalitzplot->Integral());
		dalitzplot->Draw("colz");

		int nonEmpty = 0;
		double Chi2 = 0;
///*
		TH2F *pull = new TH2F("Pull","Pull Distribution",
			dalitzplot->GetNbinsX(), dalitzplot->GetXaxis()->GetXmin(), dalitzplot->GetXaxis()->GetXmax(),
			dalitzplot->GetNbinsY(), dalitzplot->GetYaxis()->GetXmin(), dalitzplot->GetYaxis()->GetXmax());
		pull->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
		pull->GetYaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");
	
		dalitzplot->Rebin2D(drawBinM12/20.0,drawBinM13/20.0);
		dalitzData->Rebin2D(drawBinM12/20.0,drawBinM13/20.0);
		pull->Rebin2D(drawBinM12/20.0,drawBinM13/20.0);
	
		for (int i = 1; i <= dalitzplot->GetNbinsX(); ++i){
			for (int j = 1; j <= dalitzplot->GetNbinsY(); ++j){
				double dataTmp = dalitzData->GetBinContent(i,j);
				double fitTmp = dalitzplot->GetBinContent(i,j);
				if(dataTmp != 0){
	//				cout << "data bin = " << dataTmp << endl;
	//				cout << "fit bin = " << fitTmp << endl;
					double val = (fitTmp-dataTmp)/sqrt(dataTmp);
					pull->SetBinContent(i,j,val);
					Chi2 += val*val;
					nonEmpty++;
				}
	//			cout << "pull = " << pull->GetBinContent(i,j) << endl;
			}
		}
		cout << "pull: chi and chi2 non-empty bin nums = " << nonEmpty << endl;
		cout << "      chi2 = " << Chi2 << endl;
//*/


//		ChisqInfo* ChiSq = getAdaptiveChisquare(dalitzData, dalitzplot);
//		cout << "Chisquare: " <<  ChiSq->chisq << " / " << ChiSq->dof << endl;
//		TH2F *pull = ChiSq->contribPlot;
//		pull->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})");
//		pull->GetYaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");

		//Drawing data & fit comparison
		TCanvas *foo = new TCanvas("c1","c1",800,700);
		foo->Divide(2,2);
		foo->cd(1);

//		pull->Draw("colz");

		dalitzData->GetYaxis()->SetTitle("M^{2}_{K_{S}^{0}#pi^{0}} (GeV^{2}/c^{4})");
		dalitzData->GetXaxis()->SetTitle("M^{2}_{K^{+}#pi^{0}} (GeV^{2}/c^{4})    ");
	
		dalitzData->GetZaxis()->SetTitle("Events");
	
		dalitzData->SetMarkerStyle(24);
		dalitzData->SetMarkerSize(0.7);
	
		dalitzData->GetXaxis()->SetNdivisions(505);
		dalitzData->GetXaxis()->SetLabelFont(22);
		dalitzData->GetXaxis()->SetLabelOffset(0.015);
		dalitzData->GetXaxis()->SetLabelSize(0.04);
		dalitzData->GetXaxis()->SetTitleSize(0.04);
		dalitzData->GetXaxis()->SetTitleOffset(1.8);
		dalitzData->GetXaxis()->SetTitleFont(22);

		dalitzData->GetYaxis()->SetNdivisions(505);
		dalitzData->GetYaxis()->SetLabelFont(22);
		dalitzData->GetYaxis()->SetLabelOffset(0.015);
		dalitzData->GetYaxis()->SetLabelSize(0.04);
		dalitzData->GetYaxis()->SetTitleSize(0.04);
		dalitzData->GetYaxis()->SetTitleOffset(1.8);
		dalitzData->GetYaxis()->SetTitleFont(22);
	
		dalitzData->GetZaxis()->SetLabelFont(22);
		dalitzData->GetZaxis()->SetLabelSize(0.04);
		dalitzData->GetZaxis()->SetTitleSize(0.04);
		dalitzData->GetZaxis()->SetTitleFont(22);
	
		dalitzData->Draw("lego2Z");
		dalitzplot->Draw("surf3SAME");

		foo->cd(2);
		m12_data->Draw("e");
		m12_pdf_hist.Scale(m12_data->Integral()/m12_pdf_hist.Integral());
		m12_pdf_hist.Draw("Hsame");
		m12_data->Rebin(drawBinM12/50.0);
		m12_pdf_hist.Rebin(drawBinM12/50.0);

		nonEmpty = 0;
		Chi2 = 0;
		for (int i = 1; i <= m12_data->GetNbinsX(); ++i){
			double dataTmp = m12_data->GetBinContent(i);
			double fitTmp = m12_pdf_hist.GetBinContent(i);
			if(dataTmp != 0){
				double val = (fitTmp-dataTmp)/sqrt(dataTmp);
				Chi2 += val*val;
				nonEmpty++;
			}
		}
		cout << "m12: chi and chi2 non-empty bin nums = " << nonEmpty << endl;
		cout << "     chi2 = " << Chi2 << endl;


		foo->cd(3);
		m13_data->Draw("e");
		m13_pdf_hist.Scale(m13_data->Integral()/m13_pdf_hist.Integral());
		m13_pdf_hist.Draw("Hsame");
		m13_data->Rebin(drawBinM13/50.0);
		m13_pdf_hist.Rebin(drawBinM13/50.0);

		nonEmpty = 0;
		Chi2 = 0;
		for (int i = 1; i <= m13_data->GetNbinsX(); ++i){
			double dataTmp = m13_data->GetBinContent(i);
			double fitTmp = m13_pdf_hist.GetBinContent(i);
			if(dataTmp != 0){
				double val = (fitTmp-dataTmp)/sqrt(dataTmp);
				Chi2 += val*val;
				nonEmpty++;
			}
		}
		cout << "m13: chi and chi2 non-empty bin nums = " << nonEmpty << endl;
		cout << "     chi2 = " << Chi2 << endl;


	
		foo->cd(4);
		m23_data->Draw("e");
		m23_pdf_hist.Scale(m23_data->Integral()/m23_pdf_hist.Integral());
		m23_pdf_hist.Draw("Hsame");
		m23_data->Rebin(drawBinM13/50.0);
		m23_pdf_hist.Rebin(drawBinM13/50.0);

		nonEmpty = 0;
		Chi2 = 0;
		for (int i = 1; i <= m23_data->GetNbinsX(); ++i){
			double dataTmp = m23_data->GetBinContent(i);
			double fitTmp = m23_pdf_hist.GetBinContent(i);
			if(dataTmp != 0){
				double val = (fitTmp-dataTmp)/sqrt(dataTmp);
				Chi2 += val*val;
				nonEmpty++;
			}
		}
		cout << "m23: chi and chi2 non-empty bin nums = " << nonEmpty << endl;
		cout << "     chi2 = " << Chi2 << endl;



		if(!m_effPoly)
			foo->SaveAs("plots/dalitz_with_projections.C");
		else
			foo->SaveAs("plots/dalitz_with_projections_effPoly.C");

//		for (int i=0;i<nRes;i++){
//			m12_pdf_hist_res[i]->SetLineColor(i+1);
//			m12_pdf_hist_res[i]->SetLineWidth(1.2);
//			m12_pdf_hist_res[i]->Rebin(2);
//			m12_pdf_hist_res[i]->Draw("Hsame");
//
//			foo->cd(3);
//			m13_pdf_hist_res[i]->SetLineColor(i+1);
//			m13_pdf_hist_res[i]->SetLineWidth(1.2);
//			m13_pdf_hist_res[i]->Rebin(2);
//			m13_pdf_hist_res[i]->Draw("Hsame");
//
//			foo->cd(4);
//			m23_pdf_hist_res[i]->SetLineColor(i+1);
//			m23_pdf_hist_res[i]->SetLineWidth(1.2);
//			m23_pdf_hist_res[i]->Rebin(2);
//			m23_pdf_hist_res[i]->Draw("Hsame");
//		}

//		TCanvas *foo2 = new TCanvas("c2","c2",1200,400);
//		foo2->Divide(3,1);
//
//		foo2->cd(1);
//		m12_pdf_hist_res[1]->Rebin(2);
//		for (unsigned int j = 0; j < pdfValues[0].size(); ++j)
//	          m12_pdf_hist_res[1]->SetBinError(j,0);
//		m12_pdf_hist_res[1]->Scale(m12_data->Integral()/m12_pdf_hist_res[1]->Integral());
//
//		m12_data->Draw("e");
//		m12_pdf_hist_res[1]->Draw("Hsame");
//
//		foo2->cd(2);
//		m13_pdf_hist_res[1]->Rebin(2);
//		for (unsigned int j = 0; j < pdfValues[0].size(); ++j)
//	          m13_pdf_hist_res[1]->SetBinError(j,0);
//		m13_pdf_hist_res[1]->Scale(m13_data->Integral()/m13_pdf_hist_res[1]->Integral());
//
//		m13_data->Draw("e");
//		m13_pdf_hist_res[1]->Draw("Hsame");
//
//		foo2->cd(3);
//		m23_pdf_hist_res[1]->Rebin(2);
//		for (unsigned int j = 0; j < pdfValues[0].size(); ++j)
//	          m23_pdf_hist_res[1]->SetBinError(j,0);
//		m23_pdf_hist_res[1]->Scale(m23_data->Integral()/m23_pdf_hist_res[1]->Integral());
//
//		m23_data->Draw("e");
//		m23_pdf_hist_res[1]->Draw("Hsame");
//
//		foo2->SaveAs("plots/temp.C");


		m12.setNumBins(oldBins1);
		m13.setNumBins(oldBins2);
   }

	if(m_err){
		double mean[num_res];
		double rms[num_res];
		double errList[num_res];
		vector <double> fractions[num_res];
	
		for (int ii = 0; ii < num_res; ii++) mean[ii] = rms[ii] = 0;
	
		fptype mean_rate = 0;
		fptype rms_rate = 0;
		vector <fptype> rate; rate.clear();
	
		const int nSamples = 1000;
		datapdf->setRandMinuitValues(nSamples);
		for (int ii = 0; ii < nSamples; ++ii){
			datapdf->loadSample(ii);
		
			vector <vector<fptype>> fracListNew;
			fracListNew = signal->fit_fractions();
	
			double Tmp = (fracListNew[0][0]/0.333)/(fracListNew[1][1]/0.3323/0.5);
			rate.push_back(Tmp);
			mean_rate += Tmp;
			rms_rate += Tmp*Tmp;
		
			for (int i = 0; i < fracListNew.size(); ++i){
				for (int j = 0; j < fracListNew.size(); ++j){
					//cout << fracListNew[i][j] << ",  ";
					if(i == j){
						fractions[i].push_back(fracListNew[i][j]);
						mean[i] += fracListNew[i][j]; 
						rms[i] += fracListNew[i][j]*fracListNew[i][j];
					}
				}
				//cout << endl;
			}
		}
	
		//for reasonance
		TH1F* hFracs[num_res];
		TFile * froot = new TFile("plots/sigfractionHists.root", "recreate");
	
		for (int ii = 0; ii < num_res; ii++) {
			mean[ii] /= nSamples;
			rms[ii] = sqrt(rms[ii]/nSamples-mean[ii]*mean[ii]);
			sprintf(strbuffer, "hfrac_res_%d", ii);
			hFracs[ii] = new TH1F(strbuffer, "", 100, mean[ii]-4*rms[ii], mean[ii]+4*rms[ii]);
	
			for (int jj = 0; jj < nSamples; jj++)
				hFracs[ii]->Fill(fractions[ii][jj]);
			hFracs[ii]->Write();
			errList[ii] = hFracs[ii]->GetRMS();//gaus1->GetParameter(2);
		}
	
		//for K*(892)+ / Kbar*(892)0 rate
		mean_rate /= nSamples;
		rms_rate = sqrt(rms_rate/nSamples-mean_rate*mean_rate);
		TH1F* hRate = new TH1F("hRate","",100,mean_rate-4*rms_rate,mean_rate+4*rms_rate);
	
		for(int jj = 0; jj < nSamples; jj++)
			hRate->Fill(rate[jj]);
		hRate->Write();
	
		fptype fraction_rate = (fracList[0]/0.333) / (fracList[1]/0.3323/0.5);
		fptype fraction_rate_err = hRate->GetRMS();
	
		froot->Close();
		std::cout << std::endl;
		for (int ii=0;ii<num_res;ii++) 
			std::cout<<"Fit fraction for #"<<ii<<": "<<fracList[ii]<< " +- "<<errList[ii]<<std::endl;
		std::cout<<"Fit fraction for K*(892)+ / Kbar*(892)0 rate: " << fraction_rate << " +- " << fraction_rate_err << std::endl;

		if(m_test_for_pull){
			datapdf->printOriginalParams();
		}
		if(m_print_for_sys){
			cout << "tempsys" << endl;
			cout << vec_phi[0][0] << endl;
			cout << vec_phi[1][0] << endl;
			cout << fracList[0] << endl;
			cout << errList[0] << endl;
			cout << fracList[1] << endl;
			cout << errList[1] << endl;
		}
	}

	delete datapdf;
	return nll;
//	return 0;
}

void gen_test_pdf(bool use_eff = false,bool save_toy = false){
	Observable m12("m12", 0.3, 2);
	Observable m13("m13", 0.3, 2);
	EventNumber eventNumber("eventNumber");
	m12.setNumBins(BinNumsM12);
	m13.setNumBins(BinNumsM13);
//	m12.setNumBins(200);
//	m13.setNumBins(200);


	DecayInfo3 dtop0pp;
	dtop0pp.motherMass   = _mDp;
	dtop0pp.daug1Mass    = pi0Mass;
	dtop0pp.daug2Mass    = KpMass;
	dtop0pp.daug3Mass    = KsMass;
	dtop0pp.meson_radius = _radius;

	Variable K892pMass("K892p_mass", 0.89176);
	Variable K892pWidth("K892p_width", 0.0503);
	ResonancePdf *K892p = new Resonances::RBW(
		"K892p", 
		Variable("K892p_amp_real", 1.), 
		Variable("K892p_amp_imag", 0.), 
		K892pMass, 
		K892pWidth, 
		1, 
		PAIR_12);

	Variable K892zeroMass("K892zero_mass", 0.89555);
	Variable K892zeroWidth("K892zero_width", 0.0473);
	ResonancePdf *K892zero = new Resonances::RBW(
		"K892zero",
//		Variable("K892zero_amp_real", -0.3711157847882),
//		Variable("K892zero_amp_imag", 0.2194154470146),
		Variable("K892zero_amp_real", -0.3664247378405),
		Variable("K892zero_amp_imag", 0.162760084231),
		K892zeroMass,
		K892zeroWidth,
		1,
		PAIR_13);

	Variable a980pMass("a980pMass",m_a980p,0.020,0.9,1.06);
	ResonancePdf *a980p = new Resonances::FLATTE(
		"a980p",
		Variable("a980p_amp_real", 0.8655022022291),
		Variable("a980p_amp_imag", 0.9397819142879),
		a980pMass,
		Variable("g1",c_g1),
		Variable("rg2og1",c_g2og1),
		PAIR_23,
		false);

	ResonancePdf *nonr = new Resonances::NonRes(
		"nonr",
		Variable("nonr_amp_real", 1.201806574815),
		Variable("nonr_amp_imag", 3.550475936164)
		);


//     Variable K1410zeroMass("K1410zero_mass", 1.421);
//     Variable K1410zeroWidth("K1410zero_width", 0.236);
//     ResonancePdf *K1410zero = new Resonances::RBW(
//          "K1410zero",
//          Variable("K1410zero_amp_real", -1.429763770204),
//          Variable("K1410zero_amp_imag", 0.5947540558733),
//          K1410zeroMass,
//          K1410zeroWidth,
//          1,
//          PAIR_13);
//
     Variable SwaveKppi0Mass("SwaveKppi0_Mass", 1.425);
     Variable SwaveKppi0Width("SwaveKppi0_Width", 0.27);
     ResonancePdf *SwaveKppi0 = new Resonances::LASS(
          "SwaveKppi0",
          Variable("SwaveKppi0_amp_real", 1.452363736485),
          Variable("SwaveKppi0_amp_imag", -0.8115999519399),
          SwaveKppi0Mass,
          SwaveKppi0Width,
          0,
          PAIR_12);

     Variable SwaveKspi0Mass("SwaveKspi0_Mass", 1.425);
     Variable SwaveKspi0Width("SwaveKspi0_Width", 0.27);
     ResonancePdf *SwaveKspi0 = new Resonances::LASS(
          "SwaveKspi0",
          Variable("SwaveKspi0_amp_real", 2.251018780733),
          Variable("SwaveKspi0_amp_imag", 1.234646710115),
          SwaveKspi0Mass,
          SwaveKspi0Width,
          0,
          PAIR_13);


	dtop0pp.resonances.push_back(K892p);
	dtop0pp.resonances.push_back(K892zero);
//	dtop0pp.resonances.push_back(a980p);
//	dtop0pp.resonances.push_back(nonr);
//	dtop0pp.resonances.push_back(K1410zero);
	dtop0pp.resonances.push_back(SwaveKppi0);
	dtop0pp.resonances.push_back(SwaveKspi0);

	GooPdf *eff = NULL;
	if(!use_eff){
		vector<Variable> offsets;
		vector<Observable> observables;
		vector<Variable> coefficients;
		observables.push_back(m12);
		observables.push_back(m13);
		offsets.push_back(constantZero);
		offsets.push_back(constantZero);
		coefficients.push_back(constantOne);
		eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);
	}

	else{
		TFile *f = TFile::Open(eff_filename.c_str()); 
		TH2D* h2 = (TH2D*)f->Get("th2d_dalitz"); //assert(h2);
		
		int oldBins1 = m12.getNumBins();
		int oldBins2 = m13.getNumBins();
		m12.setNumBins(h2->GetNbinsX());
		m13.setNumBins(h2->GetNbinsY());
	
		BinnedDataSet *binEffData = new BinnedDataSet({m12, m13, eventNumber});
	
		int num = m12.getNumBins()*m13.getNumBins();
	
		for (int i = 0;i < num; i++){
			int iy = i/m12.getNumBins();
			int ix = i%m12.getNumBins();
			double content = h2->GetBinContent(ix+1,iy+1);
			binEffData->setBinContent(i, content);
		}
		SmoothHistogramPdf* ret = new SmoothHistogramPdf("test_eff", binEffData, Variable("smoothConst", 10));
		m12.setNumBins(oldBins1);
		m13.setNumBins(oldBins2);
		f->Close();
		eff = ret;
	}

	Amp3Body *signal = new Amp3Body("signalPDF_test", m12, m13, eventNumber, dtop0pp, eff);
	signal->setParameterConstantness(true);

	int num = 200;
	if(save_toy){
		for (int i = 0; i < num; ++i){
			if(i==Int_t(num*0.1)) cout<<"*******************************completed 10%"<<"**************************************"<<endl;
			if(i==Int_t(num*0.2)) cout<<"*******************************completed 20%"<<"**************************************"<<endl;
			if(i==Int_t(num*0.3)) cout<<"*******************************completed 30%"<<"**************************************"<<endl;
			if(i==Int_t(num*0.4)) cout<<"*******************************completed 40%"<<"**************************************"<<endl;
			if(i==Int_t(num*0.5)) cout<<"*******************************completed 50%"<<"**************************************"<<endl;
			if(i==Int_t(num*0.6)) cout<<"*******************************completed 60%"<<"**************************************"<<endl;
			if(i==Int_t(num*0.7)) cout<<"*******************************completed 70%"<<"**************************************"<<endl;
			if(i==Int_t(num*0.8)) cout<<"*******************************completed 80%"<<"**************************************"<<endl;
			if(i==Int_t(num*0.9)) cout<<"*******************************completed 90%"<<"**************************************"<<endl;
			if(i==Int_t(num*1-1)) cout<<"*******************************completed !!!"<<"**************************************"<<endl;

			saveToy(signal,m12,m13,eventNumber,i,700);
		}
	}

	UnbinnedDataSet data({m12, m13, eventNumber});
	TFile *ff = TFile::Open("data_tagAll.root");
	TTree *tree = (TTree *)ff->Get("DTag");
	assert(tree);

	fptype m_Kppi0_sq, m_Kspi0_sq, m_recoil, dE_sig;
	tree->SetBranchAddress("m_Kppi0_sq",&m_Kppi0_sq);
	tree->SetBranchAddress("m_Kspi0_sq",&m_Kspi0_sq);
	tree->SetBranchAddress("m_recoil",&m_recoil);
	tree->SetBranchAddress("dE_sig",&dE_sig);

	for(int i = 0; i < tree->GetEntries(); i++){
		tree->GetEvent(i);
		if(m_recoil<mrec_min||m_recoil>mrec_max) continue;
		if(dE_sig<dE_min||dE_sig>dE_max) continue;
		if(!cpuDalitz(m_Kppi0_sq,m_Kspi0_sq)) continue;
		m12.setValue(m_Kppi0_sq);
		m13.setValue(m_Kspi0_sq);
		eventNumber.setValue(data.getNumEvents());
		data.addEvent();
	}
	ff->Close();

	signal->setData(&data);
	signal->setDataSize(data.getNumEvents());

	signal->normalize();

	cout << fixed << setprecision(8);
	vector <vector<fptype>> fracMat;
	fracMat.clear();
	fracMat = signal->fit_fractions();
	const int num_res = fracMat.size();
	cout << "Un-diagonal elements in fraction matrix may be meaningless!!!" << endl;
	for (int i = 0; i < num_res; ++i){
		for (int j = 0; j < num_res; ++j){
			cout << fracMat[i][j] << ",  ";
		}
		cout << endl;
	}

//	ProdPdf prodpdf{"prodpdf", {signal}};
//	DalitzPlotter plotter(&prodpdf, signal);

//	UnbinnedDataSet toyMC({m12, m13, eventNumber});
//	plotter.fillDataSetMC(toyMC, 600000);

//	int oldBins1 = m12.getNumBins();
//	int oldBins2 = m13.getNumBins();
//	m12.setNumBins(drawBinM12);
//	m13.setNumBins(drawBinM13);


//	TH1F m12_sigpdf("m12_sigpdf", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
//	TH1F m13_sigpdf("m13_sigpdf", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
//	TH1F m23_sigpdf("m23_sigpdf", "", m12.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
//
//
//	fptype m12_tmp,m13_tmp,m23_tmp;
//	for (int i = 0; i < eventNumber.getValue(); ++i){
//		m12_tmp = toyMC.getValue(m12, i);
//		m12_sigpdf.Fill(m12_tmp);
//
//		m13_tmp = toyMC.getValue(m13, i);
//		m13_sigpdf.Fill(m13_tmp);
//
//		m23_tmp = cpuGetM23(m12_tmp, m13_tmp);
//		m23_sigpdf.Fill(m23_tmp);
//	}
//
//	//for out put hist
//	TCanvas *foo = new TCanvas("c1","c1",1200,350);
//	foo->Divide(3,1);
//
//	foo->cd(1);
//	m12_sigpdf.Draw("H");
//
//	foo->cd(2);
//	m13_sigpdf.Draw("H");
//
//	foo->cd(3);
//	m23_sigpdf.Draw("H");
//
//	foo->SaveAs("plots/test_pdf_plot.C");
//
//	m12.setNumBins(oldBins1);
//	m13.setNumBins(oldBins2);

}

void saveToy(Amp3Body *signal, Observable &m12, Observable &m13, EventNumber &eventNumber, int id, int num){
//	Observable m12("m12", 0.3, 2);
//	Observable m13("m13", 0.3, 2);
//	EventNumber eventNumber("eventNumber");
//	m12.setNumBins(400);
//	m13.setNumBins(400);

	ProdPdf prodpdf{"prodpdf", {signal}};
	DalitzPlotter plotter(&prodpdf, signal);
	
	UnbinnedDataSet toyMC({m12, m13, eventNumber});
	plotter.fillDataSetMC(toyMC, num);

	fptype m12_tmp,m13_tmp,m23_tmp;
	char text[100];
	sprintf(text,"./divide_root_toy_setB/toy_setB_%d.root",id);
	TFile *newfile = new TFile(text,"recreate");
	TTree *newtree = new TTree("DTag","");
	
	double m_recoil, dE_sig, m_KpKs_sq, m_Kppi0_sq, m_Kspi0_sq;
	newtree->Branch("m_recoil",&m_recoil,"m_recoil/D");
	newtree->Branch("dE_sig",&dE_sig,"dE_sig/D");
	newtree->Branch("m_Kppi0_sq",&m_Kppi0_sq,"m_Kppi0_sq/D");
	newtree->Branch("m_Kspi0_sq",&m_Kspi0_sq,"m_Kspi0_sq/D");
	newtree->Branch("m_KpKs_sq",&m_KpKs_sq,"m_KpKs_sq/D");

	for (int i = 0; i < eventNumber.getValue(); ++i){
		m12_tmp = toyMC.getValue(m12, i);
		m13_tmp = toyMC.getValue(m13, i);
		m23_tmp = cpuGetM23(m12_tmp, m13_tmp);

		m_recoil = 1.87;
		dE_sig = 0;
		m_Kppi0_sq = m12_tmp;
		m_Kspi0_sq = m13_tmp;
		m_KpKs_sq = m23_tmp;
		newtree->Fill();
	}
	newtree->Write();
	newfile->Close();
}


int main(int argc, char **argv) {
	GooFit::Application app("Dalitz example", argc, argv);
	srand((unsigned)time(NULL));
	int seed = 2333333*(double)rand()/((double)RAND_MAX);
	rnd.SetSeed(seed);

	app.add_option("-f,--filename,filename", data_filename, "File to read in", true)->check(GooFit::ExistingFile);

	GOOFIT_PARSE(app);

	GooFit::setROOTStyle();

	//Observables setup
	Observable m12("m12", 0.3, 2);
	Observable m13("m13", 0.3, 2);

	EventNumber eventNumber("eventNumber");
	m12.setNumBins(BinNumsM12);
	m13.setNumBins(BinNumsM13);

	//Prepare the data
	UnbinnedDataSet data({m12, m13, eventNumber});

	//Set up efficiency pdf
	GooPdf* eff = NULL;
	if(m_fit_use_eff){
		if(!m_effPoly)
			eff = makeHistogramPdf(data, eff_filename, "th2d_dalitz", "efficiency");
		else
			eff = fitEffPoly(data);
	}


if(m_fit_data){
	//Read in data
	getRealData(data_filename, app, data);

	//floating initial values
	vector <fptype> init_val;

	if(m_float_init){
		init_val.clear();
		for (int j = 0; j < 30; ++j)
			init_val.push_back(fRand(-5,5));

		Amp3Body *signal = makeSignalPdf(m12, m13, eventNumber, init_val, 0, eff, false);
		fptype nll = runDataFit(signal, &data, false);
		cout << "init values: " ;
		for (int k = 0; k < init_val.size(); ++k)
			cout << init_val[k] << "      ";
		cout << endl;
		cout << "NLL = " << endl << nll << endl;
	}
	else{
		//init value with min nll
		init_val.clear();
		if (!m_effPoly){
			//for smooth eff // data
			
			//Kbar(892)0
			init_val.push_back(-0.3903576304746);
			init_val.push_back(0.1297451247656);

			//K(1410)+ 
			init_val.push_back(-1.429793910011);
			init_val.push_back(0.5946039356012);

			//Kbar(1410)0
			init_val.push_back(-1.429793910011);
			init_val.push_back(0.5946039356012);

			//a_0(980)+ 
			init_val.push_back(0.86683);
			init_val.push_back(0.99428);

			//Kpi S-wave
			init_val.push_back(-1.543673413092);
			init_val.push_back(1.301284773652);

			init_val.push_back(-3.124202874666);
			init_val.push_back(-0.3446744759832);

			//nonr
			init_val.push_back(1.2332);
			init_val.push_back(3.5807);
		
			//K(1430)+
			init_val.push_back(-0.5067136006698);
			init_val.push_back(1.718314937409);

			//Kbar(1430)0
			init_val.push_back(-2.461852600875);
			init_val.push_back(2.451792324187);

			//a0(1450)+
			init_val.push_back(1.60814191448);
			init_val.push_back(-0.8384406300105);

			//for smooth eff // IO check
//			init_val.push_back(-0.249478426886);
//			init_val.push_back(-0.5757343593398);
//			init_val.push_back(0.9852732038613);
//			init_val.push_back(1.59076489141);
//			init_val.push_back(-0.960467580812);
//			init_val.push_back(2.119866925649);
//			init_val.push_back(-1.040647152845);
//			init_val.push_back(-1.192725433531);
//			init_val.push_back(-0.7938823596397);
//			init_val.push_back(2.620252723337);

		}
		else{
			//for poly eff // data
			init_val.push_back(-0.2228764621018);
			init_val.push_back(-0.5832810690462);
			init_val.push_back(0.8823176264095);
			init_val.push_back(1.637020214731);
			init_val.push_back(-1.059158275127);
			init_val.push_back(1.986955448399);
			init_val.push_back(-0.8918729448466);
			init_val.push_back(-1.194650555205);
			init_val.push_back(-0.8159530063214);
			init_val.push_back(2.625492017721);

			//for poly eff // IO check
//			init_val.push_back(-0.2975142854568);
//			init_val.push_back(-0.6277874717491);
//			init_val.push_back(0.6793008154279);
//			init_val.push_back(2.000863745694);
//			init_val.push_back(-1.050730354501);
//			init_val.push_back(2.329015298883);
//			init_val.push_back(-0.8597034302501);
//			init_val.push_back(-1.577886799133);
//			init_val.push_back(-1.732354196769);
//			init_val.push_back(2.915962244536);
		}

		Amp3Body *signal = NULL;
		if(m_fit_use_eff)
			signal = makeSignalPdf(m12, m13, eventNumber, init_val, 0, eff, false);
		else
			signal = makeSignalPdf(m12, m13, eventNumber, init_val, 0, 0, false);

		runDataFit(signal, &data, true);
	}
}
	if(m_test_pdf)
		//use_eff, save_toy
		gen_test_pdf(false,true);
	return 1;
}
