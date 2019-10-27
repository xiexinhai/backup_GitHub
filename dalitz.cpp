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

char strbuffer[1000];
const fptype _mDp       = 1.86965;
const fptype _mDp2      = _mDp * _mDp;
const fptype _mDp2inv   = 1. / _mDp2;
const fptype KsMass = 0.497611;
const fptype KpMass = 0.493677;
const fptype pi0Mass = 0.134977;
const int BinNumsM12 = 100;
const int BinNumsM13 = 100;
const int BinNumsM23 = 100;
const fptype m12_lower = (KpMass+pi0Mass)*(KpMass+pi0Mass);
const fptype m12_upper = (_mDp-KsMass)*(_mDp-KsMass);

const fptype m13_lower = (KsMass+pi0Mass)*(KsMass+pi0Mass);
const fptype m13_upper = (_mDp-KpMass)*(_mDp-KpMass);

const fptype m23_lower = (KsMass+KpMass)*(KsMass+KpMass);
const fptype m23_upper = (_mDp-pi0Mass)*(_mDp-pi0Mass);

bool m_draw_data = false;
bool m_effPoly = true;

bool m_float_init = false;
bool m_float_polyeff = false;

bool m_draw_polyeff = true;
bool m_draw_smootheff = false;

TRandom3 rnd;

void makeEffPlot(GooPdf* total, UnbinnedDataSet* data);
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
bool cpuDalitz (fptype m_12, fptype m_13, fptype bigM = _mDp, fptype dm1 = pi0Mass, fptype dm2 = KpMass, fptype dm3 = KsMass) {
//  if (m_12 > m_13) return false; // Only the upper corner is considered
  if (m_12 < pow(dm1 + dm2, 2)) return false; // This m_12 cannot exist, it's less than the square of the (1,2) particle mass.
  if (m_12 > pow(bigM - dm3, 2)) return false;   // This doesn't work either, there's no room for an at-rest 3 daughter. 
  
  // Calculate energies of 1 and 3 particles in m_12 rest frame. 
  fptype e1star = 0.5 * (m_12 - dm2*dm2 + dm1*dm1) / sqrt(m_12); 
  fptype e3star = 0.5 * (bigM*bigM - m_12 - dm3*dm3) / sqrt(m_12); 

  // Bounds for m_13 at this value of m_12.
  fptype minimum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - dm1*dm1) + sqrt(e3star*e3star - dm3*dm3), 2);
  if (m_13 < minimum) return false;
  fptype maximum = pow(e1star + e3star, 2) - pow(sqrt(e1star*e1star - dm1*dm1) - sqrt(e3star*e3star - dm3*dm3), 2);
  if (m_13 > maximum) return false;

  return true; 
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

    TH2F dalitzplot("dalitzplot",
                    "",
                    m12.getNumBins(),
                    m12.getLowerLimit(),
                    m12.getUpperLimit(),
                    m13.getNumBins(),
                    m13.getLowerLimit(),
                    m13.getUpperLimit());
	std::vector<Observable> vars;
	vars.push_back(m12);
	vars.push_back(m13);
	vars.push_back(eventNumber);

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
		if(m_recoil<1.8648||m_recoil>1.8772) continue;
		if(dE_sig<-0.03||dE_sig>0.02) continue;
		if(!cpuDalitz(m_Kppi0_sq,m_Kspi0_sq)) continue;
		m12.setValue(m_Kppi0_sq);
		m13.setValue(m_Kspi0_sq);
		eventNumber.setValue(data.getNumEvents());
		data.addEvent();
		dalitzplot.Fill(m12.getValue(), m13.getValue());
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
    dtop0pp.meson_radius = 1.5;

//	bool fixAmps = false; // Takes ~400x longer

	//K*(892)+
	Variable K892pMass("K892p_mass", 0.89176, 0.00025, 0.85, 0.93);
	Variable K892pWidth("K892p_width", 0.0503, 0.0008, 0.04, 0.06);
	ResonancePdf *K892p = new Resonances::RBW(
		"K892p", Variable("K892p_amp_real", 1), Variable("K892p_amp_imag", 0), K892pMass, K892pWidth, 1, PAIR_12);

	//K*(1410)+
	Variable K1410pMass("K1410p_mass", 1.421, 0.009, 1.2, 1.6);
	Variable K1410pWidth("K1410p_width", 0.236, 0.018, 0.1, 0.4);
	ResonancePdf *K1410p = new Resonances::RBW(
		"K1410p",
		fixAmps ? Variable("K1410p_amp_real", 0) : Variable("K1410p_amp_real", 0, 0.01, 0, 0),
		fixAmps ? Variable("K1410p_amp_imag", 0) : Variable("K1410p_amp_imag", 0, 0.01, 0, 0),
		K1410pMass,
		K1410pWidth,
		1,
		PAIR_12);

	//non resonance
	ResonancePdf *nonr = new Resonances::NonRes(
		"nonr",
		fixAmps ? Variable("nonr_amp_real", 16.82603957812) : Variable("nonr_amp_real", -1, 0.01, 0, 0),
		fixAmps ? Variable("nonr_amp_imag", 11.59309345464) : Variable("nonr_amp_imag", -1, 0.01, 0, 0));

	//Kbar*(892)0
	Variable K892zeroMass("K892zero_mass", 0.89555, 0.00020, 0.85, 0.93);
	Variable K892zeroWidth("K892zero_width", 0.0473, 0.0005, 0.04, 0.06);

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

	//Kbar*(1410)0
	Variable K1410zeroMass("K1410zero_mass", 1.421, 0.009, 1.2, 1.6);
	Variable K1410zeroWidth("K1410zero_width", 0.236, 0.018, 0.1, 0.4);

	sprintf(strbuffer, "K1410zero_amp_real_%d", i);
	Variable K1410zero_amp_real(strbuffer,init_val[2], 0.01, 0, 0);
	sprintf(strbuffer, "K1410zero_amp_imag_%d", i);
	Variable K1410zero_amp_imag(strbuffer,init_val[3], 0.01, 0, 0);
	ResonancePdf *K1410zero = new Resonances::RBW(
		"K1410zero",
		fixAmps ? Variable("K1410zero_amp_real", 1.24683082311) : K1410zero_amp_real,
		fixAmps ? Variable("K1410zero_amp_imag", 0.4483709871808) : K1410zero_amp_imag,
		K1410zeroMass,
		K1410zeroWidth,
		1,
		PAIR_13);

	//Kbar*(1430)0_2
	Variable K1430zero2Mass("K1430zero2_mass", 1.4324, 0.0013, 1.2, 1.6);
	Variable K1430zero2Width("K1430zero2_width", 0.109, 0.005, 0.1, 0.4);
	ResonancePdf *K1430zero2 = new Resonances::RBW(
		"K1430zero2",
		fixAmps ? Variable("K1430zero2_amp_real", 1.714) : Variable("K1430zero2_amp_real", 1.714, 0.01, 0, 0),
		fixAmps ? Variable("K1430zero2_amp_imag", -0.125) : Variable("K1430zero2_amp_imag", -0.125, 0.01, 0, 0),
		K1430zero2Mass,
		K1430zero2Width,
		2,
		PAIR_13);

	//a980+
	Variable a980pMass("a980pMass",0.980,0.020,0.9,1.06);
	ResonancePdf *a980p = new Resonances::FLATTE(
		"a980p",
		fixAmps ? Variable("a980p_amp_real", 80) : Variable("a980p_amp_real", 80, 0.01, 0, 0),
		fixAmps ? Variable("a980p_amp_imag", -25) : Variable("a980p_amp_imag", -25, 0.01, 0, 0),
		a980pMass,
//		Variable("g1",0.324,0.015,0,0),
//		Variable("rg2og1",1.03,0.00014,0,0),
		Variable("g1",3.24,0.015,0,0),
		Variable("rg2og1",1.03,0.00014,0,0),
		PAIR_23,
		false);

	//rho(1450)+
	Variable rho1450pMass("rho1450p_Mass", 1.465, 0.025, 1.3, 1.6);
	Variable rho1450pWidth("rho1450p_Width", 0.4, 0.06, 0.04, 0.06);
	ResonancePdf *rho1450p = new Resonances::GS(
		"rho1450p",
		fixAmps ? Variable("rho1450p_amp_real", -1.040242409497) : Variable("rho1450p_amp_real", -1.215, 0.01, 0, 0),
		fixAmps ? Variable("rho1450p_amp_imag", 1.452729026326) : Variable("rho1450p_amp_imag", -9.707, 0.01, 0, 0),
		rho1450pMass,
		rho1450pWidth,
		1,
		PAIR_23);

	//rho(1700)+
	Variable rho1700pMass("rho1700p_Mass", 1.541, 0.035, 1.4, 1.7);
	Variable rho1700pWidth("rho1700p_Width", 0.25, 0.1, 0.1, 0.4);
	ResonancePdf *rho1700p = new Resonances::GS(
		"rho1700p",
		fixAmps ? Variable("rho1700p_amp_real", -1.241) : Variable("rho1700p_amp_real", -1.241, 0.01, 0, 0),
		fixAmps ? Variable("rho1700p_amp_imag", 6.049) : Variable("rho1700p_amp_imag", 6.049, 0.01, 0, 0),
		rho1700pMass,
		rho1700pWidth,
		1,
		PAIR_23);

	//a0(1450)+
	Variable a1450pMass("a1450p_Mass", 1.474, 0.019, 1.2, 1.6);
	Variable a1450pWidth("a1450p_Width", 0.265, 0.013, 0.1, 0.3);

	sprintf(strbuffer, "a1450p_amp_real_%d", i);
	Variable a1450p_amp_real(strbuffer,init_val[4], 0.01, 0, 0);
	sprintf(strbuffer, "a1450p_amp_imag_%d", i);
	Variable a1450p_amp_imag(strbuffer,init_val[5], 0.01, 0, 0);
	ResonancePdf *a1450p = new Resonances::RBW(
		"a1450p",
		fixAmps ? Variable("a1450p_amp_real", -1.040242409497) : a1450p_amp_real,
		fixAmps ? Variable("a1450p_amp_imag", 1.452729026326) : a1450p_amp_imag,
		a1450pMass,
		a1450pWidth,
		0,
		PAIR_23);

	//SwaveKppi0
	Variable SwaveKppi0Mass("SwaveKppi0_Mass", 1.425, 0.050, 1.2, 1.6);
	Variable SwaveKppi0Width("SwaveKppi0_Width", 0.27, 0.08, 0.1, 0.3);

	sprintf(strbuffer, "SwaveKppi0_amp_real_%d", i);
	Variable SwaveKppi0_amp_real(strbuffer,init_val[6], 0.01, 0, 0);
	sprintf(strbuffer, "SwaveKppi0_amp_imag_%d", i);
	Variable SwaveKppi0_amp_imag(strbuffer,init_val[7], 0.01, 0, 0);
	ResonancePdf *SwaveKppi0 = new Resonances::LASS(
		"SwaveKppi0",
		fixAmps ? Variable("SwaveKppi0_amp_real", 1.869715774108) : SwaveKppi0_amp_real,
		fixAmps ? Variable("SwaveKppi0_amp_imag", 0.2625631972586) : SwaveKppi0_amp_imag,
		SwaveKppi0Mass,
		SwaveKppi0Width,
		0,
		PAIR_12);

	//SwaveKspi0
	Variable SwaveKspi0Mass("SwaveKspi0_Mass", 1.425, 0.050, 1.2, 1.6);
	Variable SwaveKspi0Width("SwaveKspi0_Width", 0.27, 0.08, 0.1, 0.3);

	sprintf(strbuffer, "SwaveKspi0_amp_real_%d", i);
	Variable SwaveKspi0_amp_real(strbuffer,init_val[8], 0.01, 0, 0);
	sprintf(strbuffer, "SwaveKspi0_amp_imag_%d", i);
	Variable SwaveKspi0_amp_imag(strbuffer,init_val[9], 0.01, 0, 0);
	ResonancePdf *SwaveKspi0 = new Resonances::LASS(
		"SwaveKspi0",
		fixAmps ? Variable("SwaveKspi0_amp_real", -2.61217416547) : SwaveKspi0_amp_real,
		fixAmps ? Variable("SwaveKspi0_amp_imag", -0.126650709178) : SwaveKspi0_amp_imag,
		SwaveKspi0Mass,
		SwaveKspi0Width,
		0,
		PAIR_13);

	dtop0pp.resonances.push_back(K892p);
//	dtop0pp.resonances.push_back(K1410p);
//	dtop0pp.resonances.push_back(nonr);
	dtop0pp.resonances.push_back(K892zero);
	dtop0pp.resonances.push_back(K1410zero);

//	dtop0pp.resonances.push_back(K1430zero2);

//	dtop0pp.resonances.push_back(a980p);

//	dtop0pp.resonances.push_back(rho1450p);
//	dtop0pp.resonances.push_back(rho1700p);

	dtop0pp.resonances.push_back(a1450p);
	dtop0pp.resonances.push_back(SwaveKppi0);
	dtop0pp.resonances.push_back(SwaveKspi0);

    bool fitMasses = false;

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
		int ix = i/m13.getNumBins();
		int iy = i%m13.getNumBins();
		if (ix > iy){
			int itt = iy; iy = ix; ix = itt;
		}
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
	SmoothHistogramPdf* ret = new SmoothHistogramPdf(pdfname.c_str(), binEffData, Variable("smoothConst", 0)); 

	cout << "end generate SmoothHistogramPdf" << endl;

	m12.setNumBins(oldBins1);
	m13.setNumBins(oldBins2);

	f->Close();
 	cout << "end generate efficiency pdf" << endl;
	if(m_draw_smootheff)
		makeEffPlot(ret,&data);

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
	Variable conYmax("conYmax", init_val[12]);//, 0.001, 0, 1);

	Variable decZmax("decZmax", 1.52031, 0.001, 0, 5);
	Variable conZmax("conZmax", 0.41866, 0.001, 0, 1);
	
	Variable maxDalitzX("maxDalitzX", pow(_mDp - KsMass, 2));
	TrigThresholdPdf* hiX = new TrigThresholdPdf("hiX", m12, maxDalitzX, decXmax, conXmax, true); 
	
	Variable maxDalitzY("maxDalitzY", pow(_mDp - KpMass, 2));
	TrigThresholdPdf* hiY = new TrigThresholdPdf("hiY", m13, maxDalitzY, decYmax, conYmax, true); 

//	Variable maxDalitzZ("maxDalitzZ", pow(_mDp - pi0Mass, 2));
//	TrigThresholdPdf* hiZ = new TrigThresholdPdf("hiZ", m12, m13, maxDalitzZ, decZmax, conZmax, massSum, true); 

	Variable decXmin("decXmin", 6.22596, 0.001, 0, 50);
	Variable conXmin("conXmin", 0.65621, 0.001, 0, 1);

	Variable decYmin("decYmin", 6.30722, 0.001, 0, 50);
	Variable conYmin("conYmin", 0.69527, 0.001, 0, 1);

	Variable minDalitzX("minDalitzX", pow(KpMass + pi0Mass, 2));
	TrigThresholdPdf* loX = new TrigThresholdPdf("loX", m12, minDalitzX, decXmin, conXmin, false);

	Variable minDalitzY("minDalitzY", pow(KsMass + pi0Mass, 2));
	TrigThresholdPdf* loY = new TrigThresholdPdf("loY", m13, minDalitzY, decXmin, conXmin, false);

//	Variable* decZmin = new Variable("decZmin",10.82390, 0.001, 0, 50);

//	Variable* conZmin = new Variable("conZmin", 0.31764, 0.001, 0, 1);


//	TrigThresholdPdf* loZ = new TrigThresholdPdf("loZ", m12, m13, minDalitzX, decZmin, conZmin, massSum, false);

	vector<PdfBase*> comps;
	comps.clear();
	comps.push_back(poly);
	comps.push_back(hiX);
	comps.push_back(hiY);
//	comps.push_back(hiZ);

//	comps.push_back(loX);
//	comps.push_back(loY);

//	comps.push_back(loZ);
	comps.push_back(getDPVeto(m12,m13)); 
	ProdPdf* ret = new ProdPdf("efficiency_total", comps); 
	return ret; 
}

UnbinnedDataSet *loadEffData(UnbinnedDataSet &data,std::string toyFileName = "efficiency_acceptance.root") {
	auto obs               = data.getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

//	EventNumber eventNumber("eventNumber");
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
	effpdfplot->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
	effpdfplot->GetYaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");

	TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	m12_pdf_hist.SetStats(false); 
	m12_pdf_hist.SetLineColor(kRed); 
	m12_pdf_hist.SetLineWidth(2); 

	TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m13_pdf_hist.SetStats(false); 
	m13_pdf_hist.SetLineColor(kRed); 
	m13_pdf_hist.SetLineWidth(2); 

	TH1F m23_pdf_hist("m23_pdf_hist", "", 10*m12.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
//	TH1F m23_pdf_hist("m23_pdf_hist", "", 2*m12.getNumBins(), m23_lower, m23_upper);
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
	effplot->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
	effplot->GetYaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");


	TH1F *m12_data = new TH1F("m12_data","", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	m12_data->SetMarkerStyle(20);
	m12_data->SetMarkerSize(0.6);
	m12_data->SetLineWidth(2);
	m12_data->GetYaxis()->SetTitle("Efficiency");
	m12_data->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
	m12_data->SetLineColor(1);
	
	TH1F *m13_data = new TH1F("m13_data","", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	m13_data->SetMarkerStyle(20);
	m13_data->SetMarkerSize(0.6);
	m13_data->SetLineWidth(2);
	m13_data->GetYaxis()->SetTitle("Efficiency");
	m13_data->GetXaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");
	m13_data->SetLineColor(1);
	
	TH1F *m23_data = new TH1F("m23_data","", 10*m12.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
//	TH1F *m23_data = new TH1F("m23_data","", m12.getNumBins(), m23_lower, m23_upper);
	m23_data->SetMarkerStyle(20);
	m23_data->SetMarkerSize(0.6);
	m23_data->SetLineWidth(2);
	m23_data->GetYaxis()->SetTitle("Efficiency");
	m23_data->GetXaxis()->SetTitle("M^{2}(K^{+}K_{S}^{0}) (GeV^{2}/c^{4})");
	m23_data->SetLineColor(1);


	TChain *tr_data = new TChain("DTag");tr_data->Add("efficiency_acceptance.root");
	tr_data->Project("m12_data","m_Kppi0_sq");
	tr_data->Project("m13_data","m_Kspi0_sq");
	tr_data->Project("m23_data","m_KpKs_sq");
	tr_data->Draw("m_Kspi0_sq:m_Kppi0_sq>>effplot");

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
	pull->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
	pull->GetYaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");

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
}

GooPdf* fitEffPoly(UnbinnedDataSet &data){
	auto obs               = data.getObservables();
	Observable m12         = obs.at(0);
	Observable m13         = obs.at(1);
	Observable eventNumber = obs.at(2);

//	Observable m12("m12",pow(pi0Mass+KpMass,2),pow(_mDp-KsMass,2));
//	Observable m13("m13",pow(pi0Mass+KsMass,2),pow(_mDp-KpMass,2));

	GooPdf* eff = NULL;
	UnbinnedDataSet *effdata = loadEffData(data,"efficiency_acceptance.root");
	vector <fptype> init_val;

	if(m_float_polyeff){
		init_val.clear();
		for (int j = 0; j < 9; ++j)
			init_val.push_back(fRand(-1,1));
		init_val.push_back(fRand(0,5));
		init_val.push_back(fRand(0,1));
		init_val.push_back(fRand(0,5));
		init_val.push_back(fRand(0,1));
		
	}
	else{
		init_val.clear();
		init_val.push_back(-0.09619442260917);
		init_val.push_back(0.1050472227158);
		init_val.push_back(0.0879306537382);
		init_val.push_back(0.05862313927887);
		init_val.push_back(0.02635542505541);
		init_val.push_back(0.0859644179791);
		init_val.push_back(-0.2162739116849);
		init_val.push_back(0.08941190852959);
		init_val.push_back(-0.2544651299809);

		init_val.push_back(1.809830652448);
		init_val.push_back(0.6856680565195);
		init_val.push_back(2.294445970101);
		init_val.push_back(0);
	}
	eff = makeEfficiencyPoly(m12,m13,init_val);
	eff->setData(effdata);
//	eff->setDataSize(effdata->getNumEvents());

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
		makeEffPlot(eff,effdata);
	return eff;
}

double runDataFit(Amp3Body *signal, UnbinnedDataSet *data, bool m_err, bool fixAmps = false) {
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

	//fit fractions
	vector <vector<fptype>> fracMat;
	fracMat = signal->fit_fractions();
	vector <fptype> fracList;
	fracList.clear();

	const int num_res = fracMat.size();
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
		datapdf->printCovMat();
		ProdPdf prodpdf{"prodpdf", {signal}};
	
		DalitzPlotter plotter(&prodpdf, signal);
		datapdf->printCovMat();
 
		TH2F *dalitzplot = plotter.make2D();
		dalitzplot->SetTitle("dalitzFit");
		dalitzplot->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
		dalitzplot->GetYaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");
	
		TH1F m12_pdf_hist("m12_pdf_hist", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
		m12_pdf_hist.SetStats(false); 
		m12_pdf_hist.SetLineColor(kRed); 
		m12_pdf_hist.SetLineWidth(2); 
	
		TH1F m13_pdf_hist("m13_pdf_hist", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
		m13_pdf_hist.SetStats(false); 
		m13_pdf_hist.SetLineColor(kRed); 
		m13_pdf_hist.SetLineWidth(2); 
	
		TH1F m23_pdf_hist("m23_pdf_hist", "", m12.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
		m23_pdf_hist.SetStats(false); 
		m23_pdf_hist.SetLineColor(kRed); 
		m23_pdf_hist.SetLineWidth(2); 
	
//		//for out put hist
//		TH1F m12_sigpdf("m12_sigpdf", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
//		TH1F m13_sigpdf("m13_sigpdf", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
//		TH1F m23_sigpdf("m23_sigpdf", "", m12.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));

		fptype m12_tmp,m13_tmp,m23_tmp;
	
		UnbinnedDataSet toyMC({m12, m13, eventNumber});
		plotter.fillDataSetMC(toyMC, 5000000);

		cout << "toyMC size = " << eventNumber.getValue() << endl;
		for (int i = 0; i < eventNumber.getValue(); ++i){
			m12_tmp = toyMC.getValue(m12, i);
			m12_pdf_hist.Fill(m12_tmp);
//			m12_sigpdf.Fill(m12_tmp);
	
			m13_tmp = toyMC.getValue(m13, i);
			m13_pdf_hist.Fill(m13_tmp);
//			m13_sigpdf.Fill(m13_tmp);
	
			m23_tmp = cpuGetM23(m12_tmp, m13_tmp);
			m23_pdf_hist.Fill(m23_tmp);
//			m23_sigpdf.Fill(m23_tmp);
		}

//		//for out put hist
//		m12_sigpdf.Scale(100000/m12_sigpdf.Integral());
//		m13_sigpdf.Scale(100000/m13_sigpdf.Integral());
//		m23_sigpdf.Scale(100000/m23_sigpdf.Integral());
//
//		m12_sigpdf.SaveAs("plots/m12_sigpdf.root");
//		m13_sigpdf.SaveAs("plots/m13_sigpdf.root");
//		m23_sigpdf.SaveAs("plots/m23_sigpdf.root");
	
		//for comparison
	
		TString a("Events/"); TString c(" GeV^{2}/c^{4})");
		TH1F *m12_data = new TH1F("m12_data","", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
		char b_m12[20];  sprintf(b_m12, "(%g",(m12.getUpperLimit()-m12.getLowerLimit())/m12.getNumBins());
		TString ytitle_m12 = a + b_m12 + c;
		m12_data->SetMarkerStyle(20);
		m12_data->SetMarkerSize(0.6);
		m12_data->SetLineWidth(2);
		m12_data->GetYaxis()->SetTitle(ytitle_m12);
		m12_data->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
		m12_data->SetLineColor(1);
	
	
		TH1F *m13_data = new TH1F("m13_data","", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
		char b_m13[20];  sprintf(b_m13, "(%g",(m13.getUpperLimit()-m13.getLowerLimit())/m13.getNumBins());
		TString ytitle_m13 = a + b_m13 + c;
		m13_data->SetMarkerStyle(20);
		m13_data->SetMarkerSize(0.6);
		m13_data->SetLineWidth(2);
		m13_data->GetYaxis()->SetTitle(ytitle_m13);
		m13_data->GetXaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");
		m13_data->SetLineColor(1);
	
		TH1F *m23_data = new TH1F("m23_data","", m12.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));
		char b_m23[20];  sprintf(b_m23, "(%g",(getM23UpperLimit(m12,m13)-getM23LowerLimit(m12,m13))/m12.getNumBins());
		TString ytitle_m23 = a + b_m23 + c;
		m23_data->SetMarkerStyle(20);
		m23_data->SetMarkerSize(0.6);
		m23_data->SetLineWidth(2);
		m23_data->GetYaxis()->SetTitle(ytitle_m13);
		m23_data->GetXaxis()->SetTitle("M^{2}(K^{+}K_{S}^{0}) (GeV^{2}/c^{4})");
		m23_data->SetLineColor(1);
	
		//Drawing 2D pull distribution
		TH2F *dalitzData = new TH2F ("dalitzData","dalitzData",
			dalitzplot->GetNbinsX(), dalitzplot->GetXaxis()->GetXmin(), dalitzplot->GetXaxis()->GetXmax(),
			dalitzplot->GetNbinsY(), dalitzplot->GetYaxis()->GetXmin(), dalitzplot->GetYaxis()->GetXmax());
	
	
		if(!fixAmps){
			TChain *tr_data = new TChain("DTag");tr_data->Add("data_tagAll.root");
			TCut cut_mbc = "((m_recoil>1.8648)&&(m_recoil<1.8772))";
			TCut cut_dE = "((dE_sig>-0.03)&&(dE_sig<0.02))";
	
			tr_data->Project("m12_data","m_Kppi0_sq",cut_mbc&&cut_dE);
			tr_data->Project("m13_data","m_Kspi0_sq",cut_mbc&&cut_dE);
			tr_data->Project("m23_data","m_KpKs_sq",cut_mbc&&cut_dE);
	
			tr_data->Draw("m_Kspi0_sq:m_Kppi0_sq>>dalitzData",cut_mbc&&cut_dE);
			cout << dalitzData->GetEntries() << endl;
		}
	
		TCanvas tmp("c2","c2",800,700);
		tmp.Divide(2,2);
		tmp.cd(1);
		gStyle->SetPalette(52,0);
		dalitzData->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
		dalitzData->GetYaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");
		dalitzData->Draw("colz");
	
		tmp.cd(2);
		dalitzplot->Scale(dalitzData->Integral()/dalitzplot->Integral());
		dalitzplot->Draw("colz");
	
		TH2F *pull = new TH2F("Pull","Pull Distribution",
			dalitzplot->GetNbinsX(), dalitzplot->GetXaxis()->GetXmin(), dalitzplot->GetXaxis()->GetXmax(),
			dalitzplot->GetNbinsY(), dalitzplot->GetYaxis()->GetXmin(), dalitzplot->GetYaxis()->GetXmax());
		pull->GetXaxis()->SetTitle("M^{2}(K^{+}#pi^{0}) (GeV^{2}/c^{4})");
		pull->GetYaxis()->SetTitle("M^{2}(K_{S}^{0}#pi^{0}) (GeV^{2}/c^{4})");
	
		dalitzplot->Rebin2D(5,5);
		dalitzData->Rebin2D(5,5);
		pull->Rebin2D(5,5);
	
		int nonEmpty = 0;
		double Chi2 = 0;
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
		cout << "chi and chi2 non-empty bin nums: " << nonEmpty << endl;
		cout << "chi2 = " << Chi2 << endl;
		tmp.cd(3);
		pull->Draw("colz");
//		tmp.SaveAs("plots/Original_Pull.C");

	
		//Drawing data & fit comparison
		TCanvas *foo = new TCanvas("c1","c1",800,700);
		foo->Divide(2,2);
	
		foo->cd(1);
		pull->Draw("colz");
	
		foo->cd(2);
		m12_data->Draw("e");
		m12_pdf_hist.Scale(m12_data->Integral()/m12_pdf_hist.Integral());
		m12_pdf_hist.Draw("Hsame");
//		m12_data->Rebin(2);
//		m12_pdf_hist.Rebin(2);

		foo->cd(3);
		m13_data->Draw("e");
		m13_pdf_hist.Scale(m13_data->Integral()/m13_pdf_hist.Integral());
		m13_pdf_hist.Draw("Hsame");
//		m13_data->Rebin(2);
//		m13_pdf_hist.Rebin(2);
	
		foo->cd(4);
		m23_data->Draw("e");
		m23_pdf_hist.Scale(m23_data->Integral()/m23_pdf_hist.Integral());
		m23_pdf_hist.Draw("Hsame");
//		m23_data->Rebin(2);
//		m23_pdf_hist.Rebin(2);

		if(!m_effPoly)
			foo->SaveAs("plots/dalitz_with_projections.C");
		else
			foo->SaveAs("plots/dalitz_with_projections_effPoly.C");
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
	
	//	fraction_rate *= 0.3323*0.5/0.3330;
	//	fraction_rate_err *= 0.3323*0.5/0.3330;
	
		froot->Close();
		std::cout << std::endl;
		for (int ii=0;ii<num_res;ii++) 
			std::cout<<"Fit fraction for #"<<ii<<": "<<fracList[ii]<< " +- "<<errList[ii]<<std::endl;
		std::cout<<"Fit fraction for K*(892)+ / Kbar*(892)0 rate: " << fraction_rate << " +- " << fraction_rate_err << std::endl;
	}

	delete datapdf;
	return nll;
//	return 0;
}

void gen_test_pdf(){
	Observable m12("m12", 0.3, 2);
	Observable m13("m13", 0.3, 2);
	EventNumber eventNumber("eventNumber");
	m12.setNumBins(BinNumsM12);
	m13.setNumBins(BinNumsM13);

	TH1F m12_sigpdf("m12_sigpdf", "", m12.getNumBins(), m12.getLowerLimit(), m12.getUpperLimit());
	TH1F m13_sigpdf("m13_sigpdf", "", m13.getNumBins(), m13.getLowerLimit(), m13.getUpperLimit());
	TH1F m23_sigpdf("m23_sigpdf", "", m12.getNumBins(), getM23LowerLimit(m12,m13), getM23UpperLimit(m12,m13));

	DecayInfo3 dtop0pp;
	dtop0pp.motherMass   = _mDp;
	dtop0pp.daug1Mass    = pi0Mass;
	dtop0pp.daug2Mass    = KpMass;
	dtop0pp.daug3Mass    = KsMass;
	dtop0pp.meson_radius = 1.5;

	Variable K892pMass("K892p_mass", 0.89176);
	Variable K892pWidth("K892p_width", 0.0503);
	ResonancePdf *K892p = new Resonances::RBW(
		"K892p", 
		Variable("K892p_amp_real", 1), 
		Variable("K892p_amp_imag", 0), 
		K892pMass, 
		K892pWidth, 
		1, 
		PAIR_12);

	Variable K892zeroMass("K892zero_mass", 0.89555);
	Variable K892zeroWidth("K892zero_width", 0.0473);
	ResonancePdf *K892zero = new Resonances::RBW(
		"K892zero",
		Variable("K892zero_amp_real", -0.4871474223254),
		Variable("K892zero_amp_imag", 0.001120132013781),
		K892zeroMass,
		K892zeroWidth,
		1,
		PAIR_13);

     Variable K1410zeroMass("K1410zero_mass", 1.421);
     Variable K1410zeroWidth("K1410zero_width", 0.236);
     ResonancePdf *K1410zero = new Resonances::RBW(
          "K1410zero",
          Variable("K1410zero_amp_real", -1.676234687673),
          Variable("K1410zero_amp_imag", 1.25884833851),
          K1410zeroMass,
          K1410zeroWidth,
          1,
          PAIR_13);

     Variable a1450pMass("a1450p_Mass", 1.474);
     Variable a1450pWidth("a1450p_Width", 0.265);
     ResonancePdf *a1450p = new Resonances::RBW(
          "a1450p",
          Variable("a1450p_amp_real", -1.167432844979),
          Variable("a1450p_amp_imag", 0.2174462280103),
          a1450pMass,
          a1450pWidth,
          0,
          PAIR_23);

     Variable SwaveKppi0Mass("SwaveKppi0_Mass", 1.425);
     Variable SwaveKppi0Width("SwaveKppi0_Width", 0.27);
     ResonancePdf *SwaveKppi0 = new Resonances::LASS(
          "SwaveKppi0",
          Variable("SwaveKppi0_amp_real", 5.251775650343),
          Variable("SwaveKppi0_amp_imag", 0.4439884037173),
          SwaveKppi0Mass,
          SwaveKppi0Width,
          0,
          PAIR_12);

     Variable SwaveKspi0Mass("SwaveKspi0_Mass", 1.425);
     Variable SwaveKspi0Width("SwaveKspi0_Width", 0.27);
     ResonancePdf *SwaveKspi0 = new Resonances::LASS(
          "SwaveKspi0",
          Variable("SwaveKspi0_amp_real", 0.5333076345728),
          Variable("SwaveKspi0_amp_imag", -1.568683551211),
          SwaveKspi0Mass,
          SwaveKspi0Width,
          0,
          PAIR_13);


	dtop0pp.resonances.push_back(K892p);
	dtop0pp.resonances.push_back(K892zero);
     dtop0pp.resonances.push_back(K1410zero);
     dtop0pp.resonances.push_back(a1450p);
     dtop0pp.resonances.push_back(SwaveKppi0);
     dtop0pp.resonances.push_back(SwaveKspi0);

	vector<Variable> offsets;
	vector<Observable> observables;
	vector<Variable> coefficients;
	observables.push_back(m12);
	observables.push_back(m13);
	offsets.push_back(constantZero);
	offsets.push_back(constantZero);
	coefficients.push_back(constantOne);
	GooPdf *eff = new PolynomialPdf("constantEff", observables, coefficients, offsets, 0);

	Amp3Body *signal = new Amp3Body("signalPDF_test", m12, m13, eventNumber, dtop0pp, eff);

	ProdPdf prodpdf{"prodpdf", {signal}};
	DalitzPlotter plotter(&prodpdf, signal);
	
	UnbinnedDataSet toyMC({m12, m13, eventNumber});
	plotter.fillDataSetMC(toyMC, 600000);

	fptype m12_tmp,m13_tmp,m23_tmp;
	cout << "toyMC size = " << eventNumber.getValue() << endl;

	for (int i = 0; i < eventNumber.getValue(); ++i){
		m12_tmp = toyMC.getValue(m12, i);
		m12_sigpdf.Fill(m12_tmp);

		m13_tmp = toyMC.getValue(m13, i);
		m13_sigpdf.Fill(m13_tmp);

		m23_tmp = cpuGetM23(m12_tmp, m13_tmp);
		m23_sigpdf.Fill(m23_tmp);
	}

	//for out put hist
//	m12_sigpdf.Scale(200000/m12_sigpdf.Integral());
//	m13_sigpdf.Scale(200000/m13_sigpdf.Integral());
//	m23_sigpdf.Scale(200000/m23_sigpdf.Integral());

	TCanvas *foo = new TCanvas("c1","c1",1200,350);
	foo->Divide(3,1);

	foo->cd(1);
	m12_sigpdf.Draw("");

	foo->cd(2);
	m13_sigpdf.Draw("");

	foo->cd(3);
	m23_sigpdf.Draw("");

	foo->SaveAs("plots/test_pdf_plot.C");
}

int main(int argc, char **argv) {
	GooFit::Application app("Dalitz example", argc, argv);
	srand((unsigned)time(NULL));
	int seed = 2333333*(double)rand()/((double)RAND_MAX);
	rnd.SetSeed(seed);

	std::string filename = "data_tagAll.root";
	app.add_option("-f,--filename,filename", filename, "File to read in", true)->check(GooFit::ExistingFile);

	GOOFIT_PARSE(app);

	GooFit::setROOTStyle();


    // Observables setup
	Observable m12("m12", 0.3, 2);
	Observable m13("m13", 0.3, 2);

	EventNumber eventNumber("eventNumber");
	m12.setNumBins(BinNumsM12);
	m13.setNumBins(BinNumsM13);

	// Prepare the data
	UnbinnedDataSet data({m12, m13, eventNumber});

//	Set up efficiency pdf
	GooPdf* eff = NULL;
	if(!m_effPoly)
		eff = makeHistogramPdf(data, "efficiency_acceptance.root", "th2d_dalitz", "efficiency");
	else
		eff = fitEffPoly(data);

//	Read in data
	getRealData(filename, app, data);

	//floating initial values
	vector <fptype> init_val;

	if(m_float_init){
		init_val.clear();
		for (int j = 0; j < 10; ++j)
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
			//for smooth eff
			init_val.push_back(-0.4871474953126);
			init_val.push_back(0.001120166841566);
			init_val.push_back(-1.676234590844);
			init_val.push_back(1.258848242572);
			init_val.push_back(-1.167432382186);
			init_val.push_back(0.2174461981026);
			init_val.push_back(5.251776861255);
			init_val.push_back(0.4439887755764);
			init_val.push_back(0.5333087459871);
			init_val.push_back(-1.568682267362);
		}
		else{
			//for poly eff
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
		}
		Amp3Body *signal = makeSignalPdf(m12, m13, eventNumber, init_val, 0, eff, false);
		runDataFit(signal, &data, true);
	}

//	init_val.clear();
//	init_val.push_back(-0.1567496650418);
//	init_val.push_back(-0.6644198096373);
//	init_val.push_back(1.24683082311);
//	init_val.push_back(0.4483709871808);
//	init_val.push_back(-1.040242409497);
//	init_val.push_back(1.452729026326);
//	init_val.push_back(1.86971577410);
//	init_val.push_back(0.2625631972586);
//	init_val.push_back(-2.61217416547);
//	init_val.push_back(-0.126650709178);
//	Amp3Body *signal = makeSignalPdf(m12, m13, eventNumber, init_val, 0, eff, false);
//	runDataFit(signal, &data, true);

	gen_test_pdf();
	return 1;
}
