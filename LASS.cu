#include <goofit/PDFs/physics/resonances/LASS.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/lineshapes/Lineshape.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ fpcomplex lass(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    unsigned int spin         = pc.getConstant(0);
    unsigned int cyclic_index = pc.getConstant(1);

    fptype resmass  = pc.getParameter(0);
    fptype reswidth = pc.getParameter(1);

    fptype rMassSq  = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
    fptype frFactor = 1;

    resmass *= resmass;
    // Calculate momentum of the two daughters in the resonance rest frame; note symmetry under interchange (dm1 <->
    // dm2).

//    fptype measureDaughterMoms = twoBodyCMmom(rMassSq,
//                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
//                                              (PAIR_23 == cyclic_index ? c_daug3Mass : c_daug2Mass));
//    fptype nominalDaughterMoms = twoBodyCMmom(resmass,
//                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
//                                              (PAIR_23 == cyclic_index ? c_daug3Mass : c_daug2Mass));

    fptype measureDaughterMoms = twoBodyCMmom(rMassSq,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_23 == cyclic_index ? c_daug3Mass : (PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass)));
    fptype nominalDaughterMoms = twoBodyCMmom(resmass,
                                              (PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass),
                                              (PAIR_23 == cyclic_index ? c_daug3Mass : (PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass)));

    if(0 != spin) {
        frFactor = dampingFactorSquare(nominalDaughterMoms, spin, c_meson_radius);
        frFactor /= dampingFactorSquare(measureDaughterMoms, spin, c_meson_radius);
    }

    // Implement LASS:
    /*
    fptype s = kinematics(m12, m13, _trackinfo[i]);
    fptype q = twoBodyCMmom(s, _trackinfo[i]);
    fptype m0  = _massRes[i]->getValFast();
    fptype _g0 = _gammaRes[i]->getValFast();
    int spin   = _spinRes[i];
    fptype g = runningWidthFast(s, m0, _g0, spin, _trackinfo[i], FrEval(s, m0, _trackinfo[i], spin));
    */

    fptype q = measureDaughterMoms;
    fptype g = reswidth * pow(measureDaughterMoms / nominalDaughterMoms, 2.0 * spin + 1) * frFactor / sqrt(rMassSq);
    g *= sqrt(resmass);

	//old
//    fptype _a    = 0.22357;
//    fptype _r    = -15.042;
//    fptype _R    = 1;
//    fptype _phiR = 1.10644;
//    fptype _B    = 0.614463;
//    fptype _phiB = -0.0981907;

	//new cited PhysRevD.98.112012
	fptype fac = 3.141592653/180.0;

	fptype _a    = 0.113;
//	fptype _a    = 0.113+0.006;
//	fptype _a    = 0.113-0.006;

	fptype _r    = -33.8;
//	fptype _r    = -33.8+1.8;
//	fptype _r    = -33.8-1.8;

	fptype _R    = 1;

	fptype _phiR = -109.7;
//	fptype _phiR = -109.7+2.6;
//	fptype _phiR = -109.7-2.6;
	_phiR *= fac;

	fptype _B    = 0.96;
//	fptype _B    = 0.96+0.07;
//	fptype _B    = 0.96-0.07;

	fptype _phiB = 0.1;
//	fptype _phiB = 0.1+0.3;
//	fptype _phiB = 0.1-0.3;
	_phiB *= fac;

/*
	//old cited PhysRevLett.105.081803
	fptype _a    = 0.224;
//	fptype _a    = 0.224+0.003;
//	fptype _a    = 0.224-0.003;

	fptype _r    = -15.01;
//	fptype _r    = -15.01+0.13;
//	fptype _r    = -15.01-0.13;

	fptype _R    = 1;

	fptype _phiR = 1.10;
//	fptype _phiR = 1.10+0.02;
//	fptype _phiR = 1.10-0.02;

	fptype _B    = 0.62;
//	fptype _B    = 0.62+0.04;
//	fptype _B    = 0.62-0.04;

	fptype _phiB = -0.100;
//	fptype _phiB = -0.100+0.01;
//	fptype _phiB = -0.100-0.01;

*/

    // background phase motion
    fptype cot_deltaB  = (1.0 / (_a * q)) + 0.5 * _r * q;
    fptype qcot_deltaB = (1.0 / _a) + 0.5 * _r * q * q;

    // calculate resonant part
    fpcomplex expi2deltaB = fpcomplex(qcot_deltaB, q) / fpcomplex(qcot_deltaB, -q);
    fpcomplex resT        = fpcomplex(cos(_phiR + 2 * _phiB), sin(_phiR + 2 * _phiB)) * _R;

    fpcomplex prop = fpcomplex(1, 0) / fpcomplex(resmass - rMassSq, -sqrt(resmass) * g);
    // resT *= prop*m0*_g0*m0/twoBodyCMmom(m0*m0, _trackinfo[i])*expi2deltaB;
    resT *= prop * (resmass * reswidth / nominalDaughterMoms) * expi2deltaB;

    // calculate bkg part
    resT += fpcomplex(cos(_phiB), sin(_phiB)) * _B * (cos(_phiB) + cot_deltaB * sin(_phiB)) * sqrt(rMassSq)
            / fpcomplex(qcot_deltaB, -q);

    resT *= sqrt(frFactor);
    resT *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);

    pc.incrementIndex(1, 2, 2, 0, 1);

    return resT;
}

__device__ resonance_function_ptr ptr_to_LASS = lass;

namespace Resonances {

LASS::LASS(std::string name, Variable ar, Variable ai, Variable mass, Variable width, unsigned int sp, unsigned int cyc)
    : ResonancePdf("LASS", name, ar, ai) {
    registerParameter(mass);
    registerParameter(width);

    registerConstant(sp);
    registerConstant(cyc);

    registerFunction("ptr_to_LASS", ptr_to_LASS);
}

} // namespace Resonances
} // namespace GooFit
