#include <goofit/PDFs/physics/resonances/NonResPwave.h>

#include <goofit/PDFs/ParameterContainer.h>
#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

__device__ fpcomplex nonres_pwave(fptype m12, fptype m13, fptype m23, ParameterContainer &pc) {
    unsigned int cyclic_index = pc.getConstant(0);
    unsigned int spin         = 1;

//    fptype rMassSq    = (PAIR_12 == cyclic_index ? m12 : (PAIR_13 == cyclic_index ? m13 : m23));
//    fptype rMass = sqrt(rMassSq);
//    fptype mass_daug1 = PAIR_23 == cyclic_index ? c_daug2Mass : c_daug1Mass;
//    fptype mass_daug2 = PAIR_12 == cyclic_index ? c_daug2Mass : c_daug3Mass;
//    fptype mass_daug3 = PAIR_12 == cyclic_index ? c_daug3Mass : (PAIR_23 == cyclic_index?c_daug1Mass:c_daug2Mass);

    fpcomplex result{0.0, 0.0};
    fpcomplex ret(1.0, 0.0);
    ret *= spinFactor(spin, c_motherMass, c_daug1Mass, c_daug2Mass, c_daug3Mass, m12, m13, m23, cyclic_index);

    result += ret;
    pc.incrementIndex(1, 0, 1, 0, 1);
    return result;
}

__device__ resonance_function_ptr ptr_to_NONRES_PWAVE = nonres_pwave;

namespace Resonances {

NonResPwave::NonResPwave(std::string name, Variable ar, Variable ai, unsigned int cyc)
    : ResonancePdf("NonResPwave", name, ar, ai) {
    registerConstant(cyc);
    registerFunction("ptr_to_NONRES_PWAVE", ptr_to_NONRES_PWAVE);
}

} // namespace Resonances
} // namespace GooFit
