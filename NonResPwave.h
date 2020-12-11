#pragma once

#include <goofit/PDFs/physics/resonances/Resonance.h>

namespace GooFit {

namespace Resonances {

/// Nonresonant constructor
class NonResPwave : public ResonancePdf {
  public:
    NonResPwave(std::string name, Variable ar, Variable ai, unsigned int cyc);
    ~NonResPwave() override = default;
};

} // namespace Resonances

} // namespace GooFit
