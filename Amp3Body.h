#pragma once

#include <goofit/PDFs/physics/Amp3BodyBase.h>
#include <goofit/PDFs/physics/DalitzPlotHelpers.h>

#include <goofit/detail/Complex.h>

namespace GooFit {

class SpecialResonanceIntegrator;
class SpecialResonanceCalculator;
class DalitzPlotter;

/**
A time-independent description of the Dalitz plot
as a coherent sum of resonances:

\f[
P(m^2_{12},m^2_{13};\vec\alpha) = \left|\sum\limits_i
\alpha_i B_i(m^2_{12},m^2_{13})\right|^2
\epsilon(m^2_{12},m^2_{13})
\f]

where \f$\alpha_i\f$ is a complex coefficient, \f$B_i\f$ is a resonance
parametrization (see Goofit::ResonancePdf), and \f$\epsilon\f$ is a
real-valued efficiency function. The constructor takes the
squared-mass variables \f$m_{12}\f$ and \f$m_{13}\f$, an event index (this
is used in caching), a GooFit::DecayInfo3 object which contains a `vector`
of GooFit::ResonancePdf's as well as some global information like the
mother and daughter masses, and the efficiency function.
**/

class Amp3Body : public Amp3BodyBase {
  public:
    Amp3Body(std::string n,
             Observable m12,
             Observable m13,
             EventNumber eventNumber,
             DecayInfo3 decay,
             GooPdf *eff = nullptr);
    // Note that 'efficiency' refers to anything which depends on (m12, m13) and multiplies the
    // coherent sum. The caching method requires that it be done this way or the ProdPdf
    // normalization will get *really* confused and give wrong answers.

    __host__ fptype normalize() override;
    __host__ void setDataSize(unsigned int dataSize, unsigned int evtSize = 3);
    __host__ void setForceIntegrals(bool f = true) { forceRedoIntegrals = f; }

    __host__ void populateArrays() override;

    /// Get the cached wave (device) vectors
    __host__ const thrust::device_vector<fpcomplex> &getCachedWaveNoCopy(size_t i) const { return *(cachedWaves[i]); }

    __host__ const std::vector<std::complex<fptype>> getCachedWave(size_t i) const;

    /// Sum up a cached wave
    __host__ fpcomplex sumCachedWave(size_t i) const;

    /// Get the decay info struct
    __host__ DecayInfo3 &getDecayInfo() { return decayInfo; }

    /// Calculate fit fractions (Cache should be pre-filled)
    __host__ std::vector<std::vector<fptype>> fit_fractions();

    friend DalitzPlotter;

  protected:
    DecayInfo3 decayInfo;
    Observable _m12;
    Observable _m13;
    EventNumber _eventNumber;
    fptype *dalitzNormRange;

    // Following variables are useful if masses and widths, involved in difficult BW calculation,
    // change infrequently while amplitudes, only used in adding BW results together, change rapidly.
    thrust::device_vector<fpcomplex> *cachedWaves[16]; // Caches the BW values for each event.
    fpcomplex ***integrals; // Caches the integrals of the BW waves for each combination of resonances.

    bool *redoIntegral;
    mutable bool forceRedoIntegrals;
    fptype *cachedMasses;
    fptype *cachedWidths;
    int totalEventSize;
    int cacheToUse;
    SpecialResonanceIntegrator ***integrators;
    SpecialResonanceCalculator **calculators;

    //added by xxh for fit fractions
    SpecialResonanceIntegrator ***integrators_ff;
    fpcomplex ***integrals_ff; 

    unsigned int efficiencyFunction;
};

} // namespace GooFit
