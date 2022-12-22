// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Joshua Charkow$
// $Authors: Joshua Charkow$
// --------------------------------------------------------------------------

#pragma once

// data access
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <vector>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

//#include <boost/shared_ptr.hpp>
//#include <boost/make_shared.hpp>

namespace OpenMS
{
  /** @brief A class that stores a map of theoretical spectra to be used in OpenSwathScoring
   *
   *
   *
  */
  class OPENMS_DLLAPI OpenSwathIsotopeGeneratorCacher
  {
    typedef OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;



    double massStep_; // step between adjacent isotope distributions
    double halfMassStep_; // cached value of half of the mass step
    CoarseIsotopePatternGenerator solver_; // used for computing the isotope distributions
    std::map<double, IsotopeDistribution> cachedIsotopeDistributions_;

  public:

    OpenSwathIsotopeGeneratorCacher(const Size max_isotope, const double massStep, const bool round_masses = false):
      massStep_(massStep),
      halfMassStep_(massStep / 2.0),
      solver_(max_isotope, round_masses),
      cachedIsotopeDistributions_()
    {
    };

    /// Destructor
    ~OpenSwathIsotopeGeneratorCacher();


    /** @brief Initialize cacher object
     * the scoring object
     *
     * Sets the parameters for the scoring.
     *
     * @param massStart - start mz to generate a theoretical isotope distribtuion
     * @param massEnd - end mz to generate a theoretical isotope distribtuion (non inclusive)
     * @param massStep - step between precursor mz for isotope distribution
     *
    */
    void initialize(double massStart, double massEnd, double massStep);

    double getMassStep();

    double getHalfMassStep();

    int getMaxIsotope();

    bool getRoundMasses();


    /** @brief fetches a copy of the cachedIsotopeDistribution map, useful for testing
     *
     */
    std::map<double, IsotopeDistribution> fetchCache();



    /** @brief compute and cache and isotope distribution with the corresponding mass
     * @param mass - mass to create the isotope distribution from
     */
    IsotopeDistribution addEntry(double mass);

    /** @breif gets the cached isotope distribution for the provided mass
     * if no isotope distribution exist within the halfMassStep_ then create and cache a new distribution
     *
     * @param mass - mass to fetch the isotope distribution
     */
    IsotopeDistribution get(double mass);

    /** @brief gets the theoretical spectrum corresponding with the provided m/z value
     * mz - mz of distribution to get
     * charge - charge of mz
     * mannmass - ???
     */
    //std::vector<std::pair<double, double>>  get(double mz, int charge, const double mannmass = 1.00048) const;

    /** @brief computes a missed cache but does not add corresponding mass to cache
     */
    IsotopeDistribution computeMiss(double mass) const;

    /** @brief gets form the cached isotope distribution for the provided mass
     *
     * if no isotope distribution exist within the halfMassStep_ then compute on the spot without caching
     */
    IsotopeDistribution getImmutable(double mass) const;

    /** @brief gets from cached isotope distribution for provided m/z and charge
     *  if no isotope distribution exists within the halfMassStep_ then create and cache a new distribution
     */
    std::vector<std::pair<double, double>> get(double mz, int charge, const double mannmass = 1.00048);


    /** @brief gets from cached isotope distribution for provided m/z and charge
     *  if no isotope distribution exists within the halfMassStep_ then create (but do not cache) a new distribution
     */
    std::vector<std::pair<double, double>> getImmutable(double mz, int charge, const double mannmass = 1.00048) const;
  };
}
