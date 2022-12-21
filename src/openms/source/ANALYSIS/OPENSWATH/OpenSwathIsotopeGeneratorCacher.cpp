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

#include  <OpenMS/ANALYSIS/OPENSWATH/OpenSwathIsotopeGeneratorCacher.h>
// theoretical spectrum generator
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

namespace OpenMS
{

  //
  OpenSwathIsotopeGeneratorCacher::~OpenSwathIsotopeGeneratorCacher() = default;
  // note these should be mass values rather than m/z values
  void OpenSwathIsotopeGeneratorCacher::initialize(double massStart, double massEnd, double massStep)
  {
    massStep_ = massStep;
    for (double mass=massStart; mass < massEnd; mass += massStep_)
    {
      std::cout << "Current mz: " << mass << std::endl;

      // create the theoretical distribution
      //Note: this is a rough estimate of the weight, usually the protons should be deducted first, left for backwards compatibility.
      IsotopeDistribution dist = solver_.estimateFromPeptideWeight(mass);

      cachedIsotopeDistributions_[mass] = dist;
    }
  }

  double OpenSwathIsotopeGeneratorCacher::getMassStep()
  {
    return massStep_;
  }

  double OpenSwathIsotopeGeneratorCacher::getHalfMassStep()
  {
    return halfMassStep_;
  }

  int OpenSwathIsotopeGeneratorCacher::getMaxIsotope()
  {
    return solver_.getMaxIsotope();
  }

  bool OpenSwathIsotopeGeneratorCacher::getRoundMasses()
  {
    return solver_.getRoundMasses();
  }


  IsotopeDistribution OpenSwathIsotopeGeneratorCacher::addEntry(double mass)
  {
    IsotopeDistribution dist = solver_.estimateFromPeptideWeight(mass);

    cachedIsotopeDistributions_[mass] = dist;
    return dist;
  }

  IsotopeDistribution OpenSwathIsotopeGeneratorCacher::get(double mass)
  {
    auto ele = cachedIsotopeDistributions_.upper_bound(mass);
    //auto ele = std::upper_bound(cachedIsotopeDistributions_.begin(), cachedIsotopeDistributions_.end(), mass);


    if (ele == cachedIsotopeDistributions_.end())
    {
      // element is not found so create an entry
      return addEntry(mass);
    }
    else if (std::abs(ele->first - mass) > halfMassStep_) // element found is more than halfMassStep_ away from target mass so generate a new entry
    {
      return addEntry(mass);
    }
    else // there is a current cached distribution that can be used!
    {
      auto ptrPrevious = ele--;

      if (std::abs(ptrPrevious->first - mass) < std::abs(ele->first - mass)) // check if the previous entry is closer than the first entry
      {
        return ptrPrevious->second;
      }
      else
      {
        return ele->second;
      }
    }

  }


  std::vector<std::pair<double, double>> OpenSwathIsotopeGeneratorCacher::get(double product_mz, int charge, const double mannmass)
  {
    IsotopeDistribution distribution = get(product_mz * charge);

    double currentIsotopeMz = product_mz; // mz value of current isotope

    std::vector<std::pair<double, double> > isotopes_spec;
    for (IsotopeDistribution::Iterator it = distribution.begin(); it != distribution.end(); ++it)
    {
      isotopes_spec.emplace_back(currentIsotopeMz, it->getIntensity());
      currentIsotopeMz += mannmass / charge;
    }

    return isotopes_spec;
  }
}
