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
#include <OpenMS/CONCEPT/LogStream.h>

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

  std::map<double, IsotopeDistribution> OpenSwathIsotopeGeneratorCacher::fetchCache()
  {
    return cachedIsotopeDistributions_;
  }


  IsotopeDistribution OpenSwathIsotopeGeneratorCacher::addEntry(double mass)
  {
    IsotopeDistribution dist = solver_.estimateFromPeptideWeight(mass);

    cachedIsotopeDistributions_[mass] = dist;
    return dist;
  }

  IsotopeDistribution OpenSwathIsotopeGeneratorCacher::computeMiss(double mass) const
  {
    OPENMS_LOG_DEBUG << "Cache miss :(" << std::endl;
    CoarseIsotopePatternGenerator tempSolver = solver_; //copy the isotope generator to enable const signature
    return tempSolver.estimateFromPeptideWeight(mass);
  }

  IsotopeDistribution OpenSwathIsotopeGeneratorCacher::getImmutable(double mass) const
  {
    /*
    // print the cachedIsotopeDistribution list
    for (const auto &[key, value]:cachedIsotopeDistributions_)
    {
      //std::cout << "key is: " << key << std::endl;
    }
    */

    //std::cout << "JOSH start get immutable" << std::endl;
    if (mass <= (cachedIsotopeDistributions_.begin()->first + halfMassStep_) ) // all cached greater than target (with tolerance)
    {
      //std::cout << "All cached elements strictly greater than target (with tolerance)" << std::endl;
      if ( (cachedIsotopeDistributions_.begin()->first - mass) <= halfMassStep_)
      {
        //std::cout << "first element matches";
        return cachedIsotopeDistributions_.begin()->second;
      }
      else
      {
        //std::cout << "first element does not match";
        return computeMiss(mass);
      }
    }
    else // mass is somewhere in the middle or all elements are less than mass
    {
      auto upperBound = cachedIsotopeDistributions_.upper_bound(mass);
      auto prevEle = upperBound;
      prevEle--;
      //std::cout << "ptr previous: " << prevEle->first << std::endl;
      //std::cout << "upper bound (first element not less than " << mass << ": " << upperBound->first << std::endl;

      if (upperBound == cachedIsotopeDistributions_.end())
      {
        //std::cout << "upper bound is at the end, meaning all elements are less than or equal to mass" << std::endl;
        if (mass - prevEle->first <= halfMassStep_)
        {
          //std::cout << "last element is in range" << std::endl;
          return prevEle->second;
        }
        else
        {
          //std::cout << "last element is not in range" << std::endl;
          return computeMiss(mass);
        }
      }
      else if ((mass - prevEle->first) > halfMassStep_)
      {
        //std::cout << "mass" << mass << std::endl;
        //std::cout << "Previous element is too far away (half mass step: " << (mass - prevEle->first) << ")"  << std::endl;
        if ((upperBound->first - mass) <= halfMassStep_)
        {

          //std::cout << "upper bound in range, match!" << std::endl;
          return upperBound->second;
        }
        else
        {
          //std::cout << "upper bound not in range, new entry" << std::endl;
          return computeMiss(mass);
        }
      }
      else //prevEle in range
      {
        //std::cout << "previous element is in range" << std::endl;
        if ((upperBound->first - mass) > halfMassStep_)
        {
          //std::cout << "upper bound not in range match previous element" << std::endl;
          return prevEle->second;
        }
        else // both elements in range
        {
          //std::cout << "both elements in range see which is closer" << std::endl;
          if ((upperBound->first - mass) <= (mass - prevEle->first))
          {
            //std::cout << "upper bound is closer" << std::endl;
            return upperBound->second;
          }
          else
          {
            //std::cout << "previous element is closer" << std::endl;
            return prevEle->second;
          }
        }
      }
    }
  }

  IsotopeDistribution OpenSwathIsotopeGeneratorCacher::get(double mass)
  {
    for (const auto &[key, value]:cachedIsotopeDistributions_)
    {
      //std::cout << "key is: " << key << std::endl;
    }

    if (mass <= (cachedIsotopeDistributions_.begin()->first + halfMassStep_) ) // all cached greater than target (with tolerance)
    {
      //std::cout << "All cached elements strictly greater than target (with tolerance)" << std::endl;
      if ( (cachedIsotopeDistributions_.begin()->first - mass) <= halfMassStep_)
      {
        //std::cout << "first element matches";
        return cachedIsotopeDistributions_.begin()->second;
      }
      else
      {
        //std::cout << "first element does not match";
        return addEntry(mass);
      }
    }
    else // mass is somewhere in the middle or all elements are less than mass
    {
      auto upperBound = cachedIsotopeDistributions_.upper_bound(mass);
      auto prevEle = upperBound;
      prevEle--;
      //std::cout << "ptr previous: " << prevEle->first << std::endl;
      //std::cout << "upper bound (first element not less than " << mass << ": " << upperBound->first << std::endl;

      if (upperBound == cachedIsotopeDistributions_.end())
      {
        //std::cout << "upper bound is at the end, meaning all elements are less than or equal to mass" << std::endl;
        if (mass - prevEle->first <= halfMassStep_)
        {
          //std::cout << "last element is in range" << std::endl;
          return prevEle->second;
        }
        else
        {
          //std::cout << "last element is not in range" << std::endl;
          return addEntry(mass);
        }
      }
      else if ((mass - prevEle->first) > halfMassStep_)
      {
        //std::cout << "mass" << mass << std::endl;
        //std::cout << "Previous element is too far away (half mass step: " << (mass - prevEle->first) << ")"  << std::endl;
        if ((upperBound->first - mass) <= halfMassStep_)
        {

          //std::cout << "upper bound in range, match!" << std::endl;
          return upperBound->second;
        }
        else
        {
          //std::cout << "upper bound not in range, new entry" << std::endl;
          return addEntry(mass);
        }
      }
      else //prevEle in range
      {
        //std::cout << "previous element is in range" << std::endl;
        if ((upperBound->first - mass) > halfMassStep_)
        {
          //std::cout << "upper bound not in range match previous element" << std::endl;
        }
        else // both elements in range
        {
          //std::cout << "both elements in range see which is closer" << std::endl;
          if ((upperBound->first - mass) <= (mass - prevEle->first))
          {
            //std::cout << "upper bound is closer" << std::endl;
            return upperBound->second;
          }
          else
          {
            //std::cout << "previous element is closer" << std::endl;
            return prevEle->second;
          }
        }
      }
    }
  }




/*
    if ((ele == cachedIsotopeDistributions_.end()) & (previousEle == cachedIsotopeDistributions_.end()) | (ele == cachedIsotopeDistributions_.begin()) // all cached elements are greater than target mass
    {
      //std::cout << "All elements greater than target mass, creating entry" << std::endl;
      return addEnrtry(mass);
    }
    else if ((ele == cachedIsotopeDistributions_.end()) & (previousEle != cachedIsotopeDistributions_.end())) // last element might match
    {
      //std::cout << "Last element might match" << std::endl;
      if (std::abs(previousEle->first - mass) < halfStep_)
      {
        //std::cout << "Last element matches" << std::endl;
        return previousEle->second;
      }
      else
      {
        return addEntry(mass);
      }
    }
    //else if ((ele == cachedIsotopeDistributions_.begin()) & (previousEle != cachedIsotopeDistributions_.end())) // last element might match
                                                                             //
      // at the end of the array, possible that we still have a match though if this is just slightly over the last element
      //std::cout << "At the end of the array..." << std::endl;
      if (std::abs(previousEle->first - mass) < std::abs(ele->first - mass)) // check if the previous entry is within range
      {
          double rslt = std::abs(previousEle->first - mass) < std::abs(ele->first - mass);
          //std::cout << "rslt is: " << rslt << std::endl;
          //std::cout << "Taking previous element" << std::endl;
          return ele->second;
      }
      else // previous entry not within range
      {
        double rslt = std::abs(previousEle->first - mass) < std::abs(ele->first - mass);
        //std::cout << "rslt is: " << rslt << std::endl;
        //std::cout << "No match creating entry" << std::endl;
        return addEntry(mass);
      }

    // determine if the current element or the element before it is closer
    if (std::abs(previousEle->first - mass) < std::abs(ele->first - mass)) // check if the previous entry is closer than the first entry
    {
        double rslt = std::abs(previousEle->first - mass) < std::abs(ele->first - mass);
        //std::cout << "rslt is: " << rslt << std::endl;
        //std::cout << "Taking previous element" << std::endl;
        ele = previousEle;
    }

    if (std::abs(ele->first - mass) > halfMassStep_) // element found is more than halfMassStep_ away from target mass so generate a new entry
    {
      double rslt = std::abs(ele->first - mass) > halfMassStep_;

      //std::cout << "Closest element is not in range, creating entry" << " closest element: " << ele->first << " delta: " << rslt << std::endl;
      return addEntry(mass);
    }
    else // there is a current cached distribution that can be used!
    {
      //std::cout << "Using the cache" << std::endl;

      return ele->second;
    }
    }
    */

  std::vector<std::pair<double, double>> OpenSwathIsotopeGeneratorCacher::getImmutable(double product_mz, int charge, const double mannmass) const
  {
    //std::cout << "Calling sub get immutable" << std::endl;
    IsotopeDistribution distribution = getImmutable(product_mz * charge);
    //std::cout << "done sub get immutable" << std::endl;
    //std::cout << "fetching distribution size" << std::endl;
    //std::cout << "distribition size is: " << distribution.size() << std::endl;
    double currentIsotopeMz = product_mz; // mz value of current isotope

    std::vector<std::pair<double, double> > isotopes_spec;
    for (IsotopeDistribution::Iterator it = distribution.begin(); it != distribution.end(); ++it)
    {
      isotopes_spec.emplace_back(currentIsotopeMz, it->getIntensity());
      currentIsotopeMz += mannmass / charge;
    }

    return isotopes_spec;
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
