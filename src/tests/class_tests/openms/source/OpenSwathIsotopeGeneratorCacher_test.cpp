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
// $Maintainer: Joshua Charkow $
// $Authors: Joshua Charkow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathIsotopeGeneratorCacher.h>
#include <OpenMS/CONCEPT/ClassTest.h>


// data access
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

// For creating isotopic distributions
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

using namespace std;
using namespace OpenMS;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_TEST(OpenSwathIsotopeGeneratorCacher, "$Id$")


START_SECTION(OpenSwathIsotopeGeneratorCacher(Size max_isotope, double massStep))
{
    IsotopeDistribution* nullPointer = nullptr;
    OpenSwathIsotopeGeneratorCacher* ptr = new OpenSwathIsotopeGeneratorCacher(0, 0.5);
    Size max_isotope = ptr->getMaxIsotope();
    TEST_EQUAL(max_isotope, 0)
    TEST_EQUAL(ptr->getRoundMasses(), false)
    TEST_EQUAL(ptr->getMassStep(), 0.5)
    TEST_EQUAL(ptr->getHalfMassStep(), 0.25)
    TEST_NOT_EQUAL(ptr, nullPointer)
    delete ptr;
}
END_SECTION

START_SECTION(void OpenSwathIsotopeGeneratorCacher::initialize(double massStart, double massEnd, double massStep))
{
  int maxIsotope = 2;
  double massStep = 1;
  OpenSwathIsotopeGeneratorCacher isotopeCacher(maxIsotope, massStep);
  isotopeCacher.initialize(100, 102, massStep); // since non inclusive massEnd should generate cache for 100 and 101

  std::map<double, IsotopeDistribution> cache = isotopeCacher.fetchCache();

  CoarseIsotopePatternGenerator isotopeGenerator(maxIsotope, massStep);

  // create the test cache
  std::map<double, IsotopeDistribution> test;
  test[100] = isotopeGenerator.estimateFromPeptideWeight(100);
  test[101] = isotopeGenerator.estimateFromPeptideWeight(101);


  TEST_EQUAL(cache.size(), test.size()); // both should be 2
  TEST_EQUAL(cache.size(), 2);

  // test isotopes distributions are equal
  for ( const auto &[key, testDist]: test )
  {
    TEST_EQUAL(cache[key].getContainer().size(), test[key].getContainer().size())
    for (Size i = 0; i != testDist.size(); ++i)
    {
      TEST_EQUAL(round(cache[key].getContainer()[i].getMZ()), test[key].getContainer()[i].getMZ())
      TEST_REAL_SIMILAR(cache[key].getContainer()[i].getIntensity(), test[key].getContainer()[i].getIntensity())
    }
  }
}
END_SECTION

START_SECTION(IsotopeDistribution OpenSwathIsotopeGeneratorCacher::get(double mass))
{
  // Setup
  int maxIsotope = 2;
  double massStep = 5;
  IsotopeDistribution cachedDistribution;
  IsotopeDistribution testDistribution;
  OpenSwathIsotopeGeneratorCacher isotopeCacher(maxIsotope, massStep);
  isotopeCacher.initialize(100, 120, massStep); // since non inclusive massEnd should generate cache for 100 and 101
  CoarseIsotopePatternGenerator isotopeGenerator(maxIsotope, massStep);


  // ############ Case #1 get(100) Should fetch distribution of mass 100 because it is in the spectrum ###############
  std::cout << "Case 1" << std::endl;
  cachedDistribution = isotopeCacher.get(100);
  testDistribution = isotopeGenerator.estimateFromPeptideWeight(100);

  TEST_EQUAL(testDistribution.size(), cachedDistribution.size())
  for (Size i = 0; i != testDistribution.size(); ++i)
  {
    TEST_EQUAL(round(cachedDistribution.getContainer()[i].getMZ()), testDistribution.getContainer()[i].getMZ())
    TEST_REAL_SIMILAR(cachedDistribution.getContainer()[i].getIntensity(), testDistribution.getContainer()[i].getIntensity())
  }

  // ########## Case #2 get(107.5) Should fetch distribution of mass 110 because it is the closest to the requested mass and it is less than halfStep away #################
  std::cout << "Case 2" << " get(" << 107.5 <<  ")" << std::endl;
  cachedDistribution = isotopeCacher.get(107.5);
  testDistribution = isotopeGenerator.estimateFromPeptideWeight(110);

  TEST_EQUAL(testDistribution.size(), cachedDistribution.size())
  for (Size i = 0; i != testDistribution.size(); ++i)
  {
    TEST_EQUAL(round(cachedDistribution.getContainer()[i].getMZ()), testDistribution.getContainer()[i].getMZ())
    TEST_REAL_SIMILAR(cachedDistribution.getContainer()[i].getIntensity(), testDistribution.getContainer()[i].getIntensity())
  }


  // ########## Case #3 get(103) Should fetch distribution of mass 105 because it is the closest to the requested mass (always round up if halfway of step) and it is less than halfStep away #################

  std::cout << "Case 3" << std::endl;
  cachedDistribution = isotopeCacher.get(103);
  testDistribution = isotopeGenerator.estimateFromPeptideWeight(105);

  TEST_EQUAL(testDistribution.size(), cachedDistribution.size())
  for (Size i = 0; i != testDistribution.size(); ++i)
  {
    TEST_EQUAL(round(cachedDistribution.getContainer()[i].getMZ()), testDistribution.getContainer()[i].getMZ())
    TEST_REAL_SIMILAR(cachedDistribution.getContainer()[i].getIntensity(), testDistribution.getContainer()[i].getIntensity())
  }

  // ########## Case #4 get(200.4) Should fetch distribution of mass 200.4 because there is no cached distribution close to requested mass
  std::cout << "Case 4" << std::endl;
  cachedDistribution = isotopeCacher.get(200.4);
  testDistribution = isotopeGenerator.estimateFromPeptideWeight(200.4); // don't need to recreate this

  TEST_EQUAL(testDistribution.size(), cachedDistribution.size())
  for (Size i = 0; i != testDistribution.size(); ++i)
  {
    TEST_EQUAL(round(cachedDistribution.getContainer()[i].getMZ()), testDistribution.getContainer()[i].getMZ())
    TEST_REAL_SIMILAR(cachedDistribution.getContainer()[i].getIntensity(), testDistribution.getContainer()[i].getIntensity())
  }

  // ########## Case #5 get(200.7) Should fetch distribution of mass 200.4 because this is cached (from example above) and within range
  std::cout << "Case 5" << std::endl;
  cachedDistribution = isotopeCacher.get(201);
  // testDistribution = isotopeGenerator.estimateFromPeptideWeight(200.4); // don't need to recreate this

  TEST_EQUAL(testDistribution.size(), cachedDistribution.size())
  for (Size i = 0; i != testDistribution.size(); ++i)
  {
    TEST_EQUAL(round(cachedDistribution.getContainer()[i].getMZ()), testDistribution.getContainer()[i].getMZ())
    TEST_REAL_SIMILAR(cachedDistribution.getContainer()[i].getIntensity(), testDistribution.getContainer()[i].getIntensity())
  }

  // ########## Case #6 get(90.5) Should fetch distribution of mass 90.5 because no value is close to this
  std::cout << "Case 6" << std::endl;
  cachedDistribution = isotopeCacher.get(90.5);
  testDistribution = isotopeGenerator.estimateFromPeptideWeight(90.5); // don't need to recreate this

  TEST_EQUAL(testDistribution.size(), cachedDistribution.size())
  for (Size i = 0; i != testDistribution.size(); ++i)
  {
    TEST_EQUAL(round(cachedDistribution.getContainer()[i].getMZ()), testDistribution.getContainer()[i].getMZ())
    TEST_REAL_SIMILAR(cachedDistribution.getContainer()[i].getIntensity(), testDistribution.getContainer()[i].getIntensity())
  }

  // ########## Case #7 get(89) Should fetch distribution of mass 90.5 because this is close
  std::cout << "Case 6" << std::endl;
  cachedDistribution = isotopeCacher.get(89);
  testDistribution = isotopeGenerator.estimateFromPeptideWeight(90.5); // don't need to recreate this

  TEST_EQUAL(testDistribution.size(), cachedDistribution.size())
  for (Size i = 0; i != testDistribution.size(); ++i)
  {
    TEST_EQUAL(round(cachedDistribution.getContainer()[i].getMZ()), testDistribution.getContainer()[i].getMZ())
    TEST_REAL_SIMILAR(cachedDistribution.getContainer()[i].getIntensity(), testDistribution.getContainer()[i].getIntensity())
  }





}

END_SECTION




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
