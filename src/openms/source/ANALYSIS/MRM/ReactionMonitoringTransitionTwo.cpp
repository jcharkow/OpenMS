// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransitionTwo.h>

#include <OpenMS/CONCEPT/Helpers.h>

#include <utility>

namespace OpenMS
{

  static const unsigned char DETECTING_TRANSITION_LOC = 0;
  static const unsigned char IDENTIFYING_TRANSITION_LOC = 1;
  static const unsigned char QUANTIFYING_TRANSITION_LOC = 2;

  ReactionMonitoringTransitionTwo::ReactionMonitoringTransitionTwo() :
    library_intensity_(-101),
    decoy_type_(UNKNOWN),
    precursor_mz_(0.0)
  {
    // Default is: true, false, true
    // NOTE: do not change that, the same default is implicitly assumed in TraMLHandler
    transition_flags_[DETECTING_TRANSITION_LOC] = true;
    transition_flags_[IDENTIFYING_TRANSITION_LOC] = false;
    transition_flags_[QUANTIFYING_TRANSITION_LOC] = true;
  }

  ReactionMonitoringTransitionTwo::ReactionMonitoringTransitionTwo(const ReactionMonitoringTransitionTwo & rhs) :
    name_(rhs.name_),
    peptide_ref_(rhs.peptide_ref_),
    compound_ref_(rhs.compound_ref_),
    library_intensity_(rhs.library_intensity_),
    decoy_type_(rhs.decoy_type_),
    precursor_mz_(rhs.precursor_mz_),
    product_(rhs.product_),
    transition_flags_(rhs.transition_flags_)
  {
  }

  ReactionMonitoringTransitionTwo::ReactionMonitoringTransitionTwo(ReactionMonitoringTransitionTwo && rhs) noexcept :
    name_(std::move(rhs.name_)),
    peptide_ref_(std::move(rhs.peptide_ref_)),
    compound_ref_(std::move(rhs.compound_ref_)),
    library_intensity_(std::move(rhs.library_intensity_)),
    decoy_type_(std::move(rhs.decoy_type_)),
    precursor_mz_(std::move(rhs.precursor_mz_)),
    product_(std::move(rhs.product_)),
    transition_flags_(std::move(rhs.transition_flags_))
  {
  }

  ReactionMonitoringTransitionTwo::~ReactionMonitoringTransitionTwo()
  {
  }

  ReactionMonitoringTransitionTwo & ReactionMonitoringTransitionTwo::operator=(const ReactionMonitoringTransitionTwo & rhs)
  {
    if (&rhs != this)
    {
      name_ = rhs.name_;
      peptide_ref_ = rhs.peptide_ref_;
      compound_ref_ = rhs.compound_ref_;
      precursor_mz_ = rhs.precursor_mz_;
      product_ = rhs.product_;
      library_intensity_ = rhs.library_intensity_;
      decoy_type_ = rhs.decoy_type_;
      transition_flags_ = rhs.transition_flags_;
    }
    return *this;
  }

  ReactionMonitoringTransitionTwo & ReactionMonitoringTransitionTwo::operator=(ReactionMonitoringTransitionTwo && rhs) noexcept
  {
    if (&rhs != this)
    {
      name_ = std::move(rhs.name_);
      peptide_ref_ = std::move(rhs.peptide_ref_);
      compound_ref_ = std::move(rhs.compound_ref_);
      precursor_mz_ = std::move(rhs.precursor_mz_);
      product_ = std::move(rhs.product_);
      library_intensity_ = std::move(rhs.library_intensity_);
      decoy_type_ = std::move(rhs.decoy_type_);
      transition_flags_ = std::move(rhs.transition_flags_);
    }
    return *this;
  }

  bool ReactionMonitoringTransitionTwo::operator==(const ReactionMonitoringTransitionTwo & rhs) const
  {
    return name_ == rhs.name_ &&
           peptide_ref_ == rhs.peptide_ref_ &&
           compound_ref_ == rhs.compound_ref_ &&
           precursor_mz_ == rhs.precursor_mz_ &&
           product_ == rhs.product_ &&
           library_intensity_ == rhs.library_intensity_ &&
           decoy_type_ == rhs.decoy_type_ &&
           transition_flags_ == rhs.transition_flags_;
  }

  bool ReactionMonitoringTransitionTwo::operator!=(const ReactionMonitoringTransitionTwo & rhs) const
  {
    return !(*this == rhs);
  }

  void ReactionMonitoringTransitionTwo::setName(const String & name)
  {
    name_ = name;
  }

  const String & ReactionMonitoringTransitionTwo::getName() const
  {
    return name_;
  }

  void ReactionMonitoringTransitionTwo::setNativeID(const String & name)
  {
    name_ = name;
  }

  const String & ReactionMonitoringTransitionTwo::getNativeID() const
  {
    return name_;
  }

  void ReactionMonitoringTransitionTwo::setPeptideRef(const String & peptide_ref)
  {
    peptide_ref_ = peptide_ref;
  }

  const String & ReactionMonitoringTransitionTwo::getPeptideRef() const
  {
    return peptide_ref_;
  }

  void ReactionMonitoringTransitionTwo::setCompoundRef(const String & compound_ref)
  {
    compound_ref_ = compound_ref;
  }

  const String & ReactionMonitoringTransitionTwo::getCompoundRef() const
  {
    return compound_ref_;
  }

  void ReactionMonitoringTransitionTwo::setPrecursorMZ(double mz)
  {
    precursor_mz_ = mz;
  }

  double ReactionMonitoringTransitionTwo::getPrecursorMZ() const
  {
    return precursor_mz_;
  }

  void ReactionMonitoringTransitionTwo::setProductMZ(double mz)
  {
    product_.setMZ(mz);
  }

  double ReactionMonitoringTransitionTwo::getProductMZ() const
  {
    return product_.getMZ();
  }

  int ReactionMonitoringTransitionTwo::getProductChargeState() const
  { 
    return product_.getChargeState();
  }

  bool ReactionMonitoringTransitionTwo::isProductChargeStateSet() const
  { 
    return product_.hasCharge();
  }

  void ReactionMonitoringTransitionTwo::setProduct(ReactionMonitoringTransitionTwo::Product product)
  {
    product_ = std::move(product);
  }

  const ReactionMonitoringTransitionTwo::Product & ReactionMonitoringTransitionTwo::getProduct() const
  {
    return product_;
  }

  void ReactionMonitoringTransitionTwo::updateMembers_()
  {
  }

  ReactionMonitoringTransitionTwo::DecoyTransitionType ReactionMonitoringTransitionTwo::getDecoyTransitionType() const
  {
    return decoy_type_;
  }

  void ReactionMonitoringTransitionTwo::setDecoyTransitionType(const DecoyTransitionType & d)
  {
    decoy_type_ = d;
  }

  double ReactionMonitoringTransitionTwo::getLibraryIntensity() const
  {
    return library_intensity_;
  }

  void ReactionMonitoringTransitionTwo::setLibraryIntensity(const double intensity)
  {
    library_intensity_ = intensity;
  }

  bool ReactionMonitoringTransitionTwo::isDetectingTransition() const
  {
    return transition_flags_[DETECTING_TRANSITION_LOC];
  }

  void ReactionMonitoringTransitionTwo::setDetectingTransition(bool val)
  {
    transition_flags_[DETECTING_TRANSITION_LOC] = val;
  }

  bool ReactionMonitoringTransitionTwo::isIdentifyingTransition() const
  {
    return transition_flags_[IDENTIFYING_TRANSITION_LOC];
  }

  void ReactionMonitoringTransitionTwo::setIdentifyingTransition(bool val)
  {
    transition_flags_[IDENTIFYING_TRANSITION_LOC] = val;
  }

  bool ReactionMonitoringTransitionTwo::isQuantifyingTransition() const
  {
    return transition_flags_[QUANTIFYING_TRANSITION_LOC];
  }

  void ReactionMonitoringTransitionTwo::setQuantifyingTransition(bool val)
  {
    transition_flags_[QUANTIFYING_TRANSITION_LOC] = val;
  }

} // namespace OpenMS
