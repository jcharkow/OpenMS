// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>

class QWidget;

namespace OpenMS
{
  class Annotation1DItem;

  /**
   * \brief Base class for all 1D layers, a special case of LayerData.
   *
   * 1D is a bit special because we need to remember which spectrum/chrom/IM is currently shown (there are usually many of them to choose from).
   *
   *
   */
  class OPENMS_GUI_DLLAPI LayerData1DBase : public virtual LayerDataBase
  {
  public:
    // rule of 0

    /**
     * \brief Obtain a painter which can draw the layer on a canvas
     * \return A painter
     */
    virtual std::unique_ptr<Painter1DBase> getPainter1D() const = 0;

    /**
     * \brief Given a @p partial_range for the current 1D layer (e.g. an m/z range), fill in the other
     *        dimensions (usually intensity) from all data points which are within the input range.
     * \param partial_range Range with at least one dimension populated (which is used to filter the current spectrum/chrom/...) 
     * \return Range of the data points within the input range (e.g. for spectra: m/z and intensity; or chroms: RT and intensity, etc)
     */
    virtual RangeAllType getRangeForArea(const RangeAllType partial_range) const = 0;

    /**
     * \brief Get a context menu (with lambda actions included) for this 1D layer, when a Annotation1DItem was right-clicked
     * \param annot_item The annotation item clicked on
     * \param need_repaint Reference of bool in calling function, which must know if the action requires repainting the canvas
     * \return A context menu to embed into the generic menu
     */
    virtual QMenu* getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint) = 0;


    /**
     * \brief Add a Annotation1DPeakItem to getCurrentAnnotations(). The specific type is determined by the derived class (e.g. Peak1D, ChromatogramPeak1D, etc)
     * \param peak_index Which peak should be annotated?
     * \param text Text to annotate with
     * \param color Color of the text
     * \return The item that was created
     */
    virtual Annotation1DItem* addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color) = 0;


    /// get name augmented with attributes, e.g. '*' if modified
    String getDecoratedName() const override;
        
    /// Returns a const reference to the annotations of the current spectrum (1D view)
    const Annotations1DContainer& getCurrentAnnotations() const
    {
      return annotations_1d_[current_idx_];
    }

    /// Returns a mutable reference to the annotations of the current spectrum (1D view)
    Annotations1DContainer& getCurrentAnnotations()
    {
      return annotations_1d_[current_idx_];
    }

    /// Returns a const reference to the annotations of the @p spectrum_index's spectrum (1D view)
    const Annotations1DContainer& getAnnotations(Size spectrum_index) const
    {
      return annotations_1d_[spectrum_index];
    }

    /// Returns a mutable reference to the annotations of the @p spectrum_index's spectrum (1D view)
    Annotations1DContainer& getAnnotations(Size spectrum_index)
    {
      return annotations_1d_[spectrum_index];
    }
    /// Get the index of the current spectrum (1D view)
    Size getCurrentIndex() const
    {
      return current_idx_;
    }

    /// Set the index of the current spectrum (1D view) -- and prepares annotations
    void setCurrentIndex(Size index);

    /// Does the layer have at least @p index items (e.g. spectra, chroms, etc), so a call to setCurrentIndex() is valid?
    virtual bool hasIndex(Size index) const = 0;

    /// if this layer is flipped (1d mirror view)
    bool flipped = false;

    /// Peak colors of the currently shown spectrum
    std::vector<QColor> peak_colors_1d;

  protected:
    /// Index of the current spectrum/chromatogram etc (by default, show the first one)
    Size current_idx_ = 0;

    /// Annotations of all spectra of the experiment (1D view)
    std::vector<Annotations1DContainer> annotations_1d_ = std::vector<Annotations1DContainer>(1);
  };

}// namespace OpenMS

