from Types cimport *

from RNPxlMarkerIonExtractor cimport *
# from ListUtilsIO cimport *
from PeptideIdentification cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlReport "OpenMS::RNPxlReport":
        RNPxlReport() nogil except + # compiler
        RNPxlReport(RNPxlReport &) nogil except + # compiler
        libcpp_vector[ RNPxlReportRow ] annotate(MSExperiment & spectra,
                                                 libcpp_vector[ PeptideIdentification ] & peptide_ids,
                                                 double marker_ions_tolerance) nogil except +


cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlReportRowHeader "OpenMS::RNPxlReportRowHeader":
        RNPxlReportRowHeader() nogil except + # compiler
        RNPxlReportRowHeader(RNPxlReportRowHeader &) nogil except + # compiler
        String getString(const String & separator) nogil except +

cdef extern from "<OpenMS/ANALYSIS/RNPXL/RNPxlReport.h>" namespace "OpenMS":
    
    cdef cppclass RNPxlReportRow "OpenMS::RNPxlReportRow":
        RNPxlReportRow() nogil except + # compiler
        RNPxlReportRow(RNPxlReportRow &) nogil except + # compiler

        bool no_id
        double rt
        double original_mz
        String accessions
        String RNA
        String peptide
        double best_localization_score
        String localization_scores
        String best_localization
        Int charge
        double score
        double peptide_weight
        double RNA_weight
        double xl_weight
        double abs_prec_error
        double rel_prec_error
        libcpp_map[String, libcpp_vector[ libcpp_pair[double, double] ] ] marker_ions # RNPxlMarkerIonExtractor::MarkerIonsType marker_ions
        double m_H
        double m_2H
        double m_3H
        double m_4H
        int rank

        String getString(const String & separator) nogil except +
