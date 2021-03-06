% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_import.R
\name{import_batch_tidy}
\alias{import_batch_tidy}
\title{Import single batch PLFA peak list from Excel}
\usage{
import_batch_tidy(file_path)
}
\arguments{
\item{file_path}{Full file path, including file extension for the Excel file
containing PLFA data. The input Excel file must contain Batch, DataFileName,
RetTimeSecs, MajorHeightnA, TotalPeakArea1, DisplayDelta1, and Name columns.}
}
\value{
Tibble suitable for downstream analyses in plfaR. Contains the
input columns as well as BatchDataFileName.
}
\description{
This is the first function in the plfaR workflow, loading named PLFA peak
list from a single batch for use in downstream analysis.
}
\examples{

}
