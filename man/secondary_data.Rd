% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis_2.R
\name{secondary_data}
\alias{secondary_data}
\title{Load secondary data}
\usage{
secondary_data(m, con)
}
\arguments{
\item{m}{Rmonad pipeline from stage 2}
}
\description{
Loop over all target species, produce the secondary data needed to classify
orphans. This step will take a long time (many CPU hours).
}
