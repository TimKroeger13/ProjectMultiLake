% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpeciesRichness.R
\name{speciesRichness}
\alias{speciesRichness}
\title{speciesRichness}
\usage{
speciesRichness(
  data,
  intervallBy = 100,
  NonNegative = TRUE,
  Importname1 = "",
  Importname2 = "",
  Exportname = ""
)
}
\arguments{
\item{data}{List of data generates by the Multivar function.}

\item{intervallBy}{Intervalls by to interpolate to.}

\item{NonNegative}{Creates Positive Values after Loess calculation.}

\item{Importname1}{importname 1.}

\item{Importname2}{importname 2.}

\item{Exportname}{data$Diatom$DiatomNames$}
}
\value{
Returns the same List but with new added parameters.
}
\description{
Calculates the inverse species richness for the multivar function.
}
\note{
This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
\cr Comma numbers are rounded up.
}
\author{
Tim Kroeger
}
