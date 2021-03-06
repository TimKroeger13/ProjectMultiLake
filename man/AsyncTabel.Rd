% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AsyncTabel.R
\name{AsyncTabel}
\alias{AsyncTabel}
\title{AsyncTabel}
\usage{
AsyncTabel(
  data,
  intervallBy = 100,
  Importname = "",
  Exportname = "",
  CreateVector = FALSE,
  TransformAllData = FALSE,
  vectorName,
  ValueCantBeSamlerThanZero = FALSE
)
}
\arguments{
\item{data}{List of data generates by the Multivar function.}

\item{intervallBy}{Intervalls by to interpolate to.}

\item{Importname}{importname after data$Diatom$DiatomNames$}

\item{Exportname}{data$Diatom$DiatomNames$}

\item{CreateVector}{Shuld the vector plot build on this dataset.}

\item{TransformAllData}{Z Transformas all Data.}

\item{vectorName}{Outpurname for the vector data.}

\item{ValueCantBeSamlerThanZero}{If true, values smaller than 0 become 0.}
}
\value{
Returns the same List but with new added parameters.
}
\description{
Calculates the an AsyncTable given by prior function.
}
\note{
This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
\cr Comma numbers are rounded up.
}
\author{
Tim Kroeger
}
