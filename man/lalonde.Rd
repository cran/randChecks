\name{lalonde}
\alias{lalonde}
\docType{data}
\title{Lalonde (1986) Data}
\description{
Data from Lalonde (1986).}
\usage{data(lalondeMatches)}
\format{
  The full Lalonde (1986) dataset, containing 614 units (rows) and 9 variables (columns). The columns are:
  \itemize{
    \item \code{treat}: A binary treatment variable. Equal to 1 if treated in the National Supported Work Demonstration; equal to 0 otherwise.
    \item \code{age}: age in years.
    \item \code{educ}: years of education.
    \item \code{black}: an indicator variable, equal to 1 only if the subject is black.
    \item \code{hispan}: an indicator variable, equal to 1 only if the subject is hispanic.
    \item \code{married}: an indicator variable, equal to 1 only if the subject is married.
    \item \code{nodegree}: an indicator variable, equal to 1 only if the subject does not have a degree.
    \item \code{re74}: earnings in 1974.
    \item \code{re75}: earnings in 1975.
  }
  All of the columns except \code{treat} are covariates; in these datasets, the outcome variable is not provided.
}

\examples{
data(lalonde)
lalonde
}

\references{	
	LaLonde, R. J. (1986). Evaluating the econometric evaluations of training programs with experimental data. The American Economic Review, 604-620.
}

\keyword{datasets}
\keyword{lalonde}
