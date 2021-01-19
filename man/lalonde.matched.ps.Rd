\name{lalonde.matched.ps}
\alias{lalonde.matched.ps}
\docType{data}
\title{A 1:1 Propensity Score Matched Dataset for the Lalonde (1986) Data}
\description{
A optimal 1:1 propensity score matched dataset for the Lalonde (1986) dataset.}
\usage{data(lalondeMatches)}
\format{
  370 units (rows) and 10 variables (columns). The columns are:
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
    \item \code{subclass}: The subclass denoting the pairs within the matched dataset.
  }
}

\examples{
	data(lalondeMatches)
	lalonde.matched.ps
}

\details{

The optimal 1:1 propensity score matched dataset was produced using the \code{MatchIt} \code{R} package. The propensity scores were estimated using logistic regression, where \code{treat} was the outcome and the other variables were the covariates (within no interactions included).
}

\references{	
	LaLonde, R. J. (1986). Evaluating the econometric evaluations of training programs with experimental data. The American Economic Review, 604-620.
}

\keyword{datasets}
\keyword{lalonde}
