\name{lalondeMatches}
\alias{lalondeMatches}
\docType{data}
\title{Lalonde (1986) Data and Two Matched Datasets}
\description{
Data from Lalonde (1986) and two matched datasets: One where optimal 1:1 propensity score matching was used, and one where cardinality matching was used, with the balance constraint that all standardized covariate mean differences be below 0.1.}
\usage{data(lalondeMatches)}
\format{
  Three data frames:
  \itemize{
  	\item \code{lalonde}: 614 units (rows) and 9 variables (columns). This is the full Lalonde (1986) dataset.
  	\item \code{lalonde.matched.ps}: 370 units (rows) and 10 variables (columns). This is the 1:1 propensity score matched dataset.
  	\item \code{lalonde.matched.card}: 240 units (rows) and 10 variables (columns). This is the cardinality matched dataset.
  }
  All three data frames have these 9 columns:
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
  
  Meanwhile, \code{lalonde.matched.ps} and \code{lalonde.matched.card} have one additional column, \code{subclass}, denoting the pairs for those matched datasets.
}

\examples{
data(lalondeMatches)
}

\details{

The optimal 1:1 propensity score matched dataset was produced using the \code{MatchIt} \code{R} package. The propensity scores were estimated using logistic regression, where \code{treat} was the outcome and the other variables were the covariates (within no interactions included).

The cardinality matched datset was produced using the \code{designmatch} \code{R} package.
}

\references{	
	LaLonde, R. J. (1986). Evaluating the econometric evaluations of training programs with experimental data. The American Economic Review, 604-620.
}

\keyword{datasets}
\keyword{lalonde}
