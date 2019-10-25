\name{bppInputs}
\alias{bppInputs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bppInputs
}
\description{
Generates input directories and CTL files for BPP
Test
}
\usage{
bppInputs(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wd}{
  path for working directory}
  \item{tree}{
  maximally split tree in newick format
  }
  \item{mapping}{
  txt file, one column with tip names, one column with corresponding species names
  }
  \item{priors}{
  tab delimited file, one line per combination of priors with priors in the following order: tau_alpha tau_beta theta_alpha theta_beta
  }
  \item{heredity}{
  heredity file for loci, see BPP documentation
  }
  \item{loci}{
  loci file, see BPP documentation
  }
  \item{ctl}{
  BPP template ctl file created by BPPCtlTemplate function
  }
  \item{nloci}{
  number of loci
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
