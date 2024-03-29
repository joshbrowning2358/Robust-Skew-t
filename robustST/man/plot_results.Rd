% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/plot_results.R
\name{plot_results}
\alias{plot_results}
\title{Plot results}
\usage{
plot_results(results, prefix, xi = 0, omega = 1, alpha = 0, nu = 10000)
}
\arguments{
\item{results}{A results object, as returned by runSim.}

\item{prefix}{A prefix to be added to the saved plots.}

\item{xi}{A parameter of the skew-t distribution.  This should be the
parameter used in the simulation.}

\item{omega}{A parameter of the skew-t distribution.  This should be the
parameter used in the simulation.}

\item{alpha}{A parameter of the skew-t distribution.  This should be the
parameter used in the simulation.}

\item{nu}{A parameter of the skew-t distribution.  This should be the
parameter used in the simulation.}
}
\value{
No results are returned, but plots are saved in the user's current
directory.
}
\description{
This function takes the results object as created by runSim and generates
plots summarizing the simulation.  These plots are saved in the users
current working directory.
}

