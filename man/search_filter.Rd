% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_filter.R
\name{search_filter}
\alias{search_filter}
\title{search_filter}
\usage{
search_filter(filter, type = "SIM")
}
\arguments{
\item{filter}{character vector of the Orbitrap filter string}

\item{type}{'SIM' for true wide SIM DIA data or any other value for pseudo MS1 0 volt data \(true MS2\)}
}
\value{
data.frame contraining segment search results
}
\description{
search_filter
}
