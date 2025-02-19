% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easyGO.R
\name{search_GO_by_term}
\alias{search_GO_by_term}
\title{search_GO_by_term}
\usage{
search_GO_by_term(search_term, search_field = c("TERM", "DEFINITION")[1])
}
\arguments{
\item{search_term}{character term to search for in TERM (or DEFINITION) of GO
database.}

\item{search_field}{character field to search in GO database.}
}
\value{
character vector of GO ids containing search_term.
}
\description{
search_GO_by_term
}
\examples{
search_GO_by_term("MAPK")
search_GO_by_term("MAPK", search_field = "DEFINITION")
}
