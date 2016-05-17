% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_mend_err.R
\name{parse_mend_err}
\alias{parse_mend_err}
\title{parse_mend_err: Parse output from CriMAP .loc files}
\usage{
parse_mend_err(prefile, genfile, save.mendfile = TRUE, familyPedigree)
}
\arguments{
\item{prefile}{.pre file (output from prepare function)}

\item{genfile}{.gen file}

\item{save.mendfile}{logical. Default = TRUE, saves .mnd file}

\item{familyPedigree}{data.frame containing columns ANIMAL, FATHER, MOTHER
and FAMILY. FAMILY defines family groups for crimap. IDs in ANIMAL can be
repeated within and between families if necessary.}
}
\description{
parse_mend_err: Parse output from CriMAP .loc files
}
