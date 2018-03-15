% File man/mlaure.Rd

\name{mlaure}
\alias{mlaure}
\title{R package to filter sequences, to sort headers, to perform Sylamer analysis & to plot results}
\description{
	This R package is intended to filter sequences, sort headers in order to perform seeds enrichment analysis using Sylamer and to plot results.
}
\arguments{
	This R package contains the functions needed to perform Sylamer analysis. Here I describe them:
	\item{theBiggestTheBest}{This function filters the sequences taking into account if the header is repeated into fasta file. If this is the case it will take the biggest sequence. To obtain more help visit ?theBiggestTheBest.}
	\item{sort.gNames}{It sorts the genes based on a expression contrast analysis, from the most repressed gene to the most induced gene. for more information, visit ?sort.gNames.}
	\item{plot.sylamer}{This function plots Sylamer results. For more information, visit ?plot.sylamer.}
	\item{mmus}{Mus musculus sequences stored in a DNAStringSet object. To extract the dataset execute: data(mmus)}
}
\keyword{R-package}

