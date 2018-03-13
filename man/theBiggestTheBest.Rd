% File man/theBiggestTheBest.Rd
\name{theBiggestTheBest}
\alias{theBiggestTheBest}
\title{Chooses the biggest sequence if it is duplicated in fasta file}
\description{
	Chooses the biggest sequence if there is a duplicated header in the fasta file.
}
\usage{seqs <- theBiggestTheBest(seqs=seqs, cutoff=100L)}
\arguments{
	\item{seqs}{DNAStringSet object from "Biostrings" containing the sequences}
	\item{cutoff}{An integer cutoff value to filter sequences size (pb). 50L is a default value}
	}
\seealso{
	\code{\link{Biostrings}}
	}
\examples{
	## Load Mus musculus sequences in a DNAStringSet object
	mmus <- data(mmus)
	## Using 50L as cutoff default value
	mmus <- theBiggestTheBest(seqs=mmus)
	## Using 200L as cutoff value
	mmus <- theBiggestTheBest(seqs=mmus, cutoff=200L)
	}
\keyword{filtering}
