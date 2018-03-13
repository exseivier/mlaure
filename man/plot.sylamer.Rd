% File man/plot.sylamer.Rd

\name{plot.sylamer}
\alias{plot.sylamer}
\title{Plots the Sylamer results}
\description{
	Takes the Sylamer results and microRNA seeds file and plots the log10(pvalue) of seeds enrichment based
	on the hypergeometric test for every miRNA seed, and it highlights the top 5 seeds for enrichment and depletion.
	Also, this function performs a control test. This test comprises in ploting every k-mer that does not belong to the
	microRNA seeds list. This control test is intended to check the noise of the data.
}
\usage{
	plot.sylamer(path="Path/to/Sylamer/resul/files", mirnas="miRNAsSeedFile.txt", control=FALSE)
}
\arguments{
	\item{path}{Character string. Path to sylamer output files. These files should have a *.syl extension because the function
	searches for this extension.}
	\item{mirnas}{Character string. It is the name of the microRNAs seeds.}
	\item{control}{Boolean. If it is set TRUE; plot.sylamer performs the control test.}
}
\seealso{
	\code{\link{mlore}}
}
\examples{
	plot.sylamer(path="Path/to/Sylamer/resul/files", mirnas="miRNAsSeedFile.txt", control=FALSE)
}
\keyword{plotting}
