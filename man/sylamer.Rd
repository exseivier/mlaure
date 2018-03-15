% File man/sylamer.Rd

\name{sylamer}
\alias{sylamer}
\title{Sneds sylamer task to a CPU process}
\description{
	This function take all arguments needed to run Sylamer inside R with a system call. This function performs the task in local machine and in a node in Mazorka cluster.
}

\usage{
	sylamer(headers=headers, ext=NULL, utrs=NULL, stage="local", wd=NULL)
}

\arguments{
	\item{headers}{Ths argument should be a character vector or a SimpleList S4 Object. When it is a character vector it should specify the path where the headers files are. But if it is a SimpleList, it has to contain the data imported with sort.gNames from the *.csv files containing the result table of expression contrast analysis performed with DESeq or DESeq2 R package.}
	\item{ext}{This is a character vector and specifies the extension of the headers files you want to use in sylamer analysis. (NULL as default).}
	\item{genome}{Character vector. Specifies the genome file path and name. (NULL as default).}
	\item{stage}{Character vector which specifies the machine stage: it could be "local" (as default) if you want to run sylamer locally, or "cluster" if you want to run sylamer in a cluster (with "qsub" manager)}
	\item{wd}{Character vector. Specifies the working directory where sylamer is going to run. (NULL as default).}
}

\examples{
	sylamer(headers="/Path/to/headers", ext=NULL, utrs="/Path/to/utrs.fasta", stage="local", wd="Wirking/Dir")
}

\seealso{
	\code{\link{plot.sylamer}}
}

\keyword{sylamer}
