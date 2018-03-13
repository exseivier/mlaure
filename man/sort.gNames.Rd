% File man/sort.gNames.RData
\name{sort.gNames}
\alias{sort.gNames}
\title{Sorts headers of sequences based on an expression contrast analysis}
\description{
	Sorts headers of sequences based on an expression contrast analysis using different statistics.
	"hvalue": is the honest vale, and was calculated as the quotient of log2FC / pvalue.
	"log.pvalue" & "stat": are statistics calculated as the log10(pvalue) and assigning the sign of
	the log2FC. "stat" correlates with "log.pvalue" but I do not know how this statistic was calculated :D.
	"log.fc": This is the log Fold Change obtained from the expression analysis.

	This function returns a SimpleList S4Vector R object with the datasets of the different contrast analysis
	performed in DESeq or DESeq2 R package, containing the sorted genes using the different statistics.
	Also, this function writes to files only the headers of the sorted genes. The output files have the *.hd extension.

	This function takes the path where the *.csv files are. These files are the results table of the contrast analysis
	performed with DESeq or DESeq2 R packages.
	}
\usage{
	data <- sort.gNames(path=Path/to/csv/file, sort.by="stat")
	}
\arguments{
	\item{path}{Character string. It is the path where the *.csv file are. These files are the result table fo the contrast analysis performed with
	DESeq or DESeq2 R packages}
	\item{sort.by}{Character string. it is used to select the statistics which will be used to sort the genes}
	}
\seealso{
	\code{\link{mlore}}
	}
\examples{
	data <- sort.gNames(path=Path/to/csv/file, sort.by="stat")
	}
\keyword{filtering}
