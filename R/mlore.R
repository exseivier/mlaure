#	METHODS
setGeneric("theBiggestTheBest", function(seqs, cutoff=50L) standardGeneric("theBiggestTheBest"))
setMethod("theBiggestTheBest", signature("DNAStringSet"),
	function(seqs, cutoff=50L) {
		gNames <- names(seqs)
		stopifnot(sum(grepl("[[:alpha:]]G[[:digit:]]", gNames)) == length(gNames))
			#	Selecting DNA string length higher than or equal to 50 pb
		seqs <- seqs[width(seqs) >= cutoff]
			#	Sorting DNA strings by string length
		seqs <- seqs[order(width(seqs), decreasing=T)]
			#	Picking-up no duplicated names
		seqs <- seqs[!duplicated(names(seqs))]
			#	Returning the modified DNAStringSet object
		seqs
})


setGeneric("sort.gNames", function(path, sort.by) standardGeneric("sort.gNames"))
setMethod("sort.gNames", signature("character", "character"),
	function(path, sort.by) {
		files <- list.files(path=path, pattern=".csv$")
		names <- gsub("(.*)\\..*", "\\1", files)
		data <- SimpleList()
		for(name in names){
			data[[name]] <- read.delim(paste(path, "/", name, ".csv", sep=""), header=T, sep=",", as.is=TRUE)
			if (sort.by == "hvalue") {
				data[[name]]$hvalue <- data[[name]]$log2FoldChange / data[[name]]$pvalue
				data[[name]] <- data[[name]][order(data[[name]]$hvalue),]
				data[[name]]$id <- as.character(data[[name]]$id)
				outExt <- ".hv"
			}
			if (sort.by == "log.pvalue") {
				data[[name]]$log.pvalue <- ifelse(data[[name]]$log2FoldChange < 0, log10(data[[name]]$pvalue) * -1, log10(data[[name]]$pvalue))
				data[[name]] <- data[[name]][order(data[[name]]$log.pvalue),]
				data[[name]]$id <- as.character(data[[name]]$id)
				outExt <- ".lpv"
			}
			if (sort.by == "log.fc") {
				data[[name]] <- data[[name]][order(data[[name]]$log2FoldChange),]
				outExt <- ".lfc"
			}
			if (sort.by == "stat") {
				data[[name]] <- data[[name]][order(data[[name]]$stat),]
				outExt <- ".st"
			}
			print(paste(name, ".hd", sep=""))
			write(data[[name]]$id, file=paste(name, outExt, ".hd", sep=""))
		}
		cat("Sorted headers are stored in *.hd files\n\n")
		data
})

setGeneric("plot.sylamer", function(path, mirnas, control=FALSE) standardGeneric("plot.sylamer"))
setMethod("plot.sylamer", signature("character", "character"),
	function(path, mirnas, control=FALSE) {
		syl.files <- list.files(path=path, pattern=".syl$")
		mirnas <- read.delim(mirnas, header=F, sep="\t", as.is=TRUE)
		colours <- c("blue", "red", "green", "orange", "black")
		outSeeds <- c()
		for (file in syl.files) {
			
			tmp.syl.data <- read.delim(paste(path, "/", file, sep=""), header=T, sep="\t", row.names=1)
			tmp.syl.data <- cbind("X0"=0, tmp.syl.data)
			x <- colnames(tmp.syl.data)
			x <- as.numeric(gsub("^X", "", x))
			seeds <- mirnas[,1]
			mirnas[,2] <- unlist(lapply(strsplit(mirnas[,2], split=","), function(x) x[1]))
			mirnas$V3 <- gsub("mmu-...-|\\..+", "", mirnas[,2])
			mirnas$V4 <- paste0(mirnas$V1, " (", mirnas$V3, ")")

			if (control) {
				tmp.syl.data <- tmp.syl.data[!(rownames(tmp.syl.data) %in% seeds),]
			}
			else {
				tmp.syl.data <- tmp.syl.data[rownames(tmp.syl.data) %in% seeds,]
			}

			tmp.syl.data.up <- tmp.syl.data[order(apply(tmp.syl.data[,1:(length(tmp.syl.data[1,])/2)], 1, max), decreasing=T),]
			plot(x, tmp.syl.data.up[1,], col="grey", type="l", lty=1, ylab="-log10(p-value)", xlab="Mus musculus 3'UTRs", main=file, ylim=c(-20, 20))
			
			for (i in 2:length(tmp.syl.data.up[,1])) {
				lines(x, tmp.syl.data.up[i,], col="grey", lty=1)
			}

			tmp.syl.data.up <- head(tmp.syl.data.up, n=5)
			tmp.syl.data.do <- tmp.syl.data[order(apply(tmp.syl.data[,(length(tmp.syl.data[1,])/2):length(tmp.syl.data[1,])], 1, min)),]
			tmp.syl.data.do <- head(tmp.syl.data.do, n=5)

			for(i in 1:length(tmp.syl.data.up[,1])) {
				lines(x, tmp.syl.data.up[i,], col=colours[i], lty=1)
			}
			for (i in 1:length(tmp.syl.data.do[,1])) {
				lines(x, tmp.syl.data.do[i,], col=colours[i], lty=1)
			}
			
			legend("topright", legend=mirnas[mirnas$V1 %in% rownames(tmp.syl.data.up), "V4"], col = colours, lty=1)
			legend("bottomright", legend=mirnas[mirnas$V1 %in% rownames(tmp.syl.data.do), "V4"], col = colours, lty=1)
		}
})
