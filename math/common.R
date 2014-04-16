find.files <- function(file_pattern, dir_pattern=".*", path="./") {
	unlist(lapply(dir(path=path, pattern=dir_pattern),
		function(x) {
			list.files(path=paste(path,x,sep="/"), pattern=file_pattern, full.names=TRUE)
		}
	))
}

load.files <- function(files, tag=TRUE, header=TRUE,...) {
	data.frame(do.call(rbind,lapply(files, 
		function(x) {
			if(length(grep("\\.gz$",x))) {
				y=read.table(gzfile(x), header=header, ...)
			} else {
				y=read.table(file(x), header=header, ...)
			}
			if(tag) {
				y$treatment = as.factor(sub("_","",regmatches(x,regexpr("[\\w\\d-]+_",x,perl=TRUE))))
				y$trial = as.factor(gsub("[_/]","",regmatches(x,regexpr("_\\d+/",x,perl=TRUE))))
				y$filename = x
			}
			return(y)
		}
	)))
}
