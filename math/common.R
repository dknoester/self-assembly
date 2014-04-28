library(stringr)

# find.files <- function(file_pattern=".*", path=".") {
	# list.files(recursive=TRUE, pattern=file_pattern, path=path, full.names=TRUE)
# }

find.files <- function(file_pattern, path="./", dir_pattern=".*") {
	unlist(lapply(dir(path=path, pattern=dir_pattern),
		function(x) {
			list.files(path=paste(path,x,sep="/"), pattern=file_pattern, full.names=TRUE)
		}
	))
}

load.files <- function(file_pattern, path="./", dir_pattern=".*", tag=TRUE, ...) {
	files = find.files(file_pattern, path, dir_pattern)
	data.frame(do.call(rbind,lapply(files, 
		function(x) {
			if(length(grep("\\.gz$",x))) {
				y=read.table(gzfile(x), header=TRUE, ...)
			} else {
				y=read.table(file(x), header=TRUE, ...)
			}
			if(tag) {
				m = str_match(x,"([-[:alnum:]]+)/([[:alnum:]]+)_([[:digit:]]+)/[[:alnum:]\\.]+$")
				y$filename = m[[1]]
				y$expr = as.factor(m[[2]])
				y$treatment = as.factor(m[[3]])
				y$trial = as.factor(m[[4]])
			}
			return(y)
		}
	)))
}
