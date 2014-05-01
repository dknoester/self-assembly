library(stringr)

find.files <- function(file_pattern=".*", dir_pattern=".*", path=".") {
	grep(dir_pattern, 
		list.files(recursive=TRUE, pattern=file_pattern, path=path, full.names=TRUE),
		value=TRUE)
}

# find.files <- function(file_pattern, path="./", dir_pattern=".*") {
	# unlist(lapply(dir(path=path, pattern=dir_pattern),
		# function(x) {
			# list.files(path=paste(path,x,sep="/"), pattern=file_pattern, full.names=TRUE)
		# }
	# ))
# }

filename.tags <- function(filename) {
	m = str_match(filename,"([-[:alnum:]]+)/([[:alnum:]]+)_([[:digit:]]+)/[[:alnum:]\\._]+$")
	return(list(filename=m[[1]],
		expr=as.factor(m[[2]]),
		treatment=as.factor(m[[3]]),
		trial = as.factor(m[[4]])))
}

load.files <- function(file_pattern, dir_pattern=".*", path="./", tag=TRUE, ...) {
	files = find.files(file_pattern, dir_pattern, path)
	data.frame(do.call(rbind,lapply(files, 
		function(x) {
			if(length(grep("\\.gz$",x))) {
				y=read.table(gzfile(x), header=TRUE, ...)
			} else {
				y=read.table(file(x), header=TRUE, ...)
			}
			if(tag) {
				t = filename.tags(x)
				y$filename = x
				y$expr = t$expr
				y$treatment = t$treatment
				y$trial = t$trial
			} else {
				return(y)
			}
			return(y)
		}
	)))
}
