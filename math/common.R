library(stringr)

# find.files <- function(file_pattern, dir_pattern=".*", path="./") {
	# unlist(lapply(dir(path=path, pattern=dir_pattern),
		# function(x) {
			# list.files(path=paste(path,x,sep="/"), pattern=file_pattern, full.names=TRUE)
		# }
	# ))
# }

find.files <- function(file_pattern=".*", path=".") {
	list.files(recursive=TRUE, pattern=file_pattern, path=path, full.names=TRUE)
}

# list.files(recursive=TRUE,pattern="fitness.dat.gz")

load.files <- function(files, tag=TRUE, header=TRUE,...) {
	data.frame(do.call(rbind,lapply(files, 
		function(x) {
			if(length(grep("\\.gz$",x))) {
				y=read.table(gzfile(x), header=header, ...)
			} else {
				y=read.table(file(x), header=header, ...)
			}
			if(tag) {
				m = str_match(x,"([-[:alnum:]]+)/([[:alnum:]]+)_([[:digit:]]+)/[[:alnum:]\\.]+$")
				y$filename = m[[1]]
				y$expr = as.factor(m[[2]])
				y$treatment = as.factor(m[[3]])
				y$trial = as.factor(m[[4]])
				
# # 				y$treatment = as.factor(gsub("[_/]","",regmatches(x,regexpr("/[\\w\\d-]+_",x,perl=TRUE))))
				# y$trial = as.factor(gsub("[_/]","",regmatches(x,regexpr("_\\d+/",x,perl=TRUE))))
				# y$expr = as.factor(gsub("[_/]","",regmatches(x,regexpr("/\\d+[\\w\\d-]+/",x,perl=TRUE))))
				# y$filename = x
			}
			return(y)
		}
	)))
}
