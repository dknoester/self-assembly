library(ggplot2)
source("~/research/src/self-assembly/math/common.R")
# source("~/research/src/self-assembly/math/fixup_ggplot.R")
setwd("/Users/dk/research/src/self-assembly/var")

STYLE="quick" # or "all"

quick_theme <- theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(colour = "black",size=0.75), legend.title=element_blank()) 

fitness_plot <- function(path) {
	D = rbind(load.files(find.files("fitness.dat.gz",path=path)))
	
	if(STYLE=="quick") {
		D = subset(D,update%%100==0)		
	}
	
	quartz(width=6,height=3.75)
	g = ggplot(data=D, aes(x=update, y=max_fitness)) + stat_summary(aes(color=treatment,fill=treatment),fun.data="mean_cl_boot", geom="smooth") + labs(x="Update", y="Fitness") + quick_theme + ylim(0,1) + xlim(0,10000)
	print(g)
	return(list(D,g))
}

## 010-3d-fsm-rand
#
x010 = fitness_plot("010-3d-fsm-rand")

## 009-2d-fsm-rand
#
x009 = fitness_plot("009-2d-fsm-rand")

## 008-1d-fsm-rand
#
x008 = fitness_plot("008-1d-fsm-rand")

## 007-3d-fsm-rl
#
x007 = fitness_plot("007-3d-fsm-rl")

## 006-2d-fsm-rl
#
x006 = fitness_plot("006-2d-fsm-rl")

## 005-1d-fsm-rl
#
x005 = fitness_plot("005-1d-fsm-rl")

## 004-3dd-fsm
#
x004 = fitness_plot("004-3d-fsm")

## 003-2d-fsm
#
x003 = fitness_plot("003-2d-fsm")

## 002-1d-fsm
#
x002 = fitness_plot("002-1d-fsm")

## 001-ga
#
x001 = fitness_plot("001-ga")


