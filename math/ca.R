library(ggplot2)
library(plyr)
library(stringr)
source("~/research/src/self-assembly/math/common.R")
setwd("/Users/dk/research/src/self-assembly/var")
figpath = "/Users/dk/research/src/self-assembly/doc/saso2014/figures/"

STYLE="draft" # or "final"
WIDTH=6
HEIGHT=3.75

quick_theme <- theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(colour = "black",size=0.75), legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(1,0))

# Plot the fitness of D by experiment
#
expr_fitness <- function(D) {
	if(STYLE=="draft") {
		D = subset(D,update%%100==0)		
	}
	g = ggplot(data=D, aes(x=update, y=max_fitness)) + stat_summary(aes(color=expr,fill=expr),fun.data="mean_cl_boot", geom="smooth") + labs(x="Update", y="Fitness") + quick_theme + ylim(0.5,1) + xlim(0,10000)
	return(g)
}

# Show figure x.
#
showfig <- function(x) {
	quartz(width=WIDTH,height=HEIGHT)
	print(x)
}

# Save figure x in in PDF format.
#
savefig <- function(x, name) {
	f = paste(figpath,name,".pdf",sep="")
	pdf(file=f, width=WIDTH, height=HEIGHT)
	print(x)
	dev.off()
}


# Load data for the different experiments.
x002 = load.files("fitness.dat.gz", "002-1d-fsm", "ta0")
x003 = load.files("fitness.dat.gz", "003-2d-fsm", "ta0")
x004 = load.files("fitness.dat.gz", "004-3d-fsm", "ta0")
x008 = load.files("fitness.dat.gz", "008-1d-fsm-rand", "ta0")
x009 = load.files("fitness.dat.gz", "009-2d-fsm-rand", "ta0")
x010 = load.files("fitness.dat.gz", "010-3d-fsm-rand", "ta0")

# Fitness of base experiments.
F_base_fitness = expr_fitness(rbind(x002,x003,x004))
F_base_fitness = F_base_fitness + scale_fill_discrete(breaks=c("002-1d-fsm","003-2d-fsm","004-3d-fsm"),labels=c("1D","2D","3D")) + scale_color_discrete(breaks=c("002-1d-fsm","003-2d-fsm","004-3d-fsm"),labels=c("1D","2D","3D"))
showfig(F_base_fitness)
savefig(F_base_fitness, "f001-base-fitness")

# Fitness of reinforcement learning experiments.
F_rl_fitness = expr_fitness(rbind(x008, x009, x010))
F_rl_fitness = F_rl_fitness + scale_fill_discrete(breaks=c("008-1d-fsm-rand","009-2d-fsm-rand","010-3d-fsm-rand"),labels=c("1D-rl","2D-rl","3D-rl")) + scale_color_discrete(breaks=c("008-1d-fsm-rand","009-2d-fsm-rand","010-3d-fsm-rand"),labels=c("1D-rl","2D-rl","3D-rl"))
showfig(F_rl_fitness)
savefig(F_rl_fitness, "f002-rl_fitness")

