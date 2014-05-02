library(ggplot2)
library(plyr)
library(stringr)
source("~/research/src/self-assembly/math/common.R")
setwd("/Users/dk/research/src/self-assembly/var")
figpath = "/Users/dk/research/doc/self-assembly/saso2014/non-embedded/"

STYLE="draft" # or "final"
WIDTH=6
HEIGHT=3.75

quick_theme <- theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(colour = "black",size=0.75), legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(1,0), text=element_text(family="Helvetica",size=14))


# Plot the fitness of D by experiment
#
expr_fitness <- function(D) {
	if(STYLE=="draft") {
		D = subset(D,update%%100==0)		
	}
	g = ggplot(data=D, aes(x=update, y=max_fitness)) + stat_summary(aes(color=expr,fill=expr),fun.data="mean_cl_boot", geom="smooth") + labs(x="Update", y="Fitness") + quick_theme + ylim(0.5,1)
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
savefig <- function(x, name, width=WIDTH, height=HEIGHT) {
	f = paste(figpath,name,".pdf",sep="")
	pdf(file=f, width=width, height=height)
	print(x)
	dev.off()
}

# Fitness data for the different experiments.
#
#
x002 = load.files("fitness.dat.gz", "ta0", "002-1d-fsm")
x003 = load.files("fitness.dat.gz", "ta0", "003-2d-fsm")
x004 = load.files("fitness.dat.gz", "ta0", "004-3d-fsm")
x008 = load.files("fitness.dat.gz", "ta0", "008-1d-fsm-rand")
x009 = load.files("fitness.dat.gz", "ta0", "009-2d-fsm-rand")
x010 = load.files("fitness.dat.gz", "ta0", "010-3d-fsm-rand")
x011 = load.files("fitness.dat.gz", "ta0", "011-1d-fsm-switch")
x012 = load.files("fitness.dat.gz", "ta0", "012-2d-fsm-switch")
x013 = load.files("fitness.dat.gz", "ta0", "013-3d-fsm-switch")

# Fitness of base experiments.
F_base_fitness = expr_fitness(rbind(x002,x003,x004))
F_base_fitness = F_base_fitness + scale_fill_discrete(breaks=c("002-1d-fsm","003-2d-fsm","004-3d-fsm"),labels=c("1D","2D","3D")) + scale_color_discrete(breaks=c("002-1d-fsm","003-2d-fsm","004-3d-fsm"),labels=c("1D","2D","3D"))
showfig(F_base_fitness)
savefig(F_base_fitness, "p-fitness")


# Fitness of reinforcement learning experiments.
F_rl_fitness = expr_fitness(rbind(x008, x009, x010))
F_rl_fitness = F_rl_fitness + scale_fill_discrete(breaks=c("008-1d-fsm-rand","009-2d-fsm-rand","010-3d-fsm-rand"),labels=c("1D-rl","2D-rl","3D-rl")) + scale_color_discrete(breaks=c("008-1d-fsm-rand","009-2d-fsm-rand","010-3d-fsm-rand"),labels=c("1D-rl","2D-rl","3D-rl"))
showfig(F_rl_fitness)
savefig(F_rl_fitness, "p-fitness-rl")

F_switch_fitness = expr_fitness(rbind(x011, x012, x013))
F_switch_fitness = F_switch_fitness + scale_fill_discrete(breaks=c("011-1d-fsm-switch","012-2d-fsm-switch","013-3d-fsm-switch"),labels=c("1D","2D","3D")) + scale_color_discrete(breaks=c("011-1d-fsm-switch","012-2d-fsm-switch","013-3d-fsm-switch"),labels=c("1D","2D","3D"))
showfig(F_switch_fitness)
savefig(F_switch_fitness, "p-fitness-reevolved")


# Adaptive data
#
#
A = load.files("ca_dom_adapt.dat")
A$delta_w = A$w0 - A$w1
levels(A$expr) = c("1D-rl", "2D-rl", "3D-rl")
a = ggplot(A, aes(expr,delta_w,fill=expr)) + theme(legend.position="none") + labs(x="", y="Fitness delta") + geom_boxplot() + geom_jitter() + ylim(0,1) + quick_theme
showfig(a)
savefig(a, "p-adaptive-fitness")


# Scaling data
#
#
S = load.files("ca_dom_scale.dat", "ta0")
levels(S$expr) = c("1D", "2D", "3D", "1D-rl", "2D-rl", "3D-rl")
s = ggplot(S, aes(x=factor(scale),y=w,fill=expr)) + geom_boxplot() + labs(x="Scaling factor",y="Fitness") + theme(legend.position="right") + quick_theme
showfig(s)
savefig(s, "p-scaling-fitness", width=12)