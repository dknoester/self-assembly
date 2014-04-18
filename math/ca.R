library(ggplot2)
source("~/research/src/self-assembly/math/common.R")
# source("~/research/src/self-assembly/math/fixup_ggplot.R")
setwd("/Users/dk/research/src/self-assembly/var")

STYLE="draft" # or "final"

quick_theme <- theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(colour = "black",size=0.75), legend.title=element_blank()) 


fitness_dplot <- function(D,path) {
	if(STYLE=="draft") {
		D = subset(D,update%%100==0)		
	}
	
	quartz(width=6,height=3.75)
	g = ggplot(data=D, aes(x=update, y=max_fitness)) + stat_summary(aes(color=treatment,fill=treatment),fun.data="mean_cl_boot", geom="smooth") + labs(x="Update", y="Fitness") + quick_theme + ylim(0,1) #+ xlim(0,10000)
	
	if(STYLE=="draft") {
		g =	g + ggtitle(path)
	}
	
	print(g)
	return(list(D,g))
}

fitness_plot <- function(path) {
	return(fitness_dplot(rbind(load.files(find.files("fitness.dat.gz",path=path))), path))
}

## 012-1d-fsm-rl-switch
#
# Here we're taking an RL learning system, and re-evolving it for a different objective;
# initial tests show that this is almost meaningless, as the RL signal appears to be strong enough
# to guide the system entirely.
x012 = fitness_plot("012-1d-fsm-rl-switch")

## 011-1d-fsm-switch
#
# Here we're taking a non-RL system, and re-evolving it for the opposite objective.  This addresses
# the question of how quickly can a self-adaptive system re-evolve for a new objective, vs. a 
# non-self-adaptive system.
x011 = fitness_plot("011-1d-fsm-switch")

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

## 004-3d-fsm
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


## majority, no RL
dx002 = subset(x002[[1]],treatment=="tb0")
dx002$treatment = "002-tb0"

dx003 = subset(x003[[1]],treatment=="tb0")
dx003$treatment = "003-tb0"

dx004 = subset(x004[[1]],treatment=="tb0")
dx004$treatment = "004-tb0"

D = rbind(dx002, dx003, dx004)

p001 = fitness_dplot(D,"p001")


## majority, no RL
t = "tc0"
dx002 = subset(x002[[1]],treatment==t)
dx002$treatment = "002-tb0"

dx003 = subset(x003[[1]],treatment==t)
dx003$treatment = "003-tb0"

dx004 = subset(x004[[1]],treatment==t)
dx004$treatment = "004-tb0"

D = rbind(dx002, dx003, dx004)

#p001 = fitness_dplot(D,"p001-123d-prob")
p002 = fitness_dplot(D,"p001-123d-adapt")


