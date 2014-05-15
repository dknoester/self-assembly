library(ggplot2)
library(plyr)
library(stringr)
library(extrafont)
loadfonts(quiet=TRUE)
Sys.setenv(R_GSCMD = "/opt/local/bin/gs")
source("~/research/src/self-assembly/math/common.R")
setwd("/Users/dk/research/src/self-assembly/var")
#figpath = "/Users/dk/research/doc/self-assembly/saso2014/non-embedded/"
figpath = "/Users/dk/Dropbox/goldsby-doc/saso2014/non-embedded/"


STYLE="draft" # or "final"
WIDTH=6
HEIGHT=3.75

quick_theme <- theme_bw() + theme(panel.border=element_blank(), axis.line=element_line(colour = "black",size=0.75), legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(1,0), text=element_text(size=14))


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
	quartz(width=WIDTH, height=HEIGHT)
	print(x)
}

# Save figure x in in PDF format.
#
savefig <- function(x, name, width=WIDTH, height=HEIGHT) {
	f = paste(figpath,name,".pdf",sep="")
	pdf(file=f, width=width, height=height, family="Helvetica")
	print(x)
	dev.off()
#	embed_fonts(f)#, outfile=paste(figpath,name,"-embed.pdf",sep=""))
}

# Fitness data for the different experiments.
#
#
x002 = load.files("fitness.dat.gz", "ta0", "002-1d-fsm")
x003 = load.files("fitness.dat.gz", "ta0", "003-2d-fsm")
x004 = load.files("fitness.dat.gz", "ta0", "004-3d-fsm")

# Fitness of base experiments.
F_base_fitness = expr_fitness(rbind(x002,x003,x004))
F_base_fitness = F_base_fitness + scale_fill_discrete(breaks=c("002-1d-fsm","003-2d-fsm","004-3d-fsm"),labels=c("1D","2D","3D")) + scale_color_discrete(breaks=c("002-1d-fsm","003-2d-fsm","004-3d-fsm"),labels=c("1D","2D","3D"))
showfig(F_base_fitness)
savefig(F_base_fitness, "p-fitness")

x008 = load.files("fitness.dat.gz", "ta0", "008-1d-fsm-rand")
x009 = load.files("fitness.dat.gz", "ta0", "009-2d-fsm-rand")
x010 = load.files("fitness.dat.gz", "ta0", "010-3d-fsm-rand")

# Fitness of reinforcement learning experiments.
F_rl_fitness = expr_fitness(rbind(x008, x009, x010))
F_rl_fitness = F_rl_fitness + scale_fill_discrete(breaks=c("008-1d-fsm-rand","009-2d-fsm-rand","010-3d-fsm-rand"),labels=c("1D-rl","2D-rl","3D-rl")) + scale_color_discrete(breaks=c("008-1d-fsm-rand","009-2d-fsm-rand","010-3d-fsm-rand"),labels=c("1D-rl","2D-rl","3D-rl"))
showfig(F_rl_fitness)
savefig(F_rl_fitness, "p-fitness-rl")

x011 = load.files("fitness.dat.gz", "ta0", "011-1d-fsm-switch")
x012 = load.files("fitness.dat.gz", "ta0", "012-2d-fsm-switch")
x013 = load.files("fitness.dat.gz", "ta0", "013-3d-fsm-switch")

F_switch_fitness = expr_fitness(rbind(x011, x012, x013))
F_switch_fitness = F_switch_fitness + scale_fill_discrete(breaks=c("011-1d-fsm-switch","012-2d-fsm-switch","013-3d-fsm-switch"),labels=c("1D","2D","3D")) + scale_color_discrete(breaks=c("011-1d-fsm-switch","012-2d-fsm-switch","013-3d-fsm-switch"),labels=c("1D","2D","3D"))
showfig(F_switch_fitness)
savefig(F_switch_fitness, "p-fitness-reevolved")

# 1000x dominant fitness
#
#
D = load.files("ca_dom_1000x")
D = subset(D, expr=="002-1d-fsm" | expr=="003-2d-fsm" | expr=="004-3d-fsm")
levels(D$expr) = c("1D", "2D", "3D", "1D-rl", "2D-rl", "3D-rl")
# D$delta_w = D$w1 - D$w0
d = ggplot(D, aes(x=expr,y=w1,fill=expr)) + geom_boxplot() + geom_jitter() + ylim(0,1) + quick_theme + ylab("Fitness (w)") + theme(legend.position="none") + xlab("")
showfig(d)
savefig(d, "p-1000x-fitness")

# D1 = subset(D,expr=="1D")
# D1[which.max(D1$w1),]
# individual     w0    w1                              filename expr treatment trial delta_w
# 6     494817 0.9918 0.865 .//002-1d-fsm/ta0_14/ca_dom_1000x.dat   1D       ta0    14 -0.1268

# D1 = subset(D,expr=="2D")
# D1[which.max(D1$w1),]
# individual     w0    w1                              filename expr treatment trial delta_w
# 52     498839 0.9955 0.884 .//003-2d-fsm/ta0_29/ca_dom_1000x.dat   2D       ta0    29 -0.1115

# D1 = subset(D,expr=="3D")
# D1[which.max(D1$w1),]
# individual     w0    w1                              filename expr treatment trial delta_w
# 63     333564 0.9918 0.881 .//004-3d-fsm/ta0_11/ca_dom_1000x.dat   3D       ta0    11 -0.1108

# Adaptive data
#
#
A = load.files("ca_dom_adapt.dat")
A$delta_w = A$w0 - A$w1
levels(A$expr) = c("1D-rl", "2D-rl", "3D-rl")
a = ggplot(A, aes(expr,delta_w,fill=expr)) + theme(legend.position="none") + labs(x="", y="Fitness delta") + geom_boxplot() + geom_jitter() + ylim(0,1) + quick_theme
showfig(a)
savefig(a, "p-adaptive-fitness")

# KO hidden data
K = load.files("ca_dom_ko_hidden.dat")
summary(K)
K$delta_w = K$w1 - K$w0
levels(K$expr) = c("1D", "2D", "3D")
k = ggplot(K, aes(x=expr,y=delta_w,fill=expr)) + geom_boxplot() + geom_jitter() + ylim(-1,0) + quick_theme + ylab(expression(paste(Delta, plain(w)))) + theme(legend.position="none") + xlab("")

showfig(k)
savefig(k, "p-ko-hidden")


# Scaling data
#
#
S = load.files("ca_dom_scale.dat", "ta0")
S = subset(S, expr=="002-1d-fsm" | expr=="003-2d-fsm" | expr=="004-3d-fsm")
levels(S$expr) = c("1D", "2D", "3D", "1D-rl", "2D-rl", "3D-rl")
s = ggplot(S, aes(x=factor(scale),y=w,color=expr)) + geom_point() + labs(x="Scaling factor",y="Fitness") + quick_theme + theme(legend.position="right") + geom_jitter()
showfig(s)
savefig(s, "p-scaling-fitness", width=12)

# S1 = subset(S, expr=="1D" & scale==9)
# S1[which.max(S1$w),]
# scale    w                              filename expr treatment trial
# 54     9 0.86 .//002-1d-fsm/ta0_14/ca_dom_scale.dat   1D       ta0    14

# S1 = subset(S, expr=="2D" & scale==9)
# S1[which.max(S1$w),]
# scale   w                              filename expr treatment trial
# 441     9 0.9 .//003-2d-fsm/ta0_26/ca_dom_scale.dat   2D       ta0    26

# S1 = subset(S, expr=="3D" & scale==9)
# S1[which.max(S1$w),]
# scale    w                              filename expr treatment trial
# 747     9 0.93 .//004-3d-fsm/ta0_30/ca_dom_scale.dat   3D       ta0    30

# Rule table density data
#
#
R = load.files("ca_dom_rule_density")
levels(R$expr) = c("1D")
R = merge(R, subset(S,expr=="1D" & scale==1), by="trial")
# cor(R$w, R$rho, method="spearman")
# == 0.106

P = load.files("ca_dom_sampled_rule_density")
levels(P$expr) = c("2D", "3D")
P = rbind(data.frame(w=P$w, rho=P$rho, expr=P$expr), data.frame(w=R$w, rho=R$rho, expr=R$expr.x))
r = ggplot(P, aes(x=w,y=rho,color=expr)) + geom_point() + quick_theme  + xlim(c(0,1)) + ylim(c(0,1)) + theme(legend.position="right")
showfig(r)
savefig(r, "p-rule-density")