require(lattice)
require(grid)
require(gridExtra)

source("seqLogo/pwm.R")
source("seqLogo/seqLogo.R")
source("seqLogo/AllClasses.R")
source("seqLogo/AllGenerics.R")

files = c("3223a", "3223b", "3223c", "3236a", "3236b", "3236c")
types = c("_nucMask_2_unique", "_nucMask_2_replicates")

FONTSIZE = 9
X_FONTSIZE = FONTSIZE
Y_FONTSIZE = FONTSIZE

pdf(file = "Figure4.pdf", width = 6*1.6, height=4*1.6)
grid.newpage()
pushViewport(plotViewport(c(2, 2.4, 1, 1)))
pushViewport(viewport(layout = grid.layout(4, 6)))
	vp = viewport(width=1.0, height=1.05)
	I = 0
	for(i in types)
	{
		I = I + 1
		J = 0
		for(j in files)
		{
			J = J + 1
			DATA_raw = read.table(paste("../Analysis/", j, i, ".csv", sep=""), sep=",")
			DATA = makePWM(DATA_raw)
			
			pushViewport(viewport(layout.pos.col=J, layout.pos.row=I))
				if (J == 1)
					YAXIS=TRUE
				else
					YAXIS=FALSE
				
				if (I == 2)
					XLAB=TRUE
				else
					XLAB=FALSE
				
				seqLogo(DATA, xaxis = TRUE, yaxis=YAXIS, xlab=XLAB, xfontsize=X_FONTSIZE, yfontsize=Y_FONTSIZE, ic.scale=FALSE)
				
				if (J == 3)
					grid.border(type=3, colour="black", vp=vp)
				if (J == 4)
					grid.border(type=5, colour="black", vp=vp)
			popViewport(1)
			
			if (I == 1)
				grid.text(j, y = 1.02, vp = viewport(layout.pos.col=J, layout.pos.row=I), gp=gpar(font=2, fontsize=15))
			if (J == 6)
				grid.text(ifelse(I == 1, "unique", "PCR weighted"), x = 1.02, vp = viewport(layout.pos.col=J, layout.pos.row=I), rot = 90, gp=gpar(font=2, fontsize=12))
		}
	}
	
	LARGE_SCALE_FACTOR = 1.3
	
	# plot total unique RT replicates
	DATA_raw = read.table("../Analysis/Total_unique.csv", sep=",")
	DATA = makePWM(DATA_raw)
	
	pushViewport(viewport(layout.pos.col=c(2,3), layout.pos.row=c(3,4)))
		seqLogo(DATA, xaxis = TRUE, yaxis=TRUE, xlab=TRUE, xfontsize=X_FONTSIZE*LARGE_SCALE_FACTOR, yfontsize=Y_FONTSIZE*LARGE_SCALE_FACTOR, ic.scale=FALSE, margins = c(1, 0.5, 2.4, 0.5))
	popViewport(1)
	grid.text("Total unique", y = 0.85, vp = viewport(layout.pos.col=c(2,3), layout.pos.row=c(3,4)), gp=gpar(font=2, fontsize=15))
	
	# plot total PCR replicates
	DATA_raw = read.table("../Analysis/Total_replicates.csv", sep=",")
	DATA = makePWM(DATA_raw)
	
	pushViewport(viewport(layout.pos.col=c(4,5), layout.pos.row=c(3,4)))
		seqLogo(DATA, xaxis = TRUE, yaxis=FALSE, xlab=TRUE, xfontsize=X_FONTSIZE*LARGE_SCALE_FACTOR, yfontsize=Y_FONTSIZE*LARGE_SCALE_FACTOR, ic.scale=FALSE, margins = c(1, 0.5, 2.4, 0.5))
	popViewport(1)
	grid.text("Total PCR weighted", y = 0.85, vp = viewport(layout.pos.col=c(4,5), layout.pos.row=c(3,4)), gp=gpar(font=2, fontsize=15))

popViewport()
popViewport()
dev.off()

# # Plots for CROI Slides
# pdf(file = "TotalMotifs.pdf", width = 2*4, height=2.5)
# grid.newpage()
# pushViewport(plotViewport(c(2, 3, 0, 0)))
# pushViewport(viewport(layout = grid.layout(1, 2)))
	# LARGE_SCALE_FACTOR = 1.3
	
	# # plot total unique RT replicates
	# DATA_raw = read.table("../Analysis/Total_unique.csv", sep=",")
	# DATA = makePWM(DATA_raw)
	
	# pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
		# seqLogo(DATA, xaxis = TRUE, yaxis=TRUE, xlab=TRUE, xfontsize=X_FONTSIZE*LARGE_SCALE_FACTOR, yfontsize=Y_FONTSIZE*LARGE_SCALE_FACTOR, ic.scale=FALSE, margins = c(1, 0.5, 1.5, 0))
	# popViewport(1)
	# grid.text("Total unique", y = 0.92, vp = viewport(layout.pos.col=1, layout.pos.row=1), gp=gpar(font=2, fontsize=20))
	
	# # plot total PCR replicates
	# DATA_raw = read.table("../Analysis/Total_replicates.csv", sep=",")
	# DATA = makePWM(DATA_raw)
	
	# pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
		# seqLogo(DATA, xaxis = TRUE, yaxis=FALSE, xlab=TRUE, xfontsize=X_FONTSIZE*LARGE_SCALE_FACTOR, yfontsize=Y_FONTSIZE*LARGE_SCALE_FACTOR, ic.scale=FALSE, margins = c(1, 0.5, 1.5, 0))
	# popViewport(1)
	# grid.text("Total PCR weighted", y = 0.92, vp = viewport(layout.pos.col=2, layout.pos.row=1), gp=gpar(font=2, fontsize=20))
# popViewport()
# popViewport()
# dev.off()
