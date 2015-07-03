require(plotrix)
require(RColorBrewer)
colours = brewer.pal(7, "Set1")
colours[6] = "gold"

DATA = read.table("qPCR.csv", header=TRUE, sep=",")

CONTROL_OFFSETS = c(3.65, 3.8, 3.95)
DILUTION_OFFSET = 1:3
PRIMER_RANGE = 0.45
PRIMER_OFFSETS = c(-PRIMER_RANGE/2, PRIMER_RANGE/2)
MINOR_OFFSET = seq(from=-0.13, to=0.13, length.out=4)

x_axis_ticks = c(rep(DILUTION_OFFSET, each=2) + rep(PRIMER_OFFSETS, 3), 3.8)

DILUTION=c("1_10", "1_100", "1_1000")
PRIMER=c("A", "B")
SUBTYPE=c("1", "2", "3", "4")

PCHs=c(21, 22, 23, 24)-6
CEX=c(1, 1.1, 1, 1.3)*1.2

pdf(file="Figure2.pdf", width=8, height=5)
par(oma = c(0,0,0,0))
par(mar = c(4,4,0.5,0.5))
par(mgp = c(2.5,1,0))

plot(NA, xlim=c(0.7, 4.0), ylim=c(10, 30), xaxt="n", yaxt="n", ylab=expression(C["t"]), xlab="", cex.lab=1.4)
axis(1, at=x_axis_ticks, labels=c(rep(c("3223", "3236"), 3), "Controls"), cex.axis=1.2)
axis(2, at=c(10, 15, 20, 25, 30), labels=c(10, 15, 20, 25, expression(infinity)), las=1, cex.axis=1.2)
axis.break(2, 27.5, style="zigzag")

for (i in 1:3)
{	# dilution series
	for (j in 1:2)
	{	# primer type
		for (k in 1:4)
		{	# primer subtype
			y=subset(DATA, Dilution==DILUTION[i] & Primer==PRIMER[j] & Type==SUBTYPE[k])$Ct_Mean
			points(DILUTION_OFFSET[i] + PRIMER_OFFSETS[j] + MINOR_OFFSET[k], y, pch=PCHs[k], col=colours[k], cex=CEX[k])
		}
	}
}

Control_Primer = c("A", "B", "0")
Control_Type = c("RT", "RT", "Water")
Dil_Labels = c("1:10", "1:100", "1:1000")
for (i in 1:3)
{
	mtext(Dil_Labels[i], at=DILUTION_OFFSET[i], side=1, line=2.5, font=2, cex=1.3)
	
	y=subset(DATA, Dilution=="Control" & Primer==Control_Primer[i] & Type==Control_Type[i])$Ct_Mean
	points(CONTROL_OFFSETS[i], y, col=colours[i+4], pch=16, cex=CEX[2])
}

legend("bottomright", legend=c("RT", "RT_ID", "RT_A_ID", "RT_J_ID", "Neg. Ctrl. 3223 no RT", "Neg. Ctrl. 3236 no RT", "Neg. Ctrl. Water"), col = colours, pch = c(PCHs, rep(16, 3)), pt.cex=c(CEX, rep(CEX[2], 3)))
dev.off()
