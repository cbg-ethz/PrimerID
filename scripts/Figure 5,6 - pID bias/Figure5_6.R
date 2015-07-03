RANK_CUTOFF = 100
INTERSECT_CUTOFF = 500

points_reg = 50
SEQ = unique(round(10^seq(from=0, to=log10(RANK_CUTOFF), length.out=points_reg)))

files = c("3223a", "3223b", "3223c", "3236a", "3236b", "3236c")

total_N = c(430506, 166684, 446435, 515358, 522257, 533962)
PIS=c("0.4000", "0.4246", "0.4492", "0.4738", "0.4983", "0.5229", "0.5475", "0.5721", "0.5967", "0.6213", "0.6458", "0.6704", "0.6950", "0.7196", "0.7442", "0.7688", "0.7933", "0.8179", "0.8425", "0.8671", "0.8917", "0.9163", "0.9408", "0.9654", "0.9900")
NRT=c(3000, 3014, 3028, 3042, 3056, 3070, 3084, 3098, 3113, 3127, 3142, 3156, 3171, 3185, 3200, 3215, 3230, 3245, 3260, 3275, 3290, 3305, 3321, 3336, 3351, 3367, 3382, 3398, 3414, 3430, 3445, 3461, 3477, 3493, 3510, 3526, 3542, 3559, 3575, 3592, 3608, 3625, 3642, 3658, 3675, 3692, 3709, 3727, 3744, 3761, 3779, 3796, 3814, 3831, 3849, 3867, 3885, 3903, 3921, 3939, 3957, 3975, 3994, 4012, 4031, 4049, 4068, 4087, 4106, 4125, 4144, 4163, 4182, 4202, 4221, 4241, 4260, 4280, 4300, 4320, 4340, 4360, 4380, 4400, 4420, 4441, 4461, 4482, 4503, 4524, 4544, 4565, 4587, 4608, 4629, 4651, 4672, 4694, 4715, 4737, 4759, 4781, 4803, 4825, 4848, 4870, 4893, 4915, 4938, 4961, 4984, 5007, 5030, 5053, 5077, 5100, 5124, 5147, 5171, 5195, 5219, 5243, 5268, 5292, 5316, 5341, 5366, 5391, 5415, 5440, 5466, 5491, 5516, 5542, 5567, 5593, 5619, 5645, 5671, 5697, 5724, 5750, 5777, 5804, 5830, 5857, 5884, 5912, 5939, 5966, 5994, 6022, 6050, 6078, 6106, 6134, 6162, 6191, 6219, 6248, 6277, 6306, 6335, 6365, 6394, 6424, 6453, 6483, 6513, 6543, 6574, 6604, 6635, 6665, 6696, 6727, 6758, 6789, 6821, 6852, 6884, 6916, 6948, 6980, 7012, 7045, 7077, 7110, 7143, 7176, 7209, 7242, 7276, 7310, 7343, 7377, 7412, 7446, 7480, 7515, 7550, 7584, 7620, 7655, 7690, 7726, 7762, 7797, 7833, 7870, 7906, 7943, 7979, 8016, 8053, 8091, 8128, 8166, 8203, 8241, 8279, 8318, 8356, 8395, 8434, 8473, 8512, 8551, 8591, 8631, 8670, 8711, 8751, 8791, 8832, 8873, 8914, 8955, 8997, 9038, 9080, 9122, 9164, 9207, 9249, 9292, 9335, 9378, 9421, 9465, 9509, 9553, 9597, 9641, 9686, 9731, 9776, 9821, 9866, 9912, 9958, 10004, 10050, 10097, 10143, 10190, 10237, 10285, 10332, 10380, 10428, 10476, 10525, 10573, 10622, 10671, 10721, 10770, 10820, 10870, 10921, 10971, 11022, 11073, 11124, 11175, 11227, 11279, 11331, 11384, 11436, 11489, 11542, 11596, 11649, 11703, 11757, 11812, 11866, 11921, 11976, 12032, 12087, 12143, 12199, 12256, 12313, 12369, 12427, 12484, 12542, 12600, 12658, 12717, 12776, 12835, 12894, 12954, 13014, 13074, 13134, 13195, 13256, 13317, 13379, 13441, 13503, 13565, 13628, 13691, 13754, 13818, 13882, 13946, 14011, 14076, 14141, 14206, 14272, 14338, 14404, 14471, 14538, 14605, 14672, 14740, 14808, 14877, 14946, 15015, 15084, 15154, 15224, 15295, 15365, 15436, 15508, 15579, 15652, 15724, 15797, 15870, 15943, 16017, 16091, 16165, 16240, 16315, 16391, 16466, 16543, 16619, 16696, 16773, 16851, 16929, 17007, 17086, 17165, 17244, 17324, 17404, 17484, 17565, 17647, 17728, 17810, 17893, 17975, 18058, 18142, 18226, 18310, 18395, 18480, 18565, 18651, 18738, 18824, 18911, 18999, 19087, 19175, 19264, 19353, 19442, 19532, 19622, 19713, 19804, 19896, 19988, 20080, 20173, 20267, 20360, 20454, 20549, 20644, 20740, 20836, 20932, 21029, 21126, 21224, 21322, 21420, 21520, 21619, 21719, 21819, 21920, 22022, 22124, 22226, 22329, 22432, 22536, 22640, 22745, 22850, 22956, 23062, 23168, 23276, 23383, 23491, 23600, 23709, 23819, 23929, 24040, 24151, 24263, 24375, 24488, 24601, 24715, 24829, 24944, 25059, 25175, 25291, 25408, 25526, 25644, 25763, 25882, 26001, 26122, 26242, 26364, 26486, 26608, 26731, 26855, 26979, 27104, 27229, 27355, 27482, 27609, 27737, 27865, 27994, 28123, 28253, 28384, 28515, 28647, 28780, 28913, 29046, 29181, 29316, 29451, 29588, 29724, 29862, 30000)

# 1) load simulation data
DATA = list()
for (i in 1:6)
{
	file = files[i]
	Nreads = total_N[i]
	message(file)
	
	for (p in PIS)
	{
		message(p)
		for (number_rt in NRT)
		{
			source(paste("Simulations/", file, "/Nrt_", number_rt, "_Nreads_", Nreads, "_p_", p, ".txt", sep=""))
		}
	}
}

# 2) load sequencing data
real_DATA = list()
for (i in 1:6)
{
	file = files[i]
	
	real_DATA[[i]] = read.table(file=paste("../Analysis/", file, "_nucMask_2_primers_collision_free.csv", sep=""), header=TRUE, sep=",")
}

# 3) perform OLS
best_match = list()
for (i in 1:6)
{
	file = files[i]
	Nreads = total_N[i]
	message(file)
	
	tmp_DATA = real_DATA[[i]]$Count[order(real_DATA[[i]]$Count, decreasing=TRUE)]
	
	min = 100000000000000
	for (p in PIS)
	{
		message(p)
		for (number_rt in NRT)
		{
			id = paste("Nrt_", number_rt, "_Nreads_", Nreads, "_p_", p, sep="")
			sq = sum((log(tmp_DATA[SEQ]) - log(DATA[[id]][SEQ]))^2)
			
			if (sq < min)
			{
				min = sq
				best_match[[i]] = id
			}
		}
	}
}

# 4) plot data
pdf(file="Figure6.pdf", width=3*5, height=2*5)
par(mfrow=c(2,3))
par(oma=c(1,1,1,1))
par(mgp=c(2.6,1,0))

for (i in 1:6)
{
	file = files[i]
	Nreads = total_N[i]
	message(file)
	
	tmp_DATA = real_DATA[[i]]$Count[order(real_DATA[[i]]$Count, decreasing=TRUE)]
	
	plot(1:RANK_CUTOFF, DATA[[best_match[[i]]]][1:RANK_CUTOFF], log="xy", xlab="Rank", ylab="Abundance", type="l", main=file, cex.lab=1.8, cex.axis=1.2, cex.main=2.5, ylim=c(min(tmp_DATA[1:RANK_CUTOFF]), max(tmp_DATA[1:RANK_CUTOFF])), lwd=2.5)
	
	points(1:RANK_CUTOFF, tmp_DATA[1:RANK_CUTOFF])
}
dev.off()

# 5) Venn diagrams
require("VennDiagram")
require(RColorBrewer)
colours = brewer.pal(5, "Set1")
colours[5] = "gold"

for (i in 1:6)
{
	ORDER_i = order(real_DATA[[i]]$Count, decreasing=TRUE)
	seqs_i = as.vector(real_DATA[[i]]$Primer[ORDER_i[1:INTERSECT_CUTOFF]])
	
	UNION = c()
	for (j in 1:6)
	{
		if (i != j)
		{
			ORDER_j = order(real_DATA[[j]]$Count, decreasing=TRUE)
			seqs_j = as.vector(real_DATA[[j]]$Primer[ORDER_j[1:INTERSECT_CUTOFF]])
			
			UNION = c(UNION, seqs_j)
		}
	}
	
	#message("Complement Union size ", i, " : ", length(unique(UNION)))
	message("Intersection of ", files[i], " : ", length(intersect(seqs_i, unique(UNION))))
}

# collect sets
SEQUENCES_3223 = list()
SEQUENCES_3236 = list()
SEQUENCES_ALL = list()
for (i in 1:6)
{
	file = files[i]
	ORDER = order(real_DATA[[i]]$Count, decreasing=TRUE)
	message("Processing ", file)
	message("Rank 1  :", real_DATA[[i]]$Count[ORDER[1]])
	message("Rank 100:", real_DATA[[i]]$Count[ORDER[INTERSECT_CUTOFF]])
	message("Range   :", real_DATA[[i]]$Count[ORDER[1]] - real_DATA[[i]]$Count[ORDER[INTERSECT_CUTOFF]])
	
	if (file != "3223b")
	{
		SEQUENCES_ALL[[file]] = real_DATA[[i]]$Primer[ORDER[1:INTERSECT_CUTOFF]]
	}
	
	if ((file == "3223a") || (file == "3223b") || (file == "3223c"))
	{
		SEQUENCES_3223[[file]] = real_DATA[[i]]$Primer[ORDER[1:INTERSECT_CUTOFF]]
	}
	else
	{
		SEQUENCES_3236[[file]] = real_DATA[[i]]$Primer[ORDER[1:INTERSECT_CUTOFF]]
	}
}

# plot
ALPHA = 0.45
CEX.LAB = 1.5
CEX.LAB5 = 1.7
CEX.CLASS = 2

BASE_UNIT = 2
pdf(file = "Figure5.pdf", width = 6*BASE_UNIT, height=7*BASE_UNIT)
grid.newpage()
pushViewport(plotViewport(c(0, 1, 0, 0)))
pushViewport(viewport(layout = grid.layout(7, 6)))
	pushViewport(viewport(layout.pos.col=c(1,3), layout.pos.row=c(1,3)))
		VENN_3223 = venn.diagram(SEQUENCES_3223, filename = NULL, euler.d = FALSE, scaled = FALSE, fill = colours[1:3], alpha = rep(ALPHA, 3), cex = CEX.CLASS, cat.cex = CEX.CLASS, cat.dist = c(0.055, 0.055, 0.035), margin = 0.04)
		grid.draw(VENN_3223)
		grid.text("A", x=0.01, y=0.95, gp=gpar(fontsize=35, fontface="bold"))
	popViewport(1)
	
	pushViewport(viewport(layout.pos.col=c(4,6), layout.pos.row=c(1,3)))
		VENN_3236 = venn.diagram(SEQUENCES_3236, filename = NULL, euler.d = FALSE, scaled = FALSE, fill = colours[1:3], alpha = rep(ALPHA, 3), cex = CEX.CLASS, cat.cex = CEX.CLASS, cat.dist = c(0.055, 0.055, 0.035), margin = 0.04)
		grid.draw(VENN_3236)
		grid.text("B", x=0.01, y=0.95, gp=gpar(fontsize=35, fontface="bold"))
	popViewport(1)
	
	pushViewport(viewport(layout.pos.col=c(2,5), layout.pos.row=c(4,7)))
		VENN_ALL = venn.diagram(SEQUENCES_ALL, filename = NULL, euler.d = FALSE, scaled = FALSE, fill = colours[1:5], alpha = rep(ALPHA, 5), cex = CEX.LAB5, cat.cex = CEX.CLASS, cat.dist = c(0.18, 0.23, 0.18, 0.20, 0.22))
		grid.draw(VENN_ALL)
		grid.text("C", x=0.25, y=0.95, gp=gpar(fontsize=35, fontface="bold"))
	popViewport(1)
popViewport()
popViewport()
dev.off()