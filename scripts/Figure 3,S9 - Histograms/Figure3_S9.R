files <- c("3223a", "3223b", "3223c", "3236a", "3236b", "3236c")

# 1) histogram data
DATA_abundance <- list()
DATA_pID_length <- list()
y_max <- 0
for (i in 1:6) {
    file <- files[i]
    message(file)
    
    DATA_abundance[[i]] <- read.table(paste("../Analysis/", file, "_nucMask_2_abundance_histograms.csv", sep = ""), sep = ",", header = TRUE, colClasses = "numeric")
    DATA_pID_length[[i]] <- read.table(paste("../Analysis/", file, "_nucMask_2_pID_length_histograms.csv", sep = ""), sep = ",", header = TRUE, colClasses = "numeric")
    y_max <- max(y_max, DATA_pID_length[[i]]$Count)
}

# 2) plot abundance histograms
pdf(file = "Figure3.pdf", width = 3 * 5, height = 2 * 5)
par(mfrow = c(2, 3))
par(oma = c(0, 0.2, 0, 0))
par(mgp = c(2.7, 1, 0))

for (i in 1:6) {
    file <- files[i]
    message(file)
    
    plot(DATA_abundance[[i]]$Abundance, DATA_abundance[[i]]$Count, log = "y", xaxs = "i", type = "h", lwd = 5, lend = 1, col = "gray", xlab = "Abundance", ylab = "Count", main = file, cex.lab = 1.8, cex.axis = 1.3, cex.main = 2.5, xlim = c(-2, max(DATA_abundance[[i]]$Abundance) + 10), xaxt = "n")
    axis(1, at = c(1, seq(from = 50, to = max(DATA_abundance[[i]]$Abundance) + 10, by = 50)), cex.axis = 1.3)
    
    Y_MIN <- 10^par("usr")[3]
    Y_MAX <- 10^par("usr")[4]
    
    segments(10, Y_MIN, 10, Y_MAX/3, col = "red", lwd = 1.2, lty = 2, lend = 1)
}
dev.off()

# 3) plot pID length distribution
require("RColorBrewer")
require("sfsmisc")
require("tikzDevice")

# pdf(file='FigureS9.pdf', width=10, height=5)
tikz(file = "FigureS9.tex", standAlone = FALSE, height = 3.5, width = 6.5)
colours <- brewer.pal(6, "Set1")
colours[6] <- "gold"

par(oma = c(0, 0, 0, 0))
par(mar = c(3.8, 4, 0.5, 0.5))
par(mgp = c(2.6, 1, 0))
OFFSET_FROM_MIDDLE <- 0.35
OFFSET <- (2 * OFFSET_FROM_MIDDLE)/5

plot(NA, log = "y", type = "h", lwd = 5, lend = 1, col = "black", xlab = "pID length", ylab = "Count", main = "", cex.lab = 1.8, cex.axis = 1.2, cex.main = 2.5, xlim = c(0, 20), ylim = c(1, y_max), axes = FALSE)
box()

for (i in 1:6) {
    lines(DATA_pID_length[[i]]$Length - OFFSET_FROM_MIDDLE + OFFSET * (i - 1), DATA_pID_length[[i]]$Count, type = "h", lwd = 4.5, lend = "square", col = colours[i])
    message(paste("Proportion of 10-mers in", files[i], ":", DATA_pID_length[[i]]$Count[11]/sum(DATA_pID_length[[i]]$Count)))
}

legend("topright", files, col = colours, lwd = 3)
eaxis(1, at = 0:20, f.smalltcl = 0, lab.type = "latex")
eaxis(2, lab.type = "latex") 