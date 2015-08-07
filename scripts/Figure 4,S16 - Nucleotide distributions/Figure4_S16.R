require("sfsmisc")
require("tikzDevice")

mononucs <- c("T", "G", "C", "A")
dinucs <- c("TT", "TG", "TC", "TA", "GT", "GG", "GC", "GA", "CT", "CG", "CC", "CA", "AT", "AG", "AC", "AA")

cols <- list(A = "green", C = "dodgerblue2", G = "orange", T = "red")
mono_cols <- c()

RGBcols <- list()
for (i in mononucs) {
    RGBcols[[i]] <- col2rgb(cols[[i]])
    mono_cols <- c(mono_cols, cols[[i]])
}

p_mix <- 0.8
dinuc_cols <- c()
for (j in 1:16) {
    i <- dinucs[j]
    rgb_col <- round(p_mix * RGBcols[[substr(i, start = 1, stop = 1)]] + (1 - p_mix) * RGBcols[[substr(i, start = 2, stop = 2)]])
    
    dinuc_cols[[j]] <- rgb(rgb_col[1, 1], rgb_col[2, 1], rgb_col[3, 1], max = 255)
}

## 0. load raw data
files <- c("3223a", "3223b", "3223c", "3236a", "3236b", "3236c")

DATA <- list()
for (i in files) {
    DATA[[i]] <- read.table(paste("../Analysis/", i, "_nucMask_2_pID_counts.csv", sep = ""), header = TRUE, sep = ",", colClasses = c("character", "numeric", "numeric"))
}
DATA[["Total"]] <- read.table("../Analysis/Total_pID_counts.csv", header = TRUE, sep = ",", colClasses = c("character", "numeric", "numeric"))


## 1. marginal (mono-)nucleotide barplots
freq_single_uniq <- list()
freq_single_repl <- list()

for (i in c(files, "Total")) {
    message("file: ", i)
    
    temp_DATA <- DATA[[i]]
    message("Total sum: ", sum(temp_DATA$replicates))
    
    mat_single_uniq <- mat.or.vec(nr = 4, nc = 10)
    mat_single_repl <- mat.or.vec(nr = 4, nc = 10)
    n <- dim(temp_DATA)[1]
    
    for (j in 1:10) {
        message("j: ", j)
        
        stats_uniq <- list()
        stats_repl <- list()
        
        for (k in mononucs) {
            stats_uniq[[k]] <- 0
            stats_repl[[k]] <- 0
        }
        
        for (l in 1:n) {
            if (l%%5000 == 0) {
                message(l)
            }
            
            stats_uniq[[substr(temp_DATA[l, ]$Primer, start = j, stop = j)]] <- stats_uniq[[substr(temp_DATA[l, ]$Primer, start = j, stop = j)]] + temp_DATA[l, ]$unique
            stats_repl[[substr(temp_DATA[l, ]$Primer, start = j, stop = j)]] <- stats_repl[[substr(temp_DATA[l, ]$Primer, start = j, stop = j)]] + temp_DATA[l, ]$replicates
        }
        
        for (k in 1:4) {
            mat_single_uniq[k, j] <- stats_uniq[[mononucs[k]]]
            mat_single_repl[k, j] <- stats_repl[[mononucs[k]]]
        }
        
        mat_single_uniq[, j] <- mat_single_uniq[, j]/sum(mat_single_uniq[, j])
        mat_single_repl[, j] <- mat_single_repl[, j]/sum(mat_single_repl[, j])
    }
    
    freq_single_uniq[[i]] <- mat_single_uniq
    freq_single_repl[[i]] <- mat_single_repl
}


## 2. pair-wise dinucleotide barplots
DATA_TOTAL <- DATA[["Total"]]

# loci of interest
loci <- mat.or.vec(nr = 2, nc = 9)
loci[1, ] <- 1:9
loci[2, ] <- loci[1, ] + 1

# summary statistic
n <- dim(DATA_TOTAL)[1]
freq_pairwise <- mat.or.vec(nr = 4 * 4, nc = 9 * 2)
i <- 1

for (j in 1:9) {
    message("j: ", j)
    locus <- loci[, j]
    
    message("start: ", locus[1])
    message("stop:  ", locus[2])
    
    stats_uniq <- list()
    stats_repl <- list()
    
    for (k in dinucs) {
        stats_uniq[[k]] <- 0
        stats_repl[[k]] <- 0
    }
    
    for (l in 1:n) {
        i
        if (l%%5000 == 0) {
            message(l)
        }
        
        stats_uniq[[substr(DATA_TOTAL[l, ]$Primer, start = locus[1], stop = locus[2])]] <- stats_uniq[[substr(DATA_TOTAL[l, ]$Primer, start = locus[1], stop = locus[2])]] + DATA_TOTAL[l, ]$unique
        stats_repl[[substr(DATA_TOTAL[l, ]$Primer, start = locus[1], stop = locus[2])]] <- stats_repl[[substr(DATA_TOTAL[l, ]$Primer, start = locus[1], stop = locus[2])]] + DATA_TOTAL[l, ]$replicates
    }
    
    for (k in 1:16) {
        freq_pairwise[k, i] <- stats_uniq[[dinucs[k]]]
        freq_pairwise[k, i + 1] <- stats_repl[[dinucs[k]]]
    }
    
    freq_pairwise[, i] <- freq_pairwise[, i]/sum(freq_pairwise[, i])
    freq_pairwise[, i + 1] <- freq_pairwise[, i + 1]/sum(freq_pairwise[, i + 1])
    
    i <- i + 2
}

offset <- 0.4
scale <- 1
x <- rep(c(1, 1 + offset), times = 9) + rep((0:8) * scale, each = 2)


## 3.) actual plot A.) mononucleotide block
CEX.withinbars <- 0.35
CEX.axis <- 0.9
SCALE.withinbars <- 1
SCALE.axis <- 1

pdf(file = "Figure4.pdf", width = 10, height = 6.5)
par(oma = c(2.5, 2.5, 1.2, 1.5), mar = c(1, 0.25, 1, 0.25), mgp = c(3, 0.5, 0), xpd = NA, las = 1, lty = 0, lwd = 0.1)
layout(matrix(c(1, 3, 5, 7, 9, 11, 2, 4, 6, 8, 10, 12, 16, 13, 13, 14, 14, 15, 17, 13, 13, 14, 14, 15), nrow = 4, byrow = TRUE))

for (i in c(files, "Total")) {
    # 1) unique barplots
    if (i == "Total") {
        par(mar = c(1, 0.25, 4.1, 0.25), mgp = c(3, 1, 0), las = 1, lty = 0, lwd = 0.1)
        SCALE.withinbars <- 2
        SCALE.axis <- 1.5
    }
    
    X <- barplot(freq_single_uniq[[i]], space = c(0.1), axes = FALSE, col = mono_cols)
    if (i == "3223a" || i == "Total") {
        eaxis(2, at = (0:5)/5, labels = paste((0:5) * 20, "% ", sep = ""), cex.axis = CEX.axis * SCALE.axis)
    }
    if (i == "3236c") {
        text("unique", x = par("usr")[2] * 1.08, y = (par("usr")[3] + par("usr")[4])/2, srt = 90, font = 2, cex = 1.5)
    }
    axis(1, at = X[seq(from = 1, to = 10, by = 2)], labels = seq(from = 1, to = 10, by = 2), cex.axis = CEX.axis * SCALE.axis)
    axis(1, at = X[seq(from = 2, to = 10, by = 2)], labels = seq(from = 2, to = 10, by = 2), cex.axis = CEX.axis * SCALE.axis)
    if (i == "Total") {
        mtext("Position", side = 1, line = 2.5, cex = CEX.axis * SCALE.axis/1.1)
        mtext("Total unique", side = 3, line = 0.4, font = 2, cex = 1.2)
    } else mtext(i, side = 3, line = 0.4, font = 2, cex = 1.2)
    
    for (j in 1:10) {
        low <- 0
        high <- 0
        
        for (k in 1:4) {
            low <- high
            high <- high + freq_single_uniq[[i]][k, j]
            y_target <- (low + high)/2
            text(X[j], y_target, paste(formatC(freq_single_uniq[[i]][k, j] * 100, format = "f", 1), "%", sep = ""), cex = CEX.withinbars * SCALE.withinbars)
        }
    }
    
    
    # 2) PCR replicates barplots
    X <- barplot(freq_single_repl[[i]], space = c(0.1), axes = FALSE, col = mono_cols)
    if (i == "3223a") {
        eaxis(2, at = (0:5)/5, labels = paste((0:5) * 20, "% ", sep = ""), cex.axis = CEX.axis * SCALE.axis)
    }
    if (i == "3236c") {
        text("PCR weighted", x = par("usr")[2] * 1.08, y = (par("usr")[3] + par("usr")[4])/2, srt = 90, font = 2, cex = 1.5)
    }
    axis(1, at = X[seq(from = 1, to = 10, by = 2)], labels = seq(from = 1, to = 10, by = 2), cex.axis = CEX.axis * SCALE.axis)
    axis(1, at = X[seq(from = 2, to = 10, by = 2)], labels = seq(from = 2, to = 10, by = 2), cex.axis = CEX.axis * SCALE.axis)
    mtext("Position", side = 1, line = ifelse(i == "Total", 2.5, 1.7), cex = CEX.axis * SCALE.axis/ifelse(i == "Total", 1.1, 1))
    if (i == "Total") 
        mtext("Total PCR weighted", side = 3, line = 0.4, font = 2, cex = 1.2)
    
    for (j in 1:10) {
        low <- 0
        high <- 0
        
        for (k in 1:4) {
            low <- high
            high <- high + freq_single_repl[[i]][k, j]
            y_target <- (low + high)/2
            text(X[j], y_target, paste(formatC(freq_single_repl[[i]][k, j] * 100, format = "f", 1), "%", sep = ""), cex = CEX.withinbars * SCALE.withinbars)
        }
    }
}
legend(x = par("usr")[2], y = par("usr")[4], legend = rev(mononucs), fill = rev(mono_cols), box.lwd = 1, box.lty = 0, cex = 2)
segments(grconvertX(0.5, from = "nic"), grconvertY(0.48, from = "nic"), grconvertX(0.5, from = "nic"), grconvertY(0.98, from = "nic"), lwd = 0.5, lty = 1)
dev.off()

# B.) dinucleotide plot
CEX.perc <- 0.45

tikz("FigureS16.tex", standAlone = FALSE, height = 4, width = 6.5)
par(oma = c(0, 0, 0, 0), mar = c(5.5, 3, 0.4, 3.8), xpd = TRUE, las = 1)
X <- barplot(freq_pairwise, space = c(1, 0.15), axes = FALSE, col = dinuc_cols)
legend(x = X[16] + 4.5, y = 1, legend = rev(dinucs), fill = rev(dinuc_cols))
eaxis(2, at = (0:10)/10, labels = paste("$", formatC((0:10)/10 * 100, format = "f", 0), "\\%$", sep = ""), lab.type = "latex", small.mult = 2)

x_ticks <- c()
x_tick_labels <- c()
for (i in 1:9) {
    x_ticks <- c(x_ticks, mean(c(X[2 * i - 1], X[2 * i])))
    x_tick_labels <- c(x_tick_labels, paste("\\textbf{", loci[1, i], "--", loci[2, i], "}", sep = ""))
}

text(X, par("usr")[3], labels = rep(c("unique", "PCR weighted"), 9), srt = 60, adj = c(1, 1), xpd = TRUE, cex = 0.8)
text(x_ticks, par("usr")[3] - 0.27 * (par("usr")[4] - par("usr")[3]), labels = x_tick_labels, xpd = TRUE, cex = 1.5)

for (i in 1:(2 * 9)) {
    low <- 0
    high <- 0
    for (j in 1:16) {
        low <- high
        high <- high + freq_pairwise[j, i]
        
        y_target <- (low + high)/2
        
        text(X[i], y_target, paste("$", formatC(freq_pairwise[j, i] * 100, format = "f", 1), "\\%$", sep = ""), cex = CEX.perc)
    }
}