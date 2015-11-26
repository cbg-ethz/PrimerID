files <- c("3223a", "3223b", "3223c", "3236a", "3236b", "3236c")

# 1) histogram data
DATA_hetero_Hamming_distance <- list()
max_dist = 0

y_max <- 0
for (i in 1:6) {
    file <- files[i]
    message(file)
    
    DATA_hetero_Hamming_distance[[i]] <- read.table(paste("../Analysis/", file, "_nucMask_2_min_hetero_hamming_histograms.csv", sep = ""), sep = ",", header = TRUE, colClasses = "numeric")
    max_dist = max(max_dist, min(which(DATA_hetero_Hamming_distance[[i]]$Count == 0))-1)
}


# 2) table for minimum heterozygous Hamming distance distribution
hetero_table = function() {
	message("\\begin{tabular}{r", rep("c", max_dist), "}")
	message("\\toprule")

	line = "$d_H$"
	for (i in 1:max_dist) {
		line = paste(line, " & ", i - 1, sep = "")
	}
	message(line, "\\\\")
	message("\\midrule")

	for (i in 1:6) {
		const_norm = sum(DATA_hetero_Hamming_distance[[i]]$Count)

		line = files[i]
		for (j in 1:max_dist) {
			value = DATA_hetero_Hamming_distance[[i]]$Count[j] * 100/const_norm

			if (value < 0.1) {
				if (value == 0) 
					line = paste(line, " & $0\\%$", sep = "")
				else line = paste(line, " & $\\nprounddigits{1}\\numprint{\\xintFloat[2]{", value, "}}\\%$", sep = "")
			} else {
				line = paste(line, " & $\\nprounddigits{1}\\numprint{", value, "}\\%$", sep = "")
			}
		}

		message(line, "\\\\[\\myrowskip]")
	}
	message("\\bottomrule")
	message("\\end{tabular}")
}
hetero_table()

hetero_table_abridged = function() {
	message("\\begin{tabular}{rcc}")
	message("\\toprule")

	line = " & $d_H \\leq 2$ & $d_H > 2$\\\\"
	message(line)
	message("\\midrule")

	for (i in 1:6) {
		const_norm = sum(DATA_hetero_Hamming_distance[[i]]$Count)

		line = files[i]

		value = sum(DATA_hetero_Hamming_distance[[i]]$Count[1:3]) * 100/const_norm

		if (value < 0.1) {
			if (value == 0) 
				line = paste(line, " & $0\\%$", sep = "")
			else line = paste(line, " & $\\nprounddigits{1}\\numprint{\\xintFloat[2]{", value, "}}\\%$", sep = "")
		} else {
			line = paste(line, " & $\\nprounddigits{1}\\numprint{", value, "}\\%$", sep = "")
		}

		value = sum(DATA_hetero_Hamming_distance[[i]]$Count[4:max_dist]) * 100/const_norm

		if (value < 0.1) {
			if (value == 0) 
				line = paste(line, " & $0\\%$", sep = "")
			else line = paste(line, " & $\\numprint{\\xintFloat[2]{", value, "}}\\%$", sep = "")
		} else {
			line = paste(line, " & $\\nprounddigits{1}\\numprint{", value, "}\\%$", sep = "")
		}


		message(line, "\\\\[\\myrowskip]")
	}
	message("\\bottomrule")
	message("\\end{tabular}")
}
hetero_table_abridged()
