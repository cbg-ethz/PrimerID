require(tikzDevice)
require(sfsmisc)

source("../Analysis/LogLikData.R")

tikz("HMM_LogLik.tex", standAlone = TRUE, height = 3)
par(mar = c(2.5, 6.3, 0.7, 0.0))
par(xpd = TRUE)

plot(x, y, type="n", log="x", main = "", axes=FALSE, xlab="", ylab="", xlim=c(2E-6, max(x)))
lines(x, y, type="l")

# x-axis
eaxis(1, lab.type="latex", lab.sep="cdot")
mtext("$r$", side=1, line = 1.5, cex = 1.5)

# y-axis
eaxis(2, lab.type="latex")
mtext("Log-Likelihood $\\mathcal{L}$", side=2, line = 5, cex = 1.5)

YMIN = par("usr")[3]
YMAX = par("usr")[4]

LOGLIKMAX = max(y)
LOGLIKMAX_CI = LOGLIKMAX -1.920729

lnX_CI_LOW = log(X_CI_LOW)
lnX_CI_HIGH = log(X_CI_HIGH)
lnMIDP = mean(c(lnX_CI_LOW, lnX_CI_HIGH))
MIDP = exp(lnMIDP)

XSPAN = lnX_CI_HIGH - lnX_CI_LOW
ratio = 0.28
x_left = exp(lnX_CI_LOW + ratio*XSPAN)
x_right = exp(lnX_CI_LOW + (1-ratio)*XSPAN)

segments(X_CI_LOW, YMIN, X_CI_LOW, LOGLIKMAX_CI)
segments(X_CI_HIGH, YMIN, X_CI_HIGH, LOGLIKMAX_CI)

segments(r, YMIN, r, YMIN + 0.91*(LOGLIKMAX_CI-YMIN), lty=3)
segments(r, YMIN + 1.00*(LOGLIKMAX_CI-YMIN), r, LOGLIKMAX, lty=3)

EXP = floor(log10(r))
BASE = signif(r, 3) / 10^EXP

text(r, YMIN + 0.95*(YMAX-YMIN), paste("$\\hat{r} = $", pretty10exp(r, lab.type = "latex", lab.sep = "cdot")), pos = 3, cex = 1.0)

arrows(X_CI_LOW, YMIN + 0.95*(LOGLIKMAX_CI-YMIN), x_left, YMIN + 0.95*(LOGLIKMAX_CI-YMIN), code = 1, length = 0.07)
arrows(x_right, YMIN + 0.95*(LOGLIKMAX_CI-YMIN), X_CI_HIGH, YMIN + 0.95*(LOGLIKMAX_CI-YMIN), code = 2, length = 0.07)

text(x= MIDP, y=YMIN + 0.95*(LOGLIKMAX_CI-YMIN), labels="95\\% CI")
dev.off()