###############################################################################
# CG159 Social Anxiety: Contrast Level Threshold Analysis
#   - Performed on standardised mean difference scale (SMD)
###############################################################################

# Install the nmathresh package from CRAN if not already installed
# install.packages("nmathresh")

library(nmathresh)

# The nmathresh package includes the following relevant data objects which are
# used for this analysis:
#   - SocAnx.post.summary, the posterior summaries of the variables
#   - SocAnx.post.cov, the posterior covariance matrix
#
# The raw study data is also included in the nmathresh package, which we
# partially read in for convenience below.
#
# For details of how these are created from WinBUGS output, see 
# browseVignettes("nmathresh").


## RECONSTRUCT LIKELIHOOD  ----------------------------------------------------
# First, we must reconstruct the likelihood covariance matrix.

# Create contrast design matrix X. We'll use the treatment details from the 
# study-level data for convenience.

trt.dat <- read.delim(system.file("extdata", "SocAnx_data.txt", package = "nmathresh"))[, 1:6]
K <- with(trt.dat, length(unique(c(t.1, t.2, t.3, t.4, t.5)))) -1  # Number of treatments , -1 for NA

# Work out which contrasts have data
contr.ab <- data.frame(a = c(), b = c())

for (i in 1:nrow(trt.dat)) {
   rowi <- trt.dat[i, which(!is.na(trt.dat[i, 2:6])) + 1] # non NA elements of ith row

   # get contrast from all combinations of treatments
   trtcomb <- combn(rowi, 2, function(x) sapply(x, as.numeric))

   a <- apply(trtcomb, 2, min)
   b <- apply(trtcomb, 2, max)

   # remove contrasts of treatments against themselves
   iseq <- a == b
   a <- a[!iseq]
   b <- b[!iseq]

   if (!all(iseq)) contr.ab <- rbind(contr.ab, cbind(a, b))
}

contr.ab <- unique(contr.ab[order(contr.ab$a, contr.ab$b), ])


# Contrast design matrix
X <- matrix(0, nrow = nrow(contr.ab), ncol = K - 1)
for (i in 1:nrow(X)) {
  if (contr.ab[i, "a"] > 1) X[i, contr.ab[i, "a"] - 1]  <- -1
  if (contr.ab[i, "b"] > 1)   X[i, contr.ab[i, "b"] - 1]    <- 1
}


# Reconstruct using NNLS
# Class model is used, so use the prior mean precision from the gamma distribution
prior.prec <- rep(3.9/0.35, 40)
# Other than for treatment 3
prior.prec[2] <- 0.0001

lik.cov <- recon_vcov(SocAnx.post.cov[1:(K - 1), 1:(K - 1)],  # Posterior covariance matrix
					  prior.vcov = diag(1/prior.prec),    # Prior covariance matrix
					  X = X)                              # Contrast design matrix


# Get indices of contrasts in likelihood
d.a <- d.b <- vector(length = nrow(X))
for (i in 1:nrow(X)) {
  d.a[i] <- ifelse(any(X[i, ] == -1), which(X[i, ] == -1), 0) + 1
  d.b[i] <- ifelse(any(X[i, ] == 1), which(X[i, ] == 1), 0) + 1
}

d.i <- d_ab2i(d.a, d.b, K = K)


## THRESHOLD ANALYSIS ---------------------------------------------------------
# Now we can perform threshold analysis at the contrast level.

# Get indices of basic treatment parameters d and contrasts diff in the CODA data
vnames <- sub("(.*)\\[.*","\\1", rownames(SocAnx.post.summary$statistics))
ind.d <- which(vnames == "d")
ind.sd <- which(vnames == "sd")
ind.diff <- which(vnames == "diff")
ind.delta <- which(vnames == "delta")

# Thresholds are
post.mean.d <- SocAnx.post.summary$statistics[ind.d, "Mean"]
post.cov.d <- SocAnx.post.cov[1:(K - 1), 1:(K - 1)]


thresh <- nma_thresh(mean.dk = post.mean.d,   # Posterior means of treatment effects
                     lhood = lik.cov,            # Reconstructed likelihood covariance matrix
                     post = post.cov.d,  # Posterior covariance matrix
                                                 # -- Further options below --
                     nmatype = "fixed",          # Approximate two-stage FE NMA
                     X = X,                      # Contrast design matrix
                     opt.max = FALSE)            # Lower treatment effects better


## PLOT RESULTS ---------------------------------------------------------------

# Label the contrasts and display using a forest plot, along with 95% posterior
# credible intervals as reported in CG150.1

plotdat <- data.frame(
  contr = paste0(d.b, " vs. ", d.a),
  contr.mean = SocAnx.post.summary$statistics[ind.diff[d.i], "Mean"],
  CI2.5 = SocAnx.post.summary$quantiles[ind.diff[d.i], "2.5%"],
  CI97.5 = SocAnx.post.summary$quantiles[ind.diff[d.i], "97.5%"]
)

pdf("Social Anxiety thresholds - contrast level.pdf", width = 12.5, height = 24.5)
thresh_forest(thresh,           # Threshold object produced by nma_thresh
              contr.mean,       # Means, CIs, and labels in data frame
              CI2.5, CI97.5, 
              label = contr,  
              data = plotdat,   # Data frame containing above data
			                          # -- Plotting options below here --
              label.title = "Contrast", 
			        xlab = "Standardised Mean Difference", xlim = c(-4, 3),
              y.title = "SMD", CI.title = "95% Credible Interval",
              II.cols = rgb(0.72, 0.8, 0.93),
              refline = 0, digits = 2)
dev.off()


# Only report contrasts where the threshold is less than 2 SMD.
# We'll use the orderby parameter of thresh_forest for this, creating an order
# variable with NAs for those rows to be removed, and using na.last=NA.

cutoff <- 2
plotdat$thresh.cut <- ifelse(thresh$thresholds$lo > -cutoff | thresh$thresholds$hi < cutoff,
                             1:nrow(plotdat), NA)

pdf("Social Anxiety thresholds - contrast level CUTOFF.pdf", width = 12.5, height = 10.5)
thresh_forest(thresh, contr.mean, CI2.5, CI97.5,
              label = contr, orderby = list(thresh.cut, na.last = NA), data = plotdat,
              label.title = "Contrast", xlab = "Standardised Mean Difference", xlim = c(-4, 3),
              y.title = "SMD", CI.title = "95% Credible Interval",
              II.cols = rgb(0.72, 0.8, 0.93),
              refline = 0, digits = 2)
dev.off()


# Now sorted with smallest thresholds first

absmin <- function(x) min(abs(x))
plotdat$absmin.thresh <- apply(thresh$thresholds[, c("lo", "hi")], 1, absmin)
plotdat$thresh.cut.sort <- ifelse(thresh$thresholds$lo > -cutoff | thresh$thresholds$hi < cutoff, 
                                  plotdat$absmin.thresh, NA)

pdf("Social Anxiety thresholds - contrast level CUTOFF SORT.pdf", width = 12.5, height = 10.8)
fp <- thresh_forest(thresh, contr.mean, CI2.5, CI97.5,
              label = contr, orderby = list(thresh.cut.sort, na.last = NA), data = plotdat,
              label.title = "Contrast", xlab = "Standardised Mean Difference", xlim = c(-4, 3),
              y.title = "SMD", CI.title = "95% Credible Interval",
              refline = 0, digits = 2, display = FALSE)

# Add optimal treatment caption below table
nrow_fp <- nrow(fp)
fp <- gtable::gtable_add_rows(fp, grid::unit(1, "null"))
fp <- gtable::gtable_add_grob(fp, 
                              grid::textGrob(paste0("Base-case recommended treatment is ", 
                                                    thresh$kstar, "."), 
                                             x = 0, just = "left"),
                              t = nrow_fp + 1, l = 4, clip = "off")
fp <- gtable::gtable_add_rows(fp, grid::unit(1, "null"))

# Display
grid::grid.newpage()
grid::grid.draw(fp)

dev.off()


## PSYCH TREATMENT BIAS  ------------------------------------------------------
# Here we investigate a potential common bias for all psychological treatments 
# against inactive.

# Psych treatments
psychtrts <- c(4:8, 24:36)

# Which data points compare these to an inactive trt?
psychdats <- which(contr.ab$a %in% 1:3 & contr.ab$b %in% psychtrts)

# Get U solutions by summing the individual solutions of drug data points
U.psych <- 1 / (rowSums(1 / thresh$Ukstar[,psychdats]))

# Which contrasts do the rows of Ukstar correspond to?
Ukstar.ab <- d_i2ab(1:(K * (K - 1) / 2), K)
Ukstar.ab <- Ukstar.ab[Ukstar.ab$a == thresh$kstar | Ukstar.ab$b == thresh$kstar, ]

# Thresholds are then
thresh.psych <- get.int(U.psych, thresh$kstar, 1:K, Ukstar.ab)

# Function to plot the common invariant interval with the data
plotII <- function(thresh, contr.mean, CrI.lo, CrI.hi, rowlabs, xlim, xlab, ylab, ...){

  yaxs <- length(contr.mean):1

  # split plot in two
  layout(matrix(1:2,nrow = 1), widths = c(.2, 1))

  # plot row labels on left side
  gp <- par(mar = c(5, 4, 1, 0))
  plot(rep(0,length(yaxs)), yaxs, pch = NA, ylim = c(.5,yaxs[1] + .5), ylab = "",
       xlab = "", yaxt = "n", xaxt = "n", bty = "n")
  text(0, yaxs, labels = rowlabs,xpd = NA)
  title(ylab = ylab, line = 2)

  # fake plot for setup of right side
  par(mar = c(5, 1, 1, 2))
  plot(contr.mean, yaxs, pch = NA, yaxt = "n", xlim = xlim,
       ylim = c(.5,yaxs[1] + .5), ylab = "", xlab = "", ...)
  title(xlab = xlab, line = 2)

  # reference line
  abline(v = 0, lty = 2, col = "gray")

  # combined invariant region
  polygon(rep(c(contr.mean + thresh$lo, rev(contr.mean) + thresh$hi), each = 2),
          c(rep(yaxs,each = 2) + c(.5, -.5), rep(rev(yaxs),each = 2) + c(-.5, .5)),
          col = rgb(.7, .8, .9, .7),
          border = rgb(.6, .7, .8))

  # credible intervals
  segments(y0 = yaxs, x0 = CrI.lo, x1 = CrI.hi, lend = 1)

  # contrast means
  points(contr.mean, yaxs, pch = 21, col = "black", bg = "white")

  # write invariant interval below plot
  with(thresh, title(xlab = paste0("Invariant interval about zero: ",lo.newkstar,
                                 " (",formatC(lo,digits = 2,format = "f"),", ",
                                 formatC(hi,digits = 2,format = "f"),") ",
                                 hi.newkstar), line = 3))

  par(gp)
}

# Plot
pdf("Social Anxiety thresholds - contrast level by group - psych.pdf", width = 8, height = 5.5)
plotII(thresh.psych, 
       SocAnx.post.summary$statistics[ind.diff[psychdats], "Mean"], 
       SocAnx.post.summary$quantiles[ind.diff[psychdats], "2.5%"],
       SocAnx.post.summary$quantiles[ind.diff[psychdats], "97.5%"],
       rowlabs = paste0(contr.ab[psychdats,"b"]," vs. ",contr.ab[psychdats,"a"]),
       xlim = c(-3, 2), ylab = "Psych vs. Inactive Contrasts", xlab = "SMD")
dev.off()

