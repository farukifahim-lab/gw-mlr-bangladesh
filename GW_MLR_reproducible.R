############################################################
# Interpretable MLR equations for groundwater level (GWL) in Bangladesh
# Reproducible script (relative paths; run from the repository root)
# Input : data/Final-ALL-MSD.csv   (15,000 station-month records)
# Output: outputs/*.csv, outputs/*.png
# NOTE  : cleaned from the original analysis script; please re-run in R
#         and confirm the numbers against the manuscript before archiving.
############################################################
suppressPackageStartupMessages({
  library(caTools); library(hydroGOF); library(car); library(lmtest)
})
dir.create("outputs", showWarnings = FALSE)
set.seed(123)

# ---------- 1. Load and clean ----------
dat <- read.csv("data/Final-ALL-MSD.csv", check.names = FALSE)
dat <- dat[, !(names(dat) %in% c("", "X"))]          # drop the row-index column (it is not a predictor)
names(dat) <- make.names(names(dat))                  # e.g. Mean-GWL -> Mean.GWL
# GWETPROF / GWETROOT contain the missing-value code -999 in 20 records: remove them
miss <- dat$GWETPROF < -900 | dat$GWETROOT < -900
cat("Records removed for -999 codes:", sum(miss), "\n")
dat <- dat[!miss & complete.cases(dat), ]
stopifnot(all(sapply(dat, is.numeric)))

# ---------- 2. Train-test split (70/30) ----------
split <- sample.split(dat$GWL, SplitRatio = 0.7)
train <- dat[split, ]; test <- dat[!split, ]

rmse <- function(o, p) sqrt(mean((o - p)^2))
r2   <- function(o, p) cor(o, p)^2

# ---------- 3. Full MLR baseline ----------
full <- lm(GWL ~ ., data = train)
aliased <- names(which(is.na(coef(full))))            # perfectly collinear terms (e.g. ro = sro + ssro)
if (length(aliased)) {
  cat("Aliased terms dropped:", paste(aliased, collapse = ", "), "\n")
  keep <- setdiff(names(train), c(aliased, "GWL"))
  full <- lm(reformulate(keep, "GWL"), data = train)
}
print(gof(predict(full, test), test$GWL))

# ---------- 4. Remove multicollinearity (VIF < 10) ----------
mod <- full
repeat {
  v <- vif(mod)
  if (max(v) < 10) break
  drop <- names(which.max(v)); cat("Removing", drop, "VIF =", round(max(v), 1), "\n")
  vars <- setdiff(attr(terms(mod), "term.labels"), drop)
  mod  <- lm(reformulate(vars, "GWL"), data = train)
}
summary(mod)
print(bptest(mod)); print(ncvTest(mod)); print(durbinWatsonTest(mod))

# ---------- 5. Standardized beta ranking ----------
cf <- coef(mod)[-1]
beta <- data.frame(Variable = names(cf),
                   Beta = cf * sapply(train[, names(cf)], sd) / sd(train$GWL))
beta <- beta[order(-abs(beta$Beta)), ]
write.csv(beta, "outputs/Beta_Importance_MLR.csv", row.names = FALSE)

# ---------- 6. Top-N models ----------
res <- data.frame(); preds <- list(); models <- list()
for (n in c(30, 25, 20, 15, 10, 5)) {
  n <- min(n, nrow(beta))
  vars <- beta$Variable[1:n]
  m <- lm(reformulate(vars, "GWL"), data = train)
  p <- predict(m, test)
  models[[paste0("Top", n)]] <- m; preds[[paste0("Top", n)]] <- p
  res <- rbind(res, data.frame(N = n, RMSE = rmse(test$GWL, p), R2 = r2(test$GWL, p)))
}
print(res); write.csv(res, "outputs/TopN_performance.csv", row.names = FALSE)

# ---------- 7. Equations ----------
for (k in c("Top15", "Top10", "Top5")) {
  if (is.null(models[[k]])) next
  cf <- round(coef(models[[k]]), 4)
  cat("\n", k, ": GWL = ", cf[1], paste0(ifelse(cf[-1] >= 0, " + ", " - "), abs(cf[-1]), " * ", names(cf)[-1], collapse = ""), "\n", sep = "")
}

# ---------- 8. Extremes (10th / 90th percentile of test GWL) ----------
lo <- quantile(test$GWL, 0.10); hi <- quantile(test$GWL, 0.90)
ext <- do.call(rbind, lapply(names(preds), function(k) {
  p <- preds[[k]]; o <- test$GWL
  do.call(rbind, lapply(list(Low = o <= lo, High = o >= hi), function(ix)
    data.frame(Model = k, N = sum(ix), RMSE = rmse(o[ix], p[ix]),
               Bias = mean(p[ix] - o[ix]), R2 = r2(o[ix], p[ix]))))
}))
write.csv(ext, "outputs/Extreme_performance.csv", row.names = TRUE)

# ---------- 9. Top-5 equation: evaluation on the full sample ----------
top5 <- 77.2957 + 1.0968*dat$Mean.GWL + 0.1204*dat$T2M_MAX - 0.2698*dat$stl4 -
        70.2097*dat$ssro - 2.3285*dat$swvl3
print(gof(top5, dat$GWL))
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo.txt")
