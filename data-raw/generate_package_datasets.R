## (re)generate the mice and crohn datasets for the package

# generate mice.rda ---------------------------------------------------------------------

library(coin)
library(nleqslv)
library(xtable)

dat = read.table("transformed_raw_data_STAN.csv", sep = ",", header = TRUE) ## place the file in the working directory.
Lab1 = dat[dat$lab == "Gies_Wur", ]
Lab2 = dat[dat$lab == "Man_Gas", ]
dim(Lab1)
#head(Lab1)
dim(Lab2)
#head(Lab2)
p1 = rep(NA, 29)
p1L = rep(NA, 29)
p1R = rep(NA, 29)
isL1 = rep(NA, 29)
p2 = rep(NA, 29)
p2L = rep(NA, 29)
p2R = rep(NA, 29)

isL2 = rep(NA, 29)
b = 1
for (i in 3:31) {
  Y = Lab1[, i]
  X = factor(Lab1$strain)
  print(wilcox_test(Y ~ X, distribution = "exact"))
  out1 = wilcox_test(Y ~ X, distribution = "exact")
  isL1[b] = statistic(out1) < 0 # get the direction
  p1[b]   = pvalue(out1) / 2 # get the one sided p-value
  p1L[b]  = ifelse(statistic(out1) < 0, pvalue(out1) / 2, 1 - pvalue(out1) / 2) # transform the right directed pvalues to left ones
  p1R[b]  = 1 - p1L[b] # get the right complements

  Y       = Lab2[, i]
  X       = factor(Lab2$strain)
  print(wilcox_test(Y ~ X, distribution = "exact"))
  out2 = wilcox_test(Y ~ X, distribution = "exact")
  isL2[b] = statistic(out2) < 0
  p2[b] = pvalue(out2) / 2
  p2L[b] = ifelse(statistic(out2) < 0, pvalue(out2) / 2, 1 - pvalue(out2) / 2)
  p2R[b] = 1 - p2L[b]
  b = b + 1
}

## create a dataframe to store the example for the package
mice <- data.frame(feature_name = names(Lab1)[3:31], twosided_pv1 = 2*p1, twosided_pv2 = 2*p2, dir_is_left1 = isL1, dir_is_left2 = isL2)
devtools::use_data(mice,overwrite = T)


# generate crohn.rda ---------------------------------------------------------------------

pv <- read.csv("http://www.math.tau.ac.il/~ruheller/Software/CrohnExample.csv")
crohn <- data.frame(index = pv$Index,
                    pv1 = pv$p1,
                    pv2 = pv$p2)
devtools::use_data(crohn, overwrite = T)

