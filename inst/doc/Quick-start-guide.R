## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(generalKSStat)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 50L
x <- rbeta(n, 1, 1.5)

## plot the empirical vs null CDF
plot(ecdf(x))
abline(a = 0,b = 1)

## KS stat
GKSStat(x = x, statName = "KS")

## HC stat
GKSStat(x = x, statName = "HC")

## BJ stat
GKSStat(x = x, statName = "BJ")

## -----------------------------------------------------------------------------
index <- 25L : 50L

## KS stat
GKSStat(x = x, index = index, statName = "KS")

## HC stat
GKSStat(x = x, index = index, statName = "HC")

## BJ stat
GKSStat(x = x, index = index, statName = "BJ")

## -----------------------------------------------------------------------------
indexL <- 25L : 50L
indexU <- NULL

## One-sided KS+ stat
GKSStat(x = x, indexL = indexL, indexU = NULL, statName = "KS")

## One-sided HC+ stat
GKSStat(x = x, indexL = indexL, indexU = NULL, statName = "HC")

## One-sided BJ+ stat
GKSStat(x = x, indexL = indexL, indexU = NULL, statName = "BJ")

## -----------------------------------------------------------------------------
index <- 25L : 50L

## One-sided KS+ stat
GKSStat(x = x, index = index, statName = "KS+")

## One-sided HC+ stat
GKSStat(x = x, index = index, statName = "HC+")

## One-sided BJ+ stat
GKSStat(x = x, index = index, statName = "BJ+")

## -----------------------------------------------------------------------------
sessionInfo()

