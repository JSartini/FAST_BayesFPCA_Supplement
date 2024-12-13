# Non-standard packages potentially requiring installation
if(!require(splines2)){
  install.packages("splines2")
}

if(!require(furrr)){
  install.packages("furrr")
}

if(!require(abind)){
  install.packages("abind")
}

if(!require(nleqslv)){
  install.packages("nleqslv")
}

if(!require(RiemBase)){
  install.packages("RiemBase")
}

# Calculation packages
library(tidyverse)
library(splines2)
library(rstan)
library(pracma)
library(refund)
library(furrr)
library(mgcv)
library(abind)
library(lme4)
library(nleqslv)
library(RiemBase)
