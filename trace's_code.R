#setwd("Z:/Practicum")
setwd("~/Documents/Practicum")
load("predictors.Rdata")
library(fields)
library(maps)
library(mvtnorm)
library(gstat)
library(geoR)
library(sp)
library(maptools)

long = lon-360