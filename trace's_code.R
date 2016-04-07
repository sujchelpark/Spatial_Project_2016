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

lon.new = lon-360