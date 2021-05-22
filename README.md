# EntropyProductionFromWaitingTimes
This repository provides the code from the paper "Estimating entropy production from waiting times", arXiv:2105.08681.

The scripts in this directory will reproduce the figures and tables from the paper and SI, using saved data from Data/. The Mathematica notebook AsymptoticCalc.nb shows the derivation of the asymptotic formula.

The script ExampleScript.m shows how easy it is to bound the entropy production rate from experimental measurements. It also shows how the saved data was generated, although this was performed on the MIT Supercloud cluster, taking many computational hours -- the saved data means you don't need to repeat this computation!

For convenience we include some of the data we used in the paper in Data/. These are fully referenced in the paper, and the original work should be cited if used. In most cases, the data was simply extracted from histograms directly from the article or supplementary PDF, with the exact figure referenced in Table S1. 
