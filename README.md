# HumanSignsofLife
## Created on 4/23/2021


This repo contains all of the code used in the Human Signs of Life project.

Currently the SIMULATIONS folder is the only folder
Inside is all of the Matlab code that I have been writing. 

The files that matter are 

### mcWrapper.m
  Main Simulation for ALL radar Tests
### singleMCTrials.m
  Runs a single Radar simulation for parameters given
### estFreqs.m
  Estimates the frequencies in a given chest compression signal
### getChestCompression.m
  Creates a random chest compression signal with given parameters
### getRespFilter.m
  Provides the filter coefficients for the Respiratory filter
### getHeartFilter.m
  Provides the filter coefficients for the heart rate filter
### easyFFT.m
  Function that calculates an FFT for a given input. 
