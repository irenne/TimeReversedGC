# TimeReversedGC

This package contains Matlab code for testing time-reversed Granger causality (TRGC) [1,2] as described in [1]. 

#### Files 

runSimulation.m    -  Start this to run a sample simulation and produce a plot similar to Figure 1b in [1]. 
create2dVAR        -   Function which generates a two-dimensional time series from a random autoregressive process 
timeReversedGC  -  Function which performs the following statistical tests: 
					- difference-based TRGC (based on bootstrapping) 
					- conjunction-based TRGC (based on bootstrapping) 
					- net Granger causality (based on bootstrapping) 
					- standard Granger causality  (F-test)

#### References

[1] I. Winkler, D. Panknin, D. Bartz, K.-R. MŸller, S. Haufe (2016). Validity of time reversal for testing Granger causality. IEEE Transactions on Signal Processing, 64(11), 2746-2760.

[2] S. Haufe, V.V. Nikulin, K.-R. MŸller, K. R., G. Nolte (2013). A critical assessment of connectivity measures for EEG data: a simulation study. Neuroimage, 64, 120-133.

