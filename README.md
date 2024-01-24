# BPDNWP_test

*reference manuscript, dates, authors, and copywrites

Dependencies:
- CVX Matlab package: http://cvxr.com/cvx/
- Matlab Wavelet toolbox: https://www.mathworks.com/products/wavelet.html

Processing
- Analyses were run on a computer with Intel® Xeon® Silver 4108 CPU @ 1.80 GHz, 188 GB RAM, 16 physical cores, and 32 logical cores, and running Ubuntu version 18.04.6

We recommend working through the scripts in this order:
1)	GenSimSpikeData_script
2)	ParameterSelection_script
3)	RunModel_script

Contents:
1. GenSimSpikeData_script: Script to generate and plot simulated spike rate timeseries 
2. wavelet_decomp_L2norm: Function that runs wavelet transform and calculate L2 norm
3. BPDN_fullSetup: Wrapper for functions that setup inputs to method (sampling vector, bases, etc.)
4. BPDN_setupdata: Function to format samples and timestamps for method
5. frequency_sampling: Function to define frequencies represented in DCT basis
6. DCT2_basis: Function to create DCT-II basis based on N and defined frequencies 
7. polynomial_basis: Function to create polynomial basis based on N and maximum desired polynomial degree
8. BPDN: Function to calculate basis pursuit denoising with polynomial trend 
9. BPDN_reconsig: Function to reconstruct timedomain signal based on BDPN outputs
10. ParameterSelection_script: Script to run delta parameter sweeps, plot outputs, update data with delta
11. BPDN_samplerealdata_traintest: Re-sample continuous timeseries to get sparse training sets
12. RunModel_script: Script to run BPDN on simulated data with optimal delta parameter, test output for significance, plot overall outputs
13. BPDN_wreshuffling_short: Function to calculate percentiles of BPDN DCT coefficients based on distribution of outputs from re-shuffled data 
14. BPDN_cwt2findcycles: Function to run CWT (on continuous data) and identify cycle peaks
15. BPDN_getPercentile: Function to compare output with reference distribution and compare percentiles
16. BPWPMethodsPaper_CodeWithDocumentation: Script capturing code and analysis conducted to generate each figure in the manuscript

Please contact the author for assistance if needed.
