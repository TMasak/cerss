calls_clean_all and calls_clean_all_more_columns contains raw data

qualitative_script.R performs the qualitative real data analysis, namely
- rounds up data on a common grid
- estimates the second order structure
- performs prediction (with confidence bands) for a selected surface
  -- implementation of prediction is currently vanilla, using full covariance tensor structure
Data needed for plotting purposes are saved by qualitative_script.R into TRout.RData

k_fold_serial_...R is a script run on a cluster, calculating the results for quantitative comparison in real data analysis
grid20.Rdata are the results of k_fold_serial_...

quantitative_script.R adds the results of the quantitative analysis from grid_20.Rdata into TRout.RData,
so all data needed for plotting are in TRout.RData

results.R can be used to re-produce all the plots in the paper, corresponding to the data analysis

