# Tol-Pal-model
MATLAB scripts for a model of the Tol-Pal system in E. coli. 

Requirements: Matlab with Statistics (mean, median, bootci,..), Image Processing (immse), Global Optimization (patternsearch) and Parallel Computing toolboxes. Tested on MATLAB R2020b on Ubuntu.

Analytical folder contains scripts to produce figures for the analytical toy model 

Functions:

- post_bleach.m, function using pdepe to solve the toy model given an initial bleached condition
- pre_blach.m, function using pdepe to find the steady state for the toy model

Scripts:

- dimensionless.m, file plotting the solutions to the toy model as different parameters are varied
- simulated_FRAP.m, file simulating spatialFRAP for homogenous transport and no transport

.MAT files:

- steady_state.mat, stores the steady state results for the toy model

Fitting folder contains all scripts to solve and fit to the full model

Functions:

- bindata.m, smoothen the data in column array by place it in bins and taking the mean
- effective_diff_*.m, uses fitkymo to find the effective diffusion coefficient for either Pal (non-)dividing or Pal in tolA or tolB mutants
- fitkymo.m, uses "patternsearch" from MATLAB's optimisation toolbox to find the diffusion coefficient giving the best fit to the data using Fokker-Planck diffusion
- fitkymo_pal.m, uses "patternsearch" from MATLAB's optimisation toolbox to find the best fit for the full model to the experimental data using the Pal (non-)dividing kymogrpahs and effective diffusion constants
- immsre.m, the mean square relative error between two matrices (based on MATLAB's immse.m)
- spatialFRAP_*.m, uses MATLAB's pdepe function to solve the full model for bleached initial conditions given values of the variables replicating spatialFRAP for either Pal in (non-)dividing cells or tolA, tolB mutants and TolB in (non-)dividing cells
- steady_state_*.m, uses MATLAB's pdepe function to solve the full model and find the steady state  for either Pal in (non-)dividing cells or tolA, tolB mutants and TolB in (non-)dividing cells

Scripts:

- fit_pal.m, makes use of functions to fit the experimental data to the full model
- plot_all, takes the results from fitting and plots kymographs, concentration profiles, and effective diffusion coefficients

.MAT files:

-fit_parameter.mat, stores the results of fitting

Import folder contains all of the experimental data and files to plot it

Functions:

- bindata.m, smoothen the data in column array by place it in bins and taking the mean
- fit_all_cells, runs fitkymo.m on all cells in a dataset, also calculates the median Deff for each cell and the average scaled signal from all cells
- fitkymo.m, uses "patternsearch" from MATLAB's optimisation toolbox to find the diffusion coefficient giving the best fit to the data, for comparison, uses the standard Fickian model of diffusion as well as Fokker-Planck diffusion
- immsre.m, the mean square relative error between two matrices (based on MATLAB's immse.m)
- import_data.m, reads data from excel files into a cell
- plot_data.m, plot the results for a single dataset, also plots the average kymograph and the best fit to the average
- shadederror.m, plot error bars as shading
- spatialFRAP.m, uses MATLAB's PDE solver to return a simulated kymorgraph given the intial data and diffusion coefficient, which may vary spatially
- violinplot.m, a function by Bastian Bechtold (available under a BSD 3-clause license from https://github.com/bastibe/Violinplot-Matlab)

Scripts:

- conc_cells.m, finds the correlation between concentration and cell length
- import_other_data.m, imports data from excel files and saves it in a .mat file
- main.m, calls functions to import data from an excel file, process and plot this data and save the results in a .mat file, if .mat file already exists then plots existing data
- plot_some, plots the results of data from multiple datasets, plots the average kymograph, violin plot and effective diffusion coefficient
- TolA_dist.m, loads TolA overexpression .mat file and plots data

.MAT files:

- *.mat, MATLAB data files with fluorescence data, the pixel size and the binning factor to use, as well as the produced results

