%% main.m
%% This script produces all output associated with the FRBNY DSGE model. 

%% Initialization
clear
close all
set_paths % adds necessary paths
spec_557  % sets important variables and flags

%% Estimation
% In this stage we draw from the distribution for the parameters. The modal
% parameters as well as the draws of parameters, are outputted in the /save
% folder. This program (using defaults) will take approximately 2 days to run. 
gibb_est_ant 


%% Forecasting
% Here we produce forecasts for our observable variables, one associated
% with each draw of the parameters, and saves them in the /save folder. 
% This program (in parallel, using defaults) will approximately take 45 min
% time to run.
forecast_parallel_est_ant

%% Plotting
forplot          % produces series to be plotted
plotPresentation % produces plots of series outputted from forplot, which 
                 % are saved in the /graphs folder
