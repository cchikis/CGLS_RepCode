%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: main_cgls.m
% Author: Craig A. Chikis
% Date: 01/06/2022
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -global
clear all
close all

% Fill this in per your local file setup
main_dir = "/home/cachikis/CGLS_RepCode/"; 


% You need NLopt for some of this code; enter the path where the MATLAB API is saved
% https://nlopt.readthedocs.io/en/latest/ Has instructions how to install (it is open-source)
nlopt_dir = "/home/cachikis/NLopt/mex/";
% If you don't have it just uncomment this line
% nlopt_dir = ""; 

% Set up file paths (should be self-contained from here on)
cd(main_dir)
addpath(nlopt_dir)
addpath('Code/Main/')
addpath('Code/Model/')
addpath('Code/Analysis/')
addpath('Code/Parameterizations/')

% Set up the parallel cluster 
try
	c = parcluster('local'); 
	pp = parpool(min(c.NumWorkers, 30));
catch
	warning('Parallel pool already running.')
end
st = tic;
run_paper_new(nlopt_dir); 
rt = toc(st);  

disp("Elapsed time is " + num2str(rt/60/60, '%0.1f') + " hours.") 