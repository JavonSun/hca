% example Matlab script of using HcaReml for performing heritable component
% analysis and deriving quantitative trait from a collection of sub-traits with genetic
% correlation matrix calculated by GCTA 

% before running this script, set the folder containing this script to
% be the working folder in Matlab, otherwise you many encounter "file not
% found" error

% Javon, 10/3/2015

clear
clc

addpath(sprintf('%s/..', pwd));

% read in the genetic correlation matrix (the output of GCTA)
grm = dlmread(sprintf('%s/data/gcta-caus.grm', pwd), '\t'); % in sparse format
grm = grm(:, [1 2 4]);
grm_id = dlmread(sprintf('%s/data/gcta-caus.grm.id', pwd), '\t'); % load id of individuals in the above GRM
% convert the GRM to dense format
grm = spconvert(grm);
grm = full(grm);
full_grm = grm + grm' - diag(diag(grm));

N = size(grm_id, 1);

% load the phenotypic data
phe = dataset('XLSFile', sprintf('%s/data/phenotype.csv', pwd));
d = size(phe, 2) - 2;
% note that here we have the order of individuals in phenotype file well aligned with that in GRM

% load covariates
cov = dataset('file', sprintf('%s/data/categorical-covariates.covar', pwd), ... % load categorical covariate: sex
    'Delimiter', ' ', 'VarNames', {'FID', 'IID', 'sex'}, 'ReadVarNames', false);
gender = zeros(size(cov, 1), 1);
gender(strcmp(cov.sex, 'M')) = 1;
qcov = dlmread(sprintf('%s/data/quantitative-covariates.qcovar', pwd), ' '); % load quantitative covariates
C = [gender qcov(:, 3)];

% perform heritable component analysis using HcaReml
lambda = 0; % here we arbitrarily choose lambda = 0, in practice, it can be determined via cross validation
[t, w, funval] = HcaReml(double(phe(:, 3:size(phe, 2))), C, full_grm, lambda);

