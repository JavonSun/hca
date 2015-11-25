% example Matlab script of using HcaMl for performing heritable component analysis and 
% deriving high heritable quantitative trait from a collection of sub-traits with pedigree data

% before running this script, set the folder containing this script to
% be the working folder in Matlab, otherwise you many encounter "file not
% found" error

% Javon, 10/2/2015

clear
clc

addpath(sprintf('%s/..', pwd));

% load data
ped = dataset('XLSFile', sprintf('%s/data/pedigree.csv', pwd)); % load pedigree
phe = dataset('XLSFile', sprintf('%s/data/phenotype.csv', pwd)); % load phenotype data
cov = dataset('XLSFile', sprintf('%s/data/covariates.csv', pwd)); % load data on covariates

families = get_families(double(ped)); % extracts families from the pedigree file
kmat = kinship_matrix(double(ped), families); % create the kinship matrix

% remove subjects that are not in phenotype file from families and kmats
map = match(ped.id, phe.id);
[f2, m2] = remap_mat(map, families, kmat);
k2 = m2(:, 1);

% run the analysis
lambda = 1; % here we arbitrarily choose lambda = 1, in practice, it can be determined via cross validation

% analysis without correcting covariates effect
[qt, w] = HcaMl(double(phe(:, 2:end)), zeros(0,0), ... % the output qt is the returned high heritable trait
    double(ped(match(phe.id, ped.id), 4)), f2, k2, lambda);

% analysis with correcting covariates effect
[qt_c, w_c, r, v] = HcaMl(double(phe(:, 2:end)), double(cov(:, 2:end)), ... % the output qt_c is the returned high heritable trait
    double(ped(match(phe.id, ped.id), 4)), f2, k2, lambda);



