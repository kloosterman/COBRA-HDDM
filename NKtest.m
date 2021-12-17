clear all
cd /Volumes/FB-LIP/Projects/COBRA/data/behav_data/Niels_edits
addpath(genpath(pwd))

% filename = '/Volumes/FB-LIP/Projects/COBRA/data/behav_data/Eprime_behav_data/COBRA_n-back(2012-09-27)-003-1.txt'
filename = '/Volumes/FB-LIP/Projects/COBRA/data/behav_data/Eprime_behav_data/COBRA_n-back(2012-09-27)-214-1.txt'
TR = 2;  
experiment=[];
experiment.subject = '003'
experiment.session = 1;
[block, RunningList, LevelList, OnsetLists, OnsetList]  = dz_InterpreteEprimefile(filename, TR, 'OnsetTime', '', 'TextDisplay1');
type = 'COBRA';


NK_COBRA_nback_performance_Mar28_2017