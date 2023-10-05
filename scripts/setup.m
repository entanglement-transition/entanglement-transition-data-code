%% define paths
restoredefaultpath

drtn = dir('./');
alpha_string = regexp(drtn(1).folder,'alpha(?<alpha>[0-9]+)','match');
epsilon_string = regexp(drtn(1).folder,'epsilon(?<epsilon>[0-9]+)','match');

string = drtn(1).folder;
splitString = strsplit(string, '_');
version_string  = splitString{end};

split_result = regexp(drtn(1).folder,'entanglement','split');
root_path = fullfile(split_result{1},'entanglement');
addpath(genpath(fullfile(root_path,'functions')));


%% useful global definitions
alpha_list = [38,66,76,100,200];
epsilon_list = [0,5,10,15];
diameter_list = [0.66,0.76,0.66,0.5,0.25]/0.079;
length_list = [25,50,50,50,50];
area_list = [9.5,14.5,14.5,14.5,14.5].^2;


%% plot attributes

set(0,'defaultaxesfontsize',12)
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Response dat1
clr = [166,206,227
    31,120,180
    178,223,138
    51,160,44
    251,154,153
    227,26,28
    253,191,111
    255,127,0
    202,178,214
    106,61,154]/255;

%% additional info

stack_size_list{1} = [[1251,1251,1161];
    [1251,1251,1121];
    [1251,1251,1056];
    [1251,1251,1031]];
% 66
stack_size_list{2} = [[2000,2000,1301];
    [2000,2000,1221];
    [2000,2000,1151];
    [2000,2000,1171]];

% 76
stack_size_list{3} = [[2000,2000,1051];
    [2000,2000,961];
    [2000,2000,951];
    [2000,2000,966]];
% 100
stack_size_list{4} = [[2000,2000,1191];
    [2000,2000,1111];
    [2000,2000,1101];
    [2000,2000,1101]];
% 200
stack_size_list{5} = [[2000,2000,1250];
    [2000,2000,1251];
    [2000,2000,1251];
    [2000,2000,1251]];
