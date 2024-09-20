%% Master Script for the Analysis of 2CAP Regular + Reversal Data
% Utilizes tables generated from experimental Excel sheets
% Generates MATLAB structure with structured, contextualized data
% Runs series of functions or scripts that generate outputs relevant to
% analyses at hand.
% 
% Required Analyses Script Names:
% {LIST}
% List of Analyses Completed:
% {LIST}
% 
% MDM 01/04/2022
% 
% Suggested improvements for future go below
% {Improvements}
% Bug fixes go below
% {Bugs}
%% Create core data structure variable housing relevant data

% Relevant data includes:
    % dataset name
    % approach 3D matrix
    % rat body coordinates
    % trial times
    % trial types
    % apparats coordinates
% Future data includes:
    % Spike timing matrices
    % PSTHs
    % Spike firing rates
    % Downsampled LFPs

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This section creates a structure variable utilizing behavioral tracking
% data from deeplabcut. This data is saved in .mat files from previous
% scripts and dependecies pWistDataOrganizerVer3.m and loadMAEvents.m
% loadMAEvents.m requires OpenEphys' analysis-tools from github.
% Error correction in tracking data was done previously in
% pWistDataOrganizerVer3.m. Additional requirement is the file
% masterTable.mat which houses drinking information, animal strain, session
% types, and animal weights. It is created with generateTable.m with some
% manual effort. 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Section Last Modified by MDM on 01/04/2022 
% IUPUI - Lapish Lab

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Perform basic check of dataStruct existence, create first instatiation of
% it with only folder names.

% Set paths
fpath = 'F:\dissDat';
dlcpath = [fpath filesep 'res'];

% Check and load or create
if ~exist([fpath filesep 'dataStruct.mat'],'file')
    % Create file, label data set IDs, import deeplabcut data and data
    % processed from deeplabcut from pWistDataOrganizerVer3.m
    dataPath = 'F:/dissDat/rawData';    dirList = dir([dataPath]);  folderNames = {dirList([dirList.isdir]).name}; folderNames = folderNames(~ismember(folderNames ,{'.','..'}));
    dataDirs = folderNames;
    for i = 1:length(dataDirs)
        dataStruct(i).dataSetName = folderNames(k);
        strIdx = strsplit(dataStruct(i).dataSetName{1},'-');    % *** Less than ideal ***
        load([dlcpath filesep strIdx{1} '-' strIdx{2} '-' strIdx{3} strIdx{4} '_AllData.mat'])
        dataStruct(i).approach = approach;  % 3D Approach Matrix {Trial Num x Behavior Type
        %                                   {Correct Approach, Incorrect Approach, Correct Cue, Incorrect Cue} x Binary {0,1} and Duration of Behavior 
        dataStruct(i).bodyCoords = bodyCoords;  % Rat body coordinates {frameT,SnoutTip,HeadCap,MidBack,BackTail}
        dataStruct(i).trialTimes = trialTimes;  % Trial times corrected
        dataStruct(i).trialTypes = trialTypes;  % Trial types {e.g. CS+R CS+L CS -} 
        dataStruct(i).boxCoords = boxCoords;    % Box coordinates {left sipper, right sipper, upper left corner, lower left corner, upper right corner, lower right corner}
    end
    save([fpath filesep 'dataStruct.mat'],"dataStruct")
    disp(['dataStruct.mat was created at ' fpath])
else
    % Load data set
    disp(['Data structure exists in ' fpath])
    disp('Loading data structure')
    load([fpath filesep 'dataStruct.mat'])
end

% Perform check to load masterTable.mat
if exist([fpath filesep 'masterTable.mat'],'file')
    load([fpath filesep 'masterTable.mat'])
else
    disp('Missing masterTable.mat, please revisit generateTable.m')
end

%% From here onward we can utilize separate scripts and functions to keep the bulk of this main script organized -- to an extent.
% Find and plot basic metrics on correct trials, incorrect trials, and omissions
% utilizing the following script:
