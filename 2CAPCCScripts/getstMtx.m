clear all
% Requires spikes-master (found on GitHub)
% Requires npy-matlab-master (found on GitHub)
% Requires padcat.m (found on MathWorks / MATLAB's File Site)
% set paths for where to find your data
addpath(genpath('F:\dissDat\spikes-master'))
addpath(genpath('F:\dissDat\npy-matlab-master'))
addpath(genpath('F:\dissDat\Desktop\padcat.m'))
% If navigated into data folder on MATLAB, pwd will be fine. 
% If not, set the variable myKsDir as the path to your Phy outputted data.

ksresdirpath = 'G:\17A\res';
ksdir = dir(ksresdirpath);
ksfolders = [ksdir.isdir];
allKsFolders = ksdir(ksfolders);
allKsFoldersNames = {allKsFolders(3:end).name};

% allresFolders(:) = ksdir(resFolders == 1).name;

for k = 1:length(allKsFoldersNames)

    myKsDir = [ksresdirpath filesep allKsFoldersNames{k}];




% myEventTimes = load('C:\...\data\someEventTimes.mat'); % a vector of times in seconds of some event to align to

% Loading data from kilosort/phy easily

sp = loadKSdir(myKsDir)

%sp.st are spike times in seconds
%sp.clu are cluster identities
%spikes from clusters labeled "noise" have already been omitted

if ~isempty(sp.clu)
    clusters = double(sp.clu);
    spike_times = double(sp.st);
    p_stmtx = [clusters spike_times];
    i_stmtx = unique(p_stmtx(:,1));
    id_mtx = i_stmtx';
    ir_stmtx = i_stmtx+1;
    for i = 1:numel(id_mtx)
        A = find(clusters == id_mtx(i));
        idx_mtx{i} = A;
    end
    for i = 1:numel(idx_mtx)
        stMtx{i} = spike_times(idx_mtx{i});
    end
    c2 = cellfun(@(x) [x(:)], stMtx, 'un',0);
    [stMtx,tf] = padcat(c2{:});
    writeNPY('stMtx',[myKsDir filesep 'stMtx.npy'])	%Writes stMtx variable to .npy for Python purposes if desired
    %pBin sets the bin width
    %pBin 1 = 1 second, pBin 60 = 1 minute, etc.
    %h variable is the binned spikes. i.e. how many spikes occured per time point.
    %
    pBin = 1; h=histc(stMtx,[min(stMtx(:)):pBin:max(stMtx(:))]);
    fr = h./pBin;
    
    save([myKsDir filesep 'spikes.mat'],'fr','stMtx')
end
clear sp clusters spike_times p_stmtx i_stmtx id_mtx ir_stmtx A idx_mtx tf c2 stMtx
end
% Will need to save stMtx variable from here
% Also outputs a binned firing rate variable 