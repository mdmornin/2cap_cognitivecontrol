   %%
    clear all
    addpath(genpath('analysis-tools-master'))
    spikeDir = 'F:/dissDat/res';
    loadDir = 'F:/dissDat/curated';
    rawDir = 'F:/dissDat/rawData';
    dirList = dir(rawDir);
    folderNames = {dirList([dirList.isdir]).name};
    dataDirs = folderNames(~ismember(folderNames ,{'.','..'}));

    for i = 1:length(dataDirs)
        splitName = strsplit(dataDirs{i},'-');
        spikeDirNames{i} = [splitName{1} '-' splitName{2} '-' splitName{3} splitName{4}];
    end

%%
    for iFile = 41:length(dataDirs)
        % Load timestamps
        [~,timestamps,~] = load_open_ephys_data_faster([dataDirs{iFile},filesep,'openEphysData/100_CH1.continuous']);
        % Load spike data
        load([loadDir filesep spikeDirNames{iFile} '_AllData.mat'])
        if max(abs(diff(timestamps) - diff(timestamps(1:2)))) > 10^(-8)
            disp('Unstable time stamps.')
        end

        for iNeuron = 1:size(stMtx,2)
            neuronVector = stMtx(:,iNeuron)*30000;
            neuronVector = neuronVector(~isnan(neuronVector));
            stMtx2{iNeuron} = timestamps(round(neuronVector));
            
        end
        stMtx = padcat(stMtx2{:});
        stMtx2 = [];
        save([spikeDir filesep spikeDirNames{iFile},'_AllData.mat'],'sessionInfo','trialTimes','trialTypes','approach','approachParams','boxCoords','bodyCoords','stMtx','fr')
    end