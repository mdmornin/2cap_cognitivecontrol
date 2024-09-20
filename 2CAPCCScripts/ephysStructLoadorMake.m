function ephysStruct = ephysStructLoadorMake(dataPath,checkPath,savePath,tablePath)
    
% First determine if trlStruct already exists at checkPath
% If trlStruct already exists, load it and end function. If not, continue.
if exist(checkPath,'file')
    tic;
    load(checkPath)
    toc;
    disp(['File exists and has been loaded into the workspace after ' num2str(toc) 'seconds'])
    ephysStruct = ephysStruct;
    return 
end

% 
if exist(tablePath,'file')
    load(tablePath);
else
    warning('No data table found at passed location.')
end


% Generate dataset IDs
if exist(dataPath,'dir')
    dirList = dir(dataPath);
    fileNames = {dirList.name};
    fileNames = fileNames(~ismember(fileNames, {'.','..'}));
    dataSetIDs = fileNames;
else
    warning('No data found')
    return
end

ephysStruct = struct();

for i = 1:length(dataSetIDs)
    
    load([dataPath,filesep,dataSetIDs{i}])
    ephysStruct(i).Session = dataSetIDs{i};
    ephysStruct(i).stmtx = stMtx;
    ephysStruct(i).fr = fr;

    rawDataFolder = strsplit(dataSetIDs{i},"_");
    rawDataFolder = strsplit(rawDataFolder{1},'S');
    rawDataFolder = [rawDataFolder{1} '-S' rawDataFolder{2}];
    rawDataFolder = ['rawData' filesep rawDataFolder filesep 'openEphysData'];
    

    for k = 1:64
        if exist([savePath rawDataFolder filesep '100_CH' num2str(k) '.continuous'],'file')
            rawData = load_open_ephys_data_chunked([savePath rawDataFolder filesep '100_CH' num2str(k) '.continuous'],300,420);
            rawData = downsample(rawData,30);
            noise(k) = std(rawData);
        else
            disp(['File ' savePath rawDataFolder filesep '100_CH' num2str(k) '.continuous  not found.'])
            noise(k) = NaN;
        end
    end

    [~,minNoise] = min(noise);
    rawData = load_open_ephys_data_faster([savePath rawDataFolder filesep '100_CH' num2str(minNoise) '.continuous']);
    rawData = downsample(rawData,30); 
    ephysStruct(i).LFP = rawData;
    ephysStruct(i).LFP_noise = noise;

    progress = i/length(dataSetIDs) * 100;

    disp(['Processed ' num2str(progress) '% of datasets'])
end

    if exist(savePath,'dir')
        save([savePath filesep 'ephysStruct.mat'],'ephysStruct', '-v7.3')
    else
        mkdir(savePath)
        save([savePath filesep 'ephysStruct.mat'],'ephysStruct', '-v7.3')
    end

    return

end
