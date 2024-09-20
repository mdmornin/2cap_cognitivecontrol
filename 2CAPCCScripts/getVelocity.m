function [trlStruct,velocityInfo] = getVelocity(trlStruct,figSavePath,masterTbl,flags)
% Params
trialN = 96; % First 48 trials are CS+
if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    
    for i = 1:length(trlStruct)

        if startsWith(masterTbl.SessionType{i},['Regular']) 
            for k = 1:trialN
                trlStruct(i).Velocity(:,k) = diag(pdist2(squeeze(trlStruct(i).trlnosePos(k,:,:)), squeeze(trlStruct(i).trlnosePos(k,:,:))),1)/30;
            end
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) 
            for k = 1:trialN
                trlStruct(i).Velocity(:,k) = diag(pdist2(squeeze(trlStruct(i).trlnosePos(k,:,:)), squeeze(trlStruct(i).trlnosePos(k,:,:))),1)/30;
            end
        end

    end
elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'P') 
    for i = 1:length(trlStruct)

        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'P')
            for k = 1:trialN
                trlStruct(i).Velocity(:,k) = diag(pdist2(squeeze(trlStruct(i).trlnosePos(k,:,:)), squeeze(trlStruct(i).trlnosePos(k,:,:))),1)/30;
            end
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'P')
            for k = 1:trialN
                trlStruct(i).Velocity(:,k) = diag(pdist2(squeeze(trlStruct(i).trlnosePos(k,:,:)), squeeze(trlStruct(i).trlnosePos(k,:,:))),1)/30;
            end
        end

    end
elseif strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'W') 
    for i = 1:length(trlStruct)

        if startsWith(masterTbl.SessionType{i},['Regular']) && strcmp(masterTbl.Strain{i},'W')
            for k = 1:trialN
                trlStruct(i).Velocity(:,k) = diag(pdist2(squeeze(trlStruct(i).trlnosePos(k,:,:)), squeeze(trlStruct(i).trlnosePos(k,:,:))),1)/30;
            end
        elseif startsWith(masterTbl.SessionType{i},['Reversal']) && strcmp(masterTbl.Strain{i},'W')
            for k = 1:trialN
                trlStruct(i).Velocity(:,k) = diag(pdist2(squeeze(trlStruct(i).trlnosePos(k,:,:)), squeeze(trlStruct(i).trlnosePos(k,:,:))),1)/30;
            end
        end

    end
end

velocityInfo = []; % Blank for now
end