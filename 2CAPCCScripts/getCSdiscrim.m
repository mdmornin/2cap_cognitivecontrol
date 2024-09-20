function [csDiscrim,trlStruct] = getCSdiscrim(trlStruct,masterTbl,flags)

%% Params
LeftTrials = 1:24;
RightTrials = 25:48;
disThresh = 9;
disRun = 3;
sipDescent = 10*30;
sipAscent = 18*30;
cueOn = 5*30;
    % Apporach Likelihood. Correct or Incorrect Approaches (Not combined).
    % First gather data, then plot

RegcorApproachVector = []; RevcorApproachVector = [];
RegincApproachVector = []; RevincApproachVector = [];

if strcmp(flags.SessionN,'all') && strcmp(flags.Genotype,'all')    

    for i = 1:length(trlStruct)
        if startsWith(masterTbl.SessionType{i},['Regular'])

            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));

            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            csMinApproachVector = trlStruct(i).approach(49:end,4,1);

            sortedCorrectApproach_sansCorrections = double(currentApproachVector(sortedIndex,2,1) == 0 & currentApproachVector(sortedIndex,1,1) == 1);
            sortedCorrectApproach_withCorrections = double(currentApproachVector(sortedIndex,1,1) == 1);
            sortedIncrectApproach = double(currentApproachVector(sortedIndex,2,1) == 1);
            sortedAllApproaches = double(currentApproachVector(sortedIndex,1,1) == 1) | double(currentApproachVector(sortedIndex,2,1) == 1);

            trlStruct(i).CorrectApproach_NoCorrections = sum(sortedCorrectApproach_sansCorrections);
            trlStruct(i).CorrectApproach_Corrections = sum(sortedCorrectApproach_withCorrections);
            trlStruct(i).CSMinusApproaches = sum(csMinApproachVector == 1);
            trlStruct(i).CSMinAppLH = movmean(csMinApproachVector,3);
            trlStruct(i).CorrectApproachLH = movmean(sortedCorrectApproach_sansCorrections,3);
            trlStruct(i).IncrectApproachLH = movmean(sortedIncrectApproach,3);
            trlStruct(i).CSRatio = sum(csMinApproachVector == 1) / (sum(sortedAllApproaches == 1) + sum(csMinApproachVector == 1));
            trlStruct(i).CSRatio2 = sum(csMinApproachVector == 1) / sum(sortedCorrectApproach_sansCorrections);



        elseif startsWith(masterTbl.SessionType{i},['Reversal'])
            currentTrialTimes = trlStruct(i).trialTimes;
            [~,sortedIndex] = sort(currentTrialTimes(1:48));

            currentApproachVector = trlStruct(i).approach(1:48,:,1);
            csMinApproachVector = trlStruct(i).approach(49:end,4,1);

            sortedCorrectApproach_sansCorrections = double(currentApproachVector(sortedIndex,2,1) == 1 & currentApproachVector(sortedIndex,1,1) == 0);
            sortedCorrectApproach_withCorrections = double(currentApproachVector(sortedIndex,2,1) == 1);
            sortedIncrectApproach = double(currentApproachVector(sortedIndex,1,1) == 1);
            sortedAllApproaches = double(currentApproachVector(sortedIndex,1,1) == 1) | double(currentApproachVector(sortedIndex,2,1) == 1);


            trlStruct(i).CorrectApproach_NoCorrections = sum(sortedCorrectApproach_sansCorrections);
            trlStruct(i).CorrectApproach_Corrections = sum(sortedCorrectApproach_withCorrections);
            trlStruct(i).CSMinusApproaches = sum(csMinApproachVector == 1);
            trlStruct(i).CSMinAppLH = movmean(csMinApproachVector,3);
            trlStruct(i).CorrectApproachLH = movmean(sortedCorrectApproach_sansCorrections,3);
            trlStruct(i).IncrectApproachLH = movmean(sortedIncrectApproach,3);
            trlStruct(i).CSRatio = sum(csMinApproachVector == 1) / sum(sortedAllApproaches == 1);
            trlStruct(i).CSRatio2 = sum(csMinApproachVector == 1) / sum(sortedCorrectApproach_sansCorrections);


        end

    end

    csDiscrim = [];
end