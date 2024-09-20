%% 05/25/2022
% Because I am missing the bootISI function on the new Linux, I may focus
% on the LFP Gamma Power to confirm that the approach timepoint analysis is
% broadly correct. This is a sanity check for future analyses using the
% approach timepoint. There is also some cleaning I can do of the script
% and perhaps creating functions to start cleaning up this giant script. 
%% Load and check data
clear all; close all % Refresh workspace
% PATHS
% For brain3 Linux Machine:
if ispc
    parentPath = 'F:/dissDat';
    figPath = 'F:/dissDat/figs_sepPW';
else
    parentPath = '/research/dissDat';
    figPath = '/research/dissDat/figs';
end
% LOAD
load([parentPath filesep 'ephysStruct.mat']);
load([parentPath filesep 'trlStruct.mat']);
load([parentPath filesep 'masterTable.mat']);

addpath(genpath([parentPath filesep 'analysisScripts']))
%% Set indexing based on Congruent versus Incongruent Sessions
regPIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'P');
revPIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'P');
    
regWIdx = startsWith(masterTbl.SessionType,'Regular') & startsWith(masterTbl.Strain,'W');
revWIdx = startsWith(masterTbl.SessionType,'Reversal') & startsWith(masterTbl.Strain,'W');

%% Timepoints and other parameters
numTrials = 48;


sipDescent = 6;
sipAscent = 14;
cueOn = 2;

LeftTrials = 1:24;
RightTrials = 25:48;
colors = [{'#EDB120'}, {'#7E2F8E'}, {'#0072BD'},  {'#D95319'} ];

%%
% In order to create a PSTH with approach, it will be necessary to create a
% new variable akin to 'trialTimes'.
% This new variable will be trialTimes + time of approach.
% Time of approach is in the 'approach' variable. 
% First, there needs to be a determination of correct/incorrect choices,
% and then based on that determination add the time each correct/incorrect
% choice occurs to the time of 'trialTimes'

%% Approach PSTH Creation for LFP Data
% Generally this should be similar to the spike data without the complexity
% of multiple neurons. Therefore it should simply be a TrialxSxT (Trial by Signal by Time)
% matrix. Signal is 1D, therefore we can collapse on it making a 2D matrix.
% 
% The trickiest part will be getting the timestamps to match.
%
% Pull ephys and trial data, create PSTH
    pBin = 0.1; % 33 ms bins for fr variable, match video framerate
    PSTH_Approach = [];
    if ~isfield(ephysStruct,'PSTH_Approach')
        
        for i = 1:length(trlStruct)
    
            % Load trial times from trlStruct
            trialTimes = trlStruct(i).trialTimes(1:48);
            [trialTimes, trialTimesIdx] = sort(trialTimes);
            approachVar = min(trlStruct(i).approach(:,1:2,2),[],2);
            sortedApproachVar = approachVar(trialTimesIdx);

            if startsWith(masterTbl.SessionType(i),'Regular')            
                correctIndex = trlStruct(i).approach(1:48,1,1) == 1 & trlStruct(i).approach(1:48,2,1) == 0;
                incorrectIndex = trlStruct(i).approach(1:48,2,1) == 1;
            elseif startsWith(masterTbl.SessionType(i),'Reversal')
                correctIndex = trlStruct(i).approach(1:48,1,1) == 0 & trlStruct(i).approach(1:48,2,1) == 1;
                incorrectIndex = trlStruct(i).approach(1:48,1,1) == 1;
            end

            correctIndex = correctIndex(trialTimesIdx);
            incorrectIndex = incorrectIndex(trialTimesIdx);

            % Load ephys related data
            LFP = ephysStruct(i).LFP;
            % LFP data has been previously downsampled to 1000Hz; allows
            % for analysis of data @ 500 Hz and below.

            % Find time, in seconds, that each fr index corresponds to.
            frTime = (1:length(LFP))/100;   % Denominator corresponds with pBin
    
            dlcTimestamps = trlStruct(i).bodyCoords{1};
            dlcTimestamps = dlcTimestamps(dlcTimestamps >= 0);
            dlcTimestamps = round(dlcTimestamps,2);
            trialTimes = trialTimes + sortedApproachVar;
            trialTimes = round(trialTimes,2);
   %     
        for k = 1:length(trialTimes)
            [~,tMin] = min(abs((trialTimes(k) - 2) - dlcTimestamps));
            [~,tMax] = min(abs((trialTimes(k) + 2) - dlcTimestamps)); 
            [~,tMin_fr] = min(abs((dlcTimestamps(tMin)) - frTime));
            [~,tMax_fr] = min(abs((dlcTimestamps(tMax)) - frTime));            
            trlSize(k) = size(LFP(tMin_fr:tMax_fr,:),1);

            if size(LFP(tMin_fr:tMax_fr,:),1) == 401
                LFP_PSTH(k,:,:) = LFP(tMin_fr:tMax_fr,:);
            elseif size(LFP(tMin_fr:tMax_fr,:),1) == 402
                LFP_PSTH(k,:,:) = LFP(tMin_fr:tMax_fr-1,:);
            elseif size(LFP(tMin_fr:tMax_fr,:),1) == 400
                LFP_PSTH(k,:,:) = LFP(tMin_fr:tMax_fr+1,:);          
            end

        end

        
        ephysStruct(i).LFP_PSTH = LFP_PSTH;
        ephysStruct(i).LFP_Correct = LFP_PSTH(correctIndex,:);
        ephysStruct(i).LFP_Incorrect = LFP_PSTH(incorrectIndex,:);
        
        ephysStruct(i).LFP_CorrectM = mean(ephysStruct(i).LFP_Correct);
        ephysStruct(i).LFP_IncorrectM = mean(ephysStruct(i).LFP_Incorrect);
        
        
        [~,w,t,f,ps] = spectrogram(mean(LFP_PSTH));

        
        ephysStruct(i).Correct_SpectrogramRes = [{w} {t} {f} {ps}];
        
        [~,w,t,f,ps] = spectrogram(ephysStruct(i).LFP_IncorrectM,10,0,[]);

        ephysStruct(i).Incorrect_SpectrogramRes = [{w} {t} {f} {ps}];
        
        w = []; t = []; f = []; ps = [];
        for j = 1:48
            [~,w{j},t{j},f{j},ps{j}] = spectrogram(ephysStruct(i).LFP_PSTH(j,:),10,0,[6.5:0.5:9],1000,'power');
            ephysStruct(i).allSpectrogramRes = [{w} {t} {f} {ps}];
        end

    end

    end

%% Given LFPs, create 2x2x2 analyses 
% Genotype x Session Type x Correct/Incorrect
% Take incorrect/correct trial mean per session, generate spectogram data per session,
% average spectrogram data, plot!

%% Plot theta over trials for P's W's
theta = 3;     % Manually determined range 

for j = 1:length(trlStruct)
    for i = 1:length(trialTimes)
        thetaMeans(j,i) = mean(mean(ephysStruct(j).allSpectrogramRes{1,4}{i}(:,1:40)));
    end
end

% Index out

thetaWCon = thetaMeans(regWIdx,1:15);
thetaWInc = thetaMeans(revWIdx,1:15);

thetaPCon = thetaMeans(regPIdx,1:15);
thetaPInc = thetaMeans(revPIdx,1:15);

% Plot
figure
plot(movmean(mean(thetaWCon),3),'-o','LineWidth',3)
hold on
plot(movmean(mean(thetaWInc),3),'--o','LineWidth',3)
plot(movmean(mean(thetaPCon),3),'-o','LineWidth',3)
plot(movmean(mean(thetaPInc),3),'--o','LineWidth',3)
xlabel('Trials')
ylabel('Mean Theta Power')
legend({'Congruent Wistar','Incongruent Wistar','Congruent P','Incongruent P'})
set(gca,'TickDir','out','Box','off','FontName','Arial','FontSize',20,'FontWeight','bold','LineWidth',4)
