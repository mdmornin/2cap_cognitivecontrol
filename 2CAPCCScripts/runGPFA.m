%% Format Data from Ephys Struct to Make Compatibile with GFAP

% Load ephysStruct statement
load('F:/dissDat/dataStructures/ephysStruct_1ms.mat')
% Params
sipDescent = 6*1000;
sipAscent = 14*1000;
cueOn = 1*1000;
%% Loop
for i = 1:length(ephysStruct)
    dat = [];
    for k = 1:size(ephysStruct(i).PSTH,1)

        dat(k).spikes = logical(squeeze(ephysStruct(i).PSTH(k,5001:7000,:)))';

        dat(k).trialId = k;
    end

runIdx = i;
method = 'gpfa';

% Select number of latent dimensions
xDim = 8;
% NOTE: The optimal dimensionality should be found using 
%       cross-validation (Section 2) below.

% If using a two-stage method ('fa', 'ppca', or 'pca'), select
% standard deviation (in msec) of Gaussian smoothing kernel.
kernSD = 100;
% NOTE: The optimal kernel width should be found using 
%       cross-validation (Section 2) below.

% Extract neural trajectories
try
result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD);
catch
end
% NOTE: This function does most of the heavy lifting.

% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
% NOTE: The importance of orthnormalization is described on 
 %     pp.621-622 of Yu et al., J Neurophysiol, 2009.


end
