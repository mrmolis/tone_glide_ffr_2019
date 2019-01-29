%% SECTION 0: SETTINGS & OPTIONS
inPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr_withFlat.mat'];
subsetName = 'respHorzVert_minusOnly_noArt_ynOnly_flatSubsOnly_xcorr_snr.mat';

%% SECTION 1: LOAD AND PREPARE DATA
inDir = pickOffDir(inPath);
subsetPath = [inDir '\' subsetName];
if exist(subsetPath,'file')
    load(subsetPath);
else
    load(inPath);
    upresp; % unpack fields into variables
    
    isFlat = strcmp(slope,'00');
    flatSids = unique(sid(isFlat));
    slopeLogic = ismember(slope,{'00','23'});
    sidLogic = ismember(sid,flatSids);
    dirLogic = ~strcmp(direction,'down');
    groupLogic = strcmp(group,'YN');
    polLogic = strcmp(pol,'minus');
    
    resp = resp(slopeLogic & sidLogic & dirLogic & groupLogic & polLogic...
        & ~isArt);
    save(subsetPath,'resp');
end

upresp; % update derived variables

% Clear outdated vars
clear isFlat *Logic

%% SECTION 2: TESTING
slopeList = unique(slope);
montList = unique(montType);
for i=1:numel(slopeList)
    currSlope = slopeList{i};
    slopeLogic = strcmp(slope,currSlope);
    for j=1:numel(montList)
        currMont = montList{j};
        montLogic = strcmp(montType,currMont);
        dataSnr = snr(slopeLogic & montLogic,:);
        [hSnr(i,j),pSnr(i,j),ciSnr{i,j},statsSnr(i,j)] = ...
            ttest(dataSnr(:,1),dataSnr(:,2),'tail','right');
        dataRho = rho(slopeLogic & montLogic,:);
        [hRho(i,j),pRho(i,j),ciRho{i,j},statsRho(i,j)] = ...
            ttest(dataRho(:,1),dataRho(:,2),'tail','right');
    end; clear j
end; clear i


