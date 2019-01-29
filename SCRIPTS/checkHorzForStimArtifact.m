clearvars;

%% SECTION 0: SETTINGS AND OPTIONS
dataPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respYnOnly_horzOnly_minusOnly_xcorr_snr.mat'];

 
%% SECTION 1: LOAD AND PREPARE DATA
load(dataPath);
upresp; % unpack useful subfields into variables for easy access

% Find data of interest for comparison
artSlopes = unique(slope(isArt));
artSids = unique(sid(isArt));
artDirs = unique(direction(isArt));
artPols = unique(pol(isArt));
hasArtSlope = ismember(slope, artSlopes);
hasArtSid = ismember(sid, artSids);
hasArtDir = ismember(direction, artDirs);
hasArtPol = ismember(pol, artPols);
hasArt = hasArtSlope & hasArtSid & hasArtDir & hasArtPol;

% Pare down to data of interest
resp = resp(hasArt);
upresp; % update variables to correspond to pared-down resp struct

%% SECTION 2: TESTING
% T-test for each subject to see if significant artifact present

nullMeanSnr = 0; 
nullMeanPcc = 0;
for i=1:numel(artSids)
    currSid = artSids(i);
    dataSnr = reshape(snr(isArt,:),[],1);
    [hSnr(i),pSnr(i),ciSnr{i},statsSnr(i)] = ttest(dataSnr,nullMeanSnr,...
        'tail','right');
    dataPcc = reshape(rho(isArt,:),[],1);
    [hPcc(i),pPcc(i),ciPcc{i},statsPcc(i)] = ttest(dataPcc,nullMeanPcc,...
        'tail','right');
end; clear i

%% SECTION 3: PLOTTING


%% SECTION 4: SAVING


%% SECTION 5: CLEANUP