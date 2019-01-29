clearvars;

%% SECTION 0: SETTINGS AND OPTIONS
dataPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\respHorzVert_xcorr_snr.mat'];
 
%% SECTION 1: LOAD AND PREPARE DATA
dataDir = pickOffDir(dataPath);
smallPath = [dataDir '\respHorzVert_ynOnly_minusOnly.mat'];
if ~exist(smallPath,'file')
    load(dataPath);
    upresp; % unpack useful subfields into variables for easy access
    
    % Pare down to minus and YN only
    resp = resp(strcmp(group,'YN') & strcmp(pol,'minus'));
    upresp; % update derived vars
    
    % Save for quicker access later
    save(smallPath,'resp');
else
    load(smallPath);
    upresp;
end

% Find data of interest for comparison
artSlopes = unique(slope(isArt));
artSids = unique(sid(isArt));
artDirs = unique(direction(isArt));
artMonts = unique(montType(isArt));
hasArtSlope = ismember(slope, artSlopes);
hasArtSid = ismember(sid, artSids);
hasArtDir = ismember(direction, artDirs);
hasArtMont = ismember(montType, artMonts);
hasArt = hasArtSlope & hasArtSid & hasArtDir & hasArtMont;

% Pare down to data of interest
resp = resp(hasArt);
upresp; % update derived vars

%% SECTION 2: TESTING
% T-test for each subject to see if significant artifact present

nullMeanSnr = 0; 
nullMeanPcc = 0;
for i=1:numel(artMonts)
    currMont = artMonts{i};
    isCurrMont = strcmp(montType,currMont);
    dataSnr{i} = reshape(snr(isArt & isCurrMont,:),[],1);
    [hSnr(i),pSnr(i),ciSnr{i},statsSnr(i)] = ttest(dataSnr{i},nullMeanSnr,...
        'tail','right');
    dataPcc{i} = reshape(rho(isArt & isCurrMont,:),[],1);
    [hPcc(i),pPcc(i),ciPcc{i},statsPcc(i)] = ttest(dataPcc{i},nullMeanPcc,...
        'tail','right');
end; clear i

%% SECTION 3: PLOTTING


%% SECTION 4: SAVING


%% SECTION 5: CLEANUP