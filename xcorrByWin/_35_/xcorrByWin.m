function outPath = xcorrByWin(inPath,outPath,exportInfo)

ei.functionSpec = mfilename;
ei.versionSpec = '_35_';

%% SECTION 0: USER SETTINGS
disp('SECTION 0: USER SETTINGS');
defOutDir = ['C:\Users\vhapormadseb\Desktop\tempOutput\' mfilename];
defOutName = 'resp_xcorr.mat';
defInDir = ['C:\Users\vhapormadseb\Desktop\Molis FFR datamats\'...
    'tone glides study\newProcess2017All'];
defaultFs = 20000;  %#ok<NASGU>
nWins = 3;
overlapRatio = 0;
lagMsMin = 2;
lagMsMax = 22;
winType = 'blackman';

% Get path to input data
if ~exist('inPath','var')
    msg = 'Select mat file with compiled FFR data';
    [inFile,inDir] = uigetfile(defInDir,msg);
    inPath = [inDir '\' inFile];
    loadNow = true;
    % Abort function if user cancels
    if ~inFile
        msg = 'No FFR data selected, aborting function';
        warning(msg);
        return
    end
elseif isstruct(inPath)
    resp = inPath;
    clear respPath
    loadNow = false;
else
    loadNow = true;
end

% Get save path
defOutPath = [defOutDir '\' defOutName];
if ~exist('outPath','var')
    
    [outName,outDir] = uiputfile(defOutPath,'Save output as...');
    
    if ~outName
        outPath = defOutPath;
    else
        outPath = [outDir '\' outName];
    end
    
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
end

%% SECTION 1: SETUP
disp('SECTION 1: SETUP')
% Load input data
if loadNow
    disp('Loading FFR data...');
    in = load(inPath);
    resp = in.resp;
    exportInfo = in.exportInfo;
    disp('Loaded!');
end


%% SECTION 2: PROCESS XCORR
disp('SECTION 2: PROCESS XCORR');
for i=1:numel(resp)
    % Progress display
    prog = [num2str(i) '/' num2str(numel(resp))];
    disp(['XCorr: ' prog]);
    
    % Getting started
    r = resp(i);
    fs = r.fs; %#ok<NASGU>
    spms = r.spms;
    pss = r.epochInfo.prestimSamps;
    ffr = shiftdim(r.wavAvg);
    stim = shiftdim(r.stim.wav);
    ffrTrimmed = ffr(pss+1:end);
    pad = zeros(length(ffrTrimmed)-length(stim),1);
    stimPadded = [stim;pad];
    
    % Store or load some stuff
    if ~isfield(r,'winInfo') || ~isstruct(r.winInfo)
        resp(i).winInfo = struct;
        sortfields(resp);
        r = resp(i);
    end
    
    rwi = r.winInfo;
    resp(i).winInfo.n = nWins;
    resp(i).winInfo.overlapRatio = overlapRatio;
    resp(i).winInfo.type = winType;
    stimMs = r.stim.durationMs;
    wMs = floor(stimMs/(nWins+overlapRatio-nWins*overlapRatio));
    resp(i).winInfo.durMs = wMs;
    wSamps = wMs * spms;
    resp(i).winInfo.durSamps = wSamps;
    overlapSamps = floor(wSamps * overlapRatio);
    resp(i).winInfo.overlapSamps = overlapSamps;
    overlapMs = overlapSamps * spms;
    resp(i).winInfo.overlapMs = overlapMs;
    
    
    %% 1.1) GET LAG
    lagSampsMin = lagMsMin * spms;
    lagSampsMax = lagMsMax * spms;
    inds = lagSampsMin:lagSampsMax;
    
    lagSamps = getlag(ffrTrimmed,stimPadded,inds);
    lagMs = lagSamps / spms;
    
    resp(i).lag.searchMinSamps = lagSampsMin;
    resp(i).lag.searchMaxSamps = lagSampsMax;
    resp(i).lag.searchMinMs = lagMsMin;
    resp(i).lag.searchMaxMs = lagMsMax;
    resp(i).lag.samps = lagSamps;
    resp(i).lag.ms = lagMs;
    resp(i).lag = orderfields(resp(i).lag);
    resp = orderfields(resp);
    
    %% 1.2) GET OVERALL CORRELATION DATA (NON-WINDOWED)
    % Trim FFR all the way down to stim length, accounting for lag
    stimSamps = length(stim);
    inds = (1:stimSamps) + lagSamps;
    ffr = ffrTrimmed(inds);
    
    % Normalize scale to stim
    mult = rms(stim)/rms(ffr);
    ffrScaled = ffr .* mult;
    
    % Raw correlation
    [corr,lag] = xcorr(ffrScaled,stim);
    rawCorr = abs(corr(lag==0));
    
    % Correlation coefficient (rho)
    [autocorr,lag] = xcorr(stim);
    [ref,idx] = max(abs(autocorr));
    assert(lag(idx)==0,'autocorr peak not at zero lag'); % sanity check
    rho = rawCorr / ref;
    
    % Store
    resp(i).rawCorr = rawCorr;
    resp(i).rho = rho;
    
    %% 1.3) GET WINDOWED CORRELATION DATA
    for j=1:nWins
        offset = (j-1)*(wSamps - overlapSamps);
        inds = (1:wSamps) + offset;
        ffrWin = ffr(inds);
        stimWin = stim(inds);
        
        % Normalize scale to stim
        mult = rms(stimWin)/rms(ffrWin);
        ffrWinScaled = ffrWin * mult;
        
        % Raw correlation
        [corr,lag] = xcorr(ffrWinScaled,stimWin);
        rawCorr = abs(corr(lag==0));
        
        % Correlation coefficient (rho)
        [autocorr,lag] = xcorr(stimWin);
        [ref,idx] = max(abs(autocorr));
        assert(lag(idx)==0,'autocorr peak not at zero lag'); % sanity check
        rho = rawCorr / ref;
        resp(i).byWindow(j).rho = rho;
        resp(i).byWindow(j).rawCorr = rawCorr;
    end
end
resp = sortfields(resp); %#ok<NASGU>

%% SECTION 3: EXPORT
disp('SECTION 3: EXPORT')

if exist('exportInfo','var')
    ind = numel(exportInfo) + 1;
    exportInfo(ind).functionSpec = ei.functionSpec;
    exportInfo(ind).versionSpec = ei.versionSpec;
    exportInfo(ind).timestamp = datetime;
else
    exportInfo = ei;
    exportInfo.timestamp = datetime;
end

save(outPath,'resp','exportInfo');

disp([mfilename ' completed!']);
% Double beep when done (doesn't currently work)
dos(['•&' 'exit&']);
dos(['•&' 'exit&']);

end % end main function block (xcorrByWin and nested functions)

%% HELPER FUNCTIONS
%----------------------------BEGIN HELPER FUNCTIONS------------------------
% These functions only have read-access to variables that are explicitly
% passed as input arguments, and they only have write-access to variables
% that are explicitly passed as output arguments.
%--------------------------------------------------------------------------
%% get_olag
function lag = getlag(x1,x2,inds)
% Get two-sided vector of lags (-nSamps:nSamps) and raw correlations.
[corrVec,lagVec] = xcorr(x1,x2,'none');

% Make all correlations absolute
corrVec = abs(corrVec);

% Only allow positive lag
posLogic = (lagVec>0);
corrVecPos = corrVec(posLogic);

% Further restrict to lag range allowed in settings
corrVecAllowed = corrVecPos(inds);

% Get overall lag
lagSampsMin = inds(1);
[~,ind] = max(corrVecAllowed);

% Don't let it be at range endpoints, because it's probably not true peak
offset = 0;
while true
    if ind == 1
        corrVecAllowed = corrVecAllowed(2:end);
        [~,ind] = max(corrVecAllowed);
        offset = offset + 1;
    elseif ind == numel(corrVecAllowed)
        corrVecAllowed = corrVecAllowed(1:end-1);
        [~,ind] = max(corrVecAllowed);
    else
        ind = ind + offset;
        break
    end
end
lag = ind + lagSampsMin - 1;

end % end get_olag function
%---------------------END HELPER FUNCTIONS---------------------------------
% [EOF]