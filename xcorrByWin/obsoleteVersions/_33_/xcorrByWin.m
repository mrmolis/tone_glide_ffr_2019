function xcorrByWin(respPath,outPath)
ei.functionSpec = mfilename;
ei.versionSpec = '_33_';

%% SECTION 0: USER SETTINGS
disp('SECTION 0: USER SETTINGS');
defaultSaveDir = ['C:\Users\vhapormadseb\Desktop\tempOutput\' mfilename];
defaultSaveName = 'resp_xcorr.mat';
defaultRespDir = ['C:\Users\vhapormadseb\Desktop\Molis FFR datamats\'...
    'tone glides study\newProcess2017All'];
defaultFs = 20000;  %#ok<NASGU>
nWins = 3; 
overlapRatio = 0;
lagMsMin = 4;
lagMsMax = 12;
winType = 'blackman';

% Account for any missing inputs
if ~exist('respPath','var')
    msg = 'Select mat file with compiled FFR data';
    [respFile,respDir] = uigetfile(defaultRespDir,msg);
    respPath = [respDir '\' respFile];
    loadNow = true;
    % Abort function if user cancels
    if ~respFile
        msg = 'No FFR data selected, aborting function';
        warning(msg);
        return
    end
elseif isstruct(respPath)
    in = respPath;
    clear respPath
    loadNow = false;
else
    loadNow = true;
end


defaultSavepath = [defaultSaveDir '\' defaultSaveName];
if ~exist('outPath','var')
    [outName,outDir] = uiputfile(defaultSavepath,'Save output as...'); 
    if ~outName
        warning('Output will not be saved');
    else
        if ~exist(outDir,'dir')
            mkdir(outDir);
        end
        outPath = [outDir '\' outName];
    end
end

if loadNow
    disp('Loading response data...');
    in = load(respPath);
    disp('Loaded!');
end

%% SECTION 1: PROCESS XCORR
disp('SECTION 1: PROCESS XCORR');
resp = in.resp;

for i=1:numel(resp)
    %% Progress display
    prog = [num2str(i) '/' num2str(numel(resp))];
    disp(['XCorr: ' prog]);
    
    %% Getting started
    r = resp(i);
    fs = r.fs; %#ok<NASGU>
    spms = r.spms;
    pss = r.epochInfo.prestimSamps;
    ffr = shiftdim(r.wavAvg);
    stim = shiftdim(r.stim.wav);
    ffrTrimmed = ffr(pss+1:end);
    pad = zeros(length(ffrTrimmed)-length(stim),1);
    stimPadded = [stim;pad];
    
    %% Store or load some stuff
    if ~isfield(r,'winInfo') || ~isstruct(r.winInfo)
        resp(i).winInfo = struct;
        sortfields(resp);
        r = resp(i);
    end
    
    rwi = r.winInfo;
    if ~isfield(rwi,'n')
        resp(i).winInfo.n = nWins;
    else
        nWins = rwi.n;
    end
    if ~isfield(rwi,'overlapRatio')
        resp(i).winInfo.overlapRatio = overlapRatio;
    else
        overlapRatio = rwi.overlapRatio;
    end
    if ~isfield(rwi,'type')
        resp(i).winInfo.type = winType;
    else
        winType = rwi.type;
    end
    if ~isfield(rwi,'durMs')
        stimMs = r.stim.durationMs;
        wMs = floor(stimMs/(nWins+overlapRatio-nWins*overlapRatio));
        resp(i).winInfo.durMs = wMs;
    else
        wMs = rwi.durMs;
    end
    if ~isfield(rwi,'durSamps')
        wSamps = wMs * spms;
        resp(i).winInfo.durSamps = wSamps;
    else
        wSamps = rwi.durSamps;
    end
    if ~isfield(rwi,'overlapSamps')
        overlapSamps = floor(wSamps * overlapRatio);
        resp(i).winInfo.overlapSamps = overlapSamps;
    else
        overlapSamps = rwi.overlapSamps;
    end
    if ~isfield(rwi,'overlapMs')
        overlapMs = overlapSamps * spms;
        resp(i).winInfo.overlapMs = overlapMs;
    else
        overlapMs = rwi.overlapMs; %#ok<NASGU>
    end

    %% 1.1) GET LAG
    try % don't need all these vars here, but want to make sure they exist
        lagMs = r.lag.ms; %#ok<NASGU>
        lagSamps = r.lag.samps;
        lagMsMin = r.lag.searchMsMin;
        lagMsMax = r.lag.searchMsMax;
        lagSampsMin = r.lag.searchSampsMin; %#ok<NASGU>
        lagSampsMax = r.lag.searchSampsMax; %#ok<NASGU>
    catch
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
    end
    
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
    
    %% GET WINDOWED CORRELATION DATA
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

%% 3.3) Export
if isfield(in,'exportInfo')
    exportInfo = in.exportInfo;
    ind = numel(exportInfo) + 1;
    exportInfo(ind).functionSpec = ei.functionSpec;
    exportInfo(ind).versionSpec = ei.versionSpec;
else
    exportInfo = ei; %#ok<NASGU>
end
if exist('outPath','var')
    save(outPath,'resp','exportInfo');
end

disp([mfilename ' completed!']);
% Double beep when done (doesn't currently work)
dos(['�&' 'exit&']);
dos(['�&' 'exit&']);

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

% don't let it be at range endpoints, because it's probably not true peak
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
%--------------------------------------------------------------------------
%% mergestructs
% Note: if any fields are shared between the input structs, the values in
% the second will overwrite the values in the first.
function out = mergestructs(s1,s2)

s3 = struct;
s = {s1,s2};
for K=1:2
    f = fieldnames(s{K});
    for I=1:numel(s{K})
        for J=1:numel(f)
            fld = f{J};
            s3(I).(fld) = s{K}(I).(fld);
        end; clear J
    end; clear I
end; clear K

out = s3;

end % end mergestructs
%--------------------------------END HELPER FUNCTIONS----------------------

% [EOF]