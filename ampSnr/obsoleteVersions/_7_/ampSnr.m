function ampSnr(varargin)
% ampSnr(respDir,stimDir,saveDir)
ei.functionSpec = mfilename;
ei.versionSpec = '_7_';

%% SECTION 0: USER SETTINGS
%% 0.1) Folder Settings
defaultSaveDir = ['C:\Users\vhapormadseb\Desktop\tempOutput\' mfilename];
defaultInDir = ['\\r01porhsm03\Research\NCRAR\Spectral_Dynamics_Grant\'...
    'FFR\Tone Glide Study\Methods Manuscript\xcorr_newProcess2017All'];
defaultSaveName = 'resp_snr.mat';

%% 0.2) Analysis Settings
fInterest = [350,750]; % set frequency range of interest
fStep = 1; % space between frequencies in Hz
nWins = 3; % number of desired windows
winType = 'blackman';
overlapRatio = 0;
peakSearchToleranceHz = 25;
lagMsMin = 4;
lagMsMax = 12;

%% SECTION 1: INITIALIZATION
% Parse varargin
if abs(nargin) < 2
    saveDir = '';
else
    saveDir = varargin{2};
end
if abs(nargin) < 1
    inPath = '';
else
    inPath = varargin{1};
end

% Account for any missing inputs
if isempty(inPath)
    msg = 'Select mat file with compiled FFR data';
    [inFile,inDir] = uigetfile(defaultInDir,msg);
    inPath = [inDir '\' inFile];
    loadNow = true;
    % Abort function if user cancels
    if ~inFile
        msg = 'No FFR data selected, aborting function';
        warning(msg);
        return
    end
    clear inFile inDir
elseif isstruct(inPath)
    in = inPath;
    loadNow = false;
else
    loadNow = true;
end

defaultSavePath = [defaultSaveDir '\' defaultSaveName];
if isempty(saveDir)
    if ~exist(defaultSaveDir,'dir')
        mkdir(defaultSaveDir);
    end
    [outName,outDir] = uiputfile(defaultSavePath,'Save output as...');
end
if ~outName
    warning('Output will not be saved');
else
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    outPath = [outDir '\' outName];
end

if loadNow
    disp('Loading response data...')
    load(inPath)
    disp('Loaded!');
end

%% SECTION 2: STIM DFT WINDOWED
disp('SECTION 2: STIM DFT');
if isfield(in,'resp')
    resp = in.resp;
else
    resp = in;
end
for i=1:numel(resp)
    prog = [num2str(i) '/' num2str(numel(resp))];
    disp(['Windowed stim DFT: ' prog]);
    
    % Some shorthand
    r = resp(i);
    stim = r.stim;
    durMs = stim.durationMs;
    spms = r.spms;
    
    % Store or load some stuff
    if ~isfield(r,'winInfo') || ~isstruct(r.winInfo)
        resp(i).winInfo = struct;
        resp = sortfields(resp);
        r = resp(i);
%         resp(i).winInfo.n = nWins;
%         resp(i).winInfo.overlapRatio = overlapRatio;
%         resp(i).winInfo.type = winType;
%         wMs = floor(durMs/(nWins+overlapRatio-nWins*overlapRatio));
%         resp(i).winInfo.durMs = wMs;
%         wSamps = wMs * spms;
%         resp(i).winInfo.durSamps = wSamps;
%         overlapSamps = floor(wSamps * overlapRatio);
%         resp(i).winInfo.overlapSamps = overlapSamps;
%         overlapMs = overlapSamps * spms;
%         resp(i).winInfo.overlapMs = overlapMs;
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
        wMs = floor(durMs/(nWins+overlapRatio-nWins*overlapRatio));
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
    
    % Get the DFT and PSD; manually calculate ESD (I want to see if it's
    % the same as what they use for the PSD).
    wav = stim.wav;
    w = blackman(wSamps,'periodic');
    o = overlapSamps;
    f = fInterest(1):fStep:fInterest(2);
    fs = stim.fs;
    [dft,f,t,psd] = spectrogram(wav,w,o,f,fs); 
    
    bw = struct;
    for j=1:nWins
        offset = j*(wSamps - overlapSamps);
        bw(j).ind = j;
        bw(j).freqVec = f;
        bw(j).centerTimeSecs = t(j);
        bw(j).startSampleIndex = 1 + offset;
        bw(j).endSampleIndex = wSamps + offset;
        d = dft(:,j);
        bw(j).dftRaw = d;
        bw(j).dftAbs = abs(d);
        bw(j).esd = d .* conj(d);
        bw(j).psd = psd(:,j); 
        [~,ind] = max(abs(d));
        bw(j).peakFreq = f(ind);
    end
    resp(i).stim.byWindow = sortfields(bw);
    clear bw
end
resp = sortfields(resp);

%% SECTION 3: FFR DFT WINDOWED
disp('SECTION 3: FFR DFT');
for i=1:numel(resp)
    prog = [num2str(i) '/' num2str(numel(resp))];
    disp(['Windowed FFR DFT : ' prog]);
    
    %% 3.1) FFR LAG
    % Some shorthand
    r = resp(i);
    fs = r.fs;
    spms = r.spms;
    pss = r.epochInfo.prestimSamps;
    ffr = shiftdim(r.wavAvg);
    stim = shiftdim(r.stim.wav);
    ffrTrimmed = ffr(pss+1:end);
    pad = zeros(length(ffrTrimmed)-length(stim),1);
    stimPadded = [stim;pad];
    
    % Get corr and lag vectors
    [corrVec,lagVec] = xcorr(ffrTrimmed,stimPadded,'none');
    % Make all correlations absolute & only allow positive lag
    corrVecPos = abs(corrVec(lagVec>0));
    % Further restrict to lag range allowed in settings
    lagSampsMin = lagMsMin * spms;
    lagSampsMax = lagMsMax * spms;
    inds = lagSampsMin:lagSampsMax;
    corrVecAllowed = corrVecPos(inds);
    
    % Get lag, but don't let it be at range endpoints, because it's
    % probably not a true peak
    [~,ind] = max(corrVecAllowed);
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
    
    resp(i).lag.samps = lag;
    resp(i).lag.ms = lag / spms;
    resp(i).lag.searchMinSamps = lagSampsMin;
    resp(i).lag.searchMaxSamps = lagSampsMax;
    resp(i).lag.searchMinMs = lagMsMin;
    resp(i).lag.searchMaxMs = lagMsMax;
    resp(i).lag = orderfields(resp(i).lag);
    resp = orderfields(resp);
    
    %% 3.2) FFR DFT
    rw = r.winInfo;
    wSamps = rw.durMs;
    w = blackman(wSamps); % vector with normalized window envelope
    o = floor(overlapRatio*wSamps); % overlap in number of samples
    f = fInterest(1):fStep:fInterest(2);
    
    stimSamps = length(r.stim.wav);
    inds = (1:stimSamps) + lag;
    [dft,f,t,psd] = spectrogram(ffrTrimmed(inds),w,o,f,fs);
    tol = peakSearchToleranceHz;
   
    % Prestim interval
    [dftPre,f,tPre,psdPre] = spectrogram(ffr(1:pss),w,0,f,fs); 
    tPre = tPre - pss/fs; % prestim times negative
    psi.freqVec = f;
    psi.centerTimeSecs = tPre;
    psi.startSampleIndex = 1;
    psi.endSampleIndex = pss;
    psi.dftRaw = dftPre;
    psi.dftAbs = abs(dftPre);
    psi.esd = dftPre .* conj(dftPre);
    psi.psd = psdPre;
    
    % Response windows
    for j=1:nWins
        rsbw = r.stim.byWindow(j);
        
        stimFreq = rsbw.peakFreq;
        fMin = stimFreq - tol;
        fMax = stimFreq + tol;
        fLogic = f>=fMin & f<=fMax;
        fAllowed = f(fLogic);
        
        nf = mean(abs(dftPre(fLogic)));
        bw(j).noiseFloorFreqMin = fMin;
        bw(j).noiseFloorFreqMax = fMax;
        
        offset = j*(wSamps - overlapSamps);
        bw(j).ind = j;
        bw(j).freqVec = f;
        bw(j).centerTimeSecs = t(j);
        bw(j).startSampleIndex = inds(1 + offset);
        bw(j).endSampleIndex = inds(wSamps + offset);
        d = dft(:,j);
        bw(j).dftRaw = d;
        bw(j).dftAbs = abs(d);
        bw(j).esd = d .* conj(d);
        bw(j).psd = psd(:,j);
        
        dAllowed = d(fLogic);
        [val,ind] = max(abs(dAllowed));
        snrLin = val / nf;
        snrDb = 10 * log10(snrLin);
        bw(j).peakFreq = fAllowed(ind);
        bw(j).peakSearchMin = fMin;
        bw(j).peakSearchMax = fMax;
        bw(j).peakMag = val;
        bw(j).snrLinear = snrLin;
        bw(j).snrDb = snrDb;
        bw(j).noiseFloorLevel = nf;
    end
    clear nf
    
    % Store
    resp(i).snrDb = mean(vertcat(bw.snrDb));
    resp(i).snrLinear = mean(vertcat(bw.snrLinear));
    resp(i).prestim = orderfields(psi);
    resp(i).byWindow = orderfields(bw);
    clear bw psi
    resp = orderfields(resp);
    
end; clear i % end loop thru response files

%% 3.3) Export
if isfield(in,'exportInfo')
    exportInfo = in.exportInfo;
    ind = numel(exportInfo) + 1;
    exportInfo(ind).functionSpec = ei.functionSpec;
    exportInfo(ind).versionSpec = ei.versionSpec;
else
    exportInfo = ei; %#ok<NASGU>
end

disp('Exporting updated resp file...');
save(outPath,'resp','exportInfo');
disp('Success!');
disp([mfilename ' completed!']);

end % end main function block (ampSnr and nested functions)