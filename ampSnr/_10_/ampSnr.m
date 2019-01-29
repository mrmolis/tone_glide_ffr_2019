function outPath = ampSnr(inPath,outPath,exportInfo)

ei.functionSpec = mfilename;
ei.versionSpec = '_10_';

%% SECTION 0: USER SETTINGS
disp('SECTION 0: USER SETTINGS');

%% 0.1) Folder Settings
defaultOutDir = ['C:\Users\vhapormadseb\Desktop\tempOutput\' mfilename];
defaultInDir = ['\\r01porhsm03\Research\NCRAR\Spectral_Dynamics_Grant\'...
    'FFR\Tone Glide Study\Methods Manuscript\xcorr_newProcess2017All'];
defaultOutName = 'resp_snr.mat';

%% 0.2) Analysis Settings
fInterest = [350,750]; % set frequency range of interest
fStep = 1; % space between frequencies in Hz
nWins = 3; % number of desired windows
winType = 'blackman';
overlapRatio = 0;
peakSearchToleranceHz = 25;
lagMsMin = 2;
lagMsMax = 22;

%% SECTION 1: INITIALIZATION
disp('SECTION 1: INITIALIZATION');

% Handle missing arguments
if abs(nargin) < 2
    outPath = '';
end
if abs(nargin) < 1
    inPath = '';
end

% Getting path to input data
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
    resp = inPath;
    loadNow = false;
else
    loadNow = true;
end

% Getting save path
defaultOutPath = [defaultOutDir '\' defaultOutName];
if isempty(outPath)
    if ~exist(defaultOutDir,'dir')
        mkdir(defaultOutDir);
    end
    [outName,outDir] = uiputfile(defaultOutPath,'Save output as...');
    if outName == 0
        warning('No output dir specified; default will be used');
        outPath = defaultOutPath;
    else
        outPath = [outDir '\' outName];
    end
else
    bs = find(outPath=='\');
    outDir = outPath(1:bs(end));
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
end

% Load data
if loadNow
    disp('Loading FFR data...')
    in = load(inPath);
    resp = in.resp;
    exportInfo = in.exportInfo;
    disp('Loading successful!');
end

%% SECTION 2: STIM DFT WINDOWED
disp('SECTION 2: STIM DFT WINDOWED');
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
disp('SECTION 3: FFR DFT WINDOWED');
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
    wSamps = rw.durMs * spms;
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
    
    if isfield(r,'byWindow')
        bw = r.byWindow;
    end
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
        
        offset = (j-1)*(wSamps - overlapSamps);
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
    resp(i).prestim = sortfields(psi);
    resp(i).byWindow = sortfields(bw);
    clear bw psi
    resp = sortfields(resp);
    
end; clear i % end loop thru response files

%% SECTION 4: EXPORT
disp('SECTION 4: EXPORT')
if exist('exportInfo','var')
    ind = numel(exportInfo) + 1;
    exportInfo(ind).functionSpec = ei.functionSpec;
    exportInfo(ind).versionSpec = ei.versionSpec;
    exportInfo(ind).timestamp = datetime;
else
    exportInfo = ei; 
    exportInfo.timestamp = datetime;
end

disp('Exporting updated resp file...');
save(outPath,'resp','exportInfo');
disp('Success!');
disp([mfilename ' completed!']);

end % end main function block (ampSnr and nested functions)