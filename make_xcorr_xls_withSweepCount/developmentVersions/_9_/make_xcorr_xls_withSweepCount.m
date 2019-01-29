function make_xcorr_xls_withSweepCount(respDir,stimDir,outPath,varargin)

% _9_ adds magnitude 

ei.functionSpec = mfilename;
ei.versionSpec = '_9_';
defaultInDir = ['C:\Users\vhapormadseb\Desktop\Molis FFR datamats\'...
    'tone glides study\newProcess2017All'];

%% SECTION 0: OPTIONS & SETTINGS
%% 0.1) Analysis Settings
defaultFs = 20000; % assumed sampling rate if not already specified in params
lagMsMin = 2;
lagMsMax = 22;
stimPolForComparison = 'rarefaction'; % shouldn't matter, so this is arbitrary
sweepCountStep = 50;
peakSearchToleranceHz = 25;
overlapRatio = 0;
fInterest = [350,750];
nWins = 3;
fStep = 1;

%% 0.2) Filter Settings
filterOrder = 400;
passBandLow = 300;
passBandHigh = 800;

%% SECTION 1: SETUP & INITIALIZATION
disp('SECTION 1: SETUP & INITIALIZATION');
alreadyFiltered = false;
loaded = false;
opts.saveMode = 'mainfolder';
for i=1:numel(varargin)
    if isequal(varargin{i},'-subfolder')
        opts.saveMode = 'subfolder';
    else
        warning('Unrecognized input argument, will be ignored');
    end
end

%% 1.1) Directories
% RESP
if ~exist('respDir','var')
    msg = ['Select FFR input folder OR hit ''cancel'' to load'...
        ' pre-processed data'];
    respDir = uigetdir(pwd,msg);
elseif isstruct(respDir)
    alreadyFiltered = true;
    in = respDir;
    if isfield(in,'resp')
        resp = in.resp;
    else
        error('Invalid data passed as response');
    end
    if isfield(in,'exportInfo')
        exportInfo = in.exportInfo;
    end
    loaded = true;
elseif ~exist(respDir,'dir')
    error('respDir folder path could not be found.');
end
if isequal(respDir,0)
    msg = 'Select MAT file with FFR data';
    [inFile,inDir] = uigetfile(defaultInDir,msg);
    inPath = [inDir '\' inFile];
    alreadyFiltered = true;
    if isequal(inFile,0)
        warning('User cancelled. Aborting function.')
        return
    elseif ~exist(inPath,'file')
        error('Response MAT file path could not be found.');
    end
end

% STIM
if ~exist('resp','var') || ~isfield(resp,'stim')
    if ~exist('stimDir','var') || isempty(stimDir)
        stimDir = uigetdir(pwd,'Select stim input folder');
    end
    if ~exist(stimDir,'dir')
        error('stimDir folder path could not be found');
    end
end

% OUTPUT
if ~exist('outPath','var')
    [outName,outDir] = uiputfile([pwd '\xcorr.xlsx'],...
        'Save output XLSX as...');
    outPath = [outDir outName];
else
    bs = find(outPath=='\');
    outDir = outPath(1:bs(end));
    outName = outPath(bs(end)+1:end);
end
if ~exist(outDir,'dir')
    warning('Output directory could not be found. Creating now...')
    mkdir(outDir);
    disp('Creation successful!')
end
if isequal(opts.saveMode,'subfolder')
    outDir = [outDir '\' mfilename];
    outPath = [outDir '\' outName];
    if ~exist(outDir,'dir')
        disp('Creating output subfolder...')
        mkdir(outDir,'dir');
        disp('Creation successful!')
    end
end

%% 1.2) Initialize filter
fd = fdesign.bandpass('N,Fc1,Fc2',filterOrder,passBandLow,passBandHigh,...
    defaultFs);
hfd = design(fd,'FIR'); clear fd

%% 1.3) Check for existing data so some steps can be skipped
disp('Checking for existing FFR data...');
% Check if filtered resp data already exists
alreadyStim = false;
if loaded
    disp('Data passed directly, no need to load.');
    if isfield(resp,'stim')
        alreadyStim = true;
    end
elseif alreadyFiltered
    disp('Found! Loading...');
    load(inPath);
    disp('Loaded successfully!');
    if isfield(resp,'stim')
        alreadyStim = true;
    end
else
    disp('Not found! Will continue as usual...');
    alreadyFiltered = false;
end

%% SECTION 2: STIMULUS ANALYSIS
disp('SECTION 2: STIMULUS ANALYSIS')
if alreadyStim
    disp('Skipping this section, stim data already loaded.')
else
    disp('Checking for stim directory...');
    if isequal(stimDir,0)
        warning('User did not specify. Aborting function.')
        return
    end
    %% 2.1) Listing WAV Files
    disp('2.1) Listing .wav files...');
    d = dir([stimDir '\*.wav']);
    stimList = {d.name};
    stimList = stimList(:);
    
    %% 2.2) Load WAV Files
    disp('2.2) Loading .wav files...');
    clear stim
    for i=1:numel(stimList)
        fName = stimList{i};
        disp(fName);
        [wav,fs] = audioread([stimDir '\' fName]);
        % account for stereo files
        if ~isvector(wav)
            if isequal(wav(:,1),wav(:,2))
                wav = wav(:,1);
            else
                error(['Unrecognized wav format for file ' fName]);
            end
        end
        temp = parse_stim_fname(fName);
        
        % Resample stim if samp rate doesn't match that of response
        if fs == defaultFs
            temp.wav = wav;
            temp.fs = fs;
            temp.origFs = fs;
            temp.wasResampled = false;
        else
            temp.wav = resample(wav,defaultFs,fs);
            temp.origFs = fs;
            temp.fs = defaultFs;
            temp.wasResampled = true;
        end
        
        temp.fName = fName;
        temp.fPath = [stimDir '\' fName];
        if ~exist('stim','var')
            stim = temp;
        else
            stim(i) = temp; %#ok<*AGROW>
        end
    end; clear i
end % end if

%% SECTION 3: FFR DATA PREP
disp('SECTION 3: FFR DATA PREP');

%% 3.1) Get list of files
disp('3.1) Getting list of files')
if alreadyFiltered
    disp('Skipping this section, FFR data already loaded');
else
    
    pathList = {};
    fList = {};
    subfolds = get_all_subfolders(respDir);
    for i=1:numel(subfolds)
        currDir = subfolds{i};
        temp = getfiles(currDir,'*.mat');
        fList = vertcat(fList,temp); %#ok<*AGROW>
        paths = strcat(currDir, '\', temp);
        pathList = vertcat(pathList,paths);
    end; clear i
end

%% 3.2) Parse filenames and params structs for subject and condition info
disp('3.2) Parsing filenames and params');
if alreadyFiltered
    disp('Skipping this section, data already loaded.');
else
    clear resp % just in case it exists in workspace for some reason
    for i=1:numel(fList)
        disp(['Parsing ' num2str(i) '/' num2str(numel(fList)) ': ' fList{i}]);
        temp = parse_resp_fname(fList{i});
        temp.fPath = pathList{i};
        temp.fName = fList{i};
        
        if ~exist('resp','var')
            resp = temp;
        else
            resp(i) = temp;
        end; clear temp
        
        temp = load(pathList{i},'params');
        temp.params = sortfields(temp.params);
        if  ~exist('params','var')
            params = temp.params;
        else
            params(i) = temp.params;
        end
    end; clear i temp
    if exist('paramRmList','var')
        if ~isempty(paramRmList)
            params = rmfield(params,paramRmList);
        end
    end
    resp = mergestructs(resp,params);
    clear fList
end

%% 3.3) Make averages for different sweep counts
if alreadyFiltered
    disp('Skipping this section, data already loaded');
else
    disp('3.3) Make averages for different sweep counts');
    for i=1:numel(resp)
        r = resp(i);
        prog = [num2str(i) '/' num2str(numel(resp))];
        disp([prog ' ' r.fName]);
        temp = load(r.fPath,'epochs');
        temp = cast(temp.epochs,'double');
        np = r.epochInfo.nPoints;
        pInd = find(np==size(temp));
        switch pInd
            case 1
                % no alteration needed
            case 2
                temp = transpose(temp);
            otherwise
                error('invalid dimension index (pInd)');
        end
        
        % Average waveforms at different sweep counts
        nAcc = size(temp,2);
        step = sweepCountStep;
        sweepVec = step:step:nAcc;
        sc = struct;
        for j=1:numel(sweepVec)
            n = sweepVec(j);
            disp([prog ' ' r.fName ' ' num2str(n)]);
            sc(j).count = n;
            wav = temp(:,1:n);
            sc(j).wavAvg = mean(filtfilt(hfd.Numerator,1,wav),2);
        end; clear j
        resp(i).bySweepCount = sc; % copy back into main struct
        
        % Average waveform with all available sweeps
        resp(i).wavAvg = mean(filtfilt(hfd.Numerator,1,temp),2);
        
    end; clear i
    clear temp np pInd r step sweepVec nAcc sc wav
    resp = sortfields(resp);
end

%% 3.4) Map each resp onto its corresponding stim and export updated resp
disp('3.4) Deriving additional variables');
if alreadyStim
    disp('Skipping this section, data already loaded');
else
    sd = {stim.direction};
    sp = {stim.polarity};
    ss = {stim.slope};
    sdur = vertcat(stim.durationMs);
    polLogic = shiftdim(strcmpi(stimPolForComparison,sp));
    
    for i=1:numel(resp)
        disp(resp(i).fName);
        resp(i).spms = resp(i).fs / 1000;
        resp(i).epochInfo.prestimSamps = resp(i).epochInfo.prestimMs * resp(i).spms;
        r = resp(i);
        
        % control for if stim duration was accidentally coded as string
        rsd = r.stim.durationMs;
        if ischar(rsd)
            resp(i).stim.durationMs = str2double(rsd);
        end
        clear rsd
        r = resp(i);
        
        dirLogic = shiftdim(strcmpi(r.stim.direction,sd));
        slopeLogic = shiftdim(strcmpi(r.stim.slope,ss));
        durLogic = shiftdim(r.stim.durationMs == sdur);
        comboLogic = polLogic & dirLogic & slopeLogic & durLogic;
        assert(sum(comboLogic)==1,[r.fName ': did not uniquely ID stim']);
        
        match = regexp(r.fName,'(plus|minus)','match','once');
        if ~isempty(match)
            resp(i).stim.polarity = match;
            r = resp(i);
        end
        
        % map each resp onto its corresponding stim
        temp = stim(comboLogic);
        flds = fieldnames(temp);
        flds2 = fieldnames(r.stim);
        for j=1:numel(flds)
            fld = flds{j};
            fldLogic = strcmpi(fld,flds2);
            switch sum(fldLogic)
                case 0
                    resp(i).stim.(fld) = temp.(fld);
                case 1
                    fld2 = flds2{fldLogic};
                    if isequal(r.stim.(fld2),temp.(fld));
                        if ~isequal(fld2,fld)
                            msg = [fld ', ' fld2 ': capitalization mismatch'];
                            warning(msg);
                        end
                    else
                        if isempty(r.stim.(fld2))
                            resp(i).stim.(fld2) = temp.(fld);
                            if ~isequal(fld,fld2)
                                msg = [fld ', ' fld2 ': ',...
                                    'capitalization mismatch'];
                                warning(msg);
                            end
                        else
                            switch fld
                                case 'polarity'
                                    % do nothing, doesn't matter
                                otherwise
                                    msg = [r.fName ': ' fld ': data mismatch'];
                                    error(msg);
                            end
                        end
                    end
                otherwise
                    error([r ': ' fld ': more than one field matched']);
            end
            
        end; clear j
        clear fld* msg
        
    end; clear i r sd st sp ss dirLogic typeLogic comboLogic
    
    resp = sortfields(resp);
end % end if

%% SECTION 4: XCORR ANALYSIS
disp('SECTION 4: XCORR ANALYSIS');
for i=1:numel(resp)
    r = resp(i);
    sc = r.bySweepCount;
    disp(['Analyzing ' num2str(i) '/' num2str(numel(resp)) ': ' r.fName]);
    for j=1:numel(sc)
        %% 4.1) Prepare stim and response
        wav = sc(j).wavAvg(r.epochInfo.prestimSamps+1:end);
        
        % Pad stim to match resp
        pad = zeros(length(wav)-length(r.stim.wav),1);
        stimPadded = [r.stim.wav;pad];
        
        % Not necessary to scale here because we're only interested in
        % the relative values for now.
        
        % Derive some vars needed for next chunk.
        lagSampsMin = r.spms*lagMsMin;
        lagSampsMax = r.spms*lagMsMax;
        
        % Clean up from this section
        clear scaleFactor pad idx stimCorr lagVecSamps
        
        %% 4.2) Calculate how far response lags behind stim
        inds = lagSampsMin:lagSampsMax;
        oLagSamps(j) = get_olag(wav,stimPadded,inds); % helper function
        oLagMs(j) = oLagSamps(j) / r.spms;
        
        clear inds lagSampsMin lagSampsMax
        
        %% 4.3) Calculate whole-epoch correlation coeff
        % Keep only the part of the resp that lines up with the stim
        keepInds = (1:length(r.stim.wav)) + oLagSamps(j);
        try
            wavTrimmed = wav(keepInds);
        catch
            msg = [num2str(i) ', ' num2str(j) ': ' r.fName];
            disp(msg);
            pause
        end
        
        % Rescale using stim-length inputs
        scaleFactor = rms(r.stim.wav) / rms(wavTrimmed);
        wavTrimmedScaled = wavTrimmed .* scaleFactor;
        
        % Run cross-correlation
        [corrs,lags] = xcorr(wavTrimmedScaled,r.stim.wav,'none');
        
        % Resp has already been aligned to the correct spot on the stim so
        % just take correlation at lag zero.
        corr = abs(corrs(lags==0));
        clear corrs lags
        
        % Standardize to convert from raw to coefficient
        [corrs,lags] = xcorr(r.stim.wav,r.stim.wav,'none');
        [ref,idx] = max(abs(corrs));
        assert(lags(idx)==0,...
            [r.stim.fName ' does not have autocorr peak at lag zero.']);
        oCorr(j) = corr / ref;
        
        % Store results
        sc(j).oCorr = oCorr(j);
        sc(j).oLagSamps = oLagSamps(j);
        sc(j).oLagMs = oLagMs(j);
    end; clear j
    
    % Store results
    resp(i).bySweepCount = sc;
    
    % Clean up from this section
    clear keepInds wavTrimmedScaled corrs lags ref
    
end; clear i
%% SECTION 5: NOISE FLOOR ANALYSIS
disp('SECTION 5: NOISE FLOOR ANALYSIS');
for i=1:numel(resp)
    r = resp(i);
    sc = r.bySweepCount;
    disp(['Analyzing ' num2str(i) '/' num2str(numel(resp)) ': ' r.fName]);
    for j=1:numel(sc)
        %             disp(num2str(sc(j).count));
        pre = sc(j).wavAvg(1:r.epochInfo.prestimSamps);
        nf = rms(pre);
        resp(i).bySweepCount(j).noiseFloor = nf;
    end; clear j
end; clear i

%% SECTION 6: MAGNITUDE ANALYSIS
disp('SECTION 6: MAGNITUDE ANALYSIS');
for i=1:numel(resp)
    r = resp(i);
    sc = r.bySweepCount;
    disp(['Analyzing ' num2str(i) '/' num2str(numel(resp)) ': ' r.fName]);
    tic
    for j=1:numel(sc)
        
        wav = sc(j).wavAvg(r.epochInfo.prestimSamps+1:end);
        lagSamps = sc(j).oLagSamps;
        spms = r.spms;
        fs = r.fs;
        
        % Keep only the part of the resp that lines up with the stim
        keepInds = (1:length(r.stim.wav)) + lagSamps;
        try
            wavTrimmed = wav(keepInds);
        catch
            msg = [num2str(i) ', ' num2str(j) ': ' r.fName];
            disp(msg);
            pause
        end
        
        rw = r.winInfo;
        wSamps = rw.durMs * spms;
        w = blackman(wSamps); % vector with normalized window envelope
        o = floor(overlapRatio*wSamps); % overlap in number of samples
        f = fInterest(1):fStep:fInterest(2);
        
        [dft,f,t,psd] = spectrogram(wavTrimmed,w,o,f,fs);
        tol = peakSearchToleranceHz;
        
        for k=1:nWins
            rsbw = r.stim.byWindow(k);
            stimFreq = rsbw.peakFreq;
            fMin = stimFreq - tol;
            fMax = stimFreq + tol;
            fLogic = f>=fMin & f<=fMax;
            fAllowed = f(fLogic);
            
            offset = (k-1)*(wSamps - o);
            bw(k).ind = k;
            bw(k).freqVec = f;
            bw(k).centerTimeSecs = t(k);
            bw(k).startSampleIndex = keepInds(1 + offset);
            bw(k).endSampleIndex = keepInds(wSamps + offset);
            d = dft(:,k);
            bw(k).dftRaw = d;
            bw(k).dftAbs = abs(d);
            p2 = abs(d/fs);
            bw(k).dftUv = [p2(1); 2*p2(2:end-1); p2(end)] * fs/wSamps;
            bw(k).esd = d .* conj(d);
            bw(k).psd = psd(:,k);
            
            dAllowed = bw(k).dftUv(fLogic);
            [val,ind] = max(abs(dAllowed));
            bw(k).peakFreq = fAllowed(ind);
            bw(k).peakSearchMin = fMin;
            bw(k).peakSearchMax = fMax;
            bw(k).peakUv = val;
                
        end; clear k
        
        resp(i).bySweepCount(j).byWindow = bw;
        clear bw
        resp = sortfields(resp);
        r = resp(i);
        
    end; clear j
    toc
end; clear i


%% SECTION 7: EXPORT
disp('SECTION 7: EXPORT')
%% 6.1) Export Matlab files
disp('6.1) Export Matlab files');
if exist('exportInfo','var')
    ind = numel(exportInfo)+1;
    exportInfo(ind).functionSpec = ei.functionSpec;
    exportInfo(ind).versionSpec = ei.versionSpec;
    exportInfo(ind).timestamp = datetime;
else
    exportInfo = ei;
    exportInfo.timestamp = datetime;
end
resp = sortfields(resp);
disp('Exporting response data with stim and xcorr data built in...');
outMat = ['resp_' outName(1:end-5) '.mat'];
save([outDir '\' outMat],'resp','exportInfo','-v7.3');
disp('Exported!');

%% 6.2) Export as XLSX file
disp('6.2) Export as .xlsx file');
headerRow = {...
    'ID',...
    'Number',...
    'Group',...
    'Stimulus',...
    'Slope',...
    'Direction',...
    'Duration',...
    'IsArtifact',...
    'Condition',...
    'Ear',...
    'RefElec',...
    'ActiveElec',...
    'Sweeps',...
    'CorrelationCoeff',...
    'Lag(ms)',...
    'NoiseFloor',...
    'PeakAmp_Win1(uV)',...
    'PeakAmp_Win2(uV)',...
    'PeakAmp_Win3(uV)'};

for i=1:numel(resp)
    r = resp(i);
    for j=1:numel(r.bySweepCount)
        rsc = r.bySweepCount(j);
        
        id{j,i} = r.subject.id;
        number{j,i} = r.subject.serialNumberAsDouble;
        group{j,i} = r.subject.group;
        stimulus{j,i} = r.stim.fName;
        slope{j,i} = r.stim.slope;
        direction{j,i} = r.stim.direction;
        duration{j,i} = r.stim.durationMs;
        art{j,i} = r.isArtifact;
        condition{j,i} = r.stim.polarity;
        ear{j,i} = r.stim.laterality;
        ref{j,i} = r.montage.reference;
        targ{j,i} = r.montage.channelToAnalyze;
        sweeps{j,i} = rsc.count;
        correlation{j,i} = rsc.oCorr;
        lag{j,i} = rsc.oLagMs;
        noiseFloor{j,i} = rsc.noiseFloor;
        paw1{j,i} = rsc.byWindow(1).peakUv;
        paw2{j,i} = rsc.byWindow(2).peakUv;
        paw3{j,i} = rsc.byWindow(3).peakUv;
    end; clear j
end; clear i

v = ~cellfun('isempty',id);

outCell = [id(v), number(v), group(v), stimulus(v), slope(v),...
    direction(v), duration(v), art(v), condition(v), ear(v), ref(v),...
    targ(v), sweeps(v), correlation(v), lag(v), noiseFloor(v), paw1(v),...
    paw2(v), paw3(v)];
outCell = [headerRow; outCell];
xlswrite(outPath,outCell);
disp('Success!');
disp([mfilename ' completed!']);

% End main function body
%% NESTED FUNCTIONS
%---------------------------BEGIN NESTED FUNCTIONS-------------------------
% These functions can read from and write to variables in the main function
% even if they are not explicitly passed to the nested function.
%--------------------------------------------------------------------------

% ---------------------------END NESTED FUNCTIONS--------------------------
end % end main function block

%% HELPER FUNCTIONS
%----------------------------BEGIN HELPER FUNCTIONS------------------------
% These functions only have read-access to variables that are explicitly
% passed as input arguments, and they only have write-access to variables
% that are explicitly passed as output arguments.
%--------------------------------------------------------------------------
%% parse_resp_fname
function out = parse_resp_fname(str)

out = struct;

% base type
out.stim.type = 'tone';
% subject ID
expr = '(YN||NH||HI)\d{3}';
out.subject.id = regexp(str,expr,'match','once');
% group
out.subject.group = regexp(str,'YN||NH||HI','match','once');
% serial number
sid = out.subject.id;
sn = regexp(sid,'\d{3}','match','once');
out.subject.serialNumberAsString = sn;
out.subject.serialNumberAsDouble = str2double(sn);
% presentation ear
expr = '_(R||L)_';
toks = regexp(str,expr,'tokens','once');
while iscell(toks)
    toks = toks{1};
end
out.stim.laterality = toks;
% slope
expr = '_(10||23||13)_';
toks = regexp(str,expr,'tokens','once');
while iscell(toks)
    toks = toks{1};
end
out.stim.slope = toks;
% direction
expr = 'up||down||flat';
out.stim.direction = regexp(str,expr,'match','once');
% artifact?
expr = '_art_';
if isempty(regexp(str,expr,'once'))
    out.isArtifact = false;
else
    out.isArtifact = true;
end
% reference
expr = '_(C7||A1||A2)_';
toks = regexp(str,expr,'tokens','once');
while iscell(toks) && ~isempty(toks)
    toks = toks{1};
end
if isempty(toks)
    out.montage.reference = 'A2';
else
    out.montage.reference = toks;
end
% plus/minus/rare/cond
expr = 'plus||minus||_13_\w{2}(_bc_ar)?\.mat||_46_\w{2}(_bc_ar)?\.mat';
match = regexp(str,expr,'match','once');
switch match
    case {'plus','minus'}
        out.stim.polarity = match;
    otherwise
        if ~isempty(regexp(match,'13','once'))
            out.stim.polarity = 'condensation';
        elseif ~isempty(regexp(match,'46','once'))
            out.stim.polarity = 'rarefaction';
        else
            error('unrecognized polarity');
        end
end

end % end parse_resp_fname function
%--------------------------------------------------------------------------
%% parse_stim_fname
function out = parse_stim_fname(str)

% NOTING HERE SO I DON'T FORGET
% Condensation = 1-3 = original
% Rarefaction = 4-6 = inverted

% base type
out.type = 'tone';
% parse filename
expr = ['^(?<inv>inv_)?s_(?<f1>\d{3})_(?<f2>\d{3})_(?<dur>\d{2,3})'...
    '(?:_\d+kHz)?\.wav'];
s = regexp(str,expr,'names');
% polarity
if ~isempty(s.inv)
    out.polarity = 'rarefaction';
else
    out.polarity = 'condensation';
end
% starting freq
out.startFreq = str2double(s.f1);
% ending freq
out.endFreq = str2double(s.f2);
% duration
out.durationMs = str2double(s.dur);
% direction
switch true
    case out.startFreq > out.endFreq
        out.direction = 'down';
    case out.startFreq < out.endFreq
        out.direction = 'up';
    case out.startFreq == out.endFreq
        out.direction = 'flat';
    otherwise
        error([str ': direction: something went wrong']);
end
% slope
switch out.startFreq
    case {354,707}
        out.slope = '10';
        out.slopeLong = 'Full octave';
    case {398,629}
        out.slope = '23';
        out.slopeLong = '2/3 octave';
    case {446,561}
        out.slope = '13';
        out.slopeLong = '1/3 octave';
    otherwise
        error([ str ': slope: unrecognized start freq']);
end

end % end parse_stim_fname function
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

end
%--------------------------------------------------------------------------
%% get_olag
function oLag = get_olag(x1,x2,inds)
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
oLag = ind + lagSampsMin - 1;

end % end get_olag function
%------------------------------END HELPER FUNCTIONS------------------------
% [EOF]
