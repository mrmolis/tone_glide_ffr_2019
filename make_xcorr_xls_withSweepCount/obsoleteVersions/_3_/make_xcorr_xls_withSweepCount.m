function make_xcorr_xls_withSweepCount

exportInfo.functionSpec = mfilename;
exportInfo.versionSpec = '_3_';

%% SECTION 0: OPTIONS & SETTINGS
%% 0.1) Folder Settings
% INPUT PATHS
respDir = ['C:\Users\vhapormadseb\Desktop\Molis FFR datamats\'...
    'tone glides study\newProcess2017All\CztoA2'];
stimDir = ['C:\Users\vhapormadseb\Desktop\Molis FFR stims\'...
    'tone glides study\44.1 kHz'];

% OUTPUT PATH
saveDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
    'Methods Manuscript\Matlab_Output'];

%% 0.2) Analysis Settings
defaultFs = 20000; %assumed sampling rate if not already specified in params
lagMsMin = 4;
lagMsMax = 12;
stimPolForComparison = 'rarefaction'; % shouldn't matter, so this is arbitrary
sweepCountStep = 50;

%% 0.3) Filter Settings
filterOrder = 400;
passBandLow = 300;
passBandHigh = 800;

%% SECTION 1: SETUP & INITIALIZATION
disp('SECTION 1: SETUP & INITIALIZATION');
%% 1.1) Prep output directory
exportInfo.functionSpec = mfilename;
exportInfo.versionSpec = '_3_';
outDir = [saveDir '\' mfilename];
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% 1.2) Initialize filter
fd = fdesign.bandpass('N,Fc1,Fc2',filterOrder,passBandLow,passBandHigh,...
    defaultFs);
hfd = design(fd,'FIR'); clear fd

%% 1.3) Intialize other variables
exportInfo.functionSpec = mfilename;

%% SECTION 2: STIMULUS ANALYSIS
disp('SECTION 2: STIMULUS ANALYSIS')
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

%% SECTION 3: FFR DATA PREP
disp('SECTION 3: FFR DATA PREP');
disp('Checking for existing data...');
% Check if filtered resp data already exists
if exist([outDir '\filteredRespStruct.mat'],'file')
    alreadyFiltered = true;
    disp('Found! Loading...');
else
    disp('Not found! Will continue as usual...');
    alreadyFiltered = false;
end

%% 3.1) Get list of files
if ~alreadyFiltered
    disp('3.1) Getting list of files')
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
if ~alreadyFiltered
    disp('3.2) Parsing filenames and params');
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
        temp.params = orderfields(temp.params);
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
if ~alreadyFiltered
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
        r.wavAvg = mean(filtfilt(hfd.Numerator,1,temp),2);
        
    end; clear i
    clear temp np pInd r step sweepVec nAcc sc wav
    resp = orderfields(resp);
    
    disp('Exporting filtered response data with sweep count breakdown...');
    exportInfo.timestamp = datetime;
    save([outDir '\filteredRespStruct.mat'],'resp','exportInfo');
    disp('Exported!')
else
    load([outDir '\filteredRespStruct.mat']);
    disp('Loaded!');
end

%% 3.4) Map each resp onto its corresponding stim and export updated resp
disp('3.4) Deriving additional variables');
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

resp = orderfields(resp);
exportInfo.timestamp = datetime;
disp('Exporting updated response data...');
save([outDir '\filteredRespStruct.mat'],'resp','exportInfo');
disp('Exported!');

%% SECTION 4: XCORR ANALYSIS
disp('SECTION 4: XCORR ANALYSIS');
for i=1:numel(resp)
    r = resp(i);
    sc = r.bySweepCount;
    disp(['Analyzing ' num2str(i) '/' num2str(numel(resp)) ': ' r.fName]);
    for j=1:numel(sc)
        disp(num2str(sc(j).count));
        %% 4.1) Prepare stim and response
        wav = sc(j).wavAvg(r.epochInfo.prestimSamps+1:end);
        
        % Pad stim to match resp
        pad = zeros(length(wav)-length(r.stim.wav),1);
        stimPadded = [r.stim.wav;pad];
        
        % Scale resp to rms of stim, to eliminate scale effects on raw
        % correlations.
        scaleFactor = rms(stimPadded)/rms(wav);
        wavScaled = wav .* scaleFactor;
        
        % Derive some vars needed for next chunk.
        lagSampsMin = r.spms*lagMsMin;
        lagSampsMax = r.spms*lagMsMax;
        
        % Clean up from this section
        clear scaleFactor pad idx stimCorr lagVecSamps
        
        %% 4.2) Calculate how far response lags behind stim
        inds = lagSampsMin:lagSampsMax;
        oLagSamps(j) = get_olag(wavScaled,stimPadded,inds); % helper function
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
        corr = corrs(lags==0);
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
        disp(num2str(sc(j).count));
        pre = sc(j).wavAvg(1:r.epochInfo.prestimSamps);
        nf = rms(pre);
        resp(i).bySweepCount(j).noiseFloor = nf;
    end; clear j
end; clear i

%% SECTION 6: EXPORT
disp('SECTION 6: EXPORT')
%% 6.1) Export Matlab files
disp('6.1) Export Matlab files');
exportInfo.timestamp = datetime;
resp = orderfields(resp); 
disp('Exporting response data with stim and xcorr data built in...');
save([outDir '\respWithXCorr.mat'],'resp','exportInfo');
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
    'Condition',...
    'Ear',...
    'Sweeps',...
    'CorrelationCoeff',...
    'Lag',...
    'NoiseFloor'};

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
        condition{j,i} = r.stim.polarity;
        ear{j,i} = r.stim.laterality;
        sweeps{j,i} = rsc.count;
        correlation{j,i} = rsc.oCorr;
        lag{j,i} = rsc.oLag;
        noiseFloor{j,i} = rsc.noiseFloor;
    end; clear j
end; clear i

outcell = [id(:), number(:), group(:), stimulus(:), slope(:),...
    direction(:), duration(:), condition(:), ear(:), sweeps(:),...
    correlation(:), lag(:), noiseFloor(:)];
outcell = [headerRow; outcell];
outname = 'xcorr.xlsx';
outpath = [outDir '\' outname];
xlswrite(outpath,outcell);

toc

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
expr = 'plus||minus||_13_\w{2}\.mat||_46_\w{2}\.mat';
match = regexp(str,expr,'match','once');
switch match
    case '_13.mat'
        out.stim.polarity = 'condensation';
    case '_46.mat'
        out.stim.polarity = 'rarefaction';
    otherwise
        out.stim.polarity = match;
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
%% portfields
function out = portfields(dest,src)

flds = fieldnames(src);
for I=1:numel(flds)
    dest.(flds{I}) = src.(flds{I});
end; clear j
out = dest;

end % end portfields function
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
oLag = ind + lagSampsMin;

end % end get_olag function
%------------------------------END HELPER FUNCTIONS------------------------
% [EOF]
