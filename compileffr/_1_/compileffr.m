function outPath = compileffr(respDir,stimDir,outPath)
ei.functionSpec = mfilename;
ei.versionSpec = '_1_';

%% SECTION 0: OPTIONS & SETTINGS
disp('SECTION 0: OPTIONS & SETTINGS');

% Folder preferences
defRespDir = ['C:\Users\vhapormadseb\Desktop\Molis FFR dat files\'...
    'tone glide study'];
defStimDir = ['C:\Users\vhapormadseb\Desktop\Molis FFR stims\'...
    'tone glides study\44.1 kHz'];
defOutPath = 'C:\Users\vhapormadseb\Desktop\tempOutput';

% Filter settings
filtOrder = 400;
passBandLow = 300;
passBandHigh = 800;
filtType = 'FIR';
filtShape = 'bandpass';

% Stim-to-resp mapping settings
stimPolForComparison = 'rarefaction'; % arbitrary, but need to pick s/t

%% SECTION 1: SETUP
disp('SECTION 1: SETUP')

% Add paths to helper functions
addpath('Z:\NCRAR_Share\Madsen\Matlab Code\mergestructs\_1_');

% Getting ffr input dir
disp('Getting FFR input folder...')
if ~exist('respDir','var')
    respDir = uigetdir(defRespDir,'Select FFR data folder');
    if isequal(respDir,0)
        warning('No FFR data selected. Function will abort.');
        outPath = 0;
        return
    end
end

% Getting stim input dir
disp('Getting stimulus input folder...')
if ~exist('stimDir','var')
    stimDir = uigetdir(defStimDir,'Select stimulus folder');
    if isequal(stimDir,0)
        warning('No stimulus data selected. Function will abort.');
        outPath = 0;
        return
    end
end

% Getting save path
disp('Getting savepath...')
if ~exist('outPath','var')
    [outFile,outDir] = uiputfile(defOutPath);
    if isequal(outDir,0)
        warning('No savepath entered. Default savepath will be used.');
        outPath = defOutPath;
    else
        outPath = [outDir '\' outFile];
    end
else
    bs = find(outPath=='\');
    outDir = outPath(1:bs(end));
end

% Filter setup
fd = fdesign.(filtShape)('N,Fc1,Fc2',filtOrder,passBandLow,...
    passBandHigh,defaultFs);
hfd = design(fd,filtType); clear fd
temp.denominator = 1;
temp.type = filtType;
temp.shape = filtShape;
temp.passBandLow = passBandLow;
temp.passBandHigh = passBandHigh;
temp.order = filtOrder;
hfd.params = temp;
clear temp

%% SECTION 2: COMPILE FFR DATA
disp('SECTION 2: COMPILE FFR DATA')

% Make list of stim WAV files
disp('Listing stim files (WAV)...');
d = dir([stimDir '\*.wav']);
stimList = {d.name};
stimList = stimList(:);

% Load WAV Files
disp('Loading stim files (WAV)...');
clear stim
for i=1:numel(stimList)
    prog = [num2str(i) '/' num2str(numel(stimList))];
    fName = stimList{i};
    msg = ['Loading stim ' prog ': ' fName];
    disp(msg);
    [wav,fs] = audioread([stimDir '\' fName]);
    % Account for stereo files
    if ~isvector(wav)
        if isequal(wav(:,1),wav(:,2))
            wav = wav(:,1);
        else
            error(['True multichannel stimuli are unsupported: ' fName]);
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

% List FFR files
disp('Listing FFR files (MAT)...')
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

% Parse filenames and params...
disp('Parsing FFR filenames and params structs...')
clear resp % just in case it exists in workspace for some reason
for i=1:numel(fList)
    prog = [num2str(i) '/' num2str(numel(fList))];
    disp(['Parsing FFR metadata ' prog ': ' fList{i}]);
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

% Signal averaging...
disp('Performing signal averaging...')
for i=1:numel(resp)
    r = resp(i);
    prog = [num2str(i) '/' num2str(numel(resp))];
    disp(['Averaging ' prog ': ' r.fName]);
    temp = load(r.fPath,'exportInfo');
    resp(i).exportInfo = temp.exportInfo;
    temp = load(r.fPath,'epochs');
    temp = cast(temp.epochs,'double');
    % Ensure correct orientation:
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
    temp = temp(filtfilt(hfd.Numerator,1,temp));
    % Store:
    resp(i).wavAvg = mean(temp,2);
    resp(i).filter = hfd;
    clear temp r
end; clear i

%% SECTION 3: MAP STIMULI ONTO RESPONSES
disp('SECTION 3: MAP STIMULI ONTO RESPONSES')

sd = {stim.direction};
sp = {stim.polarity};
ss = {stim.slope};
sdur = vertcat(stim.durationMs);
polLogic = shiftdim(strcmpi(stimPolForComparison,sp));

for i=1:numel(resp)
    prog = [num2str(i) '/' num2str(numel(resp))];
    disp(['Mapping ' prog ': ' r.fName]);
    resp(i).spms = resp(i).fs / 1000;
    resp(i).epochInfo.prestimSamps = resp(i).epochInfo.prestimMs * resp(i).spms;
    r = resp(i);
    
    % Control for if stim duration was accidentally coded as string:
    rsd = r.stim.durationMs;
    if ischar(rsd)
        resp(i).stim.durationMs = str2double(rsd);
    end
    clear rsd
    r = resp(i);
    
    % Find stim to map:
    dirLogic = shiftdim(strcmpi(r.stim.direction,sd));
    slopeLogic = shiftdim(strcmpi(r.stim.slope,ss));
    durLogic = shiftdim(r.stim.durationMs == sdur);
    comboLogic = polLogic & dirLogic & slopeLogic & durLogic;
    assert(sum(comboLogic)==1,[r.fName ': did not uniquely ID stim']);
    
    % Catch possible failure to categorize plus/minus properly:
    match = regexp(r.fName,'(plus|minus)','match','once');
    if ~isempty(match)
        resp(i).stim.polarity = match;
        r = resp(i);
    end
    
    % Perform mapping, while resolving/flagging potential mismatches:
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
end; clear i

%% SECTION 4: EXPORT
% Add to log
disp('Adding current function to log...')
ei.timestamp = datetime;
if exist('exportInfo','var')
    ind = numel(exportInfo) + 1; %#ok<NODEF>
    exportInfo(ind).functionSpec = ei.functionSpec;
    exportInfo(ind).versionSpec = ei.versionSpec;
    exportInfo(ind).timestamp = ei.timestamp;
else
    exportInfo = ei; %#ok<NASGU>
end
disp('Added!');

% Check if output dir exists; create if needed
disp('Checking for output dir...')
if ~exist(outDir,'dir')
    disp('Not found! Attempting to create...')
    mkdir(outDir);
    disp('Creation successful!')
else
    disp('Found!');
end

% Save
disp('Saving compiled FFR data...');
saveas(outPath,'resp','exportInfo');
disp('Save successful!');

% Announce finish
disp([mfilename ' completed!']);
% [End main function body]


end % [end main function block]

%----------------------BEGIN HELPER FUNCTIONS------------------------------
%% parse_stim_fname
function out = parse_stim_fname(str)

% NOTING HERE SO I DON'T FORGET
% Condensation = 1-3 = original
% Rarefaction = 4-6 = inverted

% Parse filename:
expr = ['^(?<inv>inv_)?s_(?<f1>\d{3})_(?<f2>\d{3})_(?<dur>\d{2,3})'...
    '(?:_\d+kHz)?\.wav'];
s = regexp(str,expr,'names');

% Stim type (always tone):
out.type = 'tone';
% Polarity:
if ~isempty(s.inv)
    out.polarity = 'rarefaction';
else
    out.polarity = 'condensation';
end
% Starting freq of glide:
out.startFreq = str2double(s.f1);
% Ending freq of glide:
out.endFreq = str2double(s.f2);
% Stimulus duration in ms:
out.durationMs = str2double(s.dur);
% Direction of glide:
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
% Slope of glide:
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
%---------------------------END HELPER FUNCTIONS---------------------------
% [EOF]