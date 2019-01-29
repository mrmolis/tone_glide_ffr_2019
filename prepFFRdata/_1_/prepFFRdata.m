function savePath = prepFFRdata(respPath,lookupDir,savePath)
% Go through an existing file of compiled FFR data and make sure
% everything's there that needs to be for subsequent analyses.  If it's
% not, look for the original .mat files in the lookup directory and fill in
% missing info using those.  ONLY DESIGNED FOR FFR TONE GLIDE STUDY
ei.functionSpec = mfilename;
ei.versionSpec = '_1_';

defRespDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
    'matFiles_newProcess2017All\summaryMats'];
defLookDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
    'matFiles_newProcess2017All\indivMats'];

if ~exist('respPath','var')
    [respName,respDir] = uigetfile(defRespDir,'Load response file');
    respPath = [respDir '\' respName];
    if isequal(respName,0)
        warning('No response file specified; function will abort');
        return
    end
end
defSavePath = respPath;

if ~exist('lookupDir','var')
    lookupDir = uigetdir(defLookDir,'Point to lookup directory');
    if isequal(lookupDir,0)
        warning('No lookup directory specified; function will abort');
    end
end

if ~exist('savePath','var')
    disp('Opening output-file dialog box...');
    [saveFile,saveDir] = uiputfile(defSavePath,'Save output');
    disp('Saving...')
    savePath = [saveDir '\' saveFile];
end

disp('Attempting to load response data...');
if isstruct(respPath)
    resp = respPath;
    clear respPath
else
    load(respPath);
    clear respFile
end
disp('Success!');

% Fields to check for
reqFlds = {'wavAvg','filterInfo'};

% Get the lookup list ready
disp('Generating lookup table...');
subfolds = get_all_subfolders(lookupDir);
fList = {};
pathList = {};
for i=1:numel(subfolds)
    fNames = shiftdim(getfiles(subfolds{i}));
    fPaths = shiftdim(strcat(subfolds{i},'\',fNames));
    fList = [fList; fNames]; %#ok<*AGROW>
    pathList = [pathList; fPaths];
end; clear i
for i=1:numel(resp)
    r = resp(i);
    if isempty(regexp(r.fName,'bc_ar','once'))
        fn = [resp(i).fName(1:end-4) '_bc_ar.mat'];
        resp(i).fName = fn;
    else
        fn = r.fName;
    end
    logic = strcmp(fn,fList);
    if sum(logic) > 1
        error('too many matches');
    elseif sum(logic) < 1
        error('no matches');
    end
    lookupPath(i) = pathList(logic); 
end

% turn non-existent fields into empty ones and pre-find lookup paths
for j=1:numel(reqFlds)
    fld = reqFlds{j};
    if ~isfield(resp,fld)
        resp(1).(fld) = []; % double check this works
    end
end

% wavAvg & filterInfo
for i=1:numel(resp)
    prog = [num2str(i) '/' num2str(numel(resp))];
    disp (['Checking ' prog '...']);
    r = resp(i);
    if isempty(r.wavAvg)
        np = r.epochInfo.nPoints;
        lp = lookupPath{i};
        resp(i).fPath = lp;
        temp = load(lp);
        epochs = cast(temp.epochs,'double'); % filtfilt needs double
        pInd = find(np==size(epochs));
        switch pInd
            case 1
                % no alteration needed
            case 2
                epochs = transpose(epochs);
            otherwise
                error('invalid dimension index (pInd)');
        end
        filterOrder = 400;
        passBandLow = 300;
        passBandHigh = 800;
        fd = fdesign.bandpass('N,Fc1,Fc2',filterOrder,passBandLow,...
            passBandHigh,r.fs);
        hfd = design(fd,'FIR'); clear fd
        resp(i).wavAvg = mean(filtfilt(hfd.Numerator,1,epochs),2);
    end
    if isempty(r.filterInfo)
        filterInfo.order = 400;
        filterInfo.passBandLow = 300;
        filterInfo.passBandHigh = 800;
        filterInfo.type = 'FIR';
        filterInfo.shape = 'bandpass';
        resp(i).filterInfo = filterInfo;
    end
end

% Order fields consistently
disp('Sorting fields...')
resp = sortfields(resp); %#ok<NASGU>

if exist('exportInfo','var')
    ind = numel(exportInfo) + 1; %#ok<NODEF>
    exportInfo(ind).functionSpec = ei.functionSpec;
    exportInfo(ind).versionSpec = ei.versionSpec;
    exportInfo(ind).timestamp = datetime;
else
    exportInfo = ei; 
    exportInfo.timestamp = datetime;
end

% Save
disp('Saving output...');
if ~exist('saveDir','var')
    if exist('savePath','var')
        bs = find(savePath=='\');
        saveDir = savePath(1:bs(end));
    else
        error('no save dir or path');
    end
end
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end
save(savePath,'resp','exportInfo');
disp('Saved!');
disp([mfilename ' completed!']);
end % main function

