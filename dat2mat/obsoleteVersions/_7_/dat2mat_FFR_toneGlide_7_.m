function dat2mat_FFR_toneGlide_7_

% import functions
addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\Matlab_Code\'...
    'artrej_by_ratio\_2_lessVerbose']);
addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\Matlab_Code\'...
    'bc_artRej_plusMinus\_3_fixedPol']);

clear variables
close all

%% SECTION 0: OPTIONS & SETTINGS
%% 0.1) Directories
inDir = 'C:\Users\vhapormadseb\Desktop\Molis FFR dat files\tone glide study\A1toA2';
outDir = [inDir '\mats'];

%% 0.2) Artifact Rejection
rejRAT = 0.1; % target rejection rate as proportion
bufferMs = 20; % ms after end of stim to include in AR range

%% SECTION 1: INITIALIZATION & SETUP
disp('SECTION 1: INITIALIZATION & SETUP');
%% 1.1) Declarations
fList=[];subject=[];params=[]; isArtifactRun=[]; slope=[]; direction=[]; %#ok<*NASGU>
epochData13=single([]);
epochData46=single([]);
refElec = []; ear=[]; durationMs=[];
testDate=[]; testTime=[]; fileType=[]; hasSweepHeaders=[]; xUnits=[];
yUnits=[]; varList=[]; prestimMs=[]; epochStartMs=[];

%% 1.2) Pick directory to process if unspecified in settings
if ~exist('inDir','var')
    inDir = uigetdir(pwd,'Choose folder to analyze.');
end

%% 1.3) Compile files
% Compile only ".dat" files
d = dir([inDir '\*.dat']);
fList = shiftdim({d.name}); clear d

% List full paths
pathList = strcat(inDir,'\',fList);

%% 1.4) Make output dir if necessary
if ~exist(outDir,'dir')
    disp('Making output folder...');
    mkdir(outDir);
    disp('Output folder created!');
end

%% SECTION 2: PROCESS FILES
disp('SECTION 2: PROCESS FILES');
for i=1:numel(fList)
    % Grab current file
    fName = fList{i};
    fPath = pathList{i};
    outpath = [outDir '\' fName(1:end-4) '.mat'];
    disp([num2str(i) '/' num2str(numel(fList)) ': ' fName]);
    if exist(outpath,'file')
        continue
    end
    
    % Parse file name
    params = parse_fname(fName);
    
    % Open file for reading
    fid = fopen([inDir '\' fName]); % use full orig filename so it loads
    str = '';

    % Get info from file headers
    % Subject
    while isempty(regexp(str,'\[Subject\]','once'))
        str = fgetl(fid);
    end
    % subject is usually blank so will not store it from here - use
    % filename instead to derive subject
    
    % Date
    while isempty(regexp(str,'\[Date\]','once'))
        str = fgetl(fid);
    end
    params.sessionInfo.testDate = regexp(str,'\d{2}/\d{2}/\d{4}',...
        'match','once');
    
    % Time
    while isempty(regexp(str,'\[Time\]','once'))
        str = fgetl(fid);
    end
    params.sessionInfo.testTime = regexp(str,'\d{2}:\d{2}:\d{2}',...
        'match','once');
    
    % Number of channels in dat file
    while isempty(regexp(str,'\[Channels\]','once'))
        str = fgetl(fid);
    end
    nChans = regexp(str,'\d+','match','once');
    nChans = str2double(nChans);
    params.sourceDatInfo.nChannelsInDat = nChans;
    
    % Sampling rate
    while isempty(regexp(str,'\[Rate\]','once'))
        str = fgetl(fid);
    end
    expr = '\d+\.?\d*';
    params.fs = str2double(regexp(str,expr,'match','once'));
    params.spms = params.fs/1000;
    
    % Neuroscan file type
    while isempty(regexp(str,'\[Type\]','once'))
        str = fgetl(fid);
    end
    params.sourceDatInfo.originFileType = regexp(str(7:end),'\S*',...
        'match','once');
    
    % Has sweep headers?
    while isempty(regexp(str,'\[Headers\]','once'))
        str = fgetl(fid);
    end
    temp = regexp(str(10:end),'\S*','match','once');
    switch temp
        case 'No'
             bool = false;
        case 'Yes'
            bool = true;
    end
    params.sourceDatInfo.hasSweepHeaders = bool;
    
    % Number of sample points
    while isempty(regexp(str,'\[Points\]','once'))
        str = fgetl(fid);
    end
    nPoints = str2double(regexp(str,'\d+',...
        'match','once'));
    params.epochInfo.nPoints = nPoints;
    
    while isempty(regexp(str,'\[Xmin\]','once'))
        str = fgetl(fid);
    end
    expr = '-?\d+\.?\d*';
    xMin = str2double(regexp(str,expr,'match','once'));
    psms = -1* xMin;
    params.epochInfo.prestimMs = psms;
    params.epochInfo.xMin = xMin;
    params.epochInfo.prestimSamps = psms * params.spms;
    
    while isempty(regexp(str,'\[Sweeps\]','once'))
        str = fgetl(fid);
    end
    expr = '\d+';
    params.epochInfo.count = str2double(regexp(str,expr,'match','once'));
    
    while isempty(regexp(str,'\[Accepted\]','once'))
        str = fgetl(fid);
    end
    expr = '\d+';
    nAcc = str2double(regexp(str,expr,'match','once'));
    params.epochInfo.nAcceptedNS = nAcc;
    
    while isempty(regexp(str,'\[Rejected\]','once'))
        str = fgetl(fid);
    end
    expr = '\d+';
    params.epochInfo.nRejectedNS = str2double(regexp(str,expr,...
        'match','once'));
    
    while isempty(regexp(str,'\[Electrode Labels\]','once'))
        str = fgetl(fid);
    end
    str = fgetl(fid);
    expr = '\[\s*(\w+)\s*\]';
    eNames = regexp(str,expr,'tokens');
    for j=1:numel(eNames)
        eNames{j} = eNames{j}{1}; %#ok<*AGROW>
    end
    params.sourceDatInfo.channelNames = eNames(:);
    
    while isempty(regexp(str,'\[Electrode XUnits\]','once'))
        str = fgetl(fid);
    end
    str = fgetl(fid);
    params.sourceDatInfo.xUnits = regexp(str,'[^ \s\]\[]*','match','once');
    
    while isempty(regexp(str,'\[Electrode YUnits\]','once'))
        str = fgetl(fid);
    end
    str = fgetl(fid);
    params.sourceDatInfo.yUnits = regexp(str,'[^ \s\]\[]*','match','once');
    
    % If multiple electrodes, only analyze the first
    fmt='%f32';
    for j=2:nChans
        fmt=[fmt '%*f32']; % * = ignore
    end; clear j
    
    % Loop through nSweeps
    for j=1:nAcc
        
        while isempty(regexp(str,'\[Epoch Data\]','once'))
            str = fgetl(fid);
        end
        str = '';
        epochs(j,:) = single(cell2mat(textscan(fid,fmt,nPoints,'CollectOutput',1)));
        
    end
    
    fclose(fid);
    
    exportInfo.timestamp = datetime;
    exportInfo.versionSpec = mfilename;
    exportInfo = orderfields(exportInfo);
    params = orderfields(params);
    
    save(outpath,'epochs','params','exportInfo');
    clear epochs params exportInfo
    
end

%% SECTION 3: BASELINE CORRECTION, ARTIFACT REJECTION & PLUS/MINUS VERSIONS
disp('SECTION 3: ARTIFACT REJECTION & PLUS/MINUS VERSIONS');

bc_artRej_plusMinus(outDir, rejRAT, bufferMs);

disp([mfilename ' completed!']);

% end main function body, put nested functions below, before "end"

%% NESTED FUNCTIONS
%---------------------------BEGIN NESTED FUNCTIONS-------------------------
% These have access to variable values from within the main function even
% if they are not passed explicitly.
%--------------------------------------------------------------------------

%----------------------------END NESTED FUNCTIONS--------------------------
end % end main function block, put helper functions below



%% HELPER FUNCTIONS
%----------------------------BEGIN HELPER FUNCTIONS------------------------
% These only have access to variables/values from main function if they are
% explicitly passed.
%--------------------------------------------------------------------------
function out = parse_fname(str)
% example format: HI007_R_10_down_120_rerefA1_13_A2

expr = ['^(?<sid>\w{2}\d{3})'...
    '_(?<ear>[LR])'...
    '_(?<slope>10|13|23)'...
    '_(?<dir>[a-z]{2,4})'...
    '_(?<dur>\d{2,3})'...
    '_?(?<art>art)?'...
    '(?:_reref)?(?<ref>\w{2})?'...
    '_(?<bins>[1346]{2})'...
    '_(?<targ>\w{2})\.dat'];
toks = regexp(str,expr,'names');
% subject ID
out.subject.id = toks.sid;
% subject group
out.subject.group = regexp(toks.sid,'\w{2}','match','once');
% subject number
sn = regexp(toks.sid,'\d{3}','match','once');
out.subject.serialNumberAsString = sn;
out.subject.serialNumberAsDouble = str2double(sn);
% presentation ear
out.stim.laterality = toks.ear;
% glide slope
out.stim.slope = toks.slope;
switch toks.slope
    case '10'
        sl = 'Full octave';
    case '23'
        sl = '2/3 octave';
    case '13'
        sl = '1/3 octave';
    otherwise
        error([str ': unrecognized slope']);
end
out.stim.slopeLong = sl;
% glide direction
out.stim.direction = toks.dir;
% stim duration in ms
out.stim.durationMs = str2double(toks.dur);
% artifact run?
if isempty(toks.art)
    out.isArtifact = false;
else
    out.isArtifact = true;
end
% reference electrode location
if isempty(toks.ref)
    out.montage.reference = 'A2';
else
    out.montage.reference = toks.ref;
end
% bins/polarity
out.stim.bins = toks.bins;
switch toks.bins
    case '13'
        pol = 'condensation';
    case '46'
        pol = 'rarefaction';
    case {'plus','minus'}
        pol = toks.bins;
    otherwise
        error([str ': unrecognized bins']);
end
out.stim.polarity = pol;
% electrode to analyze
out.montage.channelToAnalyze = toks.targ;
% other
out.stim.type = 'tone';
end

%--------------------------------------------------------------------------
%% mergeparams
function p = mergeparams(p1,p2)
flds1 = fieldnames(p1);
flds2 = fieldnames(p2);

% Copy over the uncontested ones "as is"
for i=1:numel(flds1)
    fld = flds1{i};
    nameLogic = strcmp(fld,flds2);
    switch sum(nameLogic)
        case 0
            p.(fld) = p1.(fld);
            needsResolve(i) = false;
        case 1
            if isequal(p1.(fld),p2.(fld))
                p.(fld) = p1.(fld);
                needsResolve(i) = false;
            else
                needsResolve(i) = true;
            end
        otherwise
            error('more than one field match');
    end
            
end
for i=1:numel(flds2)
    fld = flds2{i};
    nameLogic = strcmp(fld,flds1);
    switch sum(nameLogic)
        case 0
            p.(fld) = p2.(fld);
        case 1
            % will already have been caught by previous loop
        otherwise
            error('more than one field match');
    end
end

% Resolve conflicts by making two parallel fields
inds = find(needsResolve);
for i=1:numel(inds)
    ind = inds(i);
    fld = flds1{ind};
    if isstruct(p1.(fld)) && isstruct(p2.(fld))
        % Opt to split at lowest level possible. Will have to see if this
        % is kosher to have function call itself. Otherwise will need to
        % make identical function to have it call.
        p.(fld) = mergeparams(p1.(fld),p2.(fld));
    else
        fld1 = [fld '13'];
        fld2 = [fld '46'];
        p.(fld1) = p1.(fld);
        p.(fld2) = p2.(fld);
    end
end

end % end function mergeparams
%--------------------------------------------------------------------------

%-----------------------------END HELPER FUNCTIONS-------------------------
% [EOF]