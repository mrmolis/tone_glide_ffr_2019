function bc_artRej_plusMinus(inDir, rejTH, bufferMs)
% Build 5.
%
% CHANGELOG...
% BUILD 5:
% - More flexible parsing of filenames to allow for the flat-tone naming
% conventions.
% BUILD 4:
% - Uses 30 uV reject threshold instead of 10% reject rate

% import functions
addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\Matlab_Code'...
    '\artrej_by_thresh\_1_']);

close all

outDir = [inDir '\bc_ar_pm'];

% List mat files
d = dir([inDir '\*.mat']);
fList = shiftdim({d.name}); clear d
% Grab just the ones that are bins 1-3
expr = '13(?=((_\w{2}ref)?(_\w{2,3})?(_eeg)?\.mat))';
logic13 = ~cellfun('isempty',regexpc(fList,expr,'once'));
fList = fList(logic13);
% Also list full paths
pathList = strcat(inDir, '\', fList);

for i=1:numel(fList)
    fName = fList{i};
    disp([num2str(i) '/' num2str(numel(fList)) ': ' fName]);
    fPath = pathList{i};
    cond = load(fPath);
    [ind1,ind2] = regexp(fPath,expr,'once');
    rarePath = [fPath(1:ind1-1) '46' fPath(ind2+1:end)];
    bs = find(rarePath=='\');
    rareName = rarePath(bs(end)+1:end);
    if exist(rarePath,'file')
        rare = load(rarePath);
    else
        error([fName ': No matching "46" file found']);
    end
    
    % Make sure nPoints is equal
    np1 = cond.params.epochInfo.nPoints;
    np2 = rare.params.epochInfo.nPoints;
    msg = [fName ': ' np1 ' vs ' np2 ' nPoints mismatch'];
    assert(np1==np2,msg);
    minPoints = min([np1 np2]);
    
    % Make sure array is oriented properly for feeding into next part
    sz1 = size(cond.epochs);
    pLogic = sz1==np1;
    switch find(pLogic)
        case 1
            cond.epochs = transpose(cond.epochs);
        case 2
            % do nothing
        otherwise
            error(['pLogic error ' fName]);
    end
    sz2 = size(rare.epochs);
    pLogic = sz2==np2;
    switch find(pLogic)
        case 1
            rare.epochs = transpose(rare.epochs);
        case 2
            % do nothing
        otherwise
            error(['pLogic error ' fName(1:end-6) '46.mat']);
    end
    
    %% 3.1) Baseline correction & artifact rejection
    for j=1:2
        switch j
            case 1
                temp = cond;
            case 2
                temp = rare;
        end
        
        % in case it's from a batch where this wasn't added earlier...
        psms = temp.params.epochInfo.prestimMs;
        prestimSamps = psms * temp.params.spms;
        temp.params.epochInfo.prestimSamps = prestimSamps;
        stimSamps = temp.params.stim.durationMs * temp.params.spms;
        bufferSamps = temp.params.spms * bufferMs;
        
        % baseline correction (average of prestim baseline)
        disp('Performing baseline correction...');
        for k=1:size(temp.epochs,1)
%             disp(['BC: sweep ' num2str(k) ' / '...
%                 num2str(size(temp.epochs,1))]);
            psi = temp.epochs(k,1:prestimSamps);
            baseline = mean(psi);
            temp.epochs(k,:) = temp.epochs(k,:) - baseline;
        end
        
        % artifact rejection (ratio)
        disp('Performing artifact rejection...');
        [epochsAccepted,nAcc,nRej,rejThreshold] = artrej_by_thresh(...
            temp.epochs,rejTH,prestimSamps, stimSamps, bufferSamps);
        temp.epochs = epochsAccepted;
        temp.params.epochInfo.nAcceptedMatlab = nAcc;
        temp.params.epochInfo.nRejectedMatlab = nRej;
        if ~isempty(rejThreshold)
            temp.params.epochInfo.artRejThresholdAbsUv = rejThreshold;
        end
        switch j
            case 1
                cond = temp;
            case 2
                rare = temp;
        end
    end
    clear temp
    for j=1:2
        switch j
            case 1 
                temp = cond;
            case 2
                temp = rare;
        end
        flds = fieldnames(temp.params);
        for k=1:numel(flds)
            fld = flds{k};
            if isstruct(temp.params.(fld))
                temp.params.(fld) = orderfields(temp.params.(fld));
            end
        end
    end
    cond.params.stim.type = 'tone';
    rare.params.stim.type = 'tone';
    cond.params = orderfields(cond.params);
    rare.params = orderfields(rare.params);
    outName = [fName(1:end-4) '_bc_ar.mat'];
    outPath = [outDir '\' outName];
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    save(outPath,'-struct','cond');
    outName = [rareName(1:end-4) '_bc_ar.mat'];
    outPath = [outDir '\' outName];
    save(outPath,'-struct','rare');
    clear outname outpath
    
    % Take lowest number of sweeps
    sw1 = cond.params.epochInfo.nAcceptedMatlab;
    sw2 = rare.params.epochInfo.nAcceptedMatlab;
    minSweeps = min([sw1 sw2]);
    c1 = cond.params.epochInfo.count;
    c2 = rare.params.epochInfo.count;
    minCount = min([c1 c2]);
    
    %% 3.3) Plus and minus
    disp('Making plus and minus versions...');
    temp13 = cond.epochs(1:minSweeps,1:minPoints);
    temp46 = rare.epochs(1:minSweeps,1:minPoints);
    minus.epochs = (temp13 - temp46)./2; 
    plus.epochs = (temp13 + temp46)./2;  
    clear temp*
    
    params = mergeparams(cond.params,rare.params);
    params.epochInfo.nAcceptedMatlab = minSweeps;
    params.epochInfo.nRejectedMatlab = minCount - minSweeps;
    params.stim.type = 'tone';
    flds = fieldnames(params);
    for k=1:numel(flds)
        fld = flds{k};
        if isstruct(params.(fld))
            params.(fld) = orderfields(params.(fld));
        end
    end
    params = orderfields(params);
    plus.params = params;
    plus.params.stim.polarity = 'plus';
    minus.params = params;
    minus.params.stim.polarity = 'minus';
    
    % Output   
    disp('Saving...');
    [ind1,ind2] = regexp(fName,expr,'once');
    outName = [fName(1:ind1-1) 'minus' fName(ind2+1:end-4) '_bc_ar.mat'];
    outPath = [outDir '\' outName];
    save(outPath,'-struct','minus');
    outName = [fName(1:ind1-1) 'plus' fName(ind2+1:end-4) '_bc_ar.mat'];
    outPath = [outDir '\' outName];
    save(outPath,'-struct','plus');
    disp('Success!');
    
end
disp([mfilename ' completed!']);

% end main function body, put nested functions below, before true "end" of
% function

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
            needsResolve(i) = false; %#ok<*AGROW>
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
        % Opt to split at lowest level possible. Function calling itself is
        % an elegant solution to the fact that number of nested levels is
        % not known a priori.
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