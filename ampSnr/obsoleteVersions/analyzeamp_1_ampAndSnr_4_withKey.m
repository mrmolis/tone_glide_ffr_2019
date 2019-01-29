%% SECTION 0: USER SETTINGS
%% 0.1) Folder Settings
defaultStimSource='C:\Users\vhaporMadseB\Desktop\Molis FFR stims\tone glides study\20 kHz\rarefaction';
defaultRespSource='C:\Users\vhaporMadseB\Desktop\Molis FFR datamats\tone glides study\Cz to C7'; % desired default filepath (0=current working directory)
defaultSavepath='\\R01PORHSM03\Research\NCRAR\Spectral_Dynamics_Grant\FFR\Brandon\Cz to C7 montage\Amplitude-based Analysis\';
versionSpec='_1_ampAndSnr_4_withKey';

%% 0.2) Analysis Settings
fInterest =[350,750]; % set frequency range of interest
fStep=1; %space between frequencies in Hz
defaultFs = 20000; %assumed sampling rate if not already specified in params
nWins=3; %number of desired windows
overlapRatio = 0;
prestimMs = 40;
stimMs=120;
peakSearchToleranceHz=25;

%% 0.3) Filter Settings
filterOrder=400;
passBandLow=300;
passBandHigh=800;

%% 0.4) Plot Settings
spaghetti=false; % do spaghetti plots?
groupAverage=true; % do group average plots?
groupOrder={'YN','NH','HI'};
jitterAmtX=0.2;
ynMark='^'; ynColor=[0.0, 0.0,  1.0]; % blue
nhMark='s'; nhColor=[0.0, 0.5,  0.0]; % green
hiMark='o'; hiColor=[1.0, 0.0,  0.0]; % red

% How to arrange/label trial types:
ttArrange={... % how trial types should be visually laid out
    '13_down',  '13_up' ;...
    '23_down',  '23_up' ;...
    '10_down',  '10_up' }; 
ttStrings={... % text labels parallel to ttArrange
    '1/3 octave down',  '1/3 octave up'     ;...
    '2/3 octave down',  '2/3 octave up'     ;...
    'full octave down', 'full octave up'    }; 

% Number of rows/columns of subplots:
nPlotRows = 3;
nPlotCols = 2;

% Gaps between subplots (normalized units):
horizGap = 0;
vertGap = 0;

% Page margins (normalized units):
marginBottom=0.175;
marginTop=0.01;
marginLeft=0.25;
marginRight=0.25;

%% SECTION 1: INITIALIZATION
%% 1.1) Clearing Variables and Closing Figures
clear ylim xlim s f t p sn fn tn pn
close all

%% 1.2) Creating Filter
fd=fdesign.bandpass('N,Fc1,Fc2',filterOrder,passBandLow,passBandHigh,...
    defaultFs);
hfd=design(fd,'FIR');

%% 1.3) Stimulus Folder
if ~exist('stimDir','var')
    stimDir=uigetdir(defaultStimSource,'Select folder with stim files');
end

%% 1.4) FFR Folder
if ~exist('fDir','var')
    fDir=uigetdir(defaultRespSource,'Select folder with response files');
end

% FFR data directory is starting folder for save dialogue unless otherwise
% specified in settings.
if ~exist('defaultSavepath','var')
    defaultSavepath=fDir;
elseif ~defaultSavepath %#ok
    defaultSavepath=fDir;
end

%% 1.5) Output Folder
% If the starting folder for this dialogue doesn't exist, create it
if ~exist(defaultSavepath, 'dir')
    mkdir(defaultSavepath);
end

% Skip the dialogue if an output directory is already stored as a variable
if ~exist('savepath','var')
    savepath=uigetdir(defaultSavepath,'Select where to save output files');
elseif ~savepath
    savepath=uigetdir(defaultSavepath,'Select where to save output files');
end

% Automatically add subfolder for script version to destination path in
% order to limit later confusion as to which version was used to produce
% the output
outpath=[savepath '\' versionSpec];

% If the version subfolder doesn't already exist, create it
if ~exist(outpath,'dir')
    mkdir(outpath)
end

%% SECTION 2: STIMULUS ANALYSIS
%% 2.1) Listing WAV Files
d=dir([stimDir '\*.wav']);

% Count number of stims
nStims = length(d);

% Put list as part of struct that will be saved later
settings.stimList=d;

%% 2.2) Loading WAV Files
ttStim=cell(nStims,1); % Pre-allocate space
for i=1:nStims
    fName = d(i).name;
    [x,stimFs(i)]=audioread([stimDir '\' fName]); %#ok<*SAGROW>
    stimx(:,i)=x(:,1); 
    
    %% 2.3) Trial Type
    switch fName(7:9)
        case '354'
            ttStim{i}='10_up'; 
        case '398'
            ttStim{i}='23_up'; 
        case '446'
            ttStim{i}='13_up'; 
        case '561'
            ttStim{i}='13_down'; 
        case '629'
            ttStim{i}='23_down'; 
        case '707'
            ttStim{i}='10_down'; 
    end % end switch
end; clear i % end loop thru stim files
clear d fName x

%% 2.4) DFT
stimPeakFreqs = stim_analyze(stimx,nWins,'dft');

%% SECTION 3: FFR ANALYSIS
%% 3.1) Listing Files
d=dir([fDir '\*120_C7_minus.mat']);
nFiles = length(d);
settings.respList=d;

%% 3.2) Loading Existing 'Resp.Mat' File (*)
% * if available from previous run - this is purely for speed
alreadyFiltered=false;
if exist([outpath '\resp.mat'],'file')
    tempIn=load([outpath '\resp.mat']);
    resp=tempIn.resp;
    clear tempIn
    alreadyFiltered=true;
end

%% 3.3) FFR Processing
% First pre-allocate space for variables
group=cell(nFiles,1); trialType=cell(nFiles,1); oCorr=zeros(nFiles,1);...
    oLagSamps=zeros(nFiles,1); oLagMs=zeros(nFiles,1);...
    subjId=cell(nFiles,1); slope=cell(nFiles,1); direction=cell(nFiles,1);

% Begin looping through files
for i=1:nFiles
    
    %% 3.3.1) FFR Processing: Parameters
    fName=d(i).name;
    sep=find(fName=='_');
    subjId{i}=fName(1:sep(1)-1);
    group{i}=fName(1:2);
    trialType{i}=fName(sep(2)+1:sep(4)-1);
    slope{i}=fName(sep(2)+1:sep(3)-1);
    direction{i}=fName(sep(3)+1:sep(4)-1);
    stim=stimx(:,strcmp(trialType{i},ttStim));
    
    %% 3.3.2) FFR Processing: Loading & Formatting
    % Use pre-filtered data if available
    if alreadyFiltered
        currResp=resp(i).wavAvg;
        load([fDir '\' fName],'params');
    % Otherwise filter now
    else
        load([fDir '\' fName]);
        currResp=mean(filtfilt(hfd.Numerator,1,...
            cast(epochArrayAccepted,'double')),2); % xcorr requires double
        clear epochArrayAccepted % conserve memory
        resp(i).wavAvg=currResp; % prep for saving to file later
    end
    
    % Identify sampling rate (Hz)
    if any(strcmp(fieldnames(params),'fs'))
        fs(i)=params.fs; %#ok<*SAGROW>
    elseif any(strcmp(fieldnames(params),'rate'))
        fs(i)=params.rate;
    else
        fs(i)=defaultFs;
    end
    
    % Derive additional variables
    % Some needed on every loop iteration...
    spms=fs(i)/1000; % samples per millisecond
    prestimSamps=prestimMs*spms;
    waitAfterZeroLag = 4 * spms; % earliest allowed lag for peak corr
    areaOfInterest = 8 * spms; % latest allowed lag for peak corr
    respTrimmed = currResp(prestimSamps+1:end);
    pad = zeros(length(respTrimmed)-length(stim),1);
    stimPadded = [stim;pad];
    % ... others only needed once...
    if i==1
        wMs=floor(stimMs/(nWins+overlapRatio-nWins*overlapRatio)); 
        overlapMs=wMs*overlapRatio;
    end
    
    %% 3.3.3) FFR Processing: Calculate Lag
    % Set RMS of response equal to RMS of stimulus
    stimPaddedRms=rms(stimPadded);
    respTrimmedRms=rms(respTrimmed);
    respTrimmedAmped=respTrimmed.*(stimPaddedRms/respTrimmedRms);
    
    % Find autocorrelation peak that will be used to standardize the
    % cross-correlation coefficients
    [stimCorr,lag]=xcorr(stimPadded,stimPadded,'none');
    [standardFactor,idx]=max(abs(stimCorr));
    if lag(idx)~=0
        error([fName ' does not have autocorr peak at lag zero.']);
    end
    
    % Calculate how far response lags behind stim
    [r,lag]=xcorr(respTrimmedAmped,stimPadded,'none');
    r=r./standardFactor;
    lagMs=lag./spms;
    temp = abs(r((length(currResp) + waitAfterZeroLag):...
        (length(currResp) + waitAfterZeroLag + areaOfInterest))); % there is an off-by-one error here in sample selection but I kept it to be consistent with Sam
    [oCorr(i),idx]=max(temp); clear temp% off by one
    oLagSamps(i)=idx+waitAfterZeroLag; % off by one
    oLagMs(i)=lagMs(idx+length(currResp)+waitAfterZeroLag/spms); % off by one
    clear lag r lagMs idx
    
    %% 3.3.4) FFR Processing: DFT
    wSamps=wMs*spms;
    w=blackman(wSamps); % w = envelope of the window at each sample point
    o=floor(overlapRatio*wSamps); % o = overlap in number of samples
    stimSamps = stimMs*spms;
    f=fInterest(1):fInterest(2);
    tStimOn=prestimSamps+1;
    tCaptureOn=tStimOn+oLagSamps(i);
    tCaptureOff=tCaptureOn+stimSamps-1;
    [dft(:,:,i),f,t,~] =... % during stim
        spectrogram(currResp(tCaptureOn:tCaptureOff),w,o,f,fs(i)); 
    [dftPre(:,1,i),f,tn,~]=... % prestim interval
        spectrogram(currResp(tStimOn-wSamps:tStimOn-1),w,0,f,fs(i)); 
    tPre=tn-prestimMs/1000; 
end; clear i % end loop thru response files
clear tn

%% 3.4) Saving New 'Resp.Mat' File (*)
% * unless current data was already loaded from a saved Matlab file
if ~alreadyFiltered
    save([outpath '\resp.mat'],'resp');
end

%% 3.5) FFR Spectral Peaks
nResps = size(dft,3);

if exist([outpath '\ffrAmpAnalysisOut.mat'],'file')
    load([outpath '\ffrAmpAnalysisOut.mat']);
else
    tol=peakSearchToleranceHz;
    for i=1:nResps
        ttStimLogic=strcmp(trialType{i},ttStim);
        for j=1:nWins
            
            currResp=abs(squeeze(dft(:,j,i)));
            currNoise=abs(squeeze(dftPre(:,1,i)));
            stimPeak=stimPeakFreqs(ttStimLogic,j);
            
            % Generate correct frequency vector
            fStart = stimPeak - tol;
            fStop = stimPeak + tol;
            fLogic = ((f >= fStart) & (f <= fStop));
            fTemp=f(fLogic);
            clear fStart fStop
            
            % Grab corresponding part of response and noise
            tempResp=currResp(fLogic);
            tempNoise=currNoise(fLogic);
            clear currNoise currResp fLogic
            
            % Record amplitude and frequency of signal peak
            [val,ind]=max(tempResp);
            peakVal(i,j)=val;
            peakFreq(i,j)=fTemp(ind);
            clear val ind tempResp
            
            % Estimate average noise in the analysis band
            avgNoise(i,j)=mean(tempNoise);
            clear tempNoise
            
        end; clear j % end loop thru windows
    end; clear i % end loop through responses
    clear ttStimLogic tol
end % end if-else block

%% 3.6) FFR SNR
snr=10*log10(peakVal./avgNoise); % using average noise

%% 3.7) Calculating Group Averages
groupList = unique(group);          nGroups = numel(groupList);
slopeList = unique(slope);          nSlopes = numel(slopeList);
directionList = unique(direction);  nDirections = numel(directionList);
ttList=unique(trialType);           ntt = numel(ttList);

count=1;
ga=struct;
for i=1:nGroups
    currGroup=groupList{i};
    groupLogic=strcmp(group,currGroup);
    for j=1:nSlopes
        currSlope=slopeList{j};
        slopeLogic=strcmp(slope,currSlope);
        for k=1:nDirections
            currDirection=directionList{k};
            directionLogic=strcmp(direction,currDirection);
            ga(count).direction=currDirection;
            comboLogic=groupLogic&slopeLogic&directionLogic;
            ga(count).group=currGroup;
            ga(count).slope=currSlope;
            ga(count).avgNoise=mean(avgNoise(comboLogic,:));
            ga(count).peakFreq=mean(peakFreq(comboLogic,:));
            ga(count).peakVal=mean(peakVal(comboLogic,:));
            ga(count).snr=mean(snr(comboLogic,:));
            count=count+1;
        end; clear k
    end; clear j
end; clear i
clear count

%% 3.8) Saving Amplitude & SNR Analysis Output
% As MAT file
save([outpath '\ffrAmpAnalysisOut.mat'],...
    'avgNoise', 'peakFreq', 'peakVal', 'snr','ga', 'direction','slope',...
    'group','subjId','trialType');

% As XLSX files - Raw
% Individual - Raw
headerRowIndiv={'SubjectId','Group','Slope','Direction','Window1',...
    'Window2','Window3'};
rawIndiv=cell(nResps+1,numel(headerRowIndiv));
rawIndiv(1,:)=headerRowIndiv;
rawIndiv(2:nResps+1,1)=subjId(:);
rawIndiv(2:nResps+1,2)=group(:);
rawIndiv(2:nResps+1,3)=slope(:);
rawIndiv(2:nResps+1,4)=direction(:);
rawIndiv(2:nResps+1,5:7)=num2cell(peakVal(:,:));
xlswrite([outpath '\rawOut.xlsx'],rawIndiv,1,'A1');

% Group Average - Raw
headerRowGroup={'Group','Slope','Direction','Window1','Window2','Window3'};
rawGroup=cell(ntt*nGroups,numel(headerRowGroup));
rawGroup(1,:)=headerRowGroup;
for i=1:numel(ga)
    rawGroup{i+1,1}=ga(i).group;
    rawGroup{i+1,2}=ga(i).slope;
    rawGroup{i+1,3}=ga(i).direction;
    rawGroup(i+1,4:6)=num2cell(ga(i).peakVal);
end; clear i
xlswrite([outpath '\rawOut.xlsx'],rawGroup,2,'A1');

% Label worksheet tabs properly
e=actxserver('Excel.Application');
ewb=e.Workbooks.Open([outpath '\rawOut.xlsx']);
ewb.Worksheets.Item(1).Name='Indiv';
ewb.Worksheets.Item(2).Name='Avg';
ewb.Save;
ewb.Close(false);
e.Quit;

% Individual - SNR
snrIndiv=cell(nResps+1,numel(headerRowIndiv));
snrIndiv(1,:)=headerRowIndiv;
snrIndiv(2:nResps+1,1)=subjId(:);
snrIndiv(2:nResps+1,2)=group(:);
snrIndiv(2:nResps+1,3)=slope(:);
snrIndiv(2:nResps+1,4)=direction(:);
snrIndiv(2:nResps+1,5:7)=num2cell(snr(:,:));
xlswrite([outpath '\snrOut.xlsx'],snrIndiv,1,'A1');

% Group Average - SNR
snrGroup=cell(ntt*nGroups,numel(headerRowGroup));
snrGroup(1,:)=headerRowGroup;
for i=1:numel(ga)
    snrGroup{i+1,1}=ga(i).group;
    snrGroup{i+1,2}=ga(i).slope;
    snrGroup{i+1,3}=ga(i).direction;
    snrGroup(i+1,4:6)=num2cell(ga(i).snr);
end; clear i
xlswrite([outpath '\snrOut.xlsx'],snrGroup,2,'A1');

% Label worksheet tabs properly
e=actxserver('Excel.Application');
ewb=e.Workbooks.Open([outpath '\snrOut.xlsx']);
ewb.Worksheets.Item(1).Name='Indiv';
ewb.Worksheets.Item(2).Name='Avg';
ewb.Save;
ewb.Close(false);
e.Quit;

%% SECTION 4: PLOTS
m=cell(size(group)); c=cell(size(group)); 
groupMarkers={ynMark,nhMark,hiMark}; 
groupColors={ynColor,nhColor,hiColor};
for i=1:nGroups
    currGroup=groupOrder{i};
    groupLogic=strcmp(currGroup,group);
    m(groupLogic)=groupMarkers(i);
    c(groupLogic)=groupColors(i);
end; clear i % end loop through groups
clear currGroup groupLogic

quants={'Raw','SNR'};

%% 4.1) Spaghetti Plots
if spaghetti
    close all %#ok<UNRCH>
    hLin=cell(ntt,3);
    for i=1:ntt
        currType=ttList{i};
        ttLogic=strcmp(currType,trialType);
        mTemp=m(ttLogic);
        cTemp=c(ttLogic);
        for j=1:3
            switch j
                case 1
                    y=transpose(peakVal(ttLogic,:));
                    quantType='Raw';
                case 2
                    y=transpose(snr(ttLogic,:));
                    quantType='SNR';
            end
            x=transpose(repmat(1:nWins,sum(ttLogic),1));
            jitterArrayX=jitterAmtX.*rand(size(y))-0.5*jitterAmtX;
            xWithJitter = x + jitterArrayX;
            hFig=figure('Units','inches',...
                'OuterPosition',[0,0,6.5,10],... % [x1,y1,x2,y2]
                'Resize','off');
            hLin = line(xWithJitter,y);
            set(hLin,...
                {'Marker'},                     mTemp,...
                {'MarkerFaceColor'},            cTemp,...
                {'Color'},                      cTemp,...
                {'MarkerEdgeColor'},            {'k'}   );
            outname=['spag_' currType '_' currMeas '_' quantType];
            title(outname,'Interpreter','None');
            print(gcf,'-dpng','-r300',[outpath '\' outname '.png']);
            close(gcf);
        end; clear j % end loop thru absolute and SNR
        clear x xWithJitter y %jitterArrayX

        clear currMeas
        
    end; clear i % end loop through trial types
    clear cTemp currType mTemp ttLogic hLin
end % end if block (spaghetti)

%% 4.2) Group Average Plots
ttOrder=transpose(ttArrange);
ttLabels=transpose(ttStrings);
if groupAverage
    
    yUnits = 'Raw DFT Value';
    yMax=35;
    yMin=0;
    yTicks=5:5:30;
    
    for j=1:2
        switch j
            case 1
                temp=peakVal;
                quantType='Raw';
            case 2
                temp=snr;
                quantType='SNR';
                yUnits = 'SNR (dB)';
                yMax=12;
                yMin=0;
                yTicks=2:2:10;
        end % end switch
        
        hFig=figure(...
            'Units',            'inches',...
            'OuterPosition',    [0,0,6.5,10],... % [x1,y1,x2,y2]
            'Resize',           'off');
        hPlots=tight_subplot(...
            nPlotRows,      nPlotCols       ,...
            [   vertGap         horizGap    ]   ,...
            [   marginBottom,   marginTop   ]   ,...
            [   marginLeft,     marginRight ]   );
        
        for i=1:ntt
            currType = ttOrder{i};
            ttLogic = strcmp(currType,trialType);
            nSubs = sum(ttLogic);
            axes(hPlots(i)); %#ok<*LAXES>
            
            for g=1:nGroups
                currGroup = groupOrder{g};
                groupLogic = strcmp(group,currGroup);
                tempPared=temp(ttLogic&groupLogic,:);
                y = mean(tempPared,1);
                yError = std(tempPared,1)./sqrt(nSubs);
                x = 1:nWins;
                jitterX = g * 0.1 - 0.2;
                xWithJitter = x + jitterX;
                hLin(g) = line(...
                    'XData',            xWithJitter     ,...
                    'YData',            y               ,...
                    'Marker',           groupMarkers{g} ,...
                    'Color',            groupColors{g}  ,...
                    'MarkerFaceColor',  groupColors{g}  ,...
                    'MarkerEdgeColor',  'k'             );
                hold on
                hErr(g) = errorbar(xWithJitter,y,yError);
                set(hErr(g),'Color',groupColors{g});
            end; clear g
            
            xlim([0.5 3.5]);
            
            if all(ismember({'yMax','yMin'},who))
                ylim([yMin yMax]);
            end
            if exist('yTicks','var')
                set(gca,'YTick',yTicks);
            end
            
            title(ttLabels{i},... % subplot box titles
                'Interpreter',          'None'          ,...
                'Units',                'Normalized'    ,...
                'Position',             [0.5,0.85]      ,...
                'HorizontalAlignment',  'Center'        ,...
                'FontWeight',           'Normal'        ,...
                'FontSize',             11              );
            
            if mod(i,2)
                set(gca,'YTickLabel',get(gca,'YTick'));
            end
            
            if i==3
                ylabel(yUnits);
            end
            
            set(gca,'XTick',1:3);
            if ismember(i,[5,6])
                set(gca,'XTickLabel',{'0-40','40-80','80-120'});
            end
            
            grid off
            box on
            
        end; clear i % end loop thru trial types
        
        text(... % x-axis title
            'String',               'Time re: Stimulus Onset (ms)',...
            'Units',                'Normalized'                  ,...
            'Position',             [0,-0.17]                     ,...
            'HorizontalAlignment',  'Center'                      ,...
            'FontWeight',           'Normal'                      );
        
        clear yMin yMax yTicks
        outname=['groupAvg_DFT_' quantType];
        
        legend(hLin,groupOrder(:),...
            'Units',        'Normalized'                ,...
            'Location',     ([.19, .03, .62, 0.05])      ,... %[x y w h]
            'Orientation',  'Horizontal'                );
        
        print(gcf,'-dpng','-r300',[outpath '\' outname '.png']);
        close(gcf);
        
    end; clear j % end loop thru absolute and SNR
end % end if block (groupAverage)