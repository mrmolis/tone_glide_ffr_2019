% STIM_ANALYZE - updated 7/3/2017

function fPeaks=stim_analyze(stimx,nWins,choice)
% Chunk stim into NWINS equal time windows, find peak frequency in each
% window and return as FPEAKS.

%% Settings
fs = 2e4;
fLow=350;
fHigh=750;
fStep=1;
defaultChoice='dft';
defaultWins=3;

%% Rest of code

% Make frequency vector from parameters in 'Settings' section
fVec=fLow:fStep:fHigh;

% Handle missing second argument
if ~exist('nWins','var')
    nWins=defaultWins;
elseif isempty(nWins)
    nWins=defaultWins;
end

% Handle missing third argument
if ~exist('choice','var')
    choice=defaultChoice;
elseif isempty(choice)
    choice=defaultChoice;
elseif strcmp(choice,'DFT')
    choice='dft';
elseif strcmp(choice,'PSD')
    choice='psd';
end

% Count # of stims
nStims=size(stimx,2);

% Loop through stims
for i=1:nStims
    
    % Grab current stim
    x=stimx(:,i);
    
    % Calculate window length (truncates if not evently divisible by nWins)
    winSamps = floor(length(x)/nWins);
    
    % Create windowing function
    window=blackman(winSamps);
    o=0; % samples of overlap - non-overlapping here
    
    % Get DFT and PSD of stim
    [s,f,~,p] = spectrogram(x,window,o,fVec,fs);
    
    % Take absolute value of DFT to get rid of the imaginary component
    dft=abs(s);
    
    % Get peak DFT and PSD values for each window
    [vDft,~] = max(dft,[],1);
    [vPsd,~] = max(p,[],1);
    
    % Find the frequency at which each peak value occurs
    for j=1:nWins
        fPeaksDft(i,j)=f(dft(:,j)==vDft(j)); %#ok<*AGROW>
        fPeaksPsd(i,j)=f(p(:,j)==vPsd(j));
    end; clear j % wins
    
end; clear i % files

% Return appropriate peak frequencies (as per 'choice')
if strcmp(choice,'psd')
    fPeaks=fPeaksPsd;
    return
elseif strcmp(choice,'dft')
    fPeaks=fPeaksDft;
    return
else
    error('Second argument must be either ''dft'' or ''psd''');
end % end if
end % end function