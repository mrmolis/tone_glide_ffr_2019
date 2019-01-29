% calc_direction_effect

%load data
load('Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\matFiles_newProcess2017All\summaryMats\60ms\resp60ms_win40_ynOnly_minusOnly.mat')

% (UN)COMMENT TO TOGGLE LOOKING AT A SPECIFIC MONTAGE:
upresp;
chosenMont = 'vert';
montLogic = strcmp(montType,chosenMont);
resp = resp(montLogic);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% direction
stim = vertcatx(resp.stim);
direction = shiftdim({stim.direction});

% SNR
snr = vertcatx(resp.snrDb);
upSnr = snr(strcmp(direction,'up'));
downSnr = snr(strcmp(direction,'down'));
upMinusDownSnr = upSnr - downSnr;
upMinusDownSnrSd = std(upMinusDownSnr)
upMinusDownSnrAvg = mean(upMinusDownSnr)

% Correlation
rho = vertcatx(resp.rho);
upRho = rho(strcmp(direction,'up'));
downRho = rho(strcmp(direction,'down'));
upMinusDownRho = upRho - downRho;
upMinusDownRhoSd = std(upMinusDownRho)
upMinusDownRhoAvg = mean(upMinusDownRho)

%% By Window
bw = vertcatx(resp.byWindow);

% Correlation window 2 only
rhoBw = reshape(vertcatx(bw.rho),size(bw));
rhoBw = rhoBw(:,2);
upRhoBw = rhoBw(strcmp(direction,'up'));
downRhoBw = rhoBw(strcmp(direction,'down'));
upMinusDownRhoBw = upRhoBw - downRhoBw;
upMinusDownRhoBwSd = std(upMinusDownRhoBw)
upMinusDownRhoBwAvg = mean(upMinusDownRhoBw)

% SNR window 2 only
snrBw = reshape(vertcatx(bw.snrDb),size(bw));
snrBw = snrBw(:,2);
upSnrBw = snrBw(strcmp(direction,'up'));
downSnrBw = snrBw(strcmp(direction,'down'));
upMinusDownSnrBw = upSnrBw - downSnrBw;
upMinusDownSnrBwSd = std(upMinusDownSnrBw)
upMinusDownSnrBwAvg = mean(upMinusDownSnrBw)