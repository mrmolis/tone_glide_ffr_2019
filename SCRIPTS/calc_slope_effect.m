% calc_direction_effect

%load data
load('Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\matFiles_newProcess2017All\summaryMats\respYnOnly_horzVert-minusOnly-noArt_xcorr_snr.mat')

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

% Correlation by window
rhoBw = mean(reshape(vertcatx(bw.rho),size(bw)),2);
upRhoBw = rhoBw(strcmp(direction,'up'));
downRhoBw = rhoBw(strcmp(direction,'down'));
upMinusDownRhoBw = upRhoBw - downRhoBw;
upMinusDownRhoBwSd = std(upMinusDownRhoBw)
upMinusDownRhoBwAvg = mean(upMinusDownRhoBw)