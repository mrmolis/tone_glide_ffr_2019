% calc_direction_effect
close all

%load data
load('Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\matFiles_newProcess2017All\summaryMats\respYnOnly_horzVert-minusOnly-noArt_xcorr_snr.mat')

% direction
stim = vertcatx(resp.stim);
direction = shiftdim({stim.direction});

%% By Window
bw = vertcatx(resp.byWindow);

% SNR by window
snrBw = reshape(vertcatx(bw.snrDb),size(bw));
indBw = reshape(vertcatx(bw.ind),size(bw));
snrW1 = snrBw(:,1);
snrLater = mean(snrBw(:,2:end),2);
snrW1MinusLater = snrW1 - snrLater;
snrW1MinusLaterSd = std(snrW1MinusLater)
snrW1MinusLaterAvg = mean(snrW1MinusLater)

figure;
scatter(indBw(:),snrBw(:));