respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
outDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Figures_exploratory\allGroups'];
outName = 'directionByWindowInteraction';

if ~exist(outDir,'dir')
    mkdir(outDir);
end

if ~exist('resp','var')
    load(respPath);
end

close all
clear snrAvg rhoAvg
upresp;
if numel(unique(group)) < 3
    clear resp
    load(respPath);
    upresp;
end

directionList = unique(direction);
nWins = 3;
winList = {'0-40','40-80','80-120'};

for i=1:numel(directionList)
    currDirection = directionList{i};
    directionLogic = strcmpi(direction,currDirection);
    
    for j=1:nWins
        snrAvg(j,i) = mean(snr(directionLogic,j));
        rhoAvg(j,i) = mean(rho(directionLogic,j));
    end; clear j
end; clear i

% SNR
hFigSnr = figure();
x = 1:nWins;
hLine = plot(x,snrAvg);
xlabel('Time re: stimulus onset (ms)');
set(gca,'XTick',1:nWins);
set(gca,'XTickLabel',winList);
ylabel('Signal-to-noise ratio (dB)');
xlim([0.75 3.25]);
hLine(1).Marker = 'o';
hLine(2).Marker = '^';
legend('down','up');
title('Direction x Window interaction (SNR)');
saveas(gcf,[outDir '\' outName '_snr.png']);
saveas(gcf,[outDir '\' outName '_snr.fig']);

% XCORR
hFigRho = figure();
x = 1:nWins;
hLine = plot(x,rhoAvg);
xlabel('Time re: stimulus onset (ms)');
set(gca,'XTick',1:nWins);
set(gca,'XTickLabel',winList);
ylabel('Correlation coefficient (\rho)');
xlim([0.75 3.25]);
hLine(1).Marker = 'o';
hLine(2).Marker = '^';
legend('down','up');
title('Direction x Window interaction (xcorr)');
saveas(gcf,[outDir '\' outName '_xcorr.png']);
saveas(gcf,[outDir '\' outName '_xcorr.fig']);