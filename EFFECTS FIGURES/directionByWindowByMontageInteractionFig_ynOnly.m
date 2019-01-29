addpath('Z:\NCRAR_Share\Madsen\Matlab Code\upresp\')

respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
outDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Figures_exploratory\ynOnly'];
outName = 'directionByWindowByMontageInteraction';

if ~exist(outDir,'dir')
    mkdir(outDir);
end

if ~exist('resp','var')
    load(respPath);
end

close all
clear snrAvg rhoAvg
upresp;
resp = resp(strcmpi(group,'YN'));
upresp;

directionList = unique(direction);
montTypeList = unique(montType);
nWins = 3;
winList = {'0-40','40-80','80-120'};

for k=1:numel(montTypeList)
    currMontType = montTypeList{k};
    montTypeLogic = strcmpi(montType,currMontType);
    for i=1:numel(directionList)
        currDirection = directionList{i};
        directionLogic = strcmpi(direction,currDirection);
        comboLogic = directionLogic & montTypeLogic;
        ind = (k-1)*2 + i;
        
        for j=1:nWins
            snrAvg(j,ind) = mean(snr(comboLogic,j));
            rhoAvg(j,ind) = mean(rho(comboLogic,j));
        end; clear j
    end; clear i
end; clear k

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
hLine(3).Marker = 'o'; hLine(3).LineStyle = '--';
hLine(4).Marker = '^'; hLine(4).LineStyle = '--';
legend('horz down','horz up','vert down','vert up');
title('Direction x Window x Montage interaction (SNR) (YN only)');
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
ylim([0.35 0.85]);
hLine(1).Marker = 'o';
hLine(2).Marker = '^';
hLine(3).Marker = 'o'; hLine(3).LineStyle = '--';
hLine(4).Marker = '^'; hLine(4).LineStyle = '--';
legend('horz down','horz up','vert down','vert up');
title('Direction x Window x Montage interaction (xcorr) (YN only)');
saveas(gcf,[outDir '\' outName '_xcorr.png']);
saveas(gcf,[outDir '\' outName '_xcorr.fig']);