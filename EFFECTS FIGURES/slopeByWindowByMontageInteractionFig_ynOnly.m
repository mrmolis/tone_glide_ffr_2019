respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
outDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Figures_exploratory\ynOnly'];
outName = 'slopeByWindowByMontageInteraction';

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

slopeList = {'13','23','10'};
montTypeList = unique(montType);
nWins = 3;
winList = {'0-40','40-80','80-120'};

for k=1:numel(montTypeList)
    currMontType = montTypeList{k};
    montTypeLogic = strcmpi(montType,currMontType);
    for i=1:numel(slopeList)
        currSlope = slopeList{i};
        slopeLogic = strcmpi(slope,currSlope);
        comboLogic = slopeLogic & montTypeLogic;
        ind = (k-1)*3 + i;
        
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
ylim([3.5 9]);
hLine(1).Marker = 'o';
hLine(2).Marker = '^';
hLine(3).Marker = 's';
hLine(4).Marker = 'o'; hLine(4).LineStyle = '--';
hLine(5).Marker = '^'; hLine(5).LineStyle = '--';
hLine(6).Marker = 's'; hLine(6).LineStyle = '--';
legend('horz 13','horz 23','horz 10','vert 13','vert 23','vert 10');
title('Slope x Window x Montage interaction (SNR) (YN only)');
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
ylim([0.3 .95]);
hLine(1).Marker = 'o';
hLine(2).Marker = '^';
hLine(3).Marker = 's';
hLine(4).Marker = 'o'; hLine(4).LineStyle = '--';
hLine(5).Marker = '^'; hLine(5).LineStyle = '--';
hLine(6).Marker = 's'; hLine(6).LineStyle = '--';
legend('horz 13','horz 23','horz 10','vert 13','vert 23','vert 10');
title('Slope x Window x Montage interaction (xcorr) (YN only)');
saveas(gcf,[outDir '\' outName '_xcorr.png']);
saveas(gcf,[outDir '\' outName '_xcorr.fig']);