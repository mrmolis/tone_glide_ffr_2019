respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
outDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Figures_exploratory\ynOnly'];
outName = 'windowEffect';

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

winList = {'0-40','40-80','80-120'};

for i=1:numel(winList)
    
    snrAvg(i) = mean(snr(:,i));
    rhoAvg(i) = mean(rho(:,i));
    
end; clear i

x = 1:numel(winList);
[hAx,hLineSnr,hLineRho] = plotyy(x,snrAvg,x,rhoAvg);
xlabel('Time re: stimulus onset (ms)');
set(gca,'XTick',x);
set(gca,'XTickLabel',winList);
ylabel(hAx(1),'Signal-to-noise ratio (dB)');
ylabel(hAx(2),'Correlation coefficient (\rho)');
xlim(hAx(1),[0.75 3.25]);
xlim(hAx(2),[0.75 3.25]);
hLineSnr.Marker = 'o';
hLineRho.Marker = '^';
title('Window effect (YN only)');

saveas(gcf,[outDir '\' outName '.png']);
saveas(gcf,[outDir '\' outName '.fig']);