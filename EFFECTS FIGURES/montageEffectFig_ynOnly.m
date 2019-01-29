respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
outDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Figures_exploratory\ynOnly'];
outName = 'montageEffect';

if ~exist('resp','var')
    load(respPath);
end

close all
clear snrAvg rhoAvg
upresp;
resp = resp(strcmpi(group,'YN'));
upresp;

montList = unique(montType);

for i=1:numel(montList)
    currMont = montList{i};
    montLogic = strcmpi(montType,currMont);
    
    snrAvg(i) = mean(mean(snr(montLogic,:)));
    rhoAvg(i) = mean(mean(rho(montLogic,:)));
    
end; clear i

x = 1:numel(montList);
[hAx,hLineSnr,hLineRho] = plotyy(x,snrAvg,x,rhoAvg);
xlabel('Montage Type');
set(gca,'XTick',x);
set(gca,'XTickLabel',montList);
ylabel(hAx(1),'Signal-to-noise ratio (dB)');
ylabel(hAx(2),'Correlation coefficient (\rho)');
xlim(hAx(1),[0.75 2.25]);
xlim(hAx(2),[0.75 2.25]);
hLineSnr.Marker = 'o';
hLineRho.Marker = '^';
title('Montage effect');

saveas(gcf,[outDir '\' outName '.png']);
saveas(gcf,[outDir '\' outName '.fig']);
