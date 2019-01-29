respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
outDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Figures_exploratory\ynOnly'];
outName = 'slopeEffect';

if ~exist('resp','var')
    load(respPath);
end

close all
upresp;
resp = resp(strcmpi(group,'YN'));
upresp;

slopeList = {'13','23','10'};

for i=1:numel(slopeList)
    currSlope = slopeList{i};
    slopeLogic = strcmpi(slope,currSlope);
    
    snrAvg(i) = mean(mean(snr(slopeLogic,:)));
    rhoAvg(i) = mean(mean(rho(slopeLogic,:)));
    
end; clear i

x = 1:numel(slopeList);
[hAx,hLineSnr,hLineRho] = plotyy(x,snrAvg,x,rhoAvg);
xlabel('Tone-glide slope (octaves per 120 ms)');
set(gca,'XTick',x);
set(gca,'XTickLabel',{'1/3','2/3','1'});
ylabel(hAx(1),'Signal-to-noise ratio (dB)');
ylabel(hAx(2),'Correlation coefficient (\rho)');
xlim(hAx(1),[0.75 3.25]);
xlim(hAx(2),[0.75 3.25]);
hLineSnr.Marker = 'o';
hLineRho.Marker = '^';
title('Slope effect');

saveas(gcf,[outDir '\' outName '.png']);
saveas(gcf,[outDir '\' outName '.fig']);