respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
outDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Figures_exploratory\allGroups'];
outName = 'groupEffect';

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

groupList = {'YN','NH','HI'};

for i=1:numel(groupList)
    currGroup = groupList{i};
    groupLogic = strcmpi(group,currGroup);
    
    snrAvg(i) = mean(mean(snr(groupLogic,:)));
    rhoAvg(i) = mean(mean(rho(groupLogic,:)));
    
end; clear i

x = 1:numel(groupList);
[hAx,hLineSnr,hLineRho] = plotyy(x,snrAvg,x,rhoAvg);
xlabel('Subject group');
set(gca,'XTick',x);
set(gca,'XTickLabel',groupList);
ylabel(hAx(1),'Signal-to-noise ratio (dB)');
ylabel(hAx(2),'Correlation coefficient (\rho)');
xlim(hAx(1),[0.75 3.25]);
xlim(hAx(2),[0.75 3.25]);
hLineSnr.Marker = 'o';
hLineRho.Marker = '^';
title('Group effect');

saveas(gcf,[outDir '\' outName '.png']);
saveas(gcf,[outDir '\' outName '.fig']);