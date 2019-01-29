colorList = {'blue','red','green'};
groupList = unique(group);

snrFig = figure();
snrAx = axes();
hold on
rhoFig = figure();
rhoAx = axes();
hold on

for i=1:numel(groupList)
    currGroup = groupList{i};
    currColor = colorList{i};
    groupLogic = strcmpi(currGroup,group);
    
    currSnr = snr(groupLogic);
    currRho = rho(groupLogic);
    currLag = lagMs(groupLogic);
    
    scatter(snrAx,currLag,currSnr,[],currColor);
    scatter(rhoAx,currLag,currRho,[],currColor);
    
    axes(rhoAx);
    text(20,0.9-0.05*i,currGroup,...
        'Color',currColor,...
        'HorizontalAlignment','Right');
    axes(snrAx);
    text(20,14-1*i,currGroup,...
        'Color',currColor,...
        'HorizontalAlignment','Right');
    
end; clear i currGroup groupLogic currColor currSnr currRho currLag

snrAx.XLabel.String = 'Lag (ms)';
rhoAx.XLabel.String = 'Lag (ms)';

snrAx.YLabel.String = 'SNR (dB)';
rhoAx.YLabel.String = 'Correlation Coefficient (rho)';

snrAx.YLim = [-4 16];
rhoAx.YLim = [0 1];

snrAx.XLim = [2 22];
rhoAx.XLim = [2 22];

snrAx.Title.String = 'SNR by lag';
rhoAx.Title.String = 'xCorr by lag';