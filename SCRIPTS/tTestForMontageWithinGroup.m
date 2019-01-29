groupList = unique(group);
montList = unique(mont);

% Group
for i=1:numel(groupList)
    currGroup = groupList{i};
    groupLogic = strcmpi(group,currGroup);
    
    % Combined sort
    horzLogic = groupLogic & isHorz;
    vertLogic = groupLogic & isVert;
    horzLag = lagMs(horzLogic);
    vertLag = lagMs(vertLogic);
    horzSnr = snr(horzLogic);
    vertSnr = snr(vertLogic);
    horzRho = rho(horzLogic);
    vertRho = rho(vertLogic);
    
    % Lag
    [hLag(i),pLag(i)] = ttest(horzLag,vertLag); %#ok<*SAGROW>
    
    % SNR
    [hSnr(i),pSnr(i)] = ttest(horzSnr,vertSnr);
    
    % Rho
    [hRho(i),pRho(i)] = ttest(horzRho,vertRho);
    
    
end; clear i groupLogic currGroup
clear horz vert horzLogic vertLogic