% Notes to self
% MEASURES: SNR, XCORR
% SORTING VARS: slope,direction,window,montage
% OTHER: FREQ, WAVEFORM

% Vectorize
resp = resp(:);

% EXTRACT VARS...
stim = vertcatx(resp.stim);
% Slope: 
slope = shiftdim({stim.slope});
% Direction:
direction = shiftdim({stim.direction});
% SNR:
bw = vertcat(resp.byWindow);
sz = size(bw);
oSnr = vertcat(resp.snrDb);
wSnr = reshape([bw.snrDb],sz);
% XCORR
oCorr = vertcat(resp.rho);
wCorr = reshape([bw.rho],sz);
% nSubs
subject = vertcat(resp.subject);
sid = shiftdim({subject.id});
nSubs = numel(unique(sid));
% Montage
montage = vertcatx(resp.montage);
act = shiftdim({montage.channelToAnalyze});
ref = shiftdim({montage.reference});
mont = strcat(act, 'to', ref);

% Lists and counts
slopeList = unique(slope);
nSlopes = numel(slopeList);
directionList = unique(direction);
nDirections = numel(directionList);
nWins = size(bw,2);
winList = 1:nWins;
montList = unique(mont);
nMonts = numel(montList);

% So avgs will be arranged with dims ordered: MONT, WIN, DIRECTION, SLOPE
for i=1:nWins
    currWin = winList(i);
    for j=1:nSlopes
        currSlope = slopeList{j};
        slopeLogic = strcmpi(slope,currSlope);
        for k=1:nDirections
            currDirection = directionList{k};
            directionLogic = strcmpi(direction,currDirection);
            for m=1:nMonts
                currMont = montList{m};
                montLogic = strcmpi(mont,currMont);
                comboLogic = montLogic & slopeLogic & directionLogic;
                
                % wSnr
                temp = wSnr(comboLogic,currWin);
                wSnrAvg(i,j,k,m) = mean(temp);
                wSnrSem(i,j,k,m) = std(temp)/sqrt(nSubs);
                
                % wCorr
                temp = wCorr(comboLogic,currWin);
                wCorrAvg(i,j,k,m) = mean(temp);
                wCorrSem(i,j,k,m) = std(temp)/sqrt(nSubs);
                
                % oCorr
                temp = oCorr(comboLogic);
                oCorrAvg(1,j,k,m) = mean(temp);
                oCorrSem(1,j,k,m) = std(temp)/sqrt(nSubs);
                
                % oSnr
                temp = oSnr(comboLogic);
                oSnrAvg(1,j,k,m) = mean(temp);
                oSnrSem(1,j,k,m) = std(temp)/sqrt(nSubs);
                
                % key
                slopeKey{i,j,k,m} = currSlope;
                directionKey{i,j,k,m} = currDirection;
                montKey{i,j,k,m} = currMont;
                winKey(i,j,k,m) = currWin;
                oKey{1,j,k,m} = [currSlope '_' currDirection '_' currMont];
                wKey{i,j,k,m} = ['win' num2str(currWin) '_' oKey{1,j,k,m}];
            end; clear m;
        end; clear k
    end; clear j
end; clear i
