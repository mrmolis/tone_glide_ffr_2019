if ~exist('resp','var')
    load(['C:\Users\vhapormadseb\Desktop\tempOutput\lagRange_2-22ms\'...
        'respHorzVert-minusOnly-noArt_xcorr_snr.mat']);
end

clear *out

subject = vertcatx(resp.subject);
sid = shiftdim({subject.id});
group = shiftdim({subject.group});
bw = vertcatx(resp.byWindow);
snr = vertcat(bw.snrDb);
snr = reshape(snr,[],3);
rho = vertcatx(bw.rho);
rho = reshape(rho,[],3);
montage = vertcatx(resp.montage);
activeElec = shiftdim({montage.channelToAnalyze});
refElec = shiftdim({montage.reference});
mont = strcat(activeElec,'to',refElec);
montType = cell(size(mont));
montType(strcmpi(mont,'CztoC7')) = {'vert'};
montType(strcmpi(mont,'A1toA2')) = {'horz'};
stim = vertcatx(resp.stim);
slope = shiftdim({stim.slope});
direction = shiftdim({stim.direction});

measList = {'snr','rho'};
slopeList = {'13','23','10'};
directionList = unique(direction);
montTypeList = unique(montType);

ind = 0;
for j=1:numel(slopeList)
    currSlope = slopeList{j};
    slopeLogic = strcmpi(slope,currSlope);
    
    for k=1:numel(directionList)
        currDirection = directionList{k};
        directionLogic = strcmpi(direction,currDirection);
        
        for m=1:numel(montTypeList)
            currMontType = montTypeList{m};
            montTypeLogic = strcmpi(montType,currMontType);
            
            comboLogic = slopeLogic & directionLogic & montTypeLogic;
            
            for w=1:3
                currWin = num2str(w);
                
                ind = ind+1;
                vec = 1:sum(comboLogic);
                headerTag{1,ind} = ['_' currSlope '_' currDirection '_'...
                    currMontType '_win' currWin];
                sidOut(vec,ind) = sid(comboLogic);
                groupOut(vec,ind) = group(comboLogic);
                snrOut(vec,ind) = num2cell(snr(comboLogic,w));
                rhoOut(vec,ind) = num2cell(rho(comboLogic,w));
                
            end; clear w
        end; clear m
    end; clear k
end; clear j

sidHeaders = strcat('sid',headerTag);
groupHeaders = strcat('group',headerTag);
snrHeaders = strcat('snr',headerTag);
rhoHeaders = strcat('rho',headerTag);

outCell = [ sidHeaders, groupHeaders,   snrHeaders, rhoHeaders;...
            sidOut,     groupOut,       snrOut,     rhoOut];
