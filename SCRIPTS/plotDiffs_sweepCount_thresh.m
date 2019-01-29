
close all
clear variables

sourcePct = ['C:\Users\vhapormadseb\Desktop\tempOutput\CZtoA2\'...
    'pct10\resp_CztoA2_pct10_xcorr_snr.mat'];
sourceTh = ['C:\Users\vhapormadseb\Desktop\tempOutput\CZtoA2\uV30\'...
    'resp_CztoA2_uV30_xcorr_snr.mat'];

addpath('Z:\NCRAR_Share\Madsen\Matlab Code\homostruct\_1_');

fig1 = figure();
ax1 = gca;
fig2 = figure();
ax2 = gca;
fig3 = figure();
ax3 = gca;
fig4 = figure();
ax4 = gca;

for i=1:2
    switch i
        case 1
            inPath = sourcePct;
        case 2
            inPath = sourceTh;
    end
    
    in = load(inPath);
    resp = in.resp;
    
    % pare down to just minus
    try
        stim = vertcat(resp.stim);
    catch ME
        switch ME.identifier
            case 'MATLAB:catenate:structFieldBad'
                stim = shiftdim({resp.stim});
                stim = homostruct(stim);
            otherwise
                rethrow(ME);
        end
    end
    pol = shiftdim({stim.polarity});
    polLogic = strcmpi(pol,'minus');
    resp = resp(polLogic);
    stim = stim(polLogic); %#ok<NASGU>
    
    % 1) Sweep counts
    try
        epochInfo = vertcat(resp.epochInfo);
    catch ME
        switch ME.identifier
            case 'MATLAB:catenate:structFieldBad'
                epochInfo = shiftdim({resp.epochInfo});
                epochInfo = homostruct(epochInfo);
            otherwise
                rethrow(ME);
        end
    end
    nAcc{i} = vertcat(epochInfo.nAcceptedMatlab);
    axes(ax1)
    plot(ax1,nAcc{i});
    hold on
    
    % 2) Reject thresholds
    if isfield(epochInfo,'artRejThresholdAbsUv')
        th{i} = vertcat(epochInfo.artRejThresholdAbsUv);
    else
        th13 = vertcat(epochInfo.artRejThresholdAbsUv13);
        th46 = vertcat(epochInfo.artRejThresholdAbsUv46);
        th{i} = mean(horzcat(th13,th46),2);
    end
    axes(ax2)
    plot(ax2,th{i});
    hold on
    
    % 3) Overall Xcorr
    rho{i} = vertcat(resp.rho);
    axes(ax3)
    plot(ax3,rho{i});
    hold on
    
    % 4) Overall SNR
    bw = vertcat(resp.byWindow);
    sz = size(bw);
    wCorr = reshape([bw.rho],sz);
    oCorr = vertcat(resp.rho);
    snr = reshape([bw.snrDb],sz);
    snrAvg{i} = mean(snr,2);
    axes(ax4)
    plot(ax4,snrAvg{i});
    hold on
    
    clear in resp pol polLogic stim
    
end

title(ax1,'Sweep Counts')
title(ax2,'Reject Thresholds (uV)')
title(ax3,'Overall correlation')
title(ax4,'Average SNR')





