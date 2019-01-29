% SCRIPT upresp (unpack resp)

if ~exist('resp','var')
    respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
        'matFiles_newProcess2017All\summaryMats\'...
        'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
    load(respPath);
end

% Pull out substructs
stim = vertcatx(resp.stim);
sub = vertcatx(resp.subject);
montage = vertcatx(resp.montage);
bw = vertcatx(resp.byWindow);

% Pull out individual variables
slope = shiftdim({stim.slope});
direction = shiftdim({stim.direction});
group = shiftdim({sub.group});
active = shiftdim({montage.channelToAnalyze});
ref = shiftdim({montage.reference});
rho = vertcat(bw.rho); rho = reshape(rho,size(bw));
snr = vertcat(bw.snrDb); snr = reshape(snr,size(bw));
pol = shiftdim({stim.polarity});
sid = shiftdim({sub.id});
isArt = vertcatx(resp.isArtifact);

% Synthesize new variables
mont = strcat(active,'to',ref);
montType = cell(size(mont)); 
montType(strcmpi(mont,'A1toA2')) = {'horz'};
montType(strcmpi(mont,'CztoC7')) = {'vert'};