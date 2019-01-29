% Script for running all the tone glide FFR versions through the processing

% Differs from scriptRunAll in that allowed lag range is 2-22 ms (was
% previously 4-12 ms). This is based on King et al., 2016, who used 0-20 ms
% after subtracting out tube delay (we did not subtract out tube delay, so
% nudged by 2 ms).

% Written by Brandon Madsen (brandon.madsen@va.gov) for internal use only.

%%%%%%%%% EDIT THIS PART AS NEEDED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\flat\compileffr\'...
    'resp_flat_horzVert.mat'];
outDir = ['C:\Users\vhapormadseb\Desktop\tempFlatOutput\']; %#ok<NBRAK>
lookupDir = ['C:\Users\vhapormadseb\Desktop\Molis FFR datamats\'...
    'tone glides study\Flat\bc_ar_pm'];
outName = 'resp_flat.mat';
outPath = [outDir '\' outName];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create output directory if it doesn't already exist
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% Prep FFR
addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Matlab_Code\prepFFRdata\_1_']);
disp('prepFFRdata');
if ~exist(outPath,'file')
    [~] = prepFFRdata(inPath,lookupDir,outPath);
else
    disp('^^^ Skipping this, already done');
end

% XCORR ANALYSIS
addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Matlab_Code\xcorrByWin\_35_']);
inPath = outPath;
outPath = [inPath(1:end-4) '_xcorr.mat'];
disp('xcorrByWin');
if ~exist(outPath,'file')
    [~] = xcorrByWin(inPath,outPath);
else
    disp('^^^ Skipping this, already done');
end

% SNR ANALYSIS
addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Matlab_Code\ampSnr\_10_']);
inPath = outPath ;
outPath = [inPath(1:end-4) '_snr.mat'];
disp('ampSnr');
if ~exist(outPath,'file')
    [~] = ampSnr(inPath,outPath);
else
    disp('^^^ Skipping this, already done');
end

% FINISH
clear inPath lookupdir 
disp('All done!');