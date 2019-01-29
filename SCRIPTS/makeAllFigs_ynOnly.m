% Script for running all the tone glide FFR versions through the processing

% Differs from scriptRunAll in that allowed lag range is 2-22 ms (was
% previously 4-12 ms). This is based on King et al., 2016, who used 0-20 ms
% after subtracting out tube delay (we did not subtract out tube delay, so
% nudged by 2 ms).


respPath = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respHorzVert-minusOnly-noArt_xcorr_snr.mat'];
outDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
    'Figures_exploratory\ynOnly'];

%% 1) PREP

if ~exist('resp','var')
    load(respPath);
end

upresp;

%% 2) FIGS
versions = {'horz','vert'};

for i=1:numel(versions)
    
    ver = versions{i};
    currLogic = strcmpi(montType,ver) & strcmpi(group,'YN');
    in.resp = resp(currLogic);
    
    % SNR by window (fig 5?)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig5_snrByWindow\_7_']);
    outPath = [outDir '\' 'fig5_snrByWindow' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' fig5_snrByWindow']);
    fig5_snrByWindow(in.resp,outPath);
    
    % Xcorr by window (fig 6?)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig6_xcorrByWindow\_5_']);
    outPath = [outDir '\' 'fig6_xcorrByWindow' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' fig6_xcorrByWindow']);
    fig6_xcorrByWindow(in.resp,outPath);
    
    % XC/SNR Lineplots (figs 5+6 combined into one?)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'xcorrSnrLineplots\_3_']);
    outPath = [outDir '\' 'xcorrSnrLineplots' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' xcorrSnrLineplots']);
    xcorrSnrLineplots(in.resp,outPath);
    
    % Xcorr up/down diffs (fig 7)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig7_xcorrUpDownDiffs\_8_']);
    outPath = [outDir '\' 'fig7_xcorrUpDownDiffs' '_' ver '_30uV.png'];
    fig7_xcorrUpDownDiffs(in.resp,outPath);
    
    % SNR up/down diffs (fig 8)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig8_snrUpDownDiffs\_3_']);
    outPath = [outDir '\' 'fig8_snrUpDownDiffs' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' fig8_snrUpDownDiffs']);
    fig8_snrUpDownDiffs(in.resp,outPath);
    
    % Xcorr up/down diffs SINGLE PLOT (unused)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'xcorrUpDownDiffs_allInOne (unused)\_7_']);
    outPath = [outDir '\' 'xcorrUpDownDiffs_allInOne' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' xcorrUpDownDiffs_allInOne']);
    xcorrUpDownDiffs_allInOne(in.resp,outPath);
    
end; clear i

% All done!
disp('All done!');