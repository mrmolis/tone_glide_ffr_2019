% Script for running all the tone glide FFR versions through the processing

% Differs from scriptRunAll in that allowed lag range is 2-22 ms (was
% previously 4-12 ms). This is based on King et al., 2016, who used 0-20 ms
% after subtracting out tube delay (we did not subtract out tube delay, so
% nudged by 2 ms).

versions = {'A1toA2','CztoC7'};

startDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
    'matFiles_newProcess2017All\summaryMats\lagRange_4-12ms'];

for i=1:numel(versions)
    % Prep FFR
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Matlab_Code\prepFFRdata\_1_']);
    ver = versions{i};
    inPath = [startDir '\resp_' ver '_xcorr_snr_30uV.mat'];
    outDir = ['C:\Users\vhapormadseb\Desktop\tempOutput\' ver ...
        '\lagRange_2-22ms'];
    lookupDir = ['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
    'matFiles_newProcess2017All\indivMats\_' ver '\bc_ar_pm\30uV'];
    outName = ['resp_' ver '.mat'];
    outPath = [outDir '\' outName];
    disp([num2str(i) ' ' ver ' prepFFRdata']);
    if ~exist(outPath,'file')
        [~] = prepFFRdata(inPath,lookupDir,outPath);
    else
        disp('^^^ Skipping this, already done');
    end
    
    % XCORR ANALYSIS
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Matlab_Code\xcorrByWin\_35_']);
    inPath = outPath;
    outPath = [inPath(1:end-4) '_xcorr.mat'];
    disp([num2str(i) ' ' ver ' xcorrByWin']);
    if ~exist(outPath,'file')
        [~] = xcorrByWin(inPath,outPath);
    else
        disp('^^^ Skipping this, already done');
    end

    % SNR ANALYSIS
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Matlab_Code\ampSnr\_10_']);
    inPath = outPath ;
    outPath = [inPath(1:end-4) '_snr.mat'];
    disp([num2str(i) ' ' ver ' ampSnr']);
    if ~exist(outPath,'file')
        [~] = ampSnr(inPath,outPath);
    else
        disp('^^^ Skipping this, already done');
    end 
    
    % Prep for figures
    inPath = outPath;
    disp('Loading data to be used for figures...');
    in = load(inPath);
    disp('Loading successful!');
    sep = find(inPath=='\');
    inDir = inPath(1:sep(end)-1);
    
    % Results summary (fig 4???)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig4_resultsByGroup\_12_']);
    outPath = [inDir '\' 'fig4_resultsByGroup' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' fig4_resultsByGroup']);
    fig4_resultsByGroup(in.resp,outPath);
    
    % all GA waveforms (fig 4 instead?)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'waveformGrandAvgs\_6_']);
    outPath = [inDir '\' 'waveformGrandAvgs' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' waveformGrandAvgs']);
    waveformGrandAvgs(in.resp,outPath);
    
    % SNR by window (fig 5?)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig5_snrByWindow\_7_']);
    outPath = [inDir '\' 'fig5_snrByWindow' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' fig5_snrByWindow']);
    fig5_snrByWindow(in.resp,outPath);
    
    % Xcorr by window (fig 6?)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig6_xcorrByWindow\_5_']);
    outPath = [inDir '\' 'fig6_xcorrByWindow' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' fig6_xcorrByWindow']);
    fig6_xcorrByWindow(in.resp,outPath);
    
    % XC/SNR Lineplots (figs 5+6 combined into one?)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'xcorrSnrLineplots\_3_']);
    outPath = [inDir '\' 'xcorrSnrLineplots' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' xcorrSnrLineplots']);
    xcorrSnrLineplots(in.resp,outPath);
    
    % Xcorr up/down diffs (fig 7)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig7_xcorrUpDownDiffs\_8_']);
    outPath = [inDir '\' 'fig7_xcorrUpDownDiffs' '_' ver '_30uV.png'];
    fig7_xcorrUpDownDiffs(in.resp,outPath);
    
    % SNR up/down diffs (fig 8)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'fig8_snrUpDownDiffs\_3_']);
    outPath = [inDir '\' 'fig8_snrUpDownDiffs' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' fig8_snrUpDownDiffs']);
    fig8_snrUpDownDiffs(in.resp,outPath);
    
    % XC/SNR Boxplots (not used?)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'xcorrSnrBoxplots\_4_']);
    outPath = [inDir '\' 'xcorrSnrBoxplots' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' xcorrSnrBoxplots']);
    xcorrSnrBoxplots(in.resp,outPath);
    
    % Xcorr up/down diffs SINGLE PLOT (unused)
    addpath(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
        'Tone Glide Manuscript\Matlab_Code\_Figure Generation\'...
        'xcorrUpDownDiffs_allInOne (unused)\_7_']);
    outPath = [inDir '\' 'xcorrUpDownDiffs_allInOne' '_' ver '_30uV.png'];
    disp([num2str(i) ' ' ver ' xcorrUpDownDiffs_allInOne']);
    xcorrUpDownDiffs_allInOne(in.resp,outPath);
    
end; clear i

% All done!
disp('All done!');