% scriptCompareMontages

resp = load(['Z:\Spectral_Dynamics_Grant\FFR\Tone Glide Study\'...
    'matFiles_newProcess2017All\summaryMats\'...
    'respMinusAvgs_xcorr_snr_30uV.mat']);

fig5_snrByWindow_compareMontage(resp,outPath);