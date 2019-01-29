load('Z:\Spectral_Dynamics_Grant\FFR\Tone_Glide_Study\matFiles_newProcess2017All\summaryMats\respYnOnly_horzVert-minusOnly-noArt_xcorr_snr.mat')

upresp;

snrEffects = {'Extent','Montage','Window','Montage_x_Window'};
pccEffects = {'Extent','Montage','Direction_x_Window','Montage_x_Window',...
    'Direction_x_Montage_x_Window','Extent_x_Montage_x_Window'};
extent = slope;
montTypeList = unique(montType);
directionList = unique(direction);
extentList = unique(extent);
winList = 1:3;

%SNR


%PCC