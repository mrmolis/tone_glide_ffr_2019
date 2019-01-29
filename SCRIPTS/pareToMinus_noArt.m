clear variables

d = 'C:\Users\vhapormadseb\Desktop\tempOutput\A1toA2\lagRange_2-22ms\';
f = 'resp_A1toA2_xcorr_snr.mat';
p = [d f];
horz = load(p);

d = 'C:\Users\vhapormadseb\Desktop\tempOutput\CztoC7\lagRange_2-22ms\';
f = 'resp_CztoC7_xcorr_snr.mat';
p = [d f];
vert = load(p);

resp = vertcatx(horz.resp,vert.resp);
clear horz vert

d = 'C:\Users\vhapormadseb\Desktop\tempOutput\lagRange_2-22ms\';
f = 'respHorzVert_xcorr_snr.mat';
p = [d f];
save(p,'-v7.3','resp'); % v7.3 tag seems req'd to save structs > some size

stim = vertcatx(resp.stim);
pol = shiftdim({stim.polarity});
isMinus = strcmpi(pol,'minus');
isArt = vertcat(resp.isArtifact);
resp = resp(isMinus & ~isArt);
f = 'respHorzVert-minusOnly-noArt_xcorr_snr.mat';
p = [d f];
save(p,'-v7.3','resp')

clear stim pol isMinus isArt
