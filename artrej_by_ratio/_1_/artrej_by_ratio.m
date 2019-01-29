%% artrej_by_ratio

function [out,nAcc,nRej,rejThreshold] = artrej_by_ratio(in,rejRatio,...
    nTotal, prestimSamps, stimSamps, bufferSamps)
% in - Array of amplitude-over-time waveform data, with each row
% representing one sweep/epoch/trial and each column representing one time
% point relative to stimulus onset.
%
% rejRatio - Target proportion of the sweeps to be rejected, as a number
% between 0.0 and 1.0;
%
% nTotal - Total number of sweeps in the original full data set. I am
% adding this so that if for whatever reason data gets run multiple times,
% it doesn't keep paring off 10% of the remaining sweeps each time.
%
% prestimSamps - Number of samples in the prestimulus interval, so it knows
% what portion of the epoch to skip.
%
% stimSamps - Duration of stimulus in samples, to determine area of
% interest.
%
% bufferSamps - Number of extra samples to include in the artifact
% rejection analysis *after* the stimulus has turned off.

if ~ismatrix(in)
    error('First input argument (in) must be MATRIX');
end
if ~isnumeric(in)
    error('First input argument (in) must be NUMERIC matrix');
end

if nargin < 6
    error('Not enough input arguments (need 6)');
elseif nargin > 6
    warning('Excess input arguments will be ignored (only first 6 used)');
end

if ~isscalar(rejRatio)
    error('Second input argument (rejRatio) must be scalar')
end

% using new variables to define vector of indices for analysis
startInd = prestimSamps + 1;
stopInd = stimSamps + bufferSamps;
inds = startInd:stopInd;

% unused, but leaving this in for the benefit of anyone reading the code
nSamps = size(in,2);  %#ok<NASGU>

% Account for if this data has already been run before
nSweeps = size(in,1);
targetAcc = nTotal * (1-rejRatio);
if nSweeps <= targetAcc
    out = in;
    nAcc = nSweeps;
    nRej = nTotal - nSweeps;
    rejThreshold = [];
    return
end

% Main part of code
rejThreshold = max(max(in(:,inds)));
accLogic = max(abs(in(:,inds)),[],2)<rejThreshold;
nAcc = numel(accLogic);
msg = ['TH: ' num2str(rejThreshold) '; nAcc: ' num2str(nAcc) '; Pct: '...
    num2str(100*nAcc/nSweeps)];
disp(msg);
while nAcc > targetAcc
    rejThreshold = max(max(in(accLogic,inds)));
    accLogic = max(abs(in(:,inds)),[],2)<rejThreshold;
    nAcc = sum(accLogic);
    msg = ['TH: ' num2str(rejThreshold) '; nAcc: ' num2str(nAcc)...
        '; Pct: ' num2str(100*nAcc/nSweeps)];
    disp(msg);
end

% Assign remaining output variables
nRej = nSweeps - nAcc;
out = in(accLogic,:);

disp([mfilename ' completed!']);

end