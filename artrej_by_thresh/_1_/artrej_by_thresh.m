function [out,nAcc,nRej,rejThreshold] = artrej_by_thresh(in,rejThreshold,...
    prestimSamps, stimSamps, bufferSamps)
% ARTREJ_BY_THRESH
% in - Array of amplitude-over-time waveform data, with each row
% representing one sweep/epoch/trial and each column representing one time
% point relative to stimulus onset.
%
% rejThreshold - Threshold defining the maximum allowed amplitude.  Sweeps
% where the amplitude exceeds this threshold at any point in the analysis
% epoch will be rejected.
%
% prestimSamps - Number of samples in the prestimulus interval, so it knows
% what portion of the epoch to skip.
%
% stimSamps - Duration of stimulus in samples, to determine area of
% interest.
%
% bufferSamps - Number of extra samples to include in the artifact
% rejection analysis *after* the stimulus has turned off.

checkInArgs();

% using new variables to define vector of indices for analysis
startInd = prestimSamps + 1;
stopInd = stimSamps + bufferSamps;
inds = startInd:stopInd;

nSweeps = size(in,1);
% unused, but leaving this in for the benefit of anyone reading the code
nSamps = size(in,2);  %#ok<NASGU>

% Main part of code
accLogic = max(abs(in(:,inds)),[],2)<rejThreshold;
nAcc = sum(accLogic);

msg = ['TH: ' num2str(rejThreshold) '; nAcc: ' num2str(nAcc)...
        '; pctAcc: ' num2str(100*nAcc/nSweeps)];
disp(msg);

% Assign remaining output variables
nRej = nSweeps - nAcc;
out = in(accLogic,:);

disp([mfilename ' completed!']);


%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkInArgs()
        
        if ~ismatrix(in)
            error('First input argument (in) must be MATRIX');
        end
        if ~isnumeric(in)
            error('First input argument (in) must be NUMERIC matrix');
        end
        
        if nargin < 5
            error('Not enough input arguments (need 5)');
        elseif nargin > 5
            warning('Excess input arguments will be ignored (only first 5 used)');
        end
        
        if ~isscalar(rejThreshold)
            error('Second input argument (rejThreshold) must be scalar')
        end
    end % end nested function

end % end main function



% [EOF]