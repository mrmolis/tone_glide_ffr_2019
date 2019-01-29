function [ out,varargout ] = regexpc(strCell,expCell,varargin)
%REGEXPC Like regexp but for cells
%   

%% expCell conversion as needed
% Convert to cell if string
if ischar(expCell)
    expCell = {expCell};
end

% If single-element cell array, tile to match size of strCell
if iscell(expCell) && numel(expCell)==1
    expCell = repmat(expCell,size(strCell));
end

%% varargin conversion as needed
for i=1:numel(varargin)
    
    currVar = varargin{i};
    
    % Convert to cell if string
    if ischar(currVar)
        currVar = {currVar};
    end
    
    % If single-element cell array, tile to match size of strCell
    if iscell(currVar) && numel(currVar)==1
        currVar = repmat(currVar,size(strCell));
    end
    
    varargin{i} = currVar;

end; clear i currVar

%% Return output
if nargin > 2
    if nargout > 1
        [out,varargout{1:nargout-1}] = cellfun(@regexp,...
            strCell,expCell,varargin{:},...
            'UniformOutput',false);
    else
        out = cellfun(@regexp,...
            strCell,expCell,varargin{:},...
            'UniformOutput',false);
    end
else
    if nargout > 1
        [out,varargout{1:nargout-1}] = cellfun(@regexp,...
            strCell,expCell,...
            'UniformOutput',false);
    else
        out = cellfun(@regexp,...
            strCell,expCell,...
            'UniformOutput',false);
    end
end

end % end function

