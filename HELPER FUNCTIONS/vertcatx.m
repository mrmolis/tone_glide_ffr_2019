function out = vertcatx(varargin)
% Custom version of VERTCAT to handle structs with mismatched fields.
% Written by Brandon Madsen 12/28/2017 (brandon.madsen@va.gov).

try
    
    % First, see if regular vertcat can handle it:
    out = vertcat(varargin{:});
    return
    
catch ME
    
    % If regular vertcat can't handle it, check if it's our special case:
    switch ME.identifier
        case 'MATLAB:catenate:structFieldBad'
            
            % If it is the special case, handle by ensuring data is column
            % vector of cells and then passing to homostruct to finish:
            addpath('Z:\NCRAR_Share\Madsen\Matlab Code\homostruct\_1_')
            out = varargin(:);
            out = homostruct(out); % separate custom function
            return
            
        otherwise
            % If it's not the special case, rethrow the error:
            rethrow(ME);
    end
end

end % end function vertcat