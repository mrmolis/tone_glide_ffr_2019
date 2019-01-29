function out = mergestructs(s1,s2)
% MERGESTRUCTS
%
% Syntax: OUT = MERGESTRUCTS(S1,S2);
%
% Note: if any fields are shared between the input structs, the values in
% the second will overwrite the values in the first.

s3 = struct;
s = {s1,s2};
for K=1:2
    f = fieldnames(s{K});
    for I=1:numel(s{K})
        for J=1:numel(f)
            s3(I).(f{J}) = s{K}(I).(f{J});
        end; clear J
    end; clear I
end; clear K

out = s3;

end