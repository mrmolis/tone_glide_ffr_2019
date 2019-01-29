function out = sortfields(in)
% Like built-in function orderfields, except it also penetrates through to
% any depth of nesting to order fields at all levels.
%
% Written by Brandon Madsen, 12/6/2017 
% brandon.madsen@va.gov

assert(isstruct(in),'Input argument must be of type ''struct''');

s = orderfields(in);

flds = fieldnames(s);
for i=1:numel(s)
    for j=1:numel(flds)
        fld = flds{j};
        temp = s(i).(fld);
        if isstruct(temp)
            s(i).(fld) = sortfields(temp);
        end
    end
end

out = s;