function outStruct = homostruct(inCell)

for i=1:numel(inCell)
    s = inCell{i};
    f = fieldnames(s);
    n = numel(f);
    flds(1:n,i) = f(:); %#ok
end; clear i

for i=1:numel(flds)
    if ~ischar(flds{i})
        flds{i} = '';
    end
end; clear i

fldList = unique(flds);
fldList = fldList(~strcmpi(fldList,''));

for i=1:numel(inCell)
    s = inCell{i};
    for j=1:numel(fldList)
        f = fldList{j};
        if ~isfield(s,f)
            for k=1:numel(inCell{i})
                inCell{i}(k).(f) = NaN;
            end
        end
    end; clear j
end; clear i

outStruct = sortfields(vertcat(inCell{:}));