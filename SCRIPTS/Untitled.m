for i=1:numel(resp)
    if strcmp(resp(i).montage.channelToAnalyze,'A1') ...
            && strcmp(resp(i).montage.reference,'A2')
        resp(i).montage.type = 'horz';
    elseif strcmp(resp(i).montage.channelToAnalyze,'Cz') ...
            && strcmp(resp(i).montage.reference,'C7')
        resp(i).montage.type = 'vert';
    end
end