close all
lag = vertcatx(resp.lag);
lagMs = vertcat(lag.ms);

montage = vertcatx(resp.montage);
mont = strcat(...
    shiftdim({montage.channelToAnalyze}),...
    'to',...
    shiftdim({montage.reference})...
    );

isHorz = strcmpi(mont,'A1toA2');
isVert = strcmpi(mont,'CztoC7');

for i=1:2
    
    switch i
        case 1
            logic = isHorz;
            color = 'blue';
        case 2
            logic = isVert;
            color = 'red';
    end
    
    curr = lagMs(logic);
    x = 1:numel(curr);
    scatter(x,curr,color);
    hold on
    
    mu = mean(curr);
    sigma = std(curr);
    n = 30; % 30 subjects
    sem = sigma/sqrt(n);
    x = (0.25+0.5*(i-1)) * numel(curr);
    ci = sem * 2.086; % SEM * 2.086 = 95% confidence interval of mean
    errorbar(x,mu,ci,...
        'Color',color,...
        'LineWidth',3);
    plot(x,mu,...
        'Marker','o',...
        'LineWidth',3,...
        'MarkerEdgeColor',color,...
        'MarkerFaceColor','white');
end
