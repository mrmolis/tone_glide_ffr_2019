function d=just_subfolders(d)
if numel(d) < 1
    return
end
for i=1:numel(d)
    if ~d(i).isdir || ~isempty(regexp(d(i).name,'^\.(\.)?$', 'once'))
        keepLogic(i) = false;
    else
        keepLogic(i) = true; %#ok<*AGROW>
    end % end if
end % end for loop

d = d(keepLogic);

end % end function