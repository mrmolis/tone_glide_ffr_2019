% PICKOFFDIR
% Put in full path to file, get back just the directory.

function d = pickOffDir(p)

bs = find(p=='\',1,'last');
d = p(1:bs);

end % end function


% [EOF]