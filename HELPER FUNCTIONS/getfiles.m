function out = getfiles(in,varargin)
if abs(nargin)>1
    filtspec = varargin{1};
else 
    filtspec = '';
end
d = dir([in '\' filtspec]);
isDir = vertcat(d.isdir);
f = d(~isDir);
fList = {f.name};
out = fList(:);

end % end main function block