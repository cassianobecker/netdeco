function p = getRoot(varargin)

if isempty(varargin);
    up = 0;
else
    up = varargin{1};
end

fpath = mfilename('fullpath');
parts = strsplit(fpath,filesep);

p = [fullfile(parts{1:end-3-up})];

if ~ispc
   p = ['/' p];
end
    
end