function fmat = makeSavePath(fmat)

parts = strsplit(fmat,filesep);
fpath = ['/' fullfile(parts{1:end-1})];

%path = fullfile(fpath,fname);

if ~exist(fpath,'dir')
    mkdir(fpath)
end

end