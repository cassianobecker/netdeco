function p = pathPart(fname,idx)

parts = strsplit(fname,filesep);
p = fullfile(parts{1:end-idx});

end