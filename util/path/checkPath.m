function checkPath(gvpath)

s = getenv('PATH');
index = strfind(s, gvpath);
foundIt = ~isempty(index);
if ~foundIt
    setenv('PATH', [getenv('PATH') ':' gvpath]);
end

end