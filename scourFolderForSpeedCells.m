function scourFolderForSpeedCells(pathname)

x = dir(pathname);
x = x(3:end);

for i = 1:length(x)
    if strcmp(x(i).name(end-3:end), '.mat')
        isSpeedModulated(strcat(pathname,x(i).name))
    end
end

end