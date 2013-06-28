%%
close force all
diary off
fclose all;
if (exist('tmp', 'file') == 7)
    rmdir('tmp','s')
end
clear classes
