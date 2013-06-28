base_path = pwd;
data_loc = 'C:\Users\Zack\Documents\GitHub\msqc\ethaneRs\testdat';

list = ls( data_loc );
for i = 500:length(list)
    dir_name = list( i, 1:end );
    name = 'ethane';
    temp_path = strcat(data_loc, '\', dir_name)
    config_file = [temp_path, '\', name, '_config_raw.txt'];
    par_ID = fopen( config_file, 'r' );
    temp_par = fscanf( par_ID, '%f' );
    config = Fragment.defaultConfig();
    config.template = name;
    config.par = temp_par;
    config.opt = 1';
    config.calcEn = 0;
    Fragment( temp_path, config );
    i
    fclose(par_ID);
end

