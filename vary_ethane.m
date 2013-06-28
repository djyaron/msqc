data_loc = 'C:\Users\Zack\Documents\GitHub\msqc\ethaneRs\testdat';

list = ls( data_loc );
for i = 3:length(list)
    i
    dir_name = list( i, 1:end );
    temp_path = strcat(data_loc, '\', dir_name)
    vary_mol( temp_path );
end