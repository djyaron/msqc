data_loc = 'D:\testdat';
current = pwd;
list = ls( data_loc );
for i = 3:4%length(list)
    cd( current );
    i
    dir_name = list( i, 1:end );
    temp_path = strcat(data_loc, '\', dir_name)
    vary_mol( temp_path );
end