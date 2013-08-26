data_loc = 'D:\testdat';

list = ls( data_loc );
max = 5;
vary_max = 5;
frags = cell( max, vary_max );

for i = 3:max
    dir_name = list( i, 1:end );
    temp_path = strcat(data_loc, '\', dir_name);
    zip = ls( [temp_path, '\*.zip'] );
    if ~isempty( zip )
        for j = 1:vary_max
            vary = [temp_path, '\vary', num2str(j)];
            if exist( vary, 'file' ) == 7
                mat = ls( [vary, '\*.mat'] );
                if ~isempty( mat )
                    config = importdata( [vary, '\', mat(1:end) ] );
                    frags{ i, j } = Fragment( vary, config );
                    %frags{ i, j }.H2 = [];
                end
            end
        end
    end
end