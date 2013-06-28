data_loc = 'C:\Users\Zack\Documents\GitHub\msqc\ethaneRs\testdat';

list = ls( data_loc );
max = 81;
vary_max = 5;
frags = cell( max, vary_max );

for i = 3:81
    dir_name = list( i, 1:end );
    temp_path = strcat(data_loc, '\', dir_name);
    zip = ls( [temp_path, '\*.zip'] );
    if ~isempty( zip )
        for j = 1:5
            vary = [temp_path, '\vary', num2str(j)];
            if exist( vary, 'file' ) == 7
                config_file = [vary, '\config_vary_', num2str( j ), '.txt' ]
                config_ID = fopen( config_file );
                pars = fscanf( config_ID, '%f' );
                
                config = Fragment.defaultConfig();
                config.template = 'ethane'; 
                config.par = pars;
                config.basis = '6-31G';
                frags( i, j ) = Fragment( vary, config );
            end
        end
    end
end