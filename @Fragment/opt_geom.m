function opt_geom( obj, file, outfile)
%   Takes file from running Gaussian and outputs a file that stores
%   the parameters, keeping the same base for the file name
    
    file_contents = fileread( file );
    raw_nums = raw_num_list( file_contents );
    nums = process_nums( raw_nums );
    config_ID = fopen( outfile, 'w');
    for i = 1:length(nums)
        fprintf( config_ID, '%s\n', nums{i}{1} );
    end
    fclose( config_ID );
end

function list = raw_num_list( raw )
    %Pulls parameters from overall log file by regexp
    %Uses assumptions that log file will not contain any other matches
    %Is in the form X1234=-1.1
    %X can be A/B/D
    %1234 -> Any length greater than 0 of 0-9
    %Must have =, may have -
    %Any length, including 0, of 0-9, a period, any length of 0-9 again
    
    expr = ' [BAD]\d+=-?\d*\.?\d*';
    list = regexp( raw, expr, 'match');
end

function list = process_nums( raw_nums )
    %Pulls str of numbers by regexp
    %Example: it would strips the ' B1=' from ' B1=1.0002' 
    expr = '-?\d*.?\d*$';
    list = regexp( raw_nums, expr, 'match' );
end
