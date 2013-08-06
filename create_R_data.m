clear( 'all' );

path = 'C:\Users\Zack\Documents\ethaneRs';
r_groups = {
    {1}
%    {9}
%    {17}
    {8 1}
%    {7 1 1}
%    {7 8 8}
%    {6 7}
    {6 8 1}
    };

%Loops through each r group for the 3 r groups
for i = 1:length(r_groups)
    r1 = r_groups{i};
    for j = 1:i
        r2 = r_groups{j};
        for k = 1:j
            r3 = r_groups{k};
            zmat = build_ethane_z(r1, r2, r3);
            zmat.file_name_base = 'ethane';
            zmat.pars = Ethane_pars( zmat.atoms );
            temp_path = tempname([path,'\','testdat']);
            mkdir( temp_path );
            
            config = Config();
%             config.opt = 1;
            config.zmat = zmat;
            
            frag = Fragment();
            frag.runFragment( temp_path, config )
            
            disp( [i, j, k] );
        end
    end
end
