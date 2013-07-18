function initializeZipData(obj,zipFileName)
    % runs gaussian in scratch directories, and puts all relavant files
    % into zipFileName
    
%%  SET VARS  
    disp('initializing the data');
    
    newline = char(10);
    jobname = 'full';
    gjf_file = [jobname,'.gjf'];
    
%%  SET VARS FROM OBJECT
    basisSet = obj.config.basisSet;
    method   = obj.config.method;
    charge   = obj.config.charge;
    spin     = obj.config.spin;
    par      = obj.config.par;
    ctext    = obj.gaussianFile;
    dataPath = obj.dataPath;
    gaussianPath = obj.gaussianPath;
    gaussianExe  = obj.gaussianExe;
    
%%
    % Do the calculation

    tempDir = tempname([gaussianPath,'\','Scratch']);
    mkdir(tempDir);
    origdir = cd(tempDir); % save location so can move back
    fid1 = fopen(gjf_file,'w');
    fwrite(fid1, [ctext,newline,newline], 'char');
    fclose(fid1);

    setenv('GAUSS_EXEDIR', obj.gaussianPath);
    resp1 = 1;
    % counter = 0;
    % count_max = 2;

    %RUNGAUS

    movefile('temp.chk','full.chk');
    movefile('fort.32','full.f32');

    toZip = {'full.gjf','full.chk','full.f32'};

    if obj.config.opt == 1
        cd( origdir );
        obj.opt_geom( [tempDir, '\full.out'], [tempDir, '\full_opt_config.txt'] );
        cd( tempDir );
        toZip = {toZip{:},'full_opt_config.txt'};

    end

    % Do calculations with only one nucleus present at a time

    % header for the Gaussian job file (input file)
    header = ['%rwf=temp.rwf',newline,...
        '%nosave',newline,...
        '%chk=temp.chk',newline,...
        '# hf/',basisSet, newline...
        'nosymm int=noraff iop(99/6=1) ',...
        'scf=conventional',' symm=noint', newline, newline, ...
        'title', newline,newline];

    % ______________________________________________

    [n1,n2] = size(obj.H1);
    natom = obj.natom;
    obj.H1en = zeros(n1,n2,natom);

    if obj.config.calcEn == 1
        %ITERATE ATOM
    end
    zip(zipFileName,toZip);

    cd(origdir);
    % cleanup files
    %[status message messageid] = rmdir(tempDir,'s');
    status = rmdir(tempDir,'s');
    while status ~= 1
        disp('  rmdir failed. Retrying...');
        pause(0.1);
        status = rmdir(tempDir,'s');
    end
end

function iterate_atom(obj)
    tt1 = obj.templateText;
    tt2 = regexp( tt1, '\w*ATOM\d*', 'match' );
    for iatom = 1:natom
        disp(['doing calc for atom ',num2str(iatom)]);
        ctext = header;
        % SCF can have trouble converging for these single-atom calculations
        % so we add a qc keyword, to help ensure convergence
        ctext = strrep(ctext,'scf=conventional','scf=(conventional,qc)');
        % charge and spin come next
        % To keep even number of electrons (and so spin 1), add an electron
        % by increasing the charge
        atomChar = char( tt2( iatom ) );
        iEnd = regexp( atomChar, 'ATOM' );
        atomSym = atomChar( 1:(iEnd-1) );

        switch(atomSym)
            case 'H'
                tempCharge = ;
            case 'F'
                tempCharge = ;
            case 'Cl'
                tempCharge = ;
            case 'Br'
                tempCharge = ;
            case 'I'
                tempCharge = ;
            case 'N'
                tempCharge = ;
            case 'P'
                tempCharge = ;
            otherwise
                tempCharge = charge;
        end
        %         if (strcmp(atomSym,'H') || strcmp(atomSym,'F') || strcmp(atomSym,'Cl') ...
        %                 || strcmp(atomSym,'P') ...
        %                || strcmp(atomSym,'Br') || strcmp(atomSym,'I')) % should be made general, for any odd Z
        %             tempCharge = charge -1;
        %         elseif strcmp(atomSym,'N')
        %             tempCharge =
        %         else
        %             tempCharge = charge;
        %         end
        ctext = [ctext, num2str(tempCharge), ' ', num2str(2), newline];
        % For molecule specification, we first replace all ATOM# with spaces
        t1 = obj.templateText;
        % Iterate in reverse order, or replacements will not work properly
        % with more than 10 atoms.
        for jatom = natom:-1:1
            if (jatom == iatom)
                t1 = strrep(t1, ['ATOM',num2str(jatom)],' ');
            else
                t1 = strrep(t1, ['ATOM',num2str(jatom)],'-Bq');
            end
        end
        % And replace all PAR# with the parameter values
        for ipar = obj.npar:-1:1
            t1 = strrep(t1, ['PAR',num2str(ipar)], num2str(par(ipar),'%23.12f'));
        end
        ctext = [ctext, t1];

        % Do the calculation and read in data
        jobname = ['atom',num2str(iatom)];
        gjf_file = [jobname,'.gjf'];
        fid1 = fopen(gjf_file,'w');
        fwrite(fid1, [ctext,newline,newline], 'char');
        fclose(fid1);
        resp1 = 1;

        %RUNGAUS()

        movefile('fort.32',[jobname,'.f32']);
        % read in data from the polyatom output file

        toZip = {toZip{:},[jobname,'.f32']};
    end
end

function runGaus(obj)
    %Run Gaussain with .bat file and has a timeout built in (when used)
    %Code should be backwards compatible
    %Need to do something about var named terminated

    while 1
        try
            if obj.config.timeOut == -1
%%
                %This is temporary to let Fragment work while timeOut does
                %not.
                %When it is working, remove the whoel if, leaving only the
                %contents of the else statment.
                resp1 = system([gaussianPath,'\',gaussianExe,' ',gjf_file]);
%%
            else
                startTime = clock;
                timeOut = obj.config.timeOut; % seconds
                resp1 = system( [origDir, '\@Fragment\runGaus.bat ', tempDir, '\', jobname, ' &'] );
                terminated = 0;
                while exist( [tempDir, '\', jobname, '.done'], 'file' ) == 0
                    pause( 10 );
                    if timeCheck( startTime, timeOut )
                        %I don't think this TASKKILL will work... need to
                        %figure out something that works...  
                        system('TASKKILL /IM g09.exe /F');
                        terminated = 1;
                        break
                    end
                end
                disp( resp1 );
                if ( resp1 == 2057 )
                    %THIS PART DOES NOT WORK CAUSE resp1 WILL ALWAYS RETURN
                    %0 BY USING THE & IN THE SYSTEM CALL
                    disp( '  removing temporary files' );
                    delete( 'fort.6', 'gxx.d2e', 'gxx.inp', 'gxx.int', 'gxx.scr', ...
                        'temp.chk', 'temp.fch', 'temp.rwf' )
                end
            end
        catch
            disp( 'Failed, retrying...' )
            resp1 = 1;
        end
    end
end

function bool = timeCheck( start, timeOut )
    %Returns 0 or 1 for if the current time (clock) is more than timeOut 
    %seconds after start
    if timeOut == -1
        bool = 0;
    bool = etime( clock, start ) > timeOut;
end