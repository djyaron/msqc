function initializeZipData(obj,zipFileName)
    % runs gaussian in scratch directories, and puts all relavant files
    % into zipFileName
    
%%  SET VARS  
    
    jobname = 'full';
    gjf_file = [jobname,'.gjf'];
    
    basisSet = obj.config.basisSet;
    method   = obj.config.method;
    charge   = obj.config.charge;
    spin     = obj.config.spin;
    par      = obj.config.par;
    ctext    = obj.gaussianFile;
    dataPath = obj.dataPath;
    gaussianPath = obj.gaussianPath;
    gaussianExe  = obj.gaussianExe;

%%  FILE MANAGEMENT / RUN GAUS

    disp('initializing the data');
        
    % Do the calculation
    [origDir, tempDir] = mkScratchDir( gaussianPath );
    writeGjf( gjf_file, ctext );

    setenv('GAUSS_EXEDIR', obj.gaussianPath);

    runGaus( obj )
    moveFiles( obj );

%%  1 NUCLEUS CALCULATIONS

    headerObj = Header( basisSet, method );
    headerObj.link0 = {'rwf=temp.rwf' 'nosave' 'chk=temp.chk'}';
    if obj.config.opt == 1
        headerObj.route = {'opt'};
    end
    headerObj.output = {'nosymm int=noraff iop(99/6=1)' ...
        'scf=(conventional,qc)' 'symm=noint'}; %qc to ensure convergence
    header = headerObj.makeHeader();

    [n1,n2] = size(obj.H1);
    natom = obj.natom;
    obj.H1en = zeros(n1,n2,natom);

    if obj.config.calcEn == 1
        iterateAtom( obj, header );
    end
    
%%  ZIP / CLEAN UP
    
    zip(zipFileName,toZip);

    cd(origDir);
    % cleanup files
    %[status message messageid] = rmdir(tempDir,'s');
    status = rmdir(tempDir,'s');
    while status ~= 1
        disp('  rmdir failed. Retrying...');
        pause(0.1);
        status = rmdir(tempDir,'s');
    end
end

function iterateAtom(obj, header)
    tt1 = obj.templateText;
    tt2 = regexp( tt1, '\w*ATOM\d*', 'match' );
    for iatom = 1:natom
        disp(['doing calc for atom ',num2str(iatom)]);
        ctext = header;
        % To keep even number of electrons (and so spin 1), add an electron
        % by increasing the charge
        atomChar = char( tt2( iatom ) );
        iEnd = regexp( atomChar, 'ATOM' );
        atomSym = atomChar( 1:(iEnd-1) );

        switch(atomSym)
            case 'H'
                tempCharge = 1;
            case 'F'
                tempCharge = -1;
            case 'Cl'
                tempCharge = -1;
            case 'Br'
                tempCharge = -1;
            case 'I'
                tempCharge = -1;
            case 'N'
                tempCharge = 7;
            case 'P'
                tempCharge = 7;
            otherwise
                tempCharge = charge;
        end

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
        fid = fopen(gjf_file,'w');
        fwrite(fid, [ctext,newline,newline], 'char');
        fclose(fid);
        resp1 = 1;

        runGaus( boj, header );

        movefile('fort.32',[jobname,'.f32']);
        % read in data from the polyatom output file

        toZip = {toZip{:},[jobname,'.f32']};
    end
end

function runGaus(obj)
    %Run Gaussain with .bat file and has a timeout built in (when used)
    %Code should be backwards compatible
    %Need to do something about var named terminated

    try
        if obj.config.timeOut == -1
            %%
            %This is temporary to let Fragment work while timeOut does
            %not.
            %When it is working, remove the whole if, leaving only the
            %contents of the else statment.
            resp1 = system([gaussianPath,'\',gaussianExe,' ',gjf_file]);
            %%
        else
            startTime = clock;
            timeOut = obj.config.timeOut; % seconds
            system( [origDir, '\@Fragment\runGaus.bat ', tempDir, '\', jobname, ' &'] );
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

function bool = timeCheck( start, timeOut )
    %Returns 0 or 1 for if the current time (clock) is more than timeOut 
    %seconds after start
    if timeOut == -1
        bool = 0;
    else
        bool = etime( clock, start ) > timeOut;
    end
end

function [origDir, tempDir] = mkScratchDir( gaussianPath )
    tempDir = tempname([gaussianPath,'\','Scratch']);
    mkdir(tempDir);
    origDir = cd(tempDir); % save location so can move back
end

function writeGjf( gjf_file, ctext )
    fid = fopen(gjf_file,'w');
    fwrite(fid, [ctext,'\n\n'], 'char');
    fclose(fid);
end

function moveFiles( obj )
    movefile('temp.chk','full.chk');
    movefile('fort.32','full.f32');

    toZip = {'full.gjf','full.chk','full.f32'};

    if obj.config.opt == 1
        cd( origDir );
        obj.opt_geom( [tempDir, '\full.out'], [tempDir, '\full_opt_config.txt'] );
        cd( tempDir );
        toZip = {toZip{:},'full_opt_config.txt'};
    end
end