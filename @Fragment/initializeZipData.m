function initializeZipData( fragment,zipFileName )
    % runs gaussian in scratch directories, and puts all relavant files
    % into zipFileName
    
%%  SET VARS  
    
    jobname = 'full';
    gjf_file = [jobname,'.gjf'];
    
    basisSet = fragment.config.basisSet;
    method   = fragment.config.method;
    gaussianPath = fragment.gaussianPath;

%%  FILE MANAGEMENT / RUN GAUS

    disp('initializing the data');
        
    % Do the calculation
    gjf_text = buildGjf( fragment );
    writeGjf( gjf_file, gjf_text );
    [origDir, tempDir] = mkScratchDir( gaussianPath );
    movefile( [origDir,'\',gjf_file], tempDir );
%     gjf_file = [tempDir,'\',gjf_file];
    
    setenv('GAUSS_EXEDIR', fragment.gaussianPath);

    terminated = runGaus( fragment, jobname, origDir, tempDir );
    moveFiles( fragment );

%%  1 NUCLEUS CALCULATIONS

    [n1,n2] = size(fragment.H1);
    natom = fragment.natom;
    obj.H1en = zeros(n1,n2,natom);

    if fragment.config.calcEn == 1
        iterateAtom( fragment, origDir, tempDir );
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

function iterateAtom( fragment, origDir, tempDir )
    for iatom = 1:natom
        disp(['doing calc for atom ',num2str(iatom)]);

        jobname = ['atom',num2str(iatom)];
        gjf_text = buildGjf( fragment );
        writeGjf( [jobname,'.gjf'], gjf_text );

        runGaus( fragment, jobname, origDir, tempDir );

        movefile('fort.32',[jobname,'.f32']);
        % read in data from the polyatom output file

        toZip = {toZip{:},[jobname,'.f32']};
    end
end

function terminated = runGaus(fragment, jobname, origDir, tempDir)
    %Run Gaussain with .bat file and has a timeout built in (when used)
    %Code should be backwards compatible
    %Need to do something about var named terminated

    startTime = clock;
    timeOut = fragment.config.timeOut; % seconds
    system( [origDir, '\@Fragment\runGaus.bat ', tempDir, '\', jobname, ' &'] );
    terminated = 0;
    while exist( [tempDir, '\', jobname, '.done'], 'file' ) == 0
        pause( 10 );
        if timeOut ~= -1 && timeCheck( startTime, timeOut )
            %I don't think this TASKKILL will work... need to
            %figure out something that works...
            system('TASKKILL /IM g09.exe /F');
            terminated = 1;
            break
        end
    end
%   Need cleanup files if fails

%     catch
%         disp( 'Failed, retrying...' )
%         resp1 = 1;
%     end
    
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

function moveFiles( fragment )
    movefile('temp.chk','full.chk');
    movefile('fort.32','full.f32');

    toZip = {'full.gjf','full.chk','full.f32'};

    if fragment.config.opt == 1
        cd( origDir );
        obj.opt_geom( [tempDir, '\full.out'], [tempDir, '\full_opt_config.txt'] );
        cd( tempDir );
        toZip = {toZip{:},'full_opt_config.txt'};
    end
end

function gjf = buildGjf( fragment, bq, charge, spin )
    newLine = char(10);
    if nargin < 2
        bq = 0;
    end
    if nargin < 3
        charge = fragment.config.charge;
    end
    if nargin < 4
        spin = fragment.config.spin;
    end

    basisSet = fragment.config.basisSet;
    method   = fragment.config.method;
    
    
    headerObj = Header( basisSet, method, fragment.config.title );
    headerObj.link0 = {'rwf=temp.rwf' 'nosave' 'chk=temp.chk'}';
    if fragment.config.opt == 1
        headerObj.route = {'opt'};
    end
    headerObj.output = {'nosymm int=noraff iop(99/6=1)' ...
        'scf=conventional' 'symm=noint'};

    headerText = headerObj.makeHeader();
    charge_mult = [num2str(charge), ' ', num2str(spin), newLine];
    zmat_body = fragment.config.zmat.build_gjf( bq );

    gjf = [headerText, charge_mult, zmat_body];
end

function writeGjf( gjf_file, gjf_text )
    fid = fopen(gjf_file,'w');
    fwrite(fid, gjf_text);
    fclose(fid);
end