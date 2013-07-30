function varargout = elementWiseCombine(varargin)
% Compile the mex code if it hasn't been.

fname = 'elementWiseCombine';
mname = fname;
cname = [mname '.c'];

disp(' ');
disp(['Detected that the mex routine for ' fname ' is not yet built.']);
disp('Attempting to do so now ...');
disp(' ');

if( isempty(dir(cname)) )
    disp(['Cannot find the file ' fname '.c in the same directory as the']);
    disp(['file ' fname '.m. Please ensure that they are in the same']);
    disp('directory and try again. The following file was not found:');
    disp(' ');
    disp(cname);
    disp(' ');
    error(['Unable to compile ' fname '.c']);
else
    disp(['Found file ' fname '.c']);
    disp(' ');
    disp('Now attempting to compile ...');
    disp('If prompted, please press the Enter key and then select any C/C++');
    disp('compiler that is available, such as gcc or visual studio.');
    disp('If none, see http://www.mathworks.com/support/compilers/R2012b/win64.html.');
    disp(' ');
    disp(['mex(''', cname, ''')']);
    disp(' ');
    try
        mex(cname);
        disp([ fname ' mex build completed ... you may now use ' fname '.']);
        disp(' ');
    catch
        disp(' ');
        error(['Unable to compile ' fname ' ... Contact author.']);
    end
    [varargout{1:nargout}] = elementWiseCombine(varargin{:});
end
end