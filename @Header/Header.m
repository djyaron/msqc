classdef Header < handle
    %Creates header for gjf file
    %   Default header uses HF and STO-3G and 'title'
    %   Has function to output a header variable
    
    properties (SetAccess = private)
        title;      % It's a title...
        method;     % Must match methods that Gaussian recognizes
        basisSet;   % Must match basis sets that Gaussian recognizes
    end
    properties
        link0;      % Cell array of link0 commands - Append as needed
        route;      % Stuff that goes in route selection besides ...
                    % method/basis
        output;     % Stuff that goes after the route selection
        
    end
    
    methods
        function header = Header( basisSet, method, title )
            % Arg order has most frequent to change first
%%          Num of arguments handling
            if nargin < 3
                header.title = 'title';
            else
                header.title = title;
            end
            if nargin < 2
                header.method = 'hf';
            else
                header.method = method;
            end
            if nargin < 1
                header.basisSet = 'STO-3G';
            else
                header.basisSet = basisSet;
            end
%%          
            header.link0 = {};
            header.route = {};
            header.output = {};
        end
        
        function text = makeHeader( obj )
            newLine = char(10);
            head0 = obj.linkHeader();
            head1 = obj.routeHeader();
            head2 = obj.outputHeader();
            
            text = [head0, head1, head2, newLine, obj.title, newLine, newLine];
        end
        
        function text = linkHeader( obj )
            %outputs header section for link0 commands
            newLine = char(10);
            
            text = '';
            for i = 1:length( obj.link0 )
                text = [text, '%', obj.link0{i}, newLine];
            end
        end
        
        function text = routeHeader( obj )
            %outputs header section for route selection
            newLine = char(10);
            
            text = ['# ', obj.method, '/', obj.basisSet, ' '];
            for i = 1:length( obj.route )
                text = [text, obj.route{i}, ' '];
            end
            text = [text, newLine];
        end
        
        function text = outputHeader( obj )
            %outputs header section for stuff below route selection
            newLine = char(10);
            text = '';
            for i = 1:length( obj.output )
                text = [text, obj.output{i}, ' '];
            end
            text = [text, newLine];
        end
        
    end
    
end

