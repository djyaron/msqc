classdef Mixer < handle
   
   properties
      mixType % 0 for sigmoidal, 1 for linear
      par     % (1,npar) current parameters
      desc    % string description
   end
   
   methods
      function obj = Mixer(parIn,mixType)
          if (nargin < 1)
             parIn = [0];
          end
          if (nargin < 2)
             mixType = 0;
          end
          obj.par = parIn;
          obj.mixType = mixType;
      end
      function res = npar(obj)
         res = size(obj.par,2);
      end
      function res = mix(obj, v1, v2)
         x = obj.par(1);
         if (obj.mixType == 0)
            % mix objects v1 and v2, using parameter x.
            %   for x << 0, we get v1, and x>>0 we get v2, with the
            %   switch from v1 to v2 occuring mostly as x=-1..1
            c1 = (tanh(x)+1)/2.0;
            c2 = 1-c1;
            res = c2 * v1 + c1 * v2;
         else
            % want linear mix, with (v1+v2)/2 when x=0
            % res = (v1+v2)/2 + x (v2-v1)/2
            % The bounds are: res = v1 at x = -1;
            %                 res = v2 at x = 1;
            % potentially faster (since v's are matrices while x is scalar)
            res = ((1.0-x)/2.0) * v1 + ((1.0+x)/2.0) * v2;
         end
      end
   end
   
end

