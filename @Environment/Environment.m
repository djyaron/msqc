classdef Environment < handle 
   %ENVIRONMENT Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
      ncharge   % number of charges
      rho       % (1,ncharges)  value of charge
      r         % (3,ncharges)  position of charge
   end
   
   methods (Static)
      res = newCube(size,mag)
   end % static methods
   
   methods
      function plotFig(obj,nfig)
         figure(nfig);
         for ic = 1:obj.ncharge
            if (obj.rho(1,ic) < 0)
               sym = 'ro';
            else
               sym = 'bo';
            end
            if (ic == 1)
               hold off;
            else
               hold on;
            end
            plot3(obj.r(1,ic),obj.r(2,ic),obj.r(3,ic),sym);
         end
      end
      function res = compare(obj1, obj2)
         res = (obj1.ncharge == obj2.ncharge);
         if (res)
            maxdiff = max(max(abs(obj1.r - obj2.r)));
            maxdiff2 = max(max(abs(obj1.rho-obj2.rho)));
            max2 = max(maxdiff, maxdiff2);
            if (max2 > 1.0e-11)
               res = false;
            end
         end
      end
      function res = gaussianText(obj)
         newline = char(10);
         res = '';
         for ic=1:obj.ncharge
            for ix = 1:3
               res = [res, num2str(obj.r(ix,ic), '%23.12f'), ' '];
            end
            res = [res, num2str(obj.rho(ic), '%23.12f'), newline];
         end
      end
      function displace(obj, rdisp)
          if (sum(size(rdisp) == size(obj.r(:,1))) ~= 2)
              error('Environment.displace needs a 3x1 vector');
          end
          for ic = 1:obj.ncharge
              obj.r(:,ic) = obj.r(:,ic) + rdisp;
          end
      end
   end % methods
end

