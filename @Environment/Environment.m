classdef Environment < handle
   %ENVIRONMENT Environment of a Fragment.
   %   Environments support point charges and 1-4 order Gaussian fields.
   
   properties
      ncharge   % number of charges
      rho       % (1,ncharge)  value of charge
      r         % (3,ncharge)  position of charge
      nfield    % number of fields
      fieldType % specify orientation and order of multipole as [X,Y,Z]
                % multiple fields given by additional rows in matrix
                % ex. [2,1,0;0,1,2] with mags of [10,-10] -> field=XXY+10
                % field=YZZ-10
      fieldMag  % corresponding magnitudes of these fields
   end
   
   methods (Static)
      res = newCube(size,mag)
      
      res = newBox( obj, frag, cubeExtension, mag, ncharge )
      goodEnv = testEnv( obf, frag, cubeExtension, mag, ncharge, ntrial )
      res = dipoleEnvironment( frag, cubeExtension, mag, ndipole )
      res = dipCube(size, sep, mag)
      
   end % static methods
   
   methods
      function obj = Environment()
         % Return an empty environment.
         obj.ncharge = 0;
         obj.nfield = 0;
      end
      function obj = addCharge(obj,r,rho)
         obj.ncharge = obj.ncharge + 1;
         obj.rho(obj.ncharge) = rho;
         obj.r(:,obj.ncharge) = r(:);
      end
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
         % First check that we have the same number of charges.
         % Need compatibility for unset properties.
         if  size(obj1.ncharge, 1) == 0 && size(obj2.ncharge, 1) == 0
            chgEq = 1;
         elseif size(obj1.ncharge, 1) == 0 || size(obj2.ncharge, 1) == 0
            chgEq = 0;
         else
            chgEq = obj1.ncharge == obj2.ncharge;
         end
         
         % Same for fields.
         if  size(obj1.nfield, 1) == 0 && size(obj2.nfield, 1) == 0
            fldEq = 1;
         elseif size(obj1.nfield, 1) == 0 || size(obj2.nfield, 1) == 0
            fldEq = 0;
         else
            fldEq = obj1.nfield == obj2.nfield;
         end
         
         if (chgEq && fldEq)
            maxdiff = max(max(abs(obj1.r - obj2.r)));
            maxdiff2 = max(max(abs(obj1.rho-obj2.rho)));
            max2 = max(maxdiff, maxdiff2);
            if (max2 > 1.0e-11)
               chgEq = false;
            end
            if (any(any(obj1.fieldType ~= obj2.fieldType)))
               fldEq = 0;
            end
            if(max(abs(obj1.fieldMag - obj2.fieldMag)) > 1.0e-11)
               fldEq = 0;
            end
         end
         
         res = chgEq && fldEq;
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

