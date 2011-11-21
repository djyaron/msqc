classdef Aggregator < handle
   properties
      HL      % {ipar,1}   High level fragments
      LL      % {ipar,iLL} Low level fragments
      %env     % {ipar,1}   list of environments to include
   end
   methods
      function addFrags(obj, high, low1, low2)
         npar = size(obj.HL,1);
         obj.HL{npar+1,1} = high;
         obj.LL{npar+1,1} = low1;
         obj.LL{npar+1,2} = low2;
      end
      function [res, ehigh, epred] = err1(obj,par)
         ic = 0;
         for ipar = 1:size(obj.HL,1)
            t1 = obj.HL{ipar,1};
            ic = ic + t1.nenv;
         end
         res = zeros(ic,1);
         ehigh = zeros(ic,1);
         epred = zeros(ic,1);
         ic = 0;
         for ipar = 1:size(obj.HL,1)
            fhigh = obj.HL{ipar,1};
            flow1 = obj.LL{ipar,1};
            flow2 = obj.LL{ipar,2};
            for ienv = 1:fhigh.nenv
               ic = ic+1;
               ehigh(ic) = sum(sum(fhigh.partitionE1(ienv,fhigh.KE)));
               elow1 = sum(sum(flow1.partitionE1(ienv,flow1.KE)));
               elow2 = sum(sum(flow2.partitionE1(ienv,flow2.KE)));
               epred(ic) = par(1) * elow1 + par(2) * elow2;
               res(ic) = ehigh(ic)-epred(ic);
            end
         end
      end
      function [res ehigh epred elow] = err2(obj,par)
          if (size(par,2) == 4)
              par = [par 0 0 0 0];
          end
          ic = 0;
          for ipar = 1:size(obj.HL,1)%for all high level energies in object
              t1 = obj.HL{ipar,1};%figure out how many environments
              ic = ic + t1.nenv;
              %ic = ic + size(envs{ipar,1},1); %t1.nenv; %use envs as a parameter to pass in user specified parameters
              %ic = ic + size(obj.envs{ipar}); %using numel instead of size
          end
          res = zeros(ic,1);%empty matricies to hold energies
          ehigh = zeros(ic,1);
          epred = zeros(ic,1);
          elow = zeros(ic,1);
          elow2 = zeros(ic,1);
          ic = 0;
          for ipar = 1:size(obj.HL,1)
              fhigh = obj.HL{ipar,1};%pulling out high and low for particular set of parametres
              flow1 = obj.LL{ipar,1};
              flow2 = obj.LL{ipar,2};
             
              for ienv = 1:fhigh.nenv %looping over number of environments in high, should be equal to low
                  ic = ic+1; %keep count of number of environments 
                  % this is a single number, for total KE
                  highFit = fhigh.KE; %kinetic energy of basis set
                  % highFit = fhigh.H1en(:,:,1);
                  ehigh(ic) = sum(sum(fhigh.partitionE1(ienv,highFit)));%sum over partition function
                  % these are natom x natom matrices
                 
                  low1Fit = flow1.KE; %from old copy of code
                  low2Fit = flow2.KE;
                  % low1Fit = flow1.H1en(:,:,1);%KE;
                  % low2Fit = flow2.H1en(:,:,1);%KE;

                  %lowFit = flow1.H1en(:,:,1);
                  elow1 = flow1.partitionE1(ienv,low1Fit); %atom by atom, in terms of atoms
                  elow(ic)=sum(sum(elow1));
                  elow2 = flow2.partitionE1(ienv,low2Fit);
                  eres = 0.0;
                  
                  b=obj.HL{1}.Z;
                  numatom=obj.HL{1}.natom;
                  for i=1:numatom
                      for j=1:numatom
                          if (i==j) % same atom

                              if (b(i)==1) %that atom is hydrogen
                                  eres = eres + par(1) * elow1(i,j) + par(5) * elow2(i,j);
                              else % that atom is carbon, should have b(i)==6
                                  eres = eres + par(2) * elow1(i,j) + par(6) * elow2(i,j);
                              end
                          else %different atoms
                              if ((b(i)==1)&&(b(j)==1))%(both are h)
                                  eres = eres + 1 * elow1(i,j) + 1 * elow2(i,j); 
                              elseif (((b(i)==1)&&(b(j)==6))||((b(i)==6)&&(b(j)==1)))%one is c and other is h
                                  eres = eres + par(3) * elow1(i,j) + par(7) * elow2(i,j);
                              elseif ((b(i)==6)&&(b(j)==6))%(both are c)
                                  eres = eres + par(4) * elow1(i,j) + par(8) * elow2(i,j);
                              end
                          end
                                          
                          epred(ic) = eres;
                          res(ic) = ehigh(ic)-epred(ic);
                      end                     
                  end                  
               end
          end
      end
   end
end
