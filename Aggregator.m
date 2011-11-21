classdef Aggregator < handle
   properties
      HL      % (ipar,1)   High level fragments
      LL      % (ipar,iLL) Low level fragments
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
          ic = 0;
          % find total number of environments
          for ipar = 1:size(obj.HL,1)
              t1 = obj.HL{ipar,1};
              ic = ic + t1.nenv;
          end
          % create empty matrices to hold
          res = zeros(ic,1);
          ehigh = zeros(ic,1);
          epred = zeros(ic,1);
          elow = zeros(ic,1);
          ic = 0;
          for ipar = 1:size(obj.HL,1)
              fhigh = obj.HL{ipar,1};
              flow1 = obj.LL{ipar,1};
              % flow2 = obj.LL{ipar,2};
              for ienv = 1:fhigh.nenv
                  ic = ic+1;
                  % this is a single number, for total KE
                  highFit = fhigh.KE;
                  %highFit = fhigh.H1en(:,:,1);
                  ehigh(ic) = sum(sum(fhigh.partitionE1(ienv,highFit)));
                  % these are natom x natom matrices
                  lowFit = flow1.KE;
                  %lowFit = flow1.H1en(:,:,1);
                  elow1 = flow1.partitionE1(ienv,lowFit);
                  elow(ic)=sum(sum(elow1));
                  % elow2 = flow2.partitionE1(ienv,flow2.KE);
                  eres = 0.0;
                  
                  b=obj.HL{1}.Z;
                  numatom=obj.HL{1}.natom;
                  for i=1:numatom
                      for j=1:numatom
                          if (i==j) % same atom

                              if (b(i)==1) %that atom is hydrogen
                                  eres = eres + par(1) * elow1(i,j);
                              else % that atom is carbon, should have b(i)==6
                                  eres = eres + par(2) * elow1(i,j);
                              end
                          else %different atoms
                              if ((b(i)==1)&&(b(j)==1))%(both are h)
                                  eres = eres + 1 * elow1(i,j); 
                              elseif (((b(i)==1)&&(b(j)==6))||((b(i)==6)&&(b(j)==1)))%one is c and other is h
                                  eres = eres + par(3) * elow1(i,j);
                              elseif ((b(i)==6)&&(b(j)==6))%(both are c)
                                  eres = eres + par(4) * elow1(i,j);
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
