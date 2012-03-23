function res = printMixers(obj)
% primarily useful for debugging. Prints all the mixers in the model

for i = 1:size(obj.mixers,2)
   obj.mixers{i}.print;
end