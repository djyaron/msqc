function p2hf = density2p(obj,env)
% Input:
%     env:  environment number (0 for isolated fragment)
% output:
%     density2p:  (nbasis,nbasis,nbasis,nbasis) 2-particle density matrix
%                  valid only for restricted hartree fock theory

% For Restricted Hartree Fock theory, the 2-particle density is derivable 
% from the 1-particle density, as follows:

p1hf = obj.density(env)/2.0;
nb = obj.nbasis;
p2hf = zeros(nb, nb, nb, nb);
for a=1:nb
   for b=1:nb
      for c=1:nb
         for d=1:nb
            p2hf(a,b,c,d) = 2*p1hf(a,b)*p1hf(c,d) - p1hf(a,c)*p1hf(b,d);
         end
      end
   end
end


