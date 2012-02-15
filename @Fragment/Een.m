function res = Een(obj,iatom,ienv)
% energy of interaction of the molecule with the environment

res = sum(sum( obj.density(ienv).*obj.H1en(:,:,iatom) );
end
