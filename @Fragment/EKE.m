function res = EKE(obj,ienv)
% Kinetic energy in environment ienv

res = sum(sum( obj.density(ienv).*obj.KE ) );
end

