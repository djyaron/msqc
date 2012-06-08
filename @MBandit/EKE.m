function res = EKE(obj)
% Kinetic energy of molecule

res = sum(sum( obj.density.*obj.KE ) );
