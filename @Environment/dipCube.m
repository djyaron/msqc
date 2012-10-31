function res = dipCube(edge, sep, mag)
% edge: cube goes from +edge to - edge along each edge (can be vector)
% sep: distance between +q and -q in dipole (in angstroms)
% mag: magnitude of dipole moment in debye

if (length(edge) == 1)
   edge = edge*ones(3,1);
end

res = Environment;
evec = eye(3);
% two unit charges separated by 1 angstrom have a 4.8 D moment
% dipole = rho * sep * 4.8 so:
rho = mag/(sep*4.8);
for ix=-1:2:1
   for iy=-1:2:1
      for iz=-1:2:1
         pos = ...
            edge(1) * ix * evec(:,1) + edge(2) * iy * evec(:,2) ...
            + edge(3) * iz * evec(:,3);
         orient = getRandomOrientation * (sep/2);
         res.addCharge((pos+orient),rho);
         res.addCharge((pos-orient),-rho);
      end
   end
end
end
function res = getRandomOrientation()
% to avoid issues with the volume element of polar coordinates, we
% will just generate a random cartesian vector, and normalize it
r = -1 + (2 *rand(3,1));  % each element will go from -1 to +1 
res = r/norm(r);
end