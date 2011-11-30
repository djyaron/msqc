function res = newCube(size, mag)
res = Environment;
icharge = 0;
evec = eye(3);
for ix=-1:2:1
   for iy=-1:2:1
      for iz=-1:2:1
         icharge = icharge + 1;
         res.r(:,icharge) = ...
            size(1) * ix * evec(:,1) + size(2) * iy * evec(:,2) ...
            + size(3) * iz * evec(:,3);
      end
   end
end
res.ncharge = icharge;
res.rho = mag * (2 * rand(1,res.ncharge) - 1);
end
