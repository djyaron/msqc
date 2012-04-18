% 3-21G
% 
% ****
% H     0 
% S   2   1.00
%       5.4471780              0.1562850        
%       0.8245470              0.9046910        
% S   1   1.00
%       0.1831920              1.0000000    

% STO-3G
% ****
% H     0 
% S   3   1.00
%       3.42525091             0.15432897       
%       0.62391373             0.53532814       
%       0.16885540             0.44463454   

% STO-3G S orbital
d1 = [0.15432897  0.53532814  0.44463454];
e1 = [3.42525091  0.62391373  0.16885540];

% 3-21G 1
d2 =  [0.1562850   0.9046910]
e2 =  [5.4471780   0.8245470]

d3 = 1;
e3 = 0.183;

r = 0:0.1:4;

% phi(j) = sum_i d(i) * exp(-e(i)*r(j).^2)
for j = 1:size(r,2)
   phi1(j) = 0;
   for i=1:size(d1,2)
      phi1(j) = phi1(j) + d1(i) * exp(-e1(i)*r(j)^2);
   end
end
phi1 = phi1/sqrt((r.*phi1)*(r.*phi1)');

for j = 1:size(r,2)
   phi2(j) = 0;
   for i=1:size(d2,2)
      phi2(j) = phi2(j) + d2(i) * exp(-e2(i)*r(j)^2);
   end
end
phi2 = phi2/sqrt((r.*phi2)*(r.*phi2)');

for j = 1:size(r,2)
   phi3(j) = 0;
   for i=1:size(d3,2)
      phi3(j) = phi3(j) + d3(i) * exp(-e3(i)*r(j)^2);
   end
end
phi3 = phi3/sqrt((r.*phi3)*(r.*phi3)');

figure(10);
hold off;
plot(r,phi1,'k-');
hold on;
plot(r,phi2,'r-');
plot(r,phi3,'b-');

for sc = 0.5:0.1:1.5
   if (sc ~= 1)
      for j = 1:size(r,2)
         phi(j) = 0;
         for i=1:size(d1,2)
            phi(j) = phi(j) + d1(i) * exp(-e1(i)*sc^2*r(j)^2);
         end
      end
      phi = phi/sqrt((r.*phi)*(r.*phi)');
      if (sc < 1)
         plot(r,phi,'g-');
      else
         plot(r,phi,'c-');
      end
   end
end
