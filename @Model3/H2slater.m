function h2 = H2slater(F0,G1,F2)
% Generate diagonal H2 elements from Slater Condon parameters
% See INDO paper (Ridley and Zerner, Theoret. Chim Acta 32, 111 (1973) eq 3
h2 = zeros(4,4,4,4);
% fill in the non-zero elements here, based on F0 F2 and G1
G3 = G1/3;
t1 = F0 + 4*F2/25;
t2 = F0 - 2*F2/25;
t3 = 3*F2/25;
h2(1,1,1,1) = F0;
for p = 2:4
   h2(1,1,p,p) = F0;
   h2(p,p,1,1) = F0;
   h2(1,p,1,p) = G3;
   h2(1,p,p,1) = G3;
   h2(p,1,1,p) = G3;
   h2(p,1,p,1) = G3;
   h2(p,p,p,p) = t1;
   for q=(p+1):4
      h2(p,p,q,q) = t2;
      h2(q,q,p,p) = t2;
      h2(q,p,q,p) = t3;
      h2(q,p,p,q) = t3;
      h2(p,q,q,p) = t3;
      h2(p,q,p,q) = t3;
   end
end

end

