res = [];
for ifrag = 1:5
    frag = frags{4, ifrag};
    for iatom = 1
        for type = 1 % 0 is s  1 is p
            rho = diag(frag.density);
            for subtype = 1 % 1 for s, 1,2,3 for p
                ns = find((frag.basisAtom == iatom) & (frag.basisType == type) & ...
                    (frag.basisSubType == subtype));
                res(end+1) = rho(ns(2))/rho(ns(1));
            end
        end
    end
end
 
