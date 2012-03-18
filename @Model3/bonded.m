function [ bonding ] = bonded( objIn,iatom,jatom, maxbond )
% Find atoms that are close enough to be counted as bonded
if (nargin < 4)
    maxbond = 1.75;
end
obj = objIn.frag;

distAtom = zeros(obj.natom,obj.natom);
cart = obj.rcart;
distAtom(iatom,jatom) = norm(cart(:,iatom)-cart(:,jatom));


if distAtom == 0
    bonding = 0;
elseif distAtom < maxbond
    bonding = 1;
else
    bonding = 0;
end


end

