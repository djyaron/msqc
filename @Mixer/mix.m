function res = mix(obj, v0, mod, ii, jj, ienv)
% Description:
%   All aspects of mixing have been combined into a single file to boost
%   performance, although it makes for some less clear code.
%
%   Simple implementation of hybrid orbitals. Given a bond, we form a hybrid
%   orbital on atom 1, pointing towards atom 2, and a hybrid orbital on atom 
%   2 pointing towards atom 1. We then scale the interaction between these
%   orbitals, and rotate everything back to the original orientation. The
%   order of the hybrid formed (sp sp2 sp3) is based on the coordinate number
%   (number of atoms bonded to that atom, as defined in bonded array of the
%   model).
%
% Input:
%   v0:    Relevant operator data.
%   mod:   Model3 instance.
%   ii:    ilist.
%   jj:    jlist.
%   ienv:  Environment to use.
%
% Output:
%   res:   Mix result.

switch (obj.mixType)
    case (1)
        x = obj.par(1);
    case (2)
        % charge dependent mixing
        iatom = model.basisAtom(ii(1));
        ch = model.charges(iatom,ienv+1);
        x0 = obj.par(1);
        xslope = obj.par(2);
        x = x0 + xslope*ch;
    case (22)
        % charge dependent mixing quadratic
        iatom = model.basisAtom(ii(1));
        ch = model.charges(iatom,ienv+1);
        x0 = obj.par(1);
        xslope = obj.par(2);
        xquad = obj.par(3);
        x = x0 + xslope*ch + xquad*ch*ch;
    case (3)
        % bond order dependent mixing
        iatom = model.basisAtom(ii(1));
        jatom = model.basisAtom(jj(1));
        bo = model.bondOrders(iatom,jatom,ienv+1);
        x0 = obj.par(1);
        xslope = obj.par(2);
        x = x0 + xslope*(bo-1);
    case (4)
        % bond length dependent mixing
        iatom = model.basisAtom(ii(1));
        jatom = model.basisAtom(jj(1));
        Zs = zeros(1,2);
        Zs(1) = model.Z(iatom);
        Zs(2) = model.Z(jatom);
        Zs = sort(Zs);
        bl = norm(model.rcart(:,iatom)-model.rcart(:,jatom));
        if (Zs(1) == 1 && Zs(2) == 1)
            bl = bl - 0.74;
        elseif (Zs(1) == 1 && Zs(2) == 6)
            bl = bl - 1.1;
        else
            bl = bl - 1.5;
        end
        x0 = obj.par(1);
        xslope = obj.par(2);
        x = x0 + xslope*bl;
    case (5)
        % bond order and bond length dependent mixing
        iatom = model.basisAtom(ii(1));
        jatom = model.basisAtom(jj(1));
        bo = model.bondOrders(iatom,jatom,ienv+1);
        Zs = zeros(1,2);
        Zs(1) = model.Z(iatom);
        Zs(2) = model.Z(jatom);
        Zs = sort(Zs);
        bl = norm(model.rcart(:,iatom)-model.rcart(:,jatom));
        if (Zs(1) == 1 && Zs(2) == 1)
            bl = bl - 0.74;
        elseif (Zs(1) == 1 && Zs(2) == 6)
            bl = bl - 1.1;
        else
            bl = bl - 1.5;
        end
        x0 = obj.par(1);
        xslopebo = obj.par(2);
        xslopebl = obj.par(3);
        x = x0 + xslopebo*(bo-1) + xslopebl*bl;
    case (11)
        % context dependent, diagonal
        iatom = model.basisAtom(ii(1));
        xcontext = model.atomContext(iatom,ienv);
        nx = length(xcontext);
        x = obj.par(1) + sum(obj.par(2:(1+nx)).*xcontext');
    case (12)
        % context dependent, off-diagonal
        iatom = model.basisAtom(ii(1));
        jatom = model.basisAtom(jj(1));
        % disp(['iatom ',num2str(iatom),' jatom ',num2str(jatom), ...
        %       ' ienv ',num2str(ienv)]);
        xcontext = model.bondContext(iatom,jatom,ienv);
        nx = length(xcontext);
        % disp(['In mix 12 ',obj.desc,' nx ',num2str(nx),' pars ',num2str(obj.par)]);
        % if no bond, then we don't do context sensitive scaling
        if (nx == 0)
            x = obj.par(1);
        else
            x = obj.par(1) + sum(obj.par(2:(1+nx)).*xcontext');
        end
    case (32)
        % Molecular mechanics stretch function
        if (length(ii) ~= length(jj))
            error('error in MM mixer: non-square v0');
        end
        iatom = model.basisAtom(ii(1));
        jatom = model.basisAtom(jj(1));
        Zs = zeros(1,2);
        Zs(1) = model.Z(iatom);
        Zs(2) = model.Z(jatom);
        Zs = sort(Zs);
        bl = norm(model.rcart(:,iatom)-model.rcart(:,jatom));
        if (Zs(1) == 1 && Zs(2) == 1)
            bl = bl - 0.74;
        elseif (Zs(1) == 1 && Zs(2) == 6)
            bl = bl - 1.1;
        else
            bl = bl - 1.5;
        end
        x0 = obj.par(1);
        x1 = obj.par(2);
        x2 = obj.par(3);
        x = x0 + x1*bl + x2 * bl^2;
    otherwise
        error(['unknown mix type in Mixer: ',num2str(obj.mixType)]);
end

if (obj.hybrid)
    vOrig = v0;
    atom1 = mod.basisAtom(ii(1));
    atom2 = mod.basisAtom(jj(1));
    
    switch(obj.hybrid)
        case (obj.hybrid == 1)
            rot1 = getRotationSigma(mod,atom1,atom2);
            rot2 = getRotationSigma(mod,atom2,atom1);
        case (obj.hybrid == 2)
            rot1 = getRotationPi(mod,atom1);
            rot2 = getRotationPi(mod,atom2);
        otherwise
            error(['Mixer: mixGeneric: unknown hybrid ',...
                num2str(obj.hybrid)]);
    end
    % We will use a,b for hybrid orbitals and j,k for non-hybrid orbitals
    % We first determine the original H elements between the hybrid orbs:
    %     <a|H|b> = sum_j,k <a|j><j|H|k><k|b>
    %      Hhyb     =     rot' * H * rot
    v0 = rot1' * v0 * rot2;
end

switch(obj.funcType)
    case (1)
        res = ((1.0-x)/2.0) * v1 + ((1.0+x)/2.0) * v2;
    case (2)
        res = x * v0;
    case (3)
        res = x * v0 + obj.par(end) * eye(size(v0));
    case (4)
        res = ((1.0-x)/2.0) * v1 + ((1.0+x)/2.0) * v2 ...
            + obj.par(end) * eye(size(v0));
    otherwise
        error(['Mixer: mixFunction: unknown functype ',...
            num2str(obj.funcType)]);
end

% We then have to subtract the original value from the orig matrix and
% add on this modified version.
%   <j|H|k> =  <j|a><a|H|b><b|k>
%     H     =  rot * Hhyb * rot'
%     H     =  H - rot*hhyb*rot' + rot*Hmod*rot'

if (obj.hybrid)
    res = vOrig - rot1 * v0 * rot2' + rot1 * res * rot2';
end

end

function rot = getRotationSigma(mod,a1,a2)
% rot1 and rot2 is the rotation matrix to hybrid orbitals for the bond 
% between atom1 and atom2. rot = <non-hybrid|hybrid> 
%    the p orbital on atom a1 pointing along r12 (x') can be obtained by:
%       |x'> = sum_j |j><j|x'> = sum_j |j> r12_j
%     so the new p-orbital is just a linear combination of the old
%     orbitals, with the coefficients being the vector r12_j
%        <j|x'> = r12_j
%     including the s orbitals, we get.
%        rot = <non-hybrid|hybrid> = <j|a> = [c1, c2*r12]
%     for sp3 c1 = 1/2         c2 = sqrt(3)/2
%     for sp2 c1 = 1/sqrt(3)   c2 = sqrt(2/3)
%     for sp  c1 = 1/sqrt(2)   c2 = sqrt(2)
if (mod.Z(a1) == 1)
   % if hydrogen, no need to form hybrids, rot = 1x1 identity matrix
   rot = 1;
else
   coord = mod.coord(a1);
   r12 = mod.rcart(:,a2) - mod.rcart(:,a1);
   r12 = r12/norm(r12);
   if (coord == 4) % sp3 hybridization
      rot = [0.5; (sqrt(3)/2.0) * r12];
   elseif (coord == 3) % sp2 hybridization
      rot = [1.0/sqrt(3.0); sqrt(1.5) * r12];
   elseif (coord == 2)
      rot = (1.0/sqrt(2.0)) * [1; r12];
   else
      % for terminal heavy atom, assume sp3 hybrid
      rot = [0.5; (sqrt(3)/2.0) * r12];
%      error('Mixer: hybrid with coord = 1 for a non-hydrogen atom');
   end
end
end

function rot = getRotationPi(mod,a1)
% For double bonds, we use one mixer to do the sigma bond and a different
% mixer to do the pi bond. For the pi bond, we do the same thing as above, 
% but we use the non-hybridized p-orbital instead of the hybrid.
% Given three atoms bonded to a central atom (r1, r2, r3), we can get the
% axis perpendicular to the plan of these three atoms as:
%   rperp = cross(r2-r1, r3-r1).
% we then normalize this and form 
%   rot = <non-hybrid|hybrid> = <j|a> = [0, rperp]

bonded = mod.isBonded(a1,:);
coord = mod.coord(a1);
if (coord ~= 3)
   error('pi bond between atoms without 3 fold coordination');
end
[~,bondedAtoms] = find(bonded == 1);
r1 = mod.rcart(:,bondedAtoms(1));
r2 = mod.rcart(:,bondedAtoms(2));
r3 = mod.rcart(:,bondedAtoms(3));
rperp = cross(r2-r1,r3-r1);
rperp = rperp/norm(rperp);
rot = [0.0; rperp];

end

