%% Get the high level kinetic energy
ic = 0;
for ipar = 1:3
    frag = HL{ipar,1};
    for ienv = 1:5
        ic = ic + 1;
        KEpart = frag.partitionE1(ienv,frag.KE);
        HLke(ic,1) = sum(sum(KEpart));
    end
end
ndata = ic;
%% Get the low-level partitioned KE
ic = 0;
nbasis = LL{1,1}.nbasis;
natom = frag.natom;
LLke = zeros(ndata,3,natom,natom); % one choice
%LLke = cell(ndata,3);
for ipar = 1:3
    for ienv = 1:5
        ic = ic + 1;
        for iLL = 1:3
            frag = LL{ipar,1};
            KEpart = frag.partitionE1(ienv,frag.KE);
            LLke(ic,iLL,:,:) = KEpart;
            %LLke{ic,iLL} = KEpart;
        end
    end
end

%% we will use lsqnonlin
% goal is to write KE as a sum
%  HLke = sum_A,B c(type(A),type(B)) LLke(A,B)
%  type(A) = C or H, so there are 3 parameters to fit, cHH cCC cHC
%  future: include different LL Methods, so c1HH c1CC c1HC  c2HH c2CC c2HC
%          were 1,2 are the second indices of the LL cell array