function [S, H1, KE, H2, Enuc] = readpolyatom(fid1)
% Reads in matrix elements from the binary file written by Gaussian
% polyatom output routine. To get Gaussian to output this file use
% keywords: int=noraff iop(99/6=1) scf=conventional symm=noint
% Input:
%   fid1 : file handle to a binary *.f32 file (from fopen('filename')) 
% Output:
%   S    : overlap matrix between atomic orbitals
%   H1   : One-electron matrix elements (e-nuc and KE)
%   KE   : kinetic energy matrix elements
%   H2   : two-electron matrix elements H2(i,j,k,l) = (ij|kl)
%   Enuc : nuclear-nuclear repulsion energy

%fid1 = fopen('try8.f32');
%fid1 = fopen('H2_++_6.f32');
size = fread(fid1,1,'integer*4');
title = char( fread(fid1,size,'char*1'))';
size2 = fread(fid1,1,'integer*4');

size = fread(fid1,1,'integer*4');
Nbasis = fread(fid1,1,'integer*4');
junk = fread(fid1,size-4,'char*1');
size2 = fread(fid1,1,'integer*4');


skipline(fid1);
skipline(fid1);


% Nuclear-nuclear repulsion energy
size = fread(fid1,1,'integer*4');
Enuc = fread(fid1,1,'real*8');
size2 = fread(fid1,1,'integer*4');


% the copy from the read-write file of A is only Nbasis big,
% and should be Nbasis^2 big, so the orbs aren't correct
%      Do 10 I = 1, NBasis
%          Write(IntUnt) I, Eig(I)
%   10     Write(IntUnt) (A(J,I),J=1,NBasis)

%skipline(fid1);
%Eorb = zeros(Nbasis,1);
%orb  = zeros(Nbasis,Nbasis);

%for i1=1:Nbasis
%   size = fread(fid1,1,'integer*4');
%   itemp = fread(fid1,1,'integer*4');
%   Eorb(i1,1) = fread(fid1,1,'real*8');
%   size2 = fread(fid1,1,'integer*4');

%   size = fread(fid1,1,'integer*4');
%   orb(:,i1) = fread(fid1,Nbasis,'real*8');
%   size2 = fread(fid1,1,'integer*4');
%end


% Overlap matrix
size = fread(fid1,1,'integer*4');
labS = char( fread(fid1,size,'char*1'))';
size2 = fread(fid1,1,'integer*4');

S = zeros(Nbasis, Nbasis);
last = 0;
while (~(last == 1))
   size = fread(fid1,1,'integer*4');
   isize = fread(fid1,1, 'integer*4');
   last = fread(fid1,1, 'integer*4');
   ia = fread(fid1,500, 'integer*2');
   ib = fread(fid1,500, 'integer*2');
   ic = fread(fid1,500, 'integer*2');
   Sraw  = fread(fid1,500, 'real*8');
   size2 = fread(fid1,1,'integer*4');
   
   
   for i=1:isize
      S(ia(i),ib(i)) = Sraw(i);
      S(ib(i),ia(i)) = Sraw(i);
   end
end


% Kinetic energy
size = fread(fid1,1,'integer*4');
labKe = char( fread(fid1,size,'char*1'))';
size2 = fread(fid1,1,'integer*4');

KE = zeros(Nbasis, Nbasis);
last = 0;
while(~(last==1))
   size = fread(fid1,1,'integer*4');
   isize = fread(fid1,1, 'integer*4');
   last = fread(fid1,1, 'integer*4');
   ia = fread(fid1,500, 'integer*2');
   ib = fread(fid1,500, 'integer*2');
   ic = fread(fid1,500, 'integer*2');
   KEraw  = fread(fid1,500, 'real*8');
   size2 = fread(fid1,1,'integer*4');
   
   for i=1:isize
      KE(ia(i),ib(i)) = KEraw(i);
      KE(ib(i),ia(i)) = KEraw(i);
   end
end


% Hcore - KE
size = fread(fid1,1,'integer*4');
labcore = char( fread(fid1,size,'char*1'))';
size2 = fread(fid1,1,'integer*4');

tcore = zeros(Nbasis, Nbasis);
last = 0;
while(~(last==1))
   size = fread(fid1,1,'integer*4');
   isize = fread(fid1,1, 'integer*4');
   last = fread(fid1,1, 'integer*4');
   ia = fread(fid1,500, 'integer*2');
   ib = fread(fid1,500, 'integer*2');
   ic = fread(fid1,500, 'integer*2');
   core  = fread(fid1,500, 'real*8');
   size2 = fread(fid1,1,'integer*4');
%   display(['read in ', int2str(isize),' H1 elements']); 
   
   for i=1:isize
      tcore(ia(i),ib(i)) = core(i);
      tcore(ib(i),ia(i)) = core(i);
   end
end
H1 = KE + tcore;
%   disp(H1);


%% Two electron integrals
size = fread(fid1,1,'integer*4');
labh2 = char( fread(fid1,size,'char*1'))';
size2 = fread(fid1,1,'integer*4');

H2 = zeros(Nbasis,Nbasis,Nbasis,Nbasis);
last = 0;
while(~(last==1))
   size = fread(fid1,1,'integer*4');
   isize = fread(fid1,1, 'integer*4');
   last = fread(fid1,1, 'integer*4');
   ia = fread(fid1,500, 'integer*2');
   ib = fread(fid1,500, 'integer*2');
   ic = fread(fid1,500, 'integer*2');
   id = fread(fid1,500, 'integer*2');
   mu = fread(fid1,500, 'integer*2');
   iflab = fread(fid1,500, 'integer*2');
   v2  = fread(fid1,500, 'real*8');
   size2 = fread(fid1,1,'integer*4');
%   display(['read in ', int2str(isize),' H2 elements']);
   
   for i1 = 1:isize
      i=ia(i1);
      j=ib(i1);
      k=ic(i1);
      l=id(i1);
      H2(i,j,k,l) = v2(i1);
      H2(j,i,k,l) = v2(i1);
      H2(i,j,l,k) = v2(i1);
      H2(j,i,l,k) = v2(i1);
      H2(k,l,i,j) = v2(i1);
      H2(k,l,j,i) = v2(i1);
      H2(l,k,i,j) = v2(i1);
      H2(l,k,j,i) = v2(i1);
   end
end
