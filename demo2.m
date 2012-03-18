%% 
clear classes;
load('ethane4/ethaneDat.mat');
% this loads:
%   LL{ipar,ibasis} Fragment class for geometry ipar, with 
%       ibasis = 1: regular STO-3G
%       ibasis = 2: narrow STO-3G
%       ibasis = 3: diffuse STO-3G
%   HL{ipar,ibasis}  Fragment class for geometry ipar, with
%       ibasis = 1  6-31G
%       ibasis = 2  6-31G*
%       ibasis = 3  6-31G**
%  Each instance of Fragment has nenv environments in it.

% This is a model for one geometry, in nenv environments
ipar = 1;
m = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});

% Mixer(par, type): choose a function taht will mix narrow and diffuse
%    type = 0, neuron mix, 1 par
%    type = 1, linear mix, 1 par
%    type = 2, linear, charge-dependent mix, 2 par

% Current policy for selecting model is just atom identity (carbon,
% hydrogen)
mixKE = Mixer(0,1); % linear interpolation of entire KE operator
% addKEmodDiag(element # (1=H, 6=C), type (1=s 2=p), mixer)
m.addKEmodDiag(1,1,mixKE);
% The above says, for every hydrogen atom, mix the diagonal block of the
% KE matrix, using the mix-function and the parameter in mixKE

mixKE2 = Mixer(0,1);
m.addKEmodDiag(6,[1,2],mixKE2);
% The above creates another parameterized function and uses it to 
% do both the s and p orbitals of carbon. This is equivalent to:
mixKE2 = Mixer(0,1);
m.addKEmodDiag(6,1,mixKE2);
m.addKEmodDiag(6,2,mixKE2);

