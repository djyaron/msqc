%% Setup Paths
root = 'C:\Users\Alex\Programming\msqc';
exp_folder = 'test';
full_path = [ root, filesep, exp_folder ];

%% What is a reasonable value to use for the environment.
mag = [1.0 5.0 10.0 25.0];
nenv = 10;
cubSize = [6,6,6];
for imag=1:size(mag,2)
   for ienv = 1:nenv
      env{imag,ienv} = Environment.newCube(cubSize,mag(imag));
   end
end
% only run this once, or you will overwrite the env. It is commented out
% for this reason.
save( strcat( exp_folder, filesep, 'env1.mat' ),'env');
%% Generate all of the needed quantum data
load( strcat( exp_folder, filesep, 'env1.mat' ) );
config = Fragment.defaultConfig();
config.template = 'fhyde';
config.basisSet = 'sto-3g ';
frag = Fragment( full_path, config);
%% These results suggest a magnitude of 2 to 3 cause significant
% perturbations
figure(100);
i1 = 0;
plot(0,norm(frag.dipole),'ro');
for imag = 1:size(env,1)
   for ienv = 1:size(env,2)
      frag.addEnv(env{imag,ienv});
      i1 = i1+1;
      figure(100);
      hold on;
      plot(imag,norm(frag.dipoleEnv(:,i1)),'ro');
   end
end

%% Will make 200 calculations with magnitude of 2.5
clear all;
nenv = 200;
cubSize = [6,6,6];
for ienv = 1:nenv
   env{ienv,1} = Environment.newCube(cubSize,2.5);
end
% only run this once, or you will overwrite the env. It is commented out
% for this reason.
%save( strcat( exp_folder, filesep, 'env_mag25.mat' ),'env');
%%
load( strcat( exp_folder, filesep, 'env_mag25.mat' ) );
config = Fragment.defaultConfig();
config.template = 'fhyde';
config.basisSet = '6-31G**';
frag = Fragment( full_path, config);
%% These results suggest a magnitude of 2 to 3 cause significant
for ienv = 1:size(env,1)
   frag.addEnv(env{ienv,1});
end
%% High level calculations
clear classes;
load( strcat( exp_folder, filesep, 'env_mag25.mat' ) );
basisSets = {'6-31G','6-31G**','6-31++G','6-31++G**'};
config = Fragment.defaultConfig();
config.template = 'fhyde';
for ib = 1:size(basisSets,2)
   config.basisSet = basisSets{ib};
   frag = Fragment( full_path, config);
   for ienv = 1:size(env,1)
      disp(['basis ',basisSets{ib},' env ',num2str(ienv)]);
      frag.addEnv(env{ienv});
   end
end
%% Low level calculations
clear classes;
load( strcat( exp_folder, filesep, 'env_mag25.mat' ) );
config = Fragment.defaultConfig();
config.template = 'fhydeGen';
config.basisSet = 'gen';
params = {[1.0 1.0 1.0 1.0 1.0], [0.8 0.8 0.8 0.8 0.8],...
   [1.2 1.2 1.2 1.2 1.2]};
for ipar = 1:size(params,2)
   config.par = params{ipar};
   frag = Fragment( full_path, config);
   config.par = [0.8 0.8 0.8 0.8 0.8];
   for ienv = 1:size(env,1)
      disp(['ipar ',num2str(ipar),' env ',num2str(ienv)]);
      frag.addEnv(env{ienv});
   end
end
%% 
clear classes;
config = Fragment.defaultConfig();
config.template = 'fhyde';
config.basisSet = '6-31++G**';
HL = Fragment( full_path, config);
HL.loadAllEnv;
%%
config = Fragment.defaultConfig();
config.template = 'fhydeGen';
config.basisSet = 'gen';
params = {[1.0 1.0 1.0 1.0 1.0], [0.8 0.8 0.8 0.8 0.8],...
   [1.2 1.2 1.2 1.2 1.2]};
LL = cell(size(params,2),1);
for ipar = 1:size(params,2)
   config.par = params{ipar};
   LL{ipar} = Fragment( full_path, config);
   for ienv = 1:HL.nenv
      LL{ipar}.addEnv( HL.env(ienv) );
   end
end
%%
EHL = zeros(HL.nenv,1);
iatom = 2;
for ienv=1:HL.nenv
   temp = HL.partitionE1(ienv,HL.KE+HL.H1en(:,:,iatom)+HL.H1en(:,:,2));
   EHL(ienv) = temp(1,2);
end
ELL = zeros(HL.nenv,size(LL,1));
for ipar = 1:size(LL,1)
   for ienv=1:HL.nenv
      temp = LL{ipar}.partitionE1(ienv,LL{ipar}.KE+LL{ipar}.H1en(:,:,iatom) ...
         + LL{ipar}.H1en(:,:,2));
      ELL(ienv,ipar) = temp(1,2);
   end
end
%%
m2 = mean(EHL,1);
m1 = mean(ELL(:,1));
plot(ELL(:,1)-m1,EHL-m2,'r.');
hold on;
m1 = mean(ELL(:,2));
plot(ELL(:,2)-m1,EHL-m2,'b.');
hold on;
m1 = mean(ELL(:,3));
plot(ELL(:,3)-m1,EHL-m2,'g.');
%% can we write EHL = sum a_i ELL(i)
%  this becomes B(200,1) = A(200,4)*x(4)
% If A is an m-by-n matrix with m ~= n and B is a column vector 
% with m components, or a matrix with several such columns, 
% then X = A\B is the solution in the least squares sense to 
% the under- or overdetermined system of equations AX  = B. 
% In other words, X minimizes norm(A*X - B), 
B = EHL;
A = zeros(size(ELL,1),size(ELL,2)+1);
A(:,1:size(ELL,2)) = ELL;
A(:,size(ELL,2)+1) = ones(size(ELL,1),1);
x = A(1:100,:)\B(1:100,1);
B2 = A*x;
figure(2)
plot(B(1:100),B2(1:100),'b.');
hold on;
plot(B(100:200),B2(100:200),'r.');

%% Test of HF
config = Fragment.defaultConfig();
config.template = 'fhydeGen';
config.basisSet = 'gen';
params = {[1.0 1.0 1.0 1.0 1.0], [0.8 0.8 0.8 0.8 0.8],...
   [1.2 1.2 1.2 1.2 1.2]};
for ipar = 1:1
   config.par = params{ipar};
   frag = Fragment( full_path, config);
   [orb,Eorb,Ehf] = Model1.hartreeFock(frag,0);
%   for ienv = 1:HL.nenv
%      LL.addEnv( HL.env(ienv) );
%   end
end

