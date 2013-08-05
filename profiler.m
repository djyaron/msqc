%% Profile parallelTime2
clear classes;
name = 'mixer-new2';
profile clear;
profile on;
parallelTime2;
profile off;
profsave(profile('info'), ['profiles/parallelTime2/' name]);