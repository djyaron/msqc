%% Fitting multiple molecules, using makeFitme
clear classes;

train = [1:5 41:45];
test  = [6:10 46:50];

trainC{1}  = {'h2',[],'ch4',1:16,'envs',train};
testC{1} = {'h2',[],'ch4',1:16,'envs',test};
filePrefix{1} = 'ch4';

trainC{2}  = {'h2',[],'ethane',1:7,'envs',train};
testC{2} = {'h2',[],'ethane',1:7,'envs',test};
filePrefix{2} = 'c2h6';

trainC{3}  = {'h2',[],'ethylene',1:7,'envs',train};
testC{3} = {'h2',[],'ethylene',1:7,'envs',test};
filePrefix{3} = 'c2h4z2';

trainC{4}  = {'h2',[],'ch4',1:7,'ethane',1:7,'envs',train};
testC{4} = {'h2',[],'ch4',1:7,'ethane',1:7,'envs',test};
filePrefix{4} = 'ch4-c2h6';

trainC{5}  = {'h2',[],'ch4',1:7,'ethane',1:7,'ethylene',1:7,'envs',train};
testC{5} = {'h2',[],'ch4',1:7,'ethane',1:7,'ethylene',1:7,'envs',test};
filePrefix{5} = 'ch4-c2h6-c2h4';

trainC{6}  = {'h2',[],'ch4',1:16,'ethane',1:7,'envs',train};
testC{6} = {'h2',[],'ch4',1:16,'ethane',1:7,'envs',test};
filePrefix{6} = 'ch4f-c2h6';

trainC{7}  = {'h2',[],'ch4',1:16,'ethane',1:7,'ethylene',1:7,'envs',train};
testC{7} = {'h2',[],'ch4',1:16,'ethane',1:7,'ethylene',1:7,'envs',test};
filePrefix{7} = 'ch4f-c2h6-c2h4';

commonIn = {};

%%
jobs = [];
for iC = 1:7
    trainIn = trainC{iC};
    testIn = testC{iC};
    filePre = filePrefix{iC};
    jobs = [jobs batch(@dofit, 0, {trainIn, testIn, filePre, commonIn})];
end

%%
for iC = 1:7
    filePre = filePrefix{iC};
    dataDir = ['tmp/',filePre,'/'];
    diary([dataDir,'out.diary']);
    diary on;
    jobs(iC).diary
    diary off;
end

%%
for job = 1:7
    destroy(jobs(job));
end

