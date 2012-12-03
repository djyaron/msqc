%% Submit a debug job
ClusterInfo.state
ClusterInfo.clear
ClusterInfo.setQueueName('debug')
ClusterInfo.setWallTime('00:30:00')
%ClusterInfo.setUserDefinedOptions('-l licenses=MATLAB_Distrib_Comp_Engine:16')
%j=batch(@lotsOfPauses,1,{5,16},'matlabpool',8,'CaptureDiary',true);
myjob=batch(@pwd,1,{[]},'Matlabpool',2,'CurrentDirectory', ...
    '/brashear/yaron/msqc', 'CaptureDiary',true)
 
%% Check on a job
% sched = findResource();
% myjob = sched.Jobs(10);

%% Submit test job
ClusterInfo.state
ClusterInfo.clear
ClusterInfo.setQueueName('debug')
ClusterInfo.setWallTime('00:05:00')
ClusterInfo.setUserDefinedOptions('-l licenses=MATLAB_Distrib_Comp_Engine:16')
%j=batch(@lotsOfPauses,1,{5,16},'matlabpool',8,'CaptureDiary',true);
myjob=batch('nlsqtest','Matlabpool',15,'CurrentDirectory', ...
    '/brashear/yaron/nlsqtest', 'CaptureDiary',true)
 
%% Submit a real job
ClusterInfo.state
ClusterInfo.clear
ClusterInfo.setQueueName('batch')
ClusterInfo.setWallTime('12:00:00')
%ClusterInfo.setUserDefinedOptions('-l licenses=MATLAB_Distrib_Comp_Engine:16')
%j=batch(@lotsOfPauses,1,{5,16},'matlabpool',8,'CaptureDiary',true);
% myjob=batch(@context2psc,0,{},'Matlabpool',15,'CurrentDirectory', ...
%     '/brashear/yaron/msqc', 'CaptureDiary',true)
myjob=batch(@contextFactorypsc,0,{},'Matlabpool',31,'CurrentDirectory', ...
    '/brashear/yaron/msqc', 'CaptureDiary',true)
