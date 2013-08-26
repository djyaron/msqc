%%
%PSCUtils = 'D:\Users\Alex\Documents\MATLAB\MDCS_Utils_PSC_Blacklight';
PSCUtils = '/afs/andrew/usr13/acappiel/private/MDCS_Utils_PSC_Blacklight';
if isempty(strfind(path, 'MDCS_Utils_PSC_Blacklight'))
    path(path, PSCUtils);
end
ClusterInfo.state

%%
%ClusterInfo.clear
%ClusterInfo.setQueueName('debug')
%ClusterInfo.setWallTime('01:00:00')
%j = batch('multihybrid1', 'matlabpool', 15, 'CurrentDirectory', ...
%    '/usr/users/8/acappiel/dave/', 'CaptureDiary', true);

%%
ClusterInfo.clear
ClusterInfo.setQueueName('debug')
ClusterInfo.setWallTime('00:05:00')
j = batch('parallelTime2', 'matlabpool', 15, 'CurrentDirectory', ...
    '/usr/users/8/acappiel/msqc', 'CaptureDiary', true);
