%function copied from Jon's code in sharkmakeregressor.m
function [x,y]=write3Ddeconv_startTimes(file_loc,event_beg,event_end,fname,censor,noFSL)

if nargin <6
    %censor = 1;
    noFSL=0;
end
format long
x(:,1) = event_beg';
x(:,2) = event_end'-event_beg';
x=x./1000; %Convert to seconds
x(:,3) = ones(length(x),1).*censor; %originally was censor'
%write the -stim_times_FSL
if ~noFSL
    %Save to regs folder
    dlmwrite([file_loc fname '.dat'],x,'delimiter','\t','precision','%.6f')
    y=0;
else
    %write the -stim_times file
    fname = [fname '_noFSL'];
    y = x(logical(x(:,3)),1)';
    %Quick fix hack for just first ten trials troubleshoot SPMG2
    %y = y(1:10);
    dlmwrite([file_loc fname '.dat'],y,'delimiter','\t','precision','%.6f')
end
return
