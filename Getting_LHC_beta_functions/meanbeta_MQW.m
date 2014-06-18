function [avbetaxMQW,avbetayMQW,lenMQW,sbegMQW,sendMQW]=meanbeta_MQW(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type MQW.
% in input: name of the twiss file to use for this. beam=1 or 2
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


[s,betax,betay,avbetaxMQW,avbetayMQW,namesMQW,len,sbegMQW,sendMQW]=meanbetatwiss(twissfilename,'MQW');
length(namesMQW)
lenMQW=sum(sendMQW-sbegMQW)
% 48 elements. 
% avbetaxMQW=192.2401 and avbetayMQW=210.154 (in m) (injection b1)
% avbetaxMQW=192.2401 and avbetayMQW=210.154 (in m) (collision sq 2m b1)
% total length: 149.184 (m). Elias had 155m (see Mathematica notebook,
% 17-11-2009)

