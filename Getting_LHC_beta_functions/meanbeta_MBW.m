function [avbetaxMBW,avbetayMBW,lenMBW,sbegMBW,sendMBW]=meanbeta_MBW(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type MBW.
% in input: name of the twiss file to use for this. beam=1 or 2
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


[s,betax,betay,avbetaxMBW,avbetayMBW,namesMBW,len,sbegMBW,sendMBW]=meanbetatwiss(twissfilename,'MBW');
length(namesMBW)
lenMBW=sum(sendMBW-sbegMBW)
% 21 elements. 
% avbetaxMBW=138.629 and avbetayMBW=135.2122 (in m) (injection b1)
% avbetaxMBW=140.045 and avbetayMBW=136.6282 (in m) (collision sq 2m b1)
% total length: 70.6 (m). Elias had 70m (see Mathematica notebook,
% 17-11-2009)
% small check that same average beta:
%[avbetax,avbetay]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbegMBW,sendMBW)
% same result
