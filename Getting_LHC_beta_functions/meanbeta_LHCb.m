function [avbetaxLHCb,avbetayLHCb,lenLHCb,sbegLHCb,sendLHCb]=meanbeta_LHCb(twissfilename,beam)

% function to compute the average beta functions of LHCb experimental region.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the cavity
% and beginning and end (in terms of longitudinal position s) of each of them.


ringlength=26658.8832;

% IP8 position (different for b1 w.r.t b2 by ~30cm)
if (beam==1)
    sIP8=23315.37898;
else
    sIP8=23315.22662;
end

% around IP8, HOMs between 15m and 19.3 m on each side (more or less from R. Wanzenber, LHC Project Note 418
% and CERN-ATS-Note-2013-018 TECH)
sbegLHCb=[sIP8-19.3 sIP8+15]
sendLHCb=[sIP8-15 sIP8+19.3]
lenLHCb=sum(sendLHCb-sbegLHCb)

[avbetaxLHCb,avbetayLHCb]=meanbetatwiss2(twissfilename,sbegLHCb,sendLHCb);

