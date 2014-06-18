function [avbetaxCMS,avbetayCMS,lenCMS,sbegCMS,sendCMS]=meanbeta_CMS(twissfilename,beam)

% function to compute the average beta functions of CMS experimental region.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the cavity
% and beginning and end (in terms of longitudinal position s) of each of them.


ringlength=26658.8832;

% IP5 position (different for b1 w.r.t b2 by ~30cm)
if (beam==1)
    sIP5=13329.28923;
else
    sIP5=13329.59397;
end

% around IP5, HOMs between 10m and 19.3 m on each side (more or less from R. Wanzenber, LHC Project Note 418
% and CERN-ATS-Note-2013-018 TECH)
sbegCMS=[sIP5-19.3 sIP5+10]
sendCMS=[sIP5-10 sIP5+19.3]
lenCMS=sum(sendCMS-sbegCMS)

[avbetaxCMS,avbetayCMS]=meanbetatwiss2(twissfilename,sbegCMS,sendCMS);

