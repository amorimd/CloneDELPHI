function [avbetaxALICE,avbetayALICE,lenALICE,sbegALICE,sendALICE]=meanbeta_ALICE(twissfilename,beam)

% function to compute the average beta functions of ALICE experimental region.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the cavity
% and beginning and end (in terms of longitudinal position s) of each of them.


ringlength=26658.8832;

% IP2 position (different for b1 w.r.t b2 by ~30cm)
if (beam==1)
    sIP2=3332.436584;
else
    sIP2=3332.284216;
end

% around IP2, HOMs between 15m and 19.3 m on each side (from B. Salvant)
sbegALICE=[sIP2-19.3 sIP2+15]
sendALICE=[sIP2-15 sIP2+19.3]
lenALICE=sum(sendALICE-sbegALICE)

[avbetaxALICE,avbetayALICE]=meanbetatwiss2(twissfilename,sbegALICE,sendALICE);

