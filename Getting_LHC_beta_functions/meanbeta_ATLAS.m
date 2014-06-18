function [avbetaxATLAS,avbetayATLAS,lenATLAS,sbegATLAS,sendATLAS]=meanbeta_ATLAS(twissfilename,beam)

% function to compute the average beta functions of ATLAS experimental region.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the cavity
% and beginning and end (in terms of longitudinal position s) of each of them.


ringlength=26658.8832;

% around IP1, HOMs between 14.5 and 18.5m on each side (from R. Wanzenber et al, HOMs in ATLAS, CERN ATS note to be published)
sbegATLAS=[14.5 ringlength-18.5]
sendATLAS=[18.5 ringlength-14.5]
lenATLAS=sum(sendATLAS-sbegATLAS)

[avbetaxATLAS,avbetayATLAS]=meanbetatwiss2(twissfilename,sbegATLAS,sendATLAS);

