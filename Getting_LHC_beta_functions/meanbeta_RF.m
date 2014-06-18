function [avbetaxRF,avbetayRF,lenRF,sbegRF,sendRF]=meanbeta_RF(twissfilename,beam);


% function to compute the average beta functions of all RF cavities (single-harmonic).
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


%%%%
% RF cavities
%%%%

% RESULT: average beta functions of RF cavities
[s,a,b,avbetaxRF,avbetayRF,namesRF,lenRF,sbegRF,sendRF]=meanbetatwiss(twissfilename,'ACSCA');
lenRF=sum(lenRF)
