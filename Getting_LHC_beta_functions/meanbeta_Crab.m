function [avbetaxCrab,avbetayCrab,lenCrab,sbegCrab,sendCrab]=meanbeta_Crab(twissfilename,beam);


% function to compute the average beta functions of all Crab cavities.
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


%%%%
% Crab cavities
%%%%

%[s,a,b,avbetaxCrab,avbetayCrab,namesCrab,lenCrab,sbegCrab,sendCrab]=meanbetatwiss(twissfilename,'ACF');
%lenCrab=sum(lenCrab)
% previous does not work (zero length of crab cavities in twiss)
[sbegCrab,sendCrab,namesCrab,len]=extract_beamscreens(twissfilename,'ACF',6.5/3); % take average length for each cavity (actually 2.6m+2.6m+1.3m)
lenCrab=sum(sendCrab-sbegCrab)
% RESULT: average beta functions
[avbetaxCrab,avbetayCrab]=meanbetatwiss2(twissfilename,sbegCrab,sendCrab);

