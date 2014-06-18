function [avbetaxtap101_121,avbetaytap101_121,lentap101_121,sbegtap101_121,sendtap101_121]=meanbeta_taper101_121(twissfilename,beam)

% function to compute the average beta functions of all tapers between 101 
% and 121 type beam screens in IR1 & 5 triplets region, for HL-LHC.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the tapers
% and beginning and end (in terms of longitudinal position s) of each of them.

% IP5 position (different for b1 w.r.t b2 by ~30cm)
if (beam==1)
    sIP5=13329.28923;
else
    sIP5=13329.59397;
end
ringlength=26658.8832;

% BS 101: Q1 in IR1 & 5
BSlengt=9.5;sfromIP1=22.4; % length and position from IP1 of one BS of diameter 101mm (from E. Todesco, seen in 
% G. Arduini slides, 9th WP2 task leaders' meeting, 06/06/2013)
sbeg101=[sfromIP1 sIP5-sfromIP1-BSlengt sIP5+sfromIP1 ringlength-sfromIP1-BSlengt];
send101=sbeg101+BSlengt;

% tapers from BS 101 to BS 121
IClengt=0.037; % taper length for 15 deg angle (theta=0.27)
sbegtap101_121=[send101(1) sbeg101(2)-IClengt send101(3) sbeg101(4)-IClengt];
sendtap101_121=sbegtap101_121+IClengt;
lentap101_121=sum(sendtap101_121-sbegtap101_121)
[avbetaxtap101_121,avbetaytap101_121]=meanbetatwiss2(twissfilename,sbegtap101_121,sendtap101_121);
