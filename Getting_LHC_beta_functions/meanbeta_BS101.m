function [avbetax101,avbetay101,len101,sbeg101,send101]=meanbeta_BS101(twissfilename,beam)

% function to compute the average beta functions of all 101 type beam screens 
% in IR1 & 5 triplets region, for HL-LHC.
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
len101=sum(send101-sbeg101)
[avbetax101,avbetay101]=meanbetatwiss2(twissfilename,sbeg101,send101);

