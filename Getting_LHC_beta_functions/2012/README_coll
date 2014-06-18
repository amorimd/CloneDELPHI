##############################
To get the collimators beta functions (2012 lattice - one less TCTVB):
- launch 

madx_dev < thick_asbuilt_inj_ph1collimators.madx

madx_dev < thick_asbuilt_coll4TeV_sq0.6m_ph1collimators.madx
(squeeze 0.6m in IP1 and IP5 and 3m in other IPs)

madx_dev < thick_asbuilt_coll7TeV_ph1collimators.madx
(squeeze 0.55m in IP1 and IP5 and 10m in the other IPs)


- for e.g. b1 at injection, do:

cp twiss.asbuilt.b1.injection_ph1collimators.thick coll_ph1_beta_450GeV_b1.txt

suppress (with vi) 45 lines at the beginning, 1 line after column headers, and
the * at the beginning of the column headers

awk '{print $1,"\t",$5,"\t",$6}' coll_ph1_beta_450GeV_b1.txt > out

suppress with vi 2 lines with collimators from the other beam (TCTVB.4R8 and TCLIA.4L8 for B1, 
TCTVB.4L8 and TCLIA.4R2 for B2)

suppress with vi all " (use the command :%s/"//g )

with vi, add "A" to BETX and BETY (gives BETAX and BETAY)

mv out coll_ph1_beta_450GeV_b1.txt

################################
To create the file twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info (info on betas and angle)

cp twiss.asbuilt.b1.coll7tev_ph1collimators.thick twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info

in twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info, suppress (with vi) all header lines, the "*" before the 
column names, and one line after the column names. Also suppress all "RCOLLIMATOR", KEYWORD and all " (and trailing spaces)
using e.g.
:%s/"RCOLLIMATOR"              //g
:%s/"//g

awk '{print $1,$2,$3,$4,$5,$10,$11}' twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info > tmp

then in tmp, with vi add back tabs with
:%s/ /^I/g
(type "tab" and not directly ^I from a copy-paste of the above !)

mv tmp twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info

Suppress the TCLIA from the other ring (in IR8 if doing now B1, in IR2 if doing now B2)
(and also the TCTVB from the other ring...)

then merge with some info (angle) from my settings files

pr  -m -t -s\  twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info /home/nmounet/Documents/Coll_settings/coll_ph1_beta_7000GeV_sq0p55_b1_2012.txt | gawk '{print $1,$2,$3,$4,$5,$6,$7,$10}' > tmp
mv tmp twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info

then in vi suppress all collimators in IR1, IR2, IR5 and IR8 (keep IR3, IR6 and IR7), and also all TCLAs

Note: I've checked that the remaining collimators are in the correct order - same beta functions in my settings files
 and in this twiss file - to check this do instead of the above
 pr  -m -t -s\  twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info /home/nmounet/Documents/Coll_settings/coll_ph1_beta_7000GeV_sq0p55_b1_2012.txt | gawk '{print $1,$2,$3,$4,$5,$6,$7,$10,$13,$14}' > tmp

and again put back the tabs (in vi):
:%s/ /^I/g

and change the length of the primaries (0.6m - later changed to 0.95m to take artificially into account smaller number
of sigmas) and of the TCDQs (3m)
