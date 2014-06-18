# to obtain the HOMs file from the notes (in .pdf) of R. Wanzenberg (et al)

# CMS
copy and paste the longitudinal modes in e.g. tmplong.dat, if needed suppress additional
space after "-" for the first lines
# then
awk '{print $5*$7,"\t",$7,"\t",$2*1e6}' tmplong.dat > tmplong2.dat

copy and paste the transverse modes in e.g. tmptrans.dat, suppress additional
space after "-" for the first lines
# then
awk '{print $7*1e3,"\t",$6,"\t",$2*1e6}' tmptrans.dat > tmptrans2.dat

if needed, complete with modes with zero shunt impedance in the plane where there are less modes
(to have same number of lines in all files)

# assemble together the modes into a single file (that already contained BB model)
pr -m -t -s\  tmplong2.dat tmptrans2.dat tmptrans2.dat >> LHC_CMS.txt
# or (HL-LHC)
pr -m -t -s\  tmplong2.dat tmptrans2.dat tmptrans2.dat >> HLLHC_CMS.txt



# ATLAS
copy and paste the longitudinal modes in e.g. tmplong.dat
# then
awk '{print $5*$7,"\t",$7,"\t",$2*1e6}' tmplong.dat > tmplong2.dat

copy and paste the transverse modes in e.g. tmptrans.dat
# then
awk '{print $7,"\t",$6,"\t",$2*1e6}' tmptrans.dat > tmptrans2.dat

# assemble together the modes into a single file (that already contained BB model)
pr -m -t -s\  tmplong2.dat tmptrans2.dat tmptrans2.dat >> HLLHC_ATLAS.txt


# ALICE modes from Benoit
# LHC & HL-LHC (the same)
# replace all , by tabs with (in vi)
:%s/,/ /g

awk '{print $3,"\t",$2,"\t",$1*1e9}' tmp_files/ALICElongi_f_Q_R.txt > tmplong2.dat
awk '{print $3,"\t",$2,"\t",$1*1e9}' tmp_files/ALICEtrans_f_Q_R.txt > tmptrans2.dat

suppress (by hand...) all these stupid ^M
then add lines with zero modes to complete the transverse file

# assemble together the modes into a single file (that already contained columns headersl)
pr -m -t -s\  tmplong2.dat tmptrans2.dat tmptrans2.dat >> LHC_ALICE.txt

# LHCb modes from Benoit
# LHC & HL-LHC (the same)
# replace all , by tabs with (in vi)
:%s/,/ /g

awk '{print $3,"\t",$2,"\t",$1*1e9}' tmp_files/LHCblongi_f_Q_R.txt > tmplong2.dat
awk '{print $3,"\t",$2,"\t",$1*1e9}' tmp_files/LHCbtrans_f_Q_R.txt > tmptrans2.dat

#suppress in vi all these stupid ^M with
%s/^M//g
# (to get ^M type ctrl V - ctrl M)
then add lines with zero modes to complete the transverse file

# assemble together the modes into a single file (that already contained columns headersl)
pr -m -t -s\  tmplong2.dat tmptrans2.dat tmptrans2.dat >> LHC_LHCb.txt
