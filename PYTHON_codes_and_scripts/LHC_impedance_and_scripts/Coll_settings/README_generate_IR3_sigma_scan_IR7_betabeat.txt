##############################################################
# generating settings with offset from -4 to +4 sigmas in IR3

coll_settings_find_sigma.py -f collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -x 3.5 -y 3.5 -e 6500 -o collgaps_settings_sigma_mm_kept_offsetIR3_0sig.dat

# create 8 more files from -4 to +4 sigmas offset and change by hand IR3 number of sigmas

# generate new settings
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_1sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_1sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_2sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_2sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_3sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_3sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_4sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_4sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_-1sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_-1sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_-2sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_-2sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_-3sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_-3sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_-4sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_-4sig.dat

###############################################################################
# generate settings for various retractions between TCS and TCP in IR7
# (for both nominal IR3, +2 and +4sigmas in IR3)

coll_settings_find_sigma.py -f collgaps_fromRoderik_feb2014_modifNico_materialnames_mm_kept.dat -x 3.5 -y 3.5 -e 6500 -o collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat

# In IR7, TCP settings are kept the same, and retraction between TCLA and TCS
# is also kept the same. 
# Retraction between IR6 collimators and TCS in IR7 is also kept identical.

# from the file collgaps_settings_sigma_mm_kept_offsetIR3_0sig.dat (see above),
# generate 9 files with retractions between TCP and TCS in IR7 of 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5 and 6 sigmas.
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_6sig.dat
# then change by hand the number of sigmas in IR7 & IR6 in these 9 files

# do the same with 2 sigmas retraction in IR3
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_2sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2p5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_2p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_3sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3p5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_3p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_4sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4p5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_4p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5p5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_5p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_6sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_6sig.dat
# then change by hand in IR3

# do the same with 4 sigmas retraction in IR3
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_2sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2p5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_2p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_3sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3p5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_3p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_4sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4p5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_4p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5p5sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_5p5sig.dat
cp collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_6sig.dat collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_6sig.dat
# then change by hand in IR3

# generate new settings in all configurations
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_2p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_3p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_4p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_5p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_6sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_0sig_retractionIR7_6sig.dat

collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_2sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_2sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_2p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_2p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_3sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_3sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_3p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_3p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_4sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_4sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_4p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_4p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_5p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_5p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_6sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_2sig_retractionIR7_6sig.dat

collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_2sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_2sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_2p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_2p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_3sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_3sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_3p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_3p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_4sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_4sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_4p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_4p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_5p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_5p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_feb2014_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_6sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_feb2014_mm_kept_offsetIR3_4sig_retractionIR7_6sig.dat

###############################################################################
# generate settings from nominal settings in sigma, using new beta functions
colcollimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_2sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_2sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_2p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_2p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_3sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_3sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_3p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_3p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_4sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_4sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_4p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_4p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_5p5sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_5p5sig.dat
collimator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_mm_kept.dat -a 3.5 -s collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_6sig.dat -o coll.txt
mv -f coll.txt collgaps_settings_sigma_mm_kept_offsetIR3_0sig_retractionIR7_6sig.dat
limator_settings_from_sigma.py -e 6500. -r collgaps_fromRoderik_modifNico_materialnames_nominal.dat -b coll_ph1_beta_7000GeV_betabeatIR7_b1.txt -a 3.5 -s collgaps_fromRoderik_modifNico_materialnames_nominal_sigma.dat -o collgaps_settings_nominal_betabetaIR7_b1.dat

# then copy the two TCRYO from the initial file from Roderik (keep the same half-gaps -> not very correct...)
# (or better: add them directly in coll_ph1_beta_7000GeV_betabeatIR7_b1.txt beforehand)
