Crab cavities:
- initial estimate by B. Salvant thanks to R. Calaga (only broad-band resonator), in HLLHC_Crab_cavities.txt

- completely MESSED UP list of HOMs from B. Salvant (for BNL cavities - file CrabCavityModesForNico.xlsx in ~/Documents/Impedances/HLLHC/Crab_cavities/) thanks to R. Calaga, in HLLHC_Crab_cavities_BNL.txt
To copy those modes from the list (in .csv file, initially from the excel sheet of Benoit):
pr -m -t -s\  /home/nmounet/Documents/Impedances/HLLHC/Crab_cavities/CrabCavityModesForNico_modif_Nico.csv /home/nmounet/Documents/Impedances/HLLHC/Crab_cavities/CrabCavityModesForNico_modif_Nico.csv > HLLHC_Crab_cavities_BNL.txt
Then change the column headers (Rxd,..., Ryd, ...)

- more accurate list of HOMs from B. Salvant (for SLAC cavities - file RFD-cav12-hom-qext-impedance_SLAC.xlsx in ~/Documents/Impedances/HLLHC/Crab_cavities/) thanks to R. Calaga, in HLLHC_Crab_cavities_SLAC.txt

Note: in all cases, hor. and ver. cavities are treated in the same way (we divided the shunt impedance of all modes by 2 and put all the modes of both planes in each plane).
