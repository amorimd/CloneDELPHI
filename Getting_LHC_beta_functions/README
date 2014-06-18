use madx_dev (4.01.xx), it's better (02/02/2010)

From E. Laface:
- in the twiss file (if you don't change some flag when doing the twiss command)
it is always indicated the end position of the element (in column S)
- the beta function changes a lot in the triplets; with thick elements you see
only the beta function at the beginning and the end, so this is rough. You can
split them (ask Emmanuele Laface) to see the beta function on a finer mesh, it's
quite easy. It's even easier to look at the thin optics where the magnets
are already split (in 5 according to Emmanuele, but he has to check).
- at collision energy, beta function outside the IPs are matched so don't depend on 
the crossing scheme or the squeeze. but inside the IP they depend on squeeze of 
course, but also on crossing angle in principle (change path length in magnets).
- the file coll.str contains the strengths used in the magnets and other
parameters (official ones) at collision. they are the same at whatever 
collision energy considered (3.5 TeV or 7 TeV for instance).
- the file inj.str contains the strengths used in the magnets and other
parameters (official ones) at injection.
- flags:
	- on_x1, on_x2, on_x5, on_x8 : crossing angle (1 full normalized strength, 
	0 deactivated, sth in between will change the strength) in IPs. should
	be 1 (at collision and injection)
	- on_sep1, on_sep2, on_sep5, on_sep8 : separation (in the other plane) in IPs.
	should be 1 at injection and 0 at collision.
	- on_o2, on_o8, on_a2, on_a8: Massamitsu's optimization for separation and crossing at IP2 and
	IP8 (so that the beams are sligthly better centered): small offset and angle in 
	the plane of crossing angle and separation respectively. 
	Can be 0 or 1 at injection (at collision it is only for IP1 and IP5 and there it is 
	taken care of automatiquely), but usually it is deactivated (0). It is a second order effet
	anyway.
	- on_atlas and on_cms (for thin inj) or on_sol_atlas and on_sol_cms (for the rest): 
	to switch on the atlas and cms solenoids. NOTE: second order effect only (for particles
	with transverse speed), no impact on beta functions. So can be put to 0.
	- same thing for on_sol_alice : no impact so 0.
	- on_lhcb : magnet of lhcb. Should be negative at injection because only the
	negative polarity is allowed (pb of aperture and cumulation with crossing angle).
	at collision it can be positive or negative (but in practice difficult to change polarity).
	it should be of absolute value equal to 7 TeV/energy (e.g. 2 at 3.5 TeV). So, at inj
	on_lhcb=-7000/450, at collision =+-1 if 7TeV or +-2 if 3.5TeV
	- on_alice : the same as on_lhcb except that polarity can be whatever: should be then
	+7000/450 at injection and at collision +1 if 7 TeV or +2 if 3.5 TeV.
	on_alice and on_lhcb have to be considered from the time of injection (you don't want
	to change the current in those magnets during ramp).

It would be better to check the impact of all those flags on the beta functions.

In principle there is an offset of the beam due to the crossing angle (~.5mm over the full IP) 
and separation, and also the collision point in IP1 and IP5 is offset by .5mm, no offset
in IP8 (lhcb) and no collision point in IP2 (alice) because beams are separated by 5 sigmas (halo
collision) (+-2.5 sigmas for each beam).

- bv flags should be like in example2.madx : bv=1 for beam 1, bv=-1 for beam 2 (it is the direction 
of rotation)


NOTE: use "as-built" sequence: it is LHC as it is now; the other design sequence is what is foreseen
(for instance there are phase 2 collimators, etc.). but there is very low impact on optics.
also, the beam screen is on each cold magnet (plus some small amount on the sides). warm magnets are
indicated by the letter "W".

NOTE: to get the VSS markers (beginnings and ends of the beam screens near the IPs) you need to use
the file "/afs/cern.ch/eng/lhc/optics/V6.503/aperture/layoutapertures.madx" and some additional definitions
used there (see "Examples_ELaface/check-aperture-injection" and "Examples_ELaface/check-aperture-collision")


##############################
To get the collimators beta functions:
- launch 

madx_dev < thick_asbuilt_inj_ph1collimators.madx

madx_dev < thick_asbuilt_coll3.5TeV_sq1.5m_ph1collimators.madx

madx_dev < thick_asbuilt_coll7TeV_ph1collimators.madx
(squeeze 0.55m in IP1 and IP5 and 10m in the other IPs)


- for e.g. b1 at injection, do:

cp twiss.asbuilt.b1.injection_ph1collimators.thick coll_ph1_beta_450GeV_b1.txt

suppress (with vi) 45 lines at the beginning, 1 line after column headers, and
the * at the beginning of the column headers

awk '{print $1,"\t",$5,"\t",$6}' coll_ph1_beta_450GeV_b1.txt > out

suppress with vi 3 lines with collimators from the other beam (2 TCTVB and 1 TCLIA)

suppress with vi all " (use the command :%s/"//g )

with vi, add "A" to BETX and BETY (gives BETAX and BETAY)

mv out coll_ph1_beta_450GeV_b1.txt

