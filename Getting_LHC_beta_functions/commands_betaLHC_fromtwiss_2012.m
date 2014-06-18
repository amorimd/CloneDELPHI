%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes some betax and betay from LHC
% test different lattice configurations
% compare with collimator settings file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% with and without exp. magnets and crossing angle : collision 4TeV (b1)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m_noexpmag_nocross.thick');
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m.thick');
max(abs(betax2-betax)./betax)
max(abs(betay2-betay)./betay)
% 1.6e-6 difference at worst

% with and without exp. magnets and crossing angle : collision 4TeV (b2)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b2.coll4tev_sq0.6m_noexpmag_nocross.thick');
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll4tev_sq0.6m.thick');
max(abs(betax2-betax)./betax)
max(abs(betay2-betay)./betay)
% 1.3e-6 difference at worst

% with and without qprime matching : collision 4TeV (b1)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m_qprime.thick');
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m.thick');
max(abs(betax2-betax)./betax)
max(abs(betay2-betay)./betay)
% 2e-8 difference at worst

% with and without qprime matching : collision 4TeV (b2)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b2.coll4tev_sq0.6m_qprime.thick');
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll4tev_sq0.6m.thick');
max(abs(betax2-betax)./betax)
max(abs(betay2-betay)./betay)
% 7e-9 difference at worst

% with and without layout for aperture : collision 4TeV (b1)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m_nolayout.thick');
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% 1.4e-5 difference at worst

% with and without layout for aperture : collision 4TeV (b2)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b2.coll4tev_sq0.6m_nolayout.thick');
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll4tev_sq0.6m.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% 1.4e-5 difference at worst

% 2011 as-built vs 2012 as-built : injection (b1)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('twiss.asbuilt.b1.injection.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.injection.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
max(abs(interp1(s2,betax2,s)-betax)./betax)
max(abs(interp1(s2,betay2,s)-betay)./betay)
% 2e-3 max relative difference

% 2011 as-built vs 2012 as-built : injection (b2)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('twiss.asbuilt.b2.injection.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.injection.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
max(abs(interp1(s2,betax2,s)-betax)./betax)
max(abs(interp1(s2,betay2,s)-betay)./betay)
% 2e-3 max relative difference


% thick vs. thin : collision 7 TeV (b1)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev.thin');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (33% and 81%)
(avbetax2-avbetax)/avbetax
(avbetay2-avbetay)/avbetay
% also big difference on averages (4% and 10%)
figure;plot(s,betax,'xb',s2,betax2,'-og');%axis([s(end)-100 s(end) 0 500])
figure;plot(s,betay,'xr',s2,betay2,'-om');%axis([0 100 0 500])

% thick vs. thin : collision 7 TeV (b2)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev.thin');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (43% and 50%)
(avbetax2-avbetax)/avbetax
(avbetay2-avbetay)/avbetay
% not so small difference on averages (1.6% and 7.5%)
figure;plot(s,betax,'xb',s2,betax2,'-og');%axis([s(end)-100 s(end) 0 500])
figure;plot(s,betay,'xr',s2,betay2,'-om');%axis([0 100 0 500])


% thick: collision with 0.6m squeeze vs. nominal collision (0.55m) (b1)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (488% and 463%)
(avbetax2-avbetax)/avbetax
(avbetay2-avbetay)/avbetay
% small difference on averages (1% and 1%)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');

% thick: collision with 0.6m squeeze vs. nominal collision (0.55m) (b2)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b2.coll4tev_sq0.6m.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (476% and 312%)
(avbetax2-avbetax)/avbetax
(avbetay2-avbetay)/avbetay
% small difference on averages (1% and 1%)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');

% thick: collision with 1.5m squeeze (2011) vs. 0.6m squeeze (2012) (b1)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('twiss.asbuilt.b1.coll3.5tev_sq1.5m.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (621% and 620% !)
(avbetax2-avbetax)/avbetax
(avbetay2-avbetay)/avbetay
% big difference on averages (24% and 24%)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');

% thick: collision with 1.5m squeeze (2011) vs. 0.6m (2012) squeeze (b2)
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('twiss.asbuilt.b2.coll3.5tev_sq1.5m.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll4tev_sq0.6m.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (620% and 621% !)
(avbetax2-avbetax)/avbetax
(avbetay2-avbetay)/avbetay
% big difference on averages (24% and 23%)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');

% conclusion: the only important difference is between thin and thick
% optics, and (of course) between different squeezes at collision.
% The rest (magnets from experiments, crossing angle, chromaticity matching,
% layout for aperture included) doesn't matter.
% The tunes are obtained less accurately in thin wrt the real tunes,
% compared to thick (when integrating 1/beta). For this reason 
% we will use the thick optics.

%% comparison HL-LHC (0.1m) / nominal LHC (0.55m)

figpos=[50 50 1650 1150]; % figure position
gcapos=[0.1 0.1 0.85 0.8]; % axes position

% B1
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('../HLLHC_fromSFartoukh/twiss_b1_betastar10cm_7TeV_3mIP8_10mIP2');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (450% and 681%) (normal)
abs((avbetax2-avbetax)/avbetax)
abs((avbetay2-avbetay)/avbetay)
% also very big difference on averages (70% and 69%) (normal)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');

IP5=13329.28923; % IP5 position is 13329.28923 in both cases
%s3=[-80:0.1:80]+IP5;
%figure;plot(s3-IP5,interp1(s,betax,s3)./interp1(s2,betax2,s3),'b');

fig1=figure;plot(s-IP5,betax,'--b','LineWidth',4);hold on;plot(s2-IP5,betax2,'-r','LineWidth',4);
set(fig1,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('HL-LHC (\beta*=10cm)','Nominal LHC (\beta*=55cm)',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP5 [m]','FontSize',28);ylabel('\beta_x [m]','FontSize',28);
grid on;axis([-80 80 0 3.5e4])

% B2
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('../HLLHC_fromSFartoukh/twiss_b2_betastar10cm_7TeV_3mIP8_10mIP2');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (450% and 780%) (normal)
abs((avbetax2-avbetax)/avbetax)
abs((avbetay2-avbetay)/avbetay)
% also very big difference on averages (70% and 69%) (normal)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');


%% comparison HL-LHC from Stephane Fartoukh (0.1m) / HL-LHC home-made (15cm - baseline)

figpos=[50 50 1650 1150]; % figure position
gcapos=[0.1 0.1 0.85 0.8]; % axes position

% B1
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('../HLLHC_fromSFartoukh/twiss_b1_betastar10cm_7TeV_3mIP8_10mIP2');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('../HLLHC/twiss.hllhc.round.b1.coll7tev.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (450% and 681%) (normal)
abs((avbetax2-avbetax)/avbetax)
abs((avbetay2-avbetay)/avbetay)
% also very big difference on averages (70% and 69%) (normal)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');

IP5=13329.28923; % IP5 position is 13329.28923 in both cases
%s3=[-80:0.1:80]+IP5;
%figure;plot(s3-IP5,interp1(s,betax,s3)./interp1(s2,betax2,s3),'b');

fig1=figure;plot(s-IP5,betax,'--b','LineWidth',4);hold on;plot(s2-IP5,betax2,'-r','LineWidth',4);
set(fig1,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('HL-LHC (\beta*=10cm)','HL-LHC (\beta*=15cm)',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP5 [m]','FontSize',28);ylabel('\beta_x [m]','FontSize',28);
grid on;axis([-80 80 0 3.5e4])

% B2
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('../HLLHC_fromSFartoukh/twiss_b2_betastar10cm_7TeV_3mIP8_10mIP2');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('../HLLHC/twiss.hllhc.round.b2.coll7tev.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (450% and 780%) (normal)
abs((avbetax2-avbetax)/avbetax)
abs((avbetay2-avbetay)/avbetay)
% also very big difference on averages (70% and 69%) (normal)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');


%% comparison round HL-LHC (15cm) / flat HL-LHC (30/7.5cm) / nominal LHC (55cm)

figpos=[50 50 1650 1150]; % figure position
gcapos=[0.1 0.1 0.85 0.8]; % axes position

% B1
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('../HLLHC/twiss.hllhc.round.b1.coll7tev.thick');
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('../HLLHC/twiss.hllhc.flat.b1.coll7tev.thick');
[s3,betax3,betay3,avbetax3,avbetay3]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev.thick');

figure;plot(s,betax,'xb',s2,betax2,'-og',s3,betax3,'-+k');
figure;plot(s,betay,'xr',s2,betay2,'-om',s3,betay3,'-+k');

IP5=13329.28923; % IP5 position is 13329.28923 in both cases

fig1=figure;plot(s-IP5,betax/1e3,'--b','LineWidth',4);hold on;
plot(s2-IP5,betax2/1e3,'-r','LineWidth',4);plot(s3-IP5,betax3/1e3,'-.g','LineWidth',4);
set(fig1,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('HL-LHC (round - 15cm)','HL-LHC (flat - 30cm/7.5cm)','nominal LHC (55cm)',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP5 [m]','FontSize',28);ylabel('\beta_x [km]','FontSize',28);
grid on;axis([-200 200 0 3.5e1]);%axis([0 80 0 3.5e4])

fig2=figure;plot(s-IP5,betay/1e3,'--b','LineWidth',4);hold on;
plot(s2-IP5,betay2/1e3,'-r','LineWidth',4);plot(s3-IP5,betay3/1e3,'-.g','LineWidth',4);
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('HL-LHC (round - 15cm)','HL-LHC (flat - 30cm/7.5cm)','nominal LHC (55cm)',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP5 [m]','FontSize',28);ylabel('\beta_y [km]','FontSize',28);
grid on;axis([-200 200 0 5e1]);

% B2
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('../HLLHC/twiss.hllhc.round.b2.coll7tev.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('../HLLHC/twiss.hllhc.flat.b2.coll7tev.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (450% and 780%) (normal)
abs((avbetax2-avbetax)/avbetax)
abs((avbetay2-avbetay)/avbetay)
% also very big difference on averages (70% and 69%) (normal)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');

%% comparison LHC 2012 / LHC 2015 with changes in IR4 & IR8 (injection)

figpos=[50 50 1650 1150]; % figure position
gcapos=[0.1 0.1 0.85 0.8]; % axes position

% B1
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b1.injection.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.injection.2015.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (24% and 165%) (normal)
abs((avbetax2-avbetax)/avbetax)
abs((avbetay2-avbetay)/avbetay)
% very small difference on averages (<1.5%) (normal)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');

IP4=9997.005016; % IP4 position is 9997.005016 in both cases
IP8=23315.37898; % IP8 position is 23315.37898 in both cases

fig1=figure;plot(s-IP4,betax,'-b','LineWidth',4);hold on;plot(s2-IP4,betax2,'--b','LineWidth',4);
plot(s-IP4,betay,'-r','LineWidth',4);hold on;plot(s2-IP4,betay2,'--r','LineWidth',4);
set(fig1,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('\beta_x, 2012','\beta_x, 2015','\beta_y, 2012','\beta_y, 2015',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP4 [m]','FontSize',28);ylabel('\beta [m]','FontSize',28);
grid on;axis([-300 300 0 600])

fig2=figure;plot(s-IP8,betax,'-b','LineWidth',4);hold on;plot(s2-IP8,betax2,'--b','LineWidth',4);
plot(s-IP8,betay,'-r','LineWidth',4);hold on;plot(s2-IP8,betay2,'--r','LineWidth',4);
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('\beta_x, 2012','\beta_x, 2015','\beta_y, 2012','\beta_y, 2015',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP8 [m]','FontSize',28);ylabel('\beta [m]','FontSize',28);
grid on;axis([-1000 1000 0 300])

% B2
[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b2.injection.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.injection.2015.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
[mx,ix]=max(abs(interp1(s2,betax2,s)-betax)./betax)
[my,iy]=max(abs(interp1(s2,betay2,s)-betay)./betay)
% very big max relative difference (58% and 440%) (normal)
abs((avbetax2-avbetax)/avbetax)
abs((avbetay2-avbetay)/avbetay)
% very small difference on averages (<1%) (normal)
figure;plot(s,betax,'xb',s2,betax2,'-og');
figure;plot(s,betay,'xr',s2,betay2,'-om');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% octupoles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% obtaining the tune coefficients and Q''
% 3.5TeV, squeeze 1.5m (2011)
% thick: collision with 1.5m squeeze (b1)
[a2,qsec2]=octupole('twiss.asbuilt.b1.coll3.5tev_sq1.5m_oct_disp.thick')

% 4TeV, squeeze 0.6m (2011)
% thick: collision with 0.6m squeeze (b1)
[a,qsec]=octupole('2012/twiss.asbuilt.b1.coll4tev_sq0.6m_oct_disp.thick')
[ab2,qsecb2]=octupole('2012/twiss.asbuilt.b2.coll4tev_sq0.6m_oct_disp.thick')
abs(a(:)-ab2(:))./abs(a(:)) % not more than 3%
abs(qsec(:)-qsecb2(:))./abs(qsec(:)) % not more than 3%

abs(a(:)-a2(:))./abs(a(:)) % not more than 1e-7
abs(qsec(:)-qsec2(:))./abs(qsec(:)) % not more than 3e-3

% 450GeV (2012)
% thick (b1)
[ainj,qsecinj]=octupole('2012/twiss.asbuilt.b1.injection_oct_disp.thick')
abs(ainj(:)-a(:))./abs(a(:)) % not more than 2e-7
abs(qsecinj(:)-qsec(:))./abs(qsec(:)) % not more than 3e-3

% 7TeV (2012)
% thick:  collision with 0.55m squeeze (b1)
[a7,qsec7]=octupole('2012/twiss.asbuilt.b1.coll7tev_oct_disp.thick')
abs(a7(:)-a(:))./abs(a(:)) % not more than 1e-7
abs(qsec7(:)-qsec(:))./abs(qsec(:)) % not more than 1e-3

% 4TeV, squeeze 0.6m (2011)
% thick: collision with 0.6m squeeze (b2)
[ab2,qsecb2]=octupole('2012/twiss.asbuilt.b2.coll4tev_sq0.6m_oct_disp.thick')

abs(ab2(:)-a(:))./abs(a(:)) % not more than 3%
abs(qsecb2(:)-qsec(:))./abs(qsec(:)) % not more than 3%

% 450GeV (2012)
% thick (b2)
[ainj,qsecinj]=octupole('2012/twiss.asbuilt.b2.injection_oct_disp.thick')
abs(ainj(:)-a(:))./abs(a(:)) % not more than 3%
abs(qsecinj(:)-qsec(:))./abs(qsec(:)) % not more than 3%

% 7TeV (2012)
% thick:  collision with 0.55m squeeze (b2)
[a7,qsec7]=octupole('2012/twiss.asbuilt.b2.coll7tev_oct_disp.thick')
abs(a7(:)-a(:))./abs(a(:)) % not more than 3%
abs(qsec7(:)-qsec(:))./abs(qsec(:)) % not more than 3%

% conclusion: within 3% accuracy, one can take any beam and any configuration
% (energy and squeeze)

% HL-LHC 7TeV, B1 & B2, 10cm optics from Stephane Fartoukh
[a1,qsec1]=octupole_special('../HLLHC_fromSFartoukh/twiss_b1_betastar10cm_7TeV_3mIP8_10mIP2','2012/twiss.asbuilt.b1.coll7tev_oct_disp.thick')
[a2,qsec2]=octupole_special('../HLLHC_fromSFartoukh/twiss_b2_betastar10cm_7TeV_3mIP8_10mIP2','2012/twiss.asbuilt.b2.coll7tev_oct_disp.thick')
abs(a1(:)-a2(:))./abs(a1(:)) % up to 25% difference
abs(qsec1(:)-qsec2(:))./abs(qsec1(:)) % up to 2% difference

abs(a1(:)-a(:))./abs(a(:)) % big difference ! (up to factor 5.5)
abs(qsec1(:)-qsec(:))./abs(qsec(:)) % big difference ! (up to factor 2.1)
mean(a1(:)./a(:)) % 3.94
mean(qsec1(:)./qsec(:)) % 2
mean(a2(:)./ab2(:)) % 3.95
mean(qsec2(:)./qsecb2(:)) % 1.98

% HL-LHC 7TeV, B1 & B2, 15cm optics (round)
[a1,qsec1]=octupole_special('../HLLHC/twiss.hllhc.round.b1.coll7tev_oct_disp.thick','2012/twiss.asbuilt.b1.coll7tev_oct_disp.thick',14)
[a2,qsec2]=octupole_special('../HLLHC/twiss.hllhc.round.b2.coll7tev_oct_disp.thick','2012/twiss.asbuilt.b2.coll7tev_oct_disp.thick',14)
abs(a1(:)-a2(:))./abs(a1(:)) % up to 16% difference
abs(qsec1(:)-qsec2(:))./abs(qsec1(:)) % up to 2% difference

abs(a1(:)-a(:))./abs(a(:)) % big difference ! (up to factor 2.75)
abs(qsec1(:)-qsec(:))./abs(qsec(:)) % big difference ! (up to factor 1.7)
mean(a1(:)./a(:)) % 2.15
mean(qsec1(:)./qsec(:)) % 1.58
mean(a2(:)./ab2(:)) % 2.15
mean(qsec2(:)./qsecb2(:)) % 1.56

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new function to get all beta functions (and lengths)
% write everything in [twissfilename]_beta_elements.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% injection optics
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all('2012/twiss.asbuilt.b1.injection.thick','2012/twiss.asbuilt.b2.injection.thick');

% injection optics with 2015 changes (IR4 & IR8)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all('2012/twiss.asbuilt.b1.injection.2015.thick','2012/twiss.asbuilt.b2.injection.2015.thick');

% collision with 0.6m squeeze in IP1 and IP5 (3m in IP2 and 3m IP8)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all('2012/twiss.asbuilt.b1.coll4tev_sq0.6m.thick','2012/twiss.asbuilt.b2.coll4tev_sq0.6m.thick');

% collision with 0.55m squeeze in IP1 and IP5 (10m in IP2 and 10m IP8)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all('2012/twiss.asbuilt.b1.coll7tev.thick','2012/twiss.asbuilt.b2.coll7tev.thick');

% collision with 1.5m squeeze in IP1 and IP5 (10m in IP2 and 3m IP8) (2011)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all('twiss.asbuilt.b1.coll3.5tev_sq1.5m.thick','twiss.asbuilt.b2.coll3.5tev_sq1.5m.thick');

% collision with 1m squeeze in IP1 and IP5 (10m in IP2 and 3m IP8) (2011 - end)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all('twiss.asbuilt.b1.coll3.5tev_sq1m.thick','twiss.asbuilt.b2.coll3.5tev_sq1m.thick');

% HL-LHC case (optics from S. Fartoukh, June 2013, 10cm beta*)
% NOTE: in this case the devices are found from LHC twiss file (2012), only
% beta functions are HL-LHC. This is quite a dirty procedure (just for quick preliminary
% model)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all_special('../HLLHC_fromSFartoukh/twiss_b1_betastar10cm_7TeV_3mIP8_10mIP2','../HLLHC_fromSFartoukh/twiss_b2_betastar10cm_7TeV_3mIP8_10mIP2','2012/twiss.asbuilt.b1.coll7tev.thick','2012/twiss.asbuilt.b2.coll7tev.thick');

% HL-LHC case, better (twiss file from S. Fartoukh, June 2013, 10cm beta*)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all_HLLHC('../HLLHC_fromSFartoukh/twiss_b1_betastar10cm_7TeV_3mIP8_10mIP2','../HLLHC_fromSFartoukh/twiss_b2_betastar10cm_7TeV_3mIP8_10mIP2');

% HL-LHC case, better (home-made twiss with round optics, 15cm beta*)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all_HLLHC('../HLLHC/twiss.hllhc.round.b1.coll7tev.thick','../HLLHC/twiss.hllhc.round.b2.coll7tev.thick');

% HL-LHC case, better (home-made twiss with flat optics, 30cm/7.5cm beta*)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all_HLLHC('../HLLHC/twiss.hllhc.flat.b1.coll7tev.thick','../HLLHC/twiss.hllhc.flat.b2.coll7tev.thick');

% HL-LHC case, better (home-made twiss with injection optics)
[names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all_HLLHC('../HLLHC/twiss.hllhc.b1.inj.thick','../HLLHC/twiss.hllhc.b2.inj.thick');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% betas in the triplets (IR1 & 5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LHC post-LS1 case (7TeV, 0.55m squeeze in IP1 and IP5, 10m in IP2 and 10m IP8), B1
[names,betax,betay,len]=meanbeta_triplets_LHC('2012/twiss.asbuilt.b1.coll7tev.thick');

% LHC 2012 case (4TeV, 0.6m squeeze in IP1 and IP5, 3m in IP2 and IP8), B1
[names,betax,betay,len]=meanbeta_triplets_LHC('2012/twiss.asbuilt.b1.coll4tev_sq0.6m.thick');

% HL-LHC case (optics from S. Fartoukh, June 2013, 10cm beta*), B1 & B2
[names,betax,betay,len]=meanbeta_triplets_HLLHC('../HLLHC_fromSFartoukh/twiss_b1_betastar10cm_7TeV_3mIP8_10mIP2');
[names,betax,betay,len]=meanbeta_triplets_HLLHC('../HLLHC_fromSFartoukh/twiss_b2_betastar10cm_7TeV_3mIP8_10mIP2');

%% collimators for nominal LHC (0.55m)

figpos=[50 50 1650 1150]; % figure position
gcapos=[0.1 0.1 0.85 0.8]; % axes position

% B1
[s,betax,betay,avbetax,avbetay,a,b,c,d,mux,muy]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev_ph1collimators.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2,a,b,c,d,mux2,muy2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev.thick');
%[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll4tev_sq0.6m.thick');
%[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b1.injection.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
%figure;plot(s,betax,'xb',s2,betax2,'-b');
%figure;plot(s,betay,'xr',s2,betay2,'-r');

% B1: read info file and compute sigma and impedance scaling factors for the collimators
data=dlmread('2012/twiss.asbuilt.b1.coll7tev_ph1collimators.thick_info','\t',1,1);
angle=data(:,7);betaxdata=data(:,3);betaydata=data(:,4);
sdata=data(:,1);ind=zeros(1,length(sdata));
for i=1:length(sdata)
    ind(i)=find(s==sdata(i));
end
eps=2e-6/7460.52; % emittance for 7TeV
sigma=sqrt(eps*(cos(angle).^2.*betax(ind) + sin(angle).^2.*betay(ind)));
impx=betax(ind)./(sigma.^3);impy=betay(ind)./(sigma.^3);

IP3=6664.7208; % IP3 long. position
IP7=19994.1624; % IP7 long. position

% IR3
fig1=figure;subplot(311);
plot(s2-IP3,betax2,'-b','LineWidth',4);hold on;plot(s-IP3,betax,'xb','LineWidth',4,'MarkerSize',12);
plot(s2-IP3,betay2,'-r','LineWidth',4);hold on;plot(s-IP3,betay,'xr','LineWidth',4,'MarkerSize',12);
set(fig1,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[0.1 0.7 0.85 0.25]);
l=legend('\beta_x','\beta_x collimators','\beta_y','\beta_y collimators',0);set(l,'Fontsize',24);
ylabel('\beta [m]','FontSize',28);
grid on;axis([-240 260 0 400])
subplot(312);
%plot(s2-IP3,mux2,'-b','LineWidth',4);hold on;plot(s-IP3,mux,'xb','LineWidth',4,'MarkerSize',10);
%plot(s2-IP3,muy2,'-r','LineWidth',4);hold on;plot(s-IP3,muy,'xr','LineWidth',4,'MarkerSize',10);
%set(fig1,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
%set(gca,'FontSize',28,'LineWidth',2.5);%,'Position',gcapos);
%l=legend('\mu_x','\mu_x collimators','\mu_y','\mu_y collimators',0);set(l,'Fontsize',24);
%xlabel('Longitudinal position from IP3 [m]','FontSize',28);ylabel('Phase advance / 2\pi','FontSize',28);
%grid on;axis([-250 250 14.5 18])
bar(s(2:end)-1-IP3,diff(mux),1,'b');hold on;
bar(s(2:end)+1-IP3,diff(muy),1,'r');
set(fig1,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[0.1 0.4 0.85 0.25])
l=legend('\Delta \mu_x','\Delta \mu_y',0);set(l,'Fontsize',24);
ylabel('\Delta \mu / 2\pi','FontSize',28);
grid on;axis([-240 260 0 0.5])
subplot(313);
bar(data(:,1)-1-IP3,impx/1e14,1,'b');hold on;
bar(data(:,1)+1-IP3,impy/1e14,1,'r');hold on;
set(fig1,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[0.1 0.08 0.85 0.25]);
l=legend('x','y',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP3 [m]','FontSize',28);
ylabel('Impedance scaling factor [10^{14} m^{-2}]','FontSize',20);
grid on;xlim([-240 260])

% IR7
fig2=figure;subplot(311);
plot(s2-IP7,betax2,'-b','LineWidth',4);hold on;plot(s-IP7,betax,'xb','LineWidth',4,'MarkerSize',10);
plot(s2-IP7,betay2,'-r','LineWidth',4);hold on;plot(s-IP7,betay,'xr','LineWidth',4,'MarkerSize',10);
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[0.1 0.7 0.85 0.25]);
l=legend('\beta_x','\beta_x collimators','\beta_y','\beta_y collimators',0);set(l,'Fontsize',24);
ylabel('\beta [m]','FontSize',28);
grid on;axis([-250 250 0 500])
subplot(312);
bar(s(2:end)-1-IP7,diff(mux),1.,'b');hold on;
bar(s(2:end)+1-IP7,diff(muy),1.,'r');
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[0.1 0.4 0.85 0.25]);
l=legend('\Delta \mu_x','\Delta \mu_y',0);set(l,'Fontsize',24);
ylabel('\Delta \mu / 2\pi','FontSize',28);
grid on;axis([-250 250 0 0.5])
subplot(313);
bar(data(:,1)-1-IP7,impx/1e14,1,'b');hold on;
bar(data(:,1)+1-IP7,impy/1e14,1,'r');hold on;
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[0.1 0.07 0.85 0.25]);
l=legend('x','y',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP7 [m]','FontSize',28);
ylabel('Impedance scaling factor [10^{14} m^{-2}]','FontSize',20);
grid on;xlim([-250 250])

% B2
%[s,betax,betay,avbetax,avbetay]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev_ph1collimators.thick');
%ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];
%[s2,betax2,betay2,avbetax2,avbetay2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev.thick');
%ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];
%figure;plot(s,betax,'xb',s2,betax2,'-b');
%figure;plot(s,betay,'xr',s2,betay2,'-r');


%% test beta-beating in IR7

figpos=[20 50 1800 1150]; % figure position
gcapos=[0.1 0.1 0.85 0.8]; % axes position
xoff=0.07;xl=0.85/2; % offset and width of axes
betamax=800; % for y axis of beta function plots

% B1
[s,betax,betay,avbetax,avbetay,a,b,c,d,mux,muy,dx]=meanbetatwiss('2012/results_IR7_betabeating/twiss.asbuilt.b1.coll7tev_with_alpha_disp_test_betabeatir7_prematch0.2_verynice_beatings.thick');
%[s,betax,betay,avbetax,avbetay,a,b,c,d,mux,muy,dx]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev_with_alpha_disp_test_betabeatir7.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];mux(ind)=[];muy(ind)=[];dx(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2,a,b,c,d,mux2,muy2,dx2]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev_with_alpha_disp.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];mux2(ind2)=[];muy2(ind2)=[];dx2(ind2)=[];
%figure;plot(s,betax,'xb',s2,betax2,'-b');
%figure;plot(s,betay,'xr',s2,betay2,'-r');
% dispersion
%figure;plot(s,dx,'xb',s2,dx2,'-b');
%figure;plot(s,dx-dx2,'-b');
if length(s)~=length(s2)
    dx2=interp1(s2,dx2,s);
    betax2=interp1(s2,betax2,s);
    betay2=interp1(s2,betay2,s);
else
    disp('Length OK');
end
% "dispersion beating"
fig=figure;
%plot(s,abs((dx-dx2)./dx2),'.b','MarkerSize',10);
plot(s,dx-dx2,'-b','LineWidth',3);
set(fig,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
%l=legend('\Delta D_x / D_{x0}',0);set(l,'Fontsize',24);
l=legend('\Delta D_x',0);set(l,'Fontsize',24);
title('B1','FontSize',30);
xlabel('Longitudinal position [m]','FontSize',28);
ylabel('Absolute dispersion beating [m]','FontSize',28);
grid on;
%
fig=figure;
plot(s,dx,'xb',s2,dx2,'-b','LineWidth',3,'MarkerSize',10);
set(fig,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('New optics','Old optics',0);set(l,'Fontsize',24);
title('B1','FontSize',30)
xlabel('Longitudinal position [m]','FontSize',28);
ylabel('Dispersion [m]','FontSize',28);
grid on;

% beta beating everywhere
fig=figure;
plot(s,abs((betax-betax2)./betax2),'.b',s,abs((betay-betay2)./betay2),'.r','MarkerSize',10);
set(fig,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('\Delta \beta_x / \beta_{x0}','\Delta  \beta_y / \beta_{y0}',0);set(l,'Fontsize',24);
title('B1','FontSize',30)
xlabel('Longitudinal position [m]','FontSize',28);
ylabel('Relative \beta beating','FontSize',28);
grid on;

% B1, collimators
[scoll,betaxcoll,betaycoll,avbetax,avbetay,a,b,c,d,muxcoll,muycoll]=meanbetatwiss('2012/results_IR7_betabeating/twiss.asbuilt.b1.coll7tev_ph1collimators_test_betabeatir7_prematch0.2_verynice_beatings.thick');
%[scoll,betaxcoll,betaycoll,avbetax,avbetay,a,b,c,d,muxcoll,muycoll]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev_ph1collimators_test_betabeatir7.thick');
ind=find(diff(scoll)==0);scoll(ind)=[];betaxcoll(ind)=[];betaycoll(ind)=[];
[s2coll,betax2coll,betay2coll,avbetax2,avbetay2,a,b,c,d,mux2coll,muy2coll]=meanbetatwiss('2012/twiss.asbuilt.b1.coll7tev_ph1collimators.thick');
ind2=find(diff(s2coll)==0);s2coll(ind2)=[];betax2coll(ind2)=[];betay2coll(ind2)=[];
% read info file
data=dlmread('2012/twiss.asbuilt.b1.coll7tev_ph1collimatorsIR7.thick_info','\t',1,1);
angle=data(:,7);betaxdata=data(:,3);betaydata=data(:,4);
sdata=data(:,1);len=data(:,2);ind=zeros(1,length(sdata));
for i=1:length(sdata)
    ind(i)=find(scoll==sdata(i));
end
ind2=zeros(1,length(sdata));
for i=1:length(sdata)
    ind2(i)=find(s2coll==sdata(i));
end
eps=2e-6/7460.52; % emittance for 7TeV
% compute sigma and impedance scaling factors for the collimators (new optics)
sigma=sqrt(eps*(cos(angle).^2.*betaxcoll(ind) + sin(angle).^2.*betaycoll(ind)));
impx_b1=len.*betaxcoll(ind).*(cos(angle).^2+sin(angle).^2/2)./(sigma.^3);
impy_b1=len.*betaycoll(ind).*(sin(angle).^2+cos(angle).^2/2)./(sigma.^3);
% compute sigma and impedance scaling factors for the collimators (old optics)
sigma2=sqrt(eps*(cos(angle).^2.*betax2coll(ind2) + sin(angle).^2.*betay2coll(ind2)));
impx2_b1=len.*betax2coll(ind2).*(cos(angle).^2+sin(angle).^2/2)./(sigma2.^3);
impy2_b1=len.*betay2coll(ind2).*(sin(angle).^2+cos(angle).^2/2)./(sigma2.^3);


IP7=19994.1624; % IP7 long. position
% IR7
% old optics
fig2=figure;
subplot(321);
plot(s-IP7,betax2,'-b','LineWidth',4);hold on;plot(s2coll-IP7,betax2coll,'xb','LineWidth',4,'MarkerSize',10);
plot(s-IP7,betay2,'-r','LineWidth',4);hold on;plot(s2coll-IP7,betay2coll,'xr','LineWidth',4,'MarkerSize',10);
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[xoff 0.7 xl 0.25]);
l=legend('\beta_x','\beta_x collimators','\beta_y','\beta_y collimators',0);set(l,'Fontsize',24);
title('Old optics (B1)','Fontsize',30);
ylabel('\beta [m]','FontSize',28);
grid on;axis([-250 250 0 betamax])
%
subplot(323);
bar(s2coll(2:end)-1-IP7,diff(mux2coll),1.,'b');hold on;
bar(s2coll(2:end)+1-IP7,diff(muy2coll),1.,'r');
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[xoff 0.4 xl 0.25]);
l=legend('\Delta \mu_x','\Delta \mu_y',0);set(l,'Fontsize',24);
ylabel('\Delta \mu / 2\pi','FontSize',28);
grid on;axis([-250 250 0 0.5])
%
subplot(325);
bar(sdata-1-IP7,impx2_b1/1e14,1,'b');hold on;
bar(sdata+1-IP7,impy2_b1/1e14,1,'r');hold on;
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[xoff 0.07 xl 0.25]);
l=legend('x','y',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP7 [m]','FontSize',28);
ylabel('Impedance scaling factor [10^{14} m^{-2}]','FontSize',20);
grid on;axis([-250 250 0 1.])
% new optics
subplot(322);
plot(s-IP7,betax,'-b','LineWidth',4);hold on;plot(scoll-IP7,betaxcoll,'xb','LineWidth',4,'MarkerSize',10);
plot(s-IP7,betay,'-r','LineWidth',4);hold on;plot(scoll-IP7,betaycoll,'xr','LineWidth',4,'MarkerSize',10);
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[2*xoff+xl 0.7 xl 0.25]);
title('New optics (B1)','Fontsize',30);
grid on;axis([-250 250 0 betamax])
%
subplot(324);
bar(scoll(2:end)-1-IP7,diff(muxcoll),1.,'b');hold on;
bar(scoll(2:end)+1-IP7,diff(muycoll),1.,'r');
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[2*xoff+xl 0.4 xl 0.25]);
grid on;axis([-250 250 0 0.5])
%
subplot(326);
bar(sdata-1-IP7,impx_b1/1e14,1,'b');hold on;
bar(sdata+1-IP7,impy_b1/1e14,1,'r');hold on;
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[2*xoff+xl 0.07 xl 0.25]);
xlabel('Longitudinal position from IP7 [m]','FontSize',28);
grid on;axis([-250 250 0 1.])


% B2
[s,betax,betay,avbetax,avbetay,a,b,c,d,mux,muy,dx]=meanbetatwiss('2012/results_IR7_betabeating/twiss.asbuilt.b2.coll7tev_with_alpha_disp_test_betabeatir7_prematch0.2_verynice_beatings.thick');
%[s,betax,betay,avbetax,avbetay,a,b,c,d,mux,muy,dx]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev_with_alpha_disp_test_betabeatir7.thick');
ind=find(diff(s)==0);s(ind)=[];betax(ind)=[];betay(ind)=[];mux(ind)=[];muy(ind)=[];dx(ind)=[];
[s2,betax2,betay2,avbetax2,avbetay2,a,b,c,d,mux2,muy2,dx2]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev_with_alpha_disp.thick');
ind2=find(diff(s2)==0);s2(ind2)=[];betax2(ind2)=[];betay2(ind2)=[];mux2(ind2)=[];muy2(ind2)=[];dx2(ind2)=[];
%figure;plot(s,betax,'xb',s2,betax2,'-b');
%figure;plot(s,betay,'xr',s2,betay2,'-r');
% dispersion
%figure;plot(s,dx,'xb',s2,dx2,'-b');
if length(s)~=length(s2)
    dx2=interp1(s2,dx2,s);
    betax2=interp1(s2,betax2,s);
    betay2=interp1(s2,betay2,s);
else
    disp('Length OK');
end
% "dispersion beating"
fig=figure;
%plot(s,abs((dx-dx2)./dx2),'.b','MarkerSize',10);
plot(s,dx-dx2,'-b','LineWidth',3);
set(fig,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
%l=legend('\Delta D_x / D_{x0}',0);set(l,'Fontsize',24);
l=legend('\Delta D_x',0);set(l,'Fontsize',24);
title('B2','FontSize',30);
xlabel('Longitudinal position [m]','FontSize',28);
ylabel('Absolute dispersion beating [m]','FontSize',28);
grid on;
%
fig=figure;
plot(s,dx,'xb',s2,dx2,'-b','LineWidth',3,'MarkerSize',10);
set(fig,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('New optics','Old optics',0);set(l,'Fontsize',24);
title('B2','FontSize',30)
xlabel('Longitudinal position [m]','FontSize',28);
ylabel('Dispersion [m]','FontSize',28);
grid on;

% beta beating everywhere
fig=figure;
plot(s,abs((betax-betax2)./betax2),'.b',s,abs((betay-betay2)./betay2),'.r','MarkerSize',10);
set(fig,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',gcapos);
l=legend('\Delta \beta_x / \beta_{x0}','\Delta  \beta_y / \beta_{y0}',0);set(l,'Fontsize',24);
title('B2','FontSize',30)
xlabel('Longitudinal position [m]','FontSize',28);
ylabel('Relative \beta beating','FontSize',28);
grid on;

% B2, collimators
[scoll,betaxcoll,betaycoll,avbetax,avbetay,a,b,c,d,muxcoll,muycoll]=meanbetatwiss('2012/results_IR7_betabeating/twiss.asbuilt.b2.coll7tev_ph1collimators_test_betabeatir7_prematch0.2_verynice_beatings.thick');
%[scoll,betaxcoll,betaycoll,avbetax,avbetay,a,b,c,d,muxcoll,muycoll]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev_ph1collimators_test_betabeatir7.thick');
ind=find(diff(scoll)==0);scoll(ind)=[];betaxcoll(ind)=[];betaycoll(ind)=[];
[s2coll,betax2coll,betay2coll,avbetax2,avbetay2,a,b,c,d,mux2coll,muy2coll]=meanbetatwiss('2012/twiss.asbuilt.b2.coll7tev_ph1collimators.thick');
ind2=find(diff(s2coll)==0);s2coll(ind2)=[];betax2coll(ind2)=[];betay2coll(ind2)=[];
% read info file
data=dlmread('2012/twiss.asbuilt.b2.coll7tev_ph1collimatorsIR7.thick_info','\t',1,1);
angle=data(:,7);betaxdata=data(:,3);betaydata=data(:,4);
sdata=data(:,1);len=data(:,2);ind=zeros(1,length(sdata));
for i=1:length(sdata)
    ind(i)=find(scoll==sdata(i));
end
ind2=zeros(1,length(sdata));
for i=1:length(sdata)
    ind2(i)=find(s2coll==sdata(i));
end
eps=2e-6/7460.52; % emittance for 7TeV
% compute sigma and impedance scaling factors for the collimators (new optics)
sigma=sqrt(eps*(cos(angle).^2.*betaxcoll(ind) + sin(angle).^2.*betaycoll(ind)));
impx_b2=len.*betaxcoll(ind).*(cos(angle).^2+sin(angle).^2/2)./(sigma.^3);
impy_b2=len.*betaycoll(ind).*(sin(angle).^2+cos(angle).^2/2)./(sigma.^3);
% compute sigma and impedance scaling factors for the collimators (old optics)
sigma2=sqrt(eps*(cos(angle).^2.*betax2coll(ind2) + sin(angle).^2.*betay2coll(ind2)));
impx2_b2=len.*betax2coll(ind2).*(cos(angle).^2+sin(angle).^2/2)./(sigma2.^3);
impy2_b2=len.*betay2coll(ind2).*(sin(angle).^2+cos(angle).^2/2)./(sigma2.^3);


IP7=19994.1624; % IP7 long. position
% IR7
% old optics
fig2=figure;
subplot(321);
plot(s-IP7,betax2,'-b','LineWidth',4);hold on;plot(s2coll-IP7,betax2coll,'xb','LineWidth',4,'MarkerSize',10);
plot(s-IP7,betay2,'-r','LineWidth',4);hold on;plot(s2coll-IP7,betay2coll,'xr','LineWidth',4,'MarkerSize',10);
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[xoff 0.7 xl 0.25]);
l=legend('\beta_x','\beta_x collimators','\beta_y','\beta_y collimators',0);set(l,'Fontsize',24);
title('Old optics (B2)','Fontsize',30);
ylabel('\beta [m]','FontSize',28);
grid on;axis([-250 250 0 betamax])
%
subplot(323);
bar(s2coll(2:end)-1-IP7,diff(mux2coll),1.,'b');hold on;
bar(s2coll(2:end)+1-IP7,diff(muy2coll),1.,'r');
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[xoff 0.4 xl 0.25]);
l=legend('\Delta \mu_x','\Delta \mu_y',0);set(l,'Fontsize',24);
ylabel('\Delta \mu / 2\pi','FontSize',28);
grid on;axis([-250 250 0 0.5])
%
subplot(325);
bar(sdata-1-IP7,impx2_b2/1e14,1,'b');hold on;
bar(sdata+1-IP7,impy2_b2/1e14,1,'r');hold on;
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[xoff 0.07 xl 0.25]);
l=legend('x','y',0);set(l,'Fontsize',24);
xlabel('Longitudinal position from IP7 [m]','FontSize',28);
ylabel('Impedance scaling factor [10^{14} m^{-2}]','FontSize',20);
grid on;axis([-250 250 0 1.])
% new optics
subplot(322);
plot(s-IP7,betax,'-b','LineWidth',4);hold on;plot(scoll-IP7,betaxcoll,'xb','LineWidth',4,'MarkerSize',10);
plot(s-IP7,betay,'-r','LineWidth',4);hold on;plot(scoll-IP7,betaycoll,'xr','LineWidth',4,'MarkerSize',10);
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[2*xoff+xl 0.7 xl 0.25]);
title('New optics (B2)','Fontsize',30);
grid on;axis([-250 250 0 betamax])
%
subplot(324);
bar(scoll(2:end)-1-IP7,diff(muxcoll),1.,'b');hold on;
bar(scoll(2:end)+1-IP7,diff(muycoll),1.,'r');
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[2*xoff+xl 0.4 xl 0.25]);
grid on;axis([-250 250 0 0.5])
%
subplot(326);
bar(sdata-1-IP7,impx_b2/1e14,1,'b');hold on;
bar(sdata+1-IP7,impy_b2/1e14,1,'r');hold on;
set(fig2,'Color',[1 1 1],'Position',figpos,'PaperUnits','points','PaperPosition',figpos);
set(gca,'FontSize',28,'LineWidth',2.5,'Position',[2*xoff+xl 0.07 xl 0.25]);
xlabel('Longitudinal position from IP7 [m]','FontSize',28);
grid on;axis([-250 250 0 1.])

disp(['B1 new/old optics: Zx: ',num2str(sum(impx_b1)/sum(impx2_b1)),', Zy: ',num2str(sum(impy_b1)/sum(impy2_b1))]);
disp(['B2 new/old optics: Zx: ',num2str(sum(impx_b2)/sum(impx2_b2)),', Zy: ',num2str(sum(impy_b2)/sum(impy2_b2))]);
