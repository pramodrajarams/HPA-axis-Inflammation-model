
%% This code is  a part of the manuscript titled "Role of enhanced glucocorticoid receptor sensitivity
%%in inflammation in PTSD: Insights from a computational model for circadian-neuroendocrine-immune interaction".
%% The code contains the ordinary differential equations that can simulates 
%the features of HPA-axis and inflammatory pathway for different input conditions.
%% The code is developed by Dr.Pramod R. Somvanshi, Harvard J. A. Paulson School of Engineering and Applied Sciences,
%%Harvard University, Cambridge.


clc
clear all

global  n  Ki dex  R p r s0 p1 Q


p=0.0;%%%LPS injection
r=1.0;%% Change this to vary GR sensitivity
p1=0;%%%% Change this to 1 to introduce LPS for IC50 test
s0=1;%%% Change this to 0 to disconnect HPA with Inflmmation for IC50 test
%u=1;%0.1% make this 0.1 for simulationg IC50 else 1

%% Varying dexamethasone
p_span1=1;%[0,logspace(-2,2.5,20)];%[0:0.01:0.2];%%[5 10 50 100];%[5 10 20 50];% 1.5 2 2.5 3 3.5 4];10;%
%% Varying sensitivity
p_span2=1;%[0.25:0.25:2];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];
%% Varying inhibitory constant
p_span3=1;%[0.1:0.1:1.5];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];

%% Varying cytokine effect on HPA axis
p_span4=1;%[0.1:0.1:1.5];%[0.25:0.25:2];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];

%% Varying GR negative feedback
p_span5=1;%[0.1:0.1:2];%;%[0.5:0.5:4];%[0.1 0.25 0.5 0.75 1 1.5 2 2.5 3 3.5 4];


%%% The data are extarcted from Grigolet et al., 2010, Clodie et al., 2008., Copeland et al., 2005,
%%% Lauw et al., 2005, Wegner et al., 2017

%% Immune response for 0.4 ng/kg LPS%%%% 
%%The data is extracted from Grigoleit et al.,2010.
T1=[0, 1,1.5 ,2,3,4,6];
il6x=[3.3,13.2,63.7,132,83.5,30.8,7.7];
tnfx=[9,72,89,51,16.5,12.4,10];
tnfs=[2.33,7.7,10.7,6.1,6.5,3.7,2.33];
il10x=[8.5,17.9,30.8,36.8,34.6,15.4,10.7];

T2=[0,1.5,3,6];
cortx=[10,11,19,6];
Acthx=[8.0, 11.0,25.0,5.0];
T3=[0,1.5,2,3,4,5,6];
cortx1=[8,10,16,19,18,15,12];

% %%%Immune response for 2 ng/kg LPS%%%% 
% T2=[0, 1,1.5 ,2,3,4,6];
% il6x1=[6.86, 96,528,1010,761,322,20.6];
% tnfx1=[3.32,267,445,302,109,26.1,20];
% il10x1=[11.9,42.9,81,176,260,114,33.3];
% cortx2=[10.1,9.08,15,17,20,19,13];
% Acthx2=[32,31,40,75,74,35,22];



for count3=1:1:length(p_span3)
     Ki=p_span3(count3);
     
for count2=1:1:length(p_span2)
n=p_span2(count2);

 for count1=1:1:length(p_span1)
 dex=p_span1(count1);
 
 for count4=1:1:length(p_span4)
     R=p_span4(count4);
     
for count5=1:1:length(p_span5)
     Q=p_span5(count5);
 


x=[];
tspan=[0 1500];
options = odeset('RelTol',1e-10,'AbsTol',1.e-10);

%% Initial conditions
y0=[17.4896364862443,4.82352867112865,0.794555311165795,10.8943173854454,0,0,11.1479571299593,20.1160322148203,342.786973832252,24.8798181706611,31.3302385380763,0.000129677801623061,5897.54196590606,25,0.446309929849742,1.00160367257638,22.8093965060946,1.09564141359877];

[t,y] = ode15s(@HPA_INFLAM,tspan,y0,options);


CRH=y(:,1);
ACTH=y(:,2);
StARp=y(:,3);
CORT=y(:,4);
Dex1 =y(:,5);
Dex2 =y(:,6);

CORTp =y(:,7);
GR_mrna =y(:,8);
GR_prot =y(:,9);
GR_cyt =y(:,10);
GR =y(:,11);

LPS=y(:,12);
Phg=y(:,13);
Phg1= y(:,14);
TGF=y(:,15);
TNF=y(:,16);
IL10=y(:,17);
IL6=y(:,18);



V_tgf2=0.5;
Km_GR2=500*Ki;
V_acth2=11.2;
Km_tnf3=40;

V_starp2=15;
n5=1*(1+r);
n3=2*(1+r);
Ki_GR1=1.2;
V_acth1=0.9756*1.5;
n2=2;
Km_GR1=25*Ki;
nkm_Frn1=11*Ki;
vp_rp=0.57*0.49;
Kil6s=100;
n1=1*(1+r);

Flux1=(V_tgf2.*GR.^n5./((Km_GR2).^n5+GR.^n5));
Flux2=(V_acth2.*((TNF).^n2./((Km_tnf3).^n2+ (TNF).^n2)));
Fluxdiff=(V_acth1.*CRH.*((Ki_GR1.^n1./(Ki_GR1.^n1+GR.^n1)))).*(1+V_acth2.*((TNF).^n2./((Km_tnf3).^n2+ (TNF).^n2)));
Flux3=(V_starp2.*(TNF.^n2./(Km_tnf3.^n2+TNF.^n2)).*(Kil6s./(Kil6s+IL6)));
Flux4=(V_acth1.*((Ki_GR1.^n1./(Ki_GR1.^n1+GR.^n1))));
Flux5=(GR.^n3./((Km_GR1).^n3+GR.^n3));
Flux6=vp_rp.*GR;

%% Uncomment thsi to plot figure 7(A)

TGF_ss(count3,count4,count1)=mean(TGF((length(find(t<200))):end));%TGF(end);
TNF_ss(count3,count4,count1)=mean(TNF((length(find(t<200))):end));%TNF(end);
IL6_ss(count3,count4,count1)=mean(IL6((length(find(t<200))):end));%IL6(end);
IL10_ss(count3,count4,count1)=mean(IL10((length(find(t<200))):end));%IL10(end);
FRnr_ss(count3,count4,count1)=mean(GR((length(find(t<200))):end));%GR(end);
Cort_ss(count3,count4,count1)=mean(CORT((length(find(t<200))):end));%CORT(end);
FRr_ss(count3,count4,count1)=mean(GR_prot((length(find(t<200))):end));%GR_prot(end);
ACTH_ss(count3,count4,count1)=mean(ACTH((length(find(t<200))):end));%GR_prot(end);

%% Uncomment thsi to plot figure 7(B)
% TGF_ss(count5,count4,count1)=mean(TGF((length(find(t<200))):end));%TGF(end);
% TNF_ss(count5,count4,count1)=mean(TNF((length(find(t<200))):end));%TNF(end);
% IL6_ss(count5,count4,count1)=mean(IL6((length(find(t<200))):end));%IL6(end);
% IL10_ss(count5,count4,count1)=mean(IL10((length(find(t<200))):end));%IL10(end);
% FRnr_ss(count5,count4,count1)=mean(GR((length(find(t<200))):end));%GR(end);
% Cort_ss(count5,count4,count1)=mean(CORT((length(find(t<200))):end));%CORT(end);
% FRr_ss(count5,count4,count1)=mean(GR_prot((length(find(t<200))):end));%GR_prot(end);
% ACTH_ss(count5,count4,count1)=mean(ACTH((length(find(t<200))):end));%GR_prot(end);

TNF50(count1)=(1-((TNF(end)-min(TNF((length(find(t<600))):end)))./TNF(end))).*100;


 end
end
end
 end
end

%% Estimate IC50
IC50=p_span1(length(find(TNF50>max(TNF50)*0.49)));
%% Estimate percent cortisol suppression
Corsup= ((CORT(length(find(t<681.7)))-CORT(length(find(t<705.7))))./CORT(length(find(t<681.7))))*100;


%% Plotting dynamic profiles of the model variables 
% figure(1)
% subplot(3,4,1)
% plot(t,CRH); hold on;ylabel ('CRH');
% subplot(342)
% plot(t,ACTH); hold on;ylabel ('ACTH');
% subplot(3,4,3)
% plot(t,CORT); hold on;ylabel ('CORT');
% subplot(3,4,4)
% plot(t,StARp); hold on;ylabel ('Starp');
% subplot(3,4,5)
% plot(t,GR); hold on;ylabel ('GR');
% subplot(3,4,6)
% plot(t,GR_prot); hold on;ylabel ('GR_prot');
% subplot(347)
% plot(t,IL6); hold on;ylabel ('IL6');
% subplot(3,4,8)
% plot(t,TNF); hold on;ylabel ('TNF');
% subplot(3,4,9)
% %plot(t,Frb); hold on;ylabel ('C1');
% plot(t,Phg); hold on;ylabel ('Phg');
% %plot(t,cir2); hold on;ylabel ('C1');
% subplot(3,4,10)
% plot(t,TGF); hold on;ylabel ('TGF');
% subplot(3,4,11)
% plot(t,IL10); hold on;ylabel ('IL10');
% subplot(3,4,12)
% plot(t,Dex2); hold on;ylabel ('Dex');


%% Plot Figure 2(B) from the manuscript: Simulation results for validation of the 0.4 ng/kg LPS dose. 
%% Set p=0.4 in the function file

%  figure(2)
%  t1=712;
% tx=((length(find(t<711))):(length(find(t<720))));
% subplot(321)
% plot(t(tx),(TNF(tx)),'LineWidth', 2); hold on;ylabel ('TNF (pg/ml)');  xlabel ('Time (Hrs)');
% scatter(T1+t1,tnfx,'LineWidth', 2); hold on;
% xlim([711,719])
%  xticks([711 712 713 714 715 716 717 718 719])
%  xticklabels({'-1','0','1','2','3','4','5','6','7','8',})
% subplot(322)
% plot(t(tx),(IL6(tx)),'LineWidth', 2); hold on;ylabel ('IL6 (pg/ml)');  xlabel ('Time (Hrs)');
% scatter(T1+t1,il6x,'LineWidth', 2); hold on;
% xlim([711,719])
%  xticks([711 712 713 714 715 716 717 718 719])
%  xticklabels({'-1','0','1','2','3','4','5','6','7','8',})
% subplot(323)
% plot(t(tx),(IL10(tx)),'LineWidth', 2); hold on;ylabel ('IL10 (pg/ml)');  xlabel ('Time (Hrs)');
% scatter(T1+t1,il10x,'LineWidth', 2); hold on;
% xlim([711,719])
%  xticks([711 712 713 714 715 716 717 718 719])
%  xticklabels({'-1','0','1','2','3','4','5','6','7','8',})
% subplot(324)
% plot(t(tx),(CORT(tx)),'LineWidth', 2); hold on;ylabel ('Cortisol (mcg/dl)');  xlabel ('Time (Hrs)');
% scatter(T2+t1,cortx,'LineWidth', 2); hold on;
% xlim([711,719])
%  xticks([711 712 713 714 715 716 717 718 719])
%  xticklabels({'-1','0','1','2','3','4','5','6','7','8',})
% subplot(325)
% plot(t(tx),(ACTH(tx)),'LineWidth', 2); hold on;ylabel ('ACTH (mcg/dl)');  xlabel ('Time (Hrs)');
% scatter(T2+t1,Acthx,'LineWidth', 2); hold on;
% xlim([711,719])
%  xticks([711 712 713 714 715 716 717 718 719])
%  xticklabels({'-1','0','1','2','3','4','5','6','7','8',})
% subplot(326)
% plot(t(tx),(LPS(tx)),'LineWidth', 2); hold on;ylabel ('LPS'); 
% xlim([711,719])
% xticks([711 712 713 714 715 716 717 718 719])
%  xticklabels({'-1','0','1','2','3','4','5','6','7','8',})

%% Plotting Figure 3(A)from the manuscript
%% Set dx=0.5 at the time equvivalant of 11pm.
%% Vary r in function file to vary sensitivity.

% figure(3)
% tx=((length(find(t<670))):(length(find(t<725))));
% subplot(231)
% plot(t(tx),(CORT(tx))./10,'LineWidth', 2); hold on;ylabel ('Cortisol (FC)'); 
%  xlim([674,725])
%   xticks([674 680 686 692 698 704 710 716 722])
%   xticklabels({'12AM','6AM','12PM','6PM','12AM','6AM','12PM','6PM','12AM'})
% subplot(232)
% plot(t(tx),(IL6(tx)),'LineWidth', 2); hold on;ylabel ('IL6 (FC)'); 
% xlim([674,725])
% xticks([674  686  698  710  722])
% xticklabels({'12AM','12PM','12AM','12PM','12AM'})
% subplot(233)
% plot(t(tx),(TNF(tx)),'LineWidth', 2); hold on;ylabel ('TNF (FC)'); 
% xlim([674,725])
% xticks([674  686  698  710  722])
% xticklabels({'12AM','12PM','12AM','12PM','12AM'})
% subplot(234)
% plot(t(tx),(GR_prot(tx)./375),'LineWidth', 2); hold on;ylabel ('Free GR (FC)'); 
% xlim([674,725])
% xticks([674  686  698  710  722])
% xticklabels({'12AM','12PM','12AM','12PM','12AM'})
% subplot(235)
% plot(t(tx),(GR(tx)./20),'LineWidth', 2); hold on;ylabel ('Nuclear GR (FC)'); 
% xlim([674,725])
% xticks([674  686  698  710  722])
% xticklabels({'12AM','12PM','12AM','12PM','12AM'})
% subplot(236)
% plot(t(tx),100-((GR_prot(tx)-(GR(tx)))./(GR_prot(tx)))*100,'LineWidth', 2); hold on;ylabel ('% Nuclear GR (FC)'); 
% xlim([674,725])
% xticks([674  686  698  710  722])
% xticklabels({'12AM','12PM','12AM','12PM','12AM'})


%% Plotting Figure 3 (B) from the manuscript
%% Vary p_span 1 to vary dexamethasone input.
%% Vary kon to simulate change in GR affinity
%% Vary v_prot to simulate chnage in GR expression
%% Vary n5 to simulate changes in GR senstitivity
%% Vary Km_GR2 to simulate changes in GR anti-inflammatory threshold

% figure(4)
% subplot(221)
% semilogx((p_span1(1:end)./3.925e-4)*1e-9,TNF50,'LineWidth', 2);hold on;ylabel ('% TNF suppression'); xlabel ( 'Dexamethasone (AU)');
% subplot(222)
% semilogx((p_span1(1:end)./3.925e-4)*1e-9,IL650,'LineWidth', 2);hold on;
% subplot(223)
% semilogx((p_span1(1:end)./3.925e-4)*1e-9,FRnri_ss,'LineWidth', 2);hold on;
% subplot(224)
% semilogx((p_span1(1:end)./3.925e-4)*1e-9,FR_ri_ss,'LineWidth', 2);hold on;
% 
% Sens=abs((TNF50(9)-TNF50(11))/(p_span1(9)-p_span1(11)))
  

%% Plotting Figure 4 and supplementary figure S1 from the manuscript: Circadian profiles of the model variables for varying parameters  
%% Vary r and Ki in the function file

% figure(5)
% tx=((length(find(t<675))):(length(find(t<775))));
% subplot(321)
% plot(t(tx),(CORT(tx))./10,'LineWidth', 2); hold on;ylabel ('Cortisol'); 
% subplot(3,2,2)
% plot(t(tx),(ACTH(tx))./10,'LineWidth', 2); hold on;ylabel ('ACTH');
% subplot(3,2,3)
% plot(t(tx),(GR(tx))./20,'LineWidth', 2); hold on;ylabel ('Nuclear GR');
% subplot(3,2,4)
% plot(t(tx),(GR_prot(tx))./370,'LineWidth', 2); hold on;ylabel ('Free GR');
% subplot(325)
% plot(t(tx),IL6(tx),'LineWidth', 2); hold on;ylabel ('IL6');
% subplot(326)
% plot(t(tx),(TNF(tx)),'LineWidth', 2); hold on;ylabel ('TNF');

%% Plotting Figure 6 (B) from the manuscript
%% Vary r in function file and record the feedback effects as given above

% figure(6)
% subplot(231)
% plot(p_span2,Flx1_ss./(max(Flx1_ss)),'LineWidth', 2); hold on;ylabel ('GR Anti-inflamatory effect');
% plot(p_span2,Flx7_ss,'LineWidth', 2); hold on;ylabel ('GR-protein positive feedback');
% 
% subplot(232)
% plot(p_span2,1-(Flx4_ss./(max(Flx4_ss))),'LineWidth', 2); hold on;ylabel ('GR HPA inhibition effect');
% subplot(233)
% plot(p_span2,Flx3_ss./(max(Flx3_ss)),'LineWidth', 2); hold on;ylabel ('Cytokine Cortisol effect');
% subplot(234)
% plot(p_span2,Flx2_ss./(max(Flx2_ss)),'LineWidth', 2); hold on;ylabel ('Cytokine ACTH effect');
% subplot(235)
% plot(p_span2,(Flx5_ss./(max(Flx5_ss))),'LineWidth', 2); hold on;ylabel ('GR-mRNA negative feedback');
% subplot(236)
% plot(p_span2,Flx7_ss,'LineWidth', 2); hold on;ylabel ('GR-protein positive feedback');

%% Plotting Figure 6(C) from the manuscript
%% Vary r in function file and record the difference in feedback effects

%%This plots the limit cycles for the difference in feedback effect: 
figure(7)

tx=((length(find(t<700))):(length(find(t<800)))); 

plot(Flux1(tx)./0.035,(Flux3(tx)./max(Flux3(tx)))-(1-Flux4(tx)./max(Flux4(tx))),'LineWidth', 2); hold on; xlabel ('GR anti-inflammatory feedback'); ylabel ('Diffrence in Pfd and Nfd');  hold on; 


%save GCR_flux
%% Plotting Figure  7 ) from the manuscript
% outer_k=[1 2 3 4];
% 
% for kn=1:1:length(outer_k) 
%   
%       
%                 for cnt1=1:1:length(p_span2)
% 
%                                 for cnt2=1:1:length(p_span3)
%                                     
%                                     
% TNFc1(cnt2)=TNF_ss(cnt2,cnt1,outer_k(kn));
% FRnrc1(cnt2)=FRnr_ss(cnt2,cnt1,outer_k(kn));
% Cortc1(cnt2)=Cort_ss(cnt2,cnt1,outer_k(kn));
%                                 end
%                                 
%     TNFc2(cnt1,:)=TNFc1;
%     FRnrc2(cnt1,:)=FRnrc1;
%      Cortc2(cnt1,:)=Cortc1;    
%                 end
% 
%          figure(1),title('Sensitivity matrix');hold on
%         subplot(2,2,kn);surfc(p_span3,p_span2,(FRnrc2)./8.5); xlabel ('Inhibition threshold'); ylabel ('Sensitivity'); 
%         colormap hot;view([0 90]); hold on;colorbar;
%         axis([0 5 0 5 ]) 
% end
%                                     


%% Plotting Figure 7 (A) and (B) from manuscript. Need to change the X/Y lables accordingly

%% Vary p_span2 and p_span3 in for loop for figure 7A
%%  Vary p_span3 and p_span4 in for loop for figure B

%   figure(8)

%  subplot(231)
%  s1=surfc(p_span4,p_span3,Cort_ss/10,'FaceAlpha',1);zlabel ('Cortisol');ylabel ('Q'); xlabel ('R');  hold on; 
% % colormap hsv
% subplot(232)
%  surfc(p_span4,p_span3,ACTH_ss/20,'FaceAlpha',1);zlabel ('ACTH'); %ylabel ('Inhibition threshold'); xlabel ('Sensitivity');  hold on;
%  %colormap hsv
%  subplot(233)
%  surfc(p_span4,p_span3,FRr_ss/240,'FaceAlpha',1);zlabel ('Free GR receptor');% ylabel ('Inhibition threshold'); xlabel ('Sensitivity');  hold on;
%  subplot(234)
% surfc(p_span4,p_span3,FRnr_ss/20,'FaceAlpha',1);  zlabel ('Nuclear GR receptor');hold on;
% %colormap parul
% subplot(235)
%  surfc(p_span4,p_span3,TNF_ss,'FaceAlpha',1);zlabel ('TNF');% ylabel ('Inhibition threshold'); xlabel ('Sensitivity');  hold on; 
%  %colormap hsv
% subplot(236)
%  surfc(p_span4,p_span3,IL6_ss,'FaceAlpha',1);zlabel ('IL6'); %ylabel ('Inhibition threshold'); xlabel ('Sensitivity');  hold on;

%% Plotting the feedback strenghths as a function of other feedback effects
% figure(9)
%tx=((length(find(t<700))):(length(find(t<800)))); 
%  subplot(2,2,1)
%  plot(Flux1(tx),(Flux3(tx)),'LineWidth', 2); hold on; xlabel ('GR anti-inflammatory feedback'); ylabel ('Cytokine ACTH-positive feedback');  hold on; 
% subplot(2,2,2)
%  plot(Flux1(tx),(Flux3(tx)./max(Flux3(tx)))-(1-Flux4(tx)./max(Flux4(tx))),'LineWidth', 2); hold on; xlabel ('GR anti-inflammatory feedback'); ylabel ('Diffrence in Pfd and Nfd');  hold on; 
% subplot(2,2,3)
% plot((1-Flux4(tx)),(Flux2(tx)),'LineWidth', 2); hold on; xlabel ('GR HPA-negative feedback'); ylabel ('Cytokine ACTH-positive feedback');  hold on; 
%  subplot(2,2,4)
%  plot(GR(tx),TNF(tx),'LineWidth', 2); hold on; xlabel ('Nuclear GR'); ylabel ('TNF');  hold on; 


