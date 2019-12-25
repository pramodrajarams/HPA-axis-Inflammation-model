%% This code is  a part of the manuscript titled "Role of enhanced glucocorticoid receptor sensitivity
%%in inflammation in PTSD: Insights from a computational model for circadian-neuroendocrine-immune interaction".
%% The code contains the ordinary differential equations that can simulates 
%the features of HPA-axis and inflammatory pathway for different input conditions.
%% The code is developed by Dr.Pramod R. Somvanshi, Harvard J. A. Paulson School of Engineering and Applied Sciences,
%%Harvard University, Cambridge.

function df= HPA_INFLAM(t,y)


global n  Ki dex nj R Ki1  p r s0  u Q
 
CRH=y(1);
ACTH=y(2);
StARp=y(3);
CORT=y(4);
Dex1 =y(5);
Dex2 =y(6);

CORTp =y(7);
GR_mrna =y(8);
GR_prot =y(9);
GR_cyt =y(10);
GR =y(11);

LPS=y(12);
Phg=y(13);
Phg1= y(14);
TGF=y(15);
TNF=y(16);
IL10=y(17);
IL6=y(18);



%% Specify parameter perturbations here
t1=0;
te=5000;
st=0*((t>t1) & (t<te));
jt=r*((t>t1) & (t<te));
kt=0*r*((t>t1) & (t<te));
n=1*(1+jt);
nj=n;%1*(1+jt);
Ki1=Ki;%1*(1+kt);


%% Specify LPS dose here
lp=p;
It=lp*10e7*((t>712) & (t<712.1));

%% uncomment this to disconnect HPA axis and inflammation to perform IC50 test
%s0=0;


%%Specifiy dexamethasone dose here
dx=0;
dex=dx*10; %Comment this to run IC50 test
d=0;
Dex=dex*((t>696.00+d) & (t<696.25+d));


%%%%%%%%%%%%%%%% HPA-Axis%%%%%%%%%%%%%
K_strs=10;
Ki_GR1=1.2*Ki ;
n1=1*n;
n2=2;
k_crh=0.096;         
V_acth1=1.4634;%0.9756*1.5;
k_acth=0.4312;%1.15*0.1*2.5*1.5;          
k_cort=0.99;%0.44*0.9*2*1.25;     
V_cort=5.67;%0.0945*2*30;
k_dex1=1.95;
V_crh=0.0667*(1/R);
V_acth2=11.2;
Km_tnf3=40*R;
V_starp2=15;
k_starp=0.45;%0.3*1.5;
V_starp1=0.0225;%0.015*1.5;
Kil6=100;
k_dex2=0.25;
V_dex2=1.15;
Ki_GR=50;
V_GRmrna=3.4;
km_GR1=25*Ki;
kon=0.00329;
k_GRprot=0.0572 ;
krt=0.63 ;
kre=0.57;
k_GRmrna=0.1124;
V_GRprot=1.2;
fr=0.49;
tc=0.1401;
ki_tnf=1e3;
n3=2*n;
n4=2;
fac=16.65;%0.5*33.3;

%%Circadian drive from SCN clock
omega=2*pi/24;
cir=2*(1+ cos(omega*(t)));

%%Circadian cdrive from adrenal peripheral clock
cir2=(4-cir)*(Ki_GR/(Ki_GR+GR));

%% Corticotrophic releasing hormone (CRH)
dy(1)=(K_strs*(1+st)*(Ki_GR1^n1/(Ki_GR1^n1+GR^n1))*(cir)*(1+V_crh*TNF)-k_crh*CRH);

%% Adrenocorticotropic hormone (ACTH)
dy(2)= (V_acth1*CRH*(Ki_GR1^n1/(Ki_GR1^n1+GR^n1))*(1+V_acth2*((TNF)^n2/(Km_tnf3^n2+ TNF^n2)))-k_acth*ACTH) ;

%% StAR proetin
dy(3)=s0*(V_starp1*(ACTH*cir2)*(1+V_starp2*(TNF^n2/(Km_tnf3^n2+TNF^n2))*(Kil6/(Kil6+IL6)))  -k_starp*StARp);

%% Cortisol
dy(4)= (V_cort*StARp-k_cort*CORT);

%% Dexamethasone kinetics
dy(5)= Dex-k_dex1*Dex1;
dy(6)= V_dex2*Dex1-k_dex2*Dex2;

%% Peripheral cortisol
dy(7)=(1/(tc))*(CORT+(fac*Dex2)-CORTp);

%% GR mRNA
dy(8) = V_GRmrna*(1 - (GR^n3/((km_GR1)^n3+GR^n3)))- k_GRmrna*GR_mrna; 

%% GR protein
dy(9) = V_GRprot*GR_mrna + fr*kre*GR- kon*(CORTp)*GR_prot - k_GRprot*GR_prot;

%% GR-cortisol complex in cytosol
dy(10)= kon*(CORTp)*GR_prot- krt*GR_cyt*(ki_tnf^n4/(TNF^n4+ki_tnf^n4)); 

%% Nuclear GR
dy(11)= krt*GR_cyt*(ki_tnf^n4/(TNF^n4+ki_tnf^n4)) -kre*GR;


%%%%% Inflammatory pathway%%%%%%%%%

%%%%%%%Validation parametset 2 for 0.4 ng/kg @ 2 PM%%%%%%%%%%%%%

k_lps=2.7e-5;%1.35*1e-7*100*2;%0
V_phg1= 4.9956*1e7;
V_phg2=12.949;
Km_tnf1=1693.9509;
Ki_tgf1=0.00721;
Ki_IL10=7.384;%147.68*0.05;
k_phg=1.439;
V_tgf1=0.15625*1e-8;
k_tgf=0.0635;%0.03177*2;

V_tnf1=25.5194;
Km_phg1=412500;%0.075*550*1e4;
Ki_tgf2=0.143;%0.15893*0.9;
V_tnf2=106542;% 3.0*3.5514*1e4;
Km_tnf2=123.96;%0.08*1.5495*1e3;
k_tnf=1.25;
Km_il61=80;

Km_phg2=161012;% 8.0506*1e7*0.2*1e-2;
k_il10=1.6;%200*0.008;
V_il102=2.1938e3;%43875*0.05;
Km_tgf3=0.76;% 0.38*2;
V_il101=1.3374e3;%2.67480*1e6*50.0e-5;
Ki_il102=23.636;%1.1818*20*1;
V_il62=5.0e5;
Km_il6tnf=339.164;%4.8452*70;
Km_phg3=11e6;%110*10e4;
V_il61=0.55e5;
k_il6=1.625;%0.5*3.25;
tp=1.5;
tnf_b=1.25;
il6_b=1.5;
V_tgf2=0.5;
Km_GR2=500*Ki1;
n5=1*nj;
n6=4;
n7=2;
n8=2;
n9=6;
n10=2;

%% LPS
dy(12)= 1e-7*(1+It)-k_lps*LPS*Phg;

%% Phagocytes

dy(13)= (V_phg1*((1+ (V_phg2* TNF/(Km_tnf1+TNF)))* (Ki_tgf1/(Ki_tgf1+ TGF))* (Ki_IL10/(Ki_IL10+IL10)))*LPS- k_phg*Phg);

dy(14)=(1.0/(tp))*(Phg-Phg1);

%% Transformig grwth factor (TGF)
dy(15)=  (V_tgf1* Phg + (V_tgf2*GR^n5/(Km_GR2^n5+GR^n5))- k_tgf*TGF);

%% Tumour necrosis factor (TNF)
dy(16)= tnf_b+(Phg/(Km_phg1+Phg))*(V_tnf1+V_tnf2*(TNF/(Km_tnf2+TNF)))*(Ki_tgf2^n6/(Ki_tgf2^n6+TGF^n6))*(1-(IL6^n7/(Km_il61^n7+IL6^n7)))-k_tnf*(TNF);

%% Interleukin 10 (IL10)
dy(17)=((V_il101*Phg1^n8/(Phg1^n8+Km_phg2^n8))+ (V_il102* TGF^n9/(Km_tgf3^n9+ TGF^n9)) - k_il10*IL10);

%% Interleukin 6 (IL6)
dy(18)=il6_b+(Phg1/(Km_phg3+Phg1))*(V_il61+V_il62*(TNF+IL6)^n10/(Km_il6tnf^n10+(TNF+IL6)^n10))*(Ki_il102/(Ki_il102+(IL10)))-k_il6*IL6;



df=dy';

end

