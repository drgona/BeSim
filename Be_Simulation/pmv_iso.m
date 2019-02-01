function [PMVval, TCL, TS, HC] = pmv_iso(TA, TR)

CLO = 0.5;
% TA = 22;
% TR = 28;
MET = 1.2;
VEL = 0.1;
RH = 60;
FNPS = exp(16.6536 - 4030.183 / (TA + 235));
PA = RH * 10 * FNPS;
ICL = 0.155 * CLO;
M = MET * 58.15;
HC = 0;

if ICL < 0.078
    FCL = 1 + 1.29 * ICL;
else
    FCL = 1.05 + 0.645 * ICL;
end
    
HCF = 12.1 * VEL ^ 0.5;
TAA = TA + 273;
TRA = TR + 273;

TCLA = TAA + (35.5 - TA) / (3.5 * (6.45 * ICL + 0.1));
P1 = ICL * FCL;
P2 = P1 * 3.96;
P3 = P1 * 100;
P4 = P1 * TAA;
P5 = 308.7 - 0.028 * M + P2 * (TRA / 100) ^ 4;
XN = TCLA / 100;
XF = TCLA / 50;
%  XF = XN
N = 0;
EPS = 0.00015;
while abs(XN - XF) > EPS
    XF = (XF + XN) / 2;
    HCF = 12.1 * VEL ^ 0.5;
    HCN = 2.38 * abs(100 * XF - TAA) ^ 0.25;

    if HCF > HCN
        HC = HCF;
    else
        HC = HCN;
    end
    
    XN = (P5 + P4 * HC - P2 * (XF ^ 4)) / (100 + P3 * HC);
    N = N + 1;
end

TCL = 100 * XN - 273;
% 
% 
%  Skinn diff loss
% 
HL1 = 3.05 * 0.001 * (5733 - 6.99 * M - PA);
% 
% Sweat loss
% 
if M > 58.15
    HL2 = 0.42 * (M - 58.15);
else
    HL2 = 0;
end
% 
% Latent respiration loss
% 
HL3 = 1.7 * 0.00001 * M * (5867 - PA);
% 
%  Dry respiration loss
% 
HL4 = 0.0014 * M * (34 - TA);
% 
%  Radiation loss
% 
HL5 = 3.96 * FCL * (XN ^ 4 - (TRA / 100) ^ 4);
% 
%  Convection loss
% 
HL6 = FCL * HC * (TCL - TA);
% 
%  Thermal sensation to skin tran coef
% 
TS = 0.303 * exp(-0.036 * M) + 0.028;
PMVval = TS * (M - HL1 - HL2 - HL3 - HL4 - HL5 - HL6);
end