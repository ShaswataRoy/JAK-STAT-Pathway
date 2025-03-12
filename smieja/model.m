    %         full model with IRF1|STAT1 complexes,
    %         CBP proteins and IRF1|CBP complexes not modelled explicitely
    %         but included in the inertial element - up to 4th order
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function dy=model(t,y,par)

%parameters; %parametrization included in M-file 'parameters.m'
 
 NTB=1;              % NTB=1 = no translation blocking 
                    % NTB=0 - translation blocking (CHX)  
 IA=0;               % IA=1 - inactive STAT1

 block=1;             % if 0 then there is no saturation in phosphorylation; 1 otherwise

 output=0;            % output=1 - concentrations, output=0 - number of molecules

 %original parameters
 kv = par(1);
 v_irf1_transcription = par(2);
 v_stat1_transcription = par(3);
 v_stat2_transcription = par(4);
 v_lmp2_transcription = par(5);
 v_tap1_transcription = par(6);
 ktransl = par(7); 
 constant1 = par(8);
 constant2 = par(9);
 constant3 = par(10);
 constant4 = par(11);
 constant5 = par(12);
 ks1deg = par(13);
 ks1pdeg = par(14);
 ks2deg = par(15);
 ks2pdeg = par(16);
 ki1deg = par(17);
 ki1_indeg = par(18);
 kdegs1t = par(19);
 kdegs2t = par(20);
 kdegi1t = par(21);
 kdeglmp2t = par(22);
 kdegtap1t = par(23);
 kinvs1s1 = par(24); 
 kinvs1s1n = par(25);
 kinvs1s2 = par(26);
 kinvs1s2n = par(27);
 kinvpiass1s1 = par(28);
 kinvs1i1 = par(29);
 ks1s1pdeg = par(30);
 ks1s2pdeg = par(31);
 ks1i1deg = par(32);
 kdegpiass1s1 = par(33);
 ks1s1 = par(34);
 ks1s2 = par(35);
 kpiass1s1 = par(36);
 ks1i1 = par(37);
  kinacti1 = par(39);
 ks1tprod = par(40);
 ks2tprod = par(41);
 lmp2tprod = par(42);
 tap1tprod = par(43);
 es1 = par(44);
 is1 = par(45);
 es2 = par(46);
 is2 = par(47);
 is1s1 = par(48); 
 is1s2 = par(49); 
 ii1 = par(50);
 ei1 = par(51);
 ei1_in = par(52);
 ks1phos = par(53); 
 ks1dephc = par(54);
 ks2phos = par(55); 
 ks2dephc = par(56);
 ks1phos_sat = par(57);
 ks2phos_sat= par(58);
 
 %new parameters
 Ratot = par(59);
 kIFNaRa = par(60); 
 kinvIFNaRa = par(61);  
 kIFNadeg = par(62);
 Rgtot = par(63);
 kIFNgRa = par(64); 
 kinvIFNgRa = par(65);  
 kIFNgdeg = par(66); 
 Rst = par(67);
 Rf = par(68);
 kactPTP = par(69); 
 kinactPTP = par(70); 
 kPTPss = par(71); 
 kinvPTPss = par(72);
 PTPinRa = par(74);
 PTPinRg = par(75);
 ePTP = par(76);
 iPTP = par(77);
 
 
%######################################################################### 

% TR=0 IFN off, TR=1 IFN on
dy=zeros(44,1);

% equation 1 - STAT1 unphosphorylated in cytoplasm
dy(1) = - ks1deg*y(1) - Rf*y(33)/(y(33)+Rst)*ks1phos*y(1)/(1+block*ks1phos_sat*y(1)) - y(35)/(y(35)+Rst)*ks1phos*y(1)/(1+block*ks1phos_sat*y(1)) + ks1dephc*y(3) - is1*y(1) + es1*y(5) + NTB*ktransl*y(11) + 2*kinvs1s1*y(7) + kinvs1s2*y(8) + 2*kinvPTPss*y(40) + kinvPTPss*y(41);

% equation 2 - STAT2 unphosphorylated in cytoplasm
dy(2) = NTB*ktransl*y(24) - ks2deg*y(2) - Rf*y(33)/(y(33)+Rst)*ks2phos*y(2)/(1+block*ks2phos_sat*y(2)) + ks2dephc*y(4) - is2*y(2) + es2*y(6) + kinvs1s2*y(8) + kinvPTPss*y(41);

% equation 3 - STAT1 phosphorylated in cytoplasm
dy(3) = Rf*y(33)/(y(33)+Rst)*ks1phos*y(1)/(1+block*ks1phos_sat*y(1)) + y(35)/(y(35)+Rst)*ks1phos*y(1)/(1+block*ks1phos_sat*y(1)) - ks1pdeg*y(3) - ks1dephc*y(3) - 2*ks1s1*y(3)*y(3) - ks1s2*y(3)*y(4);

% equation 4 - STAT2 phosphorylated in cytoplasm
dy(4) = Rf*y(33)/(y(33)+Rst)*ks2phos*y(2)/(1+block*ks2phos_sat*y(2)) - ks2pdeg*y(4) - ks2dephc*y(4) - ks1s2*y(3)*y(4);  

% equation 5 - STAT1 unphosphorylated in nucleus
dy(5) = is1*kv*y(1) - es1*kv*y(5) - ks1deg*y(5) + kinvs1s2n*y(10) + 2*kinvs1s1n*y(9) - ks1i1*y(18)*y(5) + kinvs1i1*y(25) + 2*kinvPTPss*y(42) + kinvPTPss*y(43);

% equation 6 - STAT2 unphosphorylated in nucleus
dy(6) = is2*kv*y(2) - es2*kv*y(6) - ks2deg*y(6) + kinvs1s2n*y(10) + kinvPTPss*y(43);

% equation 7 - STAT1p|STAT1p complex in cytoplasm; 
dy(7) = ks1s1*y(3)*y(3) - is1s1*y(7) - ks1s1pdeg*y(7) - kinvs1s1*y(7) - kPTPss*y(7)*y(37); 

% equation 8 - STAT1p|STAT2p complex in cytoplasm; 
dy(8) = ks1s2*y(3)*y(4) - is1s2*y(8) - ks1s2pdeg*y(8) - kinvs1s2*y(8) - kPTPss*y(8)*y(37); 

% equation 9 - STAT1p|STAT1p complex in nucleus; 
dy(9) = is1s1*kv*y(7) - ks1s1pdeg*y(9) - kinvs1s1n*y(9) - kPTPss*y(9)*y(39);

% equation 10 - STAT1p|STAT2p complex in nucleus; 
dy(10) = is1s2*kv*y(8) - ks1s2pdeg*y(10) - kinvs1s2*y(10) - kPTPss*y(10)*y(39);

% equation 11 - STAT1 transcript     
dy(11) = NTB*ks1tprod + v_stat1_transcription*y(18) - kdegs1t*y(11); 

% ########################### hipothetic PIAS  - deleted ################################
% equation 12 - active, free PIAS
% equation 13 - PIAS|STAT1|STAT1
% equation 14 - inactive PIAS 
% ########################### hipothetic PIAS ################################

% equation 15 - IRF1 transcript with transription rate linearly dependent on TF concentration
dy(15) = v_irf1_transcription*y(9) - kdegi1t*y(15);

% equation 16 - IRF1 inactive in cytoplasm
dy(16) = - ki1_indeg*y(16) + ei1_in*y(19);

% equation 17 - IRF1 active in cytoplasm
dy(17) = NTB*ktransl*y(15) - ki1deg*y(17) - ii1*y(17) + ei1*y(18);

% equation 18 - IRF1 active in nucleus
dy(18) = kv*ii1*y(17) - kv*ei1*y(18) - ki1deg*y(18) - kinacti1*y(18)  - ks1i1*y(18)*y(5) + kinvs1i1*y(25);

% equation 19 - IRF1 inactive in nucleus
dy(19) = kinacti1*y(18) - ki1_indeg*y(19) - kv*ei1_in*y(19);

% ########################### aux STAT ################################
% equation 20 - inertial element (first - with active IRF1n as an input)
dy(20)= - constant1*y(20) + constant1*y(18);
% equation 21 - inertial element (2nd)
dy(21)= - constant1*y(21) + constant1*y(20);
% equation 22 - inertial element (3rd)
dy(22)= - constant1*y(22) + constant1*y(21);
% equation 23 - inertial element (4th - last one - it is an input for STAT1 transcription)
dy(23)= - constant1*y(23) + constant1*y(22);
% ########################### aux STAT ################################

% equation 24 - STAT2 mRNA
dy(24)= NTB*ks2tprod + NTB*v_stat2_transcription*y(22) - kdegs2t*y(24);

% equation 25 - STAT1|IRF1 complexes in nucleus
dy(25) = ks1i1*y(18)*y(5) - (ks1i1deg + kinvs1i1)*y(25);

% ########################## aux LMP - deleted ##################################
% equation 26 - inertial element for LMP2 transctription (first - with STAT1|IRF1 as an input)
% dy(26)= - constant5*y(26) + constant5*y(25);
% equation 27 - inertial element (2nd)
% dy(27)= - constant5*y(27) + constant5*y(26);
% equation 28 - inertial element (3rd)
% dy(28)= - constant5*y(28) + constant5*y(27);
% equation 29 - inertial element (4th - last one - it is an input for LMP2 transcription)
% dy(29)= - constant5*y(29) + constant5*y(28);
% ########################## aux LMP ##################################

%equation 30 - LMP2 transcript - deleted
% dy(30) = NTB*lmp2tprod + NTB*v_lmp2_transcription*y(25) - kdeglmp2t*y(30);

%equation 31 - TAP1 transcript - inertial elements neglected - time
%constant quite low anyway, as indicated by lmp2 results - deleted
% dy(31) = NTB*tap1tprod + NTB*v_tap1_transcription*y(29) - kdegtap1t*y(31);


% ########################## NEW part ##################################
%equation 32-33 - IFNa/b
dy(32) = -kIFNadeg*y(32) - kIFNaRa*y(32)*(Ratot-y(33))/(1+PTPinRa*y(37)); % extracellular IFNa/b
dy(33) = kIFNaRa*y(32)*(Ratot-y(33))/(1+PTPinRa*y(37)) - kinvIFNaRa*y(33); % active receptor IFNa/b

% %equation 34-35 - IFNgamma
dy(34) = -kIFNgdeg*y(34) - kIFNgRa*y(34)*(Rgtot-y(35))/(1+PTPinRg*y(37)); % extracellular IFNgamma
dy(35) = kIFNgRa*y(34)*(Rgtot-y(35))/(1+PTPinRg*y(37)) - kinvIFNgRa*y(35); % active receptor IFNgamma

% ########################## aux PTP (time delay) ##################################
dy(12) = - constant2*y(12) + constant3*y(33)/(y(33)+Rst) + constant4*y(35)/(y(35)+Rst);
dy(13) = - constant5*y(13) + constant1*y(12);
dy(14) = - constant5*y(14) + constant1*y(13);
% ########################## aux PTP ##################################

%equation 36-37 - inactive/active PTP in cytoplasm
dy(36) = - iPTP*y(36) + ePTP*y(38) - kactPTP*y(36)*y(14) + kinactPTP*y(37); %inactive PTP in cytoplasm
dy(37) = - iPTP*y(37) + ePTP*y(39) + kactPTP*y(36)*y(14) - kinactPTP*y(37) - kPTPss*y(7)*y(37) + kinvPTPss*y(40) - kPTPss*y(8)*y(37) + kinvPTPss*y(41);%active PTP in cytoplasm

%equation 38-39 - inactive/active PTP in nucleus
dy(38) = + kv*iPTP*y(36) - kv*ePTP*y(38); %inactive PTP in nucleus
dy(39) = + kv*iPTP*y(37) - kv*ePTP*y(39) - kPTPss*y(9)*y(39) + kinvPTPss*y(42) - kPTPss*y(10)*y(39) + kinvPTPss*y(43); %active PTP in nucleus

%equation 40-41 - STAT|STAT|PTP complexes in cytoplasm
dy(40) = kPTPss*y(7)*y(37) - kinvPTPss*y(40);  %STAT1p|STAT1p|PTP
dy(41) = kPTPss*y(8)*y(37) - kinvPTPss*y(41);  %STAT1p|STAT2p|PTP

%equation 42-43 - STAT|STAT|PTP complexes in nucleus
dy(42) = kPTPss*y(9)*y(39) - kinvPTPss*y(42);  %STAT1p|STAT1p|PTP
dy(43) = kPTPss*y(10)*y(39) - kinvPTPss*y(43);  %STAT1p|STAT2p|PTP


