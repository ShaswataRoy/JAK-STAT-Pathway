%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%   In this M-file one can introduce new parameters into the model,     %
%   also diffrent initial conditions and times characterictic for       %
%   simulation may be chosen.                                           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%############### Parametrization for system of ODE's

% Switching parameters

NTB=1;              % NTB=1 = no translation blocking 
                    % NTB=0 - translation blocking (CHX)
IA=0;               % IA=1 - inactive STAT1

block=1;            % if 0 then there is no saturation in phosphorylation; 1 otherwise

output=0;           % output=1 - concentrations, output=0 - number of molecules

trans_knock=1;      %transcription knock-down: 1-no knock-down, 0-knock-down


% basic constants and rates
%---------------------------

Av=6.022e23;                  %Avogadro number

n_vol = 250e-15;              % nuclear volume HeLa cell - similar to Gorlich at al 2002
c_vol = n_vol;              % cytoplasm volume HeLa cell - similar to Gorlich at al 2002, corrected according to pictures

kv=c_vol/n_vol;

n_genes=2;                          % number of gene copies
v_max=n_genes*.16/Av/c_vol*1e6;     % [uM/s] - maximum transcription rate (per cytoplasmic volume)
v_irf1_transcription=trans_knock*.5e-5;
v_stat1_transcription=trans_knock*.023e-6;
v_stat2_transcription=trans_knock*0;%.1e-5;
v_lmp2_transcription =trans_knock*.2e-6;
v_tap1_transcription =trans_knock*.2e-6;
ktransl = 0.36;                      % [1/s] - translation rate for all transcripts

%ratio=0;          %percentage of STAT1 made inactive in the nucleus


% % plot scaling
% %---------------------------
% if (koniec > 2)
%     time_scale=3600;
% else
%     time_scale=60;
% end
% 
if (output == 0)
    cyt_scale=Av*c_vol/1e6;
    nuc_scale=Av*n_vol/1e6;
else
    cyt_scale=1;
    nuc_scale=1;
end

% time constants for inertial elements - in seconds

time_constant1=3600*.85;
time_constant2=3600*.85;
time_constant3=3600*.85;
time_constant4=3600*.85;
time_constant5=3600*.85;

time_constant_lmp=5*60;


constant1=1/time_constant1;
constant2=constant1*2;
constant3=constant1*4;
constant4=constant1*15;
constant5=constant1/5;
constant_lmp=1/time_constant_lmp;
% t4 = 3600*1.5;               % time constant
% t_2= t4*t4;
% t_3= t_2*t4;
% t_4= t_3*t4;

% Half-times of proteins and complexes [min]:
%--------------------------------------
               % transcripts
ts1t=119;     % STAT1 mRNA
ts2t=119;     % STAT2 mRNA
tirf1t=3*60;    % IRF1 mRNA
tlmp2t = 3*60;  % LMP2 mRNA
ttap1t = 2*60;  % TAP1 mRNA

               % proteins
ts1=15*60;     % STAT1 ts1=15*60;
ts2=15*60;     % STAT2
ti1_in=20;   % inactive IRF1
ti1=20;        % active IRF1

ts1pc=20;      % STAT1p in cytoplasm (for dephosphorylation)
ts2pc=20;      % STAT2p in cytoplasm (for dephosphorylation)

               % complexes
ts1s1p_c=25;   % STAT1p|STAT1p complex in cytoplasm =25;
ts1s1p_n=25;   % STAT1p|STAT1p complex in nucleus ts1s1p_n=25;
ts1s2p_c=25;   % STAT1p|STAT2p complex in cytoplasm =25;
ts1s2p_n=25;   % STAT1p|STAT2p complex in nucleus =25;
tpiass1s1=3;   % PIAS|STAT1p|STAT1p
ts1i1 = 4;     % STAT1|IRF1 
ts1i1cyt = 30; % STAT1|IRF1 cytoplasmic


% degradation rates for single molecules [1/s]
%-----------------------------------------------------
                               % proteins
ks1deg = log(2)/ts1/60;        % STAT1 
ks1pdeg = ks1deg;           % STAT1p
ks2deg = log(2)/ts2/60;        % STAT2
ks2pdeg = ks1pdeg;          % STAT2p
ki1deg = log(2)/ti1/60;        % IRF1 active
ki1_indeg = log(2)/ti1_in/60;        % IRF1 active

                               % mRNA
kdegs1t=log(2)/ts1t/60;        % STAT1 mRNA 
kdegs2t=log(2)/ts2t/60;        % STAT2 mRNA 
kdegi1t=log(2)/tirf1t/60;      % IRF1t
kdeglmp2t = log(2)/tlmp2t/60;      % LMP2 mRNA
kdegtap1t = log(2)/ttap1t/60;      % TAP1 mRNA


% degradation (dissociation) rates for complexes (into substrates) [1/s] 
%--------------------------------------------------------------------

% for dimers of phosphorylated STATs - it represents both dephosphorylation and degradation
kinvs1s1 = log(2)/ts1s1p_c/60;      % STAT1p|STAT1p
kinvs1s1n= log(2)/ts1s1p_n/60;      % (STAT1p|STAT1p)n
kinvs1s2 = log(2)/ts1s2p_c/60;      % STAT1p|STAT2p 
kinvs1s2n= log(2)/ts1s2p_n/60;      % (STAT1p|STAT2p)n

kinvpiass1s1 = log(2)/tpiass1s1/60; % PIAS|STAT1p|STAT1p 
kinvs1i1 = log(2)/ts1i1/60;         % STAT1|IRF1
kinvs1i1cyt = log(2)/ts1i1cyt/60;   % STAT1|IRF1 cytoplasmic

% degradation rates for complexes (non-reversible degradation) [1/s]
%--------------------------------------------------------------------

ks1s1pdeg = 0.05*kinvs1s1;                % (STAT1p|STAT1p) and (STAT1p|STAT1p)n ks1s1pdeg = 0.05*kinvs1s1
ks1s2pdeg = 0.1*kinvs1s2;                % (STAT1p|STAT2p) and (STAT1p|STAT2p)n ks1s2pdeg = 0.1*kinvs1s2;   
ks1i1deg = 0.01*kinvs1i1;                % STAT1|IRF1
ks1i1cytdeg = 0.01*kinvs1i1cyt;                % STAT1|IRF1 cytoplasmic
kdegpiass1s1 = 0;

% complexes creation rates [1/(uM*s)]
%---------------------------------------------
ks1s1=4;%log(2)/tc_s1s1/60;%.6;           % STAT1p|STAT1p 
ks1s2=10;%log(2)/tc_s1s2/60;%.9;          % STAT1p|STAT2p 
kpiass1s1 = .7;                           % PIAS|STAT1p|STAT1p
ks1i1 = .01;                                % STAT1|IRF1
ks1i1cyt = .001;                             % STAT1|IRF1 cytoplasmic


%Activation of proteins
%--------------------------------------------------------------
kactivation = .7e-3;   % PIAS
tpias_inact = 20*60;
kpias_inact = log(2)/tpias_inact/60;

tacti1=1;
kacti1= log(2)/tacti1/60;
tinacti1=.25*60;
kinacti1=log(2)/tinacti1/60;


% Transcription thresholds [uM]
%---------------------------
kI1=6e-2;       % for IFN - induced IRF1 transcription
ks1=1e-3;



%Initial number of molecules
%corresponding to steady state in unstimulated cells
%-----------------------------------------------------

stat1_molecules = 1e5;                % number of STAT1 molecules in an unstimulated cell
stat2_molecules = stat1_molecules/1.5;    % number of STAT2 molecules in an unstimulated cell
pias_molecules = .5e5;                 % number of free inactive PIAS molecules in an unstimulated cell

lmpt2_molecules = 8;
tap1t_molecules = 8;

perc_s1_c = .99;                       % percentage of STAT1 molecules in cytoplasm - Assumption - 90%
perc_s2_c = .99;                      % percentage of STAT2 molecules in cytoplasm - Assumption - 99%


%Initial cellular concentrations - all [uM] (initial conditions are below them)
%corresponding to steady state in unstimulated cells
%-----------------------------------------------------

pias_ss = pias_molecules/n_vol/Av*1e6;
stat1c_ss = perc_s1_c*stat1_molecules/c_vol/Av*1e6;            % STAT1 in cytoplasm -  of molecules in cytoplasm
tot_stat1n_ss = (1-perc_s1_c)*stat1_molecules/n_vol/Av*1e6;    % total STAT1 in nucleus (including complexes)-

stat1n_ss = tot_stat1n_ss;                                     % free STAT1 in nucleus (separate variable
                                                               % since some complexes can be added in the future)

stat2c_ss = perc_s2_c*stat2_molecules/c_vol/Av*1e6;            % STAT2 in cytoplasm
stat2n_ss = (1-perc_s2_c)*stat2_molecules/n_vol/Av*1e6;        % STAT2 in nucleus

% stat1t_ss = ks1deg/ktransl*(stat1c_ss  + stat1n_ss/kv);        % STAT1 mRNA
stat1t_ss = 2.59/c_vol/Av*1e6; 
                                   
stat2t_ss = ks2deg/ktransl*(stat2c_ss  + stat2n_ss/kv);

lmp2t_ss = lmpt2_molecules/c_vol/Av*1e6;
tap1t_ss = tap1t_molecules/c_vol/Av*1e6;
                                   
% Constitutive production rates [uM/s]
%------------------------------------
ks1tprod = kdegs1t*stat1t_ss;         % STAT1 mRNA
ks2tprod = kdegs2t*stat2t_ss;         % STAT2 mRNA
lmp2tprod = kdeglmp2t*lmp2t_ss;       % LMP2 mRNA
tap1tprod = kdegtap1t*tap1t_ss;       % TAP1 mRNA

% x protein is assumed to be produced only in stimulated cells
% PIAS is not degraded (actually assumed in the steady state concerning production and degradation),
% hence no constitutive production is needed

% Nuclear import/export rates [1/s]
%------------------------------------
transport=0.5;                   % [min] %8
es1=log(2)/transport/60;       % nuclear export of unphosphorylated STAT1 (Yamada)
kes1=0;
is1=stat1n_ss/stat1c_ss*(es1 + ks1deg/kv);
                               % nuclear import of unphosphorylated STAT1
                              
es2 = es1/1;                     % nuclear export of STAT2
is2 = (es2 + ks2deg/kv)*stat2n_ss/stat2c_ss;
                               % nuclear import of STAT2

transports1s1=.5;                % [min] 0.5
is1s1=log(2)/transports1s1/60;  % nuclear import of STAT1p|STAT1p
es1s1=0;               % nuclear export of STAT1p|STAT1p
transports1s2=transports1s1;    % it is done by the same importin alpha5
is1s2=log(2)/transports1s2/60;  % nuclear import of STAT1p|STAT2p
es1s2=0;               % nuclear export of STAT1p|STAT2p

transporti1=.5;
ii1 = log(2)/transporti1/60;    %nuclear import of active IRF1
ei1 = 1e-4*ii1;

transporti1_in=2;
ei1_in = log(2)/transporti1_in/60;    %nuclear export of inactive IRF1


%############### Initial conditions in [uM] ######################

y0=zeros(1,44);
y0(1)  = stat1c_ss;              % STAT1 unphosphorylated in cytoplasm
y0(2)  = stat2c_ss;              % STAT2 unphosphorylated in cytoplasm
y0(3)  = 0;                      % STAT1 phosphorylated in cytoplasm
y0(4)  = 0;                      % STAT2 phosphorylated in cytoplasm
y0(5)  = stat1n_ss;              % STAT1 unphosphorylated in nucleus 
y0(6)  = stat2n_ss;              % STAT2 unphosphorylated in nucleus 
y0(7) = 0;                      % (STAT1p|STAT1p) in cytoplasm
y0(8) = 0;                      % (STAT1p|STAT2p) in cytoplasm
y0(9) = 0;                      % (STAT1p|STAT1p) in nucleus
y0(10) = 0;                      % (STAT1p|STAT2p) in nucleus
y0(11) = stat1t_ss;              % STAT1 transcript     
y0(12) = 0;                      % active free PIAS
y0(13) = 0;                      % (PIAS|STAT1|STAT1)
y0(14) = 0;                      % inactive PIAS
y0(15)=0;                        % IRF1 mRNA
y0(16)=0;                        % IRF1 inactive protein
y0(17)=0;                        % IRF1 active in cytoplasm
y0(18)=0;                        % IRF1 active in nucleus
y0(19)=0;                        % IRF1 inactive in nucleus
y0(20)=0;
y0(21)=0;
y0(22)=0;
y0(23)=0;
y0(24)=stat2t_ss;
y0(25)=0;
y0(30)=lmp2t_ss;
y0(31)=tap1t_ss;
y0(32)=0;
y0(33)=0;
y0(34)=0;
y0(35)=0;
y0(36)=0.01; %inactive Tc-PTP in cytoplasm
y0(37)=0;
y0(38)=0.03; % inactive Tc-PTP in nucleus
y0(39)=0;

%Phosphorylation and dephosphorylation rates [1/s] 
%-------------------------------------------
% there is no dephosphorylation of STATs in nucleus due to assumption that there are no phosphorylated monomers there
% for a rate for dimers see the section with degradation of complexes into substrates
% phosphorylation, in turn, takes place only in cytoplasm
tphos1 = 32;           % phosphorylation of STAT1
tphos2 = 30;      % phosphorylation of STAT2

ks1phos  = log(2)/tphos1/60;           	  % phosphorylation rate for STAT1 - Yamada (?)
ks1dephc = log(2)/ts1pc/60;               % dephosphorylation of STAT1 in cytoplasm  - .5e-3 in Yamada (?)
ks2phos  = log(2)/tphos2/60; 	          % phosphorylation rate for STAT2
ks2dephc = log(2)/ts1pc/60;               % dephosphorylation of STAT2 in cytoplasm

ks1phos_sat = 1e4;                   %saturation of phosphorylation due to competition for receptors
ks2phos_sat = 1e4;

ks1phos=ks1phos*(1+block*ks1phos_sat*y0(1));
ks2phos=ks2phos*(1+block*ks2phos_sat*y0(2));


 %% New parameters
 
%Receptor dynamics
Ratot = 0.1; % IFNa/b receptor total amout 
kIFNaRa = 0.0016; % IFNa/b receptor activation 
kinvIFNaRa = log(2)/30/60;  % IFNa/b receptor inactivation (halftime = 30min) 
kIFNadeg = log(2)/80/60;% IFNa/b degradation

Rgtot = 0.1; % IFNg receptor total amout 
kIFNgRa = 0.0016; % IFNg receptor activation 
kinvIFNgRa = log(2)/120/60;  % IFNg receptor inactivation (halftime = 120min)
kIFNgdeg = log(2)/132/60;% IFNg degradation

Rst = 0.0025; % Receptor stimulation factor 
Rf = 0.7; % Receptor a/b1 eficacy factor 

%Tc-PTP
kactPTP = 2.8e-05; %activation rate of Tc-PTP 
kinactPTP = 0;%0.000001; %inactivation rate of Tc-PTP 
kPTPss = .35; %Tc-PTP | STATs complexes creation rate 
kinvPTPss = log(2)/3/60; %Tc-PTP | STATs complexes dissassociation rate (halftime = 3min)
MmPTP = 0.2; % M-M constant for Tc-PTP | STATs complexes creation
PTPinRa = 7.5e3; % Tc-PTP inhibitory factor for IFNa/b receptor actiovation 
PTPinRg = 1.5e3; % Tc-PTP inhibitory factor for IFNg receptor actiovation  
ePTP = log(2)/0.6/60; % nuclear export of Tc-PTP
iPTP = log(2)/0.2/60; % nuclear import of Tc-PTP



par=[kv v_irf1_transcription v_stat1_transcription v_stat2_transcription v_lmp2_transcription v_tap1_transcription ktransl constant1 constant2 constant3 constant4 constant5 ks1deg ks1pdeg ks2deg ks2pdeg ki1deg ki1_indeg kdegs1t kdegs2t kdegi1t kdeglmp2t kdegtap1t kinvs1s1 kinvs1s1n kinvs1s2 kinvs1s2n kinvpiass1s1 kinvs1i1 ks1s1pdeg ks1s2pdeg ks1i1deg kdegpiass1s1 ks1s1 ks1s2 kpiass1s1 ks1i1 kactivation kinacti1 ks1tprod ks2tprod lmp2tprod tap1tprod es1 is1 es2 is2 is1s1 is1s2 ii1 ei1 ei1_in ks1phos ks1dephc ks2phos ks2dephc ks1phos_sat ks2phos_sat, Ratot, kIFNaRa, kinvIFNaRa, kIFNadeg, Rgtot, kIFNgRa, kinvIFNgRa, kIFNgdeg, Rst, Rf, kactPTP, kinactPTP, kPTPss, kinvPTPss, MmPTP, PTPinRa, PTPinRg, ePTP, iPTP];

