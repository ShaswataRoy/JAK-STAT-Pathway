close all;
clear all;
clc

parameters;                % parametrization of ODE's and simulation times

%############## Simulation time points ##################################

t_total=50.00*3600;             % total length of simulation

% 1st puls of IFN
t_on1=0;                        % time when 1st puls of IFN is being introduced into the system
t_off1=1*3600;                  % length of 1st puls of IFN stimulation
dose_a1=0;                      % IFNa/b 1st puls dose ng/ml
dose_g1=100;                    % IFNgamma 1st puls dose ng/ml

% 2nd puls of IFN
t_on2=7*3600;                   % time when 2nd puls of IFN is being introduced into the system 
t_off2=1*3600;                  % length of 2nd puls of IFN stimulation
dose_a2=0;                      % IFNa/b 2nd puls dose ng/ml
dose_g2=100;                    % IFNgamma 2nd puls dose ng/ml

% 3rd puls of IFN
t_on3=24*3600;                  % time when 3rd puls of IFN is being introduced into the system 
t_off3=1*3600;                  % length of 3rd puls of IFN stimulation
dose_a3=0;                      % IFNa/b 3rd puls dose ng/ml
dose_g3=0;                      % IFNgamma 3rd puls dose ng/ml

k=1000;

dt=1;                           % sampling time
T=[];
Y=[];

% ##########################################################
% ######                 1st puls IFN                    ###
% ##########################################################

t_start=t_on1;
t_stop=t_on1+t_off1;
tspan1=t_start:dt:t_stop;
y0_1=y0;
y0_1(32)=dose_a1/k;
y0_1(34)=dose_g1/k;
[T1,Y1]=ode23tb(@model,tspan1,y0_1,[],par); 

T=[T; T1];
Y=[Y; Y1];

% ##########################################################
% ######            after 1st puls IFN                   ###
% ##########################################################

t_start=T(end);
t_stop=t_on2;
tspan2=t_start:dt:t_stop; 
y0_2=Y(end,:);
y0_2(32)=0;
y0_2(34)=0;
[T2,Y2]=ode23tb(@model,tspan2,y0_2,[],par); 

T=[T; T2(2:end)];
Y=[Y; Y2(2:end,:)];

% ##########################################################
% ######                2nd puls IFN                     ###
% ##########################################################

t_start=T(end);
t_stop=t_on2+t_off2;
tspan3=t_start:dt:t_stop; 
y0_3=Y(end,:);
y0_3(32)=dose_a2/k;
y0_3(34)=dose_g2/k;
[T3,Y3]=ode23tb(@model,tspan3,y0_3,[],par); 

T=[T; T3(2:end)];
Y=[Y; Y3(2:end,:)];

% ##########################################################
% ######            after 2nd puls IFN                   ###
% ##########################################################

t_start=T(end);
t_stop=t_on3;
tspan4=t_start:dt:t_stop; 
y0_4=Y(end,:);
y0_4(32)=0;
y0_4(34)=0;
[T4,Y4]=ode23tb(@model,tspan4,y0_4,[],par); 

T=[T; T4(2:end)];
Y=[Y; Y4(2:end,:)];

% ##########################################################
% ######                3rd puls IFN                     ###
% ##########################################################

t_start=T(end);
t_stop=t_on3+t_off3;
tspan5=t_start:dt:t_stop; 
y0_5=Y(end,:);
y0_5(32)=dose_a3/k;
y0_5(34)=dose_g3/k;
[T5,Y5]=ode23tb(@model,tspan5,y0_5,[],par); 

T=[T; T5(2:end)];
Y=[Y; Y5(2:end,:)];

% ##########################################################
% ######            after 3rd puls IFN                   ###
% ##########################################################

t_start=T(end);
t_stop=t_total;
tspan6=t_start:dt:t_stop; 
y0_6=Y(end,:);
y0_6(32)=0;
y0_6(34)=0;
[T4,Y4]=ode23tb(@model,tspan6,y0_6,[],par); 

T=[T; T4(2:end)];
Y=[Y; Y4(2:end,:)];


%% Total nuclear STAT1 (phosphorylated and unphosphorylated) = 2*Y(:,9)+Y(:,10)+2*Y(:,42)+Y(:,43)+Y(:,5)
figure()
set(0,'DefaultAxesFontSize',14)
plot(T/60, (2*Y(:,9)+Y(:,10)+2*Y(:,42)+Y(:,43)+Y(:,5))*nuc_scale, 'k', 'linewidth', 3);
set(gcf,'Position',[1 1 500 240])
ylim([0 6.5*10^4])
xlim([0 900])
set(gca,'XTick',(0:200:900))
set(gca,'XTickLabels',[])
set(gca,'YTick',(0:2:7)*10^4)
set(gca,'YTickLabels',[])
grid on
% print('STAT1_time','-dtiff','-r300')

