 
% Solution for ODE IFN_model (Integration) 

% parameter values: 

p= [0, ...                % b_S
         0.08,0.013,...   % b_exp, b_imp 
         1300,0.036,...   % b_ph, b_deph 
         4680,82680,...   % k_A, k_I, 
         23400,7.3e+03,130e+03,...  % k_r, k_f, K_F  
         65,...            % b_A
         12.8, 1.0e+02,... % b_r, b_R
         2.7,10,...        % b_f, b_F
         1e-01,...         % b_a   
         0.02,...          % lambda_S
         0.03, 0.02,...    % lambda_r, lambda_R  
         0.017, 0.01,...   % lambda_f, lambda_F  
         0.006,...         % lambda_a  
         4, 3, 2, 1,...    % q, n, m, u
         6.9e-04, 0.006    % lambda_stat, basal_stat
   ];


tmax  = 500;
lcolor='k';    
options = odeset('MaxStep',1); % integration step limit 

% Solving ODE model in matlab (using a solver for stiff type problems):

sol = ode15s(@model_IFN,[0 tmax],[1000 1e+05 10 1 1 1 1 1 1],options,p); 
 
% Some graphs:

    figure(); hold on;
    plot(sol.x,(sol.y(3,:)),sprintf('%s-','b'),'LineWidth',1); hold on;
    legend('pSTAT1 oscillartory regime');
    ylabel('STAT1_p');
    
    figure(); hold on; 
    
    subplot(6,1,1); 
    plot(sol.x, sol.y(1,:),sprintf('%s-',lcolor),'LineWidth',1); hold on;
    legend('S');
    ylabel('receptor');
      
    subplot(6,1,2);
    
    plot(sol.x,(sol.y(4,:)),sprintf('%s-','k'),'LineWidth',1); hold on;
    legend('Ap_c+Ap_n');
    ylabel('STAT1_p');
    
    subplot(6,1,3);
    plot(sol.x,(sol.y(2,:)+sol.y(3,:)+sol.y(4,:)),sprintf('%s-','k'),'LineWidth',1); hold on;
    legend('A+Apc+Apn');
    ylabel('total STAT1');
 
    subplot(6,1,4);
    plot(sol.x, sol.y(5,:),sprintf('%s-','k'),'LineWidth',1);hold on;
    plot(sol.x, sol.y(7,:),sprintf('%s--','k'),'LineWidth',1); 
    ylabel('mRNA');
    legend('r','f');   

    subplot(6,1,5);
    plot(sol.x, sol.y(9,:),sprintf('%s:','k'),'LineWidth',1); hold on;
    ylabel('mRNA');
    legend('a'); 
    
    subplot(6,1,6);
    plot(sol.x, sol.y(6,:),sprintf('%s-',lcolor),'LineWidth',1); hold on; 
    plot(sol.x, sol.y(8,:),sprintf('%s--',lcolor),'LineWidth',1); 
    ylabel('downstream proteins');
    legend('R','F');
    xlabel('time (min)');    

%the end