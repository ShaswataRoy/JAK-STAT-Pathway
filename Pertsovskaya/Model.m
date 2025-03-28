%
% ODE model
%
function  dy = model_IFN(t,y,p) 
   
 % variable names:
  S = y(1);
  A = y(2);
  Ap_c = y(3);
  Ap_n = y(4);
  r = y(5);
  R = y(6);
  f = y(7);
  F = y(8);
  a = y(9);
  
  % parameter names:
  b_S = p(1);
  
  b_exp = p(2);
  b_imp = p(3);
  
  b_ph = p(4);
  b_deph = p(5);
  
  k_A = p(6);
  k_I = p(7);
  k_r = p(8);
  k_f = p(9);
  k_F = p(10);
  
  b_A = p(11);
  b_r = p(12);
  b_R = p(13);
  b_f = p(14);
  b_F = p(15);
  b_a = p(16);
  
  lambda_S = p(17);
  lambda_r = p(18);
  lambda_R = p(19);
  lambda_f = p(20);
  lambda_F = p(21);
  lambda_a = p(22);
  
  q = p(23);
  n = p(24);
  m = p(25);
  u = p(26);
   
  lambda_stat = p(27);
  basal_stat = p(28);
  
  % equations:
  dy = zeros(9,1);
  
  % dy(1)=dS/dt: S active receptors (protein)                                     

  dy(1) = b_S - lambda_S*S; 
  
  % dy(2)=dA/dt: STAT1 non-phosphorylated form (protein)
  
  dy(2) = + b_exp*Ap_n ...
          + b_deph*Ap_c ...
          + b_A*a ...
          - b_ph*S*(A/k_A)/(1+(A/k_A)+power(R/k_I,q)) ...
          - lambda_stat*A; 
  
% dy(3)=dApc/dt: STAT1 phosphorylated form (protein)
  
  dy(3) = + b_ph*S*(A/k_A)/(1+(A/k_A)+power(R/k_I,q)) ...
          - b_imp* Ap_c ...
          - b_deph* Ap_c...
          - lambda_stat*Ap_c; 
        
% dy(4)=dApn/dt: STAT1 phosphorylated form in cell nucleus (protein)
  
  dy(4) = b_imp * Ap_c ...
          - b_exp* Ap_n...
          - lambda_stat*Ap_n;
           
% dy(5)=dr/dt: SOCS1 expression (mRNA)
  
  dy(5) = + b_r * power(Ap_n/k_r,n)/(1+power(Ap_n/k_r,n)) ...
          - lambda_r * r;

% dy(6)=dR/dt: SOCS1 (protein) 
  
  dy(6) = + b_R * r ...
          - lambda_R * R;
      
% dy(7)=df/dt: IRF1 expression (mRNA)
  
  dy(7) = + b_f * power(Ap_n/k_f,m)/(1+power(Ap_n/k_f,m)) ...
          - lambda_f * f;
      
% dy(8)=dF/dt: IRF1 (protein)

  dy(8) = + b_F * f ...
          - lambda_F * F;
  
% dy(9)=da/dt: STAT1 expression (mRNA)

  dy(9) = + b_a * power(F/k_F,u)/(1+power(F/k_F,u))-lambda_a * a+ basal_stat;
 
  
end
