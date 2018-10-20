// SAMBA: Stochastic Analytical Model with a Bayesian Approach
// Based on: BCB Working Paper 239

// Code by: Paulo de Carvalho Lins and Victor Nomi
// IBRE/FGV paulo.lins@fgv.br; victor.nomi@fgv.br
// Rio de Janeiro, 14, May, 2018
// Prepared for the Macroeconomia Computável course, May 2018

//DEBUG MODE
//1. added autorregressive equation for preference shock rho_C 
//2. substituted Pie_C_bar --> pie_C_bar
%----------------------------------------------------------------
% 1. Variables
%----------------------------------------------------------------

var     y           c           c_O         c_RT        inv         g           x           m           m_C         m_I         m_X;
var     mrs;
var     k           n           y_D         y_D_C       y_D_I       y_D_X;        
var     pie_C       pie4_C      pie_A       pie_F       pie_M       pie_Y       pie_I       pie_X pie4_A pie4_F pie4_I pie4_X pie4_M;          
var     v_M         v_A         v_F         v_I         v_X;
var     mc_C        mc_I        mc_X        q_F         q_G         q_I         q_K         q_X         q_M         q_D         q_Y         r_K         w;
var     b_y_RW      l_y_RW      nx          nx_y        q           s_B_RW;   
var     b_y         pie_C_bar   pie4_C_bar  r           rr          s_B         s_y         s_y_bar;
var     tau;
var   pie_RW        q_M_RW      r_RW        y_RW        v_RW; 
var     z_A         z_B_RW      z_D         z_G         z_I         z_M         z_P         z_Q         z_R         z_W         z_M_RW     z_Z;            
var     z_C;
var     d_y_obs     d_c_obs     d_g_obs     d_inv_obs   d_x_obs     d_m_obs     d_n_obs     d_w_obs;

varexo  e_Pie_bar   e_S_bar;    
varexo  e_A         e_B         e_B_RW      e_D         e_G         e_I         e_M         e_P         e_Q         e_R         e_W         e_M_RW     e_Z;    
varexo  e_Pie_RW    e_Q_M_RW    e_R_RW      e_Y_RW      e_V_RW      e_C         e_TT;   

%----------------------------------------------------------------
% 2. Parameters
%----------------------------------------------------------------

// PREFERENCES:
parameters   gZ beta kappa varkappa sigma eta epsilon_W varphi_V_RW varphi_B_RW kappa_hat beta_hat;
gZ           = 1.01;
eta          = 1.00; 
epsilon_W    = 3.00;

//~~~~ The following are under examination - might be estimated:
varkappa     = 0.00;

//~~~~ The following will be estimated:
beta         = 0.988; // 0.9822 // RED FLAG //NEW estimation
kappa        = 0.75;  
sigma        = 2.00; // => 1.50           
varphi_V_RW  = 0.05;
varphi_B_RW  = 0.015;

// dependent parameters:
kappa_hat    = kappa/gZ;
beta_hat     = beta*(gZ)^(1-sigma);

// TECHNOLOGY:
parameters   alpha delta delta_I epsilon_I vartheta_I iota_I epsilon_C epsilon_C_P vartheta_C iota_C epsilon_X vartheta_X iota_X epsilon_M epsilon_RW vartheta_RW UIPH; 
alpha        = 0.4482;
delta        = 0.025; 
iota_I       = 0.50;
epsilon_C_P  = 11.0;
iota_C       = 0.50;
iota_X       = 0.50;
//deleted delta_Q_K

//~~~~ The following are under examination - might be estimated:
epsilon_M    = 11.0;
UIPH         = 0.06;
delta_I      = 5.00;  //=>2.0;       

//~~~~ The following will be estimated:
epsilon_I    = 1.50;
vartheta_I   = 2.50;
epsilon_C    = 1.50;     
vartheta_C   = 2.50;
epsilon_X    = 1.50;
vartheta_X   = 2.50;
epsilon_RW   = 1.50;
vartheta_RW  = 2.50;    

// POLICY:
parameters gamma_R gamma_Pie gamma_Y gamma_B gamma_S gamma_S_bar Pie_C_bar; 

Pie_C_bar = 1.011;

//~~~~ The following are under examination - might be estimated:
gamma_B   = 0.15; // 0.05 // RED FLAG // 0.3
gamma_S   = 0.50;
gamma_S_bar   = 0.35; 

//~~~~ The following will be estimated:
gamma_R   = 0.80; // => 0.7 
gamma_Pie = 1.60;         
gamma_Y   = 0.30; // 0.30 // RED FLAG //0.1


// PRICE RIGIDITY:
parameters theta_F omega_F theta_M omega_M theta_W omega_W theta_I omega_I theta_X omega_X theta_A upsilon_A varpi_P_C varpi_P_I varpi_P_X; 
theta_A    = 0.75;  

//~~~~ The following are under examination - might be estimated:
upsilon_A  = 0.60; // 0.20    
varpi_P_C  = 0.67;
varpi_P_I  = 0.22;
varpi_P_X  = 0.11;

//~~~~ The following will be estimated:
theta_F    = 0.75;  
omega_F    = 0.50;  
theta_M    = 0.75;   
omega_M    = 0.50;     
theta_W    = 0.75;  
omega_W    = 0.50; 
theta_I    = 0.75;  
omega_I    = 0.50;  
theta_X    = 0.75;  
omega_X    = 0.50; 

// STEADY STATE RELATIONS:
parameters  varpi_RT varpi_A varpi_D_I varpi_D_X;
varpi_A     = 0.30;
varpi_D_I   = 0.79;             
varpi_D_X   = 0.90;             

//~~~~ The following are under examination - might be estimated:
varpi_RT    = 0.40;
   
parameters  B_y B_y_RW Q Q_M_RW R_RW Pie_bar_RW S_B S_B_RW AF_B AF_B_RW;
B_y         = 2.00; // 1.00 / RED FLAG
B_y_RW      =-0.68;
R_RW        = 1.0074;
Pie_bar_RW  = 1.007;
S_B_RW      = 1.0142;
AF_B        = 0.55;
AF_B_RW     = 0.04;

//~~~~ The following are under examination - might be estimated:
Q           = 1.00;
Q_M_RW      = 1.05;
S_B         = 1.00;

parameters  s_C s_I s_G s_X tt; 
s_C         = 0.62;
s_I         = 0.17;
s_G         = 0.20;
s_X         = 0.13;
tt          = 0.35;

// REDUCED-FORM:
parameters  s_M s_D varpi_D_C s_M_C s_M_I s_M_X s_D_C s_D_G s_D_I s_D_X s_N s_K s_N_bar alpha_N Q_M MC_C R Z_B_RW L_y_RW;
s_M         = s_C + s_I + s_G + s_X - 1;
MC_C        = (epsilon_C_P-1)/epsilon_C_P; 
Q_M         = epsilon_M/(epsilon_M-1)*Q*Q_M_RW;
s_D         = 1 - (1 - MC_C)*s_C +(1-Q_M/Q_M_RW*(iota_C*(R_RW*S_B_RW-1)+1))*s_M + (1 - varpi_D_I)*s_I*((iota_C*(R_RW*S_B_RW-1)+1)/(iota_I*(R_RW*S_B_RW-1)+1)-1) + (1 - varpi_D_X)*s_X*((iota_C*(R_RW*S_B_RW-1)+1)/(iota_X*(R_RW*S_B_RW-1)+1)-1);
varpi_D_C   = (s_D - s_G - varpi_D_I *s_I - varpi_D_X*s_X)/(MC_C*s_C);
s_M_C       = s_C*(1 - varpi_D_C)*MC_C/(iota_C*(R_RW*S_B_RW-1)+1);
s_M_I       = s_I*(1 - varpi_D_I)/(iota_I*(R_RW*S_B_RW-1)+1);
s_M_X       = s_X*(1 - varpi_D_X)/(iota_X*(R_RW*S_B_RW-1)+1);
s_D_C       = s_C*varpi_D_C*MC_C;
s_D_G       = s_G;
s_D_I       = s_I*varpi_D_I;
s_D_X       = s_X*varpi_D_X;
s_N         = (1-alpha)*s_D;
s_K         = alpha*s_D;
s_N_bar     = 0.20*s_N;
alpha_N     = 1/(s_N/(s_N-s_N_bar));
R           = Pie_C_bar*gZ/beta_hat;
Z_B_RW      = S_B*R*Pie_bar_RW/(R_RW*Pie_C_bar);
L_y_RW      = (R_RW*S_B_RW-1)*(iota_C*s_M_C + iota_I*s_M_I + iota_X*s_M_X);

parameters  lambda_F lambda_M lambda_W lambda1_W lambda2_W lambda3_W lambda4_W lambda1_B lambda1_B_RW lambda_I lambda_X; 
lambda_F    = (1-theta_F)*(1-theta_F*beta_hat)/theta_F;
lambda_M    = (1-theta_M)*(1-theta_M*beta_hat)/theta_M;
lambda_W    = (1-theta_W)*(1-theta_W*beta_hat)/(theta_W*(1+eta*epsilon_W)*(1+beta_hat));
lambda1_W   = (1-theta_W)*(1-theta_W*beta_hat)/(theta_W*(1+eta*epsilon_W)*(1+omega_W*beta_hat));
lambda2_W   = omega_W/(1+omega_W*beta_hat);
lambda3_W   = beta_hat/(1+omega_W*beta_hat);
lambda4_W   = 1/(1+omega_W*beta_hat);
lambda1_B   = R/(gZ*Pie_C_bar);
lambda1_B_RW= Pie_C_bar*R_RW*S_B_RW/(Pie_C_bar*Pie_bar_RW*gZ);
lambda_I    = (1-theta_I)*(1-theta_I*beta_hat)/theta_I;
lambda_X    = (1-theta_X)*(1-theta_X*beta_hat)/theta_X;

parameters  varpi_RW_C varpi_RW_I varpi_RW_X;
varpi_RW_C  = iota_C*R_RW*S_B_RW/(1+iota_C*(R_RW*S_B_RW-1));
varpi_RW_I  = iota_I*R_RW*S_B_RW/(1+iota_I*(R_RW*S_B_RW-1));
varpi_RW_X  = iota_X*R_RW*S_B_RW/(1+iota_X*(R_RW*S_B_RW-1));

// SHOCKS:
parameters rho_Pie_bar rho_S_bar rho_B rho_Y_RW rho_V_RW rho_M_RW rho_R_RW rho_Pie_RW rho_Z rho_B_RW rho_I; 
parameters rho_X rho_M rho_A rho_P rho_W rho_R rho_G rho_Q rho_D rho_C rho_tau;

//~~~~ The following will be estimated:
//rho_Pie_bar = 0.95;   
//rho_S_bar   = 0.95;   
rho_Pie_bar = 0.95;   
rho_S_bar   = 0.95;   
rho_B       = 0.90;   
rho_Y_RW    = 0.80;   
rho_V_RW    = 0.50;     
rho_M_RW    = 0.50;  
rho_R_RW    = 0.90;    
rho_Pie_RW  = 0.50;  
rho_Z       = 0.80;   
rho_B_RW    = 0.90;     
rho_I       = 0.70;     
rho_X       = 0.80;    
rho_M       = 0.85;     
rho_A       = 0.50;    
rho_P       = 0.50;  
rho_W       = 0.70;   
rho_R       = 0.01; 
rho_G       = 0.80;  
rho_Q       = 0.50;   
rho_D       = 0.90;
rho_C       = 0.50; 
rho_tau     = 0.50;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model (linear);

// Households 11
c_O  = (kappa_hat/(1+kappa_hat))*c_O(-1) + (1/(1+kappa_hat))*c_O(+1) - ((1-kappa_hat)/(sigma*(1+kappa_hat)))*(r + s_B - pie_C(+1)) + ((rho_Z-kappa_hat)/(1+kappa_hat))*z_Z - (1 - rho_C)*(1-kappa_hat)*(1/(sigma*(1+kappa_hat)))*z_C; // (C.1)
c_RT = w + n - (tt/(1-tt))*tau; // (C.2)
c    = varpi_RT*c_RT + (1-varpi_RT)*c_O; // (C.3)
q    = q(+1) - ((r + s_B - pie_C(+1)) + (r_RW + s_B_RW - pie_RW(+1))) + UIPH*q(-1) + z_Q; // (C.4) Hibrida
q_K = (beta_hat/gZ)*(1-delta)*q_K(+1) + (1-(1-delta)*beta_hat/gZ)*r_K(+1) - (r + s_B - pie_C(+1)); // (C.5)
inv  = (1/(1+beta_hat))*inv(-1) + (beta_hat/(1+beta_hat))*inv(+1) + (1/(delta_I*gZ^2*(1+beta_hat)))*(q_K - q_I) - ((1-beta_hat*rho_Z)/(1+beta_hat))*z_Z + ((1-beta_hat*rho_I)/(1+beta_hat))*z_I; // (C.6) Wrong equation in the WP. Used rho_I, instead of rho_Z
k     = (1 - delta)*(1/gZ)*(k(-1) - z_Z) + (1-(1-delta)/gZ)*inv; // (C.7)
s_B_RW = - varphi_B_RW*b_y_RW + varphi_V_RW*v_RW + z_B_RW; // (C.8)
s_B      = rho_B*s_B(-1) + e_B; // (C.9)
w   = w(-1) + lambda1_W*(mrs - w ) + z_W + lambda2_W*(w(-1)-w(-2)) + lambda3_W*(w(+1)-w)-(lambda3_W+lambda4_W)*(pie_C+z_Z)+lambda4_W*(pie_C(-1)+z_Z(-1))+lambda3_W*(z_Z*rho_Z + pie_C(+1)); // (C.10)
mrs = eta*n + (sigma/(1-kappa_hat))*(c_O - kappa_hat*(c_O(-1) - z_Z));

// Firms 34
r_K = q_D + y_D - k(-1) + z_Z; // (C.11)
q_D = alpha*r_K + (1-alpha)*w - z_D; // (C.12)
n = alpha_N*(y_D + q_D - w); // (C.13)
pie_M - v_M = lambda_M*(q + q_M_RW - q_M) + beta_hat*(pie_M(+1) - v_M(+1)); // (C.14)
v_M = omega_M*pie_M(-1) + (1-omega_M)*pie_C_bar;
q_M = q_M(-1) + pie_M - pie_C; // 
mc_C = varpi_D_C*q_D + (1-varpi_D_C)*(q_M + varpi_RW_C*(r_RW + s_B_RW) + vartheta_C*((m_C - c) - (m_C(-1) - c(-1))) - z_M); // (C.15.1)
mc_I = varpi_D_I*q_D + (1-varpi_D_I)*(q_M + varpi_RW_I*(r_RW + s_B_RW) + vartheta_I*((m_I - inv) - (m_I(-1) - inv(-1))) - z_M); // (C.15.2)
mc_X = varpi_D_X*q_D + (1-varpi_D_X)*(q_M + varpi_RW_X*(r_RW + s_B_RW) + vartheta_X*((m_X - x) - (m_X(-1) - x(-1))) - z_M); // (C.15.3)
y_D_C = c - epsilon_C*(q_D - mc_C); // (C.16.1)
y_D_I = inv - epsilon_I*(q_D - mc_I); // (C.16.2)
y_D_X = x - epsilon_X*(q_D - mc_X); // (C.16.3)
m_C  = c + (epsilon_C/(1+epsilon_C*vartheta_C))*(- (q_M + varpi_RW_C*(r_RW + s_B_RW) - mc_C) + vartheta_C*(m_C(-1) - c(-1) + z_M)  ); // (C.17.1)
m_I  = inv + (epsilon_I/(1+epsilon_I*vartheta_I))*(- (q_M + varpi_RW_I*(r_RW + s_B_RW) - q_I) + vartheta_I*(m_I(-1) - inv(-1) + z_M)  ); // (C.17.2)
m_X  = x + (epsilon_X/(1+epsilon_X*vartheta_X))*(- (q_M + varpi_RW_X*(r_RW + s_B_RW) - q_X) + vartheta_X*(m_X(-1) - x(-1) + z_M) ); // (C.17.3)
pie4_I = pie_I + pie_I(-1) + pie_I(-2) + pie_I(-3);
pie_I  = v_I + lambda_I*(mc_I - q_I) + beta_hat*(pie_I(+1) - v_I(+1)) + z_P; // (C.18)
v_I = omega_I*pie_I(-1) + (1-omega_I)*pie_C_bar;
q_I = q_I(-1) + pie_I - pie_C; // (C.19)
// No pie_G e q_G (C.18) e (C.19) --> we considered q_G=1
pie_F  = v_F + lambda_F*(mc_C - q_F) + beta_hat*(pie_F(+1) - v_F(+1)) + z_P; // (C.20)
v_F = omega_F*pie_C(-1) + (1-omega_F)*pie_C_bar; // Wrong equation in the WP (V_F instead of V_C). 
q_F = q_F(-1) + pie_F - pie_C; // (C.21)
pie_A  = theta_A*pie_C_bar + (1-theta_A)*v_A; // (C.22)
v_A = (pie4_C(-1)/4 + upsilon_A*(mc_C(-1) + mc_C(-2) + mc_C(-3) + mc_C(-4))/4) + z_A; // (C.23), qsiA=1 and thetaA=1
pie4_A = pie_A + pie_A(-1) + pie_A(-2) + pie_A(-3);
pie4_C = pie_C + pie_C(-1) + pie_C(-2) + pie_C(-3); // (C.23)
pie4_F = pie_F + pie_F(-1) + pie_F(-2) + pie_F(-3);
pie_C  = varpi_A*pie_A + (1-varpi_A)*pie_F; // (C.24)
pie_X  = v_X + lambda_X*(mc_X - q_X -q ) + beta_hat*(pie_X(+1) - v_X(+1)) + z_P; // (C.25) Not include q, because of BQ.
v_X = omega_X*pie_X(-1) + (1-omega_X)*pie_C_bar;
q_X = q_X(-1) + pie_X - pie_C; // (C.26)
x    = y_RW + (epsilon_RW/(1+vartheta_RW*epsilon_RW))*(vartheta_RW*(x(-1) - y_RW(-1)) + q_X + z_M_RW); // (C.27)
pie4_X = pie_X + pie_X(-1) + pie_X(-2) + pie_X(-3);
pie4_M = pie_M + pie_M(-1) + pie_M(-2) + pie_M(-3);


// Government 7
r = gamma_R*r(-1) + (1-gamma_R)*(pie_C_bar + sigma*z_Z + gamma_Pie*(pie4_C(+4)/4 - pie4_C_bar(+4)/4) + gamma_Y*y(-1)) + z_R; // (C.28)
pie_C_bar= rho_Pie_bar*pie_C_bar(-1) + e_Pie_bar; // (C.29)
// s_y = gamma_S*s_y(-1) + gamma_S_bar*s_y_bar  + s_G*z_G; // (C.30) Not satisfing BQ
g    = 0.5*g(-1) + 0.5*(gamma_S*(s_y(-1) - s_y_bar(-1)) - gamma_B*b_y(-1)) + z_G; // Used instead of (C.30)
 s_y_bar  = rho_S_bar*s_y_bar(-1) + gamma_B*b_y(-1) + e_S_bar; // (C.31) Not satisfing BQ
//s_y_bar  = rho_S_bar*s_y_bar(-1) + e_S_bar; // Used instead of (C.31)
tau = rho_tau*tau(-1) + e_TT; // (C.32)
// g = (1/s_G)*(tau - s_y) + y + q_Y - q_G; // Strange dynamic after the e_TT shock (C.33)
s_y  = s_G*(y - g + q_Y - q_G); // Used instead of (C.33)
b_y = lambda1_B*b_y(-1) + B_y*AF_B*r  - R*s_y + lambda1_B*B_y*(y(-1) - y - pie_Y - z_Z); // (C.34)

// Resource constraints and external sector 8
y_D   = (s_D_C/s_D)*y_D_C + (s_D_G/s_D)*g + (s_D_I/s_D)*y_D_I + (s_D_X/s_D)*y_D_X; // (C.35)
m    = (s_M_C/s_M)*m_C + (s_M_I/s_M)*m_I + (s_M_X/s_M)*m_X; // (C.36)
l_y_RW = (R_RW*S_B_RW-1)*(iota_C*s_M_C*m_C + iota_I*s_M_I*m_I + iota_X*s_M_X*m_X) + R_RW*S_B_RW*(iota_C*s_M_C + iota_I*s_M_I + iota_X*s_M_X)*(r_RW + s_B_RW) + L_y_RW*(q_M - q_Y - y); // (C.37)
nx_y   = s_X*(x + q_X + q ) - s_M*(m + q + q_M_RW) - (s_X - s_M)*(y + q_Y); // (C.38) Not include q, because of BQ.
b_y_RW = lambda1_B_RW*b_y_RW(-1) + B_y_RW*AF_B_RW*(r_RW + s_B_RW) + R_RW*S_B_RW*nx_y + lambda1_B_RW*B_y_RW*(y(-1) - y - pie_Y - z_Z) - R_RW*S_B_RW*l_y_RW + lambda1_B_RW*B_y_RW*(q - q(-1) + pie_C - pie_RW); // (C.39)
y    = s_C*c + s_I*inv + s_G*g + s_X*x - s_M*m; // (C.40)
q_Y = s_I*q_I + s_G*q_G + s_X*(q_X + q) - s_M*(q + q_M_RW); // (C.41) Not include q, because of BQ.
pie_Y  = pie_C + q_Y - q_Y(-1); // (C.42)

// Rest of the World 5
y_RW    = rho_Y_RW*y_RW(-1) + e_Y_RW; // (C.43)
q_M_RW  = rho_M_RW*q_M_RW(-1) + e_Q_M_RW; // (C.44)
pie_RW  = rho_Pie_RW*pie_RW(-1) + e_Pie_RW; // (C.45)
v_RW    = rho_V_RW*v_RW(-1) + e_V_RW; // (C.46)
r_RW    = rho_R_RW*r_RW(-1) + e_R_RW; // (C.47)

// Shock - AR(1) Processes 13
z_C   = rho_C*z_C(-1) + e_C; // (C.48)
z_Q   = rho_Q*z_Q(-1) + e_Q; // (C.49)
z_B_RW= rho_B_RW*z_B_RW(-1) + e_B_RW; // (C.50)
z_D   = rho_D*z_D(-1) + e_D; // (C.51)
z_Z   = rho_Z*z_Z(-1) + e_Z; // (C.52)
z_I   = rho_I*z_I(-1) + e_I; // (C.53)
z_M   = rho_M*z_M(-1) + e_M; // (C.54)
z_M_RW= rho_M_RW*z_M_RW(-1) + e_M_RW; // (C.55)
z_W   = rho_W*z_W(-1) + e_W; // (C.56)
z_P   = rho_P*z_P(-1) + e_P; // (C.57)
z_A   = rho_A*z_A(-1) + e_A; // (C.58)
// Missing (C.59)
z_R   = rho_R*z_R(-1) + e_R; // (C.60)
z_G   = rho_G*z_G(-1) + e_G; // (C.61)

// Additional Variables not in the Working Paper 4
rr  = r - pie_C(+1); // Real interest rate
pie4_C_bar= pie_C_bar + pie_C_bar(-1) + pie_C_bar(-2) + pie_C_bar(-3);
q_G = q_D + z_P;
nx     = (s_X/(s_X-s_M))*(x + q_X) - (s_M/(s_X-s_M))*(m + q + q_M_RW);

// Observable variables
d_y_obs = y - y(-1);
d_c_obs	= c - c(-1);
d_g_obs	= g - g(-1);
d_inv_obs = inv - inv(-1);
d_x_obs	= x - x(-1);
d_m_obs	= m - m(-1);
d_n_obs	= n - n(-1);
d_w_obs	= w - w(-1);

end;

initval;
y         = 0;
c         = 0;
c_O       = 0;
c_RT      = 0;
inv       = 0; 
g         = 0;
x         = 0;
m         = 0;
m_C       = 0;
m_I       = 0;
m_X       = 0;
k         = 0;
n         = 0;
y_D       = 0;
y_D_C     = 0;
pie4_X    = 0;
pie4_I    = 0;
pie4_M    = 0;
y_D_I     = 0;
y_D_X     = 0;
pie_C     = 0;
pie4_C    = 0;
pie4_A    = 0;
pie4_F    = 0;
pie_A     = 0;
pie_F     = 0;
pie_M     = 0;
pie_Y     = 0;
pie_I     = 0;
pie_X     = 0;
mc_C      = 0;
mc_I      = 0;
m_X       = 0;
q_F       = 0;
q_G       = 0;
q_I       = 0;
q_K       = 0;
q_X       = 0;
q_M       = 0;
q_D       = 0;
q_Y       = 0;
r_K       = 0;
w         = 0;
b_y_RW    = 0;
l_y_RW    = 0;
nx        = 0;
nx_y      = 0;
q         = 0;
s_B_RW    = 0;
b_y       = 0;
pie_C_bar = 0;
r         = 0;
rr        = 0;
s_B       = 0;
s_y       = 0;
s_y_bar   = 0;
pie_RW    = 0;
q_M_RW    = 0;
r_RW      = 0;
y_RW      = 0;
v_RW      = 0;
z_A       = 0;
z_B_RW    = 0;
z_D       = 0;
z_G       = 0;
z_I       = 0;
z_M       = 0;
z_P       = 0;
z_Q       = 0;
z_R       = 0;
z_W       = 0;
z_M_RW    = 0;
z_Z       = 0;
z_C       = 0;
tau       = 0;
mrs       = 0;
v_M       = 0; 
v_A       = 0;   
v_F       = 0;  
v_I       = 0;  
v_X       = 0;

d_y_obs 	= 0;
d_c_obs		= 0;
d_g_obs		= 0;
d_inv_obs	= 0;
d_x_obs		= 0;
d_m_obs		= 0;
d_n_obs		= 0;
d_w_obs		= 0;

end;

%----------------------------------------------------------------
% 4. ~~~~ Estimation ~~~~
%----------------------------------------------------------------

varobs d_c_obs d_g_obs d_inv_obs d_x_obs d_m_obs d_w_obs s_y_bar s_y pie4_C_bar pie_A pie_F r q y_RW r_RW pie_RW q_M_RW;
//removed again: d_y_obs d_n_obs pie_C 

//removed again: d_y_obs d_n_obs pie_C 
//previous state of v3: varobs c g inv x m w s_y_bar s_y pie_F y_RW q_M_RW;

//list of renamed variables: d_y_obs d_c_obs d_g_obs d_inv_obs d_x_obs d_m_obs  d_n_obs d_w_obs
//full: y c g inv x m n w s_y_bar s_y pie4_C_bar pie_C pie_A pie_F r q y_RW r_RW pie_RW q_M_RW;
//Never again: y pie_C n 
//near full: c g inv x m w s_y_bar s_y pie4_C_bar pie_A pie_F r q y_RW r_RW pie_RW q_M_RW;
//Deleted in v3_i1: pie4_C_bar pie_A r_RW 
//Deleted in v3_i2: r q pie_RW
//Deleted in v3_i3: s_y pie_F y_RW (equal to v2)

estimated_params;
//beta,           normal_pdf,     0.99,    0.005;
//kappa,          beta_pdf,       0.85,    0.05;
sigma,          normal_pdf,     1.3,     0.25;
//varphi_V_RW,    inv_gamma_pdf,  0.05,    0.15;
//varphi_B_RW,    inv_gamma_pdf,  0.05,    0.15;
//epsilon_I,      gamma_pdf,      1,       0.5;
//vartheta_I,     gamma_pdf,      4,       1.5;
//epsilon_C,      gamma_pdf,      1,       0.5;
//vartheta_C,     gamma_pdf,      4,       1.5;
//epsilon_X,      gamma_pdf,      1,       0.5;
//vartheta_X,     gamma_pdf,      4,       1.5;
//epsilon_RW,     gamma_pdf,      1,       0.5;
//vartheta_RW,    gamma_pdf,      4,       1.5;
//gamma_R,        beta_pdf,       0.6,     0.15;
//gamma_Y,        gamma_pdf,      0.25,    0.1;
//gamma_Pie,      normal_pdf,     2,       0.35;
//theta_F,        beta_pdf,       0.65,    0.1;
//omega_F,        beta_pdf,       0.65,    0.2;
//theta_M,        beta_pdf,       0.65,    0.1;
//omega_M,        beta_pdf,       0.65,    0.2;
//theta_W,        beta_pdf,       0.75,    0.1;
//omega_W,        beta_pdf,       0.65,    0.15;
//theta_I,        beta_pdf,       0.65,    0.1;
//omega_I,        beta_pdf,       0.65,    0.2;
//theta_X,        beta_pdf,       0.65,    0.1;
//omega_X,        beta_pdf,       0.65,    0.2;

//rho_Pie_bar,    beta_pdf,       0.5,     0.25;
//rho_S_bar,      beta_pdf,       0.5,     0.25;
//rho_B,          beta_pdf,       0.5,     0.25;
//rho_Y_RW,       beta_pdf,       0.5,     0.25;
//rho_V_RW,       beta_pdf,       0.5,     0.25;
//rho_M_RW,       beta_pdf,       0.5,     0.25;
//rho_R_RW,       beta_pdf,       0.5,     0.25;
//rho_Pie_RW,     beta_pdf,       0.5,     0.25;
//rho_Z,          beta_pdf,       0.5,     0.25;
//rho_B_RW,       beta_pdf,       0.5,     0.25;
//rho_I,          beta_pdf,       0.5,     0.25;
//rho_M,          beta_pdf,       0.5,     0.25;
//rho_A,          beta_pdf,       0.5,     0.25;
//rho_P,          beta_pdf,       0.5,     0.25;
//rho_W,          beta_pdf,       0.5,     0.25;
//rho_Q,          beta_pdf,       0.5,     0.25;
//rho_D,          beta_pdf,       0.5,     0.25;
//rho_C,          beta_pdf,       0.5,     0.25;

//stderr e_Pie_bar,      inv_gamma_pdf,  1,       inf;
//stderr e_S_bar,        inv_gamma_pdf,  1,       inf;
//stderr e_A,            inv_gamma_pdf,  1,       inf;
//stderr e_B,            inv_gamma_pdf,  1,       inf;
//stderr e_B_RW,         inv_gamma_pdf,  1,       inf;
//stderr e_D,            inv_gamma_pdf,  1,       inf;
//stderr e_G,            inv_gamma_pdf,  1,       inf;
//stderr e_I,            inv_gamma_pdf,  1,       inf;
//stderr e_M,            inv_gamma_pdf,  1,       inf;
//stderr e_P,            inv_gamma_pdf,  1,       inf;
//stderr e_Q,            inv_gamma_pdf,  1,       inf;
//stderr e_R,            inv_gamma_pdf,  1,       inf;
//stderr e_W,            inv_gamma_pdf,  1,       inf;
//stderr e_M_RW,         inv_gamma_pdf,  1,       inf;
//stderr e_Z,            inv_gamma_pdf,  1,       inf;
//stderr e_Pie_RW,       inv_gamma_pdf,  1,       inf;
//stderr e_Q_M_RW,       inv_gamma_pdf,  1,       inf;
//stderr e_R_RW,         inv_gamma_pdf,  1,       inf;
//stderr e_Y_RW,         inv_gamma_pdf,  1,       inf;
//stderr e_V_RW,         inv_gamma_pdf,  1,       inf;
//stderr e_C,            inv_gamma_pdf,  1,       inf;
//stderr e_TT,           inv_gamma_pdf,  1,       inf;
end;

// INTRODUCING PARAMETER BOUNDS
estimated_params_bounds;
sigma,  0.25,    inf;
end;

estimated_params_init;
sigma,  0.82;
end;
%-------------------------------------------------------------
% 5. Computation
%-------------------------------------------------------------

shocks;
var e_Pie_bar;	stderr 1.00;
var e_S_bar;	stderr 1.00;
var e_A; 		stderr 1.00;
var e_B;	    stderr 1.00;
var e_B_RW;  	stderr 1.00;
var e_D; 		stderr 1.00;
var e_G; 		stderr 1.00;
var e_I; 		stderr 1.00;
var e_M; 		stderr 1.00;
var e_P;     	stderr 1.00;
var e_Q; 		stderr 1.00;
var e_R; 		stderr 1.00;
var e_W; 		stderr 1.00;
var e_M_RW;     stderr 1.00;
var e_Z; 		stderr 1.00;
var e_Pie_RW;	stderr 1.00;
var e_Q_M_RW;   stderr 1.00; 
var e_R_RW;     stderr 1.00; 
var e_Y_RW;     stderr 1.00; 
var e_V_RW;     stderr 1.00;
var e_C;        stderr 1.00;
var e_TT;       stderr 1.00;
end; 

resid(1);
//estimation(datafile=samba_data_to_dynare, mh_replic=1000, mh_nblocks=3, mh_jscale=10, mode_compute=6) y c inv;
//estimation(datafile=samba_data_to_dynare, plot_priors=0, mh_replic=5000, mh_nblocks=3, mh_jscale=1, mcmc_jumping_covariance = identity_matrix) y c inv;

estimation(datafile=samba_data_to_dynare, plot_priors=0, mh_replic=5000, mh_nblocks=3, mh_jscale=0.01, nobs=57, mcmc_jumping_covariance = identity_matrix);
//estimation(datafile=samba_data_to_dynare, plot_priors=0, mh_replic=5000, mh_nblocks=3, mh_jscale=0.5);
//estimation(datafile=samba_data_to_dynare, plot_priors=0, mh_replic=5000, mh_nblocks=3, mh_jscale=1, mcmc_jumping_covariance = identity_matrix, nobs = 64);
//estimation(datafile=samba_data_to_dynare, plot_priors=0, mh_replic=5000, mh_nblocks=3, mh_jscale=0.1) y c inv;
//estimation(datafile=samba_data_to_dynare, plot_priors=0, mh_replic=5000, mh_nblocks=3, mh_jscale=0.1, nobs = 64) y c inv;
//estimation(datafile=samba_data_to_dynare, plot_priors=0, mh_replic=2000, mh_nblocks=3, mh_jscale=0.1, mcmc_jumping_covariance = identity_matrix) y c inv;

//steady;
//check;
//stoch_simul(irf=0) s_y_bar s_y b_y g;
