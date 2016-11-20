%% Lab 2
clear all, close all, clc %#ok<*CLSCR,*CLALL,*DUALC>

%% User Input

    % System
    % 1 - Minimum Phase
    % 2 - Non-Minimum Phase
    sp = 2;

    % Task
    % 1 - Static Decoupling
    % 2 - Dynamical Decoupling
    % 3 - Glover-McFarlane
    tt = 2;
    
    % Plots
    % 0 - No Plots
    % 1 - Plots
    pp = 1;
    
%% System Initialization

if sp ~= 1 && sp ~= 2
    error('Error! Please Choose an Apprpriate System')
end

phm = pi/3;
s = tf('s');
switch sp
    case 1
        G = minphase;
        G_tf = tf(G);
        wc = 0.1;
    case 2
        G = nonminphase;
        G_tf = tf(G);
        wc = 0.02;
end

RGA = minreal(G_tf .* inv(G_tf)');
RGA_zero = evalfr(RGA,0) %#ok<*NOPTS>

W2 = eye(size(G_tf));

if tt == 1
    W1 = inv(evalfr(G_tf,0));
    Gtilde = G;
elseif tt == 2 || tt == 3
    w11 = -G_tf(2,2)/G_tf(2,1);
    w12 = -G_tf(1,2)/G_tf(1,1);
    w21 = -G_tf(2,1)/G_tf(2,2);
    w22 = -G_tf(1,1)/G_tf(1,2);
    W1 = [1 w12;w21 1];
    % W1 = [w11 1;1 w22]; W1 = (-1) * W1 * 10 * wc / (s + 10*wc);
    W1 = -minreal(ss(W1))
    Gtilde = G * W1;
    if pp == 1; figure; bode(Gtilde); title('Bode Diagram of Gtilde'); end
    isp = isproper(W1);
    if isp ~= 1; warning('W1 is not Proper'); end
else
    error('Please Choose an Apprpriate Task')
end



%% Core

% Minimum Phase
if sp == 2

[~,phG11] = bode(Gtilde(1,1),wc);
[~,phG22] = bode(Gtilde(2,2),wc);

argG11 = phG11*(pi/180);
argG22 = phG22*(pi/180);

T1 = (tan(phm-pi-argG11+pi/2))/wc;
T2 = (tan(phm-pi-argG22+pi/2))/wc;

L11 = Gtilde(1,1)*(1+(1/(s*T1)));
L22 = Gtilde(2,2)*(1+(1/(s*T2)));

[mL11,~] = bode(L11,wc);
[mL22,~] = bode(L22,wc);

K1 = (1/mL11);
K2 = (1/mL22);

end

% Non-Minimum Phase
if sp == 0

[~,phG12] = bode(Gtilde(1,2),wc);
[~,phG21] = bode(Gtilde(2,1),wc);

argG12 = phG12*(pi/180);
argG21 = phG21*(pi/180);

T1 = (tan(phm-pi-argG12+pi/2))/wc;
T2 = (tan(phm-pi-argG21+pi/2))/wc;

L12 = Gtilde(1,2)*(1+(1/(s*T1)));
L21 = Gtilde(2,1)*(1+(1/(s*T2)));

[mL12,~] = bode(L12,wc);
[mL21,~] = bode(L21,wc);

K1 = abs(1/mL12);
K2 = abs(1/mL21);

end

%% Controller Design

F1 = K1*(1+1/(s*T1));
F2 = K2*(1+1/(s*T2));
if sp == 2; Ftilde = [F1 0; 0 F2]; end
if sp == 0; Ftilde = [0 F1; F2 0]; end
        
F = W1*Ftilde;

if tt == 3 
L0 = G*W1*Ftilde;
alpha = 1.1;
[Fr,gam] = rloop(L0,alpha);
gam
F = minreal(W1*Ftilde*Fr);
end

L = minreal(G*F);


%% Singular Values

S = inv(eye(size(L))+L);
T = (inv(eye(size(L))+L))*L;
if pp == 1
figure
subplot(2,1,1)
sigma(S)
title('Singular values S')
subplot(2,1,2)
sigma(T)
title('Singular values T')
end
















