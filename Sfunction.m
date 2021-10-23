
function [sys, x0,str,ts]=Sfunction(t,x,u,flag)

% Este conjunto de lineas de codigo NO se modifica.

switch flag
   case 0
   [sys, x0,str,ts]=mdlInitializeSizes;
   case 1
   sys=mdlDerivatives(t,x,u);
   case 3
   sys=mdlOutputs(t,x,u);
   case {2,4,9}
   sys=[];
   otherwise
     error(['Unhandled flag = ',num2str(flag)]);
end

end



function [sys,x0,str,ts]=mdlInitializeSizes

% Aqui se declara el numero de ecuaciones diferenciales a integrar, de va-
% riables de entrada que ingresan al macro y de salidas que este tendra.

sizes=simsizes;
sizes.NumContStates  = 6; % Numero de ecuaciones diferenciales a integrar
sizes.NumDiscStates  = 0;  
sizes.NumOutputs     = 6; % Numero de variables de salida que tendra el 
%                           macro.
sizes.NumInputs      = 5; % Numero de variables de entrada que el macro 
%                           aceptara. T entrada, T chaqueta, Flujo alim
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;


sys=simsizes(sizes); % Expresion que permite que el macro funcione bien. No 
% se debe borrar.

% Este es el SS o condiciones iniciales que pueden editar después de la parte a).
x0 = [0.4, 0.5, 0, 6.8, 310, 210]; %X0, S0, E0, V0, Tmosto, Tchaqueta
% Dewasme: 0.4, 0.5, 0.8, 6.8
% Biener 2012: 10, 

% Estas 2 lineas de codigo NO se tocan.
str=[]; 
ts=[0 0];
end

function sys=mdlDerivatives(~,x,u)

% Variables de Estado    
X  = x(1); %[g/L] biomasa
S  = x(2); %[g/L] sustrato (glucosa)
P  = x(3); %[g/L] subproducto (etanol)
V  = x(4); %[L] volumen de cultivo
Tm = x(5); %[K] temperatura del mosto
Tj = x(6); %[K] temperatura de la chaqueta


% Entradas
Falim = u(1); %[L/h] flujo de alimentacion de glucosa
Fj    = u(2); %[L/h] flujo de refrigerante en la chaqueta
Tjin  = u(3); %[K] temperatura de la chaqueta
Tinf  = u(4); %[K] temperatura ambiente
%Tmin = u(5); %[K] temperatura mosto
Talim = u(5); %[K] temperatura alimentacion

% Parámetros
D = Falim/V; %[h] tasa de dilución
O = 0.035;   %[g/L] concentración de saturación de oxígeno
Sin = 350;   %[g/L] concentración de glucosa en alimentación

kx1 = 0.49; %coeficiente de rendimiento de biomasa 1
kx2 = 0.05; %coeficiente de rendimiento de biomasa 2
kx3 = 0.72; %coeficiente de rendimiento de biomasa 3
ks1 = 1;    %coeficiente de rendimiento de sustrato 1
ks2 = 1;    %coeficiente de rendimiento de sustrato 2
kp2 = 0.48; %coeficiente de rendimiento de subproducto 2
kp3 = 1;    %coeficiente de rendimiento de subproducto 3

mu_s   = 3.5;    %[gS/gX/h]
mu_o   = 0.256;  %[gO2/gX/h]

Ks     = 0.1;    %[gS/L]
Kip    = 10;     % Kie en dewasme
Ko     = 0.0001; %[gO2/L]
Kp     = 0.1;    % Ke en dewasme
kos    = 0.3968; %[gO2/gS] kos=ko1 para levadura
kop    = 1.104;  %[gO2/gE] ko3 en dewasme
rs     = mu_s*S/(S+Ks);
rscrit = mu_o*O*Kip/(kos*(O+Ko)*(Kip+P));

Vj     = 3 ;         %[L]              volumen chaqueta OJO, cambiar
UAj    = 0.1 ;         %[J/(ºC*h)]       VALOR INVENTADO GG
U      = 1000;       %[W/m2K]          OJO SI []>50g/L U ES MUCHO MENOR A 1000
UAe    = 0.75;       %[W/K]
rhojcpj= 4.184;      %[J/(m^3*K)]
rhocp  = 500*4184;   %[J/(m^3*K)]
mf     = Falim;      %[L/h] = [kg/h]   DUDA DE SI INCORPORAR LA DENSIDAD DE GLUCOSA (DIEGO DIJO QUE ERA AGUA)
mj     = Fj;         %[L/h] = [kg/h]
cpf    = 4184;       %[J/kgK]
cpj    = 4184;       %[J/kgK]

% Calores
qalim  = mf*cpf*(Tm-Talim); % alimentación
qj     = mj*cpj*(Tjin-Tj);  % chaqueta
qe     = UAe*(Tm-Tinf);     % environment
qs     = 0;                 % agitador
qp     = 0;                 % levadura


% Estos r son los que generan los cambios en el metabolismo (crabtree,
% overflow, etc)
r1     = min(rs, rscrit)/ks1; 
r2     = max(0, rs-rscrit)/ks2; 
r3     = max(0, kos*(rscrit-rs)*P/(kop*(P+Kp)))/kp3;


% Ecuaciones de Estado (balances de masa)
dXdt  = (kx1*r1 + kx2*r2 + kx3*r3)*X - D*X;
dSdt  = -(ks1*r1+ks2*r2)*X + D*Sin - D*S;
dPdt  = (kp2*r2-kp3*r3)*X - D*P;
dVdt  = Falim;
dTmdt = qp - qalim - qj - qe + qs;                       %(dH/rhocp)*r   
dTjdt = (Fj/Vj)*(Tj-Tjin)+ ((UAj)*(Tm-Tj)/(Vj*rhojcpj)); %temperatura que sale de la chaqueta 

sys = [dXdt, dSdt, dPdt, dVdt, dTmdt, dTjdt];
end

function sys=mdlOutputs(~,x,~)

% Aqui se presenta el vector de salidas que tendra este macro, las cuales
% se ordenan para facilitar la colocacion de los bloques:

sys = x;

end