
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
sizes.NumInputs      = 4; % Numero de variables de entrada que el macro 
%                           aceptara. T entrada, T chaqueta, Flujo alim
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;


sys=simsizes(sizes); % Expresion que permite que el macro funcione bien. No 
% se debe borrar.

% Este es el SS o condiciones iniciales que pueden editar después de la parte a).
x0 = [0.4, 0.5, 0, 6.8, 298, 303]; 
%     X0, S0, E0, V0, Tmosto, Tchaqueta
% Dewasme: 0.4, 0.5, 0.8, 6.8

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
Talim = u(4); %[K] temperatura alimentacion

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

mu_s   = 3.5;                               %[gS/gX/h]
mu_o   = 0.256;                             %[gO2/gX/h]
Ks     = 0.1;                               %[gS/L]
Kip    = 10;                                %[gP/L] Kie en dewasme
Ko     = 0.0001;                            %[gO2/L]
Kp     = 0.1;                               %[gP/L] Ke en dewasme
kos    = 0.3968;                            %[gO2/gS] kos=ko1 para levadura
kop    = 1.104;                             %[gO2/gE] ko3 en dewasme
rs     = mu_s*S/(S+Ks);                     %[gS/gX*h]
ro     = mu_o*O*Kip/((O+Ko)*(Kip+P));       %[gO2/gX*h]
rscrit = ro/kos;                            %[gS/gX*h]
ro2    = 0.1184;                            %[gO2/gX*h]
Vj     = 10;                                %[L]             
UAj    = 3750*4184*5;                       %[J/(h*K)]
rhom   = 1080;                              %[kg/m3]
rhoj   = 1000;                              %[kg/m3]  
masam  = V*rhom/1000;                       %[kg]
mf     = Falim*1.13;                        %[kg/h]
mj     = Fj;                                %[kg/h]
cpf    = 4184;                              %[J/kgK]
cpm    = 4184;                              %[J/kgK]
cpj    = 4184;                              %[J/kgK]
dH     = 91.2*10^6/180;                     %[J/kg]

% Estos r son los que generan los cambios en el metabolismo (crabtree,
% overflow, etc)
r1     = min(rs, rscrit)/ks1; 
r2     = max(0, rs-rscrit)/ks2; 
r3     = max(0, kos*(rscrit-rs)*P/(kop*(P+Kp)))/kp3;

% Calores REVISAR UNIDADES DE TODOS
qalim  = Falim*(Tm-Talim)/V;                % alimentación [K/h]
qj     = UAj*(Tm - Tj)*1000/(V*rhom*cpm);   % chaqueta [K/h]
%[/h]*K/[L/m3]
qs     = 50*3600/(masam*cpm);                   % agitador [K/h]
qp     = dH*(ks1*r1+ks2*r2)*X/(rhom*cpm);                % levadura [K/h]
%[J/kg]*[gS/gX*h]*[gS/L] / [gO2/molO2]*[kg/m3]*[J/kg*K]
% taller: (dH*r/rhocp)
% K/h

% Ecuaciones de Estado (balances de masa y energía)
dXdt  = (kx1*r1 + kx2*r2 + kx3*r3)*X - D*X;
dSdt  = -(ks1*r1+ks2*r2)*X + D*Sin - D*S;
dPdt  = (kp2*r2-kp3*r3)*X - D*P;
if V < 20
    dVdt  = Falim;
else
    dVdt = 0;
    set_param(gcs, 'SimulationCommand', 'stop');
end   
dTmdt = -qalim + qp - qj + qs;  
dTjdt = (Fj/Vj)*(Tjin-Tj)+ UAj*(Tm-Tj)/(Vj*rhoj*cpj);

sys = [dXdt, dSdt, dPdt, dVdt, dTmdt, dTjdt];
end

function sys=mdlOutputs(~,x,~)

% Aqui se presenta el vector de salidas que tendra este macro, las cuales
% se ordenan para facilitar la colocacion de los bloques:

sys = x;

end