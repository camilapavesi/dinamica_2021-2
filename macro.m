%% Funcion de transferencia global
out = sim("mysim11.slx");
samplet = 1;
step    = 17;

t = out.data1(:,1);
t = t(step/samplet:end);
X = out.data1(:,2);
X = X(step/samplet:end);

F = out.input1;
F = F(step/samplet:end);

datos1 = iddata(X,F,samplet);

sys1 = tfest(datos1,2,2);

[num1 , den1] = tfdata(sys1);
num1 = cell2mat(num1);
den1 = cell2mat(den1);
num1 = num1/den1(end);
den1 = den1/den1(end);

plot(t,X,'-k','linewidth',1.5); hold on; grid on;
T1ft = out.data1(:,3);
T1ft = T1ft(step/samplet:end);
plot(t,T1ft,'--r','linewidth',1.5); hold off
xlabel('t [hr]')
ylabel('Tasa específica de crecimiento')
legend('μ','μ_{ft}','location','best')
title('Función de Transferencia μ/Falim')

time = out.mu_perturb(:,1);
mu_ref = out.mu_perturb(:,2);
mu_perturb = out.mu_perturb(:,3);
plot(time, mu_ref,'-k','linewidth',1.5); hold on; grid on;
plot(time, mu_perturb,'--r','linewidth',1.5); hold off
xlabel('t [hr]')
ylabel('Tasa específica de crecimiento')
legend('μ referencia','μ perturbado','location','best')
title('Perturbación de μ mediante un escalón en Falim')

%% Funcion de transferencia Tm Fj
out = sim("mysim11.slx");
samplet = 1;
step3   = 15;

t3 = out.data3(:,1);
t3 = t3(step3/samplet:end);
T3 = out.data3(:,2);
T3 = T3(step3/samplet:end);

F3 = out.input3;
F3 = F3(step3/samplet:end);

datos3 = iddata(T3,F3,samplet);

sys3 = tfest(datos3,1,0,1);

[num3 , den3] = tfdata(sys3);
num3 = cell2mat(num3);
den3 = cell2mat(den3);
num3 = num3/den3(end);
den3 = den3/den3(end);

%plot(t3,T3,'-k','linewidth',1.5); hold on; grid on;
%T3ft = out.data3(:,3);
%T3ft = T3ft(step3/samplet:end);
%plot(t3,T3ft,'--r','linewidth',1.5); hold off
%xlabel('t [hr]')
%label('Temperatura mosto [K]')
%legend('Tm','Tm_{ft}','location','best')
%title('Función de Transferencia Tm/Fj')

time = out.Tm_perturb(:,1);
Tm_ref = out.Tm_perturb(:,2);
Tm_perturb = out.Tm_perturb(:,3);
plot(time, Tm_ref,'-k','linewidth',1.5); hold on; grid on;
plot(time, Tm_perturb,'--r','linewidth',1.5); hold off
xlabel('t [hr]')
ylabel('Temperatura del mosto [K]')
legend('Tm referencia','Tm perturbado','location','best')
title('Perturbación de Tm mediante un escalón en Fj')

%% grafico

out = sim("mysim11.slx");
t = out.grafico1(:,1);
X = out.grafico1(:,2);
S = out.grafico1(:,3);
P = out.grafico1(:,4);
V = out.grafico1(:,5);
TM = out.grafico1(:,6);
TJ = out.grafico1(:,7);
plot(t,X,'-b','linewidth',1.5); hold on; grid on;
plot(t,S,'-m','linewidth',1.5);
plot(t,P,'-r','linewidth',1.5);
plot(t,V,'-k','linewidth',1.5); hold off
plot(t,TM,'-r','linewidth',1.5); hold on; grid on;
plot(t,TJ,'-b','linewidth',1.5); hold off
xlabel('t[h]')
ylabel('T(K)')
legend('TM', 'TJ')
title('Simulación balance de energía con escalón de Tj t=25h')