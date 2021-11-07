%% Funcion de transferencia global
out = sim("mysim11.slx");
samplet = 1;
step    = 15;

t = out.data1(:,1);
t = t(step/samplet:end);
X = out.data1(:,2);
X = X(step/samplet:end);

F = out.input1;
F = F(step/samplet:end);

datos1 = iddata(X,F,samplet);

sys1 = tfest(datos1,2,1);

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

%% Funcion de transferencia Tm Fj
out = sim("mysim11.slx");
samplet = 1;
step3   = 13;

t3 = out.data3(:,1);
t3 = t3(step3/samplet:end);
T3 = out.data3(:,2);
T3 = T3(step3/samplet:end);

F3 = out.input3;
F3 = F3(step3/samplet:end);

datos3 = iddata(T3,F3,samplet);

sys3 = tfest(datos3,2,0);

[num3 , den3] = tfdata(sys3);
num3 = cell2mat(num3);
den3 = cell2mat(den3);
num3 = num3/den3(end);
den3 = den3/den3(end);

plot(t3,T3,'-k','linewidth',1.5); hold on; grid on;
T3ft = out.data3(:,3);
T3ft = T3ft(step3/samplet:end);
plot(t3,T3ft,'--r','linewidth',1.5); hold off
xlabel('t [hr]')
ylabel('Temperatura mosto [K]')
legend('Tm','Tm_{ft}','location','best')
title('Función de Transferencia Tm/Tj')

