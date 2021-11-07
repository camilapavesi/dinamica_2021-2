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

sys1 = tfest(datos1,3,1);

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
ylabel('Temperatura mosto [K]')
legend('TM','TM_{ft}','location','best')
title('Función de Transferencia TM/Tj')

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
legend('TM','TM_{ft}','location','best')
title('Función de Transferencia TM/Tj')

%% Funcion de transferencia Tj Fj
out = sim("mysim11.slx");
samplet = 1;
step4   = 13;

t4 = out.data4(:,1);
t4 = t4(step4/samplet:end);
T4 = out.data4(:,2);
T4 = T4(step4/samplet:end);

F4 = out.input4;
F4 = F4(step4/samplet:end);

datos4 = iddata(T4,F4,samplet);

sys4 = tfest(datos4,2,1);

[num4 , den4] = tfdata(sys4);
num4 = cell2mat(num4);
den4 = cell2mat(den4);
num4 = num4/den4(end);
den4 = den4/den4(end);

%% Funcion de transferencia Tj Falim
out = sim("mysim11.slx");
samplet = 1;
step5   = 13;

t5 = out.data5(:,1);
t5 = t5(step5/samplet:end);
T5 = out.data5(:,2);
T5 = T5(step5/samplet:end);

F5 = out.input5;
F5 = F5(step5/samplet:end);

datos5 = iddata(T5,F5,samplet);

sys5 = tfest(datos5,2,0);

[num5 , den5] = tfdata(sys5);
num5 = cell2mat(num5);
den5 = cell2mat(den5);
num5 = num5/den5(end);
den5 = den5/den5(end);