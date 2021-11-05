
out = sim("mysim11.slx");

samplet = 1;

step    = 3;

t = out.data1(:,1);
t = t(step/samplet:end);
X = out.data1(:,2);
X = X(step/samplet:end);

F = out.input1;
F = F(step/samplet:end);

datos1 = iddata(X,F,samplet);

sys1 = tfest(datos1,2,0);

[num1 , den1] = tfdata(sys1);
num1 = cell2mat(num1);
den1 = cell2mat(den1);
num1 = num1/den1(end);
den1 = den1/den1(end);

