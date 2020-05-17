m1 = 5;
m2 = 19;
m3 = 3;
m = m1+m2+m3;
f1 = 690;
f2 = 800;
f3 = 900;

w1 = f1 * 2 * pi;
w2 = f2 * 2 * pi;
w3 = f3 * 2 * pi;
z1 = 0.03;
z2 = 0.02;
z3 = 0.025;

s = tf('s');

G1 = 1/(m*s^2);
G2 = w1*w1/(s^2+2*z1*w1*s + w1*w1);
G3 = w2*w2/(s^2+2*z2*w2*s + w2*w2);
G4 = (s^2+2*z3*w3*s + w3*w3)/(w3*w3);
Gp3 = G1 * G2 * G3 * G4;
figure;bodeplot(Gp,Gp3);
%% delay factor
delayCount = 1;
s = tf('s');
delayModel = exp(-delayCount*Ts*s);
delayModel = pade(delayModel,2);
GpWithDelay3 = Gp3 * delayModel;
GpDis3 = c2d(GpWithDelay3,Ts,'zoh');
figure;
bodeplot(GpDis3,GpDis);
figure;
pzmap(GpDis3,GpDis);