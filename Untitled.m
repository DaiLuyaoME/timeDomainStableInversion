
%%
FF = tf(F);
z = tf('z',1/5000);
FF = FF * z^-1;
y = stable_inversion(FF,r,1032);
figure;plot(y);
%%
time  = Err.time;
error = Err.signals.values;
error = error(time<0.12);
time = time(time<0.12);
figure;
plot(time,error * 1e9);
