%% Code section A: plot tracking error
error1 = Err1.signals.values;
error2 = Err2.signals.values;
snapSignal = snap.signals.values;
time = snap.time;

error = [error1,error2] * 1e9;
figure;
h = plot(time,error,'linewidth',2);
h(1).DisplayName = 'error with model uncertainty';
h(2).DisplayName = 'error with accurate model';
h(2).Color = [0.9290    0.6940    0.1250];
hold on;
ratio = max(max(error)) / max(snapSignal);
plot(time,ratio * snapSignal,'linewidth',1,'displayname','scaled snap','color','r','linestyle','--');
legend show;
xlabel('time (s)');
ylabel('error (nm)');
set(gca,'fontsize',13);


