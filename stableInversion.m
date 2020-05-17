r = dis.signals.values;
figure;
plot(r);
%%
z= tf('z',1/5000);
inversePlant = 1/GpDis;
forwardOrder = numel(zero(inversePlant)) - numel(pole(inversePlant));
F =ss(inversePlant*z^-forwardOrder);
A = F.A;
[V,D] = eig(A);
[sortedD,index] = sort(abs(diag(D)));
Q = V(:,index);
D = D(:,index);
D = D(index,:);
% index = sortedD < 1;
% Qs = Q(:,index);
% Ds = sortedD(index);
% Qu = Q(:,~index);
% Du = sortedD(~index);
temp = sortedD < 1;
As = D(temp,temp);
Au = D(~temp,~temp);
B = inv(Q)* F.B;
Bs = B(temp);
Bu = B(~temp);
num = numel(r);
xs = cell(num,1);
xu = cell(num,1);
xs{1} = zeros(size(As,1),1);
xu{num} = zeros(size(Au,1),1);
%%
for i = 2:num
    xs{i} = As*xs{i-1} + Bs * r(i-1);
end
invAu = inv(Au);
for i = num-1:-1:1
   xu{i} = invAu * ( xu{i+1} - Bu * r(i) ); 
end
%%
C = F.C;
DD = F.D;
CQ = C * Q;
for i = 1:num
   y(i) = CQ * [xs{i};xu{i}] + DD * r(i); 
    
end
figure;plot(abs(y(1:500)));
%%
ffSignal = dis;
index = dis.time <= 0.16;
temp = ffSignal.signals.values;
temp = zeros(size(temp));
% temp(1:500) = abs(y(1:500));
temp(index) = real(y(index));
temp = circshift(temp,-forwardOrder);
ffSignal.signals.values = temp;
figure;plot(ffSignal.signals.values);
