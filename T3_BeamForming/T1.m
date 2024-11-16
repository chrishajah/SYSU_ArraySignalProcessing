%% MSINR准则下使用带功率的协方差公式
clc;clear;close;
N = 16;
n = 0:N-1;
f0 = 10e9;  % 雷达工作频率10GHz
c = 3e8;
lambda = c/f0;
d = lambda/2;   % 阵元距离
sita0 = 0;      % 期望信号方向
signal_e = 10^5;    % 信号能量
sita1 = deg2rad(45);
sita2 = deg2rad(-20);
v = deg2rad(-20);
sita = [sita1 sita2];     % 干扰方向
INR1 = 10^(40/10);
INR2 = 10^(40/10);
%INR3 = 10^(50/10);
%INR = [INR1 INR2 INR3];         % 干扰强度
INR = [INR1 INR2];   
A = exp(1j*2*pi*n'*sin(sita)*d/lambda);     % 干扰信号的导向矢量
alpha0 = exp(1j*2*pi*n'*sin(sita0)*d/lambda);   % 期望信号的导向矢量值
Rin = zeros(N,N);
for i = 1:length(sita)
    Rin = Rin + INR(i)*A(:,i)*A(:,i)';          % 干噪协方差矩阵
end
Rin = Rin + eye(N);                         % 加入噪声
RS = signal_e*(alpha0*alpha0');
[V,D] = eig(RS,Rin);        % 得到广义特征值与特征向量
[~,I] = max(diag(D));
w = V(:,I);              % 最大特征值所对应的特征向量就是最优权矢量
angle = -5/9*pi:0.01:5/9*pi;
A1 = exp(1j*2*pi*n'*sin(angle)*d/lambda);
CBF = db(abs(sum(A1)).^2/N/N)/2;          % 权重全为1时的增益方向图
MSINR = abs(w'*A1);
MSINR = db(MSINR/max(MSINR));

CBF(CBF < -80) = -80;
MSINR(MSINR < -80) = -80;

[~,ind1] = min(abs(angle-sita1));
[~,ind2] = min(abs(angle-sita2));
%[~,ind3] = min(abs(angle-sita3));

figure('Name','最优化波束形成');
plot(rad2deg(angle),CBF,'LineStyle','--');
hold on
plot(rad2deg(angle),MSINR);
plot(rad2deg(sita1), MSINR(ind1),'ro', 'MarkerFaceColor', 'r');
text(rad2deg(sita1), MSINR(ind1), '\leftarrow干扰1');
plot(rad2deg(sita2), MSINR(ind2),'ro', 'MarkerFaceColor', 'r');
text(rad2deg(sita2), MSINR(ind2), '\leftarrow干扰2');
%plot(rad2deg(sita3), MSINR(ind3),'ro', 'MarkerFaceColor', 'r');
%text(rad2deg(sita3), MSINR(ind3), '\leftarrow干扰3');
xlabel('角度');
ylabel('归一化功率增益');
legend('CBF','MSINR');
hold off
