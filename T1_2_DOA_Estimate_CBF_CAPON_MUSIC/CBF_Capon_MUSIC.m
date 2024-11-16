clc;close all;clear;

%DOA algorithms : CBF_Capon_MUSIC
%Author:ChrisZhou 2024.10 @ SYSU

function [doa_cbf,doa_capon,doa_music,angleIndex] = myDOA(arrayNum,snapshotLength,SNR)
%Parameters:
c=3e8;%光速
fc=1.5e9;    %雷达载频
lambda=c/fc; %波长
d=lambda/2;  %阵元间距
theta_T=[10,20].'; %目标来向角
targetNum=length(theta_T);   %目标个数
%arrayNum=16;                 %阵元个数
arrayPos = 0:d:(arrayNum-1) * d; %阵元位置
sampleRate= 48e5;            %采样率
ts = 1/sampleRate;           %采样间隔
L = snapshotLength;          %快拍数

angleIndex = asin((-512:1:512-1)/512) * 180 /pi; %角度纲量
angleIndexLength = length(angleIndex);
timeIndex = (0:1:L-1)*ts;                        %时间纲量

%Signal Genration

As = zeros(arrayNum,targetNum);  % A(theta)
S_signal = zeros(targetNum,L);   % s(n)
%分别定义A(theta) 与 s(n)
for i=1:1:targetNum
    As(:,i) = exp(-1j*(2*pi/lambda).*arrayPos.'*(sind(theta_T(i))));
    S_signal(i,:) = exp(1j*2*pi*fc*timeIndex);
    fc=fc+1000;  %假设后一个目标回波信号的频率高1000Hz
end

S_signal = As * S_signal;      % X(n) = A(theta) * S(n)
%SNR = 20;                      % 设定信噪比
S_signal = awgn(S_signal,SNR); % 添加高斯白噪声 即x(n) = x(n) + N(n)

%定义权矢量矩阵 W
W = zeros(arrayNum,angleIndexLength);
for i=1:angleIndexLength
    W(:,i)=exp(-1j*(2*pi/lambda).*arrayPos.'*sind(angleIndex(i)));
end

%信号协方差矩阵的最大似然估计（标准化后的自相关矩阵）
R_hat = (S_signal * conj(S_signal).' ) / L;
R_hatInv = inv(R_hat);
%DOA Algorithms

%CBF方法
doa_cbf=zeros(angleIndexLength,1);

%P(w)=W^H R W
for i=1:1:angleIndexLength
    doa_cbf(i)=conj(W(:,i).') * R_hat * W(:,i);
end


%Capon方法
doa_capon = zeros(angleIndexLength,1);

%P(w) = 1 / A^H R^-1 
for i=1:1:angleIndexLength
    doa_capon(i)= 1 / (conj(W(:,i).')* R_hatInv * W(:,i));
end

%MUSIC方法
doa_music = zeros(angleIndexLength,1);
[EV,D] = eig(R_hat);    %将协方差矩阵进行特征分解
EVA = diag(D)';         %提取特征值矩阵对角线，并转为一行
[EVA,I] = sort(EVA);    %将特征值从小到大排序
EV = fliplr(EV(:,I));   %对应特征矢量的排序
Un = EV(:,targetNum+1:arrayNum);
%P(w) = 1/ A*Un*Un^H*A
for i = 1:1:length(angleIndex)
    doa_music(i)=1/(conj(W(:,i).') * Un * Un' * W(:,i));
end



% ESPRIT方法
% 首先，对信号协方差矩阵进行特征分解
[EV, D] = eig(R_hat);
EVA = diag(D)';
[EVA, I] = sort(EVA); % 特征值排序
EV = fliplr(EV(:, I)); % 特征矢量排序

% 选择信号子空间和噪声子空间
Un = EV(:, 1:targetNum); % 信号子空间
Un_plus = Un(2:end, :); % 上部分
Un_minus = Un(1:end-1, :); % 下部分

% 计算phi矩阵
Phi = pinv(Un_minus) * Un_plus;

% 计算特征值
[~, D_phi] = eig(Phi);
eigenvalues = diag(D_phi);
angles = angle(eigenvalues); % 提取相位

% 计算估计的DOA
doa_esprit = asin(lambda * angles / (2 * pi * d)) * (180/pi); % 转换为角度

% 显示结果
disp('DOA Estimates using ESPRIT:');
disp(doa_esprit);


end


%%plot

%Case 1 :SNR = 20dB时，估计 R 快拍数 = 500，阵元数为 8、16、24三种情况；
[c8_cbf,c8_capon,c8_music,angleIndex] = myDOA(8,500,20);
[c16_cbf,c16_capon,c16_music] = myDOA(16,500,20);
[c24_cbf,c24_capon,c24_music] = myDOA(24,500,20);

No=1;
figure(No);
%title('xibaloma');
%theta = linspace(-90, 90, 1024);

subplot(311);
plot(angleIndex,db(c8_cbf/max(c8_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(c8_capon/max(c8_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(c8_music/max(c8_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=8 Snapshot=500 SNR = 20');
xlim([-90,90]);
legend;
grid on;

subplot(312);
plot(angleIndex,db(c16_cbf/max(c16_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(c16_capon/max(c16_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(c16_music/max(c16_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=16 Snapshot=500 SNR = 20');
xlim([-90,90]);
legend;
grid on;

subplot(313);
plot(angleIndex,db(c24_cbf/max(c24_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(c24_capon/max(c24_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(c24_music/max(c24_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=24 Snapshot=500 SNR = 20');
xlim([-90,90]);
legend;
grid on;

%Case 2 : 估计 R 快拍数 = 500，阵元数为 16，SNR 为 5dB、10dB 和 20dB 三种情况；
No=No+1;
figure(No);

[n5_cbf,n5_capon,n5_music] = myDOA(16,500,5);
[n10_cbf,n10_capon,n10_music] = myDOA(16,500,10);
[n20_cbf,n20_capon,n20_music] = myDOA(16,500,20);

subplot(311);
plot(angleIndex,db(n5_cbf/max(n5_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(n5_capon/max(n5_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(n5_music/max(n5_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=16 Snapshot=500 SNR = 5');
xlim([-90,90]);
legend;
grid on;

subplot(312);
plot(angleIndex,db(n10_cbf/max(n10_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(n10_capon/max(n10_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(n10_music/max(n10_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=16 Snapshot=500 SNR = 10');
xlim([-90,90]);
legend;
grid on;

subplot(313);
plot(angleIndex,db(n20_cbf/max(n20_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(n20_capon/max(n20_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(n20_music/max(n20_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=24 Snapshot=500 SNR = 20');
xlim([-90,90]);
legend;
grid on;

%Case 3 :SNR = 20dB时，阵元数为 16，估计 R 快拍数为：10、20、50和 100 四种情况

No=No+1;
figure(No);

[s10_cbf,s10_capon,s10_music] = myDOA(16,10,20);
[s20_cbf,s20_capon,s20_music] = myDOA(16,20,20);
[s50_cbf,s50_capon,s50_music] = myDOA(16,50,20);
[s100_cbf,s100_capon,s100_music] = myDOA(16,100,20);

subplot(411);
plot(angleIndex,db(s10_cbf/max(s10_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(s10_capon/max(s10_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(s10_music/max(s10_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=16 Snapshot=10 SNR = 20');
xlim([-90,90]);
legend;
grid on;

subplot(412);
plot(angleIndex,db(s20_cbf/max(s20_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(s20_capon/max(s20_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(s20_music/max(s20_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=16 Snapshot=20 SNR = 20');
xlim([-90,90]);
legend;
grid on;

subplot(413);
plot(angleIndex,db(s50_cbf/max(s50_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(s50_capon/max(s50_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(s50_music/max(s50_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=16 Snapshot=50 SNR = 20');
xlim([-90,90]);
legend;
grid on;

subplot(414);
plot(angleIndex,db(s100_cbf/max(s100_cbf)),'b','DisplayName', 'CBF','LineWidth',1);
hold on;
plot(angleIndex,db(s100_capon/max(s100_capon)),'r','DisplayName', 'Capon','LineWidth',1);
hold on;
plot(angleIndex,db(s100_music/max(s100_music)),'g','DisplayName', 'Music','LineWidth',1);
xlabel('Angle(°)');ylabel('dB');title('M=16 Snapshot=100 SNR = 20');
xlim([-90,90]);
legend;
grid on;