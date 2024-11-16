clc;
clear;
close all;
M=16; %阵元数
snapshots=200; %快拍数
thetas=0; %信号入射角度
d = 0.5;%阵元间距（半波长）
thetai=[45]; %干扰入射角度
n=[0:M-1]';

vs=exp(-j*2*pi*d*n*sind(thetas)); %信号方向向量
vi=exp(-j*2*pi*d*n*sind(thetai)); %干扰方向向量

f=16000; %载波频率
t=[0:1:snapshots-1]/200;

snr=10; %信噪比
inr=10; %干噪比

rvar=1;%信号功率
de=rvar*exp(j*2*pi*f*t); %期望信号

xs=vs*exp(j*2*pi*f*t); %构造有用信号

mont=100;
temp = zeros(M,mont);
for g=1:mont
    xi1=[randn(length(thetai),snapshots)+j*randn(length(thetai),snapshots)];

    % xi=sqrt(10^(inr/10))*vi*[randn(length(thetai),L)+j*randn(length(thetai),L)]; %构造干扰信号

    xi=vi*xi1; %构造干扰信号

    noise=10^(-snr/ 10)*[randn(M,snapshots)+j*randn(M,snapshots)]; %噪声

    X=xi+noise; %干扰加噪声

    A=[vs vi];% 信号接收阵列

    St=[xs;xi];%接收信号


    R=xs*xs'/snapshots; %构造协方差矩阵
    Rn=X*X'/snapshots; %噪声协方差矩阵

    [V,D]=eig(R,Rn);%广义特征值和特征向量
    [a,b]=max(diag(D));

    wop2=V(:,b);
    temp(:,g)=wop2;
end

wop3 = mean(temp,2);

sita = -90:0.01:90;
v=exp(-j*pi*n*sind(sita));
B1=abs(wop3'*v);

figure(1)
plot(sita,20*log10(B1/max(B1)),'k');
title('MSINR波束图')
xlabel("角度/degree")
ylabel("波束图/dB")
grid on
axis([-90 90 -80 0]);
hold off