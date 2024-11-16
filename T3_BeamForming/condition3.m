clc;
clear;
close all;
M = 16; % 阵元数
snapshots = 200; % 快拍数
thetas = 0; % 信号入射角度
d = 0.5; % 阵元间距（半波长）
thetai = [-20, 45]; % 干扰入射角度
n = [0:M-1]';

vs = exp(-j * 2 * pi * d * n * sind(thetas)); % 信号方向向量
vi = exp(-j * 2 * pi * d * n * sind(thetai)); % 干扰方向向量

f = 16000; % 载波频率
t = [0:1:snapshots-1] / 200;

snr = 10; % 信噪比
inr = 10; % 干噪比

rvar = 1; % 信号功率
de = rvar * exp(j * 2 * pi * f * t); % 期望信号

xs = vs * de; % 构造有用信号

num_trials = 100; % 仿真次数
B1_avg = zeros(1, 18001); % 初始化波束图累加数组
B2_avg = zeros(1, 18001); % 初始化波束图累加数组
for trial = 1:num_trials
    temp = zeros(M, 1); % 重置每次试验的权重向量累加
    for g = 1:1
        xi1 = [randn(length(thetai), snapshots) + j * randn(length(thetai), snapshots)];
        xi = vi * xi1; % 构造干扰信号
        noise = 10^(-snr / 10) * [randn(M, snapshots) + j * randn(M, snapshots)]; % 噪声

        %目标信号
        X = xs;

        % 计算协方差矩阵
        R = X * X' / snapshots;

        % 计算噪声和干扰的协方差矩阵
        Rn = (xi + noise) * (xi + noise)' / snapshots;

        % 广义特征值和特征向量
        [V, D] = eig(R, Rn);
        [a, b] = max(diag(D));

        wop2 = V(:, b);
        temp(:, g) = wop2;
    end

    wop3 = mean(temp, 2);

    % 扫描角度范围并累加波束图
    sita = -90:0.01:90;
    v = exp(-j * pi * n * sind(sita));
    B1 = abs(wop3' * v);

    % 累加每次仿真的波束图
    B1_avg = B1_avg + B1;
end


for trial = 1:num_trials
    temp = zeros(M, 1); % 重置每次试验的权重向量累加
    for g = 1:1
        xi1 = [randn(length(thetai), snapshots) + j * randn(length(thetai), snapshots)];
        xi = vi * xi1; % 构造干扰信号
        noise = 10^(-snr / 10) * [randn(M, snapshots) + j * randn(M, snapshots)]; % 噪声

        % 目标信号
        X = xs;

        % 计算协方差矩阵，包含目标信号、干扰和噪声
        R = X * X' / snapshots;

        % 计算噪声和干扰的协方差矩阵
        Rn = (xs + xi + noise) * (xs + xi + noise)' / snapshots;

        % 广义特征值和特征向量
        [V, D] = eig(R, Rn);
        [a, b] = max(diag(D));

        wop2 = V(:, b);
        temp(:, g) = wop2;
    end

    wop3 = mean(temp, 2);

    % 扫描角度范围并累加波束图
    sita = -90:0.01:90;
    v = exp(-j * pi * n * sind(sita));
    B2 = abs(wop3' * v);

    % 累加每次仿真的波束图
    B2_avg = B2_avg + B2;
end


% 取平均波束图
B2_avg = B2_avg / num_trials;

% 绘制波束图
figure(1)
plot(sita, db(B1_avg));
hold on;
plot(sita, db(B2_avg));
title('MSINR波束图 - 协方差矩阵包含目标信号（多次仿真平均）')
xlabel("角度/degree")
ylabel("波束图/dB")
grid on
axis([-90 90 -80 0]);
hold off;
