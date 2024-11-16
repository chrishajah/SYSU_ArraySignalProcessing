clc;
clear;
close all;
M = 16; % 阵元数
thetas = 0; % 信号入射角度
d = 0.5; % 阵元间距（半波长）
thetai = [5]; % 干扰入射角度
n = [0:M-1]';
vs = exp(-j * 2 * pi * d * n * sind(thetas)); % 信号方向向量
vi = exp(-j * 2 * pi * d * n * sind(thetai)); % 干扰方向向量
f = 16000; % 载波频率
snr = 10; % 信噪比
inr = 10; % 干噪比
rvar = 1; % 信号功率

snapshots_list = [50, 100, 200, 500]; % 不同的快拍数
colors = ['r', 'g', 'b', 'k']; % 用不同颜色区分不同快拍数
num_trials = 100; % 多次仿真次数

figure;
hold on;
for idx = 1:length(snapshots_list)
    snapshots = snapshots_list(idx);
    t = [0:1:snapshots-1] / 200;
    de = rvar * exp(j * 2 * pi * f * t); % 期望信号
    xs = vs * de; % 构造有用信号
    B1_avg = zeros(1, 18001); % 存储多个仿真结果的平均
    
    for trial = 1:num_trials
        xi1 = [randn(length(thetai), snapshots) + j * randn(length(thetai), snapshots)];
        xi = vi * xi1; % 构造干扰信号
        noise = 10^(-snr / 10) * [randn(M, snapshots) + j * randn(M, snapshots)]; % 噪声
        X = xs + xi + noise;
        R = X * X' / snapshots; % 包含目标信号的协方差矩阵
        Rn = (xi + noise) * (xi + noise)' / snapshots; % 干扰和噪声协方差矩阵
        [V, D] = eig(R, Rn);
        [a, b] = max(diag(D));
        wop3 = V(:, b);

        % 扫描角度范围并累加波束图
        sita = -90:0.01:90;
        v = exp(-j * pi * n * sind(sita));
        B1 = abs(wop3' * v);
        B1_avg = B1_avg + B1; % 累加每次仿真的波束图
    end
    
    % 取平均并绘制波束图
    B1_avg = B1_avg / num_trials;
    plot(sita, 20 * log10(B1_avg / max(B1_avg)), colors(idx), 'DisplayName', ['Snapshots = ' num2str(snapshots)]);
end

title('MSINR波束图 - 不同快拍数对比（多次仿真平均）')
xlabel("角度/degree")
ylabel("波束图/dB")
legend('show');
grid on;
axis([-90 90 -80 0]);
hold off;
