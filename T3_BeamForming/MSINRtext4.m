clc; clear; close all;

% 参数设置
N = 16;               % 阵元数
n = 0:N-1;            % 阵元索引
f0 = 16e9;            % 雷达工作频率16GHz
c = 3e8;              % 光速
lambda = c / f0;      % 波长
d = lambda / 2;       % 阵元间距

% 信号和干扰的方向、信噪比和干扰噪比
INR = [10, 10, 10];       % 干扰信噪比(dB)
SNR = 10;                 % 信号信噪比(dB)
theta = [-40, 20, 50];    % 干扰方向
theta0 = 0;               % 信号方向

% 设置不同快拍数
snap_values = [200, 500, 1000];
num_trials = 100;          % 进行100次蒙特卡罗仿真
angle = -5/9 * pi : 0.001 : 5/9 * pi; % 扫描角度范围

% 颜色设置
colors = ['b', 'r', 'g']; % 蓝色、红色、绿色表示不同快拍数

figure('Name', 'MSINR - Monte Carlo Average for Different Snapshots');

for s = 1:length(snap_values)
    snap = snap_values(s);  % 当前的快拍数
    F_total = zeros(1, length(angle)); % 累积方向图结果

    for trial = 1:num_trials
        % 信号生成
        signal = exp(1j * 2 * pi * f0 * (0:snap-1) / (2 * snap));
        INR = 10.^(INR / 10) / 2;

        % 生成导向向量
        A = exp(1j * 2 * pi * n' * sin(deg2rad(theta)) * d / lambda);
        alpha0 = exp(1j * 2 * pi * n' * sin(deg2rad(theta0)) * d / lambda);
        xs = sqrt(10^(SNR / 10)) * alpha0 * signal; % 信号分量

        % 干扰生成
        kj = length(theta); % 干扰数量
        rs = zeros(kj, snap);
        for i = 1:kj
            for j = 1:snap
                rs(i, j) = sqrt(INR(i)) * (randn(1) + 1j * randn(1)) * signal(j);
            end
        end

        % 噪声生成
        noise = (randn(N, snap) + 1j * randn(N, snap)) / sqrt(2);

        % 计算协方差矩阵
        Rs = 1 / snap * (xs * xs');         % 信号协方差矩阵
        J = A * rs;                          % 干扰矩阵
        Rin = Rs + 1 / snap * (J * J' + noise * noise'); % 干扰+噪声协方差矩阵

        % MSINR准则下的最优权重
        [V, D] = eig(Rs, Rin);
        [~, I] = max(diag(D));
        w = V(:, I); % 最优权重向量

        % 计算单次试验的自适应方向图
        F = w' * exp(1j * 2 * pi * n' * sin(angle) * d / lambda); 
        F_total = F_total + abs(F); % 累加方向图
    end

    % 取平均并转换为dB
    F_avg = db(F_total / num_trials / max(F_total)); % 方向图的平均值

    % 绘制平均结果
    plot(rad2deg(angle), F_avg, 'Color', colors(s), 'DisplayName', sprintf('快拍数 = %d', snap));
    hold on;
end

% 图例和标题
legend('Location', 'Best');
title('MSINR（不同快拍数下的蒙特卡洛仿真结果）');
xlabel('角度 (度)');
ylabel('增益 (dB)');
grid on;
hold off;
