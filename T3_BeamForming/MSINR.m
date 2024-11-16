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
snap = 1024;              % 快拍数
num_trials = 100;         % 进行100次蒙特卡洛仿真
angle = -1/2 * pi : 0.001 : 1/2 * pi; % 扫描角度范围

% 累积方向图结果初始化
F_total_with_signal = zeros(1, length(angle));
F_total_without_signal = zeros(1, length(angle));

for trial = 1:num_trials
    % 信号生成
    signal = exp(1j * 2 * pi * f0 * (0:snap-1) / (2 * snap));
    INR_lin = 10.^(INR / 10) / 2;

    % 生成导向向量
    A = exp(1j * 2 * pi * n' * sin(deg2rad(theta)) * d / lambda);
    alpha0 = exp(1j * 2 * pi * n' * sin(deg2rad(theta0)) * d / lambda);
    xs = sqrt(10^(SNR / 10)) * alpha0 * signal; % 信号分量

    % 干扰生成
    kj = length(theta); % 干扰数量
    rs = zeros(kj, snap);
    for i = 1:kj
        for j = 1:snap
            rs(i, j) = sqrt(INR_lin(i)) * (randn(1) + 1j * randn(1)) * signal(j);
        end
    end

    % 噪声生成
    noise = (randn(N, snap) + 1j * randn(N, snap)) / sqrt(2);

    % 计算协方差矩阵（包含信号）
    Rs = 1 / snap * (xs * xs');          % 信号协方差矩阵
    J = A * rs;                           % 干扰矩阵
    Rin_with_signal = 1 / snap * ((xs +J + noise) * (xs +J + noise)'); % 包含信号协方差矩阵

    % 计算协方差矩阵（不包含信号）
    Rin_without_signal = 1 / snap * ((J + noise) * (J + noise)'); % 不包含信号协方差矩阵

    % MSINR准则下的最优权重 - 不包含信号
    [V_without_signal, D_without_signal] = eig(Rs, Rin_without_signal);
    [~, I_without_signal] = max(diag(D_without_signal));
    w_without_signal = V_without_signal(:, I_without_signal); % 最优权重向量（不包含信号）

    % MSINR准则下的最优权重 - 包含信号
    [V_with_signal, D_with_signal] = eig(Rs, Rin_with_signal);
    [~, I_with_signal] = max(diag(D_with_signal));
    w_with_signal = V_with_signal(:, I_with_signal); % 最优权重向量（包含信号）

    % 计算单次试验的自适应方向图
    F_without_signal = w_without_signal' * exp(1j * 2 * pi * n' * sin(angle) * d / lambda);
    F_with_signal = w_with_signal' * exp(1j * 2 * pi * n' * sin(angle) * d / lambda);

    % 累加方向图
    F_total_without_signal = F_total_without_signal + abs(F_without_signal);
    F_total_with_signal = F_total_with_signal + abs(F_with_signal);
end

% 取平均并转换为dB
F_avg_without_signal = db(F_total_without_signal / num_trials); 

F_avg_with_signal = db(F_total_with_signal / num_trials);

% 绘制结果
figure('Name', 'MSINR - Signal in Covariance Matrix Comparison');
plot(rad2deg(angle), F_avg_without_signal, 'b-', 'DisplayName', '不包含信号协方差矩阵');
hold on;
plot(rad2deg(angle), F_avg_with_signal, 'r--', 'DisplayName', '包含信号协方差矩阵');

% 图例和标题
legend('Location', 'Best');
title('MSINR（包含信号协方差矩阵 vs 不包含信号协方差矩阵）');
xlabel('角度 (度)');
ylabel('增益 (dB)');
grid on;
hold off;
