clc; clear;

% 参数设置
N = 16;           % 阵元数
M = 2;           % 信号源数 (已知的信源数)
d = 0.5;         % 阵元间距 (lambda/2)
theta_actual = [30, 10]; % 实际信号DOA
snr = 20;        % 信噪比
lambda = 1;      % 信号波长
T = 200;         % 信号采样数
num_trials = 1000;  % 蒙特卡罗仿真实验次数
angles = -90:0.1:90;  % 角度范围    

%% 存储ESPRIT和MLE结果的变量
doa_histogram_esprit = zeros(1, length(angles));  % ESPRIT的二值化结果
P_mle_avg = zeros(1, length(angles));  % MLE的平均似然结果

%% 定义函数部分
% 定义ESPRIT算法的旋转不变性矩阵
function doa_est = esprit_estimation(Rxx, N, M, d)
    % 特征值分解
    [U, ~] = eig(Rxx);
    Us = U(:, N-M+1:N);  % 信号子空间
    
    % 子阵列划分
    Us1 = Us(1:N-1, :);
    Us2 = Us(2:N, :);
    
    % 旋转不变性
    Phi = pinv(Us1) * Us2;
    
    % 求解特征值，计算DOA估计
    eig_vals = eig(Phi);
    doa_est = asin(angle(eig_vals) / (2*pi*d)) * 180/pi;
end

% 定义MLE算法中的投影矩阵计算
function P_A_perp = projection_matrix(A_theta)
    P_A = A_theta * inv(A_theta' * A_theta) * A_theta';  % 投影矩阵 P_A
    P_A_perp = eye(size(P_A)) - P_A;                     % 正交投影矩阵 P_A_perp
end

% 定义MLE的似然函数
function L = mle_likelihood(theta, X, N, d, lambda, R_hat)
    A_theta = exp(1j*2*pi*d*(0:N-1)'*sin(theta*pi/180)/lambda);  % 构造导向矩阵
    P_A_perp = projection_matrix(A_theta);                      % 计算正交投影矩阵
    L = trace(P_A_perp * R_hat);                                % 似然函数的简化形式
end

%% 蒙特卡罗仿真
for trial = 1:num_trials
    % 生成导向矩阵
    A = zeros(N, M);
    for k = 1:M
        A(:, k) = exp(1j*2*pi*d*(0:N-1)'*sin(theta_actual(k)*pi/180)/lambda);
    end
    
    % 生成信号
    S = randn(M, T) + 1j*randn(M, T);  % 信号矩阵
    X = A * S;                         % 阵列接收信号
    
    % 加入噪声
    X = awgn(X, snr, 'measured');
    
    % 协方差矩阵
    Rxx = (X * X') / T;
    R_hat = Rxx;  % 对于MLE算法，协方差矩阵也需要
    
    %% ESPRIT算法
    doa_est_esprit = esprit_estimation(Rxx, N, M, d);
    
    % 将DOA估计值映射到角度范围上并更新直方图
    for k = 1:length(doa_est_esprit)
        [~, idx] = min(abs(angles - doa_est_esprit(k)));  % 找到最接近的角度索引
        doa_histogram_esprit(idx) = doa_histogram_esprit(idx) + 1; % 对应角度的计数增加
    end
    
    %% MLE算法
    P_mle = zeros(1, length(angles));
    for k = 1:length(angles)
        P_mle(k) = mle_likelihood(angles(k), X, N, d, lambda, R_hat);
    end
    
    % 累加结果，用于后续计算平均值
    P_mle_avg = P_mle_avg + P_mle;
end

%% 处理ESPRIT算法的结果 (二值化)
doa_binary_esprit = doa_histogram_esprit > 0;

%% 处理MLE算法的结果 (取平均)
P_mle_avg = P_mle_avg / num_trials;

%% 绘制对比图
figure;

% 绘制ESPRIT算法结果
subplot(2,1,1);
stem(angles, doa_binary_esprit, 'b', 'LineWidth', 2);
xlabel('DOA角度 (度)');
ylabel('检测结果 (1=有信号, 0=无信号)');
title('ESPRIT算法');
grid on;

% 绘制MLE算法结果
subplot(2,1,2);
plot(angles, 10*log10(real(1./P_mle_avg)), 'r', 'LineWidth', 2);
xlabel('DOA角度 (度)');
ylabel('平均谱强度 (dB)');
title('MLE算法');
grid on;

