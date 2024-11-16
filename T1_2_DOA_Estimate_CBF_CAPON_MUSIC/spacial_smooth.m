% 参数设置
M = 8; % 阵列元件数量
L = 3; % 子阵列数量（用于空间平滑）
K = M - L + 1; % 每个子阵列的长度
K_signal = 2; % 相干信号源数量
d = 0.5; % 相邻传感器之间的距离（以波长为单位）
angles = [80, 60]; % 信号源的真实波达方向（以度为单位）

% 生成信号
Nsnapshot = 10; % 快拍数量
S1 = exp(1j * 2 * pi * rand(1, Nsnapshot));
S = repmat(S1, K_signal, 1); % 生成相干信号源

% 阵列方向向量
A = exp(-1j * 2 * pi * d * (0:M-1).' * sind(angles));

% 生成具有相干信号的阵列数据
X = A * S;

% 增强相干性
X(:, 2:end) = X(:, 1:end-1) + 0.99 * X(:, 2:end);

% 添加噪声
SNR_db = 10; % 信噪比（以dB为单位）
SNR = 10^(SNR_db / 10);
signal_power = mean(abs(X(:)).^2);
noise_power = signal_power / SNR;
X = X + sqrt(noise_power/2) * (randn(size(X)) + 1j * randn(size(X)));

% ============ 未使用空间平滑的MUSIC算法 ============
R_no_smooth = (1/Nsnapshot) * (X * X'); % 完整阵列协方差矩阵

% 特征值分解
[Evec_no_smooth, Eval_no_smooth] = eig(R_no_smooth);
[~, idx_no_smooth] = sort(diag(Eval_no_smooth), 'descend');
Evec_no_smooth = Evec_no_smooth(:, idx_no_smooth);
En_no_smooth = Evec_no_smooth(:, K_signal+1:end); % 噪声子空间

% 使用MUSIC算法进行DOA估计（未使用空间平滑）
theta_scan = -90:0.1:90;
P_MUSIC_no_smooth = zeros(size(theta_scan));
for i = 1:length(theta_scan)
    a = exp(-1j * 2 * pi * d * (0:M-1).' * sind(theta_scan(i)));
    P_MUSIC_no_smooth(i) = 1 / (a' * (En_no_smooth * En_no_smooth') * a);
end

% 绘制空间谱（未使用空间平滑）
P_MUSIC_no_smooth = 10*log10(abs(P_MUSIC_no_smooth) / max(abs(P_MUSIC_no_smooth)));
figure;
plot(theta_scan, P_MUSIC_no_smooth, 'LineWidth', 2);
grid on;
xlabel('角度 (度)');
ylabel('空间谱 (dB)');
title('未使用空间平滑的MUSIC DOA估计');

% ============ 使用空间平滑的MUSIC算法 ============
% 空间平滑
R_smooth = zeros(K, K);
for m = 1:L
    R_smooth = R_smooth + X(m:m+K-1, :) * X(m:m+K-1, :)';
end
R_smooth = R_smooth / L;

% 特征值分解
[Evec, Eval] = eig(R_smooth);
[~, idx] = sort(diag(Eval), 'descend');
Evec = Evec(:, idx);
En = Evec(:, K_signal+1:end); % 噪声子空间

% 使用MUSIC算法进行DOA估计（使用空间平滑）
P_MUSIC = zeros(size(theta_scan));
for i = 1:length(theta_scan)
    a = exp(-1j * 2 * pi * d * (0:K-1).' * sind(theta_scan(i)));
    P_MUSIC(i) = 1 / (a' * (En * En') * a);
end

% 绘制空间谱（使用空间平滑）
P_MUSIC = 10*log10(abs(P_MUSIC) / max(abs(P_MUSIC)));
figure;
plot(theta_scan, P_MUSIC, 'LineWidth', 2);
grid on;
xlabel('角度 (度)');
ylabel('空间谱 (dB)');
title('使用空间平滑的MUSIC DOA估计');

% 显示结果（使用空间平滑）
[~, locs] = findpeaks(P_MUSIC, 'SortStr', 'descend', 'NPeaks', K_signal);
estimated_angles = theta_scan(locs);
disp('使用空间平滑的估计DOA（度）：');
disp(estimated_angles);
