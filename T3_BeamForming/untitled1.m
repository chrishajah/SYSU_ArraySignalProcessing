clc; clear; close all;

% 参数设置
angles = [10, 20]; % 目标来向角
SNR_dB = 20; % 信噪比
snapshots = 500; % 快拍数
d = 0.5; % 阵元间距（半波长）
wavelength = 1; % 假设波长为1
theta_scan = -90:0.1:90; % 扫描角范围
num_sources = length(angles); % 信号源个数
signal_power = 1; % 信号功率
noise_power = 10^(-SNR_dB / 10); % 噪声功率
S = sqrt(signal_power/2) * (randn(num_sources, snapshots) + 1j * randn(num_sources, snapshots)); % 随机信号
N = sqrt(noise_power/2) * (randn(num_elements, snapshots) + 1j * randn(num_elements, snapshots)); % 噪声
X = A * S + N; % 接收信号矩阵 (阵列接收信号)


