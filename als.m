function [q_est, R_est] = estimate_noise_covariance(Y, A, C, L, tau, N)
    % Y: 新息序列
    % A: 状态转移矩阵
    % C: 观测矩阵
    % L: Kalman增益
    % tau: 采样时间间隔
    % N: 自相关函数的阶数

    % 计算预测状态误差协方差矩阵 P 的稳态值
    % 这里假设 P 已知或通过其他方法估计得到
    % 为了简化，假设 P 是单位矩阵
    P = eye(size(A, 1));

    % 计算克罗内克积相关矩阵
    I = eye(size(A, 1));
    O = [C; C * A * L; C * A^2 * L^2]; % 这里简化为前3阶，实际应根据N计算
    Gamma = [I; -C * A * L; -C * A^2 * L^2]; % 这里简化为前3阶，实际应根据N计算

    % 构造 M 矩阵
    M = [tau, 0, 0, 0, 0, 0, 0, 0, 0;
         tau^3/3, tau^2/2, 0, tau^2/2, tau, 0, 0, 0, 0;
         tau^5/20, tau^4/8, tau^3/6, tau^4/8, tau^3/3, tau^2/2, tau^3/6, tau^2/2, tau];

    % 计算 A_Q 和 A_R
    G = eye(size(A, 1)); % 假设 G 为单位矩阵，实际应根据具体问题确定
    A_Q = (C * kron(O, I)) * ((eye(size(A, 1)^2) - kron(A, A)) \ ((G * kron(G, G)) * M));
    A_R = (C * kron(O, I)) * ((eye(size(A, 1)^2) - kron(A, A)) \ ((A * L * kron(A * L, A * L)))) + (kron(I, Gamma));

    % 计算 R(N)_s
    R_N = autocorr(Y, N); % 计算新息的自相关函数
    R_N_s = R_N(:); % 转换为列向量

    % 构造 A_LS 矩阵
    A_LS = [A_Q; A_R];

    % 最小二乘法估计 q 和 R
    theta_est = (A_LS' * A_LS) \ (A_LS' * R_N_s);
    q_est = theta_est(1:3); % 假设 q 有3个参数
    R_est = theta_est(4); % 假设 R 为标量
end

function r = autocorr(Y, N)
    % 计算新息的自相关函数
    r = zeros(N, 1);
    for n = 1:N
        r(n) = Y(1:end-n+1)' * Y(n:end) / length(Y(1:end-n+1));
    end
end