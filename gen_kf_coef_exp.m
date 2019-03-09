function [kf_coeff, kf_P] = gen_kf_coef_exp(coeff, P, tx, ty, Q, R)
% Predict    
    kf_coeff = coeff; % 3*1 coeff_{k|k-1}
    kf_P = P + Q; % 3*3 P_{k|k-1}
    
% Update
    H = [exp(-coeff(2) * tx), -coeff(1) * tx * exp(-coeff(2) * tx), 1]; % 1*3 Jacobian ()
    eps = ty - (coeff(1) * exp(-coeff(2) * tx) + coeff(3)); % 1*1 innovation
    S = H * kf_P * H' + R; % 1*1 innovation cov
    K = kf_P * H' * S ^ (-1); % 3*1 optimal kalman gain
    kf_coeff = kf_coeff + K * eps; % 3*1 updated state estimate
    kf_P = (eye(3) - K * H) * kf_P; % 3*3 updated cov estimate
end