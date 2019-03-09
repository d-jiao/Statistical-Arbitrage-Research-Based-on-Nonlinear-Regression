x = ptd_510050;
y = ptd_601668;
f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
b = [2.5, -6, 2.5];
[beta, R, ~, ~, ~, ~] = nlinfit(x, y, f, b);
beta = beta';
P = zeros(3);
Q = [604.8942, 11.1015, -614.3250; 11.1015, 0.3657, -11.3769; -614.3250, -11.3769, 624.1537];
% eye(3);%
R = 0.5771;


% 15min data 
% tx = test_ptq(:,2);
% ty = test_ptq(:,1);
% tx = test_ptq_2016(:,2);
% ty = test_ptq_2016(:,1);

% to day data
tx = test_ptd(:,2);
ty = test_ptd(:,1);
% tx = ptd_510050;
% ty = ptd_601668;

n = length(tx);
n_t = 1;
n_d = n / n_t;
  
[asset, money, cost, pnl, pos, code, n_trade, trade_record] = initialize_params();
delta = 0;
% delta = 1.4586 * 10e-2;

cut_loss = -0.01;
cut_win = 0.05;

for i = 1 : n_d
    
    temp_pos = [];
    temp_cost = [];
    for k = 1 : n_t
        open_signal = 0;
        ind = (i - 1) * n_t + k;
        temp_money = money(end);

        [beta, P] = gen_kf_coef(beta, P, tx(ind), ty(ind), Q, R);
        syms x;        
        ff = beta(1) * exp(-beta(2) * x) + beta(3);
        df = diff(ff);
        
%         ty(ind)
%         eval(subs(ff, 'x', tx(ind)))
        
        if ty(ind) <= eval(subs(ff, 'x', tx(ind))) - delta
            open_signal = 1;
            n_trade = n_trade + 1;
        end
        if open_signal == 1
            code = [code, n_trade];
            trade_record.code = [trade_record.code, n_trade];
            trade_record.open_time = [trade_record.open_time, ind];
            trade_record.close_time = [trade_record.close_time, 0];
            trade_record.end_value = [trade_record.end_value, 0];
            
            temp_pos = [temp_pos; -eval(subs(df, 'x', tx(ind))), 1];
            temp_cost = [temp_cost; tx(ind), ty(ind)];
            trade_cost = -temp_pos(end, 1) * temp_cost(end, 1) * 0.5 + temp_pos(end, 2) * temp_cost(end, 2);
            trade_record.init_value = [trade_record.init_value, trade_cost];
            temp_money = temp_money - trade_cost;
        end

        n_pos = numel(cost) / 2;
        temp_n_pos = numel(temp_cost) / 2;
                
        if n_pos ~= 0                              
            pnl = ((repmat([tx(ind), ty(ind)], n_pos, 1) - cost) .* pos * [1; 1]) ./ (cost .* pos * [-0.5; 1]);
            pnl = pnl';
        end

        stop_record = [];
        for j = 1 : n_pos
            if pnl(j) < cut_loss |  pnl(j) > cut_win
                close_cost = ((-pos(j, 1)) * cost(j, 1) * 0.5 + (-pos(j, 1)) * (cost(j, 1) - tx(ind))) + pos(j, 2) * ty(ind);
                temp_money = temp_money + close_cost;
                % initial margin + pnl
                stop_record = [stop_record, j];
                c = code(j);
                trade_record.close_time(c) = ind;
                trade_record.end_value(c) = close_cost;
            end
        end

        cost(stop_record, :) = [];
        pos(stop_record, :) = [];
        pnl(stop_record) = [];
        code(stop_record) = [];
        
        n_pos = numel(cost) / 2;

        money = [money, temp_money];
                
        temp_asset = logical(0);
        if n_pos ~= 0
            temp_asset = temp_money + sum((repmat([tx(ind), ty(ind)], n_pos, 1) - cost) .* pos * [1; 1]) + sum(cost .* pos * [-0.5; 1]);
        end

        if temp_n_pos ~= 0
            if ~temp_asset
                temp_asset = temp_money + sum((repmat([tx(ind), ty(ind)], temp_n_pos, 1) - temp_cost) .* temp_pos * [1; 1]) + sum(temp_cost .* temp_pos * [-0.5; 1]);
            else
                temp_asset = temp_asset + sum((repmat([tx(ind), ty(ind)], temp_n_pos, 1) - temp_cost) .* temp_pos * [1; 1]) + sum(temp_cost .* temp_pos * [-0.5; 1]);
            end
        end
        
        if ~temp_asset
            temp_asset = temp_money;
        end
            
        asset = [asset, temp_asset];
                
    end
    
    pos = [pos; temp_pos];
    cost = [cost; temp_cost];

end

n_pos = numel(cost) / 2;
if n_pos ~= 0
   for i = 1 : length(code)
       c = code(i);
       close_cost = ((-pos(i, 1)) * cost(i, 1) * 0.5 + (-pos(i, 1)) * (cost(i, 1) - tx(ind))) + pos(j, 2) * ty(ind);
       trade_record.close_time(c) = ind;
       trade_record.end_value(c) = close_cost;
   end
end

[rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = output(asset, trade_record, n_d, n_trade);
print_result(cut_win, cut_loss, n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r);
title('阈值为0', 'FontSize', 20)
% title('基于O-U的阈值', 'FontSize', 20)
function [asset, money, cost, pnl, pos, code, n_trade trade_record] = initialize_params()

    asset = [500];
    money = [500];
    cost = [];
    pnl = [];
    pos = [];
    code = [];
    n_trade = 0;
    trade_record = struct('code',[],'open_time',[],'close_time',[],'init_value',[],'end_value',[]);

end

function [kf_coeff, kf_P] = gen_kf_coef(coeff, P, tx, ty, Q, R)
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

function [rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = output(asset, trade_record, n_d, n_trade)

    rt = asset(end) / asset(1) - 1;
    annual_rt = rt * 252 / n_d;
    m_d = max_drawdown(asset);
    s_r = sharpe_ratio(asset);
    annual_s_r = s_r * sqrt(252);
    % plot(asset)
    plot(asset/500);
    max_win = max(trade_record.end_value(:) - trade_record.init_value(:));
    max_win_r = max((trade_record.end_value(:) - trade_record.init_value(:)) ./ trade_record.init_value(:));
    min_win = min(trade_record.end_value(:) - trade_record.init_value(:));
    min_win_r = min((trade_record.end_value(:) - trade_record.init_value(:)) ./ trade_record.init_value(:));
    avg_r = sum((trade_record.end_value(:) - trade_record.init_value(:)) ./ trade_record.init_value(:)) / n_trade;
    win_r = sum((trade_record.end_value(:) - trade_record.init_value(:)) > 0) / n_trade;
    
end

function print_result(cut_win, cut_loss, n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r)
    fprintf('止盈\t%.3f\n', cut_win);
    fprintf('止损\t%.3f\n', cut_loss);
    fprintf('收益率\t%.4f\n', rt);
    fprintf('年化收益率\t%.4f\n', annual_rt);
    fprintf('夏普比\t%.4f\n', s_r);
    fprintf('年化夏普比\t%.4f\n', annual_s_r);
    fprintf('最大回撤\t%.4f\n', m_d);
    fprintf('交易次数\t%.0f\n', n_trade);
    fprintf('最大单笔收益\t%.4f\n', max_win);
    fprintf('最大单笔收益率\t%.4f\n', max_win_r);
    fprintf('最大单笔亏损\t%.4f\n', min_win);
    fprintf('最大单笔亏损率\t%.4f\n', min_win_r);
    fprintf('平均每笔收益\t%.4f\n', avg_r);
    fprintf('交易胜率\t%.4f\n', win_r);
end

function maxdrawdown = max_drawdown(dta)

    maxdrawdown = 0;
    n = length(dta);
    for i = 1 : n - 2
        for j = i + 1 : n
            maxdrawdown = min(maxdrawdown, dta(j) - dta(i));
        end
    end
    
end

function sharperatio = sharpe_ratio(dta)

    rtd = diff(log(dta));
    sharperatio = mean(rtd) / std(rtd);
    
end