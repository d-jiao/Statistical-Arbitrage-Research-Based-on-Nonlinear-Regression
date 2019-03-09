% 设定参数
model_spec = [0, 1, 1];
param_spec = [-0.001, 0.005, 0];
data_file = ['test_ptd.mat', 'sample_ptd.mat'];
[n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = back_test(model_spec, param_spec, data_file)

function [n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = back_test(model_spec, param_spec, data_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%
%               基于非线性关系的统计套利策略回测函数
%                          VERSION 3.0
%                        焦点, 2017年5月
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% 用法:
% [n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, 
%            min_win, min_win_r, avg_r, win_r] = 
%            back_test(model_spec, param_spec)
%
% 说明:
% 函数基于论文《基于非线性回归分析的统计套利研究》介绍的模型展开回测
% 用于测试在各种参数和模型细节的设置下策略的表现
% 函数返回衡量回测效果的变量，并绘制累计收益率图像，输出结果
% =========================================================================
% 输入参数:
% model_spec                - 1 * 3矩阵，用于设定回测模型的结构
%                           - model_spec(1)：1表示指数增长模型, 0表示二次多项式模型
%                           - model_spec(2)：1表示样本内数据, 0表示样本外数据
%                           - model_spec(3)：1表示使用卡尔曼滤波, 0表示不使用
% param_spec                - 1 * 3矩阵，用于设定回测模型的参数
%                           - param_spec(1)：止盈线
%                           - param_spec(2)：止损线
%                           - param_spec(3)：入场阈值
% data_file                 - 1 * 2矩阵，用于设定回测数据存储的路径（.mat格式）
%                           - data_file(1)：样本外数据，存储的数据名为test_ptd
%                           - data_file(2)：样本内数据，存储的数据名为sample_ptd
%                           - 对于两个文件，第一列为模型因变量序列，第二列为模型因变量序列
% =========================================================================
% 输出参数:
% n_trade                   - 交易次数
% rt                        - 收益率
% annual_rt                 - 年化收益率
% m_d                       - 最大回撤
% s_r                       - 夏普比
% annual_s_r                - 年化夏普比
% max_win                   - 最大单笔收益
% max_win_r                 - 最大单笔收益率
% min_win                   - 最大单笔亏损
% min_win_r                 - 最大单笔亏损率
% avg_r                     - 平均每笔收益
% win_r                     - 交易胜率

    % 获取回测数据
    [tx ty] = get_test_data(model_spec(2), data_file);
    [x, y] = get_test_data(1, data_file);
    % 设置回测时间窗口，初始化仓位监控变量
    n = length(tx); % 数据长度 = 交易天数 * 日内盯市次数
    n_t = 1; % 日内盯市次数
    n_d = n / n_t; % 交易天数
    [asset, money, cost, pnl, pos, code, n_trade, trade_record] = initialize_params();
    % 设置回测的参数
    cut_win = param_spec(1);
    cut_loss = param_spec(2);
    delta = param_spec(3);
    % 初始化回测模型
    if model_spec(3) == 1
        [P, Q, R, beta] = df_init_spec(model_spec(1), x, y);
    else
        [ff, df] = get_model(model(1), x, y);
    end
    % 进行回测
    for i = 1 : n_d
        % 初始化每日参数，记录当日开仓情况和开仓成本
        temp_pos = [];
        temp_cost = [];
        for k = 1 : n_t
            % 初始化每次盯市的交易信号，持有现金量
            open_signal = 0;
            ind = (i - 1) * n_t + k;
            temp_money = money(end);
            % 对附加了卡尔曼滤波的策略，使用卡尔曼滤波更新参数
            if model_spec(3) == 1
                if model_spec(1) == 1
                    [beta, P] = gen_kf_coef(beta, P, tx(ind), ty(ind), Q, R);
                    syms x;        
                    ff = beta(1) * exp(-beta(2) * x) + beta(3);
                    df = diff(ff);
                else
                    [beta, P] = gen_kf_coef(beta, P, tx(ind), ty(ind), Q, R);
                    syms x;        
                    ff = beta(1) + beta(2) .* x + beta(3) .* x .^ 2;
                    df = diff(ff);
                end
            end
            % 根据入场条件检测开仓信号
            if ty(ind) <= eval(subs(ff, 'x', tx(ind))) - delta
                open_signal = 1;
                n_trade = n_trade + 1;
            end
            % 如果收到开仓信号，则将开仓的相关情况计入仓位监控
            if open_signal == 1
                code = [code, n_trade];
                temp_pos = [temp_pos; -eval(subs(df, 'x', tx(ind))), 1];
                temp_cost = [temp_cost; tx(ind), ty(ind)];
                trade_cost = -temp_pos(end, 1) * temp_cost(end, 1) * 0.5 + temp_pos(end, 2) * temp_cost(end, 2);
                temp_money = temp_money - trade_cost;
                trade_record = update_record(trade_record, n_trade, ind, 0, trade_cost, 0);
            end
            n_pos = numel(cost) / 2;
            temp_n_pos = numel(temp_cost) / 2;
            % 如果有过去一天及更早的仓位，计算其浮动盈亏
            % 由于结算是T+1的，因此对日内头寸不进行监控
            if n_pos ~= 0                              
                pnl = ((repmat([tx(ind), ty(ind)], n_pos, 1) - cost) .* pos * [1; 1]) ./ (cost .* pos * [-0.5; 1]);
                pnl = pnl';
            end
            stop_record = [];
            for j = 1 : n_pos
                % 对超出止盈止损线的仓位，记录其平仓信息，生成平仓信号
                if pnl(j) < cut_loss |  pnl(j) > cut_win
                    close_cost = ((-pos(j, 1)) * cost(j, 1) * 0.5 + (-pos(j, 1)) * (cost(j, 1) - tx(ind))) + pos(j, 2) * ty(ind);
                    temp_money = temp_money + close_cost; % initial margin + pnl
                    stop_record = [stop_record, j];
                    c = code(j);
                    trade_record.close_time(c) = ind;
                    trade_record.end_value(c) = close_cost;
                end
            end
            % 根据平仓信号更新仓位的初始成本、持仓量、盈亏、交易编号的记录
            cost(stop_record, :) = [];
            pos(stop_record, :) = [];
            pnl(stop_record) = [];
            code(stop_record) = [];
            % 更新持有的现金和当前持有资产的市值
            money = [money, temp_money];
            temp_asset = get_asset(cost, temp_cost, temp_money, pos, temp_pos, tx(ind), ty(ind));                   
            asset = [asset, temp_asset];              
        end 
        % 在每日交易结束，把临时仓位和成本记录计入总的记录
        pos = [pos; temp_pos];
        cost = [cost; temp_cost];
    end
    n_pos = numel(cost) / 2;
    % 对回测期结束时持有的仓位，强制平仓
    if n_pos ~= 0
       for i = 1 : length(code)
           c = code(i);
           close_cost = ((-pos(i, 1)) * cost(i, 1) * 0.5 + (-pos(i, 1)) * (cost(i, 1) - tx(ind))) + pos(j, 2) * ty(ind);
           trade_record.close_time(c) = ind;
           trade_record.end_value(c) = close_cost;
       end
    end
    % 计算衡量回测效果的指标
    [rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = output(asset, trade_record, n_d, n_trade);
    % 将这些指标打印到窗口中
    print_result(cut_win, cut_loss, n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r);
    % 绘制累计收益率的折线图    
    plot(asset/asset(1));
    title_name = sprintf('阈值为%.3f', param_spec(3));
    title(title_name, 'FontSize', 20);
end 

function [tx ty] = get_test_data(sig, data_file);
    % 根据模型设定读取相应的数据
    load(data_file(1))
    load(data_file(2))
    if sig == 1        
        % 样本内数据
        tx = sample_ptd(:,2);
        ty = sample_ptd(:,2);
    else
        % 样本外数据
        tx = test_ptd(:,2);
        ty = test_ptd(:,1);
    end
end

function [P, Q, R, beta] = df_init_spec(model_id, x, y)
    % 根据模型设定计算相应数据
    if model_id == 1
        % 对指数增长模型计算先验统计量
        f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
        b = [2.5, -6, 2.5];
        [beta, ~, ~, ~, ~, ~] = nlinfit(x, y, f, b);
        beta = beta';
        P = zeros(3);
        [Q, R] = get_QR_exp(x, y);
    else
        % 对二次多项式模型计算先验统计量
        [beta, ~, ~, ~, ~] = regress(y, [ones(size(x)), x, x .^ 2]);
        P = zeros(3);
        [Q, R] = get_QR_2(x, y);
    end
end

function [ff, df] = get_model(sig, x, y)
    % 根据设定获取相应模型
    if sig == 1
        [ff, df] = initialize_exp_model(x, y) % 指数增长模型
    else
        [ff, df] = initialize_quadratic_model(x, y) % 二次多项式模型
    end    
end

function [ff, df] = initialize_exp_model(x, y)
    % 初始化指数增长模型
    f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
    b = [2.5, -6, 2.5];
    [beta, R, ~, ~, ~, ~] = nlinfit(x, y, f, b);
    syms x;
    ff = beta(1) * exp(-beta(2) * x) + beta(3);
    df = diff(ff);
end

function [ff, df] = initialize_quadratic_model(x, y)
    % 初始化二次多项式模型
    [beta, ~, r, ~, ~] = regress(y, [ones(size(x)), x, x .^ 2]);
    syms x;
    ff = beta(1) + beta(2) .* x + beta(3) .* x .^ 2;
    df = diff(ff);
end

function [asset, money, cost, pnl, pos, code, n_trade trade_record] = initialize_params()
    % 初始化交易信息及交易记录簿
    asset = [100];
    money = [100];
    cost = [];
    pnl = [];
    pos = [];
    code = [];
    n_trade = 0;
    trade_record = struct('code',[],'open_time',[],'close_time',[],'init_value',[],'end_value',[]);
end

function [Q, R] = get_QR_2(x, y)
    % 计算二次多项式模型的样本协方差矩阵
    [~, ~, res, ~, ~] = regress(y, [ones(size(x)), x, x .^ 2]);
    R = std(res);
    n = length(x) / 2;
    coeff = [];
    for i = 1 : n
        [beta, ~, ~, ~, ~] = regress(y(i + 1 : n + i), [ones(n, 1), x(i + 1 : n + i), x(i + 1 : n + i) .^ 2]);
        coeff = [coeff, beta];
    end
    coeff = coeff';
    Q = cov(diff(coeff));
end

function [Q, R] = get_QR_exp(x, y)
    % 计算指数增长模型的样本协方差矩阵
    f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
    b = [2.5, -6, 2.5];
    [beta, res, ~, ~, ~, ~] = nlinfit(x, y, f, b);
    R = std(res);
    n = length(x) / 2;
    coeff = [];
    for i = 1 : n
        f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
        for j = 1 : 3
            [beta, ~, ~, ~, ~, ~] = nlinfit(x(i : i + n), y(i : i + n), f, beta);
            beta = lsqcurvefit(f, beta, x(i : i + n), y(i : i + n));
        end
        coeff = [coeff; beta];
    end
    Q = cov(diff(coeff));
end

function [kf_coeff, kf_P] = gen_kf_coef_2(coeff, P, tx, ty, Q, R)
	% 应用更新二次多项式模型的参数估计   
    % 预测
    kf_coeff = coeff; % 3*1 系数的初步预测{k|k-1}
    kf_P = P + Q; % 3*3 系数协方差阵的初步预测{k|k-1}    
	% 更新
    H = [1, tx, tx ^ 2]; % 1*3 雅可比矩阵
    eps = ty - (coeff(1) + coeff(2) * tx + coeff(3) * tx ^ 2); % 1*1 新息
    S = H * kf_P * H' + R; % 1*1 新息的协方差阵
    K = kf_P * H' * S ^ (-1); % 3*1 最优卡尔曼增益
    kf_coeff = kf_coeff + K * eps; % 3*1 更新对系数的估计
    kf_P = (eye(3) - K * H) * kf_P; % 3*3 更新对系数协方差阵的估计
end

function [kf_coeff, kf_P] = gen_kf_coef_exp(coeff, P, tx, ty, Q, R)
    % 应用更新二次多项式模型的参数估计
	% 预测   
    kf_coeff = coeff; % 3*1 系数的初步预测{k|k-1}
    kf_P = P + Q; % 3*3 系数协方差阵的初步预测{k|k-1}    
	% 更新
    H = [exp(-coeff(2) * tx), -coeff(1) * tx * exp(-coeff(2) * tx), 1]; 
    eps = ty - (coeff(1) * exp(-coeff(2) * tx) + coeff(3)); % 1*1 新息
    S = H * kf_P * H' + R; % 1*1 新息的协方差阵
    K = kf_P * H' * S ^ (-1); % 3*1 最优卡尔曼增益
    kf_coeff = kf_coeff + K * eps; % 3*1 更新对系数的估计
    kf_P = (eye(3) - K * H) * kf_P; % 3*3 更新对系数协方差阵的估计
end

function trade_record = update_record(trade_record, code, open_time, close_time, init_value, end_value)
    % 更新交易记录簿
    trade_record.code = [trade_record.code, code];
    trade_record.open_time = [trade_record.open_time, open_time];
    trade_record.close_time = [trade_record.close_time, close_time];
    trade_record.init_value = [trade_record.init_value, init_value];
    trade_record.end_value = [trade_record.end_value, end_value];
end

function temp_asset = get_asset(cost, temp_cost, temp_money, pos, temp_pos, px, py)
    % 计算持有资产的市值
    n_pos = numel(cost) / 2; % 前一交易日之前开仓的数量
    temp_n_pos = numel(temp_cost) / 2; %当前交易日的开仓数量
    temp_asset = logical(0);
    if n_pos ~= 0
        % 计算前一交易日之前仓位的市值
        temp_asset = temp_money + sum((repmat([px, py], n_pos, 1) - cost) .* pos * [1; 1]) + sum(cost .* pos * [-0.5; 1]);
    end
    if temp_n_pos ~= 0
        % 计算当前交易日开仓仓位的市值
        if ~temp_asset
            temp_asset = temp_money + sum((repmat([px, py], temp_n_pos, 1) - temp_cost) .* temp_pos * [1; 1]) + sum(temp_cost .* temp_pos * [-0.5; 1]);
        else
            temp_asset = temp_asset + sum((repmat([px, py], temp_n_pos, 1) - temp_cost) .* temp_pos * [1; 1]) + sum(temp_cost .* temp_pos * [-0.5; 1]);
        end
    end      
    if ~temp_asset
        temp_asset = temp_money;
    end
end

function [rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = output(asset, trade_record, n_d, n_trade)
    % 计算衡量回测效果的指标，含义如前所述
    rt = asset(end) / asset(1) - 1;
    annual_rt = rt * 244 / n_d;
    m_d = max_drawdown(asset);
    s_r = sharpe_ratio(asset);
    annual_s_r = s_r * sqrt(244);
    max_win = max(trade_record.end_value(:) - trade_record.init_value(:));
    max_win_r = max((trade_record.end_value(:) - trade_record.init_value(:)) ./ trade_record.init_value(:));
    min_win = min(trade_record.end_value(:) - trade_record.init_value(:));
    min_win_r = min((trade_record.end_value(:) - trade_record.init_value(:)) ./ trade_record.init_value(:));
    avg_r = sum((trade_record.end_value(:) - trade_record.init_value(:)) ./ trade_record.init_value(:)) / n_trade;
    win_r = sum((trade_record.end_value(:) - trade_record.init_value(:)) > 0) / n_trade; 
end

function print_result(cut_win, cut_loss, n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r)
    % 将相关指标输出到窗口中
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
    % 计算最大回撤
    maxdrawdown = 0;
    n = length(dta);
    for i = 1 : n - 2
        for j = i + 1 : n
            maxdrawdown = min(maxdrawdown, dta(j) - dta(i));
        end
    end
end

function sharperatio = sharpe_ratio(dta)
    % 计算夏普比
    rtd = diff(log(dta));
    sharperatio = mean(rtd) / std(rtd);
end