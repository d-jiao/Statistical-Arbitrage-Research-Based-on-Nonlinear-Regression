% �趨����
model_spec = [0, 1, 1];
param_spec = [-0.001, 0.005, 0];
data_file = ['test_ptd.mat', 'sample_ptd.mat'];
[n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = back_test(model_spec, param_spec, data_file)

function [n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = back_test(model_spec, param_spec, data_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%
%               ���ڷ����Թ�ϵ��ͳ���������Իز⺯��
%                          VERSION 3.0
%                        ����, 2017��5��
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% �÷�:
% [n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, 
%            min_win, min_win_r, avg_r, win_r] = 
%            back_test(model_spec, param_spec)
%
% ˵��:
% �����������ġ����ڷ����Իع������ͳ�������о������ܵ�ģ��չ���ز�
% ���ڲ����ڸ��ֲ�����ģ��ϸ�ڵ������²��Եı���
% �������غ����ز�Ч���ı������������ۼ�������ͼ��������
% =========================================================================
% �������:
% model_spec                - 1 * 3���������趨�ز�ģ�͵Ľṹ
%                           - model_spec(1)��1��ʾָ������ģ��, 0��ʾ���ζ���ʽģ��
%                           - model_spec(2)��1��ʾ����������, 0��ʾ����������
%                           - model_spec(3)��1��ʾʹ�ÿ������˲�, 0��ʾ��ʹ��
% param_spec                - 1 * 3���������趨�ز�ģ�͵Ĳ���
%                           - param_spec(1)��ֹӯ��
%                           - param_spec(2)��ֹ����
%                           - param_spec(3)���볡��ֵ
% data_file                 - 1 * 2���������趨�ز����ݴ洢��·����.mat��ʽ��
%                           - data_file(1)�����������ݣ��洢��������Ϊtest_ptd
%                           - data_file(2)�����������ݣ��洢��������Ϊsample_ptd
%                           - ���������ļ�����һ��Ϊģ����������У��ڶ���Ϊģ�����������
% =========================================================================
% �������:
% n_trade                   - ���״���
% rt                        - ������
% annual_rt                 - �껯������
% m_d                       - ���س�
% s_r                       - ���ձ�
% annual_s_r                - �껯���ձ�
% max_win                   - ��󵥱�����
% max_win_r                 - ��󵥱�������
% min_win                   - ��󵥱ʿ���
% min_win_r                 - ��󵥱ʿ�����
% avg_r                     - ƽ��ÿ������
% win_r                     - ����ʤ��

    % ��ȡ�ز�����
    [tx ty] = get_test_data(model_spec(2), data_file);
    [x, y] = get_test_data(1, data_file);
    % ���ûز�ʱ�䴰�ڣ���ʼ����λ��ر���
    n = length(tx); % ���ݳ��� = �������� * ���ڶ��д���
    n_t = 1; % ���ڶ��д���
    n_d = n / n_t; % ��������
    [asset, money, cost, pnl, pos, code, n_trade, trade_record] = initialize_params();
    % ���ûز�Ĳ���
    cut_win = param_spec(1);
    cut_loss = param_spec(2);
    delta = param_spec(3);
    % ��ʼ���ز�ģ��
    if model_spec(3) == 1
        [P, Q, R, beta] = df_init_spec(model_spec(1), x, y);
    else
        [ff, df] = get_model(model(1), x, y);
    end
    % ���лز�
    for i = 1 : n_d
        % ��ʼ��ÿ�ղ�������¼���տ�������Ϳ��ֳɱ�
        temp_pos = [];
        temp_cost = [];
        for k = 1 : n_t
            % ��ʼ��ÿ�ζ��еĽ����źţ������ֽ���
            open_signal = 0;
            ind = (i - 1) * n_t + k;
            temp_money = money(end);
            % �Ը����˿������˲��Ĳ��ԣ�ʹ�ÿ������˲����²���
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
            % �����볡������⿪���ź�
            if ty(ind) <= eval(subs(ff, 'x', tx(ind))) - delta
                open_signal = 1;
                n_trade = n_trade + 1;
            end
            % ����յ������źţ��򽫿��ֵ������������λ���
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
            % ����й�ȥһ�켰����Ĳ�λ�������両��ӯ��
            % ���ڽ�����T+1�ģ���˶�����ͷ�粻���м��
            if n_pos ~= 0                              
                pnl = ((repmat([tx(ind), ty(ind)], n_pos, 1) - cost) .* pos * [1; 1]) ./ (cost .* pos * [-0.5; 1]);
                pnl = pnl';
            end
            stop_record = [];
            for j = 1 : n_pos
                % �Գ���ֹӯֹ���ߵĲ�λ����¼��ƽ����Ϣ������ƽ���ź�
                if pnl(j) < cut_loss |  pnl(j) > cut_win
                    close_cost = ((-pos(j, 1)) * cost(j, 1) * 0.5 + (-pos(j, 1)) * (cost(j, 1) - tx(ind))) + pos(j, 2) * ty(ind);
                    temp_money = temp_money + close_cost; % initial margin + pnl
                    stop_record = [stop_record, j];
                    c = code(j);
                    trade_record.close_time(c) = ind;
                    trade_record.end_value(c) = close_cost;
                end
            end
            % ����ƽ���źŸ��²�λ�ĳ�ʼ�ɱ����ֲ�����ӯ�������ױ�ŵļ�¼
            cost(stop_record, :) = [];
            pos(stop_record, :) = [];
            pnl(stop_record) = [];
            code(stop_record) = [];
            % ���³��е��ֽ�͵�ǰ�����ʲ�����ֵ
            money = [money, temp_money];
            temp_asset = get_asset(cost, temp_cost, temp_money, pos, temp_pos, tx(ind), ty(ind));                   
            asset = [asset, temp_asset];              
        end 
        % ��ÿ�ս��׽���������ʱ��λ�ͳɱ���¼�����ܵļ�¼
        pos = [pos; temp_pos];
        cost = [cost; temp_cost];
    end
    n_pos = numel(cost) / 2;
    % �Իز��ڽ���ʱ���еĲ�λ��ǿ��ƽ��
    if n_pos ~= 0
       for i = 1 : length(code)
           c = code(i);
           close_cost = ((-pos(i, 1)) * cost(i, 1) * 0.5 + (-pos(i, 1)) * (cost(i, 1) - tx(ind))) + pos(j, 2) * ty(ind);
           trade_record.close_time(c) = ind;
           trade_record.end_value(c) = close_cost;
       end
    end
    % ��������ز�Ч����ָ��
    [rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r] = output(asset, trade_record, n_d, n_trade);
    % ����Щָ���ӡ��������
    print_result(cut_win, cut_loss, n_trade, rt, annual_rt, m_d, s_r, annual_s_r, max_win, max_win_r, min_win, min_win_r, avg_r, win_r);
    % �����ۼ������ʵ�����ͼ    
    plot(asset/asset(1));
    title_name = sprintf('��ֵΪ%.3f', param_spec(3));
    title(title_name, 'FontSize', 20);
end 

function [tx ty] = get_test_data(sig, data_file);
    % ����ģ���趨��ȡ��Ӧ������
    load(data_file(1))
    load(data_file(2))
    if sig == 1        
        % ����������
        tx = sample_ptd(:,2);
        ty = sample_ptd(:,2);
    else
        % ����������
        tx = test_ptd(:,2);
        ty = test_ptd(:,1);
    end
end

function [P, Q, R, beta] = df_init_spec(model_id, x, y)
    % ����ģ���趨������Ӧ����
    if model_id == 1
        % ��ָ������ģ�ͼ�������ͳ����
        f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
        b = [2.5, -6, 2.5];
        [beta, ~, ~, ~, ~, ~] = nlinfit(x, y, f, b);
        beta = beta';
        P = zeros(3);
        [Q, R] = get_QR_exp(x, y);
    else
        % �Զ��ζ���ʽģ�ͼ�������ͳ����
        [beta, ~, ~, ~, ~] = regress(y, [ones(size(x)), x, x .^ 2]);
        P = zeros(3);
        [Q, R] = get_QR_2(x, y);
    end
end

function [ff, df] = get_model(sig, x, y)
    % �����趨��ȡ��Ӧģ��
    if sig == 1
        [ff, df] = initialize_exp_model(x, y) % ָ������ģ��
    else
        [ff, df] = initialize_quadratic_model(x, y) % ���ζ���ʽģ��
    end    
end

function [ff, df] = initialize_exp_model(x, y)
    % ��ʼ��ָ������ģ��
    f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
    b = [2.5, -6, 2.5];
    [beta, R, ~, ~, ~, ~] = nlinfit(x, y, f, b);
    syms x;
    ff = beta(1) * exp(-beta(2) * x) + beta(3);
    df = diff(ff);
end

function [ff, df] = initialize_quadratic_model(x, y)
    % ��ʼ�����ζ���ʽģ��
    [beta, ~, r, ~, ~] = regress(y, [ones(size(x)), x, x .^ 2]);
    syms x;
    ff = beta(1) + beta(2) .* x + beta(3) .* x .^ 2;
    df = diff(ff);
end

function [asset, money, cost, pnl, pos, code, n_trade trade_record] = initialize_params()
    % ��ʼ��������Ϣ�����׼�¼��
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
    % ������ζ���ʽģ�͵�����Э�������
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
    % ����ָ������ģ�͵�����Э�������
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
	% Ӧ�ø��¶��ζ���ʽģ�͵Ĳ�������   
    % Ԥ��
    kf_coeff = coeff; % 3*1 ϵ���ĳ���Ԥ��{k|k-1}
    kf_P = P + Q; % 3*3 ϵ��Э������ĳ���Ԥ��{k|k-1}    
	% ����
    H = [1, tx, tx ^ 2]; % 1*3 �ſɱȾ���
    eps = ty - (coeff(1) + coeff(2) * tx + coeff(3) * tx ^ 2); % 1*1 ��Ϣ
    S = H * kf_P * H' + R; % 1*1 ��Ϣ��Э������
    K = kf_P * H' * S ^ (-1); % 3*1 ���ſ���������
    kf_coeff = kf_coeff + K * eps; % 3*1 ���¶�ϵ���Ĺ���
    kf_P = (eye(3) - K * H) * kf_P; % 3*3 ���¶�ϵ��Э������Ĺ���
end

function [kf_coeff, kf_P] = gen_kf_coef_exp(coeff, P, tx, ty, Q, R)
    % Ӧ�ø��¶��ζ���ʽģ�͵Ĳ�������
	% Ԥ��   
    kf_coeff = coeff; % 3*1 ϵ���ĳ���Ԥ��{k|k-1}
    kf_P = P + Q; % 3*3 ϵ��Э������ĳ���Ԥ��{k|k-1}    
	% ����
    H = [exp(-coeff(2) * tx), -coeff(1) * tx * exp(-coeff(2) * tx), 1]; 
    eps = ty - (coeff(1) * exp(-coeff(2) * tx) + coeff(3)); % 1*1 ��Ϣ
    S = H * kf_P * H' + R; % 1*1 ��Ϣ��Э������
    K = kf_P * H' * S ^ (-1); % 3*1 ���ſ���������
    kf_coeff = kf_coeff + K * eps; % 3*1 ���¶�ϵ���Ĺ���
    kf_P = (eye(3) - K * H) * kf_P; % 3*3 ���¶�ϵ��Э������Ĺ���
end

function trade_record = update_record(trade_record, code, open_time, close_time, init_value, end_value)
    % ���½��׼�¼��
    trade_record.code = [trade_record.code, code];
    trade_record.open_time = [trade_record.open_time, open_time];
    trade_record.close_time = [trade_record.close_time, close_time];
    trade_record.init_value = [trade_record.init_value, init_value];
    trade_record.end_value = [trade_record.end_value, end_value];
end

function temp_asset = get_asset(cost, temp_cost, temp_money, pos, temp_pos, px, py)
    % ��������ʲ�����ֵ
    n_pos = numel(cost) / 2; % ǰһ������֮ǰ���ֵ�����
    temp_n_pos = numel(temp_cost) / 2; %��ǰ�����յĿ�������
    temp_asset = logical(0);
    if n_pos ~= 0
        % ����ǰһ������֮ǰ��λ����ֵ
        temp_asset = temp_money + sum((repmat([px, py], n_pos, 1) - cost) .* pos * [1; 1]) + sum(cost .* pos * [-0.5; 1]);
    end
    if temp_n_pos ~= 0
        % ���㵱ǰ�����տ��ֲ�λ����ֵ
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
    % ��������ز�Ч����ָ�꣬������ǰ����
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
    % �����ָ�������������
    fprintf('ֹӯ\t%.3f\n', cut_win);
    fprintf('ֹ��\t%.3f\n', cut_loss);
    fprintf('������\t%.4f\n', rt);
    fprintf('�껯������\t%.4f\n', annual_rt);
    fprintf('���ձ�\t%.4f\n', s_r);
    fprintf('�껯���ձ�\t%.4f\n', annual_s_r);
    fprintf('���س�\t%.4f\n', m_d);
    fprintf('���״���\t%.0f\n', n_trade);
    fprintf('��󵥱�����\t%.4f\n', max_win);
    fprintf('��󵥱�������\t%.4f\n', max_win_r);
    fprintf('��󵥱ʿ���\t%.4f\n', min_win);
    fprintf('��󵥱ʿ�����\t%.4f\n', min_win_r);
    fprintf('ƽ��ÿ������\t%.4f\n', avg_r);
    fprintf('����ʤ��\t%.4f\n', win_r);
end

function maxdrawdown = max_drawdown(dta)
    % �������س�
    maxdrawdown = 0;
    n = length(dta);
    for i = 1 : n - 2
        for j = i + 1 : n
            maxdrawdown = min(maxdrawdown, dta(j) - dta(i));
        end
    end
end

function sharperatio = sharpe_ratio(dta)
    % �������ձ�
    rtd = diff(log(dta));
    sharperatio = mean(rtd) / std(rtd);
end