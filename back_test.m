f = @(b, x) b(1) * exp(-b(2) * x) + b(3);
x = ptd_510050;
y = ptd_601668;
b = [2.5, -6, 2.5];
[beta, R, ~, ~, ~, ~] = nlinfit(x, y, f, b);

syms x;
ff = beta(1) * exp(-beta(2) * x) + beta(3);
% f1 = @(x) f(beta, x)
df = diff(ff);

x = ptd_510050;
tx = x;%test_pthh(:,2);
ty = y;%test_pthh(:,1);
n = length(tx);
  
asset = 100;
money = [100];
cost = [];
pnl = [];
pos = [];
nums = 0;

for i = 1 : n
    if ty(i) <= eval(subs(ff, 'x', tx(i))) - std(R)
        open_signal = 1;
        nums = nums + 1;
    else
        open_signal = 0;
    end
    
    temp_money = money(end);
    
    if open_signal == 1
        temp_pos = -eval(subs(df, 'x', tx(i)));
        pos = [pos; [temp_pos, 1]];
        cost = [cost; [tx(i), ty(i)]];
        temp_money = money(end) - (temp_pos * tx(i) + ty(i));
    end
    
    n_pos = numel(cost) / 2;
    if n_pos ~= 0                              
%         pnl = sum((cost ./ repmat([tx(i), ty(i)], n_pos, 1) .* [ones(n_pos, 1), -ones(n_pos, 1)] - [ones(n_pos, 1), -ones(n_pos, 1)])')
        pnl = (sum((repmat([tx(i), ty(i)], n_pos, 1) .* pos)') - sum((cost .* pos)')) ./ abs(sum((cost .* pos .* [-0.5 .* ones(n_pos, 1), ones(n_pos, 1)])'));
    end
    
    stop_record = [];
    for j = 1 : n_pos
        if pnl(j) < -0.05 |  pnl(j) > 0.2
            temp_money = temp_money + pos(j, 1) * tx(i) + pos(j, 2) * ty(i);
            stop_record = [stop_record, j];
%         else
%             stop_record = [stop_record, 0];
        end
    end
    
    cost(stop_record, :) = [];
    pos(stop_record, :) = [];
    pnl(stop_record) = [];
    n_pos = numel(cost) / 2;
        
    money = [money, temp_money];
    if n_pos ~= 0
        asset = [asset, temp_money + sum(sum(pos .* repmat([tx(i), ty(i)], n_pos, 1)))];
    else
        asset = [asset, asset(end)];
    end
end

m_d = max_drawdown(asset)
s_r = sharpe_ratio(asset)
plot(trade_days_2016, asset(2 : end))
datetick('x', 'yyyy/mm/dd');
