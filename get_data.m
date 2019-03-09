load('ptd.mat')
load('universe.mat')
load('test_ptq.mat')
load('test_ptd.mat')
load('test_ptq_2016.mat')
load('trade_days_2016.mat')

[~, ind1] = ismember('510050.SH',universe);
[~, ind2] = ismember('510900.SH',universe);
[~, ind3] = ismember('601668.SH',universe);
[~, ind4] = ismember('000839.SZ',universe);
[~, ind5] = ismember('600028.SH',universe);

ptd_510050 = ptd(:,ind1);
ptd_510900 = ptd(:,ind2);
ptd_601668 = ptd(:,ind3);
ptd_000839 = ptd(:,ind4);
ptd_600028 = ptd(:,ind5);