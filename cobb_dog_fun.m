K = [0 : 0.1 : 10];
L = [0 : 0.1 : 10];
[x y] = meshgrid(K, L);
Q = 1.7 * K' .^ 0.6 * L .^ 0.4;
mesh(x, y, Q);
% title('������ƽ��');
xlabel('�Ͷ�');
ylabel('�ʱ�');