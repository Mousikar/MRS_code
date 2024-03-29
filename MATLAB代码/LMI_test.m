% LMI_test 时变非均匀时滞

clear;

n = 3;
tau_u = 0.5;%0.5 0.03

K1 = 2*eye(n);
K2 = 1/n*eye(n);
% follower不会和leader碰撞的网络拓扑
L = [1 -1 0;
      -1 1 0;
      0 -1 1]; % 跟随者邻接矩阵
B = diag(diag(L));
A = B-L;
d = sum(sum(A));

count=1;
A_cell = cell(d,1);

for i = 1:n
    for j = 1:n
        if A(i,j)==1
            A_celltemp = zeros(n);
            A_celltemp(i,j) = 1;
            A_cell{count} = A_celltemp;
            count = count+1;
        end
    end
end

%% 把拉普拉斯矩阵加进来的版本
setlmis([])
P = lmivar(1, [n 1]); % P是对称矩阵，阶数为n，满块矩阵
Q = lmivar(1, [n 1]);
R = lmivar(1, [n 1]);
S = lmivar(2, [n n]); % S是矩阵，阶数为n，满块矩阵

lmiterm([1 1 1 P], 1, -K1-K2*B, 's')
lmiterm([1 1 1 Q], 1, 1)
lmiterm([1 1 1 R], tau_u^2 * (-K1-K2*B)', -K1-K2*B)

lmiterm([1 1 2 S], 1, 1)

for i =1:d
    lmiterm([1 1 2+i P], 1, K2 * A_cell{i})
    lmiterm([1 1 2+i R], tau_u^2 * (-K1-K2*B)', K2 * A_cell{i})
    lmiterm([1 1 2+i S], -1/d, 1)
    lmiterm([1 1 2+i R], -1/d, 1)
end

lmiterm([1 2 2 Q], -1, 1)
lmiterm([1 2 2 R], -1, 1)

for i = 1:d
    lmiterm([1 2 2+i R], 1/d, 1)
    lmiterm([1 2 2+i -S], -1/d, 1)
end

for i = 1:d
    lmiterm([1 2+i 2+i R], -1/d, 1, 's')
    lmiterm([1 2+i 2+i S], 1/d, 1, 's')
end

for i = 1:d
    for j = i:d
        lmiterm([1 2+i 2+j R], tau_u^2*A_cell{i}.'*K2, K2*A_cell{j})
    end
end

lmiterm([-2,1,1,P],1,1)
lmiterm([2,1,1,0],0)

lmiterm([-3,1,1,Q],1,1)
lmiterm([3,1,1,0],0)

lmiterm([-5,1,1,R],1,1)
lmiterm([5,1,1,0],0)

lmiterm([-7,1,1,R],1,1)
lmiterm([-7,2,2,R],1,1)
lmiterm([-7,2,1,S],1,1)
lmiterm([7,1,1,0],0)
%%
lmis = getlmis;
[tmin, xfeas] = feasp(lmis);%,[0,0,0,0,0],-10^(-10)
%%
PP = dec2mat(lmis, xfeas, P)
QQ = dec2mat(lmis, xfeas, Q)
RR = dec2mat(lmis, xfeas, R)
SS = dec2mat(lmis, xfeas, S)
zuhe=[RR,SS.';SS,RR];
eig(zuhe)