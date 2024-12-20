%% 三个抓取点
% 定义抓取点，每列为一个抓取点
r = [0.4    -0.4    -0.4; 
    0       0.3     -0.3; 
    0       0       0];
n = size(r,2);
% 物体目标受力
Fo = [100 0 0 0 0 0].';
% 抓取矩阵
G=zeros(6,6*n);
for i = 1:n
    G(:,(i-1)*6+1:i*6) = [
        eye(3), zeros(3);
        skew(r(:,i)), eye(3)];
end
G
% Moore–Penrose inverse
G_Moore = pinv(G);
% plus_delta
G_delta = zeros(6*n,6);
for i = 1:n
    G_delta((i-1)*6+1:i*6,:) = 1/n*[
        eye(3), zeros(3);
        skew(r(:,i)).', eye(3)];
end
% free
m_star = rand(3,1);% ; % 随机分配就有内力 [1 1 1]均匀分配没有内力；因为要满足r m 求和等于0
m0 = sum(m_star);
J_star = eye(3,3);
J0 = n * J_star;
G_free = zeros(6*n,6);
for i = 1:n
    G_free((i-1)*6+1:i*6,:) = [
        m_star(i)/m0*eye(3), m_star(i)*inv(J0)*(skew(r(:,i)).');
        zeros(3), J_star/J0];
end

% G_plus = G_Moore
% G_plus = G_delta
G_plus = G_free
h = G_plus * Fo % 得到的力是平均分配的
h_ext = G_plus * G * h
h_int = (eye(6*n)-G_plus * G) * h

%% 两个抓取点
% 定义抓取点，每列为一个抓取点
r = [0.4    -0.4    ; 
    0       0     ; 
    0       0       ];
n = size(r,2);
% 物体目标受力
Fo = [10 0 0 0 0 0].';
% 抓取矩阵
G=zeros(6,6*n);
for i = 1:n
    G(:,(i-1)*6+1:i*6) = [
        eye(3), zeros(3);
        skew(r(:,i)), eye(3)];
end
G
% Moore–Penrose inverse
G_Moore = pinv(G);
% plus_delta
G_delta = zeros(6*n,6);
for i = 1:n
    G_delta((i-1)*6+1:i*6,:) = 1/n*[
        eye(3), zeros(3);
        skew(r(:,i)).', eye(3)];
end
% free 已知作用点位置，定系数，或者已知作用力大小定作用点位置都可以实现内力是0
m_star = [5 5];% rand(3,1); % 随机分配就有内力 均匀分配没有内力
m0 = sum(m_star);
J_star = eye(3,3);
J0 = n * J_star;
G_free = zeros(6*n,6);
for i = 1:n
    G_free((i-1)*6+1:i*6,:) = [
        m_star(i)/m0*eye(3), m_star(i)*inv(J0)*(skew(r(:,i)).');
        zeros(3), J_star/J0];
end

% G_plus = G_Moore
% G_plus = G_delta
G_plus = G_free
h = G_plus * Fo % 得到的力是平均分配的
h_ext = G_plus * G * h
h_int = (eye(6*n)-G_plus * G) * h


