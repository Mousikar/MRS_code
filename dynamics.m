% 简单的移动机械臂动力学 
% 3自由度的机械臂+移动平台 
% 动力学方程是：
% M(q) \ddot q + C(q,\dot q) \dot q + G(q) + d(t) = B(q) \tau + f
% 其中, 
% M(q)——对称有界正定惯性矩阵
% C(q,\dot q) \dot q——向心和 科里奥利力矩
% G(q)——重力转矩矢量
% d(t)——额外扰动
% \tau——控制输入
% B(q)——已知的变换矩阵
% f——广义约束力
% q=[qv,qa]^T,qv是移动平台vehicle的关节变量，qa是机械臂arm的关节变量
clear
%% 初始值 & 参数
g = 9.8; % m/s^2
T = 0.001; % s

x = 2; %m
y = 0;
theta = pi/2;
theta1 = pi/2; % rad
theta2 = pi/2; % 不变，为常数
theta3 = 0.1; % m

dotx = 0%.5; %m/s
doty = 0;
dottheta = 0;
dottheta1 = 0; % rad
dottheta2 = 0; % 不变，为常数
dottheta3 = 0; % m

mp = 5; %kg
m1 = 1;
m2 = 0.5;
m3 = 0.5;

Ip = 2.5; %kg m^2
I1 = 1;
I2 = 0.5;
I3 = 0.5 + m3 * theta3^2; % 迭代的时候会改变

d = 0.5; %m
l = 0.5;
r = 0.5;
l1 = 1/2;
l2 = 0.3/2;
l3 = 0.3/2;

m123=m1+m2+m3;
mp123 = mp+m123;
%% 迭代
f = zeros(6,1);
qv=[x,y,theta]';
qa=[theta1,theta2,theta3]';
q = [qv;qa];

% q = [x,y,theta,theta1,theta2,theta3]';
dotq = [dotx,doty,dottheta,dottheta1,dottheta2,dottheta3]';
ddotq = zeros(6,1);

q_his = [];
dotq_his = [];
ddotq_his = [];
q_his = [q_his,q];
dotq_his = [dotq_his,dotq];
ddotq_his = [ddotq_his,ddotq];

taul = 1;
taur = 1;
tau1 = 1;
tau2 = 1;
tau3 = 1;

for i=1:4000
    tau = [taul,taur,tau1,tau2,tau3]';
    % tau = 1 * ones(5,1);
    x = q(1); %m
    y = q(2);
    theta = q(3);
    theta1 = q(4); % rad
    theta2 = q(5);%; %pi/2 不变，为常数
    theta3 = q(6); % m
    
    dotx = dotq(1); %m/s
    doty = dotq(2);
    dottheta = dotq(3);
    dottheta1 = dotq(4); % rad
    dottheta2 =dotq(5) ; %0 不变，为常数
    dottheta3 = dotq(6); % m

    
    Av = [cos(theta), sin(theta), 0];
    I3 = 0.5 + m3 * theta3^2; % 迭代的时候会改变
    
    L3=2*l2 +l3 +theta3;% 迭代的时候会改变
    m23 = m2*l2 +m3 *L3;
    I123=I1+I2+I3+m3*L3^2;
    I23=I2+I3+m3*L3^2;
    
    Mv12 = [ m123 * d * cos(theta + theta1),...
        m123 *d *sin(theta ) + m23 *sin(theta +theta1)]';
    Mv21 = Mv12';
    Mv11=diag([mp123,mp123]);
    Mv22=Ip + I123 + m123 *d^2 +m2 * (l2^2 +2*d*l2*cos(theta1))+ m3*(L3^2 +2 *d *L3 *cos(theta1));
    
    Cv11 =0*eye(2);
    Cv12 = [-m123*d*dottheta*sin(theta )- m23 *sin(theta + theta1)*(dottheta + dottheta1), ...
        m123 * d*dottheta*cos(theta )+m23*cos(theta + theta1)*(dottheta + dottheta1)]';
    Cv21=Cv12';
    Cv22 = -2*m23*d*sin(theta1)*dottheta1;
    
    Cva1 = [-m23*sin(theta + theta1)*(dottheta + dottheta1),...
        -m23*cos(theta + theta1)*(dottheta + dottheta1),0]'; % 这里参考文献没写清楚
    Cva2 = Cva1;
    Cva3 = [-m3*cos(theta + theta1)*(dottheta + dottheta1),...
        -m3*sin(theta + theta1)*(dottheta + dottheta1),0 ]'; % 这里参考文献没写清楚
    
    Cav1 = Cva1';
    Cav2 = Cva2';
    Cav3 = [m3*cos(theta + theta1)*(dottheta + dottheta1),...
        m3*sin(theta + theta1)*(dottheta + dottheta1),m3*d*sin(theta )*dottheta1];
    
    Mva1 = [m23 * cos(theta + theta1), m23 * sin(theta + theta1), ...
        I123 + m2*(l2^2+2*d*l2*cos(theta1))+m3*(L3^2+2*d*L3*cos(theta1))]';
    Mva2=0*zeros(3,1);
    Mva3=[sin(theta + theta1), -cos(theta + theta1),0]';
    
    Mv = [Mv11, Mv12; Mv21, Mv22];
    Cv = [Cv11, Cv12; Cv21, Cv22];
    Bv = [sin(theta / r ), - cos(theta/r ), -l/r;
         -sin(theta / r ), cos(theta/r ), l/r;]';
    Gv = [0, 0, 0]';
    
    Mva = [Mva1,Mva2,Mva3];
    Mav = Mva';
    Cva = [Cva1,Cva2,Cva3];
    Cav = [Cav1;Cav2;Cav3];
    
    Ba = diag([1,1,1]);
    Ma = diag([I123,I23,m3 ]);
    Ca = diag([-m23*d*sin(theta1)*dottheta, -m23*d*sin(theta1)*dottheta, 0]);
    Ga = [0,m2*g*l2,m3*g*L3]';
    
    M = [Mv, Mva;
         Mav, Ma];
    C = [Cv, Cva;
         Cav, Ca];
    G = [Gv;
        Ga];
    B = [Bv, zeros(3);
         zeros(3,2), Ba];

    ddotq = M \ (B*tau+f-G-C*dotq);
    dotq = dotq + ddotq * T;
    q = q + dotq * T;
    q_his = [q_his,q];
    dotq_his = [dotq_his,dotq];
    ddotq_his = [ddotq_his,ddotq];

end
%% 可视化
figure
subplot(221)
plot(q_his(1,:)',q_his(2,:)')
axis equal

subplot(222)
plot(q_his(3:6,:)')
legend

subplot(223)
plot(dotq_his(3:6,:)')
legend

subplot(224)
plot(ddotq_his(3:6,:)')
legend