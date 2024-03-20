% dynamics2R
% M ddot q + H + G = tau 
clear

g = 9.8;
T = 0.001;
m1 = 0.7;
m2 = 0.5;
l1 = 0.4;
l2 = 0.28;

theta1 = 0;
theta2 = 0;

dottheta1 = 0;
dottheta2 = 0;

q = [theta1,theta2]';
dotq = [dottheta1,dottheta2]';
ddotq = zeros(2,1);
x = l1*cos(theta1) + l2*cos(theta1+theta2);
y = l1*sin(theta1) + l2*sin(theta1+theta2);

q_his = [];
dotq_his = [];
ddotq_his = [];
x_his = [];
y_his = [];
q_his = [q_his,q];
dotq_his = [dotq_his,dotq];
ddotq_his = [ddotq_his,ddotq];
x_his = [x_his,x];
y_his = [y_his,y];

xd_his = [];
yd_his = [];
%% d
step = 4000;
qd = [pi/2*ones(1,step+1);linspace(-pi/2, pi/2, step+1)];
dotqd = [0;pi / (step * T)];
ddotqd = 0;

%% control
kv = 16;
kp = 64;

for i=1:step
    % tau = [1;2]; % 毫无控制
    
    M = [(m1+m2)*l1^2 + m2*l2^2 + 2*m2*l1*l2*cos(theta2), m2*l2^2 + m2*l1*l2*cos(theta2);
        m2*l1*l2*cos(theta2)+m2*l2^2, m2*l2^2];
    H = [-m2*l1*l2*sin(theta2)*dottheta2^2 - 2*m2*l1*l2*sin(theta2)*dottheta1*dottheta2;
        m2*l1*l2*sin(theta2)*dottheta1^2];
    G = [m2*l2*g*cos(theta1+theta2) + (m1+m2)*l1*g*cos(theta1);
        m2*l2*g*cos(theta1+theta2)];

    tau = M * (ddotqd + kv * (dotqd - dotq) + kp * (qd(:,i) - q)) + H + G; % 加上了控制
    
    ddotq = M\(tau - H - G);

    dotq = dotq + ddotq * T;
    q = q + dotq * T;

    x = l1*cos(q(1)) + l2*cos(q(1)+q(2));
    y = l1*sin(q(1)) + l2*sin(q(1)+q(2));
    q_his = [q_his,q];
    dotq_his = [dotq_his,dotq];
    ddotq_his = [ddotq_his,ddotq];
    x_his = [x_his,x];
    y_his = [y_his,y];

    xd = l1*cos(qd(1,i)) + l2*cos(qd(1,i)+qd(2,i));
    yd = l1*sin(qd(1,i)) + l2*sin(qd(1,i)+qd(2,i));
    xd_his = [xd_his,xd];
    yd_his = [yd_his,yd];
end

figure
subplot(221)
plot(q_his')
hold on
plot(qd','k--')
legend('实际角度1','实际角度2','期望角度1','期望角度2')
title('关节角度')

subplot(222)
plot(x_his,y_his)
hold on
plot(xd_his,yd_his,'k--')
legend('实际轨迹','期望轨迹')
axis equal
title('末端轨迹')

subplot(223)
plot(dotq_his')
hold on
plot([dotqd(1)*ones(step,1),dotqd(2)*ones(step,1)],'k--')
legend('实际角速度1','实际角速度2','期望角速度1','期望角速度2')
title('关节角速度')

subplot(224)
plot(ddotq_his')
legend('实际角加速度1','实际角加速度2')
title('关节角加速度')