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
% d
step = 4000;
qd = [pi/2*ones(1,step+1);linspace(-pi/6, -5*pi/6, step+1)];
dotqd = [0;-4*pi / (6*step * T)];
ddotqd = 0;

%% control
kv = 16;
kp = 64;

for i=1:step
    % tau = [0;0]; % 毫无控制
    theta1 = q(1);
    theta2 = q(2);
    dottheta1 = dotq(1);
    dottheta2 = dotq(2);
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
%%
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
%% 估计参数的控制
kv = 16;
kp = 64;

params_his = [];
params = 0.01*ones(7,1);%[m2*l2^2; m2*l1*l2;m1*l1^2;m2*l1^2;m2*l2;m1*l1;m2*l1];% 0.1*ones(7,1);
for i=1:step
    theta1 = q(1);
    theta2 = q(2);
    dottheta1 = dotq(1);
    dottheta2 = dotq(2);

    J = [-l1*sin(theta1)-l2*sin(theta1+theta2), -l2*sin(theta1+theta2);
          l1*cos(theta1)+l2*cos(theta1+theta2),  l2*cos(theta1+theta2)];

    Y = [sum(ddotq), cos(theta2)*(2*ddotq(1)+ddotq(2))-sin(theta2)*dotq(2)^2-2*sin(theta2)*dotq(1)*dotq(2), ddotq(1), ddotq(1), g * cos(theta1 + theta2), g* cos(theta1), g* cos(theta1);
        sum(ddotq), cos(theta2)*ddotq(1)+sin(theta2)*dotq(1)^2, 0,0, g * cos(theta1 + theta2), 0, 0];

    tau = Y * params - kv * (dotq - dotqd) - kp * (q - qd(:,i)); % 加上了控制

    dotparams = - 0.000001 * eye(7) * Y.' * (q - qd(:,i));

    params = params+dotparams *T;

    M = [(m1+m2)*l1^2 + m2*l2^2 + 2*m2*l1*l2*cos(theta2), m2*l2^2 + m2*l1*l2*cos(theta2);
        m2*l1*l2*cos(theta2)+m2*l2^2, m2*l2^2];
    H = [-m2*l1*l2*sin(theta2)*dottheta2^2 - 2*m2*l1*l2*sin(theta2)*dottheta1*dottheta2;
        m2*l1*l2*sin(theta2)*dottheta1^2];
    G = [m2*l2*g*cos(theta1+theta2) + (m1+m2)*l1*g*cos(theta1);
        m2*l2*g*cos(theta1+theta2)];

    ddotq = inv(M)*(tau - H - G);

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
    params_his = [params_his,params];
    disp(i)
end
figure
plot(params_his')
hold on

title('params')
%% F = m a control
clear

T = 0.001;
m = 5;
q_his = [];
dotq_his = [];
ddotq_his = [];

step = 4000;
t = linspace(0, T*step, step+1);
% qd = 3*t;%sin(t);
% dotqd = 3*ones(1,length(t));%cos(t);
% ddotqd = 0*ones(1,length(t));%-sin(t);

qd = sin(t);
dotqd = cos(t);
ddotqd = -sin(t);

q = 1;
dotq=0;
ddotq=0;

kp = 64;
kv = 2*sqrt(kp);
for i=1:step+1
    M = m;
    H = 0;
    G = 0;

    tau = M * (ddotqd(i) - kv * (dotq - dotqd(i)) - kp * (q - qd(i))) ...
        + H + G;

    ddotq = M\(tau - H - G);

    % ddotq = ddotqd(i) - kv * (dotq - dotqd(i)) - kp * (q - qd(i));

    dotq = dotq + ddotq * T;
    q = q + dotq * T;

    q_his = [q_his,q];
    dotq_his = [dotq_his,dotq];
    ddotq_his = [ddotq_his,ddotq];

end
%
figure
subplot(311)
plot(q_his)
hold on
plot(qd','k--')
title('q')

subplot(312)
plot(dotq_his)
hold on
plot(dotqd','k--')
title('q')

subplot(313)
plot(ddotq_his)
hold on
plot(ddotqd','k--')
title('q')
%% F = m a estimate 好的
clear

T = 0.001;
m = 5;
q_his = [];
dotq_his = [];
%% 
ddotq_his = [];
params_his = [];

step = 4000;
t = linspace(0, T*step, step+1);
% qd = 3*t;%sin(t);
% dotqd = 3*ones(1,length(t));%cos(t);
% ddotqd = 0*ones(1,length(t));%-sin(t);

qd = sin(t);
dotqd = cos(t);
ddotqd = -sin(t);

q = 1;
dotq=0;
ddotq=0;

params = 1;

kp = 100;
kv = 20;
for i=1:step+1
    M = m;
    H = 0;
    G = 0;
    Y = ddotqd(i);
    Y = ddotq;

    tau = Y * params - kv * (dotq - dotqd(i)) - kp * (q - qd(i)); % 加上了控制
    % tau = Y * params - kv * (dotq - dotqd(i)) - kp * (q - qd(i)); % 加上了控制

    dotparams =  - 0.01 * Y.' * (200 * (dotq - dotqd(i)) + 300 * (q - qd(i)));
    % dotparams = - 0.01 * Y.' * (q - qd(i));

    params = params+dotparams *T;

    ddotq = M\(tau - H - G);

    % ddotq = ddotqd(i) - kv * (dotq - dotqd(i)) - kp * (q - qd(i));

    dotq = dotq + ddotq * T;
    q = q + dotq * T;

    q_his = [q_his,q];
    dotq_his = [dotq_his,dotq];
    ddotq_his = [ddotq_his,ddotq];
    params_his = [params_his,params];

end
%
figure
subplot(411)
plot(q_his)
hold on
plot(qd','k--')
title('q')

subplot(412)
plot(dotq_his)
hold on
plot(dotqd','k--')
title('dotq')

subplot(413)
plot(ddotq_his)
hold on
plot(ddotqd','k--')
title('ddotq')

subplot(414)
plot(params_his)
hold on
plot(m* ones(step,1),'k--')
title('params')
%% F = m a estimate 坏的
clear

T = 0.001;
m = 5;
q_his = [];
dotq_his = [];
ddotq_his = [];
params_his = [];

step = 4000;
t = linspace(0, T*step, step+1);
% qd = 3*t;%sin(t);
% dotqd = 3*ones(1,length(t));%cos(t);
% ddotqd = 0*ones(1,length(t));%-sin(t);

qd = sin(t);
dotqd = cos(t);
ddotqd = -sin(t);

q = 1;
dotq=0;
ddotq=0;

params = 3;

c = 10;
k = 20;
for i=1:step+1
    M = m;
    H = 0;
    G = 0;
    Y = ddotqd(i) - c * (dotq - dotqd(i));

    tau = params * (ddotqd(i) - k *(dotq - dotqd(i)) - c * (q - qd(i)) ); % 加上了控制
    
    Y = ddotq;
    dotparams =  2 * ( - k* (dotq - dotqd(i)) - c * (q - qd(i)) );

    params = params + dotparams * T;

    ddotq = M\(tau - H - G);

    dotq = dotq + ddotq * T;
    q = q + dotq * T;

    q_his = [q_his,q];
    dotq_his = [dotq_his,dotq];
    ddotq_his = [ddotq_his,ddotq];
    params_his = [params_his,params];

end
%
figure
subplot(411)
plot(q_his)
hold on
plot(qd','k--')
title('q')

subplot(412)
plot(dotq_his)
hold on
plot(dotqd','k--')
title('dotq')

subplot(413)
plot(ddotq_his)
hold on
plot(ddotqd','k--')
title('ddotq')

subplot(414)
plot(params_his)
hold on
plot(m* ones(step,1),'k--')
title('params')
%% F = m a estimate 最后
clear

T = 0.001;
m = 5;
q_his = [];
dotq_his = [];
ddotq_his = [];
params_his = [];

step = 4000;
t = linspace(0, T*step, step+1);
% qd = 3*t;%sin(t);
% dotqd = 3*ones(1,length(t));%cos(t);
% ddotqd = 0*ones(1,length(t));%-sin(t);

qd = sin(t);
dotqd = cos(t);
ddotqd = -sin(t);

q = 1;
dotq=0;
ddotq=0;

params = 3;

c = 10;
k = 20;
for i=1:step+1
    M = m;
    H = 0;
    G = 0;

    tau = params * (ddotqd(i) - c *(dotq - dotqd(i)) )- k * ( (dotq - dotqd(i)) + c * (q - qd(i)) ); % 加上了控制
    dotparams = 0.02 * (ddotqd(i) - c *(dotq - dotqd(i)) ) * ( (dotq - dotqd(i)) + c * (q - qd(i)) );
    % tau = params * (ddotqd(i) - k *  (dotq - dotqd(i)) - c * (q - qd(i)) ); % 加上了控制
    % dotparams = - 0.02 * ((dotq - dotqd(i)) ) * (ddotqd(i) - k *  (dotq - dotqd(i)) - c * (q - qd(i)));

    params = params + dotparams * T;

    ddotq = M\(tau - H - G);

    dotq = dotq + ddotq * T;
    q = q + dotq * T;

    q_his = [q_his,q];
    dotq_his = [dotq_his,dotq];
    ddotq_his = [ddotq_his,ddotq];
    params_his = [params_his,params];

end
%
figure
subplot(411)
plot(q_his)
hold on
plot(qd','k--')
title('q')

subplot(412)
plot(dotq_his)
hold on
plot(dotqd','k--')
title('dotq')

subplot(413)
plot(ddotq_his)
hold on
plot(ddotqd','k--')
title('ddotq')

subplot(414)
plot(params_his)
hold on
plot(m* ones(step,1),'k--')
title('params')