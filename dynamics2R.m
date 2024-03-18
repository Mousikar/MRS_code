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
kv = 8;
kp = 16;

for i=1:step
    % tau = [1;2];
    
    M = [(m1+m2)*l1^2 + m2*l2^2 + 2*m2*l1*l2*cos(theta2), m2*l2^2 + m2*l1*l2*cos(theta2);
        m2*l1*l2*cos(theta2)+m2*l2^2, m2*l2^2];
    H = [-m2*l1*l2*sin(theta2)*dottheta2^2 - 2*m2*l1*l2*sin(theta2)*dottheta1*dottheta2;
        m2*l1*l2*sin(theta2)*dottheta1^2];
    G = [m2*l2*g*cos(theta1+theta2) + (m1+m2)*l1*g*cos(theta1);
        m2*l2*g*cos(theta1+theta2)];

    tau = M * (ddotqd + kv * (dotqd - dotq) + kp * (qd(:,i) - q)) + H + G;
    
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


subplot(222)
plot(x_his,y_his)
hold on
plot(xd_his,yd_his,'k--')
legend
axis equal

subplot(223)
plot(dotq_his')
hold on
plot([dotqd(1)*ones(step,1),dotqd(2)*ones(step,1)],'k--')
legend

subplot(224)
plot(ddotq_his')
legend