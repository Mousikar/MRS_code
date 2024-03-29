%% three mobile manipulators dynamics
clear
%% 系统参数
T = 0.001; % s
g = 9.8; % m/s^2

m0 = 15; % kg
m1 = 1.5;
m2 = 1.5;

l1 = 0.3; % m
l2 = 0.3;
lf = 0.3; % 文中没给数据

rhoInf = 0.1;
rho0 = 0.5;
Tij = 3;

Kiu = 1;
Kil = 1;

kp = 3 * eye(15);
kd = 3 * eye(15);
sigma = 0.8 * eye(15);

mci = 3 * eye(2);
bci = 1 * eye(2);
kci = 3 * eye(2);

ri = 0.5 * eye(2);

iota1 = 0.8;
iota2 = 0.7;
omega1 = 40/180*pi;
omega2 = 20/180*pi;

m11 = m0+m1+m2;
m22 = m11;
%% 初值 和 初始化
% qi = [xp, yp, thetap, theta1, theta2];
xp = zeros(1,3);
yp = zeros(1,3);
thetap = zeros(1,3);
theta1 = zeros(1,3);
theta2 = zeros(1,3);

dotxp = zeros(1,3);
dotyp = zeros(1,3);
dotthetap = zeros(1,3);
dottheta1 = zeros(1,3);
dottheta2 = zeros(1,3);

ddotxp = zeros(1,3);
ddotyp = zeros(1,3);
ddotthetap = zeros(1,3);
ddottheta1 = zeros(1,3);
ddottheta2 = zeros(1,3);

Mm(:,:,3) = zeros(3);
Deltam(:,:,3) = zeros(3,1);
%% 历史变量
x1d_his = [];
y1d_his = [];
x2d_his = [];
y2d_his = [];
x3d_his = [];
y3d_his = [];
%% 迭代
for k = 1:15000
    t = k * T;
    I1 = 0.05 + 0.01 * sin(t);
    I2 = 0.05 + 0.01 * sin(t);
    I0 = 0.5 + 0.01 * sin(t); % 文中有没给数据
    if t <= 3
        x1d = 0;
        y1d = 0;
        beta = 0;
        x2d = x1d - iota1 * cos(beta + omega1);
        y2d = y1d - iota1 * sin(beta + omega1);
        x3d = x1d - iota2 * cos(omega2 - beta);
        y3d = y1d + iota2 * sin(omega2 - beta);
    elseif t > 3 && t<= 3+3.75*pi % 14.7810
        x1d = 0.8 * (t - 3);
        y1d = 2 * cos(x1d) - 2;
        beta = atan2(-2*sin(x1d),1);
        x2d = x1d - iota1 * cos(beta + omega1);
        y2d = y1d - iota1 * sin(beta + omega1);
        x3d = x1d - iota2 * cos(omega2 - beta);
        y3d = y1d + iota2 * sin(omega2 - beta);
    else
        x1d = 3 * pi;
        y1d = -4;
        beta = 0;
        x2d = x1d - iota1 * cos(beta + omega1);
        y2d = y1d - iota1 * sin(beta + omega1);
        x3d = x1d - iota2 * cos(omega2 - beta);
        y3d = y1d + iota2 * sin(omega2 - beta);
    end
    x1d_his = [x1d_his;x1d];
    y1d_his = [y1d_his;y1d];
    x2d_his = [x2d_his;x2d];
    y2d_his = [y2d_his;y2d];
    x3d_his = [x3d_his;x3d];
    y3d_his = [y3d_his;y3d];
    if t>=5 && t<=8
        fe = 5 *eye(15);
    else
        fe = 0 *eye(15);
    end
    for i = 1:3
        m33 = m2 * ( l1 * cos(theta1(i)) + l2 * cos(theta1(i)+theta2(i)) )^2 + ...
            I2 * cos(theta1(i)+theta2(i))^2 + (m0 * lf^2 + I0) + I1 * cos(theta1(i))^2;
        m13 = m0*lf*sin(thetap(i)) - m1*l1*cos(theta1(i))*sin(thetap(i)) - ...
            m2*( l1*cos(theta1(i))+l2*cos(theta1(i)+theta2(i)) )*sin(thetap(i));
        m31 = m13;
        m23 = -m0*lf*cos(thetap(i)) + m1*l1*cos(theta1(i))*cos(thetap(i)) +...
            m2*( l1*cos(theta1(i))+l2*cos(theta1(i)+theta2(i)) )*cos(thetap(i));
        m32 = m23;
        Mm(:,:,i) = [m11, 0,   m13;
              0,   m22, m23;
              m31, m32, m33];

        Deltamx = ( -m1*l1*sin(theta1(i))*cos(thetap(i)) - m2*( l1*sin(theta1(i))+l2*sin(theta1(i)...
            +theta2(i)) )*cos(thetap(i)) )*ddottheta2(i) + (-m2*l2*sin(theta1(i)+theta2(i))*cos(thetap(i)))*ddottheta2(i)...
            -( m1*l1*cos(theta1(i))*cos(thetap(i))*dottheta1(i) +m1*l1*sin(theta1(i))*sin(thetap(i))*dotthetap...
            + m2*( l1*sin(theta1(i)) + l2*sin(theta1(i)+theta2(i)) )*sin(thetap(i))*dotthetap(i)...
            - m2*( l1*cos(theta1(i)) + l2*cos(theta1(i)+theta2(i)) )*cos(thetap(i))*dottheta1(i)...
            -m2*l2*cos(theta1(i)+theta2(i))*cos(thetap(i))*dottheta2(i) )*dottheta1(i)...
            +(-m2*l2*cos(theta1(i)+theta2(i))*cos(thetap(i))*dottheta2(i) - m2*l2*cos(theta1(i)...
            +theta2(i))*cos(thetap(i))*dottheta1(i) +m2*l2*sin(theta1(i)+theta2(i))*sin(thetap(i))*dotthetap(i) )*dottheta2(i);
        Deltamy = 
        % Deltam(:,:,i) = [Deltamx, Deltamy, Deltamtheta]';
    end
end
%% 可视化
figure
plot(x1d_his,y1d_his,'r');axis equal;hold on;
plot(x2d_his,y2d_his,'g');
plot(x3d_his,y3d_his,'b');
for k = 1500:1500:15000
    plot([x1d_his(k),x2d_his(k),x3d_his(k),x1d_his(k)],...
        [y1d_his(k),y2d_his(k),y3d_his(k),y1d_his(k)])
end