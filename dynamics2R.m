% dynamics2R
% M ddot q + H + G = tau 
clear

g = 9.8;
m1 = 0.7;
m2 = 0.5;
l1 = 0.4;
l2 = 0.28;

theta1 = pi/2;
theta2 = -pi/2;

dottheta1 = 0;
dottheta2 = 0;

tau = [2;1];

M = [(m1+m2)*l1^2 + m2*l2^2 + 2*m2*l1*l2*cos(theta2), m2*l2^2 + m2*l1*l2*cos(theta2);
    m2*l1*l2*cos(theta2)+m2*l2^2, m2*l2^2];
H = [-m2*l1*l2*sin(theta2)*dottheta2^2 - 2*m2*l1*l2*sin(theta2)*dottheta1*dottheta2;
    m2*l1*l2*sin(theta2)*dottheta1^2];
G = [m2*l2*g*cos(theta1+theta2) + (m1+m2)*l1*g*cos(theta1);
    m2*l2*g*cos(theta1+theta2)];

ddotq = M\(tau - H - G)