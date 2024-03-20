% formation test 有均匀时滞

dt = 0.001;
iter = 1000*10;
q = [ 0.5; 0.7];
k = 3;
L = [1 -1; -1 1];
delay = 10;

q_his = zeros(2, delay+1+iter);
err_his = zeros(2, delay+1+iter);

for i = 1:delay
    q_his(:,i) = q;
    err_his(:,i) = ( q - ref(0) ) ;
end

for i = 1+delay:iter+delay
    t = i*dt;
    q_his(:,i) = q;
    err_his(:,i) = ( q - ref(t) ) ;
    err_his(:,i) = L * ( q - ref(t) ) ;
    q_d = ref(t);
    dotq_d = dot_ref(t);
    q_temp = q_his(:,i-delay);
    ref_temp = ref(t-delay*dt);
    dotq(1) = dotq_d(1) - k * ( q(1) - q_d(1) ) - ...
        ( ( q(1) - q_d(1) ) - ( q_temp(1) - q_temp(1) ) );
    dotq(2) = dotq_d(2) - k * ( q(2) - q_d(2) ) - ...
        ( ( q_temp(1) - q_temp(1) ) - ( q(2) - q_d(2) ) );
    % dotq = dot_ref(t) - k * ( q - ref(t) ) - L * ( q_his(:,i-delay) - ref(t-delay*dt) );
    q = q + dotq * dt;
end
t = (i+1)*dt;
q_his(:,1+i) = q;
err_his(:,1+i) = L * ( q - ref(t) ) ;

figure
plot(q_his')
figure
plot(err_his')
%%
q - ref(t)
%%
function q_d = ref(t)
q_d = [ 3*sin(t);
        2*sin(t)];
q_d = [ 3*t;
        2*t];
end

function dotq_d = dot_ref(t)
dotq_d = [ 3*cos(t);
        2*cos(t)];
dotq_d = [ 3;
        2];
end