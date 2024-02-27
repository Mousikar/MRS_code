% formation test 没有时滞

dt = 0.001;
iter = 1000*10;
q = [ 0.5; 0.7];
k = 3;
L = [1 -1; -1 1];

q_his = zeros(2, 1+iter);
err_his = zeros(2, 1+iter);

for i = 1:iter
    t = i*dt;
    q_his(:,i) = q;
    err_his(:,i) = ( q - ref(t) ) ;
    dotq = dot_ref(t) - k * ( q - ref(t) ) - L * ( q - ref(t) );
    q = q + dotq * dt;
end
t = (i+1)*dt;
q_his(:,1+i) = q;
err_his(:,1+i) = ( q - ref(t) ) ;

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