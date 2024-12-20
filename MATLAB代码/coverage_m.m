% coverage
% 区域是10m×10m的
clear
%%
n=5;
l1=-5;
l2=5;
reso=0.1;
a = l1:reso:l2; % 分辨率是1mm
b = l1:reso:l2;
% 设置5个重心
rng(3);% 随机种子88:21; 3:21; 7:13
mu = l1 + (l2-l1) *rand(2,5);
figure
plot(mu(1,:),mu(2,:),'ro')
hold on;
% 初始化机器人的位置
p = l1 + (l2-l1)/2 *rand(2,n);
plot(p(1,:),p(2,:),'bo');
axis([l1,l2,l1,l2])


phi = zeros(size(a,2),size(b,2));
for i = 1:size(a,2)
    for j = 1:size(b,2)
        q =[a(i);b(j)];
        phi(i,j)=getphi(q,mu);
    end
end

figure
[X,Y]=meshgrid(a,b);
surf(X,Y,phi)
view(2)

%% 迭代
V = zeros(size(a,2),size(b,2),size(p,2));
iter = 21;
k_v = 0.5;
col = [
    0.5358    0.0393    0.0899
    0.4611    0.5479    0.0134
    0.9900    0.6800    0.7600
    0.0000    0.6800    0.9200
    0.9800    0.7500    0.0000
    0.5955    0.2201    0.9358
    ];
figure
for k = 1:iter
    %% 计算分区
    for i = 1:size(a,2)
        for j = 1:size(b,2)
            q =[a(i);b(j)];
            dis = sum((p-q).^2);
            index = find(dis==min(dis));
            V(i,j,index) = phi(i,j);
            plot(a(i),b(j),'o','Color',col(index,:))
            hold on
        end
    end
    % 计算质心
    MV = zeros(n,1);
    CV = zeros(2,n);
    for index = 1:n
        MV(index) = sum(sum(V(:,:,index)));
    end
    for index = 1:n
        for i = 1:size(a,2)
            for j = 1:size(b,2)            
                CV(1,index) = CV(1,index) + 1/MV(index)*a(i)*V(i,j,index);
                CV(2,index) = CV(2,index) + 1/MV(index)*b(j)*V(i,j,index);
            end
        end
        plot(CV(1,index),CV(2,index),'.r','MarkerSize',20)
        hold on
        plot(p(1,:),p(2,:),'.b','MarkerSize',20);
    end
    drawnow
    % pause(1)
    hold off
    % 计算控制
    u=-k_v*(p-CV);
    p=p+u;
    disp('---------------')
    disp(k)
    disp(CV)
    disp(p)
end

%%
figure
plot(mu(1,:),mu(2,:),'ro')
hold on;
plot(p(1,:),p(2,:),'bo');
plot(CV(1,:),CV(2,:),'.r','MarkerSize',20)
axis([l1,l2,l1,l2])

% mu=[0,0];% 均值向量
% Sigma=[1 0.8;0.8 1];% 协方差矩阵
% [X,Y]=meshgrid(-3:0.1:3,-3:0.1:3);%在XOY面上，产生网格数据
% p=mvnpdf([X(:) Y(:)],mu,Sigma);%求取联合概率密度，相当于Z轴
% p=reshape(p,size(X));%将Z值对应到相应的坐标上
% 
% figure
% set(gcf,'Position',get(gcf,'Position').*[1 1 1.3 1])
% 
% subplot(2,3,[1 2 4 5])
% surf(X,Y,p),axis tight,title('二维正态分布图')
% subplot(2,3,3)
% surf(X,Y,p),view(2),axis tight,title('在XOY面上的投影')
% subplot(2,3,6)
% surf(X,Y,p),view([0 0]),axis tight,title('在XOZ面上的投影')

% 密度函数，有几个重心就有几个μ
function result = getphi(q,mu)
sigma = 0.5;
result = 0;
for i = 1:size(mu,2)
result = result + 1/sigma/sqrt(2*pi)*...
    exp(-norm(q-mu(:,i))^2/2/sigma/sigma);
end
end
