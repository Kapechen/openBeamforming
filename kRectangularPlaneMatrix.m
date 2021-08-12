clear,clc
%矩形平面阵NxM
N=8;
M=8;
%定义波长 mm和dx dy
lambda = 10;
dx=0.5*lambda;
dy=0.5*lambda;
%计算k
k=2*pi/lambda;
%角度参数设置
theta = -90:1:90;
phi = 0:1:360;
%波束调向参数设置 单位：度
theta0 = 30;
phi0 = 0;
%计算相位延迟 βx与βy
betax = -k*dx*sin(deg2rad(theta0))*cos(deg2rad(phi0));
betay = -k*dy*sin(deg2rad(theta0))*sin(deg2rad(phi0));
%权值初始化
a = ones(M,1);
b = ones(N,1);
w=zeros(M,N);
%凯瑟-贝塞尔权值
a = kaiser(M,0.5);
b = kaiser(N,0.5);
%归一化
a = mapminmax(a',0,1);
b = mapminmax(b',0,1);
%画全职图
figure(1)
stem(a);
title('矩形平面阵权值(x轴、凯瑟-贝塞尔权值)');
xlabel('阵元数量');
ylabel('归一化权值');
figure(2);
stem(b);
title('矩形平面阵权值(y轴、凯瑟-贝塞尔权值)');
xlabel('阵元数量');
ylabel('归一化权值');
%阵因子初始化
AF = zeros(length(theta),length(phi));
for i = 1:1:length(theta)
    for j=1:1:length(phi)
        for m=1:1:M
            for n=1:1:N
                w(m,n) = a(m)*b(n);
                temp = w(m,n)*exp(1j*((m-1)*(k*dx*sin(deg2rad(theta(i)))*cos(deg2rad(phi(j)))+betax)+(n-1)*(k*dy*sin(deg2rad(theta(i)))*sin(deg2rad(phi(j)))+betay)));
                AF(i,j) = AF(i,j)+temp;
            end
        end
    end
end
for i = 1:1:length(theta)
    for j=1:1:length(phi)
        AF(i,j)=abs(AF(i,j));
    end
end
figure(3)
plot(theta,AF(:,1));
title('波束调向矩形平面阵的阵因子垂直方向图(θ0=30，φ0=0)')
xlabel('θ')
ylabel('|AF|')
[X,Y,Z]=sph2cart(deg2rad(phi),(deg2rad(theta))',AF);
figure(4)
mesh(X,Y,Z);
title('波束调向矩形平面阵的三维阵因子方向图(θ0=30，φ0=0)')
xlabel('x轴')
ylabel('y轴')
zlabel('z轴')
