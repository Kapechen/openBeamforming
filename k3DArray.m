clear,clc
%三维天线阵
%x、y、z轴阵元数量初始化
N1=8;
N2=8;
N3=8;
%定义波长 cm和dx dy dz
lambda = 10;
dx = 0.5*lambda;
dy = 0.5*lambda;
dz = 0.5*lambda;
%计算k
k = 2*pi/lambda;
%分别初始化x轴、y轴、z轴权向量
wx = zeros(1,N1);
wy = zeros(1,N2);
wz = zeros(1,N3);
%凯瑟-贝塞尔权值并归一化
wx = kaiser(N1,0.5);
wy = kaiser(N2,0.5);
wz = kaiser(N3,0.5);
wx = mapminmax(wx',0,1);
wy = mapminmax(wy',0,1);
wz = mapminmax(wz',0,1);
% wx = wx';
% wy = wy';
% wz = wz';
figure(1)
stem(wx)
title('三维天线阵权值(x轴、凯瑟-贝塞尔权值)');
xlabel('阵元数量');
ylabel('归一化权值');
figure(2)
stem(wy)
title('三维天线阵权值(y轴、凯瑟-贝塞尔权值)');
xlabel('阵元数量');
ylabel('归一化权值');
figure(3)
stem(wz)
title('三维天线阵权值(z轴、凯瑟-贝塞尔权值)');
xlabel('阵元数量');
ylabel('归一化权值');
%角度参数设置
theta = -90:1:90;
phi = 0:1:360;
%波束调向参数设置
theta0 = 30;
phi0 = 0;
%调向相位延迟βx βy βz计算
betax = -k*dx*sin(deg2rad(theta0))*cos(deg2rad(phi0));
betay = -k*dy*sin(deg2rad(theta0))*sin(deg2rad(phi0));
betaz = -k*dz*cos(deg2rad(theta0));
%循环计数辅助值初始化
length1 = length(theta);
length2 = length(phi);
%初始化阵因子AF
AF = zeros(length1,length2);
for i=1:1:length1
    for j=1:1:length2
        %计算ψx，ψy ψz向量
        psix = k*dx*sin(deg2rad(theta(i)))*cos(deg2rad(phi(j)))+betax;
        psiy = k*dy*sin(deg2rad(theta(i)))*sin(deg2rad(phi(j)))+betay;
        psiz = k*dz*cos(deg2rad(theta(i)))+betaz;
        %初始化指数向量组
        ex = zeros(N1,1);
        ey = zeros(N2,1);
        ez = zeros(N3,1);
        for n1 = 1:1:N1
            ex(n1) = exp(1j*(n1-1)*psix);
        end
        for n2 = 1:1:N2
            ey(n2) = exp(1j*(n2-1)*psiy);
        end
        for n3 = 1:1:N3
            ez(n3) = exp(1j*(n3-1)*psiz);
        end
        AF(i,j) = (wx * ex) * (wy * ey) * (wz * ez);
    end
end
%取模
for i = 1:1:length1
    for j=1:1:length2
        AF(i,j)=abs(AF(i,j));
    end
end
figure(4)
plot(theta,AF(:,1));
title('波束调向三维天线阵的阵因子垂直方向图(θ0=30，φ0=0)')
xlabel('θ')
ylabel('|AF|')
%球坐标转化为直接坐标后画图
[X,Y,Z]=sph2cart(deg2rad(phi),(deg2rad(theta))',AF);
figure(5)
mesh(X,Y,Z);
title('波束调向三维天线阵的三维阵因子方向图(θ0=30，φ0=0)')
xlabel('x轴')
ylabel('y轴')
zlabel('z轴')
%剖面AF图
% figure(2)
% plot(theta,AF(:,1));

