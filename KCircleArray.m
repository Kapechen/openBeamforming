clear,clc
%圆阵且为各向同性阵元
%阵元数量
N=10;
%均匀加权
w=ones(16,1);%列向量
figure(1)
stem(w);
title('天线阵权值(均匀加权)');
xlabel('阵元数量');
ylabel('归一化权值');
%定义波长 mm
lambda = 10;
k=2*pi/lambda;
%圆阵半径
a=lambda;
%阵列相角调整角度 单位：度
theta0=45;
phi0=0;
%初始化theta phi
theta=-90:1:90;
phi=0:1:360;
length1 =length(theta);
length2 =length(phi);
%初始化阵因子
AF = zeros(length1,length2);
for i=1:1:length1
    for j=1:1:length2
        for n = 1:1:N
            phiN = (2*pi*(n-1))/N;
            tempN = w(n)*exp(-1j*k*a*(sin(deg2rad(theta(i))*cos(deg2rad(phi(j))-phiN))-sin(deg2rad(theta0))*cos(deg2rad(phi0)-phiN)));
            AF(i,j) = AF(i,j)+tempN;
        end
    end
end
for i=1:1:length1
    for j=1:1:length2
        AF(i,j)=abs(AF(i,j));
    end
end
figure(2)
plot(theta,AF(:,1))
title('波束调向圆阵的阵因子垂直方向图(θ0=45,φ0=0)')
xlabel('θ')
ylabel('|AF|')
[X,Y,Z]=sph2cart(deg2rad(phi),(deg2rad(theta))',AF);
figure(3)
mesh(X,Y,Z)
title('波束调向圆阵的三维阵因子方向图(θ0=45,φ0=0)')
xlabel('x轴');
ylabel('y轴');
zlabel('z轴');


