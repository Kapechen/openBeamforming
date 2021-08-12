clear,clc
%直线阵且为各向同性阵元
%阵元数量(偶数)
N = 16;
%汉明归一化权值
yout = mapminmax(hamming(N)',0,1);
w = yout(:,((N/2)+1):N);
figure(1)
stem(yout);
title('汉明天线阵权值');
xlabel('阵元数量');
ylabel('归一化阵权值');
%求分母
fenmu = sum(w);
%用于暂存AF的值
list = [];
%调向方位角
theta0 = 0;
%规定角度范围
for theta = -90:0.1:90
    u = 1/2*pi * (sin(deg2rad(theta))-sin(deg2rad(theta0))); %d = 0.5λ
    fenzi = 0;
    for n = 1:(N/2)
        fenzi = fenzi + w(n)*cos((2*n-1)*u); 
    end
    %计算加权后归一化阵因子
    AF = fenzi/fenmu;
    list = [list AF];
end
figure(2)
plot(-90:0.1:90,list);
title('具有汉明权值的阵因子图')
xlabel('θ');
ylabel('|AF|');
axis([-90 90 -0.5 1]);

