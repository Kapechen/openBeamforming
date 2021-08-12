clear,clc
m=16; %阵元数
theta = [0 20]; %仰角为0 方位角为20度
k=0.5; %k=d/λ = 0.5
SampleNum=500; %样本数量
L=1000;%采样精度
MeanN=0; %噪声均值
VarN=1;  %噪声方差
SNR=10;  %信噪比
INR=10;  %信干比
RVARsig = sqrt(VarN)*10^(SNR/20); %信号方差
RVARn = sqrt(VarN)*10^(INR/20);   %干扰方差
%信号生成
s = [RVARsig*exp(1i*2*pi*50*0.001*[0:SampleNum-1])
    RVARn*exp(1i*2*pi*(100*0.001*[0:SampleNum-1]+rand))];
%调向方位角
theta0=20;
%生成A矩阵
A=exp(-1i*2*pi*k*[0:m-1]'.*(sin(deg2rad(theta))-sin(deg2rad(theta0))));
%生成噪声
e=sqrt(VarN/2)*(randn(m,SampleNum)+1i*randn(m,SampleNum));
%生成ULA数据
Y=A*s+e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(1:1:500,Y(1,:))
%初始化权重系数
% %均匀权值
% w = ones(m,1);
%汉明加权
yout = mapminmax(blackman(m)',0,1);
w = yout';

%进行波束成形
beam=zeros(1,L);
%设定调向方位角theta0
theta0=0;
for i=1:L
 a=exp(-1i*2*pi*k*[0:m-1]'.*(sin(-pi/2+pi*(i-1)/L)-sin(deg2rad(theta0))));
 beam(i)=20*log10(abs(w'*a));
end
%绘图
figure(1)
stem(w);
title('16布莱克曼加权权值')
axis([1 16 0 2])
figure(2)
angle = -90:180/L:(90-180/L);
plot(angle,beam);
axis([-90 90 -50 30])
title('布莱克曼16阵元ULA方向图')
xlabel('方向角/度');
ylabel('幅度响应/dB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%