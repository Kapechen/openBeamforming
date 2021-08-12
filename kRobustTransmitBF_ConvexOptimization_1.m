%%
%Author:            Kapechen
%Date:              2021/08/09-2021/08/12
%Description:       Beamforming simulation based on convex optimization--1
%Version:           0.2
%Dependent Tools:   CVX:http://cvxr.com/cvx/
%%
clear all
%%
%%%信道状态信息矩阵初始化
%SU-Tx ULA阵元数量
SU_Tx_M = 6;
%ULA波长间距比
sp = 0.5;
%初始化SU-TX到SU-RX的导向矢量as_theta 用于计算Hs
as_theta = zeros(SU_Tx_M,1);
%设定波达角
theta_DOA_SURX = -20;
for i=1:1:SU_Tx_M
    as_theta(i,1) = exp(-1j*2*pi*sp*(i-1)*sin(deg2rad(theta_DOA_SURX)));
end
Hs = zeros(SU_Tx_M,1);
for i=1:1:SU_Tx_M
    Hs(i,1) =as_theta(i,1);
end
%%
%%%约束条件变量初始化
%定义多径总数 N>= M-1
N = 7;
%保护PU的干扰温度 单位:dB
P = -60;
A = 10^(P/10);
%估计误差
EpsilonP = sqrt(A)/2;
%PU可侦测移动范围 单位:° 以及对应SU-TX到PU的导向矢量ap_theta
theta1 = 50;
thetaN = 70;
thetaPi = linspace(theta1,thetaN,N);
ap_theta = zeros(SU_Tx_M,N);
for i = 1:1:N
    for j = 1:1:SU_Tx_M
        ap_theta(j,i) = exp(-1j*2*pi*sp*(j-1)*sin(deg2rad(thetaPi(i))));
    end
end

%%
%%%求解凸优化
%初始化权值向量w
cvx_begin
    variable w(SU_Tx_M,1) complex
    maximise(real(Hs'*w))
    subject to
        imag(Hs'*w)== 0
        for i=1:1:N
            abs(ap_theta(:,i)'*w) <= sqrt(A)-EpsilonP
        end
        norm(w,2) <= 1
cvx_end
%%
%%%求解结果验证
%计算阵因子部分
AF_theta = [-90:1:90];
length1 = length(AF_theta);
AF = zeros(length1,1);
eMatrix = zeros(SU_Tx_M,1);
for i=1:1:length1
      for j=1:1:SU_Tx_M
            psix = sp*2*pi*sin(pi+deg2rad(AF_theta(i)));
            eMatrix(j,1) = exp(1j*(j-1)*psix);
      end
     AF(i,1) = w'*eMatrix;
end
%%
%绘图部分
figure(1)
stem(abs(w));
figure(2)
%取模再转成分贝
POUT = 20*log10(abs(AF));
plot(AF_theta,POUT);
xlabel('角度(°)');
ylabel('发射功率/dB')