clear all; close all;
%%%%%%%% MUSIC for Uniform Circle Array%%%%%%%%
r=1;                                            % 半径(m)
N=16;                                           % 阵元数目
d=2*r*sin(pi/N);
 
M=1;                                            % 信源数量
gamma=2*pi/N*(0:N-1);                           % 和参考阵元的角度
fc=16000;                                       % 采样率
c=340;                                          % 声音的速度（m/s）
lambda=c/fc;
 
a_theta=[45 ];                                     % 仰角
a_phi=[60 ];                                       % 方位角
 
zeta=2*pi/lambda*r*sin(a_theta*pi/180);
A=exp(1i*zeta*cos((a_phi-gamma)*pi/180)).';     % 导向矢量
K=100;                                          % 快拍数
t=(0:K-1)/1000;
% 构建信号模型
S=sin(2*pi*fc*t);
X=A*S;
snr=10;
X1=awgn(X,snr,'measured');
% 计算协方差矩阵
R=X1*X1'/K;
% 特征值分解
[Q,D]=eig(R);
[D,I]=sort(diag(D),1,'descend');
Q=Q(:,I);
Qn=Q(:,M+1:N);                                  % 信号子空间
 
theta=0:90;
phi=0:1:360;
 
p_MUSIC=zeros(length(theta),length(phi));
 
for ii=1:length(theta)
    for iii=1:length(phi)
        zeta=2*pi/lambda*r*sin(theta(ii)*pi/180);
        A=exp(1i*zeta*cos((phi(iii)-gamma)*pi/180)).';
        p_MUSIC(ii,iii)=(1/(A'*(Qn*Qn')*A));
    end
end
mesh(phi,theta,abs(p_MUSIC))
grid on;xlabel('仰角');ylabel('方位角');zlabel('PMUSIC');title('UCA MUSIC');

