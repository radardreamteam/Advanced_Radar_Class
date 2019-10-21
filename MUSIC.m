clear all; close all;
%%%%%%%% MUSIC for Uniform Circle Array%%%%%%%%
r=1;                                            % �뾶(m)
N=16;                                           % ��Ԫ��Ŀ
d=2*r*sin(pi/N);
 
M=1;                                            % ��Դ����
gamma=2*pi/N*(0:N-1);                           % �Ͳο���Ԫ�ĽǶ�
fc=16000;                                       % ������
c=340;                                          % �������ٶȣ�m/s��
lambda=c/fc;
 
a_theta=[45 ];                                     % ����
a_phi=[60 ];                                       % ��λ��
 
zeta=2*pi/lambda*r*sin(a_theta*pi/180);
A=exp(1i*zeta*cos((a_phi-gamma)*pi/180)).';     % ����ʸ��
K=100;                                          % ������
t=(0:K-1)/1000;
% �����ź�ģ��
S=sin(2*pi*fc*t);
X=A*S;
snr=10;
X1=awgn(X,snr,'measured');
% ����Э�������
R=X1*X1'/K;
% ����ֵ�ֽ�
[Q,D]=eig(R);
[D,I]=sort(diag(D),1,'descend');
Q=Q(:,I);
Qn=Q(:,M+1:N);                                  % �ź��ӿռ�
 
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
grid on;xlabel('����');ylabel('��λ��');zlabel('PMUSIC');title('UCA MUSIC');

