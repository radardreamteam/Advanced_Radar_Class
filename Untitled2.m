clc,clear all,close all
%% �����ź�����
N=100;M=10;%�ź�������Ŀ����Ԫ����
K=2;%��Դ����
theta=[-10;40]*pi/180;
SNR=[10;20];sigma=1;
Am=sqrt(2*sigma^2*10.^(SNR/10));
% Am=[sqrt(10.^(SNR/10))];
S=Am*ones(1,N);
S(2,:)=S(2,:).*exp(1i*2*pi*rand(1,N));
for a=1:M
for b=1:K
A(a,b)=exp(-1i*(a-1)*pi*sin(theta(b)));%�� b �ж�Ӧ�Ķ��� theta(b)
end
end
V=zeros(M,N);
for m=1:M
v=wgn(1,N,0,'complex');
v=v-mean(v);
v=v/std(v);
V(m,:)=v;
end
X=A*S+V;
%% ���ý������ݹ����źŵĿռ���ؾ��� R
R=zeros(M,M);
for i=1:N
R=R+X(:,i)*X(:,i)';
end
R=R/N;%��һ��ͳ��ƽ��
%MUSIC �㷨
[VR,D]=eig(R);
D=real(D);
[B,IX]=sort(diag(D));
G=VR(:,IX(M-K:-1:1));
MUSICP=[];
for n=-pi/2:pi/180:pi/2
a=exp(-1i*[0:M-1]'*pi*sin(n));
MUSICP=[MUSICP,1/(a'*G*G'*a)];
MUSICP=real(MUSICP);end
n=length(MUSICP);
maxx=max(MUSICP);
figure,plot(-90:1:90,10*log10((MUSICP+eps)/maxx)+3.5),axis([-90,90,-60,inf]),title('MUSIC �㷨');
%RootMUSIC �㷨
syms z
pz=z.^([0:M-1]');
pz1=(z^(-1)).^([0:M-1]);
fz=z^(M-1)*pz1*G*G'*pz;
a=sym2poly(fz);
r=roots(a);
r1=abs(r);
for i=1:2*K %ÿ���ź�Դ�� K ��
[Y,I(i)]=min(abs(r1-1));
r1(I(i))=inf;
end
for i=1:2*K
theta_esti(i)=asin(-angle(r(I(i)))/pi)*180/pi;
end
%ESPRIT �㷨
S=VR(:,IX(M:-1:M-K+1));
S1=S(1:M-1,:);
S2=S(2:M,:);
fai=S1\S2;
[U_fai,V_fai]=eig(fai);
for i=1:K
ESPRITtheta_esti(i)=asin(-angle(V_fai(i,i))/pi)*180/pi;
end
%MVDR �㷨
MVDRP=[];
for n=-pi/2:pi/180:pi/2
a=exp(-1i*[0:M-1]'*pi*sin(n));
MVDRP=[MVDRP,1/(a'*inv(R)*a)];
end
n=length(MVDRP);
maxx=max(MVDRP);
figure,plot(-90:1:90,10*log10((MVDRP+eps)/maxx)+3.5),axis([-90,90,-35,inf]),title('MVDR');
%F-SAPES �㷨
P=6;%������ĿL=M+1-P;%������Ԫ��Ŀ�������� M-1
Rf=zeros(L,L);
for i=1:P
Rf=Rf+X(i:i+L-1)*X(i:i+L-1)'/N;
end
Rf=Rf/P; %����ƽ����Ŀռ���ؾ���
n1=0:P-1;
n2=0:L-1;
cc=[1 zeros(1,L-1)];
for n3=-90:.5:90
fy=exp(1i*pi*sin(n3/180*pi));
tt=[(fy.^(n1')).' zeros(1,M-P)];
Tfy=toeplitz(cc,tt);
GfTheta=1./(P^2)*Tfy*R*Tfy';
Qf=Rf-GfTheta;
aTheta=fy.^(-n2');
Wof=(Qf\aTheta)./(aTheta'*(Qf\aTheta));
sigma2sTheta(((n3+90)/.5+1))=Wof'*GfTheta*Wof;
end
maxx=max(sigma2sTheta);
figure,plot(-90:.5:90,10*log10((sigma2sTheta+eps)/maxx)+3.5),axis([-90,90,-35,inf]),title('F-SAPES')
