%clear everything
%close all; clear all;
%load the data

%chose which data you would like to use

%load('Beamform/NO_CLEAN/ECA/CFAR_files.mat');
%load('Beamform/NO_CLEAN/Wiener/CFAR_files.mat');
%load('Beamform/NO_CLEAN/NLMS/CFAR_files.mat');
%load('Beamform/NO_CLEAN/FBLMS/CFAR_files.mat');
%load('Beamform/NO_CLEAN/RLS/CFAR_files.mat');
load('myData/CFAR_files.mat');

%load('Element-by-element/ECA/CFAR_files.mat');
%load('Element-by-element/Wiener/CFAR_files.mat');
%load('Element-by-element/NLMS/CFAR_files.mat');
%load('Element-by-element/FBLMS/CFAR_files.mat');
%load('Element-by-element/RLS/CFAR_files.mat');

%load('Beamform/CLEAN/ECA/CFAR_files.mat');

%graph the data to visualize it
%figure 1 which is figure 5 from the rev7
% figure(1);
% contourf(X*1e-3,Y,rdmap_compansated');
% title("Range Doppler Map");
% ylabel("Doppler shift (Hz)");
% xlabel("Distance (Km)");

%figure 2 try to make the 3D mesh
% figure(2);%coresponds to figure(4) in rev7
% mesh(X*1e-3,Y,(rdmap_compansated')); axis tight
% title("3D mesh to get a better understanding");
% xlabel('Range (Km)')
% ylabel('Doppler shift (Hz)')

%attempt dimension reduction
r=5;
rdmap_reduced=imresize(rdmap_compansated, 1/r);%this dimension reduction works

X_OG=linspace(0,5.373801981617511e+04,2152);%current one
X_reduced=linspace(0,5.373801981617511e+04,size(rdmap_reduced,1));
X_reduced=repmat(X_reduced,[size(rdmap_reduced,2),1]);

Y_OG=linspace(500,-500,501);
Y_reduced=linspace(500,-500,size(rdmap_reduced,2));
Y_reduced=repmat(Y_reduced',[1,size(rdmap_reduced,1)]);

% figure(3);
% contourf(X_reduced*1e-3,Y_reduced,rdmap_reduced');
% title("Range Doppler Map_reduced");
% ylabel("Doppler shift (Hz)");
% xlabel("Distance (Km)");

%figure 2 try to make the 3D mesh
% figure(4);%coresponds to figure(4) in rev7
% mesh(X_reduced*1e-3,Y_reduced,(rdmap_reduced')); axis tight
% title("3D mesh to get a better understanding_reduced");
% xlabel('Range (Km)')
% ylabel('Doppler shift (Hz)')