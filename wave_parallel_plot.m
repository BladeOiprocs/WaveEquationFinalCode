close all
clear
clc

filename = "dataParallel.csv";
data = csvread(filename);
data(:,end) = [];

[L,N] = size(data);
M = N;
O = L/M;



U = zeros(M,N,O);
for k = 0:O-1
    idx1 = k*M+1;
    idx2 = (k+1)*M;
    U(:,:,k+1) = data(idx1:idx2,1:N);
end

FIG = figure();
X=0:.1:(0.1*(M-1));
Y=0:.1:(0.1*(N-1));
% [X,Y] = meshgrid(x,y);
for k=1:O
    surf(X,Y,U(:,:,k));
     view(0,90)
    zlim([-N/2 N/2]);
    drawnow
    pause(0.01);
end