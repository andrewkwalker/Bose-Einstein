%% Preamble

% Andrew Walker
% AMATH 481A - Autumn 2013
% Final Project - Bose-Einstein Condensate

close all;
clear all;

%% Initialize Variables

tspan = 0:.5:4;
n = 16;
L = 2*pi;
x = linspace(-L/2, L/2, n+1);
xspan = x(1:n);
yspan = xspan;
zspan = yspan;

[X,Y,Z] = meshgrid(xspan,yspan,zspan);

A1 = -1;
A2 = -1;
A3 = -1;

B1 = -A1;
B2 = -A2;
B3 = -A3;

% Create the K matrix

kx = (2*pi/L)*[0:(n/2-1) (-n/2):-1];
ky = kx;
kz = ky;

[KX,KY,KZ] = meshgrid(kx,ky,kz);

K = KX.^2+KY.^2+KZ.^2;

% Initial Condition psi=cos(x)cos(y)cos(z)
psi1_init = cos(X).*cos(Y).*cos(Z);

% Initial Condition psi=sin(x)sin(y)sin(z)
psi2_init = sin(X).*sin(Y).*sin(Z);

psi1_init_f = fftn(psi1_init);
psi2_init_f = fftn(psi2_init);

psi1_init_vec_f = reshape(psi1_init_f,n^3,1);
psi2_init_vec_f = reshape(psi2_init_f,n^3,1);

%% Use ODE45 to time-step

[t1,psi1_sol_vec_f] = ode45(@(t,psi_vec_f) rhsfft(t,psi_vec_f,K,A1,...
                                          A2,A3,B1,B2,B3,X,Y,Z,n),...
                          tspan,psi1_init_vec_f);

[t2,psi2_sol_vec_f] = ode45(@(t,psi_vec_f) rhsfft(t,psi_vec_f,K,A1,...
                                          A2,A3,B1,B2,B3,X,Y,Z,n),...
                          tspan,psi2_init_vec_f);

% A1 = real(psi1_sol_vec_f);
% A2 = imag(psi1_sol_vec_f);
% A3 = real(psi2_sol_vec_f);
% A4 = imag(psi2_sol_vec_f);
% 
% save A1.dat A1 -ascii
% save A2.dat A2 -ascii
% save A3.dat A3 -ascii
% save A4.dat A4 -ascii
                      
%% Visualize Solution

writer1 = VideoWriter('psi1.avi');
writer1.FrameRate = 3;
open(writer1);

for j=1:length(t1)
    psi1_sol_f = reshape(psi1_sol_vec_f(j,:),n,n,n);
    psi1_sol = ifftn(psi1_sol_f);
    figure(j);
    isosurface(real(psi1_sol));
    drawnow;
    writeVideo(writer1,getframe);
    pause(0.05);
end

close(writer1);

close all;

writer2 = VideoWriter('psi2.avi');
writer2.FrameRate = 3;
open(writer2);

for j=1:length(t2)
    psi2_sol_f = reshape(psi2_sol_vec_f(j,:),n,n,n);
    psi2_sol = ifftn(psi2_sol_f);
    figure(j+9);
    isosurface(real(psi2_sol));
    drawnow;
    writeVideo(writer2,getframe);
    pause(0.05);
end

close(writer2);
%
% writer3 = VideoWriter('psi3.avi');
% writer3.FrameRate = 5;
% open(writer3);
% 
% for j=1:length(t2)
%     psi2_sol_f = reshape(psi2_sol_vec_f(j,:),n,n,n);
%     psi2_sol = ifftn(psi2_sol_f);
%     figure(j);
%     psi2_sol_conj = zeros(n,n,n);
%     for j=1:n
%         for k=1:n
%             for l=1:n
%                 psi2_sol_conj(j,k,l)=real((psi2_sol(j,k,l))^2)+...
%                                      imag((psi2_sol(j,k,l))^2);
%             end
%         end
%     end
%     isosurface(psi2_sol_conj.*psi2_sol);
%     drawnow;
%     writeVideo(writer3,getframe);
%     pause(0.05);
% end
% 
% close(writer3);