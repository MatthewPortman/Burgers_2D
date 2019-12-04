%% Author: Matthew Portman - Due December 22
%% Program: Solving Burger's Equation in 2D

clear; close all; 
%% Definitions
% Defining the spatial step, time step for LW and FTCS as well as a
% separate step for a different implementation of the mimetic method. We
% also define some convenient parameters.
N = 500
Nt = 250
Ntt = 10;

N1 = N-1;
N11 = N+2;

% Grid max and min.
xMax = 10;
xMin = -10;

% Step width in space and time as well as convenient parameters.
dx = (xMax - xMin)/N1
dt = 0.2;
dx1 = 1/dx;
dxdx = 1/(dx*dx);

% Setting up the grids.
xgrid = xMin : dx: xMax;
tgrid = 0 : dt : Nt*dt;

% Initial Condition/Velocity distribution.
theta = pi/6;
sigma = 2;

alpha = cos(theta);
beta = sin(theta);

x0 = 0;
y0 = 0;
U1 = zeros(N,N);

IC = @(x,y) 0.5/(pi*sigma^2)*exp(-0.5/sigma^2*((x-x0)^2 + (y-y0)^2));

for j = 1:N
    for i = 1:N
        U1(i,j) = IC(xgrid(i),xgrid(j)); 
    end
end
UU = zeros(N,N,Nt);
UU(:,:,1) = U1;
UFTCS = UU;

%% Forward Time Centered Space
for n = 1:Nt-1
    for j = 2:N1
        for i = 2:N1
            Ux = 0.5*dx1*(UFTCS(i+1,j,n) - UFTCS(i-1,j,n));
            Uy = 0.5*dx1*(UFTCS(i,j+1,n) - UFTCS(i,j-1,n));
            A = alpha*UFTCS(i,j,n);
            Ax = alpha*Ux;
            Ay = alpha*Uy;
            B = beta*UFTCS(i,j,n);
            Bx = beta*Ux;
            By = beta*Uy;
              UFTCS(i,j,n+1) = UFTCS(i,j,n) - dt*(A*Ux + B*Uy) + 0.5*dt*dt*(A*(Ax*Ux) + ...
                 A*(Bx*Uy) + B*(Ay*Ux) + B*(By*Uy));
        end
    end
    % Calculating the area under the curve for every time step.
    avgFTCS(n) = norm(trapz(UFTCS(:,:,n)));
end

% Finding the standard deviation between all of those areas.
ConserveFTCS = std(avgFTCS)

% Generating the movie. This could be optimized by throwing it in the outer
% loop but that's for another time.

 Nstring = num2str(N);
 videoname = strjoin({Nstring,'_burgers_FTCS.avi'});

 Burgers = VideoWriter(videoname,'Uncompressed AVI');
 open(Burgers)
  for i = 1:Nt-1
      [X,Y] = meshgrid(xgrid,xgrid);
      drawnow
      xlabel('x'), ylabel('y'), zlabel('u');
      surf(X,Y,UU(:,:,i)');
      view(60,30);
      title(['Burgers Equation in 2D -- FTCS scheme \newline t = ' num2str(i*dt) ', N = ' num2str(N)]);
      F(i) = getframe(gcf);
      writeVideo(Burgers,F(i).cdata)
  end
 close(Burgers);

%% Lax-Wendroff -- Same format as FTCS
%% Thanks to Jakob Beran for the code/algorithm structure

for n = 1:Nt-1
    for j = 2:N1
        for i = 2:N1
            Ux = 0.5*dx1*(UU(i+1,j,n) - UU(i-1,j,n));
            Uy = 0.5*dx1*(UU(i,j+1,n) - UU(i,j-1,n));
            Uxx = dxdx*(UU(i+1,j,n) - 2*UU(i,j,n) + UU(i-1,j,n));
            Uyy = dxdx*(UU(i,j+1,n) - 2*UU(i,j,n) +  UU(i,j-1,n));
            Uxy = 0.25*dxdx*(UU(i+1,j+1,n) - UU(i-1,j+1,n) - UU(i+1,j-1,n) + ...
                UU(i-1,j-1,n));
              A = alpha*UU(i,j,n);
              Ax = alpha*Ux;
              Ay = alpha*Uy;
              B = beta*UU(i,j,n);
              Bx = beta*Ux;
              By = beta*Uy;
              UU(i,j,n+1) = UU(i,j,n) - dt*(A*Ux + B*Uy) + 0.5*dt*dt*(A*(Ax*Ux+A*Uxx) + ...
                 A*(Bx*Uy+B*Uxy) + B*(Ay*Ux+A*Uxy) + B*(By*Uy+B*Uyy));
        end
    end
    avgLW(n) = norm(trapz(UU(:,:,n)));
end

ConserveLW = std(avgLW)

Nstring = num2str(N);
videoname = strjoin({Nstring,'_burgers_LW.avi'});

set(gcf,'Renderer','zbuffer');
F(Nt) = struct('cdata',[],'colormap',[]);

Burgers = VideoWriter(videoname,'Uncompressed AVI');
open(Burgers)
 for i = 1:Nt-1
     [X,Y] = meshgrid(xgrid,xgrid);
     drawnow
     xlabel('x'), ylabel('y'), zlabel('u');
     surf(X,Y,UU(:,:,i)');
     view(60,30);
     title(['Burgers Equation in 2D -- Lax-Wendroff scheme \newline t = ' num2str(i*dt) ', N = ' num2str(N)]);
     F(i) = getframe(gcf);
     writeVideo(Burgers,F(i).cdata)
 end
close(Burgers);

%% Mimetic Method in Space; Forward in Time.
% This method combines both the space approximation via the MOLE Mimetic
% library and the forward time approximation.

addpath('../MATLAB/mole/mole_MATLAB')

% 2D Staggered grid and the length.
dx = (xMax - xMin)/N;
dt = 0.2;
xgrid = [xMin xMin+dx/2 : dx : xMax-dx/2 xMax];

D = div(2, N, dx); % 1D Mimetic divergence operator
I = interpol(N, 0.5); % 1D 2nd order interpolator

% Impose IC
U = 1/(sigma*sqrt(2*pi))*exp(-0.5/sigma^2*((xgrid-x0).^2))';

% Premultiply out of the time loop (since it doesn't change)
D = -dt*(alpha+beta)*D*I;

% Time integration loop
for i = 0:Ntt/dt
    
    % Area conservation per step.
    A(i+1) = trapz(U); 
    UUx(:,i+1) = U;
    U2 = U + U.*D*U; 
    U = U2;
end 

 ConserveM = std(A)

Nstring = num2str(N);
videoname = strjoin({Nstring,'_burgers_mimetic.avi'});

Burgers = VideoWriter(videoname,'Uncompressed AVI');
Burgers.FrameRate=15;
open(Burgers)
 for i = 0:Ntt/dt
    UU = UUx(:,i+1)*UUx(:,i+1)';
    [X,Y] = meshgrid(xgrid,xgrid);
    drawnow
    xlabel('x'), ylabel('y'), zlabel('u');
    surf(X,Y,UU);
    view(60,30);
    title(['Burgers Equation in 2D -- Forward Time and Mimetic scheme \newline t = ' num2str(i+1) ', N = ' num2str(N)]);
    F(i+1) = getframe(gcf);
    writeVideo(Burgers,F(i+1).cdata)
 end
close(Burgers)
