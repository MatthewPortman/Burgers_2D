%% Matthew Portman - Due December 19
%
%
%% Solving Burger's Equation in 2D
% In this code, we use a second order central difference approximation to
% solve the time-independent Poisson Equation in 2D.

clear; close all; 
%% Definitions
% First, we define an index to loop over as well as the number of steps to
% be used in the approximation. The loop is explicitly used for error
% calculation but is otherwise extraneous. 
%set(0,'DefaultFigureVisible','off');


k = 2;
%Nk = 100;
N = 50;
Nt = 50;
Ntt = 150;
%dt = 0.1;                         %Width of each time step

%while k >= 1 
%if k == 2
%    N = 0.5*Nk;
%else
%    N = Nk;
%end

%%
% After choosing the number of steps, we proceed. Defining some frequently 
% used terms based on N, the spatial bounds, the step width h, the grid on 
% which we are calculating the approximation and finally the equation f and
% analytical solution a_s.

N1 = N-1;
N11 = N+2;

xMax = 10;
xMin = -10;

dx = (xMax - xMin)/N1;
dt = dx;
dx1 = 1/dx;
dxdx = 1/(dx*dx);

xgrid = xMin : dx: xMax;
tgrid = 0 : dt : Nt*dt;
theta = pi/6;
sigma = 2;

alpha = cos(theta);
beta = sin(theta);

% Initial Velocity distribution/BC
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

%fig1 = figure('Renderer','zbuffer');
set(gcf,'Renderer','zbuffer');
set(gca,'NextPlot','replaceChildren');
F(Nt) = struct('cdata',[],'colormap',[]);
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

%Burgers = VideoWriter('burgers_LW.avi','Uncompressed AVI');
%open(Burgers)
% for i = 1:Nt-1
%     [X,Y] = meshgrid(xgrid,xgrid);
%     %hold on
%     drawnow
%    % figure(fig1);
%     title(['Burgers Equation in 2D -- Lax-Wendroff scheme \newline t = ' num2str(i) ', N = ' num2str(N)]);
%     xlabel('x'), ylabel('y'), zlabel('u');
%     surf(X,Y,UU(:,:,i)');
%     view(60,30);
%     F(i) = getframe(gcf);
%     writeVideo(Burgers,F(i).cdata)
% end
%close(Burgers);
% movie(F,1,30)

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
    avgFTCS(n) = norm(trapz(UFTCS(:,:,n)));
end

ConserveFTCS = std(avgFTCS)
              
%% Applying the Boundary Conditions and the Numerical Solution
% Having generated the basic L and F matrices we must now apply the
% boundary condition that u = 0 on the boundary of the surface. To do so,
% we know that F must be 0 around the boundary and that L*U should give
% those 0's. Therefore, L is calculated by assessing where F is 0 and
% applying that condition to the matrix appropriately.



%% Mimetic Method in Space; Leapfrog in Time.
% This method combines both the space approximation via the MOLE Mimetic
% library and the leapfrog time approximation to find a solution matrix U
% which is of a different form from the others.

addpath('../mole/mole_MATLAB')

% 2D Staggered grid and the length.
dx = (xMax - xMin)/N;
dt = dx;
xgrid = [xMin xMin+dx/2 : dx : xMax-dx/2 xMax];

ygrid = xgrid;
 
 II = interpol2D(N,N,1,1);
% 
% % We use the mimetic Divergence operator to discretize the divergence of
% % of U_x. We then follow the Leap Frog scheme in time as in the previous
% % section and generate a solution.
 
 D = div2D(2, N, dx, N, dx);
 DD = -dt*(alpha+beta)*D*II;
% 
U1 = zeros(N11,N11);
% 
 for j = 1:N11
     for i = 1:N11
         U1(i,j) = IC(xgrid(i),xgrid(j)); 
     end
 end

%DD(1,2) = 1;
%DD(end,end-1) = -1;
%for j =2:N+1
%  DD(j,j+1) = 1;
%  DD(j,j-1) = -1;
%end
 
U1(1,:) = 0;
U1(:,1) = 0;
U1(end,:) = 0;
U1(:,end) = 0;


% % Setting the first time step.
 UU = zeros(N11,N11,Nt);
 
 UU(:,:,1) = U1;
 UU = reshape(UU,1,N11*N11,Nt);
 
 const = (alpha+beta);%/(alpha*beta);
  for n = 1 : Nt-1
     
      UU(1,:,n+1) = UU(1,:,n)+UU(1,:,n).*(const*UU(1,:,n)*DD);
     
  end
  
  U = reshape(UU,N11,N11,Nt); 
  %for n = 1:Nt-1
   % avgM(n) = norm(trapz(UU(:,:,n)));
  %end

 ConserveM = std(avgM)

Burgers = VideoWriter('burgers_mimetic1.avi','Uncompressed AVI');
open(Burgers)
 for i = 1:Nt-1
    [X,Y] = meshgrid(xgrid,xgrid);
    %hold on
    drawnow
    %figure(fig1);
    title(['Burgers Equation in 2D -- Leapfrog and Mimetic scheme \newline t = ' num2str(i*dt) ', N = ' num2str(N)]);
    xlabel('x'), ylabel('y'), zlabel('u');
    surf(X,Y,UU(:,:,i));
    view(60,30);
    F(i) = getframe(gcf);
    writeVideo(Burgers,F(i).cdata)
end

close(Burgers)

%% Error Calculations
% The first error is calculated inside the loop to store it while the
% second is calculated outside of it. The values of the resultant matrices
% are averaged to show that the error is 0.

% if k == 2
%     err1 = U_as(2:end-1,2:end-1) - Uplot(2:end-1,2:end-1);
% end
% 
% k = k-1;
% %end
% 
% err2 = U_as(2:end-1,2:end-1) - Uplot(2:end-1,2:end-1);
%     
% avg1 = mean(mean(err1))
% avg2 = mean(mean(err2))
% 
% %% Plotting the solutions
% % Two plots are generated separately for the numerical and analytical 
% % solutions since an overlay would be impossible to distinguish. 
% 
% [X,Y] = meshgrid(xGrid,xGrid);
% 
% surf(X,Y,Uplot)
% 
% grid on
% xlabel('x')
% ylabel('y')
% zlabel('u(x,y)')
% title('2D Poisson Problem: Numerical Solution')
% 
% figure()
% 
% surf(X,Y,U_as)
% grid on
% xlabel('x')
% ylabel('y')
% zlabel('u(x,y)')
% title('2D Poisson Problem: Analytical Solution')

%% Closing Remarks
% The solution was remarkably easy to calculate once the L matrix was
% determined. I would imagine that this technique is more efficient
% than a loop up to a certain size N after which the MatLab operator
% "\" would no longer be useful. I would be interested to see how this
% calculation is performed using mimetic methods after discussing the
% solution briefly in class.
