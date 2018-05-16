% Marker-and-Cell (MAC) Finite Volume Solver for the steady incompressible
% Stokes equations in 2 dimensions. Domain is a
% rectangle, with zero velocity BC
%
% Author: Gustavo Buscaglia
% Last Modification: 12/6/2017

%Problem definition
global N1 N2
Lx = 30.0;
Ly = 30.0;
N1 = 100; %number of cells in the domain along x
N2 = 100;  %number of cells in the domain along y
%Material properties
mu = 1/30;
%Default mesh
aux = Lx/N1;
Nx = N1; %total cells along x
X = [0:aux:Lx]; 
aux = Ly/N2;
Ny = N2; %total cells along y
Y = [0:aux:Ly];
%Cell sizes
dx = X(2:Nx+1)-X(1:Nx);
dy = Y(2:Ny+1)-Y(1:Ny);
%Cell center coordinates
Xh = 0.5*(X(2:Nx+1)+X(1:Nx));
Yh = 0.5*(Y(2:Ny+1)+Y(1:Ny));
dxh = Xh(2:Nx)-Xh(1:Nx-1);
dyh = Yh(2:Ny)-Yh(1:Ny-1);
%Arrays
dimU = (Nx+1)*Ny;
dimV = Nx*(Ny+1);
dimP = Nx*Ny;
%Global matrix and rhs
nunk = dimU + dimV + dimP;
Ag = spalloc(nunk,nunk,27*nunk);
rhs = zeros(nunk,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%X-momentum equation ((Nx+1)*Ny control volumes)
%%The leftmost and rightmost unknowns are always 0
ij2ng =@(N,i,j) i + (j-1)*N;


for j=1:Ny
   nguP = ij2ng(Nx+1,1,j);
   Ag(nguP,nguP)=1; rhs(nguP)=0;
   nguP=ij2ng(Nx+1,Nx+1,j);
   Ag(nguP,nguP)=1; rhs(nguP)=0;
end

%%Now the inner rows and columns
for i=2:Nx
 for j=1:Ny
  nguP=ij2ng(Nx+1,i,j);
  xP=X(i);yP=Yh(j);
%%  rhs(nguP)     =rhs(nguP)     + force(xP,yP,1)*dxh(i-1)*dy(j);
  if ((i==floor(Nx/2))&&(j==floor(Ny/2))) 
   rhs(nguP)     =rhs(nguP)     + 1;
  end
  %%% W boundary (integral of -p+2*muW*du/dx)
    muW=mu;
    nguWW=ij2ng(Nx+1,i-1,j); ngpW=ij2ng(Nx,i-1,j)+dimU+dimV;
    Ag(nguP,nguP)=Ag(nguP,nguP)+2*muW*dy(j)/dx(i-1);
    Ag(nguP,nguWW)=Ag(nguP,nguWW)-2*muW*dy(j)/dx(i-1);
    Ag(nguP,ngpW)=-dy(j);
  %%% E boundary (integral of p-2*muE*du/dx)
    muE=mu;
    nguEE=ij2ng(Nx+1,i+1,j); ngpE=ij2ng(Nx,i,j)+dimU+dimV;
    Ag(nguP,nguP)=Ag(nguP,nguP)+2*muE*dy(j)/dx(i);
    Ag(nguP,nguEE)=Ag(nguP,nguEE)-2*muE*dy(j)/dx(i);
    Ag(nguP,ngpE)=dy(j);
  %%% N boundary (integral of -mu(dv/dx))
    muN=mu;
    nguNN=ij2ng(Nx+1,i,j+1);ngvNW=ij2ng(Nx,i-1,j+1)+dimU;ngvNE=ij2ng(Nx,i,j+1)+dimU;
    Ag(nguP,ngvNE)=Ag(nguP,ngvNE)-muN*dxh(i-1)/dxh(i-1);
    Ag(nguP,ngvNW)=Ag(nguP,ngvNW)+muN*dxh(i-1)/dxh(i-1);
  %%% N boundary (integral of -mu(du/dy))
    if (j ~= Ny) %% north cell boundary fully fluid
     Ag(nguP,nguP) =Ag(nguP,nguP) +muN*dxh(i-1)/dyh(j);
     Ag(nguP,nguNN)=Ag(nguP,nguNN)-muN*dxh(i-1)/dyh(j);
    else %%N cells both boundary
     Ag(nguP,nguP) =Ag(nguP,nguP) +muN*dxh(i-1)/(dy(j)/2);
%%     rhs(nguP)     =rhs(nguP)     +...
%%      muN*dxh(i-1)/(dy(j)/2)*(0+0)/2; zero u at the wall assumed
    end 
    %%% S boundary (integral of mu(dv/dx))
    muS=mu;
    nguSS=ij2ng(Nx+1,i,j-1);ngvSW=ij2ng(Nx,i-1,j)+dimU;ngvSE=ij2ng(Nx,i,j)+dimU;
    Ag(nguP,ngvSE)=Ag(nguP,ngvSE)+muS*dxh(i-1)/dxh(i-1);
    Ag(nguP,ngvSW)=Ag(nguP,ngvSW)-muS*dxh(i-1)/dxh(i-1);
    %%% S boundary (integral of mu(du/dy))
    if (j~=1) %% south cell boundary fully fluid
     Ag(nguP,nguP) =Ag(nguP,nguP) +muS*dxh(i-1)/dyh(j-1);
     Ag(nguP,nguSS)=Ag(nguP,nguSS)-muS*dxh(i-1)/dyh(j-1);
    else  %%S cells both boundary
     Ag(nguP,nguP) =Ag(nguP,nguP) +muS*dxh(i-1)/(dy(j)/2);
%     rhs(nguP)     =rhs(nguP)     +...
%      muS*dxh(i-1)/(dy(j)/2)*(0+0)/2; zero u at the wall assumed
    end 
 end
end
%%%%%%%%%%%%%%%End of X-Momentum equations
%%%%%%%%%%%%%%%Now Y-Momentum
%%%Y-momentum equation (Nx*(Ny+1) control volumes)
%%The northmost and southmost unknowns are always zero
for i=1:Nx
   ngvP=ij2ng(Nx,i,1)+dimU;
   Ag(ngvP,ngvP)=1; rhs(ngvP)=0;
   ngvP=ij2ng(Nx,i,Ny+1)+dimU;
   Ag(ngvP,ngvP)=1; rhs(ngvP)=0;
end
%%Now the inner rows and columns
for i=1:Nx
 for j=2:Ny
  ngvP=ij2ng(Nx,i,j)+dimU;
  xP=Xh(i);yP=Y(j);
%%  rhs(ngvP)     =rhs(ngvP)     + force(xP,yP,2)*dx(i)*dyh(j-1);
  if (i==floor(Nx/2)&&j==floor(Ny/2))
   rhs(ngvP)     =rhs(ngvP)     + 0;
  end
  %%% S boundary (integral of -p+2*muS*dv/dy)
    muS=mu;
    ngvSS=ij2ng(Nx,i,j-1)+dimU; ngpS=ij2ng(Nx,i,j-1)+dimU+dimV;
    Ag(ngvP,ngvP)=Ag(ngvP,ngvP)+2*muS*dx(i)/dy(j-1);
    Ag(ngvP,ngvSS)=Ag(ngvP,ngvSS)-2*muS*dx(i)/dy(j-1);
    Ag(ngvP,ngpS)=-dx(i);
  %%% N boundary (integral of p-2*muN*dv/dy)
    muN=mu;
    ngvNN=ij2ng(Nx,i,j+1)+dimU; ngpN=ij2ng(Nx,i,j)+dimU+dimV;
    Ag(ngvP,ngvP)=Ag(ngvP,ngvP)+2*muN*dx(i)/dy(j);
    Ag(ngvP,ngvNN)=Ag(ngvP,ngvNN)-2*muN*dx(i)/dy(j);
    Ag(ngvP,ngpN)=dx(i);
  %%% E boundary (integral of -mu(du/dy))
    muE=mu;
    ngvEE=ij2ng(Nx,i+1,j)+dimU;nguNE=ij2ng(Nx+1,i+1,j);nguSE=ij2ng(Nx+1,i+1,j-1);
    Ag(ngvP,nguNE)=Ag(ngvP,nguNE)-muE*dyh(j-1)/dyh(j-1);
    Ag(ngvP,nguSE)=Ag(ngvP,ngvSE)+muE*dyh(j-1)/dyh(j-1);
  %%% E boundary (integral of -mu(dv/dx))
    if (i~=Nx) %% east cell boundary fully fluid
     Ag(ngvP,ngvP) =Ag(ngvP,ngvP) +muE*dyh(j-1)/dxh(i);
     Ag(ngvP,ngvEE)=Ag(ngvP,ngvEE)-muE*dyh(j-1)/dxh(i);
    else %% E cells both boundary
     Ag(ngvP,ngvP) =Ag(ngvP,ngvP) +muE*dyh(j-1)/(dx(i)/2);
  %%   rhs(ngvP)     =rhs(ngvP)     +...
  %%    muE*dyh(j-1)/(dx(i)/2)*(0+0)/2;
    end 
    %%% W boundary (integral of mu(du/dy))
    muW=mu;
    ngvWW=ij2ng(Nx,i-1,j)+dimU;nguNW=ij2ng(Nx+1,i,j);nguSW=ij2ng(Nx+1,i,j-1);
    Ag(ngvP,nguNW)=Ag(ngvP,nguNW)+muW*dyh(j-1)/dyh(j-1);
    Ag(ngvP,nguSW)=Ag(ngvP,nguSW)-muW*dyh(j-1)/dyh(j-1);
    %%% W boundary (integral of mu(dv/dx))
    if (i~=1) %% W cell boundary fully fluid
     Ag(ngvP,ngvP) =Ag(ngvP,ngvP) +muW*dyh(j-1)/dxh(i-1);
     Ag(ngvP,ngvWW)=Ag(ngvP,ngvWW)-muW*dyh(j-1)/dxh(i-1);
    else %% W cells both boundary
     Ag(ngvP,ngvP) =Ag(ngvP,ngvP) +muW*dyh(j-1)/(dx(i)/2);
%%     rhs(ngvP)     =rhs(ngvP)     +...
%%      muW*dyh(j-1)/(dx(i)/2)*(0+0)/2;
    end 
 end
end
%%%%%%%%%%%%%%%End of Y-Momentum equations
%%%%%%%%%%%%%%%Beginning of Mass equations
%%%Mass equation (Nx*Ny control volumes)
for i=1:Nx
 for j=1:Ny
  ngpP = ij2ng(Nx,i,j)+dimU+dimV;
  xP = Xh(i);yP=Yh(j);
  if (i+j == 2) %% cell with imposed pressure
   Ag(ngpP,ngpP)=1; rhs(ngpP)=0;
  else    %% pressure unknown is alive (fluid cell)
   %%% S boundary (integral of -v)
   ngvS=ij2ng(Nx,i,j)+dimU;
   Ag(ngpP,ngvS)=Ag(ngpP,ngvS)-dx(i);
   %%% N boundary (integral of v)
   ngvN=ij2ng(Nx,i,j+1)+dimU;
   Ag(ngpP,ngvN)=Ag(ngpP,ngvN)+dx(i);
   %%% E boundary (integral of u)
   nguE=ij2ng(Nx+1,i+1,j);
   Ag(ngpP,nguE)=Ag(ngpP,nguE)+dy(j);
   %%% W boundary (integral of -u)
   nguW=ij2ng(Nx+1,i,j);
   Ag(ngpP,nguW)=Ag(ngpP,nguW)-dy(j);
  end
 end
end
%%%%%%%%%%%%%%%End of Mass equations
xxx = Ag \ rhs;
%%%%%%%%%%Cell-average velocities
Uav=zeros(Nx,Ny);
Vav=zeros(Nx,Ny);
P=zeros(Nx,Ny);
U=zeros(Nx+1,Ny);
V=zeros(Nx,Ny+1);
for i=1:Nx+1
 for j=1:Ny
  U(i,j)=xxx(ij2ng(Nx+1,i,j));
 end
end
for i=1:Nx
 for j=1:Ny+1
  V(i,j)=xxx(ij2ng(Nx,i,j)+dimU);
 end
end
for i=1:Nx
 for j=1:Ny
   Uav(i,j)=0.5*(xxx(ij2ng(Nx+1,i,j))+xxx(ij2ng(Nx+1,i+1,j)));
   Vav(i,j)=0.5*(xxx(ij2ng(Nx,i,j)+dimU)+xxx(ij2ng(Nx,i,j+1)+dimU));
   P(i,j)=xxx(ij2ng(Nx,i,j)+dimU+dimV);
 end
end

%%%%%%%%%%%%%%%
figure(1);
h=quiver(Xh,Yh,Uav',Vav',3);
set(h,'linewidth',2);

figure(2);
contourf(X,Yh,U');

figure(3);
contourf(Xh,Y,V');

%%contourf(Xh(2:Nx-1),Y(2:Ny),V(2:Nx-1,2:Ny)')   %%only active ones
figure(4);
contourf(Xh(2:Nx-1),Yh(2:Ny-1),P(2:Nx-1,2:Ny-1)');
