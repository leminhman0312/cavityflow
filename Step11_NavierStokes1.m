%Step 10 :  Poisson Equation
 
% 
 
clear all
nx=5; ny=5; nt=100; nit=50;
xmin=0; xmax=2; 
ymin=0; ymax=2;
dx = (xmax-xmin)/(nx-1); 
dy=(ymax-ymin)/(ny-1);
x=linspace(0,2,nx);
y=linspace(0,2,ny);
[Y,X]=meshgrid(y,x);
 
rho=1;
nu=0.1;
dt=0.001;
 
 
% Init
u=zeros(ny,nx);
v=zeros(ny,nx);
p=zeros(ny,nx);
b=zeros(ny,nx);
%Pressure Field
%Square Brackets of Poissons Equation
   
 
 
for it = 1:nt+1
        for i=2:(nx-1)
            for j=2:(ny-1)
                b(i,j)=rho*(1/dt*((u(i+1,j)-u(i-1,j))/(2*dx)+(v(i,j+1)-v(i,j-1))/(2*dy))-((u(i+1,j)-u(i-1,j))/(2*dx))^2-2*((u(i,j+1)-u(i,j-1))/(2*dy)*(v(i+1,j)-v(i-1,j))/(2*dx))-((v(i,j+1)-v(i,j-1))/(2*dy))^2);
            end
        end
    for iit=1:nit+1
        pn=p;
        for i=2:(nx-1)
            for j=2:(ny-1) 
            p(i,j)=((pn(i+1,j)+pn(i-1,j))*dy^2+(pn(i,j+1)+pn(i,j-1))*dx^2)/(2*(dx^2+dy^2))-dx^2*dy^2/(2*(dx^2+dy^2))*b(i,j);
            end
        end
        p(:,ny) =p(:,ny-1); %%dp/dy = 0 at y = 2
        p(1,:) = p(2,:);    %%dp/dy = 0 at y = 0
        p(:,1)=p(:,2);      %%dp/dx = 0 at x = 0
        p(:,ny)=0;          %%dp = 0 at y = 2
 
    end
     
    un = u;
    vn = v;
     
    for j=2:nx-1
        for i=2:ny-1        
        %Velocity Field
        u(i,j) = un(i,j)-un(i,j)*dt/dx*(un(i,j)-un(i-1,j))-vn(i,j)*dt/dy*(un(i,j)-un(i,j-1))-dt/(2*rho*dx)*(p(i+1,j)-p(i-1,j))+nu*(dt/dx^2*(un(i+1,j)-2*un(i,j)+un(i-1,j))+ (dt/dy^2*(un(i,j+1)-2*un(i,j)+un(i,j-1))));
        v(i,j) = vn(i,j)-un(i,j)*dt/dx*(vn(i,j)-vn(i-1,j))-vn(i,j)*dt/dy*(vn(i,j)-vn(i,j-1))-dt/(2*rho*dy)*(p(i,j+1)-p(i,j-1))+nu*(dt/dx^2*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))+ (dt/dy^2*(vn(i,j+1)-2*vn(i,j)+vn(i,j-1))));
        end
    end
        u(1,:)=0;
        u(:,1)=0;
        u(nx,:)=0;
        u(:,ny)=1;
        v(1,:)=0;
        v(ny,:)=0;
        v(:,1)=0;
        v(nx,:)=0;
 
 end
contourf(x,y,p.','linewidth',1)
hold on
quiver(x,y,u.',v.',2)
xlabel('x')
ylabel('y')
colorbar