clear
%clf;

%Polar grid

%Longitude
xnum=101;
xmin=0;
xmax=600e+3;
W=xmax-xmin;
dx=W/(xnum-1);
x=xmin:dx:xmax;

%Radial distance
ynum=51;
ymin=0e+3;
ymax=400e+3;
H=ymax-ymin;
dy=H/(ynum-1);
y=ymin:dy:ymax;

%Define constants
kx=1; %horiz. wave number
vx0=1e-9; %reference horiz. velocity

ky=1; %vertical wave number
vy0=vx0*2*kx/ky*H/W; %reference vertical velocity

%Define matrices and vectors
x2d=zeros(xnum,ynum); y2d=x2d; vy2d=x2d; vx2d=x2d; div2d=x2d;
vx=zeros(xnum*ynum,1); vy=vx;

%Compute V,P,T
for i1=1:xnum
    for i2=1:ynum

        
        gi = i2 + (i1-1)*ynum;

        sinx = sin(kx*2*pi*x(i1)/W); siny = sin(ky*pi*y(i2)/H); 
        cosx = cos(kx*2*pi*x(i1)/W); cosy = cos(ky*pi*y(i2)/H);

        vx(gi)  = -vx0*sinx*cosy;
        vy(gi)  =  vy0*cosx*siny;

        %2D matrices
        vx2d(i1,i2) = vx(gi);
        vy2d(i1,i2) = vy(gi);
        div2d(i1,i2) = -vx0*cosx*cosy*kx*2*pi/W + vy0*cosx*cosy*ky*pi/H;

        %Cartesian coordinates in km
        x2d(i1,i2)=x(i1)/1000;
        y2d(i1,i2)=y(i2)/1000;

    end
end

figure(10)
subplot(221)
pcolor(x2d,y2d,(vx2d.^2 + vy2d.^2).^0.5/vy0)
shading interp
colorbar('southoutside')
hold on
contour(x2d,y2d,(vx2d.^2 + vy2d.^2).^0.5/vy0,20,'k','LineWidth',0.5);
hold off
axis image
title('Scaled velocity magnitude')

subplot(222)
pcolor(x2d,y2d,vy2d./vy0)
shading interp
colorbar('southoutside')
hold on
contour(x2d,y2d,vy2d./vy0,20,'k','LineWidth',0.5);
hold off
axis image
title('Scaled radial Velocity')

subplot(223)
pcolor(x2d,y2d,vx2d./vy0)
shading interp
colorbar('southoutside')
hold on
contour(x2d,y2d,vx2d./vy0,20,'k','LineWidth',0.5);
hold off
axis image
title('Scaled tangential Velocity')

subplot(224)
pcolor(x2d,y2d,div2d)
shading interp
colorbar('southoutside')
axis image
title('Divergence')

%Save hdf5 file
fname='vtp0001.h5';
dt=1e+5*365.25*86400;
t=[dt 0];
hdf5write(fname,'/Nodes/V1',vx);
hdf5write(fname,'/Nodes/V2',vy, 'WriteMode', 'append');
%hdf5write(fname,'/Nodes/Tk',Tk, 'WriteMode', 'append');
%hdf5write(fname,'/Nodes/P',Pa, 'WriteMode', 'append');
h5writeatt(fname,'/','Time',t);