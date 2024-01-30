clear
%clf;

%Polar grid

%Longitude
phinum=101;
phimin=0;
phimax=2*pi;
W=phimax-phimin;
dphi=W/(phinum-1);
phi=0:dphi:2*pi;

%Radial distance
rnum=51;
rmin=3504e+3;
rmax=6371e+3;
H=rmax-rmin;
dR=H/(rnum-1);
r=rmin:dR:rmax;

%Define constants
kphi=5; %long. wave number

kr=1; %radial wave number
vr0=1e-9; %reference radial velocity

T0=298; %surface T
T1 = 4000; %bottom T
dTmax = 1000; %max. amplitude of T perturbation

P0=0.1e+6; %surface P
P1=135e+9; %bottom P

%Define matrices and vectors
x2d=zeros(phinum,rnum); y2d=x2d; Tk2d=x2d; Pa2d=x2d; vr2d=x2d; vphi2d=x2d; div2d=x2d;
vphi=zeros(phinum*rnum,1); vr=vphi; Tk=vr; Pa=vr;

%Compute V,P,T
for i1=1:phinum
    for i2=1:rnum
        
        gi = i2 + (i1-1)*rnum;
        
        %Normalized radial distance
        R = (r(i2)-rmin)/H;
        
        sinp = sin(kphi*2*pi*phi(i1)/W); cosp =cos(kphi*2*pi*phi(i1)/W);
        sinr = sin(kr*pi*R); cosr = cos(kr*pi*R);
        
        vphi0 = vr0*2/kphi*(sinr/cosr + (r(i2))/H*kr*pi);
        vphi(gi)= -vphi0*sinp*cosr;
        vr(gi)  =    vr0*cosp*sinr;
        
        %Define P-T conditions
        dT = 0;
        if vr0 ~= 0
            dT = dTmax*vr(gi)/vr0; %perturbation dependent on radial velocity
        end
        Tk(gi) = T0 + (rmax-r(i2))/H*T1 + dT;
        Pa(gi) = P0 + (rmax-r(i2))/H*P1;
        
        %2D matrices
        vr2d(i1,i2) = vr(gi);
        vphi2d(i1,i2) = vphi(gi);
        div2d(i1,i2) = vr(gi) + (r(i2))/H*vr0*cosp*cosr*kr*pi - vphi0*cosp*kphi/2*cosr;
        
        Tk2d(i1,i2)=Tk(gi);
        Pa2d(i1,i2)=Pa(gi);
        
        %Cartesian coordinates in km
        x2d(i1,i2)=cos(phi(i1))*r(i2)/1000;
        y2d(i1,i2)=sin(phi(i1))*r(i2)/1000;
        
    end
end

figure(10)
subplot(2,3,1)
pcolor(x2d,y2d,(vr2d.^2 + vphi2d.^2).^0.5/vr0)
shading interp
colorbar('southoutside')
hold on
contour(x2d,y2d,(vr2d.^2 + vphi2d.^2).^0.5/vr0,20,'k','LineWidth',0.5);
hold off
axis image square
title('Scaled velocity magnitude')

subplot(2,3,2)
pcolor(x2d,y2d,vr2d./vr0)
shading interp
colorbar('southoutside')
hold on
contour(x2d,y2d,vr2d./vr0,20,'k','LineWidth',0.5);
hold off
axis image square
title('Scaled radial Velocity')

subplot(2,3,3)
pcolor(x2d,y2d,vphi2d./vr0)
shading interp
colorbar('southoutside')
hold on
contour(x2d,y2d,vphi2d./vr0,20,'k','LineWidth',0.5);
hold off
axis image square
title('Scaled tangential Velocity')

subplot(234)
pcolor(x2d,y2d,div2d)
shading interp
colorbar('southoutside')
axis image square
title('Divergence')

subplot(235)
pcolor(x2d,y2d,Tk2d)
shading interp
colorbar('southoutside')
hold on
[c,h]=contour(x2d,y2d,(Tk2d-T0)./(T1-T0),20,'k','LineWidth',0.5);
hold off
axis image square
title('Temperature')

subplot(236)
pcolor(x2d,y2d,(Pa2d-P0)./(P1-P0))
shading interp
colorbar('southoutside')
axis image square
title('Pressure')


%Save hdf5 file
fname='vtp0001.h5';
dt=1e+5*365.25*86400;
t=[dt 0];
hdf5write(fname,'/Nodes/V1',vphi);
hdf5write(fname,'/Nodes/V2',vr, 'WriteMode', 'append');
hdf5write(fname,'/Nodes/Tk',Tk, 'WriteMode', 'append');
hdf5write(fname,'/Nodes/P',Pa, 'WriteMode', 'append');
h5writeatt(fname,'/','Time',t);
