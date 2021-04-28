%Make vtp*.h5 files for gloabl spherical domain with YIN YANG grids

clear

%YIN YANG grids
%Axis 1: longitude
phinum= 101;
dphi=270/(phinum-1-2);
phi=-135-dphi:dphi:135+dphi;

%Axis 2: radial distance
rnum= 21;
dr=1/(rnum-1-2);
r=1.2222222:dr:2.2222222;

%Axis 3: colatitude
colatnum= 31;
dcolat=90/(colatnum-1-2);
colat=45-dcolat:dcolat:135+dcolat;

%Convert long. and colat. in radians
deg2rad=pi/180;
phi=phi.*deg2rad;
colat=colat.*deg2rad;

%Plot YIN YANG grids
plotyinyang(phi,r(10),colat,phinum,colatnum)

%Number of nodes
nx1= phinum  ;%number of nodes along axis 1
nx2= rnum    ;%number of nodes along axis 2
nx3= colatnum;%number of nodes along axis 3
nodenum=nx1*nx2*nx3;


%File infos
finit = 1; %Initial file number
fstp  = 1; %Step file number
fend  = 1; %Final file number
for f = finit:fstp:fend
    
    %A) Load here the Eulerian fields of your geodynamic model
    V1model=...; velocity along the axis 1
    V2model=...; velocity along the axis 2
    V3model=...; velocity along the axis 3
    Tkmodel=...; Temperature
    Pmodel=... ; Pressure
    Fdmodel=...; Fraction of deformation accommodated by disl. creep (adimensional)
        
    time=...; %total elapsed time
    dt=...  ; %elapsed time from previous output file
        
    %B) Convert your arrays into 1D D-REX_M format
    V1=zeros(nodenum*2,1); V2=V1; V3=V1; Tk=V1; P=V1; Fd=V1;
    for yy = 1:2
        for i1 = 1:nx1
            for i2 = 1:nx2
                for i3 = 1:nx3
                    
                    %Global index for D-REX_M 1D array
                    i4 = nodenum*(yy-1) + i2 + (i1-1)*nx2 + (i3-1)*nx1*nx2;
                    
                    if (yy==1)
                        %Coord. of YIN grid in the spherical grid
                        pp = phi(i1);
                        cc = colat(i3);
                    else
                        %Coord. of YANG grid in the spherical grid
                        [pp,cc]=yin2yang(phi(i1),colat(i3));
                    end
                    rr=r(i2);
                    
                    %Interpolate here the eulerian fields of the geodynamic model
                    % to the long. (pp), radial (rr) and colat. (cc) coordinates.
                    V1model_interpolated = ...;
                    
                    V1(i4) = V1model_interpolated; % in m/s or adimensional
                    V2(i4) = V2model_interpolated; % in m/s or adimensional
                    V3(i4) = V3model_interpolated; % in m/s or adimensional
                    Tk(i4) = Tkmodel_interpolated; % in Kelvin
                    P(i4)  =  Pmodel_interpolated; % in Pa
                    Fd(i4) = Fdmodel_interpolated; % adimensional, 0?Fd?1
                    
                end
            end
        end
    end
    
    t=[dt time]; % array with time infos (in sec or adimensional)
    
    %C)Save arrays in HDF5 format files
    fname=['vtp', num2str(f,'%.4d'),'.h5'];
    hdf5write(fname,'/Nodes/V1',V1);
    hdf5write(fname,'/Nodes/V2',V2, 'WriteMode', 'append');
    hdf5write(fname,'/Nodes/V3',V3, 'WriteMode', 'append');
    hdf5write(fname,'/Nodes/Tk',Tk, 'WriteMode', 'append');
    hdf5write(fname,'/Nodes/P',P, 'WriteMode', 'append');
    hdf5write(fname,'/Nodes/Fd',Fd, 'WriteMode', 'append');
    h5writeatt(fname,'/','Time',t);
    
end



function plotyinyang(phi,R,colat,phinum,colatnum)

%YIN grid at given radius
X=zeros(phinum,colatnum); Y=X; Z=X;
C=ones(phinum,colatnum);
%YIN grid
for p=1:phinum
    for c=1:colatnum
        X(p,c)=R*sin(colat(c))*cos(phi(p));
        Y(p,c)=R*sin(colat(c))*sin(phi(p));
        Z(p,c)=R*cos(colat(c));
    end
end

%YANG grid at given radius
X1=zeros(phinum,colatnum); Y1=X1; Z1=X1;
for p=1:phinum
    for c=1:colatnum
        
        %Convert YANG grid coordinates into YIN coord. system
        [pp,cc]=yin2yang(phi(p),colat(c));
        
        X1(p,c)=R*sin(cc)*cos(pp);
        Y1(p,c)=R*sin(cc)*sin(pp);
        Z1(p,c)=R*cos(cc);
    end
end

%Plot YIN YANG grids
figure(1)
surf(X,Y,Z,C)
hold on
surf(X1,Y1,Z1,C)
hold off
end


function [pp,cc]=yin2yang(p,c)
sinl=sin(p);
cosl=cos(p);
sinc=sin(c);
cosc=cos(c);
cosmc=sinl*sinc;
if(cosmc<-1.0); cosmc=-1.0; end
if(cosmc> 1.0); cosmc= 1.0; end

%New long. and colat.
pp=atan2(cosc,-sinc*cosl);
cc=acos(cosmc);
end