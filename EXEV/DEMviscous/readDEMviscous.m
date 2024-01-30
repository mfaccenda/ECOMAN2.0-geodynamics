%Read DEM
read=1;
if read
    clear all
    read=1;
end

close all;

path = './';

logscale=1;
printmod=0;

if read

    fname = [path,'viscoustensordem.h5'];

    nx=h5read(fname,'/gridnum');
    etacontr=h5read(fname,'/etacontr');
    Vol=h5read(fname,'/Vol');
    a1a2=h5read(fname,'/esa12');
    a2a3=h5read(fname,'/esa23');
    C11=h5read(fname,'/C11');
    C22=h5read(fname,'/C22');
    C33=h5read(fname,'/C33');
    C44=h5read(fname,'/C44');
    C55=h5read(fname,'/C55');
    C66=h5read(fname,'/C66');
    C12=h5read(fname,'/C12');
    C13=h5read(fname,'/C13');
    C23=h5read(fname,'/C23');
    
    gridnum = h5read(fname,'/gridnum');

    etacontrnum = gridnum(1);
    volnum = gridnum(2);
    a1a2num = gridnum(3);
    a2a3num = gridnum(4);

    %Strain rate components
    eps_zz=1; eps_xx=-0.5; eps_yy=-0.5; eps_xz=1;

    %Deviatoric normal stress during uniaxial compression parallel to z axis
    sign=C33.*eps_zz + C13*eps_xx + C23*eps_yy;

    %Average shear stress in the vertical plane
    sigsv=eps_xz*(C44+C55)./2;

    %Average shear stress in the horizontal plane
    sigsh=eps_xz*C66;

    N = 0.125*(C11+C22) -0.25*C12 + 0.5*C66;
    L = 0.5*(C44+C55);

    G = 0.5.*(C55-C44);

end

%Choose viscosity contrast to plot
etac = 0.1;
etaidx = find(etacontr>=etac,1);
if etaidx > etacontrnum
    etaidx = etacontrnum;
end
if etaidx < 1
    disp('Etacontr index < 1. Increase viscosity contrast to be plotted')
    return;
end

%Choose volume fraction to plot
vol = 0.20;
Vstp = 0.05; %same as Volsavestep
volidx = floor(vol/Vstp)+1;
if volidx > volnum
    volidx = volnum;
end
if volidx < 1
    disp('Volume index < 1. Increase volume fraction to be plotted')
    return;
end


%2D Plots
f = figure('OuterPosition',get(0,'ScreenSize').*[1 1 .5 0.9],'Visible','off');
figure(f)
x=zeros(a1a2num,a2a3num); y=x; z=x;
x(:,:)=10.^a2a3(:,:);
y(:,:)=10.^a1a2(:,:);

x1 = x; y1 = y;
%Set values in the upper rigth triangle equal to diagonal values
for i=1:a1a2num
    x(i,i+1:end)=x1(i,i);
    y(i,i+1:end)=y1(i,i);
end

%Deviatoric normal stress during uniaxial compression parallel to z axis
subplot(3,2,1)
z(:,:)=sign(etaidx,volidx,:,:);
plotdem(x,y,z,a1a2num,a2a3num,'{\sigma}_n^*',logscale)

subplot(3,2,3)
z(:,:)=sigsv(etaidx,volidx,:,:);
plotdem(x,y,z,a1a2num,a2a3num,'{\sigma}_{sv}',logscale)

subplot(3,2,4)
z(:,:)=sigsh(etaidx,volidx,:,:);
plotdem(x,y,z,a1a2num,a2a3num,'{\sigma}_{sh}',logscale)

%Normal vs vertical shear viscosity ratio
subplot(3,2,2)
z(:,:)=sign(etaidx,volidx,:,:)./sigsv(etaidx,volidx,:,:);
plotdem(x,y,z,a1a2num,a2a3num,'{\sigma}_n/{\sigma}_{sv}',logscale)

%Radial viscous anisotropy
subplot(3,2,5)
%z(:,:)=sigsh(etaidx,volidx,:,:)./sigsv(etaidx,volidx,:,:);
%plotdem(x,y,z,a1a2num,a2a3num,'{\sigma}_{sh}/{\sigma}_{sv}',logscale)
z(:,:)=N(etaidx,volidx,:,:)./L(etaidx,volidx,:,:);
plotdem(x,y,z,a1a2num,a2a3num,'Radial anisotropy',logscale)

%Azimuthal viscous anisotropy
subplot(3,2,6)
%z(:,:)=(C44(etaidx,volidx,:,:)./C55(etaidx,volidx,:,:)-1).*100;
%plotdem(x,y,z,a1a2num,a2a3num,'C44/C55-1 (%)',logscale)
z(:,:)=G(etaidx,volidx,:,:);
plotdem(x,y,z,a1a2num,a2a3num,'Azimuthal anisotropy',logscale)

if printmod
    filename=strcat('viscousanisotropy',num2str(vol),'.jpg');
    print('-djpeg', '-r300',filename);
    movefile(filename,path);
end
