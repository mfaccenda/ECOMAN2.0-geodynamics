%Read DEM
read=1;
if read
    clear all
    read=1;
end

close all;

path = './';

logscale=1;
printmod=1;

if read

    fname = [path,'elastictensordem.h5'];

    C11=h5read(fname,'/C11');
    C22=h5read(fname,'/C22');
    C33=h5read(fname,'/C33');
    C44=h5read(fname,'/C44');
    C55=h5read(fname,'/C55');
    C66=h5read(fname,'/C66');
    C12=h5read(fname,'/C12');
    C13=h5read(fname,'/C13');
    C23=h5read(fname,'/C23');
    a1a2=h5read(fname,'/esa12');
    a2a3=h5read(fname,'/esa23');
    Vol=h5read(fname,'/Vol');
    ro=h5read(fname,'/ro');

    gridnum = h5read(fname,'/gridnum');

    volnum = gridnum(1);
    a1a2num = gridnum(2);
    a2a3num = gridnum(3);

    %P-wave anisotropy
    Vpfast = (C11+C22)/2;
    Vpslow = C33;
    Vpiso = (Vpfast + Vpslow)./2;
    dVpiso = (Vpiso(:,:,:)./Vpiso(:,1,1)-1)*100;
    fp = (Vpfast - Vpslow)./Vpiso*100/2;

    M = 1e+9*(3/15*(C11+C22+C33)+2/15*(C12+C13+C23)+4/15*(C44+C55+C66));
    G = 1e+9*(1/15*(C11+C22+C33)-1/15*(C12+C13+C23)+1/5*(C44+C55+C66));

end

%Choose volume fraction to plot
Vstp = 0.005;
vol = 0.05; %
volidx = floor(vol/Vstp)+1;
if volidx > volnum
    volidx = volnum;
end
if volidx < 1
    disp('Volume index < 1. Increase volume fraction to be plotted')
    return;
end

%2D Plots
f = figure('OuterPosition',get(0,'ScreenSize').*[1 1 .4 0.9],'Visible','off');
figure(f)
x=zeros(a1a2num,a2a3num); y=x; z=x;
x(:,:)=10.^a2a3(:,:);
y(:,:)=10.^a1a2(:,:);

x1 = x; y1 = y;
%Set values in the upper rigth triangle equal to diagonal values
for i=1:a1a2num
    x(i,i+1:end)=x(i,i);
    y(i,i+1:end)=y(i,i);
end

nexttile
z(:,:)=fp(volidx,:,:);
plotdem(x,y,z,a1a2num,a2a3num,'Vp anis (%)',logscale)

nexttile
z(:,:)=dVpiso(volidx,:,:);
plotdem(x,y,z,a1a2num,a2a3num,'dVp (%)',logscale)

nexttile
z(:,:)=(M(volidx,:,:)./ro(volidx)).^0.5;
plotdem(x,y,z,a1a2num,a2a3num,'Vp (m/s)',logscale)

nexttile
z(:,:)=(G(volidx,:,:)./ro(volidx)).^0.5;
plotdem(x,y,z,a1a2num,a2a3num,'Vs (m/s)',logscale)

nexttile
z(:,:)=(M(volidx,:,:)./G(volidx,:,:)).^0.5;
plotdem(x,y,z,a1a2num,a2a3num,'Vp/Vs',logscale)

sgtitle(['Volume fraction ',num2str(vol*100),' %'],'fontsize',24,'fontweight','bold','fontangle','italic');

if printmod
    filename=strcat('dem_elastic',num2str(vol),'.jpg');
    print('-djpeg', '-r300',filename);
    movefile(filename,path);
end