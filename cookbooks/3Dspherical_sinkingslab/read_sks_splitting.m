clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                Set input parameters                  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../../SKS-SPLIT/viz/')

%Path to input file split.*.staat
splitpath='out/splitting/';
%Path to output file 
outpath='out';

filenumber='0020'; %set number for input split.*.staat and output files, string with 4 digits

%Coordinate system 
cartspher  = 2; % 1 (cartesian) vs 2 (polar/spherical)
dimensions = 3; % 2 (2D) vs 3 (3D)

Surface = 6371e+3; %Set vertical coordinate of surface where to plot the seismic stations

%Set coordinates of additional reference seismic station with 1 sec of delay time
X1ref=90/180*pi;
X2ref=6371e+3;
X3ref=78/180*pi;

printmod=1; % Save (>0) or not (0) the imahe with the splitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

splitname=strcat(splitpath,'split.0.',filenumber,'.sstat'); %set name for input split.*.staat file
outname=strcat('sks.timestep.',filenumber,'.h5'); %set name for Paraview output files
M=dlmread(splitname);

np=size(M,1); %number of seismic stations

xyz=zeros(3,np+1); %seismic stations coordinates
xyz0 = xyz;
sks=zeros(3,np+1); %fast azimuths coordinates

xyz0(1,1:np)  = M(:,1);
xyz0(2,1:np)  = Surface;
xyz0(3,1:np)  = M(:,2);
mean_fazi     = M(:,3);
%std_dev_fazi = M(:,4);
mean_dt       = M(:,5);
%std_dev_dt   = M(:,6);

%Convert fazi into unit vector components
for i=1:np
    if mean_dt(i)
        sks(1,i)=cosd(mean_fazi(i));
        sks(3,i)=sind(mean_fazi(i));
    end
end

%Set position, fast azimuth and delay time (1 sec) for reference station
xyz0(1,np+1)=X1ref;
xyz0(2,np+1)=X2ref;
xyz0(3,np+1)=X3ref;
sks(1,np+1)=1;
mean_dt(np+1)=1;

if(cartspher == 1)
    xyz = xyz0;
else
    %Convert spherical coordinates to cartesian coordinates
    for i=1:np+1
        xyz(1,i)=xyz0(2,i)*sin(xyz0(3,i))*cos(xyz0(1,i));
        xyz(2,i)=xyz0(2,i)*sin(xyz0(3,i))*sin(xyz0(1,i));
        xyz(3,i)=xyz0(2,i)*cos(xyz0(3,i));
    end
end

%Rotate
if(cartspher == 2)
    
    acs=zeros(3,3);
    SKS=zeros(3,np+1);
    for i=1:np+1
        %Z-X-Z: Euler angles in radians
        %Rotate forward with respect to phi (X3s(s)) and theta (X1s(s)) of the seismic station
        %First rotate around X to set Z vertical, then around Z.
        
        phi1 = xyz0(1,i)-pi*3/2; theta = 0; phi2 = 0;
        if(dimensions == 3)
            theta = +xyz0(3,i) - pi/2;
        end
        
        %Transform Euler angles into direction cosine matrix (Z-X-Z)
        acs(1,1)=cos(phi2)*cos(phi1)-cos(theta)*sin(phi1)*sin(phi2);
        acs(2,1)=cos(phi2)*sin(phi1)+cos(theta)*cos(phi1)*sin(phi2);
        acs(3,1)=sin(phi2)*sin(theta);
        
        acs(1,2)=-sin(phi2)*cos(phi1)-cos(theta)*sin(phi1)*cos(phi2);
        acs(2,2)=-sin(phi2)*sin(phi1)+cos(theta)*cos(phi1)*cos(phi2);
        acs(3,2)=cos(phi2)*sin(theta);
        
        acs(1,3)=sin(theta)*sin(phi1);
        acs(2,3)=-sin(theta)*cos(phi1);
        acs(3,3)=cos(theta);
        
        
        for j=1:3
            for k=1:3
                SKS(j,i)=SKS(j,i) + acs(j,k)*sks(k,i);
            end
        end
        
    end
    
    sks = SKS;
    
end

%Scaling measuraments by the delay time
figure(1)
if(cartspher == 1)
    
    sks(1,:)=sks(1,:).*mean_dt';
    sks(3,:)=sks(3,:).*mean_dt';
    h=quiver(-xyz(1,:),xyz(3,:),-sks(1,:),sks(3,:),0.4,'r','LineWidth',1);
    h.ShowArrowHead = 'off';
    hold on
    h=quiver(-xyz(1,:),xyz(3,:),sks(1,:),-sks(3,:),0.4,'r','LineWidth',1);
    h.ShowArrowHead = 'off';
    grid on
    hold off
    axis image
    
else
    
    sks(1,:)=sks(1,:).*mean_dt';
    sks(2,:)=sks(2,:).*mean_dt';
    sks(3,:)=sks(3,:).*mean_dt';
    [x,y,z]=sphere;
    surf(x*X2ref,y*X2ref,z*X2ref)
    %shading interp
    colormap([1 1 1])
    
    view_elev=mean(xyz0(3,:))/pi*180-90;
    if mean(xyz0(3,:))/pi*180>90
        view_elev=-view_elev;
    end
    view_az=mean(xyz0(1,:))/pi*180+90;
    view(view_az,view_elev)
    hold on
    h=quiver3(xyz(1,:),xyz(2,:),xyz(3,:),sks(1,:),sks(2,:),sks(3,:),0.4,'r','LineWidth',1);
    h.ShowArrowHead = 'off';
    h=quiver3(xyz(1,:),xyz(2,:),xyz(3,:),-sks(1,:),-sks(2,:),-sks(3,:),0.4,'r','LineWidth',1);
    h.ShowArrowHead = 'off';
    xlabel('X')
    ylabel('Y')
    hold off
    axis image
    camzoom(2)
    
end

if printmod
    filename=strcat('sks',filenumber,'.png');
    print('-dpng', '-r300',filename);
    movefile(filename,outpath);
end

hdf5write(outname,'/xyz', xyz);
hdf5write(outname,'/sks', sks,'WriteMode','append');
hdf5write(outname,'/dt', mean_dt,'WriteMode','append');
movefile(outname,outpath);

npstr=num2str(np+1);

nl = sprintf('\n'); % new line
stringa=['<?xml version="1.0" ?>', nl, ...
    '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">',nl, nl, ...
    '<Domain>', nl, ...
    '   <Grid Name="materialSwarm" >',nl...
    '      <Time Value="0.000000" />', nl, ...
    '         <Topology Type="POLYVERTEX" NodesPerElement="',npstr,'"> </Topology>',nl,...
    '         <Geometry Type="XYZ">',nl,...
    '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',npstr,32,'3">',outname,':/xyz </DataItem>',nl,...
    '         </Geometry>',nl,nl];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  field_name,file_name,is_vector,dimensions,timestep
ff1=add_swarm_to_XDMF('sks',strcat(outname,':/sks'),1,npstr,0);
ff2=add_swarm_to_XDMF('dt',strcat(outname,':/dt'),0,npstr,0);
ff=[ff1,nl,ff2,nl];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scriptfile=[stringa,nl,ff,nl, '  </Grid>',nl,'</Domain>',nl,'</Xdmf>'];
xmffilename=strcat('XDMF.SKS.',filenumber,'.xmf');
fid = fopen(xmffilename, 'w');
fwrite(fid,scriptfile);
fclose(fid);

movefile(xmffilename,outpath);
