%clear;
close all;
clc;

global printmod figname input_dir output_dir fname0 filename strain stp_strain strainnum
global vp vs1 vs2 pp ps1 ps2
global cs ss odf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                Set input parameters                  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Path to D-REX_S output files directory
input_dir = 'Atype_/';
%Output directory, where images and videos are saved
output_dir = input_dir;
%First part of the D-REX_S output file name (usually same as output_dir)
fname0 = 'Atype_';

%Min/Step/Max strain as filenumber of D-REX_S output files
min_strain = 0.1;
stp_strain = 0.1;
max_strain = 1.0;

%Choose what to plot/save (No = false; Yes = true)
plot_Cij_LPO = true; %Activate plotting tensors and LPO. Set to false if want to make only movies of images saved in a previous run

plot_Voigt = true; %Plot Vp, dVs due to aggregate LPO, Voigt average
plot_Reuss = true; %Plot Vp, dVs due to aggregate LPO, Reuss average
plot_Mixed = true; %Plot Vp, dVs due to aggregate LPO, mixed Voigt/Reuss average
plot_Phase1LPO = true; %Plot LPO pole figures of main anisotropic phase
plot_Phase2LPO = true; %Plot LPO pole figures of enstatite (only for upper mantel aggr.)

makevideo = true; %Make .mp4 movie of the  LPO evolution

plot_azirad = true; % Plot azimuthal, radial anisotropy as a function of strain

plot_dec = true; % Plot components of elastic tensor as a function of strain

plot_SPO   = false; %Plot Vp, dVs when an SPO model is superimposed on the LPO (i.e., spomod > 0)

printmod = true; %Save the MatLab figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strain=min_strain:stp_strain:max_strain;
strainnum = int16((max_strain-min_strain)/stp_strain+1);

if ~exist(output_dir,'dir')
    system(['mkdir ',output_dir]);
end
%system(['cd ',curdir]);
if plot_Phase1LPO
    MI = zeros(strainnum,1); JI = MI;
end
if plot_Phase2LPO
    MI2 = zeros(strainnum,1); JI2 = MI2;
end

for stp_strain = 1:strainnum

    if plot_Cij_LPO

        figname0 = [fname0,num2str(strain(stp_strain),'%5.1f')];

        filename = strcat(input_dir,[figname0,'.h5']);

        %Read aggregate density
        rho = h5read(filename,'/Density');
        %Read aggregate rocktype
        rocktype = h5read(filename,'/Rocktype');

        %Make sure that LPO of 2nd phase is plotted only for upper mantle aggregates
        if rocktype ~= 1
            plot_Phase2LPO = 0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  Read and plot aggregate elastic tensor %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Set triclin symmetry for the sample
        cs = crystalSymmetry('-1');
        ss = specimenSymmetry('-1');
        
        if plot_Voigt

            %Read output elastic tensor
            C = h5read(filename,'/Voigt');

            %Plotting with MTex
            C=stiffnessTensor(C.*1e+9,cs,'density',rho);
            [vp,vs1,vs2,pp,ps1,ps2] = velocity(C);
            
            figname = ['VpdVs_Voigt_',figname0,'.png'];
            plotfigures(1)

        end

        if plot_Reuss

            %Read output elastic tensor
            C = h5read(filename,'/Reuss');

            %Plotting with MTex
            C=stiffnessTensor(C.*1e+9,cs,'density',rho);
            [vp,vs1,vs2,pp,ps1,ps2] = velocity(C);
            
            figname = ['VpdVs_Reuss_',figname0,'.png'];
            plotfigures(1)

        end

        if plot_Mixed

            %Read output elastic tensor
            C = h5read(filename,'/Mixed');

            %Plotting with MTex
            C=stiffnessTensor(C.*1e+9,cs,'density',rho);
            [vp,vs1,vs2,pp,ps1,ps2] = velocity(C);
            
            figname = ['VpdVs_Mixed_',figname0,'.png'];
            plotfigures(1)

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  Read and plot phase 1 LPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if plot_Phase1LPO

            %Need to set symmetry to orthorhombic otherwise does not plot antipodal mode
            cs = crystalSymmetry('mmm');%'mmm',[4.7646,10.2296,5.9942],'x||a','z||c','mineral','Olivine');
            ss = specimenSymmetry('mmm');%'mmm',[4.7646,10.2296,5.9942],'x||a','z||c','mineral','Olivine');

            acs = h5read(filename,'/acs');
            odfdrex = h5read(filename,'/odf');

            acs=permute(acs,[3,2,1]);
            numgrains=size(acs,3);

            %Orientation function
            ori   = orientation('matrix',acs,{cs},{ss});

            %Weights
            w = odfdrex.*numgrains;

            odf = calcDensity(ori,'weights',w,'halfwidth',10/180*pi);

            MI(stp_strain) = calcMIndex(odf);
            JI(stp_strain) = norm(odf).^2;

            figname = ['LPO_Phase1_',figname0,'.png'];
            plotfigures(2)

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  Read and plot Enstatite LPO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plot_Phase2LPO

            %Need to set symmetry to orthorhombic otherwise does not plot antipodal mode
            cs = crystalSymmetry('mmm');%'mmm',[4.7646,10.2296,5.9942],'x||a','z||c','mineral','Olivine');
            ss = specimenSymmetry('mmm');%'mmm',[4.7646,10.2296,5.9942],'x||a','z||c','mineral','Olivine');

            acs = h5read(filename,'/acs_ens');
            odfdrex = h5read(filename,'/odf_ens');

            acs=permute(acs,[3,2,1]);
            numgrains=size(acs,3);

            %Orientation function
            ori   = orientation('matrix',acs,{cs},{ss});

            %Weights
            w = odfdrex.*numgrains;

            odf = calcDensity(ori,'weights',w,'halfwidth',10/180*pi);

            MI2(stp_strain) = calcMIndex(odf);
            JI2(stp_strain) = norm(odf).^2;

            figname = ['LPO_Phase2_',figname0,'.png'];
            plotfigures(3)

        end

    end

end

if plot_Cij_LPO && strainnum > 1

    if plot_Phase1LPO
        save('MI_Phase1LPO.mat','MI')
        movefile('MI_Phase1LPO.mat',output_dir);

        save('JI_Phase1LPO.mat','JI')
        movefile('JI_Phase1LPO.mat',output_dir);

        f = figure('Name','M-index J-index','NumberTitle','off');
        figure(f)
        plot(strain, MI,'LineWidth',1.5)

        hold on
        plot(strain, JI./100,'LineWidth',1.5)
        axis([0 max_strain 0 0.5])
        ylabel('M-index / J-Index')
        xlabel('\gamma strain ')
        legend('M-index','J-index/100')
        hold off

        figname=['MI_JI_1_',fname0,'.png'];
        if printmod
            print('-dpng', '-r150',[output_dir,figname]);
        end

    end

    if plot_Phase2LPO
        save('MI_Phase2LPO.mat','MI2')
        movefile('MI_Phase2LPO.mat',output_dir);

        save('JI_Phase2LPO.mat','JI2')
        movefile('JI_Phase2LPO.mat',output_dir);

        f = figure('Name','M-index J-index','NumberTitle','off');
        figure(f)
        plot(strain, MI2,'LineWidth',1.5)

        hold on
        plot(strain, JI2./100,'LineWidth',1.5)
        axis([0 max_strain 0 0.5])
        ylabel('M-index / J-Index')
        xlabel('\gamma strain ')
        legend('M-index','J-index/100')
        hold off

        figname=['MI_JI_2_',fname0,'.png'];
        if printmod
            print('-dpng', '-r150',[output_dir,figname]);
        end
    end

end

if plot_azirad

    figname = 'AziRad.png';
    plotazirad(plot_Voigt,plot_Reuss,plot_Mixed)

end

if plot_dec

    figname = 'dec.png';
    plotdec()

end

if makevideo

    dirfile             =   dir([output_dir,'*.png']);
    NFile               =   length(dirfile);

    if plot_Voigt

        % prepare to write to the video
        filename=['VpdVs_Voigt_',fname0,'.mp4'];
        vidObj              =   VideoWriter(filename,"MPEG-4");
        vidObj.FrameRate    =   4; % frames default is 30
        open(vidObj);

        for stp_strain = 1:strainnum
            fname           =   ['VpdVs_Voigt_',fname0,num2str(strain(stp_strain),'%5.1f'),'.png'];
            data_file       =   fullfile(output_dir, fname);
            [X,MAP]         =   imread(data_file);
            f               =   im2frame(uint8(X),MAP);
            writeVideo(vidObj,f);
        end
        close(vidObj);
        movefile(filename,output_dir);

    end

    if plot_Reuss

        % prepare to write to the video
        filename=['VpdVs_Reuss_',fname0,'.mp4'];
        vidObj              =   VideoWriter(filename,"MPEG-4");
        vidObj.FrameRate    =   4; % frames default is 30
        open(vidObj);

        for stp_strain = 1:strainnum
            fname           =   ['VpdVs_Reuss_',fname0,num2str(strain(stp_strain),'%5.1f'),'.png'];
            data_file       =   fullfile(output_dir, fname);
            [X,MAP]         =   imread(data_file);
            f               =   im2frame(uint8(X),MAP);
            writeVideo(vidObj,f);
        end
        close(vidObj);
        movefile(filename,output_dir);

    end
    if plot_Mixed

        % prepare to write to the video
        filename=['VpdVs_Mixed_',fname0,'.mp4'];
        vidObj              =   VideoWriter(filename,"MPEG-4");
        vidObj.FrameRate    =   4; % frames default is 30
        open(vidObj);

        for stp_strain = 1:strainnum
            fname           =   ['VpdVs_Mixed_',fname0,num2str(strain(stp_strain),'%5.1f'),'.png'];
            data_file       =   fullfile(output_dir, fname);
            [X,MAP]         =   imread(data_file);
            f               =   im2frame(uint8(X),MAP);
            writeVideo(vidObj,f);
        end
        close(vidObj);
        movefile(filename,output_dir);

    end

    if plot_Phase1LPO

        % prepare to write to the video
        filename=['LPO_Phase1_',fname0,'.mp4'];
        vidObj              =   VideoWriter(filename,"MPEG-4");
        vidObj.FrameRate    =   4; % frames default is 30
        open(vidObj);

        for stp_strain = 1:strainnum
            fname           =   ['LPO_Phase1_',fname0,num2str(strain(stp_strain),'%5.1f'),'.png'];
            data_file       =   fullfile(output_dir, fname);
            [X,MAP]         =   imread(data_file);
            f               =   im2frame(uint8(X),MAP);
            writeVideo(vidObj,f);
        end
        close(vidObj);
        movefile(filename,output_dir);

    end

    if plot_Phase2LPO

        % prepare to write to the video
        filename=['LPO_Phase2_',fname0,'.mp4'];
        vidObj              =   VideoWriter(filename,"MPEG-4");
        vidObj.FrameRate    =   4; % frames default is 30
        open(vidObj);

        for stp_strain = 1:strainnum
            fname           =   ['LPO_Phase2_',fname0,num2str(strain(stp_strain),'%5.1f'),'.png'];
            data_file       =   fullfile(output_dir, fname);
            [X,MAP]         =   imread(data_file);
            f               =   im2frame(uint8(X),MAP);
            writeVideo(vidObj,f);
        end
        close(vidObj);
        movefile(filename,output_dir);

    end

end


if plot_SPO > 0

    filename = strcat(input_dir,fname0,'SPO.h5');

    %Read aggregate density
    rho = h5read(filename,'/Density');
    %Read spomod
    spomod = h5read(filename,'/spomod');

    %Set triclinic symmetry for the sample
    cs = crystalSymmetry('-1');
    ss = specimenSymmetry('-1');

    %Read output elastic tensor
    C = h5read(filename,'/Mixed');

    %Plotting with MTex
    C=stiffnessTensor(C.*1e+9,cs,'density',rho);
    [vp,vs1,vs2,pp,ps1,ps2] = velocity(C);
    %Rotate such that x is EW, y is vertical and z is NS
    %r = rotation('Euler',0*degree,90*degree,0*degree);
    %C=rotate(C,r);

    figname = ['VpdVs_SPO=',num2str(spomod),'.png'];
    plotfigures(1)

end