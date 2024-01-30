function plotazirad(plot_Voigt,plot_Reuss,plot_Mixed)

global printmod input_dir output_dir fname0 figname filename strain strainnum

Rs = zeros(strainnum,1); Rp = Rs; As = Rs; Ap = Rs; ETAeff = Rs;

v = [0      0.4470 0.7410;
     0.8500 0.3250 0.0980;
     0.9290 0.6940 0.1250;
     0.4940 0.1840 0.5560;
     0.4660 0.6740 0.1880;
     0.3010 0.7450 0.9330;
     0.6350 0.0780 0.1840];

if plot_Voigt == 0 && plot_Reuss == 0 && plot_Mixed == 0
    disp('NO AZIRAD plot: set plot_Voigt, plot_Reuss or plot_Mixed > 0')
    return
end

f = figure('Name','Azimuthal and Radial anisotropy','NumberTitle','off'); clf;

if plot_Voigt
    
    for stp_strain = 1:strainnum
        
        filename = [input_dir,fname0,num2str(strain(stp_strain),'%5.1f'),'.h5'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  Read stiff matrix and plot strength of anisotropies %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %Read output elastic tensor
        Cstiff = h5read(filename,'/Voigt');
        
        A=3/8*(Cstiff(1,1)+Cstiff(3,3))+1/4*Cstiff(1,3)+0.5*Cstiff(5,5);
        C=Cstiff(2,2);
        L=0.5*(Cstiff(4,4)+Cstiff(6,6));
        N=1/8*(Cstiff(1,1)+Cstiff(3,3))-1/4*Cstiff(1,3)+0.5*Cstiff(5,5);
        F=0.5*(Cstiff(1,2)+Cstiff(2,3));
        
        
        %calculate strength of azim and rad anisotropy
        Rs(stp_strain)=((N/L)^0.5 - 1)*100;
        Rp(stp_strain)=((A/C)^0.5 - 1)*100;
        As(stp_strain)=((Cstiff(6,6)/Cstiff(4,4))^0.5 - 1)*100;
        Ap(stp_strain)=((Cstiff(1,1)/Cstiff(3,3))^0.5 - 1)*100;
        ETAeff(stp_strain)=F/(A-2*L);
        
    end
    
    
    subplot(2,1,1)
    plot(strain,Rs,'-','Color',v(1,:),'LineWidth',1.5);
    hold on
    plot(strain,Rp,'-','Color',v(2,:),'LineWidth',1.5);
    plot(strain,As,'-','Color',v(3,:),'LineWidth',1.5);
    plot(strain,Ap,'-','Color',v(4,:),'LineWidth',1.5);
    hold off
    ylabel('Rs Rp As Ap (%)')
    xlabel('\gamma strain ')
    legend('Rs', 'Rp', 'As', 'Ap','Location','NorthWest')
    subplot(2,1,2)
    plot(strain,ETAeff,'k-','LineWidth',1.5)
    legend('ETA effective','Location','NorthEast')
    ylabel('\eta')
    xlabel('\gamma strain ')
    
end

if plot_Reuss
    
    for stp_strain = 1:strainnum
        
        filename = [input_dir,fname0,num2str(strain(stp_strain),'%5.1f'),'.h5'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  Read stiff matrix and plot strength of anisotropies %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %Read output elastic tensor
        Cstiff = h5read(filename,'/Reuss');
        
        
        A=3/8*(Cstiff(1,1)+Cstiff(3,3))+1/4*Cstiff(1,3)+0.5*Cstiff(5,5);
        C=Cstiff(2,2);
        L=0.5*(Cstiff(4,4)+Cstiff(6,6));
        N=1/8*(Cstiff(1,1)+Cstiff(3,3))-1/4*Cstiff(1,3)+0.5*Cstiff(5,5);
        F=0.5*(Cstiff(1,2)+Cstiff(2,3));
        
        
        %calculate strength of azim and rad anisotropy
        Rs(stp_strain)=((N/L)^0.5 - 1)*100;
        Rp(stp_strain)=((A/C)^0.5 - 1)*100;
        As(stp_strain)=((Cstiff(6,6)/Cstiff(4,4))^0.5 - 1)*100;
        Ap(stp_strain)=((Cstiff(1,1)/Cstiff(3,3))^0.5 - 1)*100;
        ETAeff(stp_strain)=F/(A-2*L);
        
    end
    
    
    subplot(2,1,1)
    hold on
    plot(strain,Rs,'--','Color',v(1,:),'LineWidth',1.5);
    plot(strain,Rp,'--','Color',v(2,:),'LineWidth',1.5);
    plot(strain,As,'--','Color',v(3,:),'LineWidth',1.5);
    plot(strain,Ap,'--','Color',v(4,:),'LineWidth',1.5);
    hold off
    ylabel('Rs Rp As Ap (%)')
    xlabel('\gamma strain ')
    legend('Rs', 'Rp', 'As', 'Ap','Location','NorthWest')
    subplot(2,1,2)
    hold on
    plot(strain,ETAeff,'--k','LineWidth',1.5)
    hold off
    legend('ETA effective','Location','NorthEast')
    ylabel('\eta')
    xlabel('\gamma strain ')
    
end

if plot_Mixed
    
    for stp_strain = 1:strainnum
        
        filename = [input_dir,fname0,num2str(strain(stp_strain),'%5.1f'),'.h5'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%  Read stiff matrix and plot strength of anisotropies %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %Read output elastic tensor
        Cstiff = h5read(filename,'/Mixed');
        
        
        A=3/8*(Cstiff(1,1)+Cstiff(3,3))+1/4*Cstiff(1,3)+0.5*Cstiff(5,5);
        C=Cstiff(2,2);
        L=0.5*(Cstiff(4,4)+Cstiff(6,6));
        N=1/8*(Cstiff(1,1)+Cstiff(3,3))-1/4*Cstiff(1,3)+0.5*Cstiff(5,5);
        F=0.5*(Cstiff(1,2)+Cstiff(2,3));
        
        
        %calculate strength of azim and rad anisotropy
        Rs(stp_strain)=((N/L)^0.5 - 1)*100;
        Rp(stp_strain)=((A/C)^0.5 - 1)*100;
        As(stp_strain)=((Cstiff(6,6)/Cstiff(4,4))^0.5 - 1)*100;
        Ap(stp_strain)=((Cstiff(1,1)/Cstiff(3,3))^0.5 - 1)*100;
        ETAeff(stp_strain)=F/(A-2*L);
        
    end
    
    
    subplot(2,1,1)
    hold on
    plot(strain,Rs,'+','Color',v(1,:),'LineWidth',1.5);
    plot(strain,Rp,'+','Color',v(2,:),'LineWidth',1.5);
    plot(strain,As,'+','Color',v(3,:),'LineWidth',1.5);
    plot(strain,Ap,'+','Color',v(4,:),'LineWidth',1.5);
    hold off
    ylabel('Rs Rp As Ap (%)')
    xlabel('\gamma strain ')
    legend('Rs', 'Rp', 'As', 'Ap','Location','NorthWest')
    subplot(2,1,2)
    hold on
    plot(strain,ETAeff,'+k','LineWidth',1.5)
    hold off
    legend('ETA effective','Location','NorthEast')
    ylabel('\eta')
    xlabel('\gamma strain ')
    
end

if printmod
    print('-dpng', '-r150',[output_dir,figname]);
end
