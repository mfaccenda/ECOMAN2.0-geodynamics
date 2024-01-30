function plotfigures(fignum)

global printmod figname output_dir strain stp_strain
global vp vs1 vs2 pp ps1 ps2
global ss cs odf

SS = get(0,'ScreenSize');

if fignum==1
    
    f = figure('Name','Vp,Vs1,dVs',...% left corner (x,y)    | of screen size (x,y)
            'NumberTitle','off','OuterPosition',[ 0.2*SS(3)   0.6*SS(4)  0.9*SS(3)   0.7*SS(4)],'Visible','off');

    % plotting convention - plot a-axis to east
    plota2east
    
    % set colour map to seismic color map : blue2redColorMap
    %setMTEXpref('defaultColorMap',blue2redColorMap)
    setMTEXpref('defaultColorMap',red2blueColorMap);
    setMTEXpref('xAxisDirection','east')
    setMTEXpref('FontStyle','Arial');
    
    % some options
    blackMarker = {'Marker','s','MarkerSize',10,'antipodal',...
        'MarkerEdgeColor','white','MarkerFaceColor','black','doNotDraw'};
    whiteMarker = {'Marker','o','MarkerSize',10,'antipodal',...
        'MarkerEdgeColor','black','MarkerFaceColor','white','doNotDraw'};
    
    % some global options for the titles
    %titleOpt = {'FontSize',getMTEXpref('FontSize'),'visible','on'}; %{'FontSize',15};
    titleOpt = {'visible','on','color','k','FontSize',20};
    
    % Setup multiplot
    % define plot size [origin X,Y,Width,Height]
    mtexFig = mtexFigure('position',[0 0 1000 500],'Name','VpVs1dVs','NumberTitle','off');
    
    % set up spacing between subplots default is 10 pixel
    mtexFig.innerPlotSpacing = 15;
    
    % Standard Seismic plot with 8 subplots in 3 by 3 matrix
    %
    % Plot matrix layout
    %        1 Vp        2 AVs      3 AVs + S1 polarizations
    %
    %**************************************************************************
    % Vp : Plot P-wave velocity (km/s)
    %**************************************************************************
    
    % Plot P-wave velocity (km/s)
    plot(vp,'contourf','complete','upper','contourf')
    mtexTitle('Vp (km/s)',titleOpt{:})
    mtexColorbar%('FontWeight','bold')
    
    
    % extrema
    [maxVp, maxVpPos] = max(vp);
    [minVp, minVpPos] = min(vp);
    
    % percentage anisotropy
    AVp = 200*(maxVp-minVp) / (maxVp+minVp);
    
    % mark maximum with black square and minimum with white circle
    hold on
    plot(maxVpPos,blackMarker{:})
    plot(minVpPos,whiteMarker{:})
    hold off
    
    % subTitle
    xlabel(['Vp Anisotropy = ',num2str(AVp,'%6.1f')],titleOpt{:})
    
    % create a new axis
    nextAxis
    
    % Plot S1-wave velocity (km/s) 
    plot(vs1,'contourf','complete','upper','contourf');
    mtexTitle('Vs1 (km/s)',titleOpt{:})
    mtexColorbar

    % Max/min vs1
    [maxvs1,maxvs1Pos] = max(vs1);
    [minvs1,minvs1Pos] = min(vs1);
   
    %xlabel(['Max Vs1 = ',num2str(maxvs1,'%6.1f')],titleOpt{:})
    
    % mark maximum with black square and minimum with white circle
    hold on
    plot(maxvs1Pos,blackMarker{:})
    plot(minvs1Pos,whiteMarker{:})
    hold off
    
    % mark crystal axes
    text([xvector,yvector,zvector],{'X ','Y ','Z'},...
        'backgroundcolor','w','doNotDraw');
    
    % create a new axis
    nextAxis
    
    %Plot S-wave anisotropy (percent)
    AVs = 200*(vs1-vs2)./(vs1+vs2);

    % Max/min percentage anisotropy
    [maxAVs,maxAVsPos] = max(AVs);
    [minAVs,minAVsPos] = min(AVs);

    mtexTitle('AVs (%), Vs1 polarization',titleOpt{:})
    
    plot(AVs,'contourf','complete','upper','contourf');
    mtexColorbar
    
    xlabel(['\gamma = ',num2str(strain(stp_strain),'%5.1f')],titleOpt{:})
    
    hold on
    plot(ps1,'linewidth',2,'color','black','doNotDraw')
    plot(maxAVsPos,blackMarker{:})
    plot(minAVsPos,whiteMarker{:})
    hold off
    
    f.Position = [123   831   1273   419];

    mtexTitle('AVs (%) + Vs1 polarization',titleOpt{:})

    if printmod>0
        print('-dpng', '-r300',[output_dir,figname]);
    end
    
    close(figure(fignum))
    
end

if fignum == 2 || fignum == 3
    
    if fignum == 2
        f = figure('Name','LPO Phase1',...% left corner (x,y)    | of screen size (x,y)
            'NumberTitle','off','OuterPosition',[ 0.2*SS(3)   0.6*SS(4)  0.9*SS(3)   0.7*SS(4)],'Visible','off');
    end
    if fignum == 3
        f = figure('Name','LPO Phase2',...% left corner (x,y)    | of screen size (x,y)
            'NumberTitle','off','OuterPosition',[ 0.2*SS(3)   0.6*SS(4)  0.9*SS(3)   0.7*SS(4)],'Visible','off');
    end

    %titleOpt = {'visible','on','color','k','FontSize',44};
    titleOpt = {'visible','on','color','k','FontSize',20,'FontStyle','Arial'};

    setMTEXpref('defaultColorMap',WhiteJetColorMap);
    setMTEXpref('xAxisDirection','east')
    setMTEXpref('FontStyle','Arial');

    plotPDF(odf,[Miller(1,0,0,cs)],'antipodal','contourf');
    %caxis([0 5.0])
    mtexTitle('[100]',titleOpt{:})
    mtexColorbar

    % create a new axis
    nextAxis
    plotPDF(odf,[Miller(0,1,0,cs)],'antipodal','contourf');
    %caxis([0 5.0])
    mtexTitle('[010]',titleOpt{:})
    mtexColorbar

    % create a new axis
    nextAxis
    plotPDF(odf,[Miller(0,0,1,cs)],'antipodal','contourf');
    %caxis([0 5.0])
    mtexTitle('[001]',titleOpt{:})
    mtexColorbar
    %annotate([xvector,yvector,zvector],'label',{'x','y','z'},'backgroundcolor','w');
    xlabel(['\gamma = ',num2str(strain(stp_strain),'%5.1f')],'FontSize',20)
    
    f.Position = [123   831   1273   419];

    mtexTitle('[001]',titleOpt{:})
    
    if printmod>0
        print('-dpng', '-r300',[output_dir,figname]);
    end
    
    close(figure(fignum))
    
end