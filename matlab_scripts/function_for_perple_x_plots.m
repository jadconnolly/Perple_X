function [] = function_for_perple_x_plots (x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,LineStyle,LineWidth,Marker,FontSize,titl)
% Generic function to make 2- and 3-d plots from Perple_X tab format files.
%                                               JADC, 5/2011.
%                                                                
% Modifications:
%
%   3-d surface options moved to allow vector export of contour plots. square
%   axes applied to all plots. 
%                                               Philippe Goncalves, 2/2012. 

fig1 = figure(1);

if nvar == 1 % two cases: 1d - table -> 2d plot
    
    [kvar, ok] = listdlg('PromptString','Select the INDEPENDENT (X) variable:','ListSize',[240 500],'SelectionMode','single','ListString',dnames{1});
    if ok == 0, errordlg(['You did not choose a variable, I quit!']), end;
    
    if mvar > 2,
        [dvar, ok] = listdlg('PromptString','Select the DEPENDENT (Y) variables:','ListSize',[240 500],'ListString',dnames{1});
        if ok == 0, errordlg(['You did not choose a variable, I quit!']), end;
        [n m] = size(dvar);
        jvar = n*m;
    else   % mvar must be 2
        dvar = 2;
        jvar = 1;
    end
    
    hold all
    
    for i = 1:jvar,plot(a(kvar,1:nrow),a(dvar(i),1:nrow),'LineStyle',LineStyle,'LineWidth',LineWidth,'Marker',Marker),end
    
    legend(dnames{1}{dvar},'Location','EastOutside'); 
    %axis square; 
    %legend HIDE
    axis tight; xlabel(dnames{1}{kvar},'FontSize',FontSize); title(titl,'FontSize',FontSize);
    
elseif nvar == 2 % 2d - table -> 2/3d plot
    
    amin = min(a(:)); amax = max(a(:)); disp(['Grid data range is ',num2str(amin),' - >',num2str(amax)])
    
    choice = questdlg('Select plot style','Plot Style','3D Surface','Auto-Contour','Contour','3D Surface');
    
    switch choice;
        
        case '3D Surface';
            surf(x,y,a); d2 = 0; light; shading interp; lighting phong; zlabel(zname); 
            colorbar; %uncomment for colorbar
        case 'Auto-Contour';
            [C,h]=contour(x,y,a); clabel(C,h); d2 = 1;
        case 'Contour';
             prompt = {'Minimum contour:','Maximum contour:','Contour interval:'};
             dlg = 'Contour specification';
             num_lines = 1;
             da = (amax-amin)/11;
             def = {num2str(amin+da/2),num2str(amax-da/2),num2str(da)};
             c = inputdlg(prompt,dlg,num_lines,def);
             helpdlg('Carefully select contours for labeling. When done, press RETURN while the Graph window is the active window.');
             contours = [str2num(c{1}):str2num(c{3}):str2num(c{2})];
%            contours = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15, 20, 30, 40]
           % contours = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02,  0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 20, 30, 40]
            [C,h]= contour(x,y,a,contours);clabel(C,h,'manual'); d2 = 1;
 
    end
    
    if d2 == 1,
        if strcmp(titl,' ')
            titl = zname;
        else
            titl = [titl ', ' zname];
        end
    end
    
    axis square; axis tight; xlabel(xname); ylabel(yname); title(titl);
    
end

