function [] = function_for_perple_x_plots (x,y,a,symb,xname,yname,zname, ...
    nvar,mvar,nrow,dnames,LineStyle,LineWidth,Marker,FontSize,titl,type)
%
% Generic function to make 2- and 3-d plots from Perple_X tab format files.
%                                               JADC, 5/2011.
%
% Modifications:
%
%   3-d surface options moved to allow vector export of contour plots.
%   square axes applied to all plots.
%                                               Philippe Goncalves, 2/2012.
%
%   added point file plots and type flag to indicate plot type
%     1 - 2d tab file plot
%     2 - 1d tab file plot
%     3 - 2d pts file plot
%                                               JADC, 2/2025

fig1 = figure(1);

choices2d = {'3D Surface','3D Auto-Contour','3D Manual-Contour','2D Color-Filled Contour'};

if type == 2 % three cases according to "type"

    [kvar, ok] = listdlg('PromptString','Select the INDEPENDENT (X) variable:','ListSize',[240 500],'SelectionMode','single','ListString',dnames{1});
    if ok == 0, errordlg('You did not choose a variable, I quit!'), end

    if mvar > 2
        [dvar, ok] = listdlg('PromptString','Select the DEPENDENT (Y) variables:','ListSize',[240 500],'ListString',dnames{1});
        if ok == 0, errordlg('You did not choose a variable, I quit!'), end
        [n, m] = size(dvar);
        jvar = n*m;
    else   % mvar must be 2
        dvar = 2;
        jvar = 1;
    end

    hold on

    for i = 1:jvar,plot(a(kvar,1:nrow),a(dvar(i),1:nrow),'LineStyle',LineStyle,'LineWidth',LineWidth,'Marker',Marker),end

    legend(dnames{1}{dvar},'Location','EastOutside');
    %axis square;
    %legend HIDE
    axis tight; xlabel(dnames{1}{kvar},'FontSize',FontSize); title(titl,'FontSize',FontSize);

elseif type == 1 % 2d - table -> 2/3d plot

    amin = min(a(:)); amax = max(a(:)); disp(['Grid data range is ',num2str(amin),' - >',num2str(amax)])

    [idx, tf] = listdlg('ListString', choices2d,...
        'SelectionMode', 'Single', 'PromptString', 'Select plot style', 'Initialvalue', 1,'Name', 'Make choice');

    if ~tf, display('You did not choose a plot type, I quit!'),return,end

    if idx == 1

        surf(x,y,a); d2 = 0; light; shading interp; lighting phong; zlabel(zname);
        colorbar; %uncomment for colorbar

    elseif idx == 2

        [C,h]=contour(x,y,a); clabel(C,h); d2 = 1;

    elseif idx == 3

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

    elseif idx ==4

        [C,h]=contourf(x,y,a); clabel(C,h); d2 = 1;

    end

    if d2 == 1
        if strcmp(titl,' ')
            titl = zname;
        else
            titl = [titl ', ' zname];
        end
    end

    axis square; axis tight; xlabel(xname); ylabel(yname); title(titl);

elseif type == 3 % 2d pts file plots
    %                                                symbol list
    sc(1:12) = ['b','r','b','r','b','r','b','r','b','r','b','r'];
    mt(1:12) = ['o','o','s','s','v','v','o','o','s','s','v','v'];
    %                                                 first figure
    hold on
    isym = 0;

    isym1 = 0;
    isym2 = 0;
    isym3 = 0;
    isym4 = 0;
    %                               group points by symbol
    indx(1:5,1:nrow) = 0;
    %                               best central model symbol = 1
    j = find(symb==1);
    jct(1) = size(j,1);
    indx(1,1:size(j,1)) = j;
    %                               perturbed model symbol = 2
    j = find(symb==2);
    jct(2) = size(j,1);
    indx(2,1:size(j,1)) = j;
    %                               central model symbol = 3
    j = find(symb==3);
    jct(3) = size(j,1);
    indx(3,1:size(j,1)) = j;

    if indx(2,1) > 0
    %                               best perturbed model
        j = find(min(a(indx(2,1:jct(2)))))
        if j ~= indx(1,1)
        jct(4) = 1;
        indx(4,1) = j;
        end
    end

    j = find(symb>=4);
    indx(5,1:size(j,1)) = j;

    if indx(5,1) > 0
        errordlg ('invalid symbol code in *.pts file, valid codes [1,4]');
    end

    hold on

    iplot = 0;
    jplot = 0;

    for i = 1:4

        if indx(i,1) > 0

            iplot = iplot + 1;

            if i == 1

                plot3 (x(indx(i,1)),y(indx(i,1)),a(indx(i,1)),['o'], ...
                    'MarkerFaceColor',"#EDB120",...
                    'MarkerEdgeColor','k',...
                    'MarkerSize',20);

                strg = "Best central model" + ...
                    newline + ...
                    ['  ',xname,' = ',num2str(x(indx(i,1)),'%10.1f')] + ...
                    newline + ...
                    ['  ',yname,' = ',num2str(y(indx(i,1)),'%10.1f')];

                t(iplot) = cellstr(strg);

            elseif i == 2 || i == 3

                isymb = symb(indx(i,1))
                fc = sc(isymb);
                sz = 6;

                plot3 (x(indx(i,1:jct(i))),y(indx(i,1:jct(i))),...
                    a(indx(i,1:jct(i))),[fc,mt(isymb)], ...
                    'MarkerFaceColor',fc,'MarkerSize',sz);

                if isymb == 3
                    t(iplot) = cellstr('Central model tries');
                else
                    t(iplot) = cellstr('Central model perturbations');
                end

                jplot = iplot;

            elseif i == 4

                plot3 (x(indx(i,1)),y(indx(i,1)),a(indx(i,1)),['o'], ...
                    'MarkerFaceColor',"#0072BD",...
                    'MarkerEdgeColor','none',...
                    'MarkerSize',12);

                strg = "Best perturbed central model" + ...
                    newline + ...
                    ['  ',xname,' = ',num2str(x(indx(i,1)),'%10.1f')] + ...
                    newline + ...
                    ['  ',yname,' = ',num2str(y(indx(i,1)),'%10.1f')];

                t(iplot) = cellstr(strg);

            end
        end
    end

    hold off

    choice = questdlg('Fit a surface to the points','Spline fit','Yes','No','Yes');

    switch choice
        case 'Yes'
            hold on
            ok = 0;

        case 'No'
            ok = 1;
    end

    if ok == 0

        % grid the data
        % ngrid = int8(max(size(x))/3); % Adjust the number of points as needed
        ngrid = 10;
        [xgrid, ygrid] = ndgrid(linspace(min(x),max(x),ngrid),linspace(min(y),max(y),ngrid));
        % linear, cubic, natural, nearest
        method = 'nearest';
        agrid = griddata(x, y, a, xgrid, ygrid, method);

        amin = min(agrid(:)); amax = max(agrid(:));
        disp(['Grid data range is ',num2str(amin),' - >',num2str(amax)])

        [idx, tf] = listdlg('ListString', choices2d,...
            'SelectionMode', 'Single', 'PromptString', ...
            'Select plot style', 'Initialvalue', 1,'Name', 'Make choice');

        if ~tf, display('You did not choose a plot type, I quit!'),return,end

        if idx == 1

            surf(xgrid,ygrid,agrid); d2 = 0; light; shading interp; lighting phong; zlabel(zname);
            colorbar; %uncomment for colorbar

        elseif idx == 2

            [C,h]=contour(xgrid,ygrid,agrid); clabel(C,h); d2 = 1;

        elseif idx == 3

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
            [C,h]= contour(xgrid,ygrid,agrid,contours);clabel(C,h,'manual'); d2 = 1;

        elseif idx ==4

            [C,h]=contourf(xgrid,ygrid,agrid); clabel(C,h); d2 = 1;

        end

    end

    choice = questdlg('Plot error ellipse','Covariance','Yes','No','Yes');

    switch choice
        case 'Yes'
            hold on
            [avg,covariance] = function_to_plot_covariance_ellipse([x y]);

            if jplot ~= 0
                strg = string(t(jplot)) + newline + ...
                    ['  \mu_{',xname,'} = ',num2str(avg(1),'%10.1f'),...
                ' \pm ',num2str(sqrt(covariance(1,1)),'%10.1f')] + ...
                newline + ['  \mu_{',yname,'} = ',num2str(avg(2),'%10.1f'),...
                ' \pm ',num2str(sqrt(covariance(2,2)),'%10.1f')];
                t(jplot) = cellstr(strg);
            end

        case 'No'
            ok = 1;
    end

    if iplot > 0, legend (t); end
    xlabel(xname);
    ylabel(yname);
    zlabel(zname);
    axis square;
    box;
    %legend('boxoff')

end

end