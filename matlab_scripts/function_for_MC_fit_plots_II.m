function [] = function_for_MC_fit_plots_II ( ...
    x,y,a,symb,fit,xname,yname,zname,sname, ...
    nrow,LineStyle,LineWidth,Marker,FontSize,fmt, file, ...
    icov, xy, xz, yz, minx, maxx, miny, maxy, bayes)
%
% manipulate and plot the data in *_perturbed.pts
%
%                                               JADC, 8/2025

figure(1);

%                                                symbol list
sc(1:12) = ['b','r','b','r','b','r','b','r','b','r','b','r'];
mt(1:12) = ['o','o','s','s','v','v','o','o','s','s','v','v'];
%                                                 first figure
hold on

jct(1:4) = 0;
%                               group points by symbol
indx(1:5,1:nrow) = 0;
%                               best central model symbol = 1
j = find(symb==1);
jct(1) = size(j,1);
indx(1,1:size(j,1)) = j;
aopt = a(j(1));
%                               perturbed model symbol = 2 if Misfit
%                                                      = 5 if Bayes
j = find(symb==2+bayes);
jct(2) = size(j,1);
indx(2,1:size(j,1)) = j;
if jct(2) == 0
    errordlg('No successful perturbations to plot')
    error('The *_perturbed.pts contains no results, I quit')
end
%                               to prevent the use min/max values in
%                               statistical analysis set them to NaNs
j = find(symb == 9 ); x(j) = NaN; y(j) = NaN; a(j) = NaN;
j = find(symb == 10); x(j) = NaN; y(j) = NaN; a(j) = NaN;
%                               check for central model symbols

%                               find min/max x-y values
j = find(symb~=0);
k = find(x(j) == min(x(j)));
if x(k(1)) < minx, minx = min(x(j)); end
k = find(x(j) == max(x(j)));
if x(k(1)) > maxx, maxx = max(x(j)); end
k = find(y(j) == min(y(j)));
if y(k(1)) < miny, miny = min(y(j)); end
k = find(y(j) == max(y(j)));
if y(k(1)) > maxy, maxy = max(y(j)); end

jct(4) = 0;

if indx(2,1) > 0 && indx(1,1) > 0
    %                           if the best perturbed model is better than
    %                           the best central model, save it:
    if bayes == 0
        if min(a(indx(2,1:jct(2)))) < a(indx(1,1))
            jct(4) = 1;
            indx(4,1) = find(a(:) == min(a(indx(2,1:jct(2)))));
        end
    else
        if max(a(indx(2,1:jct(2)))) > a(indx(1,1))
            jct(4) = 1;
            indx(4,1) = find(a(:) == max(a(indx(2,1:jct(2)))));
        end
    end
end

hold on

xy2 = subplot (2,5,[4,5,9,10]);

hold on

if jct(2) > 0

    isymb = symb(indx(2,1));
    %                           set edge color
    ec = sc(isymb);
    sz = 6;
    %                           save the indices
    ind = indx(2,1:jct(2));
    %                           make plot of all models that fit
    jnd = ind(find(fit(ind)==1));

    if ~(isempty(jnd))
        h1a = plot3 (xy2,x(jnd),y(jnd),a(jnd),[mt(isymb)], ...
            'MarkerFaceColor',ec,...
            'MarkerEdgeColor',ec,...
            'MarkerSize',sz);
    end

    t(2) = cellstr('\bf{Perturbations}');

        %                           make plot of all models that don't fit
    jnd = ind(find(fit(ind)==0));

    if ~(isempty(jnd))    
       h1b = plot3 (xy2,x(jnd),y(jnd),a(jnd),[mt(isymb)], ...
            'MarkerFaceColor','none',...
            'MarkerEdgeColor',ec,...
            'MarkerSize',sz);
    end

end

if jct(4) > 0

    ec = "#0072BD";

    if fit(indx(4,1)) == 1, fc = ec; lw = 1; else, fc = 'none'; lw = 2; end

    h2 = plot3 (xy2,x(indx(4,1)),y(indx(4,1)),a(indx(4,1)),'o', ...
        'MarkerFaceColor',fc,...
        'MarkerEdgeColor',ec,...
        'MarkerSize',25,'LineWidth',lw);

    if bayes == 3, bopt = max(a); else bopt = min(a); end

    strg = "\bf{Best overall}"  + newline + xname + " = " + ...
        num2str(x(indx(4,1)),5) + newline + yname + " = " + ...
        num2str(y(indx(4,1)),5) + newline + zname + " = " + ...
        num2str(bopt,fmt);

    t(4) = cellstr(strg);

end

view(xy2,0,90) % XY

xlabel(xy2,xname);
ylabel(xy2,yname);
zlabel(xy2,zname);
%zscale(xy2,'log') 2024
xy2.ZScale = 'log';
axis(xy2,"square");
box(xy2);

if icov ~= 0

    hold on

    [hcov,~,covariance] = function_to_plot_covariance_ellipse([x y],aopt,2);

    strg = "\sigma_{" + sname(1) + "}\rm = \pm" + num2str(sqrt(covariance(1,1)),fmt) ...
        + newline + ...
        "\sigma_{" + sname(2) + "}\rm = \pm" + num2str(sqrt(covariance(2,2)),fmt) ...
        + newline + ...
        "\epsilon_{\rm" + sname(3) + "} = \pm" + num2str(std(a,'omitnan'),fmt);

    t(5) = cellstr(strg);

end

if jct(1) > 0

    ec = "#EDB120";

    if fit(indx(1,1)) == 1, fc = ec; ec = 'k'; lw = 1; else, fc = 'none'; lw = 2; end

    h0 = plot3 (xy2,x(indx(1,1)),y(indx(1,1)),a(indx(1,1)),'o', ...
        'MarkerFaceColor',fc,...
        'MarkerEdgeColor',ec,...
        'MarkerSize',25,'LineWidth',lw);

    strg = "\bf{Best central model}" + ...
        newline + xname + " = " + num2str(x(indx(1,1)),5) + ...
        newline + yname + " = " + num2str(y(indx(1,1)),5) + ...
        newline + zname + " = " + num2str(min(a),fmt);

    t(1) = cellstr(strg);

end

if icov == 0                     % no covariance, this ain't gonna happen
    if jct(4) == 1
        if exist('h1a','var')
        legend (xy2,[h1a h2],[t(2) t(4)]);
        else
        legend (xy2,[h1b h2],[t(2) t(4)]);
        end
    else
        if exist('h1a','var')
        legend (xy2,h1a,[t(2)]);
        else 
        legend (xy2,h1b,[t(2)]);
        end
    end
else
    if jct(4) == 1
        if exist('h1a','var')
        legend (xy2,[h1a h2 hcov],[t(2) t(4) t(5)]);
        else
        legend (xy2,[h1b h2 hcov],[t(2) t(4) t(5)]);
        end 
    else
        if exist('h1a','var')
        legend (xy2,[h1a hcov],[t(2) t(5)]);
        else
        legend (xy2,[h1b hcov],[t(2) t(5)]);
        end 
    end
end

legend(xy2,'Location','northwest'); 

title (xy2,'Perturbation Analysis')

choice = questdlg('Set the same X-Y limits for all plots?', ...
    'X-Y limits','Yes','No','Yes');

switch choice
    case 'Yes'
        ok = 0;

    case 'No'
        ok = 1;
end

if ok == 0

    prompt = {'Lower X-axis limit:','Upper X-axis limit:'};
    lines = 1; options.Resize='on';
    def = {num2str(minx), num2str(maxx)};
    answer = inputdlg(prompt,'X limits',lines,def,options);

    minx = str2num(answer{1});
    maxx = str2num(answer{2});

    prompt = {'Lower Y-axis limit:','Upper Y-axis limit:'};
    lines = 1; options.Resize='on';
    def = {num2str(miny), num2str(maxy)};
    answer = inputdlg(prompt,'Y limits',lines,def,options);

    miny = str2double(answer{1});
    maxy = str2double(answer{2});

    xr = [minx,maxx];
    yr = [miny,maxy];

    xlim(xy,xr);xlim(xz,xr);xlim(yz,xr);xlim(xy2,xr);
    ylim(xy,yr);ylim(xz,yr);ylim(yz,yr);ylim(xy2,yr);

end

% stripped file name for plot output

file = strrep (file,"_perturbed.pts","");

% to use my screen settings, set my_settings ~= 1

my_settings = 0;

if my_settings == 0

    choice = questdlg('Add a title?','Title','Yes','No','No');

    switch choice

        case 'Yes'

            answer = inputdlg('Enter the title:');
            name = string(answer{1});

        case 'No'

            name = "Click on axes to rearrange, rotate, etc.";

    end

    sgtitle(name)

    % msg = "Output plots (" + file + ".png, " + file + ".pdf)?";
    % choice = questdlg(msg,'Plot','Yes','No','Yes');
    % 
    % switch choice
    % 
    %     case 'Yes'

            exportgraphics(gcf,strcat(file,'.png'),'ContentType','image','Resolution',300);
            exportgraphics(gcf,strcat(file,'.pdf'),'ContentType','vector','BackgroundColor','none');

    % end

else

    %                             CUSTOMIZED SETTINGS:
    %                             reposition the xz/yz plots
    %                             this is one off and specifically
    %                             for my screen
    xz.Position(2) = 0.49;
    yz.Position(2) = 0.2;

    name = "Unfiltered bl_Chi problem results";

    % normally matlab uses sgtitle(name) but takes up too much space
    % so i use annotation here instead
    
    annotation('textbox',...
        [0.29475 0.792439862542954 0.42478125 0.0329896907216494],...
        'String',name,...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',15,...
        'FitBoxToText','on');

    exportgraphics(gcf,strcat(file,'.png'),'ContentType','image','Resolution',300);
    exportgraphics(gcf,strcat(file,'.pdf'),'ContentType','vector','BackgroundColor','none');


end

