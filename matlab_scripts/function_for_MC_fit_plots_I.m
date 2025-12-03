function [] = function_for_MC_fit_plots_I (x,y,a,symb,fit,fitonly, ...
    xname,yname,zname,sname, ...
    nrow,xvar,yvar,zvar,file,path, ...
    LineStyle,LineWidth,Marker,FontSize)
%
% manipulate and plot the data in *_central.pts
%
%                                               JADC, 8/2025

figure(1);
%                                                symbol list
sc(1:12) = ['b','r','b','r','b','r','b','r','b','r','b','r'];
mt(1:12) = ['o','o','s','s','v','v','d','d','s','s','v','v'];
%                                                 format specification for legend text:
fmt = '%.2e';

%                                                 first figure
hold on

maxx = 1d99;
minx = -1d99;
maxy = maxx;
miny = minx;

icov = 0;

jct(1:4) = 0;
%                               check if it's misfit or bayes:
bayes = 0;
if size(find(symb==5),1) ~= 0 || size(find(symb==6),1) ~= 0, bayes = 3; end

if bayes == 3
    zname = "Bayes";
    sname(3) = zname;
end
%                               find the central model results
idx = find(symb==3+bayes);

if ~isempty(idx)
%                               must be doing a central pts file
    if isempty(find(symb==1, 1))
%                               but there's no best central model
%                               because mc_fit was terminated early
%                               make one if central models are present
%                               so make one
        if bayes == 3, [aopt, k] = max(a(idx)); else [aopt, k] = min(a(idx)); end
        symb(k) = 1;

    end

end
%                               group points by symbol
indx(1:5,1:nrow) = 0;
%                               best central model symbol = 1
j = find(symb==1);
jct(1) = size(j,1);
indx(1,1:size(j,1)) = j;
if bayes == 3, aopt = max(a(j)); else aopt = min(a(j)); end
%                               perturbed model symbol = 2 if Misfit
%                                                      = 5 if Bayes
j = find(symb==2+bayes);
jct(2) = size(j,1);
indx(2,1:size(j,1)) = j;
%                               central model symbol = 3 if misfit
%                                                    = 6 if Bayes
j = find(symb==3+bayes);
jct(3) = size(j,1);
indx(3,1:size(j,1)) = j;
%                               this code should be general with respect to
%                               central vs perturbed *.pts files except for
%                               the following test, delete either or both
%                               to make a generic *.pts plotting program
if jct(2) > 0 && jct(3) > 0
    errordlg ('symbol types 2|5 & 3|6 cannot both be present in a pts file')
    error ('I quit')
else if jct(2) > 0
    errordlg ('You did not select a *_central.pts file')
    error ('You did not select a *_central.pts file, I quit!')
end
%                               mc_fit saves the min/max parameter
%                               ranges as pts with symbol 9 & 10 use
%                               these to locate MPP coordinate
j = find(symb == 9);
k = find(symb == 10);
MPPX = (x(j) + x(k)) / 2; MPPY = (y(j) + y(k)) / 2; MPPZ = (a(j) + a(k)) / 2;
%                               find min/max x-y values, these will
%                               include the range min/max values from
%                               the above
j = find(symb~=0);
k = find(x(j) == min(x(j))); minx = x(k(1));
k = find(x(j) == max(x(j))); maxx = x(k(1));
k = find(y(j) == min(y(j))); miny = y(k(1));
k = find(y(j) == max(y(j))); maxy = y(k(1));
%                               to prevent the use min/max values in
%                               statistical analysis set them to NaNs
j = find(symb == 9 ); x(j) = NaN; y(j) = NaN; a(j) = NaN;
j = find(symb == 10); x(j) = NaN; y(j) = NaN; a(j) = NaN;

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

jplot = 0;

xy = subplot (2,5,[1,2,6,7]);

hold on
%                               coded so there will now always be
%                               a best central model, plot it:
ec = "#EDB120";

if fit(indx(1,1)) == 1, fc = ec; lw = 1; ec = 'k'; 
else, fc = 'none'; lw = 2; end

h0 = plot3 (xy,x(indx(1,1)),y(indx(1,1)),a(indx(1,1)),'o', ...
    'MarkerFaceColor',fc,...
    'MarkerEdgeColor',ec,...
    'MarkerSize',25, ...
    'LineWidth',lw);

strg = "\bf{Best central model}" + ...
    newline + xname + " = " + num2str(x(indx(1,1)),5) + ...
    newline + yname + " = " + num2str(y(indx(1,1)),5) + ...
    newline + zname + " = " + num2str(aopt,fmt);

t(1) = cellstr(strg);


if indx(2,1) > 0 || indx(3,1) > 0
    %                           assume only one symbol type (2|3)
    if indx(2,1) > 0, i = 2; else, i = 3; end
    %                           set the symbol code
    isymb = symb(indx(i,1));
    %                           save the indices
    ind = indx(i,1:jct(i));
    %                           set the edge color
    ec = sc(isymb);

    sz = 6;
    %                           make plot of all models that fit
    jnd = ind(find(fit(ind)==1));

    if ~isempty(jnd)

        h1a = plot3 (xy,x(jnd),y(jnd),a(jnd),[mt(isymb)], ...
            'MarkerFaceColor',ec,...
            'MarkerEdgeColor',ec,...
            'MarkerSize',sz);
    end

    %                           make plot of all models that don't fit
    jnd = ind(find(fit(ind)==0));

    if ~isempty(jnd)

        h1b = plot3 (xy,x(jnd),y(jnd),a(jnd),[mt(isymb)], ...
            'MarkerFaceColor','none',...
            'MarkerEdgeColor',ec,...
            'MarkerSize',sz);
    end

    if i == 3
        t(i) = cellstr('\bf{Central model tries}');
        jplot = i;
    else
        t(i) = cellstr('\bf{Perturbations}');
        jplot = i;
    end

end

if indx(4,1) > 0

    ec = "#0072BD";

    if fit(indx(4,1)) == 1, fc = ec; lw = 1; else, fc = 'none'; lw = 2; end

    h2 = plot3 (xy,x(indx(4,1)),y(indx(4,1)),a(indx(4,1)),'o', ...
        'MarkerFaceColor',fc,...
        'MarkerEdgeColor',ec,...
        'MarkerSize',25, ...
        'LineWidth',lw);

    if bayes == 3, bopt = max(a); else bopt = min(a); end

    strg = "\bf{Best overall}"     + newline + xname + " = " + ...
        num2str(x(indx(4,1)),5) + newline + yname + " = " + ...
        num2str(y(indx(4,1)),5) + newline + zname + " = " + ...
        num2str(bopt,fmt);

    t(4) = cellstr(strg);

end


view(xy,0,90) % XY

xlabel(xy,xname);
ylabel(xy,yname);
zlabel(xy,zname);
%zscale(xy,'log')
axis(xy,"square");
box(xy);
%                               to prevent the use min/max values in
%                               statistical analysis set them to NaNs
j = find(symb == 9 ); x(j) = NaN; y(j) = NaN; a(j) = NaN;
j = find(symb == 10); x(j) = NaN; y(j) = NaN; a(j) = NaN;

% choice = questdlg('Plot error ellipse?','Covariance','Yes','No','Yes');
%
% switch choice
%
%     case 'Yes'

hold on

[hcov,~,covariance] = function_to_plot_covariance_ellipse([x y],aopt,jplot);

strg = "\sigma_{" + sname(1) + "}\rm = \pm" + num2str(sqrt(covariance(1,1)),fmt) ...
    + newline + ...
    "\sigma_{" + sname(2) + "}\rm = \pm" + num2str(sqrt(covariance(2,2)),fmt);

t(5) = cellstr(strg);

icov = 5;

% end
%                               plot the MPP coordinate
hmpp = plot3 (xy,MPPX,MPPY,MPPZ,'d', ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',15);
%                               copy all plotted points from xy into xz and
%                               yz
%for ax = [xz yz]
%    copyobj(allchild(xy),ax);
%end
%                               text for legend
strg = "\bf{MPP}" + newline + xname + " = " + num2str(MPPX,5) + ...
    newline + yname + " = " + num2str(MPPY,5);

t(6) = cellstr(strg);

if jplot > 0   % allow for the eventuallity of no data

    if icov == 0   % no covariance
        if jct(4) == 1 &&  jplot > 0
            legend (xy,[h0,h1a,h2,hmpp],[t(1) t(jplot) t(4) t(6)]);
        else
            legend (xy,[t(1) t(jplot) t(6)]);
        end
    else
        if jct(4) == 1
            if exist('h1a','var')
            legend (xy,[h0 h1a h2 hcov hmpp],[t(1) t(jplot) t(4) t(icov) t(6)]);
            else 
            legend (xy,[h0 h1b h2 hcov hmpp],[t(1) t(jplot) t(4) t(icov) t(6)]);
            end
        else
            if exist('h1a','var')
            legend (xy,[h0 h1a hcov hmpp],[t(1) t(jplot) t(icov) t(6)]);
                        else 
            legend (xy,[h0 h1b hcov hmpp],[t(1) t(jplot) t(icov) t(6)]);
            end
        end
    end

end

legend ('Location','northwest')
title (xy,'Central Model Analysis')
xy.ZScale = 'log';

xz = subplot (2,5,3);

if exist('h0','var'), copyobj(h0, xz), end
if exist('h1a','var'), copyobj(h1a, xz), end
if exist('h1b','var'), copyobj(h1b, xz), end
if exist('h2','var'), copyobj(h2, xz), end

view(xz,0,0)

xlabel(xz,xname);
ylabel(xz,yname);
zlabel(xz,zname);
xz.ZScale = 'log';
axis(xz,"square");
box(xz);

yz = subplot (2,5,8);

if exist('h0','var'), copyobj(h0, yz), end
if exist('h1a','var'), copyobj(h1a, yz), end
if exist('h1b','var'), copyobj(h1b, yz), end
if exist('h2','var'), copyobj(h2, yz), end

view(yz,90,0)

xlabel(yz,xname);
ylabel(yz,yname);
zlabel(yz,zname);
yz.ZScale = 'log';
axis(yz,"square");
box(yz);
%                             open the perturbation result file:
first = 1;

[x,y,a,symb,fit,~,~,~,~,~,nrow,~,~,~,~,~] = function_to_get_MC_fit_plot_file ( ...
    '*_perturbed.pts',xvar,yvar,zvar, fitonly, xname, yname, zname, sname, file, path, first);

function_for_MC_fit_plots_II (x, y, a, symb, fit, xname, yname, zname, sname, nrow, ...
    LineStyle, LineWidth, Marker, FontSize, fmt, file, ...
    icov, xy, xz, yz, minx, maxx, miny, maxy, bayes);

end
