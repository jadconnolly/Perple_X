% load saved figures
fig1   = hgload('/path');
fig2   = hgload('/path');
%  Prepare 'subplot'
figure
h1=subplot(1,1,1);
h2=subplot(1,1,1);
copyobj(allchild(get(fig1,'CurrentAxes')),h1)
hold on
copyobj(allchild(get(fig2,'CurrentAxes')),h1)



Fe^{3+}/Fe^{2+} Plagioclase-rich domain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a point
plot3 ([],[],[],['o'], ...
    'MarkerFaceColor',"#EDB120",...
    'MarkerEdgeColor','k',...
    'MarkerSize',25);


%--------------------------------------------------------


fig1   = hgload('ag_assemblage.fig');
fig2   = hgload('ag_bulk.fig')

fig1   = hgload('ag_assemblage_simple.fig');
fig2   = hgload('ag_bulk_simple.fig');

figure (4)
h1=subplot(1,2,1);
h2=subplot(1,2,2);
copyobj(allchild(get(fig1,'CurrentAxes')),h1)
copyobj(allchild(get(fig2,'CurrentAxes')),h2)
xlim(h2,[500,1000])
box(h2)
legend(h1)
legend(h2)

figure (4)
h1=subplot(1,1,1);
h2=subplot(1,1,1);
copyobj(allchild(get(fig1,'CurrentAxes')),h1)
hold on
copyobj(allchild(get(fig2,'CurrentAxes')),h1)

legend
box
axis square
xlabel('\itT\rm, K');
ylabel('\itP\rm, bar');
zlabel('Score');

subplot

%----------------------------------------------------------
% 3 3d sectionscla plot

% make original figure handle = f0

figure(1)
h1 = subplot (2,4,[1,2,5,6])
h2 = subplot (2,4,[3,4])
h3 = subplot (2,4,[7,8])
copyobj(allchild(get(f0,'CurrentAxes')),h1)
copyobj(allchild(get(f0,'CurrentAxes')),h2)
copyobj(allchild(get(f0,'CurrentAxes')),h3)
view(h2,0,0)
view(h3,90,0)
legend(h1)
xlabel(h1,xname);
xlabel(h2,xname);
ylabel(h1,yname);
ylabel(h3,yname);
zlabel(h2,zname);
zlabel(h3,zname);
zscale(h2,"log");
zscale(h3,"log");
%vk30
%xr = [573,973]
%yr = [1e4,3e4]
%bl synthetic
xr = [800,1000]
yr = [1e3,15e3]
xlim(h1,xr);xlim(h2,xr);xlim(h3,xr)
ylim(h1,yr);ylim(h2,yr);ylim(h3,yr)
box(h1); box(h2); box(h3);
axis(h1,"square") %'auto','auto x', etc, 'tight','normal'






%----------------------------------------------------------
% 2 3d sections 

% make original figure handle = f0

f0 = figure(1)
figure(2)
h1 = subplot (1,2,1)
h2 = subplot (1,2,2)
copyobj(allchild(get(f0,'CurrentAxes')),h1)
copyobj(allchild(get(f0,'CurrentAxes')),h2)
view(h1,0,0)
view(h2,90,0)
%legend(h1)
xlabel(h1,xname);
ylabel(h2,yname);
zlabel(h1,zname);
zscale(h1,"linear");
zscale(h2,"linear");
xr = [800,1000]
yr = [5e3,15e3]
zr = [0.0,0.2]
xlim(h1,xr);xlim(h2,xr);
ylim(h1,yr);ylim(h2,yr);
zlim(h1,zr);zlim(h2,zr);
box(h1); box(h2); 
%axis(h1,"square") %'auto','auto x', etc, 'tight','normal'
a12 =subplot(122);
a12.Position
a12.Position(2) = 0.16; a12.Position(1) = 0.48;
a12.ZTickLabel = []

a11 =subplot(121);
a11.Position
a11.Position(2) = 0.16;
a11.XTick = [850,900,950]
%--------------------------------------------------------
fig1   = hgload('AB.fig');
fig2   = hgload('CD.fig')
fig3   = hgload('EF.fig');
fig4   = hgload('GH.fig');

figure (5)
h1=subplot(2,2,1);
h2=subplot(2,2,2);
h3=subplot(2,2,3);
h4=subplot(2,2,4);
copyobj(allchild(get(fig1,'CurrentAxes')),h1)
copyobj(allchild(get(fig2,'CurrentAxes')),h2)
copyobj(allchild(get(fig3,'CurrentAxes')),h3)
copyobj(allchild(get(fig4,'CurrentAxes')),h4)



xr = [750,1050]
yr = [3e3,11e3]
xlim(xr)
ylim(yr)
box;
%--------------------------------------------------------
% 3 error plots

clf
cla

fig1   = hgload('bl_total.fig');
fig2   = hgload('bl_thermo.fig')
fig3   = hgload ('bl_analytical.fig')

figure (1)
h1=subplot(1,3,1);
h2=subplot(1,3,2);
h3=subplot(1,3,3);
copyobj(allchild(get(fig1,'CurrentAxes')),h1)
copyobj(allchild(get(fig2,'CurrentAxes')),h2)
copyobj(allchild(get(fig3,'CurrentAxes')),h3)

%vk30
%xr = [573,973]
%yr = [1e4,3e4]
%bl synthetic
xr = [700,1100]
yr = [1e3,15e3]
xlim(h1,xr);xlim(h2,xr);xlim(h3,xr)
ylim(h1,yr);ylim(h2,yr);ylim(h3,yr)
axis (h1, 'square');axis (h2, 'square');axis (h3, 'square');
xlabel(h1,'\itT\rm, K');
xlabel(h2,'\itT\rm, K');
xlabel(h3,'\itT\rm, K');
ylabel(h1,'\itP\rm, bar');
box(h1);box(h2);box(h3)
legend(h1)
legend(h2)
legend(h3)
title (h1,'Total Uncertainty')
title (h2,'Thermodynamic Uncertainty')
title (h3,'Analytical Uncertainty')
% xt = 725;
% yt = 8000;
% text (h1,xt,yt,['\mu_{\itT\rm} = 1009 \pm 40 K',newline,'\mu_{\itP\rm} = 10.3 \pm 1.9 kbar'],FontSize=14)
% text (h2,xt,yt,['\mu_{\itT\rm} = 1014 \pm 33 K',newline,'\mu_{\itP\rm} = 10.7 \pm 0.9 kbar'],FontSize=14)
% text (h3,xt,yt,['\mu_{\itT\rm} = 1003 \pm 24 K',newline,'\mu_{\itP\rm} = 10464 \pm 1.0 kbar'],FontSize=14)
yticklabels(h3,[""])
yticklabels(h2,[""])
%--------------------------------------------------------
% 3 error plots

cla
clf

fig1   = hgload('vk30_matrix_total.fig');
fig2   = hgload('vk30_matrix_thermo.fig')
fig3   = hgload ('vk30_matrix_analytic.fig')

figure (1)
h1=subplot(1,3,1);
h2=subplot(1,3,2);
h3=subplot(1,3,3);
copyobj(allchild(get(fig1,'CurrentAxes')),h1)
copyobj(allchild(get(fig2,'CurrentAxes')),h2)
copyobj(allchild(get(fig3,'CurrentAxes')),h3)

%vk30
% xr = [573,973]
% yr = [1e4,3e4]
%bl synthetic
xr = [700,1100]
yr = [1e3,15e3]
xlim(h1,xr);xlim(h2,xr);xlim(h3,xr)
ylim(h1,yr);ylim(h2,yr);ylim(h3,yr)
axis (h1, 'square');axis (h2, 'square');axis (h3, 'square');
xlabel(h1,'\itT\rm, K');
xlabel(h2,'\itT\rm, K');
xlabel(h3,'\itT\rm, K');
ylabel(h1,'\itP\rm, bar');
box(h1);box(h2);box(h3)
legend(h1)
legend(h2)
legend(h3)
title (h1,'Total Uncertainty')
title (h2,'Thermodynamic Uncertainty')
title (h3,'Analytical Uncertainty')
% xt = 573;
% yt = 28000;
%text (h1,xt,yt,['\mu_{\itT\rm} = 771 \pm 27 K',newline,'\mu_{\itP\rm} = 16.0 \pm 2.4 kbar'],FontSize=14)
%text (h2,xt,yt,['\mu_{\itT\rm} = 777 \pm 26 K',newline,'\mu_{\itP\rm} = 16.8 \pm 0.8 kbar'],FontSize=14)
%text (h3,xt,yt,['\mu_{\itT\rm} = 807 \pm 14 K',newline,'\mu_{\itP\rm} = 17.3 \pm 1.0 kbar'],FontSize=14)
yticklabels(h3,[""])
yticklabels(h2,[""])
%-------------------------------------------------------------
                 strg = ''
                 strg = string(strg) + ['\sigma = \pm',num2str(sqrt(covariance(1,1)))] + ...
                      newline + ['\sigma = \pm',num2str(sqrt(covariance(2,2)))];

                    
                                 strg = string(t(jplot)) + newline + ...
                    ['\sigma =\pm',num2str(sqrt(covariance(1,1)))] + ...
                    newline + ['\pm ',num2str(sqrt(covariance(2,2)))];

                                 strg = string (['\sigma =\pm',num2str(sqrt(covariance(1,1)))])
                                 strg = string(strg) + newline ...
                                     + ['\sigma =\pm',num2str(sqrt(covariance(2,2)))];
%-------------------------------------------------------------
figure (4)
h1=subplot(1,1,1);
h2=subplot(1,1,1);
copyobj(allchild(get(fig1,'CurrentAxes')),h1)
hold on
copyobj(allchild(get(fig2,'CurrentAxes')),h1)

legend
box
axis square
xlabel('\itT\rm, K');
ylabel('\itP\rm, bar');
zlabel('Score');

subplot



ps = find(y < 3e4 & y > 1e4)

ipts = 0 
for i = 1, 699
    if x(ps(i)) <= 973 & x(ps(i)) >= 573, 
        ipts = ipts + 1
        inds(ipts) = ps(i);
    end
end

