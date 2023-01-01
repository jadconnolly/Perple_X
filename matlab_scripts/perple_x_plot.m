
% MatLab demo script to plot Perple_X tab and ctr format files.
% JADC March 12, 2011

%clear all;

%clf(fig1);

% c(1) = 'r'
% c(2) = 'g'
% c(3) = 'b'
% c(4) = 'c'
% c(5) = 'm'
% c(6) = 'y'
% c(7) = 'k'
%
% for i = 1:7:70
%     i
%     for j = 1:1:7
%         j
%         i + j - 1
%         x(i+j-1) = char (c(j),'-')
%     end
% end

%set (gca, 'PlotBoxAspectRatio', [1 0.5 1],...
%    'FontSize',14.,...
%    'LineWidth',1.0)

LineWidth = 1.0
FontSize = 14.0
LineSpec = '-'
Marker = 'none'

[x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,titl] = function_to_get_perple_x_ss_file; %open the Perple_X file

function_for_perple_x_plots (x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,LineSpec,LineWidth,Marker,FontSize,titl);

