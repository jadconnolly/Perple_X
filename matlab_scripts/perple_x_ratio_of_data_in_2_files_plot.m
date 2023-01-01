
% MatLab script to plot the ratio of properties in two Perple_X 2d tab format files.
% JADC March 30, 2011

clear all; clf;

LineWidth = 1.0
FontSize = 14.0
LineSpec = '-'
Marker = 'none'
%                                 NUMERATOR
helpdlg('The next prompts select the file and variable to be used for the numerator of the ratio','Numerator Specification');
[x,y,a,xname,yname,z1name,nvar,mvar,nrow,dnames,titl] = function_to_get_perple_x_file; %open the Perple_X file
if nvar == 1 ,errordlg(['This script is only for 3d data, I quit!']), end;
%                                 DENOMINATOR
helpdlg('The next prompts select the file and variable to be used for the denominator of the ratio','Denominator Specification');
[x,y,b,xname,yname,z2name,nvar,mvar,nrow,dnames,titl] = function_to_get_perple_x_file; %open the Perple_X file
if nvar == 1 ,errordlg(['This script is only for 3d data, I quit!']), end;
%                                 make the ratio, could check sizes first. 
a = a./b;
zname = [z1name,'/',z2name];

function_for_perple_x_plots (x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,LineSpec,LineWidth,Marker,FontSize,titl);





