% MatLab demo script to plot MC_fit *.pts format files, the script prompts
% for the *_central.pts and automatically reads the corresponding 
% *_peturbed.pts file.

% to customize the arrangement of the final subplots see my_settings
% near the end of function_for_MC_plots_II.m

% JADC Aug 27, 2025

clf

LineWidth = 1.0;
FontSize = 14.0;
LineSpec = '-';
Marker = 'none';
%                                       open the *_central.pts file
first = 0;xvar=0;yvar=0;zvar=0;fitonly=0;file="";path="";xname="";
yname="";zname="";sname="";

[x,y,a,symb,fit,fitonly,xname,yname,zname,sname,nrow,xvar,yvar,zvar, ...
    file,path] = function_to_get_MC_fit_plot_file ('*_central.pts',xvar, ...
    yvar, zvar, fitonly, xname, yname, zname, sname, file, path, first); 

file = strrep (file,"central","perturbed");

function_for_MC_fit_plots_I (x,y,a,symb,fit,fitonly,xname,yname,zname,sname,nrow, ...
                             xvar,yvar,zvar,file,path, ...
                             LineSpec,LineWidth,Marker,FontSize);

