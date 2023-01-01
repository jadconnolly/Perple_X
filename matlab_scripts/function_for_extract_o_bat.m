function [] = function_for_extract_o_bat (x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,titl)
%Generic function to make 2- and 3-d plots from Perple_X tab format files
%

figure(1);

if nvar == 1 % two cases: 1d - table -> 2d plot
    
    errordlg(['oy Nimrod, this function is for 2d tables, I quit!'])
    
elseif nvar == 2 % 2d - table -> 2/3d plot
    
    amin = min(a(:)); amax = max(a(:));
    prompt = {[zname,' range is ',num2str(amin),' -> ',num2str(amax),'. Specify the contour value:']};dlg = 'Contour specification';
    num_lines = 1; def = {num2str((amin+amax)/2)};
    c = inputdlg(prompt,dlg,num_lines,def);
    conts = str2num(c{1});
    % for some reason the explicit vector of contour levels cannot
    % be of size 1 (because contourc interprets it as the number of levels)
    hold on
    C = contourc(x,y,a,[conts 1e99]);
    plot(C(1,2:C(2,1)),C(2,2:C(2,1)))
    xlabel(xname);ylabel(yname)
    
    [s1 s2] = size(C);

    if C(2,1)+1 <  s2,   
        helpdlg(['WARNING! the contour is not continuous'])
    end 
    
    
    filename = ['contour=',c{1},'.txt'];
    helpdlg(['output will be written to file: ',filename])
    fid = fopen(filename, 'wt');
    fprintf(fid, '%g %g\n', C(3:s1*s2))
    fclose(fid)
    

    
end

