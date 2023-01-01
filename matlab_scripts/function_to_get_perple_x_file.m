function [x,y,a,xname,yname,zname,nvar,mvar,nrow,dnames,titl] = function_to_get_perple_x_file

% MatLab script to read Perple_X tab files see:
%    perplex.ethz.ch/faq/Perple_X_tab_file_format.txt
% for format details.

% JADC March 26, 2011

% if nvar = 2: a 2d table at evenly spaced increments of x & y the 
% numbers of x-nodes and y-nodes are inc(1) and inc(2) and on return
% a(inc(2),inc(1)) is an arrary containing the value of a dependent
% property selected in this function by the user

% if nvar = 1: assumes a 1d table with nrow arbitrarily spaced rows
% consisting of mvar (>1) properties. nrow is computed by dividing the size
% of the data array a by mvar. This format requires reduntant information
% if the table is regularly spaced.

% see get_perple_x_file_with_regular_1d_grid for a mode efficient data
% format.

ok = 0;

while ok == 0;
    
    [filename, pathname] = uigetfile('*.tab', 'Select a Perple_X tab file');

    data_file=fullfile(pathname, filename);

    fid = fopen(data_file, 'rt');
    
    fmt = fgetl(fid); % read revision tag
    
    if strcmp(fmt,'|6.6.6')     % valid revision
        
        ok = 1;
        
        titl  = fgetl(fid); % title
        nvar  = fscanf(fid, '%f', 1); % number of independent variables
        
        for i = 1:nvar                % independent variables
                                      % the problem with using fscanf here
                                      % is that an error will result if the
                                      % variable name exceeds the length
                                      % specified by the format, there's
                                      % gotta be a better way to do this.
            vname(i,:) = fscanf(fid, '%9c', 1);
            vmin(i)    = fscanf(fid, '%f', 1);
            dv(i)      = fscanf(fid, '%f', 1);
            inc(i)     = fscanf(fid, '%f', 1);
            v(i,1:inc(i))   = vmin(i):dv(i):vmin(i)+(inc(i)-1)*dv(i);
        end
        
        mvar  = fscanf(fid, '%f', 1); % number of dependent variables
        dnames = textscan(fid,'%s',mvar); % dependent variable names
        
        fclose(fid);
        
        a = textread(data_file, '%f','headerlines',4*(nvar+1)+1); % read the numeric data
        
        x = v(1,1:inc(1));  % primary variable values, use also for size in 1d
        xname = vname(1,:); % primary variable name
        nrow = inc(1);      % number of rows
        
        if nvar == 2 % 2d table
            
            if mvar == 1
                
                dvar = mvar;
                
            else
                
                [dvar, ok] = listdlg('PromptString','Select the dependent variable:','ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
                if ok == 0, errordlg(['You did not choose a variable, I quit!']), end;
                
            end
            
            zname = dnames{1}{dvar};
            
            a = reshape(a,mvar,inc(1),inc(2));
            a = reshape(a(dvar,1:inc(1),1:inc(2)),inc(1),inc(2));
            a = rot90(a); 
            a = flipud(a);
                        % filter the a-values
            dlg_title = ['Filter dialog...'];
            prompt = {['Minimum ',zname,' value to be plotted:'],['Maximum ',zname,' value to be plotted:']}; 
            lines = 1;options.Resize='on';
            def = {num2str(min(a(:))), num2str(max(a(:))) };
            answer = inputdlg(prompt,dlg_title,lines,def,options); 
            i = find(a(:)<str2num(answer{1}));
            a(i) = NaN;
            i = find(a(:)>str2num(answer{2}));
            a(i) = NaN;
          
            
            y = v(2,1:inc(2));
            yname = vname(2,:);
            
        elseif nvar == 1 % 1d table
            
            [m n] = size(a); nrow = m*n/mvar;
            a = reshape(a,mvar,nrow);
            yname = xname; % assign values to unused output arguments
            zname = xname;
            y = x;
            
        else
            
            errordlg(['The input data is ',nvar,'-dimensional, this script is only configured for 1-2d']);
            
        end
        
    else
        
        choice = questdlg('Invalid file format, try again?','File error','Yes','No','Yes');
        
        switch choice;
            case 'Yes';
                ok = 0;
                
            case 'No';
                ok = 1;
                break;
        end
        
    end % end while
    
end

end