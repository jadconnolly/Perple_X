function [x,y,z,symb,xname,yname,zname,nvar,mvar,nrow,dnames,titl,type] ...
    = function_to_get_perple_x_ss_file

% MatLab script to read Perple_X tab files see:
%    perplex.ethz.ch/faq/Perple_X_tab_file_format.txt
% for format details.

% JADC March 26, 2011

% modified to allow arbitrary plotting of 3d data, requires spread sheet
% format. Old version is function_to_get_perple_x_file.

% JADC April 27, 2013

% modified to read *.pts files and return a type flag:
%    type = 1 - 2d tab file
%    type = 2 - 1d tab file
%    type = 3 - 2d pts file

% JADC February 24, 2025

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
symb = 0;

while ok == 0

    [filename, pathname, indx] = uigetfile({'*.tab;*.pts',' ';'*.*','All Files (*.*)'}, ...
        'Select a Perple_X tab or pts file');


    if indx == 0, errordlg('You did not choose a tab or pts file, I quit!'), break, end

    data_file=fullfile(pathname, filename);

    fid = fopen(data_file, 'rt');

    fmt = fgetl(fid); % read revision tag

    if strcmp(fmt,'|6.6.6')     % valid tab file revision

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

            type = 1;

            if mvar <= 2 || strtrim(dnames{1}{1}) ~= strtrim(vname(1,:))
                errordlg('The input data is not in spreadsheet format, set spreadsheet keyword to T and rerun WERAMI'); %not in ss format, write error
            end


            [xvar, ok] = listdlg('PromptString','Select the X-axis variable:','ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
            if ok == 0, errordlg('You did not choose a variable, I quit!'), break, end
            [yvar, ok] = listdlg('PromptString','Select the Y-axis variable:','ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
            if ok == 0, errordlg('You did not choose a variable, I quit!'), break, end
            [zvar, ok] = listdlg('PromptString','Select the Z-axis variable:','ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
            if ok == 0, errordlg('You did not choose a variable, I quit!'), break, end

            xname = dnames{1}{xvar};
            yname = dnames{1}{yvar};
            zname = dnames{1}{zvar};

            a = reshape(a,mvar,inc(1),inc(2));

            x = reshape(a(xvar,1:inc(1),1:inc(2)),inc(1),inc(2));
            x = rot90(x);
            x = flipud(x);
            y = reshape(a(yvar,1:inc(1),1:inc(2)),inc(1),inc(2));
            y = rot90(y);
            y = flipud(y);
            z = reshape(a(zvar,1:inc(1),1:inc(2)),inc(1),inc(2));
            z = rot90(z);
            z = flipud(z);

        elseif nvar == 1 % 1d table

            type = 2;

            [m , n] = size(a); nrow = m*n/mvar;
            z = reshape(a,mvar,nrow);
            yname = xname; % assign values to unused output arguments
            zname = xname;
            y = x;

        else

            errordlg(['The input data is ',nvar,'-dimensional, this script is only configured for 1-2d']);

        end

    else % *.pts file

        type = 3;

        T = readtable(data_file,'FileType','delimitedtext');
        ans = T.Properties.DimensionNames;
        A = table2array(T);
        ans = size (A);
        nrow = ans(1);
        cols = ans(2);
        %                                    name the variables
        dnames{1}{1} = 'symbol code';
        for i = 2:cols-1;dnames{1}{i} = ['parameter_',num2str(i-1)];end
        dnames{1}{cols} = 'Score';

        %                                    variable name hack
        dnames{1}{3} = '\itT\rm, K'
        dnames{1}{2} = '\itP\rm, bar'
        %                                    select the variables
        [xvar, ok] = listdlg('PromptString','Select the X-axis variable:','ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
        if ok == 0, errordlg('You did not choose a variable, I quit!'), end
        [yvar, ok] = listdlg('PromptString','Select the Y-axis variable:','ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
        if ok == 0, errordlg('You did not choose a variable, I quit!'), end
        [zvar, ok] = listdlg('PromptString','Select the Z-axis variable:','ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
        if ok == 0, errordlg('You did not choose a variable, I quit!'), end

        xname = dnames{1}{xvar};
        yname = dnames{1}{yvar};
        zname = dnames{1}{zvar};

        symb = A(:,1);
        x = A(:,xvar);
        y = A(:,yvar);
        z = A(:,zvar);
        nvar = 3;

        mvar = 3;

        titl  = 'point file (*.pts) plot'; % title

    % else
    % 
    %     choice = questdlg('Invalid file format, try again?','File error','Yes','No','Yes');
    % 
    %     switch choice
    %         case 'Yes'
    %             ok = 0;
    % 
    %         case 'No'
    %             ok = 1;
    %             break;
    %     end

    end % end if

end % end while

if mvar > 1

    choice = questdlg('Restrict range of data to be plotted?','Filter','Yes','No','No');

    switch choice
        case 'Yes'
            ok = 0;

        case 'No'
            ok = 1;
    end

         i = find(z(:)<0);
         z(i) = NaN;
         i = find(z(:)>1);
         z(i) = NaN;

    if ok == 0

        % filter the x-values
        dlg_title = 'Filter dialog...';
        prompt = {['Minimum ',xname,' value to be plotted:'],['Maximum ',xname,' value to be plotted:']};
        lines = 1;options.Resize='on';
        def = {num2str(min(x(:))), num2str(max(x(:))) };
        answer = inputdlg(prompt,dlg_title,lines,def,options);
        i = find(x(:)<str2num(answer{1}));
        x(i) = NaN;
        i = find(x(:)>str2num(answer{2}));
        x(i) = NaN;

        % filter the y-values
        dlg_title = 'Filter dialog...';
        prompt = {['Minimum ',yname,' value to be plotted:'],['Maximum ',yname,' value to be plotted:']};
        lines = 1;options.Resize='on';
        def = {num2str(min(y(:))), num2str(max(y(:))) };
        answer = inputdlg(prompt,dlg_title,lines,def,options);
        i = find(y(:)<str2num(answer{1}));
        y(i) = NaN;
        i = find(y(:)>str2num(answer{2}));
        y(i) = NaN;

        % filter the z-values
        dlg_title = 'Filter dialog...';
        prompt = {['Minimum ',zname,' value to be plotted:'],['Maximum ',zname,' value to be plotted:']};
        lines = 1;options.Resize='on';
        def = {num2str(min(z(:))), num2str(max(z(:))) };
        answer = inputdlg(prompt,dlg_title,lines,def,options);
        i = find(z(:)<str2num(answer{1}));
        z(i) = NaN;
        i = find(z(:)>str2num(answer{2}));
        z(i) = NaN;

    end

end