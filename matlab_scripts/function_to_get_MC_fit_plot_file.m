function [x,y,z,symb,fit,fitonly,xname,yname,zname,sname,nrow,xvar,yvar,zvar,file,path] ...
    = function_to_get_MC_fit_plot_file (type,xvar,yvar,zvar,fitonly, ...
    xname,yname,zname,sname,file,path,first)

% function to read *.pts file for MC_fit_plot script, this script
% is essentially the same as function_to_get_perple_x_ss except that
% type specifies a filter value (*_central.pts or *_perturbed.pts)
% and the tab file capabilities has been cut.

if first == 0

    [file, path, indx] = uigetfile({type,' ';'*.*', ...
        'All Files (*.*)'}, ['Select an MC_fit ',type,' file']);

    if indx == 0, error(['You did not choose a',type,' pts file, I quit!']), end

end

data_file = fullfile(path, file);

fid = fopen (data_file, 'rt');

T = readtable (data_file,'FileType','delimitedtext','PreserveVariableNames',true);

fclose (fid);

colnames = string(T.Properties.VariableNames);

ans = T.Properties.DimensionNames;
A = table2array(T);
ans = size (A);
nrow = ans(1);
cols = ans(2);

% the MC_fit pts files use a "999" to indicate the end of the parameter
% data, after the 999 the additional values are misfit, bayes score,
% and bulk composition

cols = find(A(1,1:cols)==999) - 1;
%                                    check if misfit or bayes
i = find(A(:,1)==3);
bayes = 3;
if ~isempty(i), bayes = 0; end
if bayes == 0, zname = 'Misfit'; else, zname = 'Bayes'; end
%                                    variable name hack
dnames{1}{1} = '\itP\rm, bar';
dnames{1}{2} = '\itT\rm, K';
%                                    for the rest use perplex names
if cols > 5
    for i = 5:cols - 1
        dnames{1}{i-2} = colnames(i);
    end
end

if first == 0
    % select the variables
    ok = 0;

    [xvar, ok] = listdlg('PromptString','Select the X-axis variable:', ...
        'ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
    if ok == 0, error('You did not choose a variable, I quit!'), end
    [yvar, ok] = listdlg('PromptString','Select the Y-axis variable:', ...
        'ListSize',[200 400],'SelectionMode','single','ListString',dnames{1});
    if ok == 0, error('You did not choose a variable, I quit!'), end

    xname = dnames{1}{xvar};
    yname = dnames{1}{yvar};
    % make subscript names:
    if xvar == 2
        sname(1) = "\itT";
    elseif xvar == 1
        sname(1) = "\itP";
    else
        sname(1) = "\it p_{\rm" + string(xvar) + "}";
    end

    if yvar == 2
        sname(2) = "\itT";
    elseif yvar == 1
        sname(2) = "\itP";
    else
        sname(2) = "\it p_{\rm" + string(xvar) + "}";
    end

    xvar = xvar + 2;
    yvar = yvar + 2;
    zvar = cols;

    sname(3) = zname;

end

symb = A(:,1);
fit  = A(:,2);
x = A(:,xvar);
y = A(:,yvar);
z = A(:,zvar);

if first == 0

    text = 'the central';
    % check if the models are mixed,
    % some fit data, others don't
    all = ismember(symb, [2 3 5 6]); sall = sum(all);
    ins = ismember(symb, [2 3 5 6]) & fit == 1; sins = sum(ins);
    outs = ismember(symb, [2 3 5 6]) & fit == 0; souts = sum(outs);

    if sall ~= sins && sall ~= souts
        % mixed
        choice = questdlg([num2str(sins),' out of the ',num2str(sall), ...
            ' central model results fit all the analytical',newline, ...
            'data within its uncertainty. Remove the models that do not fit?'], ...
            'Fit or not','Yes','No','Yes');

        switch choice

            case 'Yes'

                fitonly = 1;
                idx = find(outs);
                x(idx) = NaN; y(idx) = NaN; z(idx) = NaN;

            case 'No'

                fitonly = 0;
        end

    end

    choice = questdlg(['Filter the ',text,' models?'],'Filter','Yes','No','No');

    switch choice

        case 'Yes'

            dlg_title = 'Filter dialog...'; options.Resize='on';
            prompt = {['Minimum ',zname,' score:'],['Maximum ',zname,' score:']};
            def = {num2str(min(z(:))), num2str(max(z(:))) };
            answer = inputdlg(prompt,dlg_title,1,def,options);

            lb = 0.999*str2double(answer{1});
            ub = 1.001*str2double(answer{2});

            idx = find(z(:) < lb | z(:) > ub);

            x(idx) = NaN; y(idx) = NaN; z(idx) = NaN;

    end

else

    if fitonly == 1
        idx = find(ismember(symb, [2 3 5 6]) & fit == 0);
        x(idx) = NaN; y(idx) = NaN; z(idx) = NaN;
    end

end

