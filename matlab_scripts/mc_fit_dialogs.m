function varargout = mc_fit_dialogs(key, varargin)
% MC_FIT_DIALOGS Centralized dialogs for mc_fit plotting.
%
% Dialog families
%   Filter by data fit
%       ask_fit_filter_single
%           Appears in single-file workflows when a fit-based filter is
%           requested for one file only. In mc_fit_plot the perturbation
%           models are asked first; if declined, the central models are then
%           asked separately when their file is opened.
%
%   Filter by misfit value
%       ask_manual_misfit_filter
%           Appears when a manual misfit-value filter is available.
%       ask_manual_misfit_limits
%           Appears immediately after acceptance of the manual misfit-value
%           filter and collects the lower/upper ln(Misfit) or ln(Bayes) limits.
%
% Other dialogs
%       ask_apply_xyz_limits
%       ask_xyz_limits
%       ask_xy_limits
%       select_plot_file
%       select_axis_variable
%       ask_small_sample_continue
%
% Model type strings (central / perturbation) are embedded in dialog text.

switch lower(key)
    case 'ask_fit_filter_single'
        nfit = varargin{1};
        ntrial = varargin{2};
        model_text = varargin{3};

        ttl = 'Filter by data fit';
        msg = sprintf([ ...
            '%d out of %d %s models fit all analytical data within uncertainty. ' ...
            'Exclude the remaining %s models?' newline, newline,...
            'Answer NO if the analytical uncertainties are not well constrained.'], ...
            nfit, ntrial, model_text, model_text);
        choice = questdlg(msg, ttl, 'Yes', 'No', 'No');
        varargout{1} = strcmp(choice, 'Yes');

    case 'ask_manual_misfit_filter'
        model_text = varargin{1};
        ttl = 'Filter by misfit value';
        msg = sprintf('Apply a manual misfit-value filter to %s models?', model_text);
        choice = questdlg(msg, ttl, 'Yes', 'No', 'No');
        varargout{1} = strcmp(choice, 'Yes');

    case 'ask_manual_misfit_limits'
        model_text = varargin{1};
        zlabel = char(varargin{2});
        zmin = varargin{3};
        zmax = varargin{4};

        ttl = 'Filter by misfit value';
        options.Resize = 'on';
        prompt = { ...
            sprintf('Minimum %s value for %s models:', zlabel, model_text), ...
            sprintf('Maximum %s value for %s models:', zlabel, model_text)};
        def = {num2str(zmin), num2str(zmax)};
        answer = inputdlg(prompt, ttl, 1, def, options);
        if isempty(answer)
            varargout{1} = [];
            varargout{2} = [];
        else
            varargout{1} = str2double(answer{1});
            varargout{2} = str2double(answer{2});
        end

    case 'ask_apply_xyz_limits'
        ttl = 'X-Y-Z limits';
        if nargin > 1 && ~isempty(varargin{1})
            msg = char(varargin{1});
        else
            msg = 'Set the same X-Y-Z limits for all plots?';
        end
        choice = questdlg(msg, ttl, 'Yes', 'No', 'Yes');
        varargout{1} = strcmp(choice, 'Yes');

    case 'ask_xyz_limits'
        xr0 = varargin{1};
        yr0 = varargin{2};
        zr0 = varargin{3};
        options.Resize = 'on';
        answer = inputdlg({'Lower X-axis limit:','Upper X-axis limit:'}, ...
            'X limits', 1, {num2str(xr0(1)), num2str(xr0(2))}, options);
        if isempty(answer), varargout{1}=[]; varargout{2}=[]; varargout{3}=[]; return, end
        xr = [str2double(answer{1}), str2double(answer{2})];

        answer = inputdlg({'Lower Y-axis limit:','Upper Y-axis limit:'}, ...
            'Y limits', 1, {num2str(yr0(1)), num2str(yr0(2))}, options);
        if isempty(answer), varargout{1}=[]; varargout{2}=[]; varargout{3}=[]; return, end
        yr = [str2double(answer{1}), str2double(answer{2})];

        answer = inputdlg({'Lower Z-axis limit:','Upper Z-axis limit:'}, ...
            'Z limits', 1, {num2str(zr0(1)), num2str(zr0(2))}, options);
        if isempty(answer), varargout{1}=[]; varargout{2}=[]; varargout{3}=[]; return, end
        zr = [str2double(answer{1}), str2double(answer{2})];

        varargout{1} = xr;
        varargout{2} = yr;
        varargout{3} = zr;

    case 'ask_xy_limits'
        xr0 = varargin{1};
        yr0 = varargin{2};
        options.Resize = 'on';
        answer = inputdlg({'Lower X-axis limit:','Upper X-axis limit:'}, ...
            'X limits', 1, {num2str(xr0(1)), num2str(xr0(2))}, options);
        if isempty(answer), varargout{1}=[]; varargout{2}=[]; return, end
        xr = [str2double(answer{1}), str2double(answer{2})];

        answer = inputdlg({'Lower Y-axis limit:','Upper Y-axis limit:'}, ...
            'Y limits', 1, {num2str(yr0(1)), num2str(yr0(2))}, options);
        if isempty(answer), varargout{1}=[]; varargout{2}=[]; return, end
        yr = [str2double(answer{1}), str2double(answer{2})];

        varargout{1} = xr;
        varargout{2} = yr;

    case 'select_plot_file'
        type = varargin{1};
        [file, path, indx] = uigetfile({type,' ';'*.*','All Files (*.*)'}, ...
            ['Select an mc_fit ',type,' file']);
        varargout{1} = file;
        varargout{2} = path;
        varargout{3} = indx;

    case 'select_axis_variable'
        which_axis = varargin{1};
        dnames = varargin{2};
        prompt = sprintf('Select the %s-axis variable:', upper(which_axis));
        [pick, ok] = listdlg('PromptString',prompt, ...
            'ListSize',[200 400], 'SelectionMode','single', 'ListString',dnames);
        varargout{1} = pick;
        varargout{2} = ok;

    case 'ask_small_sample_continue'
        msg = varargin{1};
        choice = questdlg(msg, 'Data scarcity warning', 'Continue', 'Abort', 'Abort');
        varargout{1} = strcmp(choice, 'Continue');

    otherwise
        error('Unknown dialog key: %s', key)
end
end
