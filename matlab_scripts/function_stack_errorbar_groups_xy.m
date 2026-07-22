function function_stack_errorbar_groups_xy(hblue, hred, sigx_blue, sigy_blue, sigx_red, sigy_red)
% Stack x and y error-bar components independently.
% If the perturbation (red) bar is smaller in a given direction, put it in front.
% Otherwise leave the imprecision (blue) bar in front.

if isempty(hblue) || isempty(hred)
    return
end

if sigx_red < sigx_blue
    local_stack_front(local_collect_group(hred, 'x'));
else
    local_stack_front(local_collect_group(hblue, 'x'));
end

if sigy_red < sigy_blue
    local_stack_front(local_collect_group(hred, 'y'));
else
    local_stack_front(local_collect_group(hblue, 'y'));
end
end

function h = local_collect_group(hs, which)
    h = gobjects(0,1);
    if isempty(hs)
        return
    end

    if isstruct(hs)
        switch lower(which)
            case 'x'
                h = local_concat_fields(hs, {'xbar','xcap_y','xcap_z'});
            case 'y'
                h = local_concat_fields(hs, {'ybar','ycap_x','ycap_z'});
        end
    else
        try
            hs = hs(:);
            h = hs(isgraphics(hs));
        catch
            h = gobjects(0,1);
        end
    end
end

function h = local_concat_fields(s, names)
    h = gobjects(0,1);
    for i = 1:numel(names)
        if isfield(s, names{i})
            hi = s.(names{i});
            try
                hi = hi(:);
                h = [h; hi(isgraphics(hi))]; %#ok<AGROW>
            catch
            end
        end
    end
end

function local_stack_front(h)
    if isempty(h)
        return
    end
    try
        uistack(h, 'top')
    catch
        % Ignore stacking failures on older MATLAB versions/graphics states.
    end
end
