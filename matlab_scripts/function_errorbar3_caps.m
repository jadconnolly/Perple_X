function h = function_errorbar3_caps(ax,x,y,z,sigx,sigy,epsz,varargin)
%ERRORBAR3_CAPS Draw 3D error bars with caps.
%
%   h = ERRORBAR3_CAPS(ax,x,y,z,sigx,sigy,epsz)
%   h = ERRORBAR3_CAPS(ax,x,y,z,sigx,sigy,epsz,Name,Value,...)
%
% Draws 3D error bars centered at (x,y,z), with half-width errors:
%   x-direction: +/- sigx
%   y-direction: +/- sigy
%   z-direction: +/- epsz
%
% Caps are drawn with lengths equal to 10% of the orthogonal errors:
%   x-bar caps lie in the yz-plane:
%       cap along y = 0.1*sigy
%       cap along z = 0.1*epsz
%   y-bar caps lie in the xz-plane:
%       cap along x = 0.1*sigx
%       cap along z = 0.1*epsz
%   z-bar caps lie in the xy-plane:
%       cap along x = 0.1*sigx
%       cap along y = 0.1*sigy
%
% INPUTS
%   ax    - axes handle; pass [] to use gca
%   x,y,z - coordinates (scalars or equal-sized arrays)
%   sigx  - x error half-widths
%   sigy  - y error half-widths
%   epsz  - z error half-widths
%
% NAME-VALUE OPTIONS
%   'Color'     - line color, default [0 0 0]
%   'LineWidth' - line width, default 0.75
%   'LineStyle' - line style, default '-'
%   'cap_fraction' - cap size relative to orthogonal error, default 0.1
%
% OUTPUT
%   h is a struct with fields:
%       h.xbar, h.ybar, h.zbar
%       h.xcap_y, h.xcap_z
%       h.ycap_x, h.ycap_z
%       h.zcap_x, h.zcap_y
%
% If epsz(i) = 0 for a point, no z-direction error bar or z-related caps
% are drawn for that point.
%
% EXAMPLE
%   figure
%   ax = axes;
%   hold(ax,'on')
%   grid(ax,'on')
%   view(ax,3)
%   plot3(ax,1,2,3,'ko','MarkerFaceColor','k')
%   h = function_errorbar3_caps(ax,1,2,3,0.2,0.3,0.4,'Color','b','LineWidth',1.2);

    if nargin < 1 || isempty(ax)
        ax = gca;
    end

    if ~ishandle(ax) || ~strcmp(get(ax,'Type'),'axes')
        error('First argument must be an axes handle or [].');
    end

    % Parse name-value options
    p = inputParser;
    p.FunctionName = 'function_errorbar3_caps';
    addParameter(p,'Color',[0 0 0]);
    addParameter(p,'LineWidth',0.75,@(v) isnumeric(v) && isscalar(v) && v > 0);
    addParameter(p,'LineStyle','-',@(s) ischar(s) || (isstring(s) && isscalar(s)));
    addParameter(p,'cap_fraction',0.1,@(v) isnumeric(v) && isscalar(v) && v >= 0);
    parse(p,varargin{:});

    col = p.Results.Color;
    lw  = p.Results.LineWidth;
    ls  = char(p.Results.LineStyle);
    cap_fraction = p.Results.cap_fraction;

    % Force column vectors and check sizes
    x    = x(:);
    y    = y(:);
    z    = z(:);
    sigx = sigx(:);
    sigy = sigy(:);
    epsz = epsz(:);

    n = numel(x);

    if numel(y)    ~= n || numel(z)    ~= n || ...
       numel(sigx) ~= n || numel(sigy) ~= n || numel(epsz) ~= n
        error('x, y, z, sigx, sigy, and epsz must all have the same number of elements.');
    end

    if any(sigx < 0) || any(sigy < 0) || any(epsz < 0)
        error('sigx, sigy, and epsz must be non-negative.');
    end

    lineopts = {'Color',col,'LineWidth',lw,'LineStyle',ls,'HandleVisibility','off'};
    lineoptz = {'Color',col,'LineWidth',lw,'LineStyle',ls,'HandleVisibility','off'};

    hold_state = ishold(ax);
    hold(ax,'on')

    % Preallocate graphics handles
    h.xbar   = gobjects(n,1);
    h.ybar   = gobjects(n,1);
    h.zbar   = gobjects(n,1);
    h.xcap_y = gobjects(2*n,1);
    h.xcap_z = gobjects(2*n,1);
    h.ycap_x = gobjects(2*n,1);
    h.ycap_z = gobjects(2*n,1);
    h.zcap_x = gobjects(2*n,1);
    h.zcap_y = gobjects(2*n,1);

    for i = 1:n
        xi = x(i);
        yi = y(i);
        zi = z(i);

        sx = sigx(i);
        sy = sigy(i);
        ez = epsz(i);

        % Cap half-lengths = cap_fraction of orthogonal errors
        cx = cap_fraction * sx;
        cy = cap_fraction * sy;
        cz = cap_fraction * ez;

        % Main bars
        h.xbar(i) = line(ax,[xi-sx, xi+sx],[yi, yi],[zi, zi], lineopts{:});
        h.ybar(i) = line(ax,[xi, xi],[yi-sy, yi+sy],[zi, zi], lineopts{:});

        if ez > 0
            h.zbar(i) = line(ax,[xi, xi],[yi, yi],[zi-ez, zi+ez], lineoptz{:});
        end

        % X-bar caps at both ends:
        % caps perpendicular to x, drawn along y and z
        h.xcap_y(2*i-1) = line(ax,[xi-sx, xi-sx],[yi-cy, yi+cy],[zi, zi], lineopts{:});
        h.xcap_y(2*i)   = line(ax,[xi+sx, xi+sx],[yi-cy, yi+cy],[zi, zi], lineopts{:});

        if ez > 0
            h.xcap_z(2*i-1) = line(ax,[xi-sx, xi-sx],[yi, yi],[zi-cz, zi+cz], lineopts{:});
            h.xcap_z(2*i)   = line(ax,[xi+sx, xi+sx],[yi, yi],[zi-cz, zi+cz], lineopts{:});
        end

        % Y-bar caps at both ends:
        % caps perpendicular to y, drawn along x and z
        h.ycap_x(2*i-1) = line(ax,[xi-cx, xi+cx],[yi-sy, yi-sy],[zi, zi], lineopts{:});
        h.ycap_x(2*i)   = line(ax,[xi-cx, xi+cx],[yi+sy, yi+sy],[zi, zi], lineopts{:});

        if ez > 0
            h.ycap_z(2*i-1) = line(ax,[xi, xi],[yi-sy, yi-sy],[zi-cz, zi+cz], lineopts{:});
            h.ycap_z(2*i)   = line(ax,[xi, xi],[yi+sy, yi+sy],[zi-cz, zi+cz], lineopts{:});
        end

        % Z-bar caps at both ends:
        % caps perpendicular to z, drawn along x and y
        if ez > 0
            h.zcap_x(2*i-1) = line(ax,[xi-cx, xi+cx],[yi, yi],[zi-ez, zi-ez], lineoptz{:});
            h.zcap_x(2*i)   = line(ax,[xi-cx, xi+cx],[yi, yi],[zi+ez, zi+ez], lineoptz{:});

            h.zcap_y(2*i-1) = line(ax,[xi, xi],[yi-cy, yi+cy],[zi-ez, zi-ez], lineoptz{:});
            h.zcap_y(2*i)   = line(ax,[xi, xi],[yi-cy, yi+cy],[zi+ez, zi+ez], lineoptz{:});
        end
    end

    if ~hold_state
        hold(ax,'off')
    end
end