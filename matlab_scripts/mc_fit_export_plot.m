function mc_fit_export_plot(fig_handle, file_root, outPDF)
% mc_fit_export_plot export final plot as either png bitmap or vector pdf.
% outPDF = false -> png (default)
% outPDF = true  -> pdf

if nargin < 3 || isempty(outPDF)
    outPDF = false;
end

if outPDF
    exportgraphics(fig_handle, strcat(file_root,'.pdf'), ...
        'ContentType','vector', 'BackgroundColor','none');
else
    exportgraphics(fig_handle, strcat(file_root,'.png'), ...
        'ContentType','image', 'Resolution',300);
end
end
