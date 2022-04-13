% Export a plot as eps and as pdf.
%
% Uses matlab/octave function to write eps, and convert_eps_2_pdf for
% pdf.
%
% input:
% ------
% path: string, path where to store the files.
% name: string, name, without extension, under which to store the files.
% varargin: ????
%
% output:
% -------
% none
%
% side effects:
% -------------
% Creates the two files.
function export_plot( path, name, varargin )

    plot_device = '-depsc2';
    ext = '.eps';
    plot_name = [path, name, ext];
    print(plot_device, plot_name);
    convert_eps_2_pdf(plot_name);

end
