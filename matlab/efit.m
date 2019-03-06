classdef efit < handle
%classdef efit < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class serves as an interface for the EFIT file format G-EQDSK in
% MATLAB. You can create an instance and fill the properties with your data
% and then write the file. A test script "efit_test.m" should be given.
% Functionality: only write of files supported by now.
%##########################################################################
% description of G-EQDSK file format:
% -------------------------------------------------------------------------
% G EQDSK FORMAT
% Briefly, a right-handed cylindrical coordinate system (R, φ, Ζ) is used.
% The G EQDSK provides information on the pressure, poloidal current
% function, q profile on a uniform flux grid from the magnetic axis to
% the plasma boundary and the poloidal flux function on the rectangular
% computation grid. Information on the plasma boundary and the surrounding
% limiter contour in also provided.
%##########################################################################
% properties:
%--------------------------------------------------------------------------
% *) path
% *) CASE, idum, nw, nh
% *) rdim, zdim, rleft, zmid, rmaxis, zmaxis, simag, sibry, rcentr, bcentr,
% current, xdum
% *) fpol, pres, ffprim, pprime, psirz, qpsi
% *) nbbbs, limitr, rbbbs, zbbbs, rlim, zlim
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = efit(path, nw, nh)
% *) function fillempty(obj)
% *) function write(obj)
%##########################################################################

%author: Philipp Ulbl
%created: 08.02.2019
%modified: 05.03.2019

    properties (Access = public)
        path = {}; %path of file without name

        CASE = {'EFIT', '', ' 00000', '  0 ms', '', ''}; %Identification character string

        idum = 3; %dummy ??? set to 3 (because seen in many files).
        nw;       %Number of horizontal R grid points
        nh;       %Number of vertical Z grid points

        rdim;     %Horizontal dimension in meter of computational box
        zdim;     %Vertical dimension in meter of computational box
        rleft;    %Minimum R in meter of rectangular computational box
        zmid;     %Z of center of computational box in meter
        rmaxis;   %R of magnetic axis in meter
        zmaxis;   %Z of magnetic axis in meter
        simag;    %poloidal flux at magnetic axis in Weber /rad
        sibry;    %poloidal flux at the plasma boundary in Weber /rad
        rcentr;   %R in meter of vacuum toroidal magnetic field BCENTR
        bcentr;   %Vacuum toroidal magnetic field in Tesla at RCENTR
        current;  %Plasma current in Ampere
        xdum = 0; %dummy ??? set to zero.

        fpol;     %Poloidal current function in m-T, F = RB T on flux grid
        pres;     %Plasma pressure in nt / m^2 on uniform flux grid
        ffprim;   %FF'(psi) in (mT) / (Weber /rad) on uniform flux grid
        pprime;   %P'(psi) in (nt /m^2 ) / (Weber /rad) on uniform flux grid
        psirz;    %Poloidal flux in Weber / rad on the rectangular grid points
        qpsi;     %q values on uniform flux grid from axis to boundary

        nbbbs;    %Number of boundary points
        limitr;   %Number of limiter points

        rbbbs;    %R of boundary points in meter
        zbbbs;    %Z of boundary points in meter
        rlim;     %R of surrounding limiter contour in meter
        zlim;     %Z of surrounding limiter contour in meter
    end
    properties (Access = public)
        fname;
    end

    methods (Access = public)

        function obj = efit(path, nw, nh)
            %##############################################################
            %function obj = efit(path, nw, nh)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of class efit. Path of efit file and the number
            % of rows and columns of the array-like properties is needed.
            %##############################################################
            % path  ... location of the efit file (read from or write to)
            % nw    ... number of rows (width, horizontal grid points)
            % nh    ... number of cols (height, vertical grid points)
            %##############################################################

            obj.CASE{2} = datestr(datetime('today'), 'mm/dd/yy');

            obj.path = path;

            obj.nw = nw;
            obj.nh = nh;
        end

        function fillempty(obj)
            %##############################################################
            %function fillempty(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % fills all empty properties with zeros of the appropriate size.
            %##############################################################

            %fill empty scalar variables
            if isempty(obj.rdim), obj.rdim = 0; end
            if isempty(obj.zdim), obj.zdim = 0; end
            if isempty(obj.rleft), obj.rleft = 0; end
            if isempty(obj.zmid), obj.zmid = 0; end
            if isempty(obj.rmaxis), obj.rmaxis = 0; end
            if isempty(obj.zmaxis), obj.zmaxis = 0; end
            if isempty(obj.simag), obj.simag = 0; end
            if isempty(obj.sibry), obj.sibry = 0; end
            if isempty(obj.rcentr), obj.rcentr = 0; end
            if isempty(obj.bcentr), obj.bcentr = 0; end
            if isempty(obj.current), obj.current = 0; end

            %fill empty scalar variables
            if isempty(obj.fpol), obj.fpol = zeros(1, obj.nw); end
            if isempty(obj.pres), obj.pres = zeros(1, obj.nw); end
            if isempty(obj.ffprim), obj.ffprim = zeros(1, obj.nw); end
            if isempty(obj.pprime), obj.pprime = zeros(1, obj.nw); end
            if isempty(obj.qpsi), obj.qpsi = zeros(1, obj.nw); end

            %fill empty matrix variables
            if isempty(obj.psirz), obj.psirz = zeros(obj.nh, obj.nw); end

            %fill bbbs and limitr variables
            if isempty(obj.nbbbs), obj.nbbbs = 1; end
            if isempty(obj.limitr), obj.limitr = 1; end
            if isempty(obj.rbbbs), obj.rbbbs = zeros(1, obj.nbbbs); end
            if isempty(obj.zbbbs), obj.zbbbs = zeros(1, obj.nbbbs); end
            if isempty(obj.rlim), obj.rlim = zeros(1, obj.limitr); end
            if isempty(obj.zlim), obj.zlim = zeros(1, obj.limitr); end
        end

        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % writes all properties into a file in the specified location
            % (given by "path"). Checks the dimensions of the properties
            % beforehand.
            %##############################################################

            obj.check();

            dat = erase(obj.CASE{2}, '/');
            %generate filename without spaces
            fn = ['g', dat, '.', erase(obj.CASE{3}, ' '), ...
                            '.', erase(obj.CASE{5}, ' '), ...
                            '.', erase(obj.CASE{6}, ' ')];
            obj.fname = [obj.path, fn];

            %check if file exists and delete
            if isfile(obj.fname)
                delete(obj.fname);
            end

            %write 1st block: CASE strings, nw, nh and other scalars
            %write 1st line
            obj.write2000(obj.CASE, [obj.idum, obj.nw, obj.nh]);
            %write 2nd line
            obj.write2020([obj.rdim, obj.zdim, obj.rcentr, obj.rleft, obj.zmid]);
            %write 3rd line
            obj.write2020([obj.rmaxis, obj.zmaxis, obj.simag, obj.sibry, obj.bcentr]);
            %write 4th line
            obj.write2020([obj.current, obj.simag, obj.xdum, obj.rmaxis, obj.xdum]);
            %write 5th line
            obj.write2020([obj.zmaxis, obj.xdum, obj.sibry, obj.xdum, obj.xdum]);

            %write 2nd block: fpol, pres, ffprim, pprime, psirz, qpsi
            obj.write2020(obj.fpol);
            obj.write2020(obj.pres);
            obj.write2020(obj.ffprim);
            obj.write2020(obj.pprime);
            obj.write2020(obj.psirz');
            obj.write2020(obj.qpsi);

            %write last block with bbbs and lim data
            obj.write2022([obj.nbbbs, obj.limitr]);
            obj.write2020([obj.rbbbs; obj.zbbbs;]);
            obj.write2020([obj.rlim; obj.zlim;]);

        end

    end

    methods (Access = private)

        function check(obj)

            %check if nw, nh are empty
            if any(isempty([obj.nw, obj.nh]))
                error('nw an nh must not be empty!')
            end

            %check dimensions of vector-properties
            tochk = {obj.fpol, obj.pres, obj.ffprim, obj.pprime, obj.qpsi};
            for k = 1:numel(tochk)
                if numel(tochk{k}) ~= obj.nw
                    error('dim of "fpol, pres, ffprim, pprime or qpsi" does not match nw')
                end
            end

            %check dimensions of matrix-properties
            if any(size(obj.psirz) ~= [obj.nw, obj.nh])
                error('dim of "psirz" does not match nw, nh')
            end

        end
        function write2000(obj, varstr, varint)
            %2000 format (6a8,3i4)
            fmt = '%s%4i%4i%4i';

            fid = fopen(obj.fname, 'a');

            %make full string out of cell array of strings
            str = '';
            for k = 1:numel(varstr)
                str = [str, sprintf('%-8s', varstr{k})];
            end

            fprintf(fid, fmt, str, varint);
            fprintf(fid, '\n');

            fclose(fid);

        end
        function write2020(obj, var)
            %2020 format (5e16.9)
            fmt = '%16.9e';

            fid = fopen(obj.fname, 'a');

            %vector valued
            if any(size(var)==1)
                for k = 1:5:numel(var) %5-spacing to match format
                    fprintf(fid, fmt, var(k:(k + min(4, numel(var)-k))));
                    fprintf(fid, '\n');
                end

            %matrix valued
            elseif all(size(var)>1)
                for k = 1:5:numel(var) %inner-loop: 5-spacing to match format
                    fprintf(fid, fmt, var(k:(k + min(4, numel(var)-k))));
                    fprintf(fid, '\n');
                end

            %scalars, cells, etc.
            else
                error('var not printable.')
            end

            fclose(fid);

        end
        function write2022(obj, var)
            %2022 format (2i5)
            fmt = '%5i';

            fid = fopen(obj.fname, 'a');

            fprintf(fid, fmt, var);
            fprintf(fid, '\n');

            fclose(fid);
        end
    end
end