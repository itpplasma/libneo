classdef efit < handle
%classdef efit < handle
%##########################################################################
% description of class:
%--------------------------------------------------------------------------
% This class serves as an interface for the EFIT file format G-EQDSK in
% MATLAB. You can create an instance and fill the properties with your data
% and then write the file. A test script "efit_example_read.m" and
% "efit_example_write.m" should be given.
% Functionality: write and read files is supported by now.
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
% *) fname
% *) CASE, idum, nw, nh
% *) rdim, zdim, rleft, zmid, rmaxis, zmaxis, simag, sibry, rcentr, bcentr,
%    current, xdum
% *) fpol, pres, ffprim, pprime, psirz, qpsi
% *) nbbbs, limitr, rbbbs, zbbbs, rlim, zlim
%##########################################################################
% methods:
%--------------------------------------------------------------------------
% *) function obj = efit(fname, nw, nh)
% *) function buildcase(obj, name, date, shot, time)
% *) function fillempty(obj)
% *) function write(obj)
% *) function read(obj)
% *) function plot1d(obj, var)
% *) function plot2d(obj)
% *) function plotsummary(obj)
%##########################################################################

%author: Philipp Ulbl
%created: 08.02.2019

    properties (Access = public)
        fname = {}; %name of file with path
        
        CASE = {'EFIT    ', '', '', ' #000000', '  0000ms', ''}; %Identification character string

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
    properties (Access = private)

    end

    methods (Access = public)

        function obj = efit(fname, nw, nh)
            %##############################################################
            %function obj = efit(fname, nw, nh)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % constructor of class efit. Path of efit file and the number
            % of rows and columns of the array-like properties is needed.
            %
            % write efit: fname is only the path. nw, nh must be set.
            % read efit: fname is path and name. nw, nh can be [].
            %##############################################################
            % fname ... filename of the efit file. can be only a path 
            %           (with \ at the end) for writing or path+name for
            %           reading.
            % nw    ... number of rows (width, horizontal grid points)
            % nh    ... number of cols (height, vertical grid points)
            %##############################################################

            obj.CASE{3} = datestr(datetime('today'), 'mm/dd/yy');

            obj.fname = fname;
            
            obj.nw = nw;
            obj.nh = nh;
        end
        
        function buildcase(obj, name, date, shot, time)
            %##############################################################
            %function buildcase(obj, name, date, shot, time)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Fills the character string CASE in the head of the file.
            % Parameters to this function can have either a value or can be
            % empty [].
            %##############################################################
            % name  ... name of the file. string with max 8-chars
            % date  ... date of the file/shot. string with max 8-chars
            % shot  ... shot number. max 6-digits
            % time  ... time of shot in ms. max 4-digits
            %##############################################################
            
            %check if variables are in input or not and write if yes
            if ~isempty(name)
                if numel(name) <= 8
                    obj.CASE{1} = name;
                else
                    error(['name <', name, '> should not exceed 8 chars.']);
                end
            end
            
            if ~isempty(date)
                if numel(name) <= 8
                    %8char date string
                    obj.CASE{2} = date;
                else
                    error(['date <', date, '> should not exceed 8 chars.']);
                end
            end
            
            if ~isempty(shot)
                %empty space + # followed by 6digit shot number
                obj.CASE{4} = [' #', sprintf('%06d', shot)];
            end
            
            if ~isempty(time)
                %8char string made of 4digit time in ms (zero padded),
                %padded with empty spaces
                obj.CASE{5} = sprintf('%8s', [sprintf('%04d', time), 'ms']);
            end
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
            if isempty(obj.rdim),   obj.rdim = 0; end
            if isempty(obj.zdim),   obj.zdim = 0; end
            if isempty(obj.rleft),  obj.rleft = 0; end
            if isempty(obj.zmid),   obj.zmid = 0; end
            if isempty(obj.rmaxis), obj.rmaxis = 0; end
            if isempty(obj.zmaxis), obj.zmaxis = 0; end
            if isempty(obj.simag),  obj.simag = 0; end
            if isempty(obj.sibry),  obj.sibry = 0; end
            if isempty(obj.rcentr), obj.rcentr = 0; end
            if isempty(obj.bcentr), obj.bcentr = 0; end
            if isempty(obj.current),obj.current = 0; end

            %fill empty scalar variables
            if isempty(obj.fpol),   obj.fpol = zeros(1, obj.nw); end
            if isempty(obj.pres),   obj.pres = zeros(1, obj.nw); end
            if isempty(obj.ffprim), obj.ffprim = zeros(1, obj.nw); end
            if isempty(obj.pprime), obj.pprime = zeros(1, obj.nw); end
            if isempty(obj.qpsi),   obj.qpsi = zeros(1, obj.nw); end

            %fill empty matrix variables
            if isempty(obj.psirz),  obj.psirz = zeros(obj.nh, obj.nw); end

            %fill bbbs and limitr variables
            if isempty(obj.nbbbs),  obj.nbbbs = 1; end
            if isempty(obj.limitr), obj.limitr = 1; end
            if isempty(obj.rbbbs),  obj.rbbbs = zeros(1, obj.nbbbs); end
            if isempty(obj.zbbbs),  obj.zbbbs = zeros(1, obj.nbbbs); end
            if isempty(obj.rlim),   obj.rlim = zeros(1, obj.limitr); end
            if isempty(obj.zlim),   obj.zlim = zeros(1, obj.limitr); end
        end

        function write(obj)
            %##############################################################
            %function write(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Writes all properties into a file in the specified location
            % (given by "fname"). Checks the dimensions of the properties
            % beforehand.
            %##############################################################

            obj.check();
            
            %check if file, then extract path
            if isfile(obj.fname)
                [obj.fname, ~, ~] = fileparts(obj.fname);
                obj.fname = [obj.fname, '/'];
            end
            %check if folder exists, otherwise create
            if ~isfolder(obj.fname)
                mkdir(obj.fname);
            end
            
            %generate filename without spaces
            name = ['g', erase(erase(obj.CASE{4}, ' '), '#'), '.',...
                             erase(erase(obj.CASE{5}, ' '), 'ms'),'.',...
                             erase(obj.CASE{1}, ' ')];
            
            %put together path+name
            obj.fname = [obj.fname, name];
            
            %check if file exists and delete if yes
            if isfile(obj.fname)
                delete(obj.fname);
                warning(['existing file <', obj.fname ,'> overwritten.']);
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

        function read(obj)
            %##############################################################
            %function read(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Reads all properties from a file in the specified location
            % (given by "fname").
            %##############################################################

            %check if file exists
            if ~isfile(obj.fname)
                error(['file <', obj.fname, '> does not exist.'])
            end
            
            %line counter - is increased each time a line is read
            ll = 0;
            
            %read 1st block: CASE strings, nw, nh and 4 lines of scalars
            [ll, obj.CASE, var] = read2000(obj, ll);
            [obj.idum, obj.nw, obj.nh] = deal(var{:});
            
            [ll, var] = read2020(obj, ll, 'line');
            [obj.rdim, obj.zdim, obj.rcentr, obj.rleft, obj.zmid]      = deal(var{:});
            [ll, var] = read2020(obj, ll, 'line');
            [obj.rmaxis, obj.zmaxis, obj.simag, obj.sibry, obj.bcentr] = deal(var{:});
            [ll, var] = read2020(obj, ll, 'line');
            [obj.current, obj.simag, obj.xdum, obj.rmaxis, obj.xdum]   = deal(var{:});
            [ll, var] = read2020(obj, ll, 'line');
            [obj.zmaxis, obj.xdum, obj.sibry, obj.xdum, obj.xdum]      = deal(var{:});
            
            %read 2nd block: fpol, pres, ffprim, pprime, psirz, qpsi
            [ll, obj.fpol]   = read2020(obj, ll, 'vector');
            [ll, obj.pres]   = read2020(obj, ll, 'vector');
            [ll, obj.ffprim] = read2020(obj, ll, 'vector');
            [ll, obj.pprime] = read2020(obj, ll, 'vector');
            [ll, obj.psirz]  = read2020(obj, ll, 'matrix', [obj.nw, obj.nh]);
            [ll, obj.qpsi]   = read2020(obj, ll, 'vector');
            
            %write last block with bbbs and lim data
            [ll, obj.nbbbs, obj.limitr] = read2022(obj, ll);
            [ll, var] = read2020(obj, ll, 'matrix', [obj.nbbbs, 2]);
            obj.rbbbs = var(:, 1);
            obj.zbbbs = var(:, 1);
            [~, var]  = read2020(obj, ll, 'matrix', [obj.limitr, 2]);
            obj.rlim  = var(:, 1);
            obj.zlim  = var(:, 1);
        end
        
        function plot1d(obj, var)
            %##############################################################
            %function plot1d(obj, var)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Plots any 1d variable specified by var.
            %##############################################################
            % var   ... name of the 1d variable as a string
            %           possibilities: f, p, ff, pp, psi, q
            %##############################################################

            switch var
                case 'f'
                    y = obj.fpol; l = 'm T';
                case 'p'
                    y = obj.pres; l = 'N m^{-2}';
                case 'ff'
                    y = obj.ffprim; l = 'm T';
                case 'pp'
                    y = obj.pprime; l = 'N m^{-2} Wb^{-1}';
                case 'psi'
                    %index of the middle row/col of psi
                    ind = ceil(size(psi, 1)/2);
                    y = obj.psirz(ind, ind:end);
                    l = 'Wb (= T m^2)';
                case 'q'
                    y = obj.qpsi;
                    l = '1';
            end
            
            x = linspace(obj.rcentr, obj.rleft + obj.rdim, obj.nw);

            plot(x, y, '-r', 'DisplayName', var);
            xlabel('r / m')
            ylabel([var, ' / ', l])
            legend
        end
        
        function plot2d(obj)
            %##############################################################
            %function plot2d(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Plots the flux as 2d contour lines.
            %##############################################################
            
            x = linspace(obj.rleft, obj.rleft + obj.rdim, obj.nw);
            y = linspace(obj.zmid-obj.zdim, obj.zmid+obj.zdim, obj.nh);

            [X, Y] = meshgrid(x, y);
            
            contour(X, Y, obj.psirz')
            xlabel('r / m')
            ylabel('z / m')
            c = colorbar;
            ylabel(c, '\psi / Wb (= T m^2)')
            axis image
        end
        
        function plotsummary(obj)
            %##############################################################
            %function plotsummary(obj)
            %##############################################################
            % description:
            %--------------------------------------------------------------
            % Creates an OMFIT-style summary plot of the EFIT-file.
            %##############################################################
            
            %plot psi as 2d on the left
            subplot(1, 2, 1)
            obj.plot2d();
            
            %plot 4 profiles on the right
            
            %plot p and q in 1st row
            subplot(2, 4, 3)
            obj.plot1d('p');
            subplot(2, 4, 4)
            obj.plot1d('q');
            
            %plot pp, f in 2nd row
            subplot(2, 4, 7)
            obj.plot1d('pp');
            subplot(2, 4, 8)
            obj.plot1d('ff');
        end
    end

    methods (Access = public)

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

            %open file with permission to write. append data
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

            %open file with permission to write. append data
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
                error('var not writeable.')
            end

            fclose(fid);

        end
        function write2022(obj, var)
            %2022 format (2i5)
            fmt = '%5i';

            %open file with permission to write. append data
            fid = fopen(obj.fname, 'a');

            fprintf(fid, fmt, var);
            fprintf(fid, '\n');

            fclose(fid);
        end
        
        function [ll, varstr, varint] = read2000(obj, ll)
            %2000 format (6a8,3i4)
            
            %open file with permission to read
            fid = fopen(obj.fname, 'r');
            
            %read single line
            str = fgetl(fid);
            
            varstr = cell(1, 6);
            %extract 6*8characters for the field CASE
            for k=1:numel(obj.CASE)
                varstr{k} = str(((k-1)*8+1):(k*8));
            end
            
            %extract 3 integers
            str = str(49:end);
            varint = strsplit(strtrim(str), ' ');
            varint = cellfun(@str2num, varint, 'UniformOutput', false);
            
            ll = ll + 1;
            
            fclose(fid);
        end
        function [ll, var]            = read2020(obj, ll, type, siz)
            %2020 format (5e16.9)
            
            %open file with permission to read
            fid = fopen(obj.fname, 'r');
            
            %skip already read lines
            for k=1:ll
                fgets(fid);
            end
            
            %read single line with 5 scalars
            if strcmp(type, 'line')
                %get single line, read in exponential format and convert to
                %cell array
                var = num2cell(sscanf(fgetl(fid), '%e'));
                ll = ll + 1;
 
            %read a vector
            elseif strcmp(type, 'vector')
                var = [];
                %calculate the number of rows (5 cols)
                nlines = ceil(obj.nw/5);
                %append new elements to existing vector
                for k=1:nlines
                    var = [var; sscanf(fgetl(fid), '%e')];
                end
                ll = ll + nlines;
                
            %read a matrix
            elseif strcmp(type, 'matrix')
                var = zeros(siz);
                %calculate the number of rows (5 cols)
                nlines = ceil(numel(var)/5);
                %fill successivly 5 elements until the number of elements
                %of the matrix is reached
                for k=1:nlines
                    var(((k-1)*5+1):min(k*5, numel(var))) = ...
                        sscanf(fgetl(fid), '%e');
                end
                ll = ll + nlines;
                
            else
                error(['type ', type, ' unknown.']);
            end
            
            fclose(fid);
        end
        function [ll, var1, var2]     = read2022(obj, ll)
            %2022 format (2i5)
            
            %open file with permission to write. append data
            fid = fopen(obj.fname, 'r');
            %skip already read lines
            for k=1:ll
                fgets(fid);
            end
            
            %extract 2 integers
            str  = strsplit(strtrim(fgetl(fid)), ' ');
            var1 = str2double(str{1}); %str2double faster than str2num
            var2 = str2double(str{2});
            
            ll = ll + 1;
            
            fclose(fid);
        end
    end
end