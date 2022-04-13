close all

%% Option:
% 1 ... Passing particle
% 2 ... Trapped particle
option = 1;

%% Libneo
addpath('~/GITHUB/libneo/matlab/EFIT')

%% Full orbit

switch option
    case 1
        full_orbit_rphiz = load('full_orbit_plot_rphiz_passing.dat');
    case 2
        full_orbit_rphiz = load('full_orbit_plot_rphiz_trapped.dat');
end
full_orbit_xyz = full_orbit_rphiz;
full_orbit_xyz(:,1) = full_orbit_rphiz(:,1).*cos(full_orbit_rphiz(:,2));
full_orbit_xyz(:,2) = full_orbit_rphiz(:,1).*sin(full_orbit_rphiz(:,2));

n_points_full = numel(full_orbit_rphiz(:,1));

off_set = 30;
R_max = max(full_orbit_rphiz(:,1)) + off_set;
Z_max = max(full_orbit_rphiz(:,3)) + off_set;
Z_min = min(full_orbit_rphiz(:,3)) - 1.6.*off_set;

%% Poincaré cut
switch option
    case 1
        poincare_rphiz = load('poincare_plot_phi_0_rphiz_passing.dat');
    case 2
        poincare_rphiz = load('poincare_plot_phi_0_rphiz_trapped.dat');
end
poincare_xyz = poincare_rphiz;
poincare_xyz(:,1) = poincare_rphiz(:,1).*cos(poincare_rphiz(:,2));
poincare_xyz(:,2) = poincare_rphiz(:,1).*sin(poincare_rphiz(:,2));



%% Equilibrium Data for Torus

e = efit('g26884.4300', [], []);
e.read();

phi = linspace(0, 1.5*pi, 100);
r = e.rbbbs .* 100;
z = e.zbbbs .* 100;

[R, PHI] = meshgrid(r, phi);
XB = R .* cos(PHI);
YB = R .* sin(PHI);
ZB = repmat(z, numel(phi), 1);


%% Plot settings
col = 0.7 .* [0, 1, 1];

% Number of pushings through tetrahedra
n_limit = 8000;

% Number of pushing line behind particle
switch option
    case 1
        n_show_points = 20;
    case 2
        n_show_points = 300;
end

%Logical for Poincaré cuts appearing when crossing phi=0 plane
L1 = ismember(full_orbit_rphiz([1:n_limit],:),poincare_rphiz);
L1 = L1(:,1);

% Number of Poincaré plots with slow velocity
switch option
    case 1
        n_slow = 3;
        n_slow_2 = 11;
    case 2
        n_slow = 2;
        n_slow_2 = 11;
end


%% Figure with plots
figure

% line behind particle
p0 = plot3(nan,nan,nan,'r-');
p0.LineWidth = 2.0;

hold on
% Particle
switch option
    case 1
        particle_name = 'Passing particle';
    case 2
        particle_name = 'Trapped particle';
end
p1 = scatter3(nan,nan,nan,'r','DisplayName', particle_name);
p1.MarkerFaceColor = 'r';


% Poincaré cut
p2 = scatter3(nan,nan,nan,'b','DisplayName','Poincaré cut');
p2.MarkerFaceColor = 'b';
p2.SizeData = 6.0;

% Torus
sb = surf(XB, YB, ZB, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Plasma boundary');
patch(reducepatch(surf2patch(sb), 1e-6), 'EdgeColor', 'none', 'HandleVisibility', 'Off')
set(sb,'FaceLighting','gouraud')
light


% Limits of plot
xlim([-R_max,R_max]);
ylim([-R_max,R_max]);
zlim([Z_min,Z_max]);

% View of 3D plot
az = 24.2848;
el = 60.7185;
view(az,el);

axis off


%% Plot Plane

threshold = 0; 
zp = [Z_min, Z_max];
xp = [min(poincare_rphiz(:,1)) - off_set, max(poincare_rphiz(:,1)) + off_set];
% Use the axes x and Y limits to find the co-ordinates for the patch
x1 = [ xp(1) xp(2) xp(2) xp(1)];
z1 = [ zp(1) zp(1) zp(2) zp(2)];y1 = ones(1,numel(z1))* threshold; 
v = patch(x1,y1,z1, 'b');
set(v,'facealpha',0.2);
set(v,'edgealpha',0.2);
set(gcf,'renderer','opengl') ;


%% Text
switch option
    case 1
        t = text(poincare_xyz(1,1)+35,poincare_xyz(1,2)-28,poincare_xyz(1,3)-10,'$\varphi=0$','Interpreter','latex');
    case 2
        t = text(poincare_xyz(1,1)+1,poincare_xyz(1,2)-60,poincare_xyz(1,3)-65,'$\varphi=0$','Interpreter','latex');
end
t.FontSize = 15;
t.Rotation = -14;
t.FontWeight = 'bold';


%% Video options
switch option
    case 1
        video_file_name = 'orbit_passing.mp4';
    case 2
        video_file_name = 'orbit_trapped.mp4';
end
v = VideoWriter(video_file_name,'MPEG-4');
v.Quality = 100;
v.FrameRate = 80;
open(v);

%shading interp

%% Loop for orbit
lh = legend([sb,p1,p2]);
lh.Position = [0.2279 0.7798 0.2009 0.1452];
lh.FontSize = 12;

counter_poincare = 0;
for i = 1:n_limit
    
   if (counter_poincare >= n_slow) && ~L1(i)
        if mod(i,8) ~=0, continue; end
        n_show_points = 50;
   end
   
   if (counter_poincare >= n_slow_2) && ~L1(i)
        if mod(i,16) ~=0, continue; end
        n_show_points = 50;
   end
   
   p1.XData = full_orbit_xyz(i,1);
   p1.YData = full_orbit_xyz(i,2);
   p1.ZData = full_orbit_xyz(i,3);
   
   % line behind particle
   if i <= n_show_points
      p0.XData = full_orbit_xyz([1:i],1);
      p0.YData = full_orbit_xyz([1:i],2);
      p0.ZData = full_orbit_xyz([1:i],3); 
   else
      p0.XData = full_orbit_xyz([i-n_show_points:i],1);
      p0.YData = full_orbit_xyz([i-n_show_points:i],2);
      p0.ZData = full_orbit_xyz([i-n_show_points:i],3);
   end
   
   %Plot Poincaré plots when crossing phi=0 plane
   if L1(i)
       p2.XData = full_orbit_xyz(L1(1:i),1);
       p2.YData = full_orbit_xyz(L1(1:i),2);
       p2.ZData = full_orbit_xyz(L1(1:i),3);
       
       counter_poincare = counter_poincare + 1;
   end

   frame = getframe(gcf);
   writeVideo(v,frame);
% drawnow
end

close(v)

    