%% | Truman Abbe | Utah State University | truman.abbe23@gmail.com | June 6th, 2026 |
%% Plotting Script for AttitudeSim001_Results.mat
%% Quaternion RK45 Validation Against CSSWE Euler Angle RK45

clear; close all; clc;

%% -------------------- LOAD RESULTS --------------------------------------
data_filename = fullfile(pwd, 'AttitudeSim001_Results.mat');
S = load(data_filename);

T = S.T;
T_days = S.T_days;
Y = S.Y;
pointing_error = S.pointing_error;
B_hyst_history = S.B_hyst_history;
H_body_history = S.H_body_history;

fprintf('Loaded results from: %s\n', data_filename);

%% -------------------- PLOT SETTINGS -------------------------------------
singlePos = [0.13 0.18 0.82 0.74];

set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultLegendFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 16);
set(groot, 'defaultTextFontSize', 16);
set(groot, 'defaultLegendFontSize', 16);

c_wx   = [0 0.25 0.85];    % blue
c_wy   = [0.85 0.1 0.1];   % red
c_wz   = [0.93 0.69 0.13]; % yellow
c_wmag = [0 0 0];          % black
c_beta = [0 0 0];          % black

x_ticks_days = 0:1:14;
idx = 1:size(Y,1);

%% -------------------- HYSTERESIS LOOP: Y ROD ----------------------------
fig_loop_y = figure('Name','Hysteresis Loop - Y Rod','NumberTitle','off','Color','w');
set(fig_loop_y,'Units','inches');
set(fig_loop_y,'Position',[1 1 7.2 3.8]);
set(fig_loop_y,'PaperPositionMode','auto');

axy = axes('Parent',fig_loop_y,'Position',singlePos);
plot(axy, H_body_history(2,:), B_hyst_history(2,:), 'LineWidth', 1.4);

xlabel(axy,'H [A/m]','FontSize',18);
ylabel(axy,'B_{hyst} [T]','FontSize',18);
title(axy,'Y Rod','FontSize',18);
grid(axy,'on');
box(axy,'on');

axy.LineWidth = 1.2;
axy.TickDir   = 'in';
axy.Layer     = 'top';
axy.FontSize  = 16;

exportgraphics(fig_loop_y,'CSSWE_Hysteresis_Loop_Y-Rod.png','Resolution',600);

%% -------------------- HYSTERESIS LOOP: Z ROD ----------------------------
fig_loop_z = figure('Name','Hysteresis Loop - Z Rod','NumberTitle','off','Color','w');
set(fig_loop_z,'Units','inches');
set(fig_loop_z,'Position',[1 1 7.2 3.8]);
set(fig_loop_z,'PaperPositionMode','auto');

axz = axes('Parent',fig_loop_z,'Position',singlePos);
plot(axz, H_body_history(3,:), B_hyst_history(3,:), 'LineWidth', 1.4);

xlabel(axz,'H [A/m]','FontSize',18);
ylabel(axz,'B_{hyst} [T]','FontSize',18);
title(axz,'Z Rod','FontSize',18);
grid(axz,'on');
box(axz,'on');

axz.LineWidth = 1.2;
axz.TickDir   = 'in';
axz.Layer     = 'top';
axz.FontSize  = 16;

exportgraphics(fig_loop_z,'CSSWE_Hysteresis_Loop_Z-Rod.png','Resolution',600);

%% -------------------- ANGULAR VELOCITY ----------------------------------
fig_vel = figure('Name','Angular Velocity Damping Profile',...
                 'NumberTitle','off','Color','w');
set(fig_vel,'Units','inches');
set(fig_vel,'Position',[1 1 7.2 3.8]);
set(fig_vel,'PaperPositionMode','auto');

axv = axes('Parent',fig_vel,'Position',singlePos);

Y_deg = Y(:,5:7) * (180/pi);
w_total_deg = sqrt(sum(Y_deg.^2,2));

plot(axv,T_days(idx),Y_deg(idx,1),'Color',c_wx,'LineWidth',1.4); hold(axv,'on');
plot(axv,T_days(idx),Y_deg(idx,2),'Color',c_wy,'LineWidth',1.4);
plot(axv,T_days(idx),Y_deg(idx,3),'Color',c_wz,'LineWidth',1.4);
plot(axv,T_days(idx),w_total_deg(idx),'Color',c_wmag,'LineWidth',2.4);
yline(axv,0,'-k','HandleVisibility','off','LineWidth',1);

xlabel(axv,'Time [days]','FontSize',18);
ylabel(axv,'Angular Velocity [deg/sec]','FontSize',18);
legend(axv,{'$\omega_x$','$\omega_y$','$\omega_z$','$|\omega|$'},...
       'Interpreter','latex',...
       'Location','northeast',...
       'Box','on',...
       'FontSize',16);

grid(axv,'on');
box(axv,'on');
xlim(axv,[0 14]);
xticks(axv,x_ticks_days);
ylim(axv,[-10 15]);

axv.LineWidth = 1.2;
axv.TickDir   = 'in';
axv.Layer     = 'top';
axv.FontSize  = 16;

exportgraphics(fig_vel,'CSSWE_Angular_Velocity.png','Resolution',600);

%% -------------------- POINTING ERROR ------------------------------------
fig_point = figure('Name','Pointing Error',...
                   'NumberTitle','off','Color','w');
set(fig_point,'Units','inches');
set(fig_point,'Position',[1 1 7.2 3.8]);
set(fig_point,'PaperPositionMode','auto');

axp = axes('Parent',fig_point,'Position',singlePos);
plot(axp,T_days(idx),pointing_error(idx),'Color',c_beta,'LineWidth',1.4);

xlabel(axp,'Time [days]','FontSize',18);
ylabel(axp,'$\beta$ [degrees]','Interpreter','latex','FontSize',18);

grid(axp,'on');
box(axp,'on');
xlim(axp,[0 14]);
xticks(axp,x_ticks_days);
ylim(axp,[0 180]);

axp.LineWidth = 1.2;
axp.TickDir   = 'in';
axp.Layer     = 'top';
axp.FontSize  = 16;

exportgraphics(fig_point,'CSSWE_Pointing_Error.png','Resolution',600);

fprintf('Finished generating figures.\n');