% ¦FLOW VELOCITIES, PEAK DISCHARGES, AND MAXIMUM FLOW HEIGHTS OF FLUID
% FLOW: A COMPARISION OF THREE DIFFERENT APROACHES USING BOULDER SIZE FOR 
% PALEOHYDROLOGIC RECONSTRUCTIONS IN A FLUVIAL CHANNEL¦
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% This Matlab® script calculates and plots...
%
%
%       (1) peak fluid-flow velocities
%
%
% of a flood flow for a predefined range of...
%
%           - boulder diameters of allochthonous, fluvial transported
%               boulders and
%           - a estimated boulder rock material density
%
% by applying approaches after Costa (1983), Clarke (1996), and
% Alexander & Cooker(2016). These approaches are following the incipient
% motion principle: maximum transported grain sizes represent maximum flow
% conditions in a fluvial channel.
%
% Furthermore this script utilizes the Rosenwinkel et al. (2017)'s code
% manningseq.m, which applies the empirical Gauckler–Manning–Strickler
% formula to an...
%
%           - input topographic cross-section
%               (default or input optional) and
%           - channel bed slope (default or input optional)
%
%                           representing a fluvial channel, to calculate...
%
%
% (2) peak discharges and
%
% (3) maximum flow height
%
%
% able to sustain transport for the predefined range of boulder diameters.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Velocities employed to the Gauckler–Manning–Strickler formula are
% defined as cross-sectional average flow velocities.
% Be aware that the fluid flow velocities derived by the three different 
% approaches mentioned above are not implicitly defined as cross-sectional
% average flow velocities! Costa (1983) and Clarke (1996), use an upscaling
% from bed or critical velocity by a factor 1.2 to estimate "average
% velocity". Alexander & Cooker (2016) define their velocity as fluid
% velocity incident on one boulder (conversion is omitted for the
% computation that is set in this code here).
%
% Turbulent, Newtonian fluid flow behaviour is assumed for the calculations
% below. Calculations are not applicable for non-Newtonian, plastic fluids
% that perform laminar flow due to the establishment of shear strength in
% the fluid material.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% The function 'manningseq.m' required for this script uses itself the
% function 'getdistance2.m', which can be requested by the author, and the
% function 'fminsearchbnd.m' by John D'Errico, which is available on the
% Mathworks® File Exchange
% (http://www.mathworks.com/matlabcentral/fileexchange/8277).
%
% Look into read_me.docx file for further details and reference.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Input:
%
% - optional: topographic cross-section (d_eg, z_eg) and
%   channel bed slope (S_eg)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Output:
%
% (_) visualization of input topographic cross-sectional profile and flow
%     height of a turbulent, Newtonian fluid flow just sustaining transport
%     for a 3 m diameter boulder, following Clarke (1996):
%
% (1) peak cross-sectional fluid-flow velocities of Costa (1983), Clarke
%     (1996), and Alexander & Cooker(2016) approaches (V_eg) for a given
%     boulder diameter range (D)
%
% (2) peak discharges (Q_eg) from Costa (1983), Clarke (1996), and
%     Alexander & Cooker(2016) approaches for a given boulder diameter
%     range (D) and input topographic profile (incl. channel bed slope)
%     applied
%
% (3) maximum flow height (h_eg) from Costa (1983), Clarke (1996), and
%     Alexander & Cooker(2016) approaches for a given boulder diameter
%     range (D) and input topographic profile data (incl. channel bed
%     slope) applied
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Author: Marius Huber, ETH Zurich 2018 (marius.huber@posteo.net)
%
% Have fun!
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Let's start!

clear all
clc
close all

%myFilePath = 'C:\Users\Marius\Documents\ETH\Masterarbeit';
%addpath(myFilePath)
%addpath(genpath(myFilePath))

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% set input parameters

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

D       = 0:0.1:40;   %range of boulder diameters [m]

rho_m   = 2700;       %rock material density [kg/m^3]

rho_f   = 1200;       %flood water density [kg/m^3]

g       = 9.81;       %gravitational acceleration [m/s^2]

n       = 0.04;       %Manning's roughness coefficient for mountain streams
                      % (Chow, 1959), [s/(m^(1/3))]


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Command window conversations, choice of example topographoic cross-secion
%   input

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

disp(...
    'Do you want to integrate your own example topographic profile data?');
str = input('y/n[n]\n','s');

while(1)
if str == 'y'
    fprintf('\n');
    fprintf('Mention: distance and height value arrays for topographic\n');
    fprintf('profile input need to have the same length!\n');
    input('(press enter)\n','s');
    
    fprintf('Type in distance values first, seperated by spaces:\n');
    d_eg = input('(omit space at the end)\n','s');
    d_eg = strsplit(d_eg, ' '); d_eg = str2double(d_eg);
    
    fprintf('Type in elevation values now, seperated by spaces:\n');
    z_eg = input('(omit space at the end)\n','s');
    z_eg = strsplit(z_eg, ' '); z_eg = str2double(z_eg);
    
    fprintf('Finally, type in  a channel slope value [rad],[m/m]:\n');
    S_eg = input('(omit space at the end)\n','s');
    S_eg = str2double(S_eg);
    
    break
elseif str == 'n'                 %default example topographic profile data
    d_eg = [  0,  5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65,  70];
    z_eg = [110, 90, 70, 50, 30, 20, 10,  5, 10, 20, 30, 50, 70, 90, 110];
    S_eg = 0.03;
    break
elseif isempty(str)
    str = 'n';
else
    disp('Please choose "y" (yes) or "n" (no)!');
    return
    
    
end
end


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Calculations

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


V_eg.Costa   = empiricalCosta1983(D);
V_eg.Clarke  = Clarke1996(D,S_eg,rho_m,rho_f,g);
V_eg.Clarke(1) = 0;
V_eg.Alexander = AlexanderCooker2016(D,rho_m,rho_f,g);
V_eg.Alexander(1) = 0;

Q_eg.Costa   = zeros(1,length(D));
h_eg.Costa   = zeros(1,length(D));
for i= 1:length(Q_eg.Costa)
    [Q_eg.Costa(i),h_eg.Costa(i)] = ...
        manningseq(V_eg.Costa(i),S_eg,n,d_eg(:),z_eg(:),false);
end

Q_eg.Clarke   = zeros(1,length(D));
h_eg.Clarke   = zeros(1,length(D));
for i = 1:length(Q_eg.Clarke)
    [Q_eg.Clarke(i),h_eg.Clarke(i)] = ...
        manningseq(V_eg.Clarke(i),S_eg,n,d_eg(:),z_eg(:),false);
end

Q_eg.Alexander   = zeros(1,length(D));
h_eg.Alexander  = zeros(1,length(D));
for i = 1:length(Q_eg.Alexander)
    [Q_eg.Alexander(i),h_eg.Alexander(i)] = ...
        manningseq(V_eg.Alexander(i),S_eg,n,d_eg(:),z_eg(:),false);
end

clearvars i

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Plotting 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Visualization: Topo-profile and flow height of a turbulent, Newtonian
% fluid flow just sustaining transport for a 3 m diameter boulder,
% following Alexander and Cooker (2016)

figure(1)

compare = manningseq(V_eg.Alexander(31),S_eg,n,d_eg,z_eg);
title({'topo + flow height of flow sustaining 3 m diameter boulder';...
    'following Alexander&Cooker (2016)'});

clearvars compare

%---

% Peak fluid flow velocities vs. boulder diameter for the three different
%   approaches:

figure(2)

plot(D,V_eg.Costa,'m','LineWidth',2,'DisplayName','empirical Costa 1983');
hold on
plot(D,V_eg.Clarke,'c','LineWidth',2,'DisplayName','Clarke 1996');
plot(D,V_eg.Alexander,'g','LineWidth',2,'DisplayName',...
    'Alexander&Cooker 2016');

title('Peak fluid flow velocities vs. boulder diameter');
xlabel('boulder diameter [m]'); ylabel('flow velocity [m/s]');
legend('show','Location','southeast');
xlim([0 40]); %ylim([0 25]);


%---

% Peak discharge vs. boulder diameter for the three different approaches:

figure(3)

semilogy(D,Q_eg.Costa,'m','LineWidth',2,'DisplayName', ...
    'empirical Costa 1983'); hold on
semilogy(D,Q_eg.Clarke,'c','LineWidth',2,'DisplayName','Clarke 1996');
semilogy(D,Q_eg.Alexander,'g','LineWidth',2,'DisplayName', ...
    'Alexander&Cooker 2016');

title('Peak discharge vs. boulder diameter');
xlabel('boulder diameter [m]'); ylabel('peak discharge [m^3/s]');
legend('show','Location','southeast');
xlim([0 40]); %ylim([10^-1 10^5]);


%---

% Maximum flow height vs. boulder diameter for the three different
%   approaches:

figure(4)

plot(D,h_eg.Costa,'m','LineWidth',2,'DisplayName','empirical Costa 1983');
hold on
plot(D,h_eg.Clarke,'c','LineWidth',2,'DisplayName','Clarke 1996');
plot(D,h_eg.Alexander,'g','LineWidth',2,'DisplayName',...
    'Alexander&Cooker 2016');

title('Maximum flow height vs. boulder diameter');
xlabel('boulder diameter [m]'); ylabel('maximum flow height [m]');
legend('show','Location','southeast');
xlim([0 40]); %ylim([0 25]);


%% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

% References:

% Alexander, J. and Cooker, M.J., 2016. Moving boulders in flash floods and
%   estimating flow conditions using boulders in ancient deposits. 
%   Sedimentology, 63(6): 1582-1595.
% Clarke, A.O., 1996. Estimating probable maximum floods in the Upper Santa
%   Ana Basin, Southern California, from stream boulder size. Environmental
%   & Engineering Geoscience, 2(2): 165-182.
% Costa, J.E., 1983. Paleohydraulic reconstruction of flash-flood peaks 
%   from boulder deposits in the Colorado Front Range. Geological Society 
%   of America Bulletin, 94(8): 986-1004.
% Rosenwinkel, S., Landgraf, A., Schwanghart, W., Volkmer, F., Dzhumabaeva,
%   A., Merchel, S., Rugel, G., Preusser, F. and Korup, O., 2017. Late 
%   Pleistocene outburst floods from Issyk Kul, Kyrgyzstan? Earth Surface 
%   Processes and Landforms.


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Copyright (C) 2021 Marius Huber, full license notice in the read-me file

