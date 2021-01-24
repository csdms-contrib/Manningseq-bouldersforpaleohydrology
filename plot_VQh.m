% ¦PALEOHYDROLOGICAL PEAK DISCHARGE CALCULATION FROM BOULDER SIZE¦
% ¦PLOTTING¦
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
%
% This Matlab® script plots...
%
% - peak fluid-flow velocities (V) and
%
% - peak discharges (Q) and 
%
% - maximum flow height (h)
%
% for a fluid flow transporting allochthonous, fluvial transported boulders
% (diameters and densities) in a channel 
%
% after using the... 
%
% allochthonous_boulders_for_paleohydrology.m
%
% Matlab® script, which computes these values with the help of
% Rosenwinkel et al. (2017)'s manningseq.m code (application of the
% empirical Gauckler–Manning–Strickler formula).
% Input peak cross-sectional fluid-flow velocities were calculated with
% approaches after Costa (1983), Clarke (1996), and Alexander & Cooker
% (2016) utilizing input topographic profiles. These approaches are
% following the incipient motion principle: maximum transported grain sizes
% represent maximum flow conditions in a fluvial channel. Look into
% read_me.docx file and allochthonous_boulders_for_paleohydrology.m script
% for further details and reference.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Author: Marius Huber, ETH Zurich 2018 (marius.huber@posteo.net)
%
% Have fun!

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Plotting peak fluid-flow velocities (V)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

figure(1)

%it's a line plot
plot(boulders.diameters,V.Costa,'r+',...
    'DisplayName','empirical Costa 1983'); hold on  
plot(boulders.diameters,V.Clarke_arithmean,'r.',...
    'DisplayName','Clarke 1996');
plot(boulders.diameters,V.Alexander,'r*',...
    'DisplayName','Alexander&Cooker 2016');

% ---

%xlim([0 20]);          %adjust x-axis limits if needed
%ylim([0 40]);          %adjust y-axis limits if needed


title('peak fluid-flow velocity estimates');
xlabel('boulder diameter [m]');
ylabel('peak fluid-flow velocity [m/s]');
legend('show','Location','southeast');

% ---

%add horizontal lines for better comparision

line([0 40],[5 5],'HandleVisibility','off');
line([0 40],[10 10],'HandleVisibility','off');

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Plotting peak discharges (Q)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

figure(2)

%it's a semilogarithmic plot
semilogy(boulders.diameters,Q.Costa_arithmean,'k+',...
    'DisplayName','empirical Costa 1983'); hold on
semilogy(boulders.diameters,Q.Clarke_arithmean,'k.'...
    ,'DisplayName','Clarke 1996'); 
semilogy(boulders.diameters,Q.Alexander_arithmean,'k*',...
    'DisplayName','Alexander&Cooker 2016');

% ---

xlim([0 40]);           %adjust x-axis limits if needed
ylim([10^3 10^6]);      %adjust y-axis limits if needed


title('Peak discharge estimates');
xlabel('boulder diameter [m]'); 
ylabel('peak discharge [m^3/s]');
legend('show','Location','southeast');

% ---

%add horizontal lines for better comparision

line([0 40],[10^4 10^4],'HandleVisibility','off');
line([0 40],[10^5 10^5],'HandleVisibility','off');

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Plotting maximum flow height (h)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

figure(3)

%it's a line plot
plot(boulders.diameters,h.Costa_arithmean,'b+',...
    'DisplayName','empirical Costa 1983'); hold on  
plot(boulders.diameters,h.Clarke_arithmean,'b.',...
    'DisplayName','Clarke 1996');
plot(boulders.diameters,h.Alexander_arithmean,'b*',...
    'DisplayName','Alexander&Cooker 2016');

% ---

%xlim([0 20]);          %adjust x-axis limits if needed
%ylim([0 40]);          %adjust y-axis limits if needed


title('Maximum flow height estimates');
xlabel('boulder diameter [m]');
ylabel('maximum flow height [m]');
legend('show','Location','southeast');

% ---

%add horizontal lines for better comparision

line([0 40],[10 10],'HandleVisibility','off');
line([0 40],[20 20],'HandleVisibility','off');

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