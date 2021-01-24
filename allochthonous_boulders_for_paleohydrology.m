% ¦PALEOHYDROLOGICAL FLOOD FLOW RECONSTRUCTION FROM BOULDER SIZE IN A
% FLUVIAL CHANNEL USING INPUT DATA FROM SURVEYED BOULDERS AND ASSOCIATED
% CHANNEL REACHES¦
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Use       ¦  plot_Q_and_h.m  ¦            script for plotting!
%
% Use       ¦  compare_VQh.m  ¦             script for better...
%
% comparision of the three different approaches: Costa (1983), Clarke
% (1996), and Alexander & Cooker(2016) and visualization.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% This Matlab® script calculates...
%
%
%       (1) peak fluid-flow velocities
%
%
% of a flood flow for a input set of...
%
%           - boulder diameters of allochthonous, fluvial transported 
%               boulders and
%           - associated estimated boulder rock material densities
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
%           - input topographic cross-sections and
%           - associated channel bed slopes
%
% expressing an approximation of a fluvial channel geometry,
% to calculate...
%
%
% (2) peak discharges and
%
% (3) maximum flow height
%
%
% able to sustain transport for the input set of boulder
% diameters and densities.
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
% - excel® table 
%                   ¦  boulders_import.xlsx  ¦
%   
%   with input set of boulder diameters and associated estimated boulder
%       rock material densities
%
% - excel® table
%                   ¦  topoprofile_import.xlsx  ¦
%
%   with input topographic cross-sections (d, z) associated and channel bed
%   slope (S)
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Output:
%
% - excel® table
%                   ¦  paleohydrology_results_from_boulders.xlsx  ¦
%
%
% (1) peak cross-sectional fluid-flow velocities of Costa (1983), Clarke
%     (1996), and Alexander & Cooker(2016) approaches (V) for the input
%     set of boulder diameters/densities
%
% (2) peak discharges (Q_eg) from Costa (1983), Clarke (1996), and
%     Alexander & Cooker(2016) approaches for the input set of
%     boulder diameters/densities and corresponding input topographic
%     profile data (incl. channel bed slope) applied
%
% (3) maximum flow height (h_eg) from Costa (1983), Clarke (1996), and
%     Alexander & Cooker(2016) approaches for the input set of
%     boulder diameters/densities and corresponding input topographic
%     profile data (incl. channel bed slope) applied
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
% define basic parameters

rho_f   = 1500;       % flood water density [kg/m^3]
g       = 9.81;       % gravitational acceleration [m/s^2]

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%for manningseq.m function, 
n       = 0.04;       % Manning's roughness coefficient for mountain 
                      %         streams (Chow, 1959), [s/(m^(1/3))]

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% input data

% ¦boulder diameters in [m] (maximum diameters) and¦
% ¦     density of boulder rock material [kg/m^3]  ¦

% input excel spreadsheet

filenames.boulders = 'boulders_import_centralNepal.xlsx'; 
    % adjust to appropriate filename!
    
[boulders_num, boulders_txt] = xlsread(filenames.boulders);     %read file
boulders.names      = boulders_txt(4:end,1);
boulders.diameters  = boulders_num(1:end,1);
boulders.densities  = boulders_num(1:end,2);
boulders.topos      = boulders_txt(4:end,4);

clearvars boulders_num boulders_txt

bouldercount= length(boulders.topos);

for i= 1:bouldercount
    boulders.topos_seperate{i,1} = ...
        strsplit(boulders.topos{i,1},',');
end
clearvars i

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% ¦topographic profiles¦

% input excel spreadsheet

filenames.topos = 'topoprofile_import_centralNepal.xlsx';
    % adjust to appropriate filename!

[topos_num, topos_txt] = xlsread(filenames.topos);               %read file
topocount = length(topos_txt(5,1:end))./3;
for i = 1:topocount %#ok<BDSCI>
topos.names{i}          = topos_txt(3,(i-1)*3+1);
topos.coordinates{i}    = topos_txt(4,(i-1)*3+1);
topos.distances{i}      = topos_num(1:end,(i-1)*3+1);
topos.elevations{i}      = topos_num(1:end,(i-1)*3+2);
topos.slopes{i}         = topos_num(1,(i-1)*3+3);
end
clearvars i topos_num topos_txt

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% velocity calculations

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Empirical approach after Costa (1983)
% (see function 'empiricalCosta1983.m' for more details).


V.Costa  = empiricalCosta1983(boulders.diameters);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Force balance approach followed by Clarke (1996)
% (see function 'Clarke1996.m' for more details).

% calculation loop

for i = 1:bouldercount
    for j = 1:length(boulders.topos_seperate{i,1})
    for k = 1:length(topos.names)
      if  strcmp(string(boulders.topos_seperate{i,1}{1,j}), string(k))
            V.Clarke(i,j) = Clarke1996(boulders.diameters(i),...
            topos.slopes{1,k},boulders.densities(i),rho_f,g);
      end
    end
    end
end

% Arithmetic mean from slope values of different profiles used for 
% calculation

for i = 1:bouldercount
    for j = 1:length(boulders.topos_seperate{i,1})
        V.Clarke_arithmean(i,1) = mean(V.Clarke(i,1:j));
    end
end

clearvars i j k


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Force balance approach after Alexander & Cooker (2016)
% (see function 'AlexanderCooker2016' for more details).


V.Alexander = ...
    AlexanderCooker2016(boulders.diameters,boulders.densities,rho_f,g);


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% discharge calculations

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Preparation of topo data for calculation loops below 
% (cut NANs in distance and elevation arrays).

for i = 1:topocount %#ok<BDSCI>
    topos.distances{1,i}(isnan(topos.distances{1,i}))=[];
    topos.elevations{1,i}(isnan(topos.elevations{1,i}))=[];
end

clearvars i

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Empirical approach after Costa (1983)


% calculation loop

for i = 1:bouldercount
    for j = 1:length(boulders.topos_seperate{i,1})
    for k = 1:length(topos.names)
      if  strcmp(string(boulders.topos_seperate{i,1}{1,j}), string(k))
    
[Q.Costa(i,j),h.Costa(i,j)] = ...
manningseq(V.Costa(i,1),topos.slopes{1,k},n,topos.distances{1,k},...
topos.elevations{1,k},false);
        
      end
    end
    end
end

% Arithmetic mean from values of different profiles used for calculation

for i = 1:bouldercount
    for j = 1:length(boulders.topos_seperate{i,1})
        Q.Costa_arithmean(i,1) = mean(Q.Costa(i,1:j));
        h.Costa_arithmean(i,1) = mean(h.Costa(i,1:j));
    end
end

clearvars i j k


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Force balance approach followed by Clarke (1996)


% calculation loop

for i = 1:bouldercount
    for j = 1:length(boulders.topos_seperate{i,1})
    for k = 1:length(topos.names)
      if  strcmp(string(boulders.topos_seperate{i,1}{1,j}), string(k))
    
[Q.Clarke(i,j),h.Clarke(i,j)] = ...
manningseq(V.Clarke(i,j),topos.slopes{1,k},n,topos.distances{1,k},...
topos.elevations{1,k},false);
        
      end
    end
    end
end

% Arithmetic mean from values of different profiles used for calculation

for i = 1:bouldercount
    for j = 1:length(boulders.topos_seperate{i,1})
     Q.Clarke_arithmean(i,1) = mean(Q.Clarke(i,1:j));
     h.Clarke_arithmean(i,1) = mean(h.Clarke(i,1:j));
    end
end

clearvars i j k


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Force balance approach after Alexander & Cooker (2016)


% calculation loop

for i = 1:bouldercount
    for j = 1:length(boulders.topos_seperate{i,1})
    for k = 1:length(topos.names)
      if  strcmp(string(boulders.topos_seperate{i,1}{1,j}), string(k))
    
[Q.Alexander(i,j),h.Alexander(i,j)] = ...
manningseq(V.Alexander(i,1),topos.slopes{1,k},n,topos.distances{1,k},...
topos.elevations{1,k},false);
        
      end
    end
    end
end

% Arithmetic mean from values of different profiles used for calculation

for i = 1:bouldercount
    for j = 1:length(boulders.topos_seperate{i,1})
        Q.Alexander_arithmean(i,1) = mean(Q.Alexander(i,1:j));
        h.Alexander_arithmean(i,1) = mean(h.Alexander(i,1:j));
    end
end

clearvars i j k
 

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% export data

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Compiling output table

column_headings = {...
'sample name', ...
'boulder diameter [m]', ...
'density of boulder rock material [kg/m^3]', ...
'topoprofiles used for calculation', ...
'flow velocity  [m/s] (averaged), Costa (1983)', ...
'flow velocity  [m/s] (averaged), Clarke (1996)', ...
'flow velocity  [m/s], Alexander&Cooker (2016)', ...
'flow discharge [m^3/s], Costa (1983)', ...
'flow discharge [m^3/s], Clarke (1996)', ...
'flow discharge [m^3/s], Alexander&Cooker (2016)', ...
'flow height [m], Costa (1983)', ...
'flow height [m], Clarke (1996)', ...
'flow height [m], Alexander&Cooker (2016)'};

export_data = zeros(bouldercount+1, length(column_headings));
export_data = num2cell(export_data);

export_data(1,:)      = column_headings;
export_data(2:end,1)  = boulders.names;
export_data(2:end,2)  = num2cell(boulders.diameters);
export_data(2:end,3)  = num2cell(boulders.densities);
export_data(2:end,4)  = boulders.topos;
export_data(2:end,5)  = num2cell(V.Costa);
export_data(2:end,6)  = num2cell(V.Clarke_arithmean);
export_data(2:end,7)  = num2cell(V.Alexander);
export_data(2:end,8)  = num2cell(Q.Costa_arithmean);
export_data(2:end,9)  = num2cell(Q.Clarke_arithmean);
export_data(2:end,10) = num2cell(Q.Alexander_arithmean);
export_data(2:end,11) = num2cell(h.Costa_arithmean);
export_data(2:end,12) = num2cell(h.Clarke_arithmean);
export_data(2:end,13) = num2cell(h.Alexander_arithmean);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Export data to Excel® spreadsheet file

filenames.export = 'paleohydrology_results_from_boulders.xlsx';
xlswrite(filenames.export,export_data,1)


%% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

% References:

% Alexander, J. and Cooker, M.J., 2016. Moving boulders in flash floods and
%   estimating flow conditions using boulders in ancient deposits. 
%   Sedimentology, 63(6): 1582-1595.
% Chow, V.T., 1959. Open channel hydraulics. McGraw-Hill Book Company, 
%   Inc; New York. 
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

