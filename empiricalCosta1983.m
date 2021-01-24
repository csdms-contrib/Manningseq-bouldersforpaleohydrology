function V_avg = empiricalCosta1983(D_i)


% This Matlab-function calculates flow velocities from boulder sizes with 
% empirical approach after Costa (1983).
%
%
%  input:
%
%           D_i      intermediate boulder diameter in [m]
%
%  output:
%
%         V_avg      average velocity of stream flow
%
%
%
% Author: Marius Huber (marius.huber@posteo.net)
% Date: 24. October, 2018
%--------------------------------------------------------------------------

%adjustment of D-value

D_i_mm = D_i.*1000;  %now its [mm], like in the regression plot,
                     %    Figure 3, Costa (1983)

%coefficients

a = 0.27; %for particles coarser than 500 mm
b = 0.40; %for particles coarser than 500 mm

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%the empirical power function

V_avg = (D_i_mm.^b).*a;

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%plotting

%figure(1)
%plot(D,V_avg);

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
end

% Reference:

% Costa, J.E., 1983. Paleohydraulic reconstruction of flash-flood peaks 
%   from boulder deposits in the Colorado Front Range. Geological Society 
%   of America Bulletin, 94(8): 986-1004.

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Copyright (C) 2021 Marius Huber, full license notice in the read-me file

