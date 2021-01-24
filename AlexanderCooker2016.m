function V = AlexanderCooker2016(D,rho_m,rho_f,g)


% This Matlab-function calculates flow velocities from boulder diameters 
% after the approach of Alexander & Cooker (2016), which is a force balance
% approach applied to the initiation of boulder movement in a turbulent and
% Newtonian fluid flow with an additional impulsive force generated 
% by unsteady flow (incipient motion).
%
% F = F_d + F_a - F_f
%
% The force F is the sum of the drag F_d and the impulsive force F_a and
% minus the frictional force F_f between the boulder and the bed.
%
% input arguments:
%
%           D                  diameter of a spherical boulder in [m]
%         rho_m                boulder density [kg/m^3]
%         rho_f                fluid density [kg/m^3]
%           g                  standard acceleration due to gravity [m/s^2]
%        
% output arguments:
%         
%         V                    time-averaged fluid velocity
%
%
% Author: Marius Huber (marius.huber@posteo.net)
% Date: 24. October, 2018
%--------------------------------------------------------------------------

%set parameters

lambda  = 0.4;           %coefficient of friction
k       = 0.5;           %coefficient of shape
a       = 0.5;           %impulsive acceleration [m/s^2], 
%                                   rate of change in time of the fluid 
%                                   velocity relative to the boulder

% ¦Lambda is a dimensionless parameter describing resistance to flow  
% ¦exerted by an obstacle on the bed (see French, 1971). k is a  
% ¦dimensionless constant that depends on the shape of the boulder; 
% ¦0.5 for a sphere, 1 for a cylinder with vertically oriented axis.
                    
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                     

Vol = [(D.^3).*(pi/6)];  %volume of the boulder 
                         %      (we assume a spherical boulder)
A = [(D./2).^2].*pi;     %cross-sectional area of the boulder in plane 
                         %      perpendicular to the flow direction (we
                         %      assume a spherical boulder)
L = Vol./A;              %length associated with the boulder and its 
                         %      orientation L = volume / cross-sectional 
                         %      area

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

V = [(L.*g*2).*(((rho_m/rho_f)-1)*lambda-(k*(a/g)))].^(0.5);

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%plotting

%figure(1)
%plot(D,V_avg);

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
end

% Reference:

% Alexander, J. and Cooker, M.J., 2016. Moving boulders in flash floods and
%   estimating flow conditions using boulders in ancient deposits. 
%   Sedimentology, 63(6): 1582-1595.

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Copyright (C) 2021 Marius Huber, full license notice in the read-me file

