function V_avg = Clarke1996(D_n,S,rho_m,rho_f,g)


% This Matlab-function calculates flow velocities from boulder sizes
% following Clarke (1996), who uses a force balance approach for incipient
% motion after Costa (1983). Incipient motion of boulders under turbulent 
% and Newtonian fluid-flow is presupposed. 
%
% F_c = F_l + F_d
%
% Critical force F_c is equal to the sum of the lift force F_l and drag
% force F_d operating on the boulder.
%
%
% input arguments:
%
%          D_n      nominal diameter in [m], 
%                       equal to (three major axes)^(1/3)
%           S       bed slope angle [rad]
%         rho_m     boulder density [kg/m^3]
%         rho_f     fluid density [kg/m^3]
%           g       gravitational acceleration [m/s^2]
%
%  output arguments:
%
%         V_avg     average velocity of stream flow
%
%
% Author: Marius Huber (marius.huber@posteo.net)
% Date: 24. October, 2018
%--------------------------------------------------------------------------

%set coefficients

mu_s = 0.625;  %coefficent of sliding friction for a cubic boulder
mu_r = 0.225;  %coefficent of rolling friction for a spherical boulder

C_l_c = 0.178;  %lift coefficient for a cubic boulder
C_l_s = 0.20;   %lift coefficient for a spherical boulder
C_d_c = 1.18;   %drag coefficient for a cubic boulder
C_d_s = 0.20;   %drag coefficient for a spherical boulder

%--------------------------------------------------------------------------

%mass of a cubic boulder
M_c = (D_n.^3).*rho_m;

%mass of a spherical boulder
M_s = [(D_n.^3).*(pi/6)].*rho_m;

%resisting force for cubic boulder
F_r_c = M_c.*[(rho_m-rho_f)/rho_m].*g.*(((cos(S))*mu_s)-(sin(S)));

%resisting force for spherical boulder
F_r_s = M_s.*[(rho_m-rho_f)/rho_m].*g.*(((cos(S))*mu_r)-(sin(S)));

%critical force for cubic boulder
F_c_c = F_r_c;

%critical force for spherical boulder
F_c_s = F_r_s;

%drag force for cubic boulder
F_d_c = (F_c_c.*C_d_c)./(C_l_c+C_d_c);

%drag force for spherical boulder
F_d_s = (F_c_s.*C_d_s)./(C_l_s+C_d_s);

%cross-sectional area of the boulder for a cubic boulder
A_c = D_n.^2;

%cross-sectional area of the boulder for a spherical boulder
A_s = ((D_n./2).^2).*pi;

%critical velocity for a cubic boulder
V_c_c = [(((F_d_c./C_d_c)./rho_f).*2)./A_c].^(0.5);

%critical velocity for a spherical boulder
V_c_s = [(((F_d_s./C_d_s)./rho_f).*2)./A_s].^(0.5);

%average critical velocity for both forms
V_critical = (V_c_c+V_c_s)./2;

%average velocity of stream flow (multiplication with 1.2)
V_avg = V_critical.*1.2; 
%adjusted to approximate mean velocities (Baker 1973, p. 27)

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%plotting
%figure(1)
%plot(D,V_avg);

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
end

% Reference:

% Clarke, A.O., 1996. Estimating probable maximum floods in the Upper Santa
%   Ana Basin, Southern California, from stream boulder size. Environmental
%   & Engineering Geoscience, 2(2): 165-182.
% Costa, J.E., 1983. Paleohydraulic reconstruction of flash-flood peaks 
%   from boulder deposits in the Colorado Front Range. Geological Society 
%   of America Bulletin, 94(8): 986-1004.

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Copyright (C) 2021 Marius Huber, full license notice in the read-me file

