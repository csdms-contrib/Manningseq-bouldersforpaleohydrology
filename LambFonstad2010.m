function [Q,V,A] = LambFonstad2010(D,S,rho_m,rho_f,g,phi,w,k_s)


% This Matlab-function calculates discharge after Lamb and Fonstad 2010, 
% which goes back to Parker (1991) and Lamb et al. (2008) as used in 
% Lang et al. (2013).
%
% input arguments:
%
% D      % maximum diameter [m], 
%           (actually intermediate axis, mean grain size)
% S      % slope [m/m]
% rho_m  % density of rock material [kg7m^3]
% rho_f  % density of rock material [kg7m^3]
% g      % standard acceleration due to gravity [m/s^2]
% phi    % slope of valley flanks [rad]
% w      % valley bottom width [m]
% k_s    % bed roughness length scale ks = 0.1 to 1 m
%
% output arguments:
%
% Q      % peak discharge in [m^3/s]
% V      % depth averaged velocity
% A      % cross-sectional area of the flow
%
%
% after Data Repository of Lang et al. (2013):
%
% ''We followed the approach of Lamb and Fonstad (2010) where
%
% Q = 8.1*A*[tau_b/rho_f]^(1/2)*[h/k_s]^(1/6)            (1)
%
% h is the flow depth and A is the cross sectional area of the flow, which
% we modeled as a trapezoidal valley
%
% A = (h^2/tan(phi)) + w*h                                  (2)
%
% w is the flat bottom width. Bed shear stress tau_b is
%
% tau_b = rho_f*g*h_r*S                                  (3)
%
% h_r is the hydraulic radius, closely approximated by mean depth h_mean
%
% h_r = [A*sin(phi)]/[2*h+w]  ~= h_mean                   (4)
%
% Using tau_b we solve for the intermediate axis length of a median block 
% size D using the relation:
%
% tau_c = 0.15*S^0.25                                   (5)
%
% for the critical stress for incipient motion from Lamb et al. (2008) and
% citations therein
%
% D = tau_b /[tau_c*g*(rho_m-rho_f)]                      (6)   .''
%
% Author: Marius Huber (marius.huber@posteo.net)
% Date: 24. October, 2018
%--------------------------------------------------------------------------

% calculate tau_b with (5) and (6):

tau_c = 0.15*S^0.25;
tau_b = D.*[tau_c*g*(rho_m-rho_f)];

% calculate hydraulic radius with (3)

h_r = tau_b./(rho_f*g*S);

% solving quadratic formula after inserting (2) into (4)
% x_1 = [(-1)*b + (b^2-4*a*c)^(1/2)] / 2*a  (this is the relevant value)
% x_2 = [(-1)*b - (b^2-4*a*c)^(1/2)] / 2*a  (not relevant)

a = sin(phi)/tan(phi);
b = h_r.*(-2) + sin(phi)*w;           %sin(phi)*w-2*h_r;
c = h_r.*(-1)*w;

h = [(b.*(-1))+[(b.^2)-(c.*4*a)].^(1/2)]/2*a;

% calculating the cross sectional area of the flow

A = ((h.^2)./tan(phi))+h.*w;

% calculating discharge

Q = 8.1*A.*[tau_b./rho_f].^(1/2).*[h./k_s].^(1/6);

V = Q./A;

%plotting

%figure (1)
%plot(D,V);
%semilogx(Q,D)


end

% References:

% Lamb, M.P., Dietrich, W.E., Aciego, S.M., DePaolo, D.J. and Manga, M., 
%   2008. Formation of Box Canyon, Idaho, by megaflood: Implications for 
%   seepage erosion on Earth and Mars. Science, 320(5879): 1067-1070.
% Lamb, M.P. and Fonstad, M.A., 2010. Rapid formation of a modern bedrock 
%   canyon by a single flood event. Nature Geoscience, 3(7): 477-481.
% Lang, K.A., Huntington, K.W. and Montgomery, D.R., 2013. Erosion of the 
%   Tsangpo Gorge by megafloods, eastern Himalaya. Geology, 41(9): 
%   1003-1006.
% Parker, G., 1991. Selective sorting and abrasion of river gravel. II: 
%   Applications. Journal of Hydraulic Engineering, 117(2): 150-171.

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Copyright (C) 2021 Marius Huber, full license notice in the read-me file


