function [U,D]=backslip_elastic3D(m,xloc)

%%INPUTS:
%fault model (standard Okada paramterization -- must be column vector):
% m(1) = length (km)
% m(2) = width (km)
% m(3) = depth to down-dip edge (km)
% m(4) = dip (degrees)
% m(5) = strike (degrees)
% m(6) = east position of down-dip edge (km)
% m(7) = north position (km)
% m(8) = strike-slip (m)
% m(9) = dip-slip (m)
% m(10) = tensile (m)
% xloc = 2xn matrix of n station coordinates (km) -- first row x, second row y
%%OUTPUTS:
%U is a 3xn matrix of VELOCITIES (first row east, second row north, third row up) 


%elastic part for back-slip
[Ue,D,S,flag]=disloc3d(m',xloc,1,.25);%
U=-Ue;  %minus sign because it is backslip
D=-D;
