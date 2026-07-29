function [dudn] = normal_derivative(obPoint, elements, normal)

% normal_derivative calculates the normal derivative of u at obPoint
% Inputs:
%   obPoint: Coordinates of the observation point [x_ob, y_ob]
%   elements: Coordinates of the element point [x_el, y_el]
%   normal: Normal vector components [n_x, n_y]
% Outputs:
%   dudn: Normal derivative of u at obPoint

% Calculate 2*pi*r^2 where r is the distance between obPoint and elements
r_2_pi = 2 * pi * sum((obPoint - elements).^2);

% Calculate partial derivatives of u with respect to x and y
% u = (1/2*pi) * ln(1/r)
% du/dx = -(1/2*pi) * (x/r^2)
% du/dy = -(1/2*pi) * (y/r^2)
dudx = (obPoint(1) - elements(1)) / r_2_pi;
dudy = (obPoint(2) - elements(2)) / r_2_pi;

% Calculate the normal derivative dudn = n_x * dudx + n_y * dudy
dudn = normal(1) * dudx + normal(2) * dudy;

