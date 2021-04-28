function [lat2,lon2] = reckon(lat1,lon1,delta,az)
% RECKON: Compute coordinates along great circle given a starting position,
% distance, and bearing.
%
% INPUT
%    lat1: Array of starting latitudes (degrees)
%    lon1: Array of starting longitudes (degrees)
%   delta: Array of great circle arc lengths (degrees)
%      az: Array of bearings measured clockwise-positive from north (degrees)
%
% OUTPUT
%    lat2: Latitude along great circle (degrees)
%    lon2: Longitude along great circle (degrees)
%
% B. VanderBeek (OCT-2020)
%
% Created to reproduce results of function of same name in Matlab's Mapping
% toolbox which is not always included in the Matlab installation.
%

% Convert to radians
lat1  = lat1*pi./180;
lon1  = lon1*pi./180;
delta = delta*pi./180;
az    = az*pi./180;

% Compute coordinates along great-circle
lat2 = asin(sin(lat1).*cos(delta) + cos(lat1).*sin(delta).*cos(az));
lon2 = lon1 + atan2(sin(az).*sin(delta).*cos(lat1),cos(delta) - (sin(lat1).*sin(lat2)));

% Convert back to degrees
lat2 = lat2*180./pi;
lon2 = lon2*180./pi;
