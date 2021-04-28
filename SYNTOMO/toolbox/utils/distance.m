function [d,az] = distance(lat1,lon1,lat2,lon2)
% DISTANCE: Compute distance and bearing between geographic coordinate
% pairs.
%
% INPUT
%    lat1: Array of starting latitudes (degrees)
%    lon1: Array of starting longitudes (degrees)
%    lat2: Array of ending latitudes (degrees)
%    lon3: Array of ending longitudes (degrees)
%
% OUTPUT
%       d: Arc length between geographic coordinate pairs (degrees)
%      az: Bearing of second geographic coordinate with respect to first 
%          measured clockwise-positive from North (degrees)
%
% B. VanderBeek (OCT-2020)
%
% Created to reproduce results of function of same name in Matlab's Mapping
% toolbox which is not always included in the Matlab installation.
%

% Convert to radians
lat1 = lat1*pi./180;
lon1 = lon1*pi./180;
lat2 = lat2*pi./180;
lon2 = lon2*pi./180;

% Arc distance
d = acos(cos(lat1).*cos(lat2).*cos(lon1-lon2) + sin(lat1).*sin(lat2));

% Bearing
az = atan2(sin(lon2-lon1).*cos(lat2),...
    cos(lat1).*sin(lat2) - sin(lat1).*cos(lat2).*cos(lon2-lon1));

% Return degrees
d  = d*180./pi;
az = az*180./pi;
