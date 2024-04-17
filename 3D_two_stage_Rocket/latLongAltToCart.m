function [Rx, Ry, Rz] = latLongAltToCart(lat, long, alt, Re)
    % Convert latitude, longitude, and altitude to Cartesian coordinates
    % Inputs:
    %   lat: Latitude in degrees
    %   long: Longitude in degrees
    %   alt: Altitude in kilometers
    %   Re: Radius of the sphere in kilometers
    % Outputs:
    %   Rx: x-coordinate
    %   Ry: y-coordinate
    %   Rz: z-coordinate
    
    
    % Calculate Cartesian coordinates
    Rx = (Re + alt) * cos(lat) * cos(long);
    Ry = (Re + alt) * cos(lat) * sin(long);
    Rz = (Re + alt) * sin(lat);
end
