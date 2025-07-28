function [pxl_x, pxl_y] = ang2px(angle,distance,cfg)
%% Calculate number of px given settings
% Input     
%           angle:     Visual angle in ° (e.g. 90 (1/2*pi))
%           distance:  Distance to the screen, unit (m, cm, mm) corresponds to output
%           cfg:       Screen parameters
% Output
%           pxl_x:      Number of pixels horizontal
%           pxl_y:      Number of pixels vertical
%%
switch nargin 
    case 3
    otherwise 
    error('3 Input arguments required')
end 

if ~isfield(cfg, 'xlen') || ~isfield(cfg, 'ylen') || ~isfield(cfg, 'hPxl') || ~isfield(cfg, 'vPxl')
    error('Missing screen parameter')
end 

if length(angle)==1
    ob_size = 2*distance*tan(deg2rad(angle)/2); % the size of the object
    
    pxl_x = round(ob_size/(cfg.xlen / cfg.hPxl));
    pxl_y = round(ob_size/(cfg.ylen / cfg.vPxl));
    
elseif length(angle)==2
    ob_sizex = 2*distance*tan(deg2rad(angle(1))/2); % the size of the object X
    ob_sizey = 2*distance*tan(deg2rad(angle(2))/2); % the size of the object X
    
    pxl_x = round(ob_sizex/(cfg.xlen / cfg.hPxl));
    pxl_y = round(ob_sizey/(cfg.ylen / cfg.vPxl));
    
end
