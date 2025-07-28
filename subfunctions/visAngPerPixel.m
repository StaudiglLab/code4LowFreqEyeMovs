function DPP = visAngPerPixel(screenXorYCm, distance, screenXorYPix)

%% compute visual angle per pixel 
% DPP = visAngPerPixel(screenXorYCm, distance, screenXorYPix)

visAng = 2*atand(0.5*screenXorYCm/distance); % visual angle of the screen length in degrees 
DPP = visAng/screenXorYPix;  % visual angle per pixel 


end