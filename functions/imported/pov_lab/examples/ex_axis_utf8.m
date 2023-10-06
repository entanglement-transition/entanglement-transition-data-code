%% 'axis' method usage example
% Common setup for all examples

ex_setup
% Start timer

tic
% Povlab object

pl = povlab( povray_version,...
             povray_path, ...
             povray_out_dir);
% Camera

pl.include_begin('camera');
    pl.camera('angle', 35, 'location', [12 14 8], 'look_at', [0 0 0.5]);
pl.include_end();
% Llights

pl.include_begin('lights');
    pl.light('location', [10 10 7],     'color', [3 3 3], 'shadowless', true);
    pl.light('location', [-10 10 30],   'color', [3 3 3], 'shadowless', true);
pl.include_end();
% Floor

pl.include_begin('floor');
    % Texture
    finish = "phong 1 reflection {0.02 metallic 0.8}";
    tex_grey   = pl.declare("tex_grey",   pl.texture('pigment_odd', [0.01 0.01 0.01], 'pigment_even', [0.5 0.5 0.5], 'finish', finish));
    
    % Plane
    pl.plane('normal', [0,0,1], 'limits', [-4 -4 4 4], 'texture', tex_grey);
pl.include_end();
% Default axis

pl.scene_begin('scene_file', 'axis.pov', 'image_file', 'axis.png');
    pl.include("camera");
    pl.include("lights");
    pl.include("floor");
    
    pl.axis('length', [5 5 4]);
pl.scene_end();
% Render and display

image = pl.render();
imshow(image);
% *Axis of custom color and thickness*  


pl.scene_begin('scene_file', 'axis_colored.pov', 'image_file', 'axis_colored.png');
    pl.include("camera");
    pl.include("lights");
    pl.include("floor");
    
    % Texture
    tex_gray   = pl.declare("tex_gray",   pl.texture('pigment', [0.2 0.2 0.2]));
    tex_yellow = pl.declare("tex_yellow", pl.texture('pigment', [0.6 0.6 0.0]));
    tex_pink   = pl.declare("tex_pink",   pl.texture('pigment', [0.8 0.2 0.3]*1.2));
    tex_green  = pl.declare("tex_green",  pl.texture('pigment', [0.0 0.3 0.0]));

    pl.axis('length', [5 5 4], 'radius', 0.07, 'tex_common', tex_gray, 'tex_x', tex_yellow, 'tex_y', tex_pink, 'tex_z', tex_green);
pl.scene_end();
% Render and display

image = pl.render();
imshow(image);
% Elapsed time

toc