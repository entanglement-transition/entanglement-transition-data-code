% Total cleanup
clear povlab;
close all;

%
% Set version and path to POV-Ray, that installed in your system, output directory for scene files and rendered images
%
povray_version = "3.7";

% Init outout directory
povray_out_dir = get_root_dir() + "/examples/out";

if isunix
%     povray_path = '"/usr/local/bin/povray +A -W1920 -H1080 -L'+ povray_out_dir + '"';
    povray_path = '"/usr/local/bin/povray +A -W640 -H360 -L'+ povray_out_dir + '"';
elseif ispc
    % Windows version of Povray is not supporting image rendering parametrs in command line, so please
    % set up Povray from menu 'Render/Edit Settings' before using Povlab
    povray_path = "C:/Program Files/POV-Ray/v3.7/bin/pvengine64.exe";
    % povray_path = "C:/Program Files/POV-Ray/v3.8-beta/bin/pvengine64.exe";
elseif ismac
    % Was not tested
    povray_path = '"/usr/local/bin/povray +A -W800 -H600 -L'+ povray_out_dir + '"';
else
    disp('Unknown platform');
end

disp('Setup complete');