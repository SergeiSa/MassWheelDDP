close all; clear; clc; %clear classes; %clc;

%Create ground link
Ground = SRD_get_Link_Ground();

index = 1;
Link_1 = SRD_get_Link(...
    'Order', index, ...
    'Name', 'Link_1', ...
    'RelativeBase', zeros(3, 1), ...
    'RelativeFollower', [0;0;-0.38], ...
    'RelativeCoM', [0;0;-0.3134], ...
    'Mass', 0.67, ...
    'Inertia', diag([0, 0.003792, 0]), ...
    'ToDisplay', true, ...
    'Color', 'r', ...
    'StlPath', []);

Joint_1 = SRD_get_Joint_PivotY(...
    'Name', 'Joint_1', ...
    'ChildLink',  Link_1, ...
    'ParentLink', Ground, ...
    'ParentFollowerNumber', 1, ...
    'UsedGeneralizedCoordinates', 1, ...
    'UsedControlInputs', [], ...
    'DefaultJointOrientation', eye(3));


index = index + 1;
Link_2 = SRD_get_Link(...
    'Order', index, ...
    'Name', 'Link_2', ...
    'RelativeBase', zeros(3, 1), ...
    'RelativeFollower', [0;0;0.5], ...
    'RelativeCoM', [0;0;0], ...
    'Mass', 0, ...
    'Inertia', diag([0, 0.0017, 0]), ...%0.0017
    'ToDisplay', true, ...
    'Color', 'r', ...
    'StlPath', []);

Joint_2 = SRD_get_Joint_PivotY(...
    'Name', 'Joint_1', ...
    'ChildLink',  Link_2, ...
    'ParentLink', Link_1, ...
    'ParentFollowerNumber', 1, ...
    'UsedGeneralizedCoordinates', 2, ...
    'UsedControlInputs', 1, ...
    'DefaultJointOrientation', eye(3));

LinkArray = [Ground; Link_1; Link_2]; %Create array of links
SRD_save(LinkArray, 'LinkArray');

% InitialPosition = [0, 0, pi/10, -pi/10-2*pi/30]; %Define initial position of the robot
InitialPosition = [pi, -pi/2]; %Define initial position of the robot
SRD_save(InitialPosition, 'InitialPosition');

Chain = SRD_Chain(LinkArray);
SRD_save(Chain, 'Chain');
%Chain.Update(InitialPosition)


SRD_SetupVisuals('AxisLimits', [-1; 1; -1; 1; -1; 1], ...
    'ViewAngle', [-37.5, 30], ...
    'ToDrawMeshes', false, ...
    'Animation_ToUseGrid', true, ...
    'Animation_ToUseGridMinor', true, ...
    'DrawRobot_Default_RobotColor', [0.6 0.3 0], ...
    'DrawRobot_Default_EdgeAlpha', 0.3, ...
    'DrawRobot_Default_FaceAlpha', 1, ...
    'DrawRobot_Default_LineWidth', 0.5, ...
    'DrawRobot_STL_FaceColor', [0.8 0.8 1.0], ...
    'DrawRobot_STL_EdgeColor', 'none', ...
    'DrawRobot_STL_FaceLighting', 'gouraud', ...
    'DrawRobot_STL_AmbientStrength', 0.15, ...
    'DrawRobot_STL_camlight', 'headlight', ...
    'DrawRobot_STL_material', 'metal', ... %shiny dull metal
    'ToDrawFrames', false, ...
    'DrawRobot_Frame_Scale', 0.2, ...
    'DrawRobot_Frame_LineWidth', 1);

DrawRobot_function = SRD_DrawRobot_get_function('DrawRobot_Type', 'Default', ... %'Default' or 'STL' or 'Custom'
    'DrawRobot_Custom_handle', [], ...
    'Function_Type', 'DrawGivenPosition', ... %'DrawGivenPosition' or 'DrawInitialPosition'  or 'DrawCurrentPosition'
    'Chain', Chain);

DrawRobot_function(InitialPosition, [])
SRD__make_default_scene('Default')
            