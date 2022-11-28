%{
Author: Amreen Imrit
Project: MIE301 Redesign Project
Description: Animates the path and linkage motion of a drafting pantograph at various scales, designed to work with any drawing
given to the code. Also shows the force analysis on each link at the various scales and drawing. 
This is designed to work with anydrawing, but the inpur has to be array with the various points making up a drawing, and points should be stored
in the array sequentially
%}

%Setting up
close all; % closes all figures
clear all; % clears all variables from memory
clc;       % clears all calculations from the Matlab workspace


xmin = -10; % leftmost window edge
xmax = 70;  % rightmost window edge
ymin = -40; % bottom window edge
ymax = 40;  % top window edge
%}

% Defining the major points on the pantograph, these will be used moving forward
%{
   A = fixed point
   B = pivot on the link attached to the groud
   C = top central pivot
   D = bottom central pivot (i.e input)
   E = pivot on the output link
   F = output point
   
   Pantograph rule ish
   BD = CE = AB
   BC = DE = EF
   AC = CF

   link 2 is AC
   link 3 is CF
   link 4 is DE
   link 5 is BD
%}

R2 = 33;    % Length of link 2 The link in real life is 34 cm. 33 cm is used in the analysis so that 0.5 cm is left on each side to add the pivots
R3 = 33;    %Length of link 3
%scale_vector = [2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 9, 10];
scale_vector = [2, 2.25, 2.5, 2.75, 3, 3.5, 4, 5, 6, 7, 8, 10];

%Variables to change
scale = scale_vector(10)                         %CHANGE the scale to try different scale factor

drawing_type = {'Circle', 'Square', 'Heart'};
original_drawing = 1;                           %CHANGE the number (1 - circle, 2 - rectangle, 3- heart) to try different drawings

yA = 0;                                         %y coordinate of fixed point
xA = 0;                                         %x coordinate of fixed point

BD = get_shortest_link_length(R2, R3, scale);   %length of shortest link

%Creating the x and y vector of point D so that it traces a circle
x_origin = sqrt(2)*R2/scale;                    %center of the input drawing. Value can be changes if you want the center to be somewhere else.
y_origin = 0;                                   %y coordinate of the center of the input drawing

%x_origin = 20 %IGNORE
%y_origin = 5 % IGNORE

[xD, yD] = input_coordinates (original_drawing, x_origin, y_origin, 60); %obtains coordinate of input **to change size of input drawing, need to change it in the fucntion. size set right now so that it work for every scale factor

%animating
%
figure (1);
set(1,'WindowStyle','Docked')

for i = 1:length(xD)
    hold off
    [theta5, theta2, theta3] = get_angles_for_point_postition(R2, BD, 0, 0, xD(i), yD(i)); %get the thetas of the some links at the specific position of D
    [xF(i), yF(i)] = draw_pantograph(R3, BD, 0, 0, xD(i), yD(i), theta5, theta2, theta3, false); %draws the pantograph at the specific position of D
    
    input_path = plot (xD, yD, 'm', 'LineWidth', 1);    %Draws the input path
    output_path = plot (xF, yF, 'k', 'LineWidth', 1);   %Draws the output path
    
    %sends path backwards so that it looks like the pantograph is on top of the path
    uistack(input_path,'bottom')    
    uistack(output_path,'bottom')
    
    %Labelling the graph
    xlabel('x (cm)', 'fontsize', 15);   % axis label
    ylabel('y (cm)', 'fontsize', 15);   % axil label 
    title_graph = strcat("Draws ", drawing_type(original_drawing), " with Scale Factor ", num2str(scale));
    title(title_graph);
    axis equal; 
    grid on;
    axis( [xmin xmax ymin ymax] );
    pause(0.05);
end
hold on;
plot (33*sqrt(2), 0, 'ro')
legend (" Output Path", " Input Path", " Link 2", " Link 3", " Link 4", "Link 5"); %USE only when unistack is on(i.e. uncommented)
%legend (" Link 2", " Link 3", " Link 4", " Link 5", " Input Path", " Output Path");

input_dist = calculate_total_distance_travelled(xD, yD)
output_dist = calculate_total_distance_travelled(xF, yF)
distance_travelled_ratio = output_dist/input_dist
%}

%
clear scale BD x_origin y_origin

figure(2);
set(2,'WindowStyle','Docked');
%scale_vector = [2.5, 2.75, 3, 3.5, 4, 5, 6];
friction_mag = 0.0065*0.2*9.81
for i=1:length(scale_vector)
%for i = 1:1
    scale = scale_vector(i); %scale 6, things looks weird, Weird in general something wrong
    BD = get_shortest_link_length(R2, R3, scale);
    x_origin = sqrt(2)*R2/scale;     
    y_origin = 0; 
    [xD_FA, yD_FA] = input_coordinates (original_drawing, x_origin, y_origin, 360);

    [F_input, thetaF, x, y, friction_dir_test, F_input_test] = perform_force_analysis(xD_FA, yD_FA, R2, BD);
    thetaF = thetaF*180/pi(); %Angle between the horizontal and the line joining point F to the center of the drawing[point (33*sqrt(2), 0)]
    F_input(1)
    hold on;
    plot (thetaF, F_input, 'linewidth', 2)
    %scatter (thetaF, F_input, 5, 'filled')
    xlabel('Angle of Output (deg)', 'fontsize', 12);   % axis label
    ylabel('Force at input (N)', 'fontsize', 12);   % axil label 
    title_graph2 = strcat ("Force Required at Input to Overcome Friction at Output (", drawing_type(original_drawing), ")"); 
    title(title_graph2);
    grid on;
    axis( [0 380 0.02 0.14] )
end
legend ("Scale: 2", "Scale: 2.25", "Scale: 2.5", "Scale: 2.75", "Scale: 3", "Scale: 3.5", "Scale: 4", "Scale: 5", "Scale: 6", "Scale: 7", "Scale: 8", "Scale: 10");
%legend ("Scale: 2.5", "Scale: 2.75", "Scale: 3", "Scale: 3.5", "Scale: 4", "Scale: 5", "Scale: 6")
%}

%--------------------------FORCE ANALYSIS FUNCTIONS------------------------
function [mag_F_input, theta_drawing, X, Y, friction_dir_test, F_input_test]=perform_force_analysis(xD, yD, R2, BD)
    for i=1:length(xD)
        [theta5, theta2, theta3] = get_angles_for_point_postition(R2, BD, 0, 0, xD(i), yD(i)); %get the thetas of the some links at the specific position of D
        [X(i, :), Y(i, :)] = draw_pantograph(R2, BD, 0, 0, xD(i), yD(i), theta5, theta2, theta3, true); %gets the coordinate of all points on the pantograph at any point in time
    end

    for i = 2:length(xD)
    %for i = 2:2
        %i = 10+1;
        % for link 3
        friction_mag = 0.0065*0.2*9.81; %calculate magnitude of friction force
        friction_dir = (-1/sqrt((X(i, 6) - X(i-1, 6))^2 + (Y(i, 6) - Y(i-1, 6))^2))*[X(i, 6) - X(i-1, 6), Y(i, 6) - Y(i-1, 6)]; %gets slope of line tangent to what is being drawn = LOA of kinetic friction
        %friction_dir_test(:,i) = friction_dir;
        %thetaF(i-1) = atan2(friction_dir(2), friction_dir(1));
        theta_drawing(i-1) = mod(atan2(Y(i, 6) - 0, X(i, 6) - 33*sqrt(2)), 2*pi());
        
        slope_friction =  friction_dir(2)/friction_dir(1); %gets LOA of friciton force
        friction_dir_test(i) = slope_friction;
        slope_F43 = (Y(i, 4) - Y(i, 5))/(X(i, 4) - X(i, 5)); %gets slope of line DE = line of action (LOA) F43
        slope_F23 = find_slope_and_force_magnitude (slope_friction, X(i,6), Y(i,6), slope_F43, X(i, 4), Y(i, 4), X(i, 3), Y(i, 3), false); %obtain slope from point C to to K3
        [F23, F43] = find_slope_and_force_magnitude (slope_F43, friction_mag*friction_dir(1), friction_mag*friction_dir(2), slope_F23, 0, 0, 0, 0, true);
        %
        %for link 2
        slope_F52 = (Y(i, 4) - Y(i, 2))/(X(i, 4) - X(i, 2)); %slope of F52 is the same as slope of line BD
        slope_F12 = find_slope_and_force_magnitude(slope_F23, X(i, 3), Y(i, 3), slope_F52, X(i, 4), Y(i, 4), 0, 0, false);
        [F12, F52] = find_slope_and_force_magnitude (slope_F52, -F23(1), -F23(2), slope_F12, 0, 0, 0, 0, true); %x value for 52 is -23 because we are suppose to use F32 = -F23

        %for link 5
        F_input = F43 + F52; %F_unknown + F25 + F45 = 0 and we know F45 = -F43 & F25 = -F52
        mag_F_input(i-1) = sqrt(F_input(2)^2+F_input(1)^2);
        F_input_test(i) =  F_input(2)/F_input(1);
        %}
    end
    [theta_drawing, order] = sort(theta_drawing); %ordering theta_drawing so that it goes from 0 to 360 in ascending order
    mag_F_input = mag_F_input(order); %reordering F-input based on how theta_drawing was ordered
end

function [tail_vector, head_vector] = find_slope_and_force_magnitude (slope1, x1, y1, slope2, x2, y2, x3, y3, force_polygon)
    %{
    Purpose: finds the point of intersection of 2 lines, and return either the
    slope of a third line that intersect the other lines at the same point or
    when force_polygon = true, it returns the vector representing the 2 unknown
    forces

    Fomula:
    Uses idea of force polygon to get the vector of the unknown forces in a
    force polygon, and it uses system if linear equation technique
    (elimination) to solve for point of intersection

    Parameters:
    slope1 -> slope of first line (when force_polygon = true slope2 is slope of force attached to head of known force)
    x1 -> x coordinate of a point on line 1
    y1 -> matching y coordinate of the point on line 1
    slope2 -> slope of second line (when force_polygon = true slope2 is slope of force attached to tail of known force)
    x2 -> x coordinate of a point on line 2
    y2 -> matching y coordinate of the point on line 2
    x3 -> x coordinate point of the third line whose slope we are trying to
    find (only used when force_polyon = false)
    y3 -> matching y coordinate point of the third line whose slope we are trying to
    find (only used when force_polyon = false)
    force_polygon -> bolean, true indicates we must return vector of unknown
    forces, false indicate to return only the slope of third line

    Returns: 
    when force_polygon == true
    tail_vector -> vector of force that is attached to the tail of the
    known force (force whose mag and dir we did not know at the beginning)
    head_vector -> vector of force that is attached to the head of the know
    force (force who mag we knew but not the dir)

    When force_polygon == false
    tail_vector -> return slope of third line
    head_vector -> not used
    %}

    b1 = y1-slope1*x1;
    b2 = y2-slope2*x2;

    if (slope1 == inf || slope1 == -inf)     %if line 1 is vertical, this handles it
        commonpoint_x = x1;
        commonpoint_y = slope2*commonpoint_x + b2;
    elseif (slope2 == inf || slope2 == -inf) %if line 2 is vertical, this hanfels it
        commonpoint_x = x2;
        commonpoint_y = slope1*commonpoint_x + b1;
    else                                     %line are not vertical, so normal elimination will work
        commonpoint_x = (b2-b1)/(slope1-slope2);
        commonpoint_y = slope1*commonpoint_x + b1;
    end 
    
    if (force_polygon == true) 
        tail_vector = [x2 - commonpoint_x, y2 - commonpoint_y]; %vector of force attached to the tail of known force
        head_vector = [commonpoint_x - x1, commonpoint_y - y1]; %vector of force attached to the head of known force
    else
        tail_vector = (y3 - commonpoint_y)/(x3 - commonpoint_x); %slope of unknown line (LOA of some force for a 3 force polygon)
    end
end

%--------------------------FUNCTIONS--------------------------------------
function BD = get_shortest_link_length(R2, R3, scale) 
    %{
    Purpose: gets the link length of the links that becomes smaller

    Fomula:
    It calculats the hypothernuse formed when link R2 and R3 are in a
    90 degree configuration, divides the hypotenuse by the scale factor
    giving us where the central pivot should be at when the pantgraph
    is in the 90 degree configuration. It then multiplies this value by
    cos 45 to get the length of BD

    Parameters:
    R2 -> Length of link 2 (or AC)
    R3 -> Length of link 3 (or CF)
    Scale -> Scale factor of the pantograph

    Returns: 
    BD -> the length of the link that connects the bottom central pivot
    to the link which is attahced to the ground
    %}
    
    hyp = sqrt(R2^2+R3^2);  %calculates hypotenuse
    x = hyp/scale;          % calculates distance from A to D when in the 90 degree configuration
    BD = x*cos(pi()/4);     % since the traingle ABD is an isoceles triagle, lenght of AD can be calculated
end

function [xD, yD] = input_coordinates (original_drawing, x_origin, y_origin, increment)
    %{
    Purpose: returns array of the x and y coordinate of the input drawing

    Parameters:
    original_drawing -> value range from 1-3, each represent which drawing
    we are doing
    x_origin -> center of the drawing (exept for heart is's the pointy thing on a heart)
    y_origin -> center of drawing
    increment -> how many points we want the drawing to have, lower
    increments makes animating faster, higher increment results in more
    accurate drawing and force analysis

    Returns: 
    xD -> array of input x coordinate
    yD -> array of input y coordinate
    %}

    %Input drawing parameters
    if original_drawing == 1
        r_circle = 1.5; %radius of circle
        %r_circle = 8
        [xD, yD] = circle_input_coordinates(r_circle, y_origin, x_origin, increment);     %obtain x,y of input path
    elseif original_drawing == 2
        len = 4; %dimension of rectangle
        width = 1; %dimension of rectangle
        %len = 3; %dimension of rectangle
        %width = 19; %dimension of rectangle
        [xD, yD] = square_input_coordinates (x_origin, y_origin, len, width, increment);  %obtain x,y of input path
    elseif original_drawing == 3
        size = 0.95; %size of heart
        %size = 3.5;
        [xD, yD] = heart_input_coordinates (x_origin, y_origin, size, increment);         %obtain x,y of input path
    end
end

function [theta5, theta2, theta3] = get_angles_for_point_postition(R2, BD, xA, yA, xD, yD)
    %{
    Purpose: computes theta5, theta2, theta3 based on the link lenghts
    and the position of point D and point A

    Parameters:
    R2 -> Length of R2 (or AC)
    BD -> legth of the smallest link
    xA -> x poistion of fixed point
    yA -> y position of fixed point
    xD -> x position of bottom central pivot
    yD -> y position of bottom central pivot

    Returns: 
    theta5 -> CCW angle from horizontal to link 5 (BD)
    theta2 -> CCW angle from horizontal to link 2 (AB, AC)
    theta3 -> ccw angle from horizontal to link 3 (CF, EF)
    %}

    %Finding Theta 5
    AD = sqrt((xA-xD)^2+(yA-yD)^2);                         %finds the length of AD
    angle_AD_to_BD = acos(AD/(2*BD));                       %using cosine law, the angle between AD and BD can be found
    CCW_angle_horizontal_to_AD = atan2(yA-yD, xA-xD);       %using slope af line FD and the fact the tangent of the angle of a line = its slope, the CCW angle can be found
    theta5 = CCW_angle_horizontal_to_AD - angle_AD_to_BD;   %knowing that CW_angle_horizontal_to_FD includes is theta5 + angle_AD_to_BD, theta 5 can be found
    
    %Finding theta2
    xB = xD + BD*cos(theta5);                               %find x coordinate of  B
    yB = yD + BD*sin(theta5);                               %find y coordinate of B
    theta2 = atan2 (yB - yA, xB - xA);                      %from the slope of line AB, theta 2 can be found
    
    %Finding theta3
    angle_AC_to_CF = pi() - 2*angle_AD_to_BD;               %all angles in a traingle sum to 180 gives angle FA to BD, then using the idea that all angles on a line sum to 180 gives angle BA to AD, finally using the rule of angles in a parallelogram angle FO to BO is found to = angle BA to AD
    theta3 = angle_AC_to_CF + theta2;                       %using trigonometry, theta 3 can be found
    
end

function [xReturn, yReturn] = draw_pantograph (R3, BD, xA, yA, xD, yD, theta5, theta2, theta3, force_analysis)
    %{
    Purpose: Draw the pantographs links at any position of the central
    pivot

    Parameters:
    R3 -> Length of R3 (or FB)
    BD -> legth of the smallest link
    xA -> x poistion of fixed point
    yA -> y position of fixed point
    xD -> x position of bottom central pivot
    yD -> y position of bottom central pivot
    theta5 -> CCW angle from horizontal to link 5 (BD)
    theta2 -> CCW angle from horizontal to link 2 (AB, BC)
    theta3 -> ccw angle from horizontal to link 3 (CF, EF)
    force_analysis -> boolean, says whether this is force analysis or not

    Returns: 
    xF -> x position of the output
    yF -> y position of the output

    if force analysis is true then
    xF -> array of position for x's of all points (order: A, B, C, D, E, F)
    yF -> array of position for y's of all points
    %}

    %calculate position of A using theta5 and the length of BD
    xB = xD + BD*cos(theta5); 
    yB = yD + BD*sin(theta5);
    
    %calculate position of B based using theta2 and the length of link 3 (= link 2)
    xC = xA + R3*cos(theta2);
    yC = yA + R3*sin(theta2);
    
    %calculate position of C based using theta3 and the length of AD (=BC)
    xE = xC - BD*cos(theta3);
    yE = yC - BD*sin(theta3);
    
    %calculate position of C based using theta3 and the link 3
    xF = xC - R3*cos(theta3);
    yF = yC - R3*sin(theta3);
        
    if (force_analysis == true)
        xReturn = [xA, xB, xC, xD, xE, xF];
        yReturn = [yA, yB, yC, yD, yE, yF];
    else
        xReturn = xF;
        yReturn = yF;
    end
    
    %Drawing the links
    if (force_analysis == false)
        plot([xA, xC], [yA, yC], 'Color','r','LineWidth',4)                             %draws link 2
        hold on;
        plot([xA, xB], [yA, yB], 'Color','r','LineWidth',4, 'HandleVisibility','off')   %draws line AB
        plot([xC, xF], [yC, yF], 'Color','b','LineWidth',4)                             %draws link 3
        plot([xC, xE], [yC, yE], 'Color','b','LineWidth',4, 'HandleVisibility','off')   %draws line CE
        plot([xD, xE], [yD, yE], 'Color','c','LineWidth',4)                             %draws link 4
        plot([xD, xB], [yD, yB], 'Color','g','LineWidth',4)                             %draws link 5

        %Add points + point Labels
        %draw base pivot for link2
        ground_size = 2.5;                                   
        plot([xA, ground_size], [yA, -ground_size], 'b', 'HandleVisibility', 'off');       
        plot([xA, -ground_size], [yA, -ground_size], 'b', 'HandleVisibility', 'off');   

        %draw circle for each pivot
        draw_pivot (xA, yA, 'A');
        draw_pivot (xB, yB, 'B');
        draw_pivot (xC, yC, 'C');
        draw_pivot (xE, yE, 'E');
        draw_pivot (xD, yD, 'D');
        draw_pivot (xF, yF, 'F');
    end
end

function draw_pivot (x, y, label)
    plot (x, y, 'ko', 'MarkerFaceColor','w', 'HandleVisibility','off');
    text(x+1.5, y, label, 'color', 'k', 'fontsize', 7);
end

function [xD, yD] = circle_input_coordinates(r_circle, y_origin, x_origin, increment)
   %{
    Purpose: creates 2 vectors, one for the x coordinates of point D and
    the other for the y coordinates of point D so as to draw a circle

    Parameters:
    r_circle -> Radius of the circle
    y-origin -> the y coordinate of the center of the circle
    x_origin -> the x coordinate of the center of the circle

    Returns: 
    xD -> x position of the central pivot (point D)
    yD -> y position of the central pivot (point D)
    %}

    theta_circle = linspace(0, 2*pi(), increment);
    yD = y_origin + r_circle*sin(theta_circle);
    xD = x_origin + r_circle*cos(theta_circle);
end

function [xD, yD] = square_input_coordinates (x_center, y_center, L, W, increment)
   %{
    Purpose: creates 2 vectors, one for the x coordinates of point D and
    the other for the y coordinates of point D so as to draw a square

    Parameters:
    x_center -> the x coordinate of the center of the circle
    y-center -> the y coordinate of the center of the circle
    L -> lenght of the square
    W -> width of the square

    Returns: 
    xD -> x position of the central pivot (point D)
    yD -> y position of the central pivot (point D)
    %}

    leftmost = x_center - W/2;      %the leftmost point of the square
    rightmost = x_center + W/2;     % the rightmodt point of the square
    topmost = y_center + L/2;       %the topmost point of the square
    bottommost = y_center - L/2;    %the bottommost point of the square
    increment = round(increment/(2*L+2*W));

    % Draw the bottom horizontal line
    line1X = linspace(leftmost, rightmost, increment*W); 
    line1Y = bottommost + zeros(1, length(line1X));

    %draws the right vertical line
    line2Y = linspace(bottommost, topmost,increment*L);
    line2X = rightmost + zeros(1, length(line2Y));
    
    %draws the top horizontal line
    line3X = linspace(rightmost, leftmost, increment*W);
    line3Y = topmost + zeros(1, length(line3X));

    %draws the left vertical line
    line4Y = linspace(topmost, bottommost, increment*L);
    line4X = leftmost + zeros(1, length(line4Y));

    % concatenates all previous so as to create the xD and yD vector
    xD = [line1X, line2X, line3X, line4X];
    yD = [line1Y, line2Y, line3Y, line4Y];
end

function [xD, yD] = heart_input_coordinates (x_center, y_center, size, increment)
    t = linspace(-pi,pi, increment);
    xD = x_center + size*t .* sin( pi * .872*sin(t)./t);
    yD = y_center+1 -size*abs(t) .* cos(pi * sin(t)./t);
end

function distance = calculate_total_distance_travelled (x, y)
    %{
    Purpose: calculated the total distance travelled by point O and point D

    Parameters:
    x -> vector containing all x coordinates
    y -> vector containing all y coordinates 

    Returns: 
    distance -> the total distnace travelled by point decribed by vector
    x,y
    %}

    distance = sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2);
    distance = sum (distance, 'all');
end
