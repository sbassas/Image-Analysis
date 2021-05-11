%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pie Wedging Algorithm, August 28 2020, Séraphin Bassas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: The purpose of this algorithm is to determine the relationship
% between FRET ratio intensity as a function of the azimuthal angle of the
% bin in which a pixel is located. The algorithm takes the center of a cell
% and uses it as the middle of a circle which we separate into bins. Each
% bin is a slice of the pie/circle and the sum of the intensity of the
% pixels within each bin is then computed, as the intensity of the pixel is
% semi-quantitatively representative of the FRET-ration intensity at that
% location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREDIT: adapted from a matlab script for "pie wedging" given to me by Dr.
% Amy Sutton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE: This program takes as input two .tif images, the first being the 
% FRET sensor intensity image of a cell, and the second being the FRET-
% ratio image of the same cell. The sensor intensity image is segmented and
% its outline is used as a mask to overlay on the FRET-ratio image for 
% subsequent computation of FRET-ratio values. The masked image is then 
% divided into an even number of "slices of a pie" that resemble pie wedges
% so that the intensity of the cell can be interpreted at different angles 
% from the center. Speaking of, the center of the "pie" is taken as the 
% centroid of the cell in the sensor intensity image. The results are then 
% displayed in the form of an array called "Results" which features a first
% column with the sum of the pixel intensity values found within each pie 
% wedge (separated by pie wedge), and a second column with the orientation,
% in degrees, of the pie wedge. The orientation is determined in the 
% following way: East = 0 ; North = 90 ; West = 180 ; South = 270
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCLAIMER: This script assumes that the number of pie wedges to be
% created is even; the script will not work as intended for an odd number
% of wedges.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up
clear
clc

%% loading and processing image
%{
I = imread('/Users/seraphinbassas/Desktop/sensor-intensity.tif');
FretIm = imread('/Users/seraphinbassas/Desktop/fret-ratio.tif');
%}
%% alternative way to load images
export_fig = @export_fig.export_fig;

% opening the Sensor Intensity file and reading it to the workspace
[f, PathName] = uigetfile({'*.tif','Tagged Image File Format (*.tif)'}, 'Select Sensor Intensity file', mfilename);
SelectedFile = fullfile(PathName,f);
I = imread(SelectedFile);

% opening the Fret Ratio file and reading it to the workspace
[f2, PathName2] = uigetfile({'*.tif','Tagged Image File Format (*.tif)'}, 'Select Fret Ratio file', mfilename);
SelectedFile2 = fullfile(PathName2,f2);
FretIm = imread(SelectedFile2);

% throwing a warning and end if no file is selected in the dialog above
if ~f, warning('No file selected.'); return; end
if ~f2, warning('No file selected.'); return; end
% filename = [pathname f];
% save('output.mat', 'filename', '-mat'); %saves the variable containing the 
%path to the image in a output.mat file 

warning('off', 'Images:initSize:adjustingMag');

%% processing
%Sensor image
BW = imfill(imbinarize(imgaussfilt(I,3),0.1), 'holes');
ObjectInfo = regionprops('table', BW, 'centroid', 'MaxFeretProperties');
[MaxFeret,Index1] = max(ObjectInfo{:,2});

% H = bwferet(BW, 'MaxFeretProperties');
% [MaxFeret,Index] = max(H{:,1});


%% visualizing Feret Objects
%{
Temp = cell2mat(ObjectInfo{Index1,4});
points = zeros(2,2);
points(1,1) = Temp(1,1); points(1,2) = Temp(1,2); points(2,1) = Temp(2,1); points(2,2) = Temp(2,2);
figure, imshow(BW), hold on
drawline('Position', points)
hold off
%}

%% erosion, substraction, and creation of the mask

StructElmnt = strel('diamond', 2);
ErodedIm = imerode(BW, StructElmnt);
ErodedIm = imerode(ErodedIm,StructElmnt);
M = imsubtract(BW, ErodedIm);
FretIm(~M) = 0;

%% setttings

slice = 1;
num_slices = 44;  %this has to be an even number
Center = ObjectInfo{Index1,1};
radius_cell = ObjectInfo{Index1,2};
Results = zeros(num_slices,2);
CM = colormap(jet(60));

%% explanation
% num_slices dictates how many slices the pie is divided in
% center_line is a horizontal line that extends right from the
% center of the cell
% sign_x determines in which direction of the x-step
%% pie wedging

if slice == 1
    %first find starting line (horizontally to right of center point, 0pi angle)
    
    pie_angle = 360 / num_slices;
    slice_angle_bi = NaN(1,num_slices);
    
    p(1,1) = Center(1) + radius_cell;
    p(1,2) = Center(2); %center_pt_xy(2) - radius_cell;
    central_line = vertcat(Center,p);
    % imshow(I);
    % hold on
    % fig2 = gcf;
    % roi1 = drawline('Position',central_line);
    
    %% fings going on
    
    
    %determine slope of starting line (should be 0 for the horizontal line)
    %first pie slice will start pie_angle/2 to left of central_line
    dx = central_line(2,1) - central_line(1,1);
    dy = central_line(2,2) - central_line(1,2);
    slope = dy/dx;
    
    slice_line{1} = vertcat(Center,p); %temporary
    slice_endpt{1}= p;
    %delete(roi1);
    %roi3=drawline('Position',slice_line{1});
    %delete(roi3);
    
    imshow(I);
    hold on
    l = 0; num_angles = NaN;
    for s = 1:num_slices %remember origin is top left corner, not near center
        slice_angle_bi(1,s) = (s-1)* pie_angle;
        if slice_angle_bi(1,s) > 180
            slice_angle_bi(1,s) = -(num_slices - s+1)*pie_angle;
            if l == 0
                num_angles = s-1;
                l = 1;
            end
        end
        %% determining the sign of x for the direction of x-step
        
        if s == 1
            sign_x = -1; %central_line always horizontal towards +x direction, first step to left
            angle = pie_angle/2;
            i = 1;
        elseif s > 1
            %sign of x-step depends on quadrant
            sign_x = -sign(slice_line{s-1}(1,2)-slice_line{s-1}(2,2)); %origin is not near center
            angle = pie_angle;
            i = s-1;
        end
        
        %% defining the slices of the pie taking into consideration sign switches
        
        
        %find a temporary point that defines a triangle starting from the
        %center point that makes a right angle with the previous line.
        if abs(slope) > 0.001
            slope_perp = -1/slope; %slope of perpendicular line to line 1
            intercept_perp = slice_line{i}(2,2) - slope_perp*slice_line{i}(2,1); %using slope_perp and pt p on line 1
            
            %{
    %draw a line (plot) for perpendicular line (debug help)
    new_x = slice_line{i}(2,1)+50*sign_x;
    line_x = linspace(slice_line{i}(2,1),new_x,5);
    line_y = linspace(slice_line{i}(2,2),(slope_perp*new_x+intercept_perp),5);
    line_z = m1*ones(size(line_x));
    plot3(line_x,line_y,line_z,'-','LineWidth',3,'Color',[0 0 0]);
            %}
            
        %ROI debug
        %draw a line (plot) for perpendicular line (debug help)
        xp = slice_line{i}(2,1) + 50 * sign_x;
        yp = slope_perp * xp + intercept_perp;
        x2 = slice_line{i}(2,1);
        y2 = slice_line{i}(2,2);
        
        %roi_center = drawpoint('Position',[slice_line{i}(1,1) slice_line{i}(1,2)]);
        roi_line = drawline('Position',[slice_line{i}(1,1) slice_line{i}(1,2); slice_line{i}(2,1) slice_line{i}(2,2)]);
        roi_p = drawline('Position',[x2 y2; xp yp]);
        %delete(roi_center);
        delete(roi_p); delete(roi_line);
            
            
            %find intersection of perpendicular line with line 2
            length1 = sqrt((slice_line{i}(2,1)-slice_line{i}(1,1))^2 + (slice_line{i}(2,2)-slice_line{i}(1,2))^2); %length original line segment
            length2 = length1 / cosd(angle);
            %this is the length of the segment perpendicular to line 1 and
            %intersecting line 2
            length_perp = length2 * sind(angle);
            %this is the point on line 2, which will become the next slice
            %line, that intersects with the circle outside the scope of imshow
            x3 = sign_x*length_perp/sqrt(1+slope_perp^2)+slice_line{i}(2,1);
            y3 = slope_perp*x3 + intercept_perp;
            
            %for debugging, but all these points cannot be seen in imshow
            %because they are outside the scale of the imshow window
            roi2 = drawpoint('Position',[x3 y3]);
            delete(roi2);
            
        elseif abs(slope)< 0.001 %switch x and y axes
            if s == 1
                sign_x = 1; %first step underneath horizontal 0pi line, origin top left corner
            elseif s > 1
                sign_x = sign(slice_line{i}(2,1));
            end
            slope = dx / dy;
            slope_perp = -1/slope; %slope of perpendicular line to line 1
            intercept_perp = slice_line{i}(2,1) - slope_perp*slice_line{i}(2,2); %using slope_perp and pt p on line 1
            
            %{
            %ROI debug
            %draw a line (plot) for perpendicular line (debug help)
            yp = slice_line{i}(2,2)+50*sign_x;
            x2 = slice_line{i}(2,1); %plotting must be without x/y axes flipped
            y2 = slice_line{i}(2,2);
            xp = slope_perp*yp+intercept_perp;
            %roi_center = drawpoint('Position',[slice_line{i}(1,1) slice_line{i}(1,2)]);
            %roi_x2y2 = drawpoint('Position',[y2 x2]);
            roi_line = drawline('Position',[slice_line{i}(1,1) slice_line{i}(1,2); x2 y2]);
            roi_p = drawline('Position',[x2 y2; xp yp]);
            %delete(roi_center); delete(roi_x2y2);
            delete(roi_line); delete(roi_p);
            %}
            
            %find intersection of perpendicular line with line 2
            length1 = sqrt((slice_line{i}(2,2)-slice_line{i}(1,2))^2 + (slice_line{i}(2,1)-slice_line{i}(1,1))^2); %length original line segment
            length2 = length1 / cosd(angle);
            length_perp = length2 * sind(angle);
            %the x-step must be in the direction of the circle, depends on quadrant
            y3 = sign_x*length_perp/sqrt(1+slope_perp^2)+slice_line{i}(2,2);
            x3 = slope_perp*y3 + intercept_perp;
            
            %         %for debugging
            %         roi2 = drawpoint('Position',[x3 y3]);
            %         delete(roi2);
        end
        %% creating pie mask and finding intesnity values within each slice
        
        %creating a clone so each iteration doesn't modify FretIm
        %irreversibly
        Clone = FretIm; 
        %extracting the dimensions of the image to create a mask of
        %identical size for poly2mask function below
        Dimensions = size(FretIm);
        %extracting the vertices of the slice at this iteration
        pt(1,1) = Center(1); pt(1,2) = Center(2);
        pt(2,1) = p(1,1); pt(2,2) = p(1,2);
        pt(3,1) = x3; pt(3,2) = y3;
        %creating a mask out of the above points 
        SliceMask = poly2mask(pt(:,1),pt(:,2),Dimensions(1),Dimensions(2));
        %applying mask to FretIm that is already masked with sensor outline
        Clone(~SliceMask) = 0;
        %finding all non zero values within slice
        [row,col,Storage] = find(Clone);
        %summing all intensity values
        if ~(isempty(Storage)),Results(s,1) = sum(Storage); end
        %vectorizing the central line and the new slice line to calculate
        %orientation
        u = [(central_line(2,1) -central_line(1,1)) (central_line(2,2)-central_line(1,2))];
        v = [(pt(3,1) -pt(1,1)) (pt(3,2)-pt(1,2))];
        %finding the angle between the vectors in degrees
        CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
        Degrees = real(acosd(CosTheta));
        if s > (num_slices/2), Degrees = 180 + (180-Degrees); end
        Results(s,2) = Degrees;
        
        %visualizing how each part of the cell outline falls into which
        %slice, but produces as many figures as there are slices
        %         figure
        %         imshow(Clone,CM)
        
        %% updating variables for next iteration
        
        %now, find line 2 using the center_pt_xy, and pt [x3 y3] same as above
        dx = x3 - central_line(1,1);
        dy = y3 - central_line(1,2);
        slope = dy/dx; %new slope value for next iteration
        intercept = y3 - slope*x3;
        
        %find delta_x and delta_y for a given radius on that line
        delta_x = sign(dx) * sqrt(radius_cell^2 / (slope^2+1) );
        %delta_x = sign_x*sqrt(radius_cell ^2/(slope^2+1));
        delta_y = slope*delta_x;
        p(1,1) = delta_x + Center(1);
        p(1,2) = delta_y + Center(2);
        
        slice_line{s} = vertcat(Center,p);
        slice_endpt{s}= p;
        
        %{
    %draw a line (plot) - watch out for origin location on tiff images
    line_x = linspace(slice_line{s}(1,1),slice_line{s}(2,1),5);
    line_y = linspace(slice_line{s}(1,2),slice_line{s}(2,2),5);
    %line_z = m1*ones(size(line_x));
    %plot3(line_x, line_y,line_z,'-','LineWidth',3,'Color',[0 0.4470 0.7410]);
    plot(line_x, line_y,'-','LineWidth',1,'Color',[0 0.4470 0.7410]);
    %hold on
    %view(2);
        %}
        
        
        %ROI all slices visualization
                if s > 1
                    vertices = [slice_line{s-1}(1,:); slice_line{s-1}(2,:); slice_line{s}(2,:)];
                    roi_slice = drawpolygon('Position',vertices,'LineWidth',0.5); %draws closed shape
                    if s == num_slices
                        vertices = [slice_line{s}(1,:); slice_line{s}(2,:); slice_line{1}(2,:)];
                        roi_slice = drawpolygon('Position',vertices,'LineWidth',0.5); %draws closed shape
                    end
                end
        
    end
    
    %% to visualize a single slice whose index you can retrieve from the results variable
    %{
    %ROI visualization of single slice
    f2 = figure;
    imshow(FretIm, CM)
    
    slice_index = 16; %change this parameter according to the slice
    %you would like to visualize
    vertices = [slice_line{slice_index-1}(1,:); slice_line{slice_index-1}(2,:); slice_line{slice_index}(2,:)];
    roi_slice = drawpolygon('Position',vertices,'LineWidth',0.5); %draws closed shape
    if slice_index == num_slices
        vertices = [slice_line{slice_index}(1,:); slice_line{slice_index}(2,:); slice_line{1}(2,:)];
        roi_slice = drawpolygon('Position',vertices,'LineWidth',0.5); %draws closed shape
    end
    %}
    %{
    save figure showing radial slices
    headerdisp = '';
    saveas(fig2,strcat(folder, '/t', num2str(t), '_Displacement_',num2str(num_slices),'pieslices.jpg'));
    saveas(fig1,strcat(folder, '/t', num2str(t), '_deltaD_pieslices.svg'));
    saveas(fig1,strcat(folder, '/t', num2str(t), '_deltaD_pieslices.fig'));
    
    export_fig(fig2,[folder, '/t' num2str(t) '_Displacement_',num2str(num_slices),'pieslices.tif'],'-transparent','-r150'); %'-r150' increases pixel res
    close(fig2);
    %}
    
    slice = 0;

end

