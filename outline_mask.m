%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mask of Cell Outline, June 2020, Séraphin Bassas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: The purpose of this script, given an image of a cell, is to 
% create a mask of the outline of the cell border, and explore the 
% relationship between FRET ratio intensity as a function of distance from 
% an arbitrary point on the cell outline.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREDIT: N/A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE:  This program takes as input two .tif images, the first being the 
% FRET sensor intensity image of a cell, and the second being the FRET-
% ratio image of the same cell. The sensor intensity image is segmented and
% its outline is used as a mask to overlay on the FRET-ratio image for 
% subsequent computation of FRET-ratio values. The left most non zero pixel
% that results from the substraction of the cell outline mask from the
% FRET-ratio image is used as the arbitrary point of origin. From this
% point, all non zero pixel intensity values are classified in a TOP and a
% BOTTOM array based on whether they lie above or below the point of origin
% latitudinally. Pixel intensity values corresponding to FRET-ratio 
% intensity values are then plotted as a function of distance from the 
% point of origin.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCLAIMER: This script assumes that the input images will be FRET images
% and thus the segmentation steps are based on the expected fluoresence
% intensity of this imaging technique. Any images with differing
% fluorescence instensity peaks should change the value inputed to the
% imbinarize function in this script. The script also assumes that no
% manipulation of the matrix obtained from the "find" MATLAB function be
% performed before classification in top and bottom matrices, as these
% matrices rely on the fact that find scans pixels longitudinally first and
% then latitudinally .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear
clear
clc

%% Basic Functions .tiff setup
%This function saves a figure or single axes to one or more vector and/or 
%bitmap file formats
export_fig = @export_fig.export_fig;

% opening the Sensor Intensity file and reading it to the workspace
if ~ispc; menu('Please select Sensor Intensity file','OK'); end
[f, PathName] = uigetfile({'*.tif','Tagged Image File Format (*.tif)'}, 'Select Sensor Intensity file', mfilename);
SelectedFile = fullfile(PathName,f);
SensorIm = imread(SelectedFile);

% opening the Fret Ratio file and reading it to the workspace
if ~ispc; menu('Please select Fret Ratio file','OK'); end
[f2, PathName2] = uigetfile({'*.tif','Tagged Image File Format (*.tif)'}, 'Select Fret Ratio file', mfilename);
SelectedFile2 = fullfile(PathName2,f2);
FretIm = imread(SelectedFile2);

% throwing a warning and end if no file is selected in the dialog above
if ~f, warning('No file selected.'); return; end
if ~f2, warning('No file selected.'); return; end
% filename = [pathname f];
% path to the image in a output.mat file 
% save('output.mat', 'filename', '-mat'); %saves the variable containing the 


warning('off', 'Images:initSize:adjustingMag');

%% filtering and thresholding
FilteredIm = imgaussfilt(SensorIm,3);% smoothing with guassian filter

% EITHER

% ThreshVal = graythresh(FilteredIm); %auto thresholding using Otsu
% BinaryIm = imbinarize(FilteredIm, ThreshVal); %passing Otsu value as limit

% OR
BinaryIm = imbinarize(FilteredIm, 0.1);% thresholding image 
% 0.3 for Ajinkya & 0.1 for LP

IIT=imfill(BinaryIm,'holes');

%% erosion and substraction

StructElmnt = strel('diamond', 2);
ErodedIm = imerode(IIT, StructElmnt);
ErodedIm = imerode(ErodedIm,StructElmnt);
M = imsubtract(IIT, ErodedIm);

%% applying the mask

FretIm(~M)=0; %changes all the pixel values in FretIm that do not belong to
              %the mask to 0 value
              
%% finding non zero coordinates of an image
% find gives all the non zero elements in the image matrix
% it stores the row in Y, the column in X, and the value in Val
[Y, X ,Val] = find(FretIm);

%% visualization

f1 = figure('Name', 'Fret Image Processing Steps');
subplot(2,2,1)
imshow(SensorIm)
title('Original Fret Image')
subplot(2,2,2)
imshow(IIT)
title('Binarized')
subplot(2,2,3)
imshow(ErodedIm)
title('Eroded')
subplot(2,2,4)
imshow(M)
title('Mask')

f3=figure('Name', 'Cell Border Intensity Profile');
subplot(1,2,1)
CM = colormap(jet(50)); %changes the colormap display for imshow
imshow(FretIm, CM) %display what's left to FretIm once the mask is applied
subplot(1,2,2)

%here I create a scatter plot with the x-axis being the column, or 
%x-coordinate of a non zero value in original FretIm, and on the y-axis 
%its corresponding value
scatter(X,Val,'x')
hold on
plot(X(1), Y(1))
hold off

%% distance

%saving the first point that was identified as reference for distance
InitX = X(1);
InitY = Y(1);
%initializing an array with double variables
T = {double(0),double(0)};
B = {double(0),double(0)};

%creating a variable that tells the program if a pixel is on the top or 
%bottom of the cell
Boolean = true;
% looping through all the points that belong to the border and pushing them
% to either array T (for Top) or B (for Bottom) which allows me to get two
% arrays with points in ascending distance from the point of origin, along
% with their intensities
for k = 2:length(X) %starts at 2 because point 1 is the origin
    
    if Y(k) < Y(k-1) %in the data, this condition represents whether the 
        %subsequent series of points is part of the bottom or top
        
        Boolean = ~Boolean; %every time there is a change, boolean changes
    end
    if Boolean %here I compute the distanc between the point and the origin
        T = vertcat(T, { double(sqrt( (X(k)-InitX)^(2) + (Y(k)-InitY)^(2) )), double(Val(k))});
    else
        B = vertcat(B, { double(sqrt( (X(k)-InitX)^(2) + (Y(k)-InitY)^(2) )), double(Val(k))});
    end
end

%converting T and B to regular Matlab arrays
T = cell2mat(T);
B = cell2mat(B);

%visualizing the data in a plot representing intensity as function of 
%distance from the origin
f4= figure('Name', 'Border Intensity as a function of distance from the left most point');
subplot(1,2,1)
scatter( T(:,1), T(:,2), 25, 'x')
title('Top half of cell')
subplot(1,2,2)
scatter( B(:,1), B(:,2), 25, 'x')
title('Bottom half of cell')



