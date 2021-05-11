%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concentric Circles, June 2020, Séraphin Bassas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: The purpose of this script is to represent the relationship
% between the FRET intensity ratio of a protein expressed as a function of
% distance from the leading edge of the cell containing said protein. This
% infromation is extracted from a time-lapse movie of a cell on which frame
% by frame TFM analysis has been carried out.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREDIT: N/A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE: This script takes as input two stacks of same thickness. It asks
% for the stack of Fret images (.lif files) in the first stack, and the TFM
% counterpart of each FRET image at the corresponding index of the second
% stack (.tif files). The script then loops through each stack and for each
% Fret image, it virtually separates the cell into different rings of
% increasing radius. It does so using concentric circles which have as
% their center the leading edge of the cell, which is assumed to be the
% point of highest traction on its counterpart TFM image. A customizable
% amount of rings can be used, each with a radius spaced linearly between
% the smallest ring (having radius 50 pixels) and the largest ring (having
% radius the max feret diameter of the binarized image of the cell * 1.5).
% The output of this script is a text file containing the total pixel
% intensity values encompassed within each concentric circle which can be
% used for further analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCLAIMER: This script assumes that the leading edge of a cell can be
% identified as the point of highest traction on the output of TFM analysis
% performed on the same image of said cell.
% If ever the code does not work for a new set of images, the
% problem might be the color mask function that is at the end of the
% script. In this event, either remove or comment out the functions at the
% bottom of the script, open the 'colorThresholder' app with your new set
% of TFM images and create a thresholder that leaves only the highest
% traction areas, and export the function to the bottom of this script,
% make sure the name of the function corresponds to the one being called in
% the script and you're good to go
%% clean up
clear
clc

%% Basic Functions .tiff setup
%loads the Leica Image Format to matlab
lif_load = @ci_loadLif.ci_loadLif;

%This function saves a figure or single axes to one or more vector and/or
%bitmap file formats
export_fig = @export_fig.export_fig;
% load('path.mat', 'filename')
%prompting the user to choose a file
[f, pathname] = uigetfile({'*.tif','Tagged Image File Format (*.tif)'}, 'Select TIF file', mfilename);
SelectedFile = fullfile(pathname,f);
if ~f, warning('No file selected.'); return; end

% filename = [pathname f];
%save('output.mat', 'filename', '-mat');
warning('off', 'Images:initSize:adjustingMag');

%prompting user to choos the corrseponding TFM image
[f2, PathName2] = uigetfile({'*.tif','Tagged Image File Format (*.tif)'}, 'Select Fret Ratio file', mfilename);
SelectedFile2 = fullfile(PathName2,f2);

if ~f2, warning('No file selected.'); return; end

%% checking stack length matches
% extracting length of Fret Stack and then TFM stack
Info1 = imfinfo(SelectedFile);
Info2 = imfinfo(SelectedFile2);

if length(Info1) ~= length(Info2), warning('Stacks do not have corresponding size');return;end

%% opening a .txt file to output the sum of pixel intensities in each circle
MyResults = fopen("ConcentricCircleResults.txt", "w");

%% reading the stacks and incrementing through the slices

FretStack = imread(SelectedFile);
TFMStack = imread(SelectedFile2);
disp('Starting analysis')

for frame = 1 : length(Info1)
    
    FretIm = imread(SelectedFile, frame);
    FretIm = rgb2gray(FretIm); %converts to 8-bit image
    TractionIm = imread(SelectedFile2, frame);
    
    fprintf('Now processing image %2.0f ... \n', frame);
    
    %% thresholding Fret image
    
    %retrieving auto thresholding value of filtered FRET image
    ThreshVal = graythresh(imgaussfilt(FretIm,5));
    % creating the desired binary image
    IIT= imfill(imbinarize(imgaussfilt(FretIm,5), ThreshVal), 'holes');
    
    
    %% initializing the variables to build the circles
    
    %extracting centroid and feret diameter measurements from the binarized im
    %The options mean that this function outputs the centroid and feret
    %properties to a table
    FretObjects = regionprops('table', IIT,'centroid', 'MaxFeretProperties');
    
    %creating variables for the loop
    NumCircles = 7;
    
    % [MaxFeretFret,Index1] = max(FeretDiameterTableFret{:,1});
    [MaxFeretFret,Index1] = max(FretObjects{:,2});
    FretCenter = FretObjects{Index1,1};
    Radius = linspace(50, MaxFeretFret*1.5  , NumCircles);
    
    %% processing the TFM image
    
    %calling the function at the bottom of the script to create a mask of
    %traction image
    [BW, ImMask] = createMask(TractionIm);
    %finding the centroid of each object identified in the binary mask
    TFMobjects = regionprops('table',BW, 'centroid', 'MaxFeretProperties');
    f4 = figure;
    imshowpair(TractionIm,ImMask,'montage')
    
    %finding the maximum feret diameter object identified in the TFM object
    %which is supposedly the leading edge object
    [MaxFeretTFM, Index2] = max(TFMobjects{:,2});
    %setting the center of the concentric circles to the center of the TFM
    %object
    TFMcenter = TFMobjects{Index2,1};
    %Modifiable center of the circle
    % Center = [(TFMcenter(1) +FretCenter(1))/2, (TFMcenter(2) + FretCenter(2))/2];
    Center = TFMcenter;
    
    %% using masks to create the concentric rings and measure them
    
    %updating the index of the image that is being analyzed in output file
    fprintf(MyResults,'The following measurements were extracted from image %.0f of the stack :\n', frame);
    fprintf(MyResults, 'Total sum of pixel intensity values (0-255) :\n');
    %creating a loop for each circle and its mask
    %the loop determines how many circles are made creating 5 equally spaced-
    %out circles between 50 pixels of radius to, the max Feret diameter
    CM = colormap(jet(70));
    f1 = figure('Name', 'Different Rings of the image');
    for k = 1 : length(Radius)
        %creating a duplicate of the FRET image at each circle creation
        CloneIm = FretIm;
        %creating an array that allows to compute the circle's edge
        Theta = 0:pi/50:2*pi;
        %these compute the x and y coordinates of the circle using above theta
        xunit = Radius(k) * cos(Theta) + Center(1);
        yunit = Radius(k) * sin(Theta) + Center(2);
        %creating a mask that takes everything inside the circle that was just
        %created as foreground and zeros everything outside of it to create an
        %image of the same dimensions as FretIm
        FretSize = size(FretIm);
        ImMaskL = poly2mask(xunit,yunit,FretSize(1),FretSize(2));
        
        %if the circle drawn is the smallest, then it sets to 0 all the pixels
        %outside the radius of the smallest circle
        if Radius(k) == Radius(1)
            CloneIm(~ImMaskL) = 0;
            CloneIm(~IIT) = 0;
            %visualization of the circle
            subplot(2,4,k)
            imshow(CloneIm,CM)
            %else the pixels outside the circle just drawn are set to 0 and so are
            %the ones within the circle that is one radius smaller
        else
            %once again this takes all pixels from the fret image that are not
            %in the cicrcle and zeros them
            CloneIm(~ImMaskL)=0;
            %here we change the radius of the circle to the previous radius
            %size, we trace the circle, and create the circle mask
            xunit = Radius(k-1) * cos(Theta) + Center(1);
            yunit = Radius(k-1) * sin(Theta) + Center(2);
            ImMaskS = poly2mask(xunit,yunit,FretSize(1),FretSize(2));
            %this sets everything from the fret image that is inside the circle
            %or radius one size smaller to 0, that way everything outisde the
            %big circle and inside the small one have zero values
            CloneIm(ImMaskS) = 0;
            CloneIm(~IIT) = 0;
            
            subplot(2,4,k)
            imshow(CloneIm,CM)
        end
        
        
        %here, the images are left with only the ROI's having non-zero pixel
        %intensities, so I find all those pixels, and then I find the total sum
        %of intensities in the ROI
        [X1,Y1,Val] = find(medfilt2(CloneIm));
        Total = sum(Val);
        
        %Here I output the data to an output file called
        %'ConcentricCircleResults.txt'
        fprintf(MyResults,'\tWithin Circle %.0f\n', k);
        fprintf(MyResults, '\t%.0f\n', Total);
        
        %here I create an overflow bin where i measure the total intensity of
        %everything that falls outside the biggest ring
        %checking that we are at the last iteration of the circle
        if k == length(Radius)
            CloneIm2 = FretIm;
            CloneIm2(ImMaskL) = 0;
            CloneIm2(~IIT) = 0;
            subplot(2,4,k+1)
            imshow(CloneIm2,CM)
            %finding all non zero values outside the largest circle
            [X2,Y2,Val2] = find(medfilt2(CloneIm2));
            Total2 = sum(Val2);
            fprintf(MyResults,'\nIn the overflow bin\n');
            fprintf(MyResults, '\t%.0f\n', Total2 );
        end
        
    end
    
    fprintf(MyResults, 'Processing of slice %.0f complete.\n\n', frame);
end

disp('End of stack reached')
fclose(MyResults);
disp('Program done')




% **************************** FUNCTIONS BELOW ****************************
% *************************************************************************

%% function generated and exported from the colorThresholder matlab app

function [BW,maskedRGBImage] = createMask(RGB)
%createMask  Threshold RGB image using auto-generated code from colorThresholder app.
%  [BW,MASKEDRGBIMAGE] = createMask(RGB) thresholds image RGB using
%  auto-generated code from the colorThresholder app. The colorspace and
%  range for each channel of the colorspace were set within the app. The
%  segmentation mask is returned in BW, and a composite of the mask and
%  original RGB images is returned in maskedRGBImage.

% Auto-generated by colorThresholder app on 12-Aug-2020
%------------------------------------------------------


% Convert RGB image to chosen color space
I = RGB;

% Define thresholds for channel 1 based on histogram settings
channel1Min = 0.000;
channel1Max = 255.000;

% Define thresholds for channel 2 based on histogram settings
channel2Min = 0.000;
channel2Max = 255.000;

% Define thresholds for channel 3 based on histogram settings
channel3Min = 0.000;
channel3Max = 255.000;

% Create mask based on chosen histogram thresholds
sliderBW = (I(:,:,1) >= channel1Min ) & (I(:,:,1) <= channel1Max) & ...
    (I(:,:,2) >= channel2Min ) & (I(:,:,2) <= channel2Max) & ...
    (I(:,:,3) >= channel3Min ) & (I(:,:,3) <= channel3Max);

% Create mask based on selected regions of interest on point cloud projection
I = double(I);
[m,n,~] = size(I);
polyBW = false([m,n]);
I = reshape(I,[m*n 3]);

% Project 3D data into 2D projected view from current camera view point within app
J = rotateColorSpace(I);

% Apply polygons drawn on point cloud in app
polyBW = applyPolygons(J,polyBW);

% Combine both masks
BW = sliderBW & polyBW;

% Initialize output masked image based on input image.
maskedRGBImage = RGB;

% Set background pixels where BW is false to zero.
maskedRGBImage(repmat(~BW,[1 1 3])) = 0;

end

function J = rotateColorSpace(I)

% Translate the data to the mean of the current image within app
shiftVec = [6.981305 31.350006 216.665937];
I = I - shiftVec;
I = [I ones(size(I,1),1)]';

% Apply transformation matrix
tMat = [-0.001626 -0.001338 0.000000 0.699473;
    0.001206 -0.001803 0.000000 0.103623;
    0.000000 0.000000 -0.002308 9.160254;
    0.000000 0.000000 0.000000 1.000000];

J = (tMat*I)';
end

function polyBW = applyPolygons(J,polyBW)

% Define each manually generated ROI
hPoints(1).data = [0.484197 0.481865;
    0.233390 0.486250;
    0.205522 0.350302;
    0.238963 0.240667;
    0.470263 0.196813;
    0.629108 0.245052;
    0.617961 0.442396];

% Iteratively apply each ROI
for ii = 1:length(hPoints)
    if size(hPoints(ii).data,1) > 2
        in = inpolygon(J(:,1),J(:,2),hPoints(ii).data(:,1),hPoints(ii).data(:,2));
        in = reshape(in,size(polyBW));
        polyBW = polyBW | in;
    end
end

end