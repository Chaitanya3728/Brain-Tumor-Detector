clc;
close all;
clear all;

%% Input image
s=imread('brain tumor.jfif');                              % Read the input image
s=imresize(s, [256, 256]);                                 % Resize image to 256x256 pixels
figure; imshow(s); title('Input image','FontSize',10);     % Display the resized image

%% Filtering
%Set the number of iterations, time step, and diffusion coefficient
num_iter = 10;
delta_t = 1/7;
kappa = 15;
option = 2;
disp('Preprocessing image please wait . . .');             % Display a message for image preprocessing 
inp = anisodiff(s,num_iter,delta_t,kappa,option);          % Apply anisotropic diffusion filter to the image
inp = uint8(inp);                                          % Convert the filtered image to uint8 data type
inp=imresize(inp,[256,256]);                               % Resize image to 256x256 pixels
% Convert RGB image to grayscale if it has multiple color channels
if size(inp,3)>1
    inp=rgb2gray(inp);
end
figure; imshow(inp); title('Filtered image','FontSize',10); % Display the preprocessed image

%% Thresholding
if size(s, 3) == 3                                         % Check if the image is RGB
    grayImage = rgb2gray(s);                               % If yes, convert it to grayscale using rgb2gray
else
    grayImage = s;                                         
end
threshold = 0.8;                                           % Set the threshold value to 0.8
binaryImage = imbinarize(grayImage, threshold);            % Convert the grayscale image to binary image using imbinarize

sout=imresize(inp,[256,256]);           % Resize image to 256x256 pixels
t0=mean(s(:));                          % Compute the mean value
th=t0+((max(inp(:))+min(inp(:)))./2);   % Compute a threshold value
for i=1:1:size(inp,1)                   % Loop through the rows of input image
    for j=1:1:size(inp,2)               % Loop through the columns of input image
        if inp(i,j)>th                  % If pixel value > threshold value
            sout(i,j)=1;                % Set the corresponding pixel value in the output image to 1
        else                            
            sout(i,j)=0;                % Set the corresponding pixel value in the output image to 0
        end
    end
end

%% Morphological Operation
label=bwlabel(sout);                    % Labeling connected components in the binary image
stats=regionprops(logical(sout),label,'Solidity','Area','BoundingBox');     % Calculating properties of the labeled regions
density=[stats.Solidity];               % Extracting solidity values 
area=[stats.Area];                      % Extracting area values
high_dense_area=density>0.7;            % Finding regions with high solidity values
max_area=max(area(high_dense_area));    % Finding maximum area among the high solidity regions
tumor_label=find(area==max_area);       % Finding label of the region with maximum area
tumor=ismember(label,tumor_label);      % Creating a binary image of the tumor region

if max_area>200                         % Checking if the tumor area > threshold
   figure; imshow(tumor); title('tumor alone','FontSize',10);
else
    h = msgbox('No Tumor detected !!','status');  % Displaying a message box if no tumor is found
    %disp('no tumor');
    return;
end
            
%% Bounding box
box = stats(tumor_label);                         % Get the bounding box coordinates for the tumor
wantedBox = box.BoundingBox;                      % Saving coordinates in a variable
figure, imshow(inp); title('Bounding Box','FontSize',10); % display
hold on;                                          % Hold the current plot
rectangle('Position',wantedBox,'EdgeColor','y');  % Draw a rectangle on the image using the bounding box coordinates
hold off;                                         % Release the current plot


%% Getting Tumor Outline - image filling, eroding, subtracting
%erosion the walls by a few pixels
dilationAmount = 5;                               % sets the dilation amount for morphological operations
rad = floor(dilationAmount);                      % calculates the radius for morphological operations
[r,c] = size(tumor);                              % getting size of the input image
filledImage = imfill(tumor, 'holes');             % fills the holes in the input image

% for each pixel in the image, creates a square structuring element of
% radius 'rad' centered at the pixel, and erodes the filled image with
% that structuring element to get the eroded output image
for i=1:r
   for j=1:c
       x1=i-rad;
       x2=i+rad;
       y1=j-rad;
       y2=j+rad;
       if x1<1
           x1=1;
       end
       if x2>r
           x2=r;
       end
       if y1<1
           y1=1;
       end
       if y2>c
           y2=c;
       end
       erodedImage(i,j) = min(min(filledImage(x1:x2,y1:y2)));
   end
end
figure, imshow(erodedImage); title('eroded image','FontSize',10); % display eroded image

%% Subtracting eroded image from original BW image
tumorOutline=tumor;                                 % assigning binary tumor image
tumorOutline(erodedImage)=0;                        % Remove eroded portion from tumor outline
figure; imshow(tumorOutline); title('Tumor Outline','FontSize',10); % displaying tumor outlined image


%% Inserting the outline in filtered image in red color
rgb = inp(:,:,[1 1 1]);    % converts the input image to RGB format

red = rgb(:,:,1);          % extracts red channel
red(tumorOutline)=255;     % sets pixel values of tumor outline to white in the red channel
green = rgb(:,:,2);        % extracts green channel
green(tumorOutline)=0;     % sets pixel values of tumor outline to black in the green channel
blue = rgb(:,:,3);         % extracts blue channel
blue(tumorOutline)=0;      % sets pixel values of tumor outline to black in the blue channel

tumorOutlineInserted(:,:,1) = red;   % assigns the modified red channel to 1st layer
tumorOutlineInserted(:,:,2) = green; % assigns the modified green channel to 2nd layer
tumorOutlineInserted(:,:,3) = blue;  % assigns the modified blue channel to 3rd third layer

figure, imshow(tumorOutlineInserted); title('Detected Tumor','FontSize',10); % displaying tumor detected image

%% Display Altogether
figure,
subplot(2,3,1);imshow(s);title('Input image','FontSize',10); 
subplot(2,3,2);imshow(inp);title('Filtered image','FontSize',10);
subplot(2,3,3);imshow(inp);title('Bounding Box','FontSize',10); 
hold on;rectangle('Position',wantedBox,'EdgeColor','y');hold off; % Superimposing a rectangle on the filtered image at the wanted box location
subplot(2,3,4);imshow(tumorOutlineInserted);title('Detected Tumor','FontSize',10);
subplot(2,3,5);imshow(tumor);title('tumor alone','FontSize',10);
subplot(2,3,6); imshow(binaryImage); title('Thresholded Image','Fontsize' ,10);

%% Conclusion
%In this project, we performed the segmentation of a brain tumor image 
% using image processing techniques. 
%The input image is read and then resized to 256x256 pixels. 
%If the image is RGB, it is converted to grayscale. 
%Then, an anisotropic diffusion filter is applied to the image to 
% reduce noise and enhance edges. 
%The filtered image is thresholded to convert it to a binary image. 
%The binary image is then subjected to morphological operations to 
% label connected components and calculate their properties. 
%The region with the highest area is selected as the tumor region. 
%If the tumor area is greater than a threshold of 200 pixels, the 
% tumor alone is displayed along with its bounding box. 
%The code then erodes the tumor region to get a smoother outline and 
% subtracts it from the original binary image to obtain the tumor outline. 
%The tumor outline is then displayed. 
%If no tumor is detected, message box is displayed saying 'no tumor detected'.