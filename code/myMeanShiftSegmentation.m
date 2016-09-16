function myMeanShiftSegmentation( sigmacolor, sigmaspacial, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

inpImg = '..\data\baboonColor.png';
original = imread(inpImg);
channelCount = size(original ,3);
num_rows_original = size(original, 1);
num_cols_original = size (original, 2);
subsampled = zeros(floor(num_cols_original/5)+1, floor(num_rows_original/5)+1, channelCount);
subsampled = uint8(subsampled);
for i = 1:channelCount
    original1 = original(:, :, i);
    original1 = double(original1);
    gfilter = fspecial('gaussian', 5, 1);
    I1 = conv2(original1, gfilter, 'same');
    I1 = uint8(I1);
    I1 = I1(1:5:end, 1:5:end);
    subsampled(:, :, i) = I1;
end;
num_rows = size(subsampled, 1);
num_cols = size (subsampled, 2);

segmented = subsampled;

newCord = zeros(num_rows, num_cols, 2);

Imin = ones(channelCount);
Imin = Imin*255;
Imax = zeros(channelCount);
Imintotal = 255;
Imaxtotal = 0;
for y=1:num_rows
    for x=1:num_cols
        newCord(x,y,1) = x;
        newCord(x,y,2) = y;
        for i = 1:channelCount
            if(subsampled(x,y,i)>Imax(i)) Imax(i)=subsampled(x,y,i); end;
            if(subsampled(x,y,i)<Imin(i)) Imin(i) = subsampled(x,y,i); end;
        end;
    end;
end;

gspacial = fspecial('gaussian', (2*num_rows +1), sigmaspacial);

gcolor1 = gausswin((2*(Imax(1)-Imin(1))+1), (1/sigmacolor)); %fspecial('gaussian', (Imax(i)-Imin(i)), sigmacolor);
gcolor2 = gausswin((2*(Imax(2)-Imin(2))+1), (1/sigmacolor)); %fspecial('gaussian', (Imax(i)-Imin(i)), sigmacolor);
gcolor3 = gausswin((2*(Imax(3)-Imin(3))+1), (1/sigmacolor)); %fspecial('gaussian', (Imax(i)-Imin(i)), sigmacolor);

for n1=1:n
    %disp([' n is ',num2str(n1)]);
    for y=1:num_rows
        for x=1:num_cols
            %disp([' x is ',num2str(x), ' and y is ',num2str(y)]);
            %(x,y) is current pixel
            numx = 0;
            numy = 0;
            wsum = 0;
            wcolorFinal =1;
            wcolor = zeros(channelCount);
            for y1=1:num_rows
                for x1 = 1:num_cols
                    wspacial = gspacial( ((x1 - x)+(num_rows+1)) ,((y1 - y)+(num_rows+1)) );
                    %wcolor(x1, y1, i) = gspacial ( ( (subsampled(x1, y1, i) - subsampled(x,y, i))+((Imax(i)-Imin(i))/2) ), () )
                    %disp(newCord(x1,y1,1));
                    wcolor(1) = gcolor1( (subsampled(newCord(x1,y1,1), newCord(x1,y1,2), 1) - subsampled(newCord(x,y,1),newCord(x,y,2), 1)+((Imax(1)-Imin(1))+1)) ) ;
                    wcolor(2) = gcolor2( (subsampled(newCord(x1,y1,1), newCord(x1,y1,2), 2) - subsampled(newCord(x,y,1),newCord(x,y,2), 2)+((Imax(2)-Imin(2))+1)) ) ;
                    wcolor(3) = gcolor3( (subsampled(newCord(x1,y1,1), newCord(x1,y1,2), 3) - subsampled(newCord(x,y,1),newCord(x,y,2), 3)+((Imax(3)-Imin(3))+1)) ) ;
                    for i=1:channelCount
                        wcolorFinal = wcolorFinal*wcolor(i);
                    end;
                    numx = numx + (x1)*wspacial*wcolorFinal;
                    numy = numy + (y1)*wspacial*wcolorFinal;
                    wsum = wsum + wspacial*wcolorFinal;
                end;
            end;
            newCord(x,y,1) = uint8(numx/wsum);
            newCord(x, y, 2) = uint8(numy/wsum);
        end;
    end;
end;
for y=1:num_rows
    for x=1:num_cols
        for i=1:channelCount
            segmented(x, y, i) = subsampled(newCord(x,y,1), newCord(x,y,2), i);
        end;
    end;
end;
%{
%Mean-shift segmentation
%create a matrix Z such that 1st column is X, 2nd column is Y and 3rd
%column is intensity at X,Y
i=1;
for y=1:num_rows
    for x=1:num_cols
        z(i,1) = x;
        z(i,2) = y;
        z(i,3) = subsampled(x,y, 1);
        z(i,4) = subsampled(x,y, 2);
        z(i,5) = subsampled(x,y, 3);
        i = i+1;
    end;
end;

%create a function for comparing distance between x and each row of Z,
%where x is a row and Z is a matrix
% change this function to find nearest neighbor based on distance (col 1
% and col2) and intensity(col3) using sigmaspacial and sigmacolor. I am not
% able to do this
w = [0.4; 0.6];% weight matrix
mydist = @(x,Z)sqrt((bsxfun(@minus,x,Z).^2)*w);

%following statement returns a matrix IDX with 256*256 rows and 8 columns
%where values in ith row denotes the index of 8 points closest to point at
%index i using our distance function 'mydist'
IDX = knnsearch(double(segmented(:,:,1)), double(segmented(:,:,1)), 'Distance', mydist, 'k', 8);

%go on updating matrix IDX to get k rows at the last using the mean of
%nearest point(segmentation). I am not able to do this.
%}
disp(' Tuned parameter values are: ');
disp([' Sigmacolor is ',num2str(sigmacolor)]);
disp([' Sigmaspacial is ',num2str(sigmaspacial)]);
disp([' Number of iterations is ',num2str(n)]);

figure(1);
subplot(1, 2, 1);
imshow (subsampled); % phantom is a popular test image
title('Original subsampled');
daspect ([1 1 1]);
axis tight;

subplot(1, 2, 2);
imshow (segmented);
title('Segmented Image');
daspect ([1 1 1]);
axis tight;

imwrite(original,'..\images\baboonColorOriginal.png');
imwrite(segmented,'..\images\baboonColorSegmented.png');

end

