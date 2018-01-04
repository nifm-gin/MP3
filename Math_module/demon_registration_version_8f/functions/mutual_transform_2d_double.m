function [I1_TF,I2_TF]=mutual_transform_2d_double(I1,I2,I1_moved,I2_moved,range,HistogramBins,Hkernel, SampleDistance)
% This function transform one image modality into another.
% 
% [I1_TF,I2_TF]=mutual_transform_2d_double(I1,I2,I1_moved,I2_moved,range,HistogramBins,Hkernel, SampleDistance)
%
% Function is written by D. Kroon University of Twente(March 2009)


% Loop through all locations on which the histogram data is sampled
I1_TF=zeros(size(I1));
I2_TF=zeros(size(I2));
        
Histograms_sizes=ceil(size(I1)/SampleDistance);

halfKernel=(size(Hkernel,1)-1)/2;
halfSample=(SampleDistance-1)/2;

for hist_x=0:(Histograms_sizes(1)-1)
    for hist_y=0:(Histograms_sizes(2)-1)
        % Calculate the image coordinates of a location histogram
        image_x_center=hist_x*SampleDistance+halfSample;
        image_y_center=hist_y*SampleDistance+halfSample;

        % Check boundary conditions
        image_x_center(image_x_center>(size(I1,1)-1)) = size(I1,1)-1;
        image_y_center(image_y_center>(size(I1,2)-1)) = size(I1,2)-1;
        
        % Calculate the part of the kernel inside the image
        kernel_x_start=halfKernel-image_x_center; 
        kernel_y_start=halfKernel-image_y_center; 
        kernel_x_end=size(I1,1)-1+halfKernel-image_x_center; 
        kernel_y_end=size(I1,2)-1+halfKernel-image_y_center; 
       
        % Kernel Boundary conditions
        kernel_x_start(kernel_x_start<0)=0;
        kernel_y_start(kernel_y_start<0)=0;
        kernel_x_end(kernel_x_end>(size(Hkernel,1)-1))=(size(Hkernel,1)-1);
        kernel_y_end(kernel_y_end>(size(Hkernel,2)-1))=(size(Hkernel,2)-1);
                 
        % Image Region which uses this kernel
        image_x_start=image_x_center-halfSample; 
        image_x_end=image_x_center+halfSample; 
        image_y_start=image_y_center-halfSample; 
        image_y_end=image_y_center+halfSample; 
        
        % Boundary conditions on the image region
        image_x_start(image_x_start<0)=0;
        image_x_end(image_x_end>(size(I1,1)-1))=size(I1,1)-1;
        image_y_start(image_y_start<0)=0;
        image_y_end(image_y_end>(size(I1,2)-1))=size(I1,2)-1;
        
        % Make an histogram for I1 on I2 transformed and visa versa
        Histogram_I1=zeros([HistogramBins HistogramBins]);
        Histogram_I2=zeros([HistogramBins HistogramBins]);

        % Loop through the kernel
        for kernel_x=kernel_x_start:kernel_x_end,
            for kernel_y=kernel_y_start:kernel_y_end,
                 % The current image coordinates
                 image_x=image_x_center+kernel_x-halfKernel;
                 image_y=image_y_center+kernel_y-halfKernel;
                 
                 % The current histogram location
                 hist2d_x=round((HistogramBins-1)/(range(2)-range(1))*(I1(image_x+1,image_y+1)-range(1)));
                 hist2d_y=round((HistogramBins-1)/(range(2)-range(1))*(I2_moved(image_x+1,image_y+1)-range(1)));
                 
                 % Update the 2D histogram with the kernel value
                 Histogram_I1(hist2d_x+1,hist2d_y+1)=Histogram_I1(hist2d_x+1,hist2d_y+1)+Hkernel(kernel_x+1,kernel_y+1);
                 
                 % The current histogram location
                 hist2d_x=round((HistogramBins-1)/(range(2)-range(1))*(I1_moved(image_x+1,image_y+1)-range(1)));
                 hist2d_y=round((HistogramBins-1)/(range(2)-range(1))*(I2(image_x+1,image_y+1)-range(1)));
                 
                 % Update the 2D histogram with the kernel value
                 Histogram_I2(hist2d_x+1,hist2d_y+1)=Histogram_I2(hist2d_x+1,hist2d_y+1)+Hkernel(kernel_x+1,kernel_y+1);
            end
        end

        % Loop throught the image region which uses this kernel/histograms
        for image_x=image_x_start:image_x_end
            for image_y=image_y_start:image_y_end;
                
                   hist2d_x=round((HistogramBins-1)/(range(2)-range(1))*(I1(image_x+1,image_y+1)-range(1)));
                   [temp,hist2d_y]=max(Histogram_I1(hist2d_x+1,:)); hist2d_y=hist2d_y-1;
                   I1_TF(image_x+1,image_y+1)=hist2d_y/((HistogramBins-1)/(range(2)-range(1)))+range(1);

                   hist2d_y=round((HistogramBins-1)/(range(2)-range(1))*(I2(image_x+1,image_y+1)-range(1)));
                   [temp,hist2d_x]=max(Histogram_I2(:,hist2d_y+1)); hist2d_x=hist2d_x-1;
                   I2_TF(image_x+1,image_y+1)=hist2d_x/((HistogramBins-1)/(range(2)-range(1)))+range(1);
            end
        end        
        

    end
end











