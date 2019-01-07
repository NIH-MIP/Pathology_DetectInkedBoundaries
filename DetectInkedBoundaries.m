%% DetectInkedBoundaries
%
% Stephanie Harmon
% stephanie.harmon@nih.gov
% Leidos Biomedical Research
% Molecular Imaging Branch, National Cancer Institute
% National Institutes of Health
% December 2018
% 
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%% DESCRIPTION
%
% This program takes in digital pathology (.czi/.svs) files with inked
% markings, automatically detects boundaries of ink, provides user
% interface for editing of detected boundaries, and writes to annotation
% file (.cz/.xml). 
%
% (Optional) Additionally, users can provide corresponding un-inked
% digital imaging, which will be automatically registered to inked imaging
% and annotations saved
%
% (Optional - NEW 1/7/2019) users can elect to manually outline inked borders, 
% which are then automatically registered to un-inked imaging (if provided) and
% annotations automatically saved. This may come in use if the 'auto'
% option fails to accurately detect boundaries
%
% Annotations will be saved to the highest resolution image contained in
% the digital stack (i.e. 40x if full digital image) provided at input.
%
% This program utilizes other freely available toolboxes, see dependencies
%
% General workflow:
%       1. read in image stack, find lowest magnification ratio (smallest
%           sampling of full-res image)
%       2. Option 'auto' ink detection/outlining
%               2a. deconvolve RGB H&E image using Khan et al Classifier
%               2b. Make a mask of tissue sample using H channel
%               2c. Use k-means clustering to detect inked markings from 
%                   remaining channels within tissue sample
%               2d. Grow-shrink morhological operations to determine ROIs
%               2e. Initiate user interface for accept/reject of all 
%                   proposed ROIs
%               2f. Initiate user interface to allow editing of accepted ROIs
%          Option 'manual' ink outlining
%               3a. user prompted for number of ROIs
%               3b. initiate user interface for manual outlining 
%               3c. after each ROI is outlined, user should close figure
%       3. Correlation-based image registration of inked and un-inked
%       4. Write annotation file 
%           (.xml if supplied .svs or .cz if supplied .czi)
%
%% DEPENDENCIES
%
% - local installation of MATLAB 2018b
%
% - MATLAB Image Processing Toolbox
%
% - Bio-Formats toolbox for MATLAB
%     https://docs.openmicroscopy.org/bio-formats/5.3.4/users/matlab/index.html
%
% - color decon toolbox from Khan et al 
%     https://github.com/lun5/color-deconvolution/tree/master/stain_normalisation_toolbox
%
% - polygon decimate function from 
%     https://www.mathworks.com/matlabcentral/fileexchange/34639-decimate-polygon   
%
%% USAGE
%   
%   DetectInkedBoundaries('--Marked','/example/path/to/marked.svs','--Unmarked', '/example/path/to/unmarked.svs','--Method','auto')
%
%   input:
%           1. (required) digital image in CZI or SVS format with inked markings
%           2. (optional) digital image of specimen with removed markings,
%                   if provided, will register low res digital images and register digital
%                   markings to 'clean' specimen
%           3. (optional) method of annotation. default = 'auto'. users who
%                   wish to outline images themselves should use 'manual'
%
%% TIPS AND TRICKS 
%
%  - Points can be added and deleted from ROIs after accept/reject stage, this stage is meant to exclude 
%       false positives arising from other inked notes/markings (i.e. "EPE", arrows, etc).
%
%  - If several proposed ROIs are part of the same region, simply add more points to join ROIs as part of
%       Workflow Step 7. A closing operation is performed that will merge all overlapping ROIs.
%
%  - If an ROI is open-ended (i.e. not closed circle), the tool will encompass only ink. To expand to other 
%       regions, all add and/or drag points to encompass the full ROI area
%
%% main
function [] = DetectInkedBoundaries(varargin)
   
    % DIRECTORY
    dir_id = find(strcmpi(varargin,'--Marked'));  
    if(size(dir_id) > 0)
        czi_m = varargin{dir_id+1};
        [filepath_m,name_m,ext_m] = fileparts(czi_m);
        if(~exist(czi_m)) % check to make sure it exists
            disp(['ERROR: file path does not exist']);
            return
        end
    else % user must provide parent directory
        disp('ERROR: required to input marked image path');
        return
    end
    
    % FILE TYPE
    typ_id = find(strcmpi(varargin,'--Unmarked'));
    if(size(typ_id) > 0)
        imgs_in = 2;
        %save path here
        czi_u = varargin{typ_id+1};
        [filepath_u,name_u,ext_u] = fileparts(czi_u);
    else
        imgs_in = 1;
    end 
    
    % MANUAL VS AUTO CONTOURS
    meth_id = find(strcmpi(varargin,'--Method'));
    if(size(meth_id) > 0)
        method_in = varargin{meth_id+1};
        %save path here
        if(strcmpi(method_in,'manual'))
            method = 'manual';
        else
            method = 'auto';
        end  
    else
        method = 'auto';
    end 
    

 % read in inked (marked) image : referred to with subscript "m" throughout remainder of code
    % grab original image - detemine if its a full stack or if its a scout 
        disp ('reading in original inked image')
        mData = bfGetReader(czi_m);
        mMeta = mData.getMetadataStore(); 
        imgLevels_m = mMeta.getImageCount;
        
    % if SCOUT (single img stack)
        if(imgLevels_m == 1)
            convertFactor = mData.getBitsPerPixel;
            RefLevel_m = 1;
            mIMG = readinIMG(mData, RefLevel_m);
            % sometimes the scout images are not correct base, we expect unit16
            if(convertFactor ~= 8) 
                mIMG = uint8(mIMG./(2^convertFactor/2^8));
            end
            mlevels = 'single';
            levelInfo_m = {};
            levelInfo_m{1,1} = eval(mMeta.getPixelsSizeX(0));
            levelInfo_m{1,2} = eval(mMeta.getPixelsSizeY(0));
            levelInfo_m{1,3} = 1; % wratio
            levelInfo_m{1,4} = 1; % hratio
            levelInfo_m{1,5} = 1; 
            max_ratio = 1;
        else
            levelInfo_m = {};
            levelInfo_m{1,1} = eval(mMeta.getPixelsSizeX(0));
            levelInfo_m{1,2} = eval(mMeta.getPixelsSizeY(0));
            max_ratio = 1;
            % identify magnification structure within image file
            for i=1:imgLevels_m
                levelInfo_m{i,1} = eval(mMeta.getPixelsSizeX(i-1));
                levelInfo_m{i,2} = eval(mMeta.getPixelsSizeY(i-1));
                levelInfo_m{i,3} = levelInfo_m{i,1}/levelInfo_m{1,1}; % wratio
                levelInfo_m{i,4} = levelInfo_m{i,2}/levelInfo_m{1,2}; % hratio
                if(1/str2num(sprintf('%0.4f',levelInfo_m{i,3})) == 1/str2num(sprintf('%0.4f',levelInfo_m{i,4}))) % sub of base
                    levelInfo_m{i,5} = int2str(round(1/str2num(sprintf('%0.4f',levelInfo_m{i,3})))); % mag ratio
                    if(str2num(levelInfo_m{i,5}) > max_ratio)
                        max_ratio = levelInfo_m{i,5};
                    end
                else
                    levelInfo_m{i,5} = 'NA'; %scout image or label image
                end
            end
            % we will determine RefLevel_m later, depending whether or not
            % input includes unmarked image
            mlevels = 'multi';
        end
        
        
        switch imgs_in
            case 1
                switch mlevels
                    case 'single'
                        if(strcmpi(method,'manual'))
                            boundaryLocation = manualContour(mIMG);
                        end
                        if(strcmpi(method,'auto'))
                            boundaryLocation = processIMG(mIMG);
                        end
                        boundaries_m = boundaryLocation; %no scaling necessary
                        writeXML(boundaries_m,[filepath_m filesep name_m], ext_m) 
                        %save
                    case 'multi'
                        % select lowest level (highest ratio)
                        RefLevel_m = find(strcmpi(levelInfo_m(:,5),max_ratio));
                        mIMG = readinIMG(mMeta, mData, RefLevel_m);
                        if(strcmpi(method,'manual'))
                            boundaryLocation = manualContour(mIMG);
                        end
                        if(strcmpi(method,'auto'))
                            boundaryLocation = processIMG(mIMG);
                        end
                        %scale and save
                        boundaries_m = scalebounds(boundaryLocation, levelInfo_m(RefLevel_m,:));
                        writeXML(boundaries_m,[filepath_m filesep name_m], ext_m) 
                end
            case 2
                switch mlevels
                    case 'single'
                        [RefLevel_m, RefLevel_u, uData, uMeta, levelInfo_u] = parseMatchedImage(czi_u, mlevels, levelInfo_m);
                        uIMG = readinIMG(uData, RefLevel_u);
                        if(strcmpi(method,'manual'))
                            boundaryLocation = manualContour(mIMG);
                        end
                        if(strcmpi(method,'auto'))
                            boundaryLocation = processIMG(mIMG);
                        end
                        boundaries_reg = registerImages(mIMG,uIMG,boundaryLocation,mlevels);
                        boundaries_m = boundaryLocation; %no scaling necessary
                        boundaries_u = scalebounds(boundaries_reg, levelInfo_u(RefLevel_u,:));
                        writeXML(boundaries_m,[filepath_m filesep name_m], ext_m) 
                        writeXML(boundaries_u,[filepath_u filesep name_u], ext_u) 
                    case 'multi'
                        [RefLevel_m, RefLevel_u, uData, uMeta, levelInfo_u] = parseMatchedImage(czi_u, mlevels, levelInfo_m);
                        mIMG = readinIMG(mData, RefLevel_m);
                        uIMG = readinIMG(uData, RefLevel_u);
                        if(strcmpi(method,'manual'))
                            boundaryLocation = manualContour(mIMG);
                        end
                        if(strcmpi(method,'auto'))
                            boundaryLocation = processIMG(mIMG);
                        end
                        boundaries_reg = registerImages(mIMG,uIMG,boundaryLocation,mlevels);
                        boundaries_m = scalebounds(boundaryLocation, levelInfo_m(RefLevel_m,:));
                        boundaries_u = scalebounds(boundaries_reg, levelInfo_u(RefLevel_u,:));
                        writeXML(boundaries_m,[filepath_m filesep name_m], ext_m) 
                        writeXML(boundaries_u,[filepath_u filesep name_u], ext_u) 
                end
        end
        
end        

%%  if an un-inked image is provided, read in. referred to with subscript "u" throughout remainder of code
function [RefLevel_m, RefLevel_u, uData, uMeta, levelInfo_u] = parseMatchedImage(czi_u, mlevels, levelInfo_m)
    % grab unmarked image - detemine if its a full stack or if its (throw error if scout)
        disp ('reading in un-inked image')
        uData = bfGetReader(czi_u);
        uMeta = uData.getMetadataStore(); 
        imgLevels_u = uMeta.getImageCount;
        
        if(imgLevels_u == 1)
            % add message box here eventually warning user they selected a
            % file that only has 1 level, making registration difficult
            % and/or if its a scout then cant map to full res
            answer = questdlg('The provided un-inked image only contains a single channel. Running this program on scout will not map annotations to full resolution image', ...
            'WARNING', ...
            'Proceed, this is not a scout image','Quit','Proceed, this is not a scout image');            
            switch answer
                case 'Proceed, this is not a scout image'
                    levelInfo_u = {};
                    levelInfo_u{1,1} = eval(uMeta.getPixelsSizeX(0));
                    levelInfo_u{1,2} = eval(uMeta.getPixelsSizeY(0));
                    RefLevel_u = 1;
                case 'Quit'
                    return
            end
        else
            levelInfo_u = {};
            levelInfo_u{1,1} = eval(uMeta.getPixelsSizeX(0));
            levelInfo_u{1,2} = eval(uMeta.getPixelsSizeY(0));

            for i=1:imgLevels_u
                levelInfo_u{i,1} = eval(uMeta.getPixelsSizeX(i-1));
                levelInfo_u{i,2} = eval(uMeta.getPixelsSizeY(i-1));
                levelInfo_u{i,3} = levelInfo_u{i,1}/levelInfo_u{1,1}; % wratio
                levelInfo_u{i,4} = levelInfo_u{i,2}/levelInfo_u{1,2}; % hratio
                if(1/str2num(sprintf('%0.4f',levelInfo_u{i,3})) == 1/str2num(sprintf('%0.4f',levelInfo_u{i,4}))) % sub of base
                    levelInfo_u{i,5} = int2str(round(1/str2num(sprintf('%0.4f',levelInfo_u{i,3})))); % mag ratio
                    if(size(find(strcmpi(levelInfo_u{i,5}, levelInfo_m(:,5))),1)>0)
                        levelInfo_u{i,6} = 'Y'; %has matching mag ratio to marked image
                    end
                else
                    levelInfo_u{i,5} = 'NA';
                end
            end

        
            switch mlevels

                % if marked image only has 1 level, find unmarked of most similar size
                % else if both are multi-layer, find matching layers
                case 'single'
                    RefLevel_m = 1;
                    matchOpt_u = levelInfo_u(~strcmpi(levelInfo_u(:,5),'NA'),:);
                    [valu,locu] = min(abs(1-cellfun(@(x) mod(x/levelInfo_m{1,1},levelInfo_m{1,1}), matchOpt_u(:,1))));
                    [RefLevel_u, valR] = find(cellfun(@(x) isequal(x,matchOpt_u{locu,1}),levelInfo_u(:,1)));

                case 'multi' %work on the highest level possible to avoid overloading the system
                    matchOpt_u = levelInfo_u(strcmpi(levelInfo_u(:,6),'Y'),:);
                    RefLevel_u = find(strcmpi(levelInfo_u(:,5),num2str(max(cellfun(@str2num, matchOpt_u(:,5))))));
                    RefLevel_m = find(strcmpi(levelInfo_m(:,5),num2str(max(cellfun(@str2num, matchOpt_u(:,5))))));
            end

        end
end


%% read in images at pre-specified level
%uses bioformats toolbox
function [lvlIMG] = readinIMG(i_data, i_level)
    i_data.setSeries(i_level-1);
    lvl1 = bfGetPlane(i_data,1); %R
    lvl2 = bfGetPlane(i_data,2); %G
    lvl3 = bfGetPlane(i_data,3); %B
    lvlIMG(:,:,1) = lvl1; 
    lvlIMG(:,:,2) = lvl2; 
    lvlIMG(:,:,3) = lvl3;
end
    

%% color decon and contour generation
function [boundaryLocation] = processIMG(mIMG)         

    %first deconvolve using random forest classifier to identify stains
    %uses Khan et al toolbox
    [M, Labels] = EstUsingSCD(mIMG);
    stains = Deconvolve(mIMG, M, 0);
    [s1, s2, s3] = PseudoColourStains(stains, M);

    %define large tissue area to avoid smaller usless points
    tiss_bw = rgb2gray(s1);
    tiss_label = imsegkmeans(tiss_bw,2,'NumAttempts',3);
    tiss_mask = zeros(size(tiss_bw));
    if(mean(tiss_bw(tiss_label == 1)) > mean(tiss_bw(tiss_label == 2)))
        tiss_mask(tiss_label == 2) = 1;
    else
        tiss_mask(tiss_label == 1) = 1;
    end
    tiss_mask = imclearborder(tiss_mask); % clear markings along borders
    tiss_mask = bwareaopen(tiss_mask,1000); 
    tiss_mask = imfill(tiss_mask,'holes'); % fill holes to connect regions 
    se = strel('sphere',5); 
    tiss_mask = imdilate(tiss_mask,se); 
    tiss_mask = imfill(tiss_mask,'holes'); % fill holes to connect regions 
    tiss_mask = imerode(tiss_mask,se); % shrink with same condition
    
    % markings are typically darker than h&e and dominate 2nd and 34d
    % channels, we will agregate from both
    I = rgb2gray(s2); %s2
    pixel_labels = imsegkmeans(I,2,'NumAttempts',3);
    mask_init = zeros(size(I)); % create initial mask - this will be choppy
    mask_init(pixel_labels>1) = 1;
    
    % take blue too and then join
    I = rgb2gray(s3); %s2
    pixel_labels = imsegkmeans(I,2,'NumAttempts',3);
    mask_init(pixel_labels>1) = 1;
    
    %grow/shrink to smooth edges
    mask_init = imclearborder(mask_init); % clear markings along borders
    mask_init = bwareaopen(mask_init,50); % find regions >50 pixels
    mask_init = imfill(mask_init,'holes'); % fill holes to connect regions
    se = strel('sphere',5);  % artifically grow, otherwise too many sections
    mask_init = imdilate(mask_init,se);
    mask_init = imclearborder(mask_init); % clear markings along borders 
    mask_init = imfill(mask_init,'holes'); % fill holes again
    mask_init = imclose(mask_init, se); % fill with same shape as dilation
    se = strel('sphere',5);  % now shrink
    mask_init = imerode(mask_init,se); % shrink with same condition
    boundaries = bwboundaries(mask_init); % identify boundaries of objects
    
    %only keep ones that are touching the rough estimate of tissue outline
    % this is to avoid the user going through too many
    boundaries_init = {};
    nb = 1;
    for j = 1:size(boundaries,1)
        b = boundaries{j};
        if(length(find(ismember(sub2ind(size(mask_init),b(:,1),b(:,2)), find(tiss_mask>0)))>0)>0)
            boundaries_init{nb,1} = b;
            nb = nb + 1;
        end
    end
        
    % select annotations for final contours
    boundaries_final = {};
    nb = 1;
    f = figure;
    for k=1:size(boundaries_init,1)
       b = boundaries_init{k};
       imshow(mIMG)
       hold on;
       plot(b(:,2),b(:,1),'y','LineWidth',3);
       answer = questdlg('Accept annotation?', ...
            'ROI Selection', ...
            'Yes','No','Yes');
        % Handle response
        switch answer
            case 'Yes'
                boundaries_final{nb,1} = b;
                nb = nb+1;
            case 'No'
        end
    end
    close(f)
    
    %allow freehand editing
        f = figure;
        imshow(mIMG, []);
        for ind = 1:numel(boundaries_final)
            % Convert to x,y order.
            pos = boundaries_final{ind};
            pos = fliplr(pos);
            %decimate to 5% of points
            [pos_out,i_rem,CI]=DecimatePoly(pos,[0.10 2],false);
            %only display half of points
            pos_pts = false(size(pos_out,1),1);
            pos_pts(1:round(0.25*size(pos_out,1))) = true;
            pos_pts = pos_pts(randperm(size(pos_pts,1)));
            % Create a freehand ROI.
             drawfreehand('Position', pos_out, 'Waypoints', pos_pts, 'Smoothing',0);
        end
        
        set(gcf, 'CloseRequestFcn', ' uiresume(gcbf); set(gcf,''Visible'',''Off'');')
        uiwait(gcf)

        hfhs = findobj(gcf, 'Type', 'images.roi.Freehand');
     
    %find final boundaries, merging overlapping regions    
        boundaryLocation = [];
        mask_all = zeros(size(mIMG,1),size(mIMG,2));
        for ind = 1:numel(hfhs)
            b_pts = hfhs(ind).Position;
            mask_b = poly2mask(b_pts(:,1),b_pts(:,2),size(mIMG,1),size(mIMG,2));
            mask_all(find(mask_b>0)) = 1;
        end
        
        se = strel('sphere',5);  % artifically grow, otherwise too many sections
        mask_all = imdilate(mask_all,se);
        mask_all = imfill(mask_all,'holes'); % fill holes again
        se = strel('line',25,45);  %line object to bring detected edges from outside of ink to inside
        mask_all = imerode(mask_all,se); %shrink to fall within boundaries
        
        boundaryLocation = bwboundaries(mask_all); % identify boundaries of objects
        close(f)

end

%% color decon and contour generation
function [boundaryLocation] = manualContour(mIMG) 

    num_rois = inputdlg('Enter number of ROIs:','Input',[1 25],{'1'});
    N_roi = str2num(num_rois{1});
    if(N_roi > 1)
        g = msgbox('close image window after each ROI');
    end

    %allow freehand editing
    f = figure;
    imshow(mIMG, []);
    
    all_pts = [];    
    for roi_i = 1:N_roi 
        drawfreehand;
        
        set(gcf, 'CloseRequestFcn', ' uiresume(gcbf); set(gcf,''Visible'',''Off'');')
        uiwait(gcf)

        hfhs = findobj(gcf, 'Type', 'images.roi.Freehand');
        i_pts = hfhs.Position;
        all_pts = cat(1,all_pts,{i_pts});
        clear hfhs
    end
    
    %find final boundaries, merging overlapping regions    
        boundaryLocation = [];
        mask_all = zeros(size(mIMG,1),size(mIMG,2));
        for ind = 1:size(all_pts,1)
            if(size(all_pts{ind})>0)
                b_pts = all_pts{ind};
                mask_b = poly2mask(b_pts(:,1),b_pts(:,2),size(mIMG,1),size(mIMG,2));
                mask_all(find(mask_b>0)) = 1;
            end
        end

        boundaryLocation = bwboundaries(mask_all); % identify boundaries of objects
        close(f)
        close(g)
end



%% registration
%Register grayscale images

function [boundaries_reg] = registerImages(mIMG,uIMG,boundaryLocation, mlevels)
 
    % Convert RGB images to grayscale
    uIMG_reg = rgb2gray(uIMG);
    mIMG_reg = rgb2gray(mIMG);

    % Default spatial referencing objects
    fixedRefObj = imref2d(size(uIMG_reg));
    movingRefObj = imref2d(size(mIMG_reg));

    % Default spatial referencing objects
    fixedRefObj = imref2d(size(uIMG_reg));
    movingRefObj = imref2d(size(mIMG_reg));    

    switch mlevels
        case 'single'
            % Phase correlation with similarity
            tform = imregcorr(mIMG_reg,movingRefObj,uIMG_reg,fixedRefObj,'transformtype','similarity','Window',true);
            MOVINGREG.Transformation = tform;
            MOVINGREG.RegisteredImage = imwarp(mIMG_reg, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
        
        case 'multi'
            % Phase correlation with translation
            tform = imregcorr(mIMG_reg,movingRefObj,uIMG_reg,fixedRefObj,'transformtype','translation','Window',true);
            MOVINGREG.Transformation = tform;
            MOVINGREG.RegisteredImage = imwarp(mIMG_reg, movingRefObj, tform, 'OutputView', fixedRefObj, 'SmoothEdges', true);
    end

    % Nonrigid registration
    [MOVINGREG.DisplacementField,MOVINGREG.RegisteredImage] = imregdemons(MOVINGREG.RegisteredImage,uIMG_reg,100,'AccumulatedFieldSmoothing',1.0,'PyramidLevels',3);

    % Store spatial referencing object
    MOVINGREG.SpatialRefObj = fixedRefObj;

    % apply registration to mask from boundaries
    mask_final = zeros(size(mIMG(:,:,1)));
    for k=1:size(boundaryLocation,1)
       b = boundaryLocation{k};
       mask_b = poly2mask(b(:,2),b(:,1),size(mIMG,1),size(mIMG,2));
       mask_final = mask_final + mask_b;
    end

    outputView = imref2d(size(uIMG));
    mask_reg = imwarp(mask_final, MOVINGREG.Transformation, 'OutputView', outputView);

    % make bwboundaries again
    boundaries_reg = bwboundaries(mask_reg);
      

end

%% scale to hgihest resolution provided in image stack
function [boundaries_scaled] = scalebounds(boundaryInfo, levelInfo)

    wratio = 1/levelInfo{1,3};
    hratio = 1/levelInfo{1,4};
    boundaries_scaled = {};
    for i = 1:length(boundaryInfo)
        b = boundaryInfo{i};
        b_new = zeros(size(b));
        b_new(:,1) = b(:,1).*wratio;
        b_new(:,2) = b(:,2).*hratio;
        boundaries_scaled{i,1} = b_new; 
        clear b_new
    end
    
end
%% write boundary points to annotation files
% default/dummy files are provided for .czi and .svs image formats
% each roi is given a unique region identifier 
 function [] = writeXML(boundaryLocation,saveName, saveType) 
    loadFilePath = which('DetectInkedBoundaries');
    [master_path,filename] = fileparts(loadFilePath);
    switch saveType
        case '.czi'
            load([master_path filesep 'cz_master.mat']);
            
            fileID = fopen([saveName '.txt'],'w');

            for mm = 1:size(cz_master.header,1)
                fprintf(fileID,'%s\r\n',cz_master.header(mm));
            end

            for jj = 1:length(boundaryLocation)
                %edit ROI id
                    roi_str = cz_master.roi;
                    roi_loc = find(contains(roi_str,'Bezier Id'));
                    roi_str(roi_loc) = strrep(roi_str(roi_loc),'NULL',int2str(jj));
                %edit ROI points
                    pts_loc = find(contains(roi_str,'<Points>'));
                    roi_pts = '';
                    for i = 1:length(boundaryLocation{jj})
                        pts_i = sprintf('%.0f,%.0f ',boundaryLocation{jj}(i,2),boundaryLocation{jj}(i,1));
                        roi_pts = [roi_pts pts_i];
                    end
                    roi_str(pts_loc) = strrep(roi_str(pts_loc),'NULL',roi_pts);
                %send to file    
                for mm = 1:size(roi_str,1)
                    fprintf(fileID,'%s\r\n',roi_str(mm));
                end
            end

            for mm = 1:size(cz_master.end,1)
                fprintf(fileID,'%s\r\n',cz_master.end(mm));
            end

            fclose(fileID);
            copyfile([saveName '.txt'],[saveName '.cz']);
    
        case '.svs'
            load([master_path filesep 'svs_master.mat']);
            
            fileID = fopen([saveName '.txt'],'w');  
            
            for mm = 1:size(svs_master.header,1)
                fprintf(fileID,'%s\r\n',svs_master.header(mm));
            end
            
            for jj = 1:length(boundaryLocation)
                %edit ROI id in header
                    roi_str = svs_master.roi_header;
                    roi_loc = find(contains(roi_str,'Region Id'));
                    roi_str(roi_loc) = strrep(roi_str(roi_loc),'NULL',int2str(jj));
                    for mm = 1:size(roi_str,1)
                        fprintf(fileID,'%s\r\n',roi_str(mm));
                    end
                %edit ROI points
                    for i = 1:length(boundaryLocation{jj})
                        ver_str = svs_master.roi_vertex;
                        ver_str = strrep(ver_str,'NUX',sprintf('%.0f',boundaryLocation{jj}(i,2)));
                        ver_str = strrep(ver_str,'NUY',sprintf('%.0f',boundaryLocation{jj}(i,1)));
                        fprintf(fileID,'%s\r\n',ver_str);
                    end
                %send ROI footer    
                    roi_str = svs_master.roi_footer;
                    for mm = 1:size(roi_str,1)
                        fprintf(fileID,'%s\r\n',roi_str(mm));
                    end
            end
            
            for mm = 1:size(svs_master.end,1)
                fprintf(fileID,'%s\r\n',svs_master.end(mm));
            end

            fclose(fileID);
            copyfile([saveName '.txt'],[saveName '.xml']);
            
    end
    delete([saveName '.txt']);
 end