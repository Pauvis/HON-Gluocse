%%%%%% Analysis workflow:
% cell_f_extraction_2P_doublecheckplanes_24012022a.m
% ExportPreprocesseddata_2PGlusensstim08022022.m
% Combine_2Pdata0802.m
% Analysis_combined_data_2P.m
% Locate2Pcells.m
% Script for 2P image stack processing after FIJI macro:
% 1) Open motioncorrected tiff stacks one by one and compare to "Stack"
% 2) Load existing labelled ROI data files ("ROI_list.xlsx" & "RoiSet.zip")
% 3) Get fluorescence values per ROI
%% Load Files. Files are to be organised as:
% Working directory
% Animal ID
% Experiments
% files with "...bin.tif", "ROI_list.xlsx" & "RoiSet.zip"
clear all
clc
plane_num = 6;
ignore_plane = 1;
% scenario 1: frame order inverted as per Mahesh system: plane_i [1,2,3,4,5,61] == plane_ii [6,5,4,3,2,1]
scenario1_list = [6,5,4,3,2,1];
% scenario1_list = [5,4,3,2,1,6];
% register_tiffs = 0; % if == 1 will motion correct, else just compare to tiff stacks already there

fid = 'C:\Users\pauliusv\Desktop\Test folder';
fid = 'P:\Alexander\Summer 2019\To_do\Sensory_stim_exp';
fps=5.148800329523 %INSS volume MAKE SURE TO UPDATE EVERY TIME!!!
halo1=6;
halo2=2;
% Get a list of all files and folders in this folder.
files = dir(fid);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..
% if ~exist('all_data','var')
%     all_data = [];
% end

for id_i = 1:length(subFolderNames)
    %     if isempty(all_data)
    %         all_data(1).IDs = subFolderNames{id_i};
    %     elseif  isempty(find(strcmpi({all_data.IDs},subFolderNames(id_i))))
    %         all_data(end+1).IDs = subFolderNames{id_i};
    %     end
    if ~strcmpi(subFolderNames{id_i},'Graphs')
        % find expnames:
        fid_exps = fullfile(fid,subFolderNames{id_i});
        files = dir(fid_exps);
        dirFlags = [files.isdir];
        subFolders = files(dirFlags);
        subFolderNames_exp = {subFolders(3:end).name};

        %     idx = find(strcmpi({all_data.IDs},subFolderNames(id_i)));

        for exp_i = 1:length(subFolderNames_exp)
            fid_exp = fullfile(fid_exps,subFolderNames_exp{exp_i});
            if exist(fullfile(fid_exp,'datanew.mat')) == 2 %&& 1==2
                disp('Exp done already');
            else
                %%
                disp(strcat('Running exp:',fid_exp));
                tif_stack_names = dir(fullfile(fid_exp,'*bin*.tif'));
                idx_to_del = find(contains({tif_stack_names.name},'STD')|contains({tif_stack_names.name},'bin.tif'));
                STD_stack = tif_stack_names(find(contains({tif_stack_names.name},'Stack')));
                tif_stack_names(idx_to_del) = [];

                %       all_data(idx).exps.(subFolderNames_exp{exp_i}).exp_name = tif_stack_name.name;

                % Load ROIs
                ROI_list=xlsread([fid_exp,'\ROI_list.xlsx']);
                [cvsROIs] = ReadImageJROI([fid_exp,'\RoiSet.zip']);
                % Conts
                CONTSfull=cvsROIs;
                for i=1:size(CONTSfull,2)
                    CONTSfull{1,i}=cvsROIs{1,i}.mnCoordinates;
                end

                % Read STD tiff stack
                tic;
                STD_stack_slices = tiff_read(fullfile(fid_exp,STD_stack.name),0);
                scenarionums = []; clear data CONTS cellconts npconts cell_f np_f mat np_temp outsx outsy

                for plane_i = 1:plane_num
                    if ignore_plane == scenario1_list(plane_i)
                    else
                        % Get ROI contours
                        %%% select ROIs from cont file
                        tempt=ROI_list(:,end-plane_i+1);
                        tempt(isnan(tempt))=[];
                        CONTS=CONTSfull(tempt);
                        %                   STD_stack_t_name = strcat(strrep(strrep(STD_stack.name,'.tif','_'),'Stack','STD'),num2str(plane_i),'.tif');
                        %                   %%
                        %                   STD_stack_slice_t = tiff_read(fullfile(fid_exp,STD_stack_t_name));
                        plane_stack_t_name = strcat(strrep(strrep(STD_stack.name,'.tif','_'),'Stack',''),num2str(plane_i),'.tif');
                        if exist(fullfile(fid_exp,plane_stack_t_name)) == 2
                        elseif exist(fullfile(fid_exp,strrep(plane_stack_t_name,'6.tif','61.tif'))) ==2
                            plane_stack_t_name = fullfile(fid_exp,strrep(plane_stack_t_name,'6.tif','61.tif'));
                        end
                        %%
                        plane_stack_slice_t = tiff_read(fullfile(fid_exp,plane_stack_t_name),0);
                        %Make STD projection
                        std_projection =  std(plane_stack_slice_t,[],3);


                        %% get all the pixels that contain a contour element
                        pixelsY=size(plane_stack_slice_t,1);
                        pixelsX=size(plane_stack_slice_t,2);
                        mat = zeros(pixelsY,pixelsX);
                        matx = zeros(pixelsY,pixelsX);
                        clear ROIs
                        for c = 1:length(CONTS)
                            ROIs(c).mat=zeros(pixelsY,pixelsX);
                            temp = CONTS{c};
                            Xind = ones(pixelsY,1)*(1:pixelsX);
                            Yind = (1:pixelsY)'*ones(1,pixelsX);
                            [IN ON] = inpolygon(Xind,Yind,temp(:,1),temp(:,2));
                            mat(find(IN))=1;%
                            matx(find(IN))=c;
                            ROIs(c).mat(find(matx==c))=1;
                        end


                        %
                        [c1,r1]=find(mat==0);
                        [c2,r2]=find(mat==1);
                        nppx=fliplr([r1,c1]);
                        cellpx=fliplr([r2,c2]);


                        % get the distance to all np pixels for each cell
                        clear np_CONTS;
                        %                 for i=1:length(CONTS)
                        %                     cellcenter=flip(floor(mean(CONTS{i})));
                        %                     pixdist=zeros(1,length(nppx));
                        %                     for j=1:length(nppx);
                        %                         clear temp;
                        %                         temp=distmat([cellcenter;nppx(j,:)]);
                        %                         pixdist(j)=temp(1,2);
                        %                     end
                        %                 end

                        %% alternatively get np conts by following cell outlines
                        cellconts=zeros(pixelsY,pixelsX);
                        npconts=zeros(pixelsY,pixelsX);
                        for i=1:length(CONTS);

                            temp=fliplr(ceil(CONTS{i}));
                            %ind1=find(temp(:,1)<0);
                            %temp=temp+1;
                            for j=1:length(temp);
                                cellconts(temp(j,1),temp(j,2))=i;
                            end
                            %find nearest nppx to cont
                            p(1,2)=0;
                            for j=1:length(temp);
                                pixel=temp(j,:);
                                for s=-halo1:halo1 %adjust here for halo size
                                    a1=find(nppx(:,1)==pixel(1)+s);
                                    a2=find(nppx(:,2)==pixel(2)+s);
                                    a3=[nppx(a1,:);nppx(a2,:)];
                                    sss=1;
                                    a4(1,2)=0;
                                    for ss=1:length(a3)
                                        if ((pixel(1)-a3(ss,1))^2+(pixel(2)-a3(ss,2))^2)^0.5 <halo1 %adjust here for halo size
                                            a4(sss,:)=a3(ss,:);
                                            sss=sss+1;

                                        end
                                    end
                                    p((size(p,1)+1):(size(p,1)+size(a4,1)),:)=a4;

                                    clear a1 a2 a3 a4
                                end
                            end

                            A=p;
                            A(find(A(:,1)==0),:)=[];
                            A(find(A(:,2)==0),:)=[];
                            for f=1:size(A,1)
                                npconts(A(f,1),A(f,2))=i;

                            end
                            clear A p

                            clear temp
                            [I(1,:),I(2,:)]=find(npconts==i);
                            np_CONTS(i)={I'};
                            clear I
                        end

                        %exclusion loop
                        for i=1:length(CONTS);

                            temp=fliplr(ceil(CONTS{i}));
                            %ind1=find(temp(:,1)<0);
                            %temp=temp+1;
                            for j=1:length(temp);
                                cellconts(temp(j,1),temp(j,2))=i;
                            end
                            %find nearest nppx to cont
                            p(1,2)=0;
                            for j=1:length(temp);
                                pixel=temp(j,:);
                                for s=-halo2:halo2 %set exclusion
                                    a1=find(nppx(:,1)==pixel(1)+s);
                                    a2=find(nppx(:,2)==pixel(2)+s);
                                    a3=[nppx(a1,:);nppx(a2,:)];
                                    sss=1;
                                    a4(1,2)=0;
                                    for ss=1:length(a3)
                                        if ((pixel(1)-a3(ss,1))^2+(pixel(2)-a3(ss,2))^2)^0.5 <halo2 %set exclusion
                                            a4(sss,:)=a3(ss,:);
                                            sss=sss+1;

                                        end
                                    end
                                    p((size(p,1)+1):(size(p,1)+size(a4,1)),:)=a4;

                                    clear a1 a2 a3 a4
                                end
                            end

                            A=p;
                            A(A(:,1)==0,:)=[];

                            for f=1:size(A,1)
                                npconts(A(f,1),A(f,2))=0;

                            end
                            clear A p

                            clear temp
                            [I(1,:),I(2,:)]=find(npconts==i);
                            np_CONTS(i)={I'};
                            clear I
                        end
                        %
                        %                 figure;
                        %                 imagesc(cellconts);
                        %                 figure;
                        %                 imagesc(npconts);

                        %now get the ROI mats for NPs
                        for c = 1:length(np_CONTS)
                            NPs(c).mat=zeros(pixelsY,pixelsX);
                            NPs(c).mat(find(npconts==c))=1;
                        end

                        %% extract the data from the main movie matrix
                        framen=size(plane_stack_slice_t,3);
                        cell_f=zeros(framen,length(CONTS));
                        np_f=zeros(framen,length(CONTS));
                        for i=1:length(CONTS);
                            % cell_temp=fliplr(floor(CONTS{i}));
                            % outsx=find(cell_temp(:,1)>pixelsY);
                            % outsy=find(cell_temp(:,2)>pixelsX);
                            % cell_temp(outsx,1)=pixelsY;
                            % cell_temp(outsy,2)=pixelsX;

                            %  [cell_temp(1,:),cell_temp(2,:)]=find(ROIs(i).mat==1);

                            np_temp=floor(np_CONTS{i});
                            outsx=find(np_temp(:,1)>pixelsY);
                            outsy=find(np_temp(:,2)>pixelsX);
                            np_temp(outsx,1)=pixelsY;
                            np_temp(outsy,2)=pixelsX;

                            for j=1:framen;
                                cell_f(j,i)=sum(sum(plane_stack_slice_t(:,:,j).*ROIs(i).mat))/length(find(ROIs(i).mat==1)); %multiplies by the ROI mat and sums then divides to get arithmetic mean.
                                np_f(j,i)=sum(sum(plane_stack_slice_t(:,:,j).*NPs(i).mat))/length(find(NPs(i).mat==1));
                                %         cell_f(j,i) = mean(mean(mean(Y(cell_temp(:,1),cell_temp(:,2),j))));
                                %         np_f(j,i) = mean(mean(mean(Y(np_temp(:,1),np_temp(:,2),j))));
                            end
                        end


                        %write
                        assignin('base','cell_f',cell_f);
                        assignin('base','np_f',np_f);
                        assignin('base','npconts',npconts);
                        assignin('base','cellconts',cellconts);
                        assignin('base','mat',mat);
                        %% plot cells inside np areas

                        both_conts=mat*(max(max(npconts)))+npconts;
                        %                 figure;
                        % plot to check if they match
                        figure('units','normalized','outerposition',[0 0 0.75 0.75]);
                        set(gcf,'color','w');

                        subplot(2,2,1); imagesc(zscore(std_projection(40:end-40,40:end-40)));subplot(2,2,2);imagesc(zscore(STD_stack_slices(40:end-40,40:end-40,scenario1_list(plane_i))))
                        mat_t = mean(cat(3,zscore(STD_stack_slices(:,:,scenario1_list(plane_i)))*max(max(both_conts)),both_conts),3);
                        subplot(2,2,3);imagesc(mat_t(40:end-40,40:end-40));
                        subplot(2,2,4);
                        imagesc((cell_f-np_f)');
                        saveas(gcf,fullfile(fid_exp,num2str(plane_i)), 'png');
                        close all;

                        %                 imagesc(both_conts);

                        %

                        %% plot
                        %                 figure;

                        %rename and clear
                        data(plane_i).CONTS=CONTS;
                        data(plane_i).cellconts=cellconts;
                        data(plane_i).npconts=npconts;
                        data(plane_i).cell_f=cell_f;
                        data(plane_i).np_f=np_f;
                    end
                end
                toc;
                save(fullfile(fid_exp,'datanew'),'data');
            end
        end
    end
end

function Y = tif_read(name)

info = imfinfo(name);

numFramesStr = regexp(info.ImageDescription, 'images=(\d*)', 'tokens');
numFrames = str2double(numFramesStr{1}{1});

% Use low-level File I/O to read the file
fp = fopen(name , 'rb');
% The StripOffsets field provides the offset to the first strip. Based on
% the INFO for this file, each image consists of 1 strip.
%     fseek(fp, info.StripOffsets, 'bof');
% Assume that the image is 16-bit per pixel and is stored in big-endian format.
% Also assume that the images are stored one after the other.
% For instance, read the first 100 frames

Y = zeros(info(1).Height,info(1).Width,numFrames);
for cnt = 1:numFrames
    Y(:,:,cnt) = fread(fp, [info.Width info.Height], 'uint32', 0, 'ieee-be')';
end
fclose(fp);


info = imfinfo(name);
T = numel(info);

d1 = info(1).Height;
d2 = info(1).Width;

Y = zeros(d1,d2,numFrames);
for t = 1:numFrames
    temp = imread(name, t, 'Info',info);
    Y(:,:,t) = temp(1:end,1:end);
end

end

function Y = tiff_read(name,T)

info = imfinfo(name);
if T == 0
    T = numel(info);
end

d1 = info(1).Height;
d2 = info(1).Width;
Y = zeros(d1,d2,T);
for t = 1:T
    temp = imread(name, t, 'Info',info);
    Y(:,:,t) = temp(1:end,1:end);
end
end