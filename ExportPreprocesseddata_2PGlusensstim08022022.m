%% 2P glucose +sensory stim analysis draft
%%%%%% Analysis workflow:
% cell_f_extraction_2P_doublecheckplanes_24012022a.m
% ExportPreprocesseddata_2PGlusensstim08022022.m
% Combine_2Pdata0802.m
% Analysis_combined_data_2P.m
% Locate2Pcells.m
clear all
clc
fid = 'P:\Alexander\Summer 2019\To_do\Sensory_stim_exp';
animal_list = {'077','082'};
animal_list2 = {'126b','152f','G35','G36','m143a'};
Available_animal_list = getsubfolder_names_fun(fid);

idx_to_del = find(strcmpi(Available_animal_list,'Graphs')==1);
Available_animal_list(:,idx_to_del) = [];

% iterate animals
for animal_i = 1:length(Available_animal_list)
    exp_list = getsubfolder_names_fun(fullfile(fid,Available_animal_list{animal_i}));
    % iterate experiments for a given animal
    for exp_i = 1:length(exp_list)
        fid_exp = fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i});
        if exist(fullfile(fid_exp,'preprocessed_data.mat')) == 2 && 1==2
            disp(strcat('Exp done already: ',Available_animal_list{animal_i},'-',exp_list{exp_i}));
        else
            disp(strcat('Running exp:',fid_exp));
            
%         disp(exp_i);
        %determine if glucose or saline experiment
        if contains(exp_list{exp_i},'glu')
            exp_type = 'glu';
        elseif contains(exp_list{exp_i},'sal')
            exp_type = 'sal';
        end
        % load cell fluorescence data
        load(fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i},'datanew.mat'));

        %% Here combine cells that are the same between layers TODO!!
        clear preprocessed_data
        % Get m%atrices of contours 
        pixelsY=size(data(1).cellconts,1);
        pixelsX=size(data(1).cellconts,2);
        
        clear ROIs
        for plane_i = 1:length(data) % for each plane
            ROIs(plane_i).mat=zeros(pixelsY,pixelsX);
            plane_conts = data(plane_i).CONTS;
            for cell_i = 1:length(plane_conts)
                mat = zeros(pixelsY,pixelsX);
                matx = zeros(pixelsY,pixelsX);
                temp = plane_conts{cell_i};
                Xind = ones(pixelsY,1)*(1:pixelsX);
                Yind = (1:pixelsY)'*ones(1,pixelsX);
                [IN ON] = inpolygon(Xind,Yind,temp(:,1),temp(:,2));
                mat(find(IN))=1;%
      
                preprocessed_data(plane_i).ROIsmat{cell_i}=mat;
            end
        end
        % for each cell look for overlapping cells in a downstream planes
        cell_number = 0;
        for plane_i = 1:length(preprocessed_data) % for each plane except the last one
            preprocessed_data(plane_i).cell_ID = zeros(length(preprocessed_data(plane_i).ROIsmat),1);
        end
        for plane_i = 1:length(preprocessed_data)-1 % for each plane except the last one
            for cell_i = 1:length(preprocessed_data(plane_i).ROIsmat) % for each cell in that plane
                cell_to_take = preprocessed_data(plane_i).ROIsmat{cell_i};
 
                cell_similar_to = [];
                for cell_vs = 1:length(preprocessed_data(plane_i+1).ROIsmat) % next plane
                    cell_to_compare_to = preprocessed_data(plane_i+1).ROIsmat{cell_vs};
                    cell_to_take_vec = reshape(cell_to_take,1,[]);
                    cell_to_compare_to_vec = reshape(cell_to_compare_to,1,[]);
                    
                    idx1 = find(cell_to_take_vec == 1);
                    idx2 = find(cell_to_compare_to_vec == 1);
                    Overlap_percentages = [length(find(cell_to_compare_to_vec(idx1)==1))/length(idx1),length(find(cell_to_take_vec(idx2)==1))/length(idx2)];
                    cell_to_take_f = smoothdata(zscore(data(plane_i).cell_f(:,cell_i)),'movmedian',10);
                    cell_to_compare_to_f = smoothdata(zscore(data(plane_i+1).cell_f(:,cell_vs)),'movmedian',10);
                    recorded_framelen = min([length(cell_to_take_f),length(cell_to_compare_to_f)]);
                    cell_to_take_f = cell_to_take_f(1:recorded_framelen);cell_to_compare_to_f = cell_to_compare_to_f(1:recorded_framelen);
                    [C,lags] = xcorr(cell_to_take_f,cell_to_compare_to_f,'coeff');
                    [max_c,lag_at_c] = max(C);
                    lag_at_c = lags(lag_at_c);
                    if any(Overlap_percentages>0.05) && abs(lag_at_c)<10 && max_c >0.9
%                          figure;
%                          subplot(2,1,1);
%                          t =imfuse(repmat(cell_to_take,1,1,3),repmat(cell_to_compare_to,1,1,3));
%                          imagesc(t)
%                          subplot(2,1,2);
%                          plot(cell_to_take_f); hold on; plot(cell_to_compare_to_f);
%                          title(strjoin({'C =',num2str(max_c),'Lag =',num2str(lag_at_c),'Overlap =',num2str(max(Overlap_percentages))}));
%                         suptitle(strjoin({'Plane',num2str(plane_i),', cell,',num2str(cell_i),'&',num2str(cell_vs)}));
                        cell_similar_to(end+1) = cell_vs;
                    end
                end
                % Was the cell already found to be similar to smth else?
               
                if preprocessed_data(plane_i).cell_ID(cell_i) == 0
                    cell_number = cell_number+1; % add another cell
                    preprocessed_data(plane_i).cell_ID(cell_i) = cell_number;
                    if ~isempty(cell_similar_to)
                        
                        preprocessed_data(plane_i+1).cell_ID(cell_similar_to) = cell_number;
                    end
                else % Cell was already labelled as smth
                    if ~isempty(cell_similar_to)
                        preprocessed_data(plane_i+1).cell_ID(cell_similar_to) = preprocessed_data(plane_i).cell_ID(cell_i);
                    end
                end
            end
        end
        
        % fill the last plane
         for cell_i = 1:length(preprocessed_data(end).ROIsmat) % for each cell in that plane
             if preprocessed_data(end).cell_ID(cell_i) == 0
                 cell_number = cell_number+1; % add another cell
                 preprocessed_data(end).cell_ID(cell_i) = cell_number;
             end
         end
        %%
%         % concatenate a list of cells
%         cell_f_mat = [];
%         for data_i = 1:5 % for each slice
%             cell_f_mat = [cell_f_mat,data(data_i).cell_f(recorded_framelen)];
%         end
%         plot(zscore(cell_f_mat));
        
        % load stim parameters
        stim_file_name = getstimfile_names_fun(fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i}));
        load(fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i},stim_file_name));
        
        % load Frame info
if any(strcmpi(animal_list,Available_animal_list{animal_i}))
        filename = fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i},'FRAME.lvm');
        frame_info = readtable(filename,'FileType','text', 'HeaderLines', 22); % This works for new ID files
     
        frames = zeros(length(frame_info.FRAME),1);
        frames(find(frame_info.FRAME>2.5)) = 1;
        frame_trigs = find(diff(frames)==1)+1; 
        
        stim_trigs = zeros(length(frame_info.FRAME),1);
        fldnm_t = fieldnames(frame_info);
        fldnm_t = fldnm_t(3:end-4);
        for fldnm_i = 1:length(fldnm_t)
%             plot(frame_info.(fldnm_t{fldnm_i}));hold on;
            stim_trigs(find(frame_info.(fldnm_t{fldnm_i})>2.2)) = 1;
        end
         stim_trigs = find(diff(stim_trigs)==1)+1;
         
         stim_frames_posthoc = [];
         for stim_i = 1:length(stim_trigs)
             t =  (frame_trigs - stim_trigs(stim_i) > 0);
             idx = find(diff(t) == 1);
             stim_frames_posthoc(stim_i) = (idx );
         end

        % Check if this was running for too long ie another exp was started
        find_exp_end = min([find(diff(frame_trigs)>100);length(frame_trigs)]);
        
        
        fps = (frame_trigs(find_exp_end)-frame_trigs(1))/length(frame_trigs(1:find_exp_end));
        fps_per_plane = fps/6;
        Time_x = linspace(1,size(data(1).cell_f,1)/fps_per_plane,size(data(1).cell_f,1))/60;
        frames_per_plane = (length(frame_trigs(1:find_exp_end))-6) /6;

        % Get IP timing
        allhigh = zeros(length(frame_info.clitter),1);
        IP_timing = allhigh;
        allhigh(find(frame_info.clitter>2.0)) = 1;
        IP_timing(find(frame_info.clitter>1.0)) = 1;
        IP_timing = IP_timing-allhigh;
        IP_timing_sum = movsum(IP_timing,100);
        idx = find(IP_timing_sum == 100);
        % plot(IP_timing); hold on; plot(IP_timing_sum); xline(idx(1)-50);  xline(idx(end));
        IP_timing = [idx(1)-50,idx(end)];
      
        IP_timing_frame = [];
        for stim_i = 1:length(IP_timing)
            t =  (frame_trigs - IP_timing(stim_i) > 0);
            idx = find(diff(t) == 1);
            IP_timing_frame(stim_i) = (idx);
        end
        
       
        % Get running info
        filename = fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i},'loco.lvm');
        running_info = readtable(filename,'FileType','text', 'HeaderLines', 23); % This works for new ID files
        
        running_vec = zeros(length(running_info.Voltage_1),1);
        running_vec2 = running_vec;
        running_vec(find(running_info.Voltage_1>3)) = 1;
        
        running_vec2(find(running_info.Voltage_1<-1)) = 1;
        running_vec = running_vec+running_vec2;
        
        running_vec = (abs([0;diff(running_vec(frame_trigs(1):frame_trigs(find_exp_end)))]));
        
        running_X = linspace(1,length(running_vec)/(1000*60),length(running_vec));
        abs_run_vec = movsum(running_vec,[1000,0]);
        abs_run_vec = abs_run_vec * pi * (7*1*2) / 1440;
        IP_timing_minutes = running_X(IP_timing);
        abs_run_vec = downsample(abs_run_vec,100);
        running_X = downsample(running_X,100);

        %%
%         % Plot all cells
%         plot(Time_x,(smoothdata(zscore(cell_f_mat,[],1),'movmedian',round(fps))')); hold on
%         plot(Time_x,mean(zscore(cell_f_mat,[],1),2),'LineWidth',3);hold on; xline(IP_timing_frame(end)/(fps*60));
        
        stim_fldnms = fieldnames(parameters);
        idx_todel = find(contains(stim_fldnms,'stim')==0 | contains(stim_fldnms,'list')==1);
        stim_fldnms(idx_todel) = [];
        %%
        % Get a stim timing mat
        Stims_mat_duringrec = [];  Stims_mat_duringrec_v = []; stims_post = stim_fldnms(find(cellfun(@length,{stim_fldnms{:}}) == 6));
        for stim_i = 1:length(stims_post)
            Stims_mat_duringrec(stim_i,:) = parameters.(stims_post{stim_i});
            Stims_mat_duringrec_v = [Stims_mat_duringrec_v,Stims_mat_duringrec(stim_i,:)];
        end
%         stim_detection_offsets = (sort(Stims_mat_duringrec_v)-sort(stim_frames_posthoc));
        stim_frames_to_use = Stims_mat_duringrec_v;
else
    filename = fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i},'FRAME.lvm');
        frame_info = readtable(filename,'FileType','text', 'HeaderLines', 22); % This works for new ID files
     
        frames = zeros(length(frame_info.Voltage),1);
        stim_trigs = frames;
        frames(find(frame_info.Voltage>2.2)) = 1;
        frame_trigs = find(diff(frames)==1)+1; 
        stim_trigs(find(frame_info.Voltage_1>3)) = 1;
        stim_trigs = find(diff(stim_trigs)==1)+1; 
        
        for stim_i = 1:length(stim_trigs)
        end
        
        stim_frames_posthoc = [];
        for stim_i = 1:length(stim_trigs)
            t =  (frame_trigs - stim_trigs(stim_i) > 0);
            idx = find(diff(t) == 1);
            stim_frames_posthoc(stim_i) = (idx );
        end

        % Check if this was running for too long ie another exp was started
        find_exp_end = min([find(diff(frame_trigs)>100);length(frame_trigs)]);
        
        
        fps = (frame_trigs(find_exp_end)-frame_trigs(1))/length(frame_trigs(1:find_exp_end));
        fps_per_plane = fps/6;
        Time_x = linspace(1,size(data(1).cell_f,1)/fps_per_plane,size(data(1).cell_f,1))/60;
         frames_per_plane = (length(frame_trigs(1:find_exp_end))-6) /6;

        % Get IP timing
        IP_timing_frame = parameters.injection_frame;
        IP_timing_minutes = [(IP_timing_frame/ fps -15 )/60 ,(IP_timing_frame/ fps +15 )/60 ];
%         IP_timing_minutes = IP_timing_frame/ fps +10 )/60 ;
        
        % Get running info
        filename = fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i},'loco.lvm');
        running_info = readtable(filename,'FileType','text', 'HeaderLines', 23); % This works for new ID files
        
        running_vec = zeros(length(running_info.Voltage),1);
        running_vec(find(running_info.Voltage>-0.31)) = 1;

        
        running_vec = (abs([0;diff(running_vec(frame_trigs(1):frame_trigs(find_exp_end)))]));
        
        running_X = linspace(1,length(running_vec)/(1000*60),length(running_vec));
        abs_run_vec = movsum(running_vec,[1000,0]);
        abs_run_vec = abs_run_vec * pi * (7*1*2) / 20;
        abs_run_vec = downsample(abs_run_vec,100);
        running_X = downsample(running_X,100);
        
        
 
        %%
%         % Plot all cells
%         plot(Time_x,(smoothdata(zscore(cell_f_mat,[],1),'movmedian',round(fps))')); hold on
%         plot(Time_x,mean(zscore(cell_f_mat,[],1),2),'LineWidth',3);hold on; xline(IP_timing_frame(end)/(fps*60));
        
        stim_fldnms = fieldnames(parameters);
        idx_todel = find(contains(stim_fldnms,'stim')==0 | contains(stim_fldnms,'list')==1);
        stim_fldnms(idx_todel) = [];
        %%
        % Get a stim timing mat
         % Get a stim timing mat
        Stims_mat_duringrec = [];  Stims_mat_duringrec_v = []; stims_post = stim_fldnms(find(cellfun(@length,{stim_fldnms{:}}) == 6));
        for stim_i = 1:length(stim_fldnms)
            Stims_mat_duringrec(stim_i,:) = parameters.(stim_fldnms{stim_i});
            Stims_mat_duringrec_v = [Stims_mat_duringrec_v,Stims_mat_duringrec(stim_i,:)];
        end
%         stim_detection_offsets = (sort(Stims_mat_duringrec_v)-sort(stim_frames_posthoc));
        stim_frames_to_use = stim_frames_posthoc;
end
%%
% Plot to check

% plot(running_X,abs_run_vec); hold on;
% xline(IP_timing_minutes(end));
% plot(Time_x,zscore(mean(data(1).cell_f,2,'omitnan')));
% for stim_i = 1:length(stim_frames_posthoc)
%     xline(Time_x(round(stim_frames_posthoc(stim_i)/6)));
% end
%%

            % copy over relevant fields
            for data_i = 1:length(data)
                preprocessed_data(data_i).cell_f = data(data_i).cell_f;
                preprocessed_data(data_i).np_f = data(data_i).np_f;
            end

            preprocessed_data(1).fps = fps;
            preprocessed_data(1).Time_x = Time_x';
            preprocessed_data(1).Parameters = parameters;
            preprocessed_data(1).Stims_mat_duringrec = Stims_mat_duringrec;
            preprocessed_data(1).stim_frames_posthoc = stim_frames_posthoc;
            preprocessed_data(1).stim_frames_to_use = stim_frames_to_use;
%             preprocessed_data(1).stim_detection_offsets = stim_detection_offsets;
            preprocessed_data(1).running_X = running_X;
            preprocessed_data(1).abs_run_vec = abs_run_vec;
            preprocessed_data(1).running_vec = running_vec;
            preprocessed_data(1).IP_timing_minutes = IP_timing_minutes ;
            preprocessed_data(1).find_exp_end = running_X(end);
            preprocessed_data(1).frames_per_plane = frames_per_plane;
            preprocessed_data(1).frames_detected = length(frame_trigs(1:find_exp_end));
            preprocessed_data(1).frame_trigs = frame_trigs(1:find_exp_end);
            
%             preprocessed_data(1).exp_end = frame_trigs(1:find_exp_end);
            save(fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i},'preprocessed_data.mat'),'preprocessed_data', '-v7.3');
       clearvars -except exp_list fid Available_animal_list animal_list Available_animal_list exp_i  animal_i   %deletes all variables except X in workspace
        end
    end
end


function stim_file_name = getstimfile_names_fun(fid)
    mat_files = dir(fullfile(fid,'*.mat'));
    file_idx = find(contains({mat_files.name},'stim') == 1 & contains({mat_files.name},'bin') == 0);
    stim_file_name = mat_files(file_idx).name;
end
function subFolderNames = getsubfolder_names_fun(fid)
    % Get a list of all files and folders in this folder.
    files = dir(fid);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags); % A structure with extra info.
    % Get only the folder names into a cell array.
    subFolderNames = {subFolders(3:end).name};
end