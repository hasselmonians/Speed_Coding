clear
close all
clc

% Parameters
min_R = 0.75; % minimal absolute R value a statistically significantly speed modulated cell has to display in order to be considered a speed cell
min_speed_cells_per_session = 4; % minimum number of speed cells defined by having min_R and statistically significantly speed modulated

%% create dialog box to choose experimental group
[animal_list,group_ID] = dialog_box_for_selecting_groups;
[textfile] = dialog_box_for_selecting_conditions;
%% select individual animals and sessions
% animal_list{1} = 'ArchTChAT_10';
% group_ID = 'ArchTChAT';
% textfile = '\temp_files.txt';

% counter for all cells
counter = 0; % initialize counter for all sessions containing

% go to data folder
animal_folder = 'D:\Dropbox (hasselmonians)\hdannenb\UnitRecordingData\';
cd(animal_folder)

for a = 1:length(animal_list)
    cd(strcat(animal_folder,animal_list{a}))
    if exist(strcat(pwd,textfile),'file') ==  2
        fileID = fopen(textfile);
        C = textscan(fileID,'%s');
        for session = 1:length(C{1})
            tic
            disp(['loading ',animal_list{a},', session ',C{1}{session}])
            load(C{1}{session})
            disp('session loaded')
            toc
            [~,session_name,~] = fileparts(root.name);
            %% for sessions with speed cells: pull out path to session file and speed cells within file
            tic
            disp('finding speed modulated cells in this session')
            cell_no = nan(length(speedMod.Non.sig_R),2); % initialize matrix for storage for cells matching the speed cells criteria within one session
            for i = 1:length(speedMod.Non.sig_R) % for loop over all cells active during the baseline session
                if speedMod.Non.sig_R(i) == 1 && abs(speedMod.Non.R(i)) > min_R % if cell showed statistically significant speed tuning
                    cell_no(i,:) = [root.cells(i,1) root.cells(i,2)];
                end
            end
            cell_counts = sum(~isnan(cell_no(:,1)));
            if cell_counts >= min_speed_cells_per_session
                counter = counter + 1;
                Pointer{counter,1}.session_name = [pwd,'\',session_name];
                Pointer{counter,1}.cel = cell_no;
            end
            toc
            clear cell_no cell_counts
        end
    end
end

% save result in animal parent folder
if exist('Pointer','var')
    cd(animal_folder)
    savename = ['Speed_modulated_cells','_',textfile(2:end-4),'_',group_ID];
    save(savename,'Pointer')
end

