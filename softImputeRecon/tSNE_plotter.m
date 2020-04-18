%% Setup for plots
rose_gold = [246/255,170/255,180/255];
sea_foam = [0 0.925 0.7];
burgundy = [0.502,0,0.125];
gold = [255/255 215/255 0];

color_matrix = [0 0.925 0.7;...
                        246/255,170/255,180/255;...
                        255/255 215/255 0;...
                        0.502,0,0.125;...
                        190/255,186/255,218/255;...
                        128/255,177/255,211/255;...
                        253/255,180/255,98/255;...
                        179/255,222/255,105/255;...
                        217/255,217/255,217/255;...
                        188/255,128/255,189/255]; 
                        
make_figs = true;
save_figs = false;
formats = {'png','epsc','fig'};

%% messing around with filepaths
% put the path of your github repo here. End it with a slash
zoom_path = '../';

pathname_in = [zoom_path,'SplatGenData/'];
pathname_in = PathSlashCorrector(pathname_in);

pathname_out = SubfolderMaker(zoom_path,'Figures/');
pathname_out = PathSlashCorrector(pathname_out);

%% Data sets upon which to run code
nGroups = 10; %[2,5,10];
nCells = 10000;% [1000,10000];
nGenes =  1000;%[5000, 1000]
methods = {'_sI','_ALRA','_merged','_sI_thresh','_ALRA_hack','_truth','_data'};

method_names = {'softImpute','ALRA','merged','softImpute thresh','ALRA zeros','Ground Truth','Observed'};
%% Load tSNE and operate on htem
for i = nGroups
    for j = 1:length(nCells)
        for k = 1:length(methods)
            dataset_name = char([ num2str(i), '_groups_', num2str(nCells(j)), '_cells_', num2str(nGenes(j)), '_genes/']);
            filepath_in = char([pathname_in,dataset_name]);
            filename_in = char([filepath_in,'tSNE',methods{k},'.csv']);
            tSNE_vals = csvread(filename_in,1,0);
            
            %ideas:
            % put a clustering method in here, color according to
            % cluster
            % load in the true clustering, color according to true cluster
            % (YES?this seems like better idea).
            % or do both and compare.
            
            filename_group = char([filepath_in,'group_data.csv']);
            group_list = csvread(filename_group);
            
            filepath_out = SubfolderMaker(pathname_out,dataset_name);
            if make_figs
                h= figure;
                hold on
                for el = 1:i
                    color_mixed = color_matrix(el,:);
                    gp_inds = (group_list == el);
                    plot(tSNE_vals(gp_inds,1),tSNE_vals(gp_inds,2),'.','MarkerSize',12,'Color',color_mixed)
                end
                set(gca,'FontName','Helvetica Neue','FontSize',20,'FontWeight','Bold')
                set(gca,'TickDir','out');
                title(char(['2 dim t-SNE of ',method_names{k},', colored by truth']))
                
                % save plots in various figure formats
                if save_figs
                    for em = 1:length(formats)
                        save_name = [filepath_out,'tSNE_vis',methods{k}];
                        saveas(h,save_name,char(formats(em)));
                    end
                else 
%                     pause
                end
%                 close all
            end
            
        end
    end
end