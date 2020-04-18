%% Setup for plots
rose_gold = [246/255,170/255,180/255];
sea_foam = [0 0.925 0.7];
burgundy = [0.502,0,0.125];
gold = [255/255 215/255 0];
make_figs = true;
save_figs = true;
formats = {'png','epsc','fig'};

%% messing around with filepaths
% put the path of your github repo here. End it with a slash
zoom_path = '../';

pathname_in = [zoom_path,'SplatGenData/'];
pathname_in = PathSlashCorrector(pathname_in);

pathname_out = SubfolderMaker(zoom_path,'Figures/');
pathname_out = PathSlashCorrector(pathname_out);

%% Data sets upon which to run code
nGroups = [2,5,10];
nCells = [1000,10000];
nGenes = [5000, 1000]
methods = {'_sI','_ALRA','_merged','_sI_thresh','_ALRA_hack'};

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
                h= figure(1);
                hold on
                for el = 1:i
                    gp_inds = (group_list == el);
                    plot(tSNE_vals(gp_inds,1),tSNE_vals(gp_inds,2),'o','MarkerSize',10)
                end
                set(gca,'FontName','Helvetica Neue','FontSize',20,'FontWeight','Bold')
                set(gca,'TickDir','out');
                title(char(['2-dim t-SNE with true color, method',methods{k}(2:end)]))
                
                % save plots in various figure formats
                if save_figs
                    for em = 1:length(formats)
                        save_name = [filepath_out,'tSNE_vis'];
                        saveas(h,save_name,char(formats(em)));
                    end
                    close all
                end
            end
            
        end
    end
end