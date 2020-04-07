function Name = SubfolderMaker(FilePath,SubFolder)

%% Name of new folder 
Name = [FilePath,SubFolder];

%% Correct Name's folder delimiter based on  OS
if ismac
    Name(Name == '\') = '/'; 
elseif isunix
    Name(Name == '\') = '/';
elseif ispc
    Name(Name == '/') = '\'; 
else
    disp('Platform not supported')
end

%% Make new folder 
FolderMade = 0;
if ~(exist(Name,'file') == 7)
    mkdir(Name);
    FolderMade = 1;  
else
    FolderMade = 1;
end
end

