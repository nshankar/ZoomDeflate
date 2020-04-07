function Name = PathSlashCorrector(Name)
% correct slashes in file paths 
if ismac
    Name(Name == '\') = '/'; 
elseif isunix
    Name(Name == '\') = '/';
elseif ispc
    Name(Name == '/') = '\'; 
else
    disp('Platform not supported')
end
end 