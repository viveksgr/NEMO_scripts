function var_set = variable_extract(dirs,matname,varname,mat_)
% Compiles variables named varname stored in file matname in directories
% given by cell array dirs into array var_set
% mat_(logical) = true implies var_set is a multidimensional array. Else a
% cell array. Only works if varname is a row vector of equal size.

if nargin<4
    mat_ = true;
end
ndir = length(dirs);
var_set = cell(1,ndir);
for ii = 1:ndir
    temp = load(fullfile(dirs{ii},matname),varname);
    eval(sprintf('var_set{ii}=temp.%s;',varname)) 
end

if mat_
    if size(var_set{1},1)>1
        var_set = horzcat(var_set{:});
    else
        var_set = vertcat(var_set{:});
    end
end
