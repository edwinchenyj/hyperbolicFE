function dirname = sim_directory(arg_cell)
dirname = '';
for i = 1:length(arg_cell)
    dirname = strcat(dirname,num2str(arg_cell{i}),'__');
    if i == length(arg_cell)
        % remove the last character
        dirname = dirname(1:end-1);
    end
end