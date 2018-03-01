function output = parse_parameter(funcname,parsed,unparsed)
if isempty(unparsed)
    % end of recursion, remove the last comma
    parsed = parsed(1:end-1);
    disp([funcname,'(',parsed,')']);
    output = eval([funcname,'(',parsed,')']);
else
    parsed = strcat(parsed,'''',num2str(unparsed{1}),'''',',');
    
    num = length(unparsed{2});
    output = cell(1,num);
    for i_count = 1:length(unparsed{2})
        % parse the list for the parameter
        if iscell(unparsed{2})
            % if the options are strings
            temp_parsed = strcat(parsed,'''',num2str(unparsed{2}{i_count}),'''',',');
        else
            % if the options are numbers
            temp_parsed = strcat(parsed,num2str(unparsed{2}(i_count)),',');
        end
        % recursion
        output{i_count} = parse_parameter(funcname,temp_parsed,unparsed(3:end));
    end
    % unfolding cell level. ugly.
    if num == 1
        output = output{1};
    end
end


end