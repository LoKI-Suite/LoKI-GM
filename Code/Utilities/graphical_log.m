% Script to plot the errors of convergence cycles from log file

% Create a simple graphical interface to prompt for the file path
[file, path] = uigetfile('*.txt', 'Select log file');
if isequal(file, 0)
    disp('No files were selected. Exit.');
    return;
end
file_path = fullfile(path, file);
% Read the contents of the selected file
content = read_file(file_path);
% Parse the contents of the file
parsed_content = parse_content(content);
pressure_iterations = parsed_content.pressure_iterations;
pressure_errors = parsed_content.pressure_errors;
neutrality_iterations = parsed_content.neutrality_iterations;
neutrality_errors = parsed_content.neutrality_errors;
global_iterations = parsed_content.global_iterations;
global_errors = parsed_content.global_errors;
reduced_fields = parsed_content.reduced_fields;

% Plot the errors
plot_errors(pressure_iterations, pressure_errors, neutrality_iterations, neutrality_errors, global_iterations, global_errors);

function content = read_file(file_path)
    % Read all the lines in the file and return them as a single cell
    fid = fopen(file_path, 'r');
    content = textscan(fid, '%s', 'Delimiter', '\n');
    content = content{1};
    fclose(fid);
end

function parsed_content = parse_content(content)
    % Initialize lists to store the different types of errors and reduced fields
    pressure_iterations = [];
    pressure_errors = [];
    neutrality_iterations = [];
    neutrality_errors = [];
    global_iterations = [];
    global_errors = [];
    reduced_fields = [];
    % Go through every line of the file's contents
    for i = 1:length(content)
        line = content{i};
        % If the line contains "New pressure cycle iteration," extract the pressure error value
        if contains(line, 'New pressure cycle iteration')
            pressure_iterations(end+1) = str2double(extractBetween(line, '(', ')'));
            pressure_errors(end+1) = str2double(extractAfter(line, 'error = '));
        % If the line contains "New neutrality cycle iteration," extract the neutrality error value
        elseif contains(line, 'New neutrality cycle iteration')
            neutrality_iterations(end+1) = str2double(extractBetween(line, '(', ')'));
            neutrality_errors(end+1) = str2double(extractAfter(line, 'error = '));
        % If the line contains "New global cycle iteration," extract the global error value
        elseif contains(line, 'New global cycle iteration')
            global_iterations(end+1) = str2double(extractBetween(line, '(', ')'));
            global_errors(end+1) = str2double(extractAfter(line, 'error = '));
        % If the line contains "Updated reduced field," extract the value of the reduced electric field
        elseif contains(line, 'Updated reduced field')
            reduced_field_value = str2double(extractBetween(line, '(', ' '));
            reduced_fields(end+1) = reduced_field_value;
        end
    end
    % Return a structure containing the lists of errors and reduced electric fields
    parsed_content = struct('pressure_iterations', pressure_iterations, ...
                            'pressure_errors', pressure_errors, ...
                            'neutrality_iterations', neutrality_iterations, ...
                            'neutrality_errors', neutrality_errors, ...
                            'global_iterations', global_iterations, ...
                            'global_errors', global_errors, ...
                            'reduced_fields', reduced_fields);
end

function plot_errors(pressure_iterations, pressure_errors, neutrality_iterations, neutrality_errors, global_iterations, global_errors)
    figure;
    subplot(1, 3, 1);
    plot_error(global_iterations, global_errors, 'Global Errors', 'Global Cycle Errors');
    subplot(1, 3, 2);
    plot_error(neutrality_iterations, neutrality_errors, 'Neutrality Errors', 'Neutrality Cycle Errors');
    subplot(1, 3, 3);
    plot_error(pressure_iterations, pressure_errors, 'Pressure Errors', 'Pressure Cycle Errors');
end

function plot_error(iterations, errors, label, titleStr)
    errors = errors(:);
    pos_indices = find(errors >= 0);
    neg_indices = find(errors < 0);
    
    hold on;
    if ~isempty(pos_indices)
        plot(iterations(pos_indices), abs(errors(pos_indices)), 'o', 'DisplayName', [label ' (Positive)'], 'Color', 'blue');
    end
    if ~isempty(neg_indices)
        plot(iterations(neg_indices), abs(errors(neg_indices)), 'o', 'DisplayName', [label ' (Negative)'], 'Color', 'red');
    end
    plot(iterations, abs(errors), '-', 'Color', 'black', 'HandleVisibility', 'off');  % Line connecting all points
    hold off;
    
    xlabel('Iteration');
    ylabel('Error Value');
    title(titleStr);
    set(gca, 'YScale', 'log');
    legend;
    grid on;
    
    % Check if there is no data for a particular cycle and add a label "empty" if needed
    if isempty(errors)
        text(0.5, 0.5, 'empty', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'normalized');
    end
end