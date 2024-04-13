%% clear previous work

clc; 
clear;
close all;

%% Read in data from the .mat file

load Employ_Data.mat

%% Example of a matrix from Employ_Data

% Define column names for each month
column_names = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', ...
    'Sep', 'Oct', 'Nov', 'Dec'};

% Define row names for each year
row_names = {'2000', '2001', '2002', '2003', '2004', '2005', '2006', ...
    '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015',...
    '2016', '2017', '2018', '2019', '2020', '2021', '2022', '2023', '2024'};

% Convert the matrix to a table with specified row and column names
employment_table = array2table(total_employment, 'VariableNames', ...
    column_names, 'RowNames', row_names);

% Display the table
fprintf('(Example) Total Employment Table: \n\n')
disp(employment_table);

%% Create proportion matrices

% for the Dec. to Jan. first column 
dec_1999_employed = 134523;
dec_1999_unemployed = 5653; 
dec_1999_notinlf = 68655;

% calculate the proportions using the calculate_proportion function (see
% documentation below)
p_employed_to_employed = calculate_proportions(total_employment, ...
    employed_to_employed, dec_1999_employed);
p_unemployed_to_employed = calculate_proportions(total_unemployment, ...
    unemployed_to_employed, dec_1999_unemployed);
p_notinlf_to_employed = calculate_proportions(not_in_lf, ...
    not_in_lf_to_employed, dec_1999_notinlf);
p_employed_to_unemployed = calculate_proportions(total_employment, ...
    employed_to_unemployed, dec_1999_employed);
p_unemployed_to_unemployed = calculate_proportions(total_unemployment, ...
    unemployed_to_unemployed, dec_1999_unemployed);
p_notinlf_to_unemployed = calculate_proportions(not_in_lf, ...
    not_in_lf_to_unemployed, dec_1999_notinlf);

p_employed_to_notinlf = (1 - p_employed_to_employed) - p_employed_to_unemployed;
p_unemployed_to_notinlf = (1- p_unemployed_to_employed) - p_unemployed_to_unemployed;


p_notinlf_to_notinlf = (1 - p_notinlf_to_employed) - p_notinlf_to_unemployed;

% set the NA values to 0 so that way we can calculate the sum and average
% later on
p_employed_to_notinlf(25, 2:12) = 0;
p_unemployed_to_notinlf(25, 2:12) = 0;
p_notinlf_to_notinlf(25, 2:12) = 0;

%% transition probabilities

% calculate the average proportions
% call the function calc_average (see the functions section below)

avg_p_employed_to_employed = calc_average(p_employed_to_employed(:));
avg_p_unemployed_to_employed = calc_average(p_unemployed_to_employed(:));
avg_p_notinlf_to_employed = calc_average(p_notinlf_to_employed(:));

avg_p_employed_to_unemployed = calc_average(p_employed_to_unemployed(:));
avg_p_unemployed_to_unemployed = calc_average(p_unemployed_to_unemployed(:));
avg_p_notinlf_to_unemployed = calc_average(p_notinlf_to_unemployed(:));

avg_p_employed_to_notinlf = calc_average(p_employed_to_notinlf(:));
avg_p_unemployed_to_notinlf = calc_average(p_unemployed_to_notinlf(:));
avg_p_notinlf_to_notinlf = calc_average(p_notinlf_to_notinlf(:));

% create the transition matrices using all the average proportions
% calculated above. 

fprintf("Transition Matrices\n")

long_run_smooth = [
    avg_p_employed_to_employed, avg_p_employed_to_unemployed, avg_p_employed_to_notinlf;
    avg_p_unemployed_to_employed, avg_p_unemployed_to_unemployed, avg_p_unemployed_to_notinlf;
    avg_p_notinlf_to_employed, avg_p_notinlf_to_unemployed, avg_p_notinlf_to_notinlf
]

deep_recession = [
    calc_deep_rec(p_employed_to_employed), calc_deep_rec(p_employed_to_unemployed), ...
    calc_deep_rec(p_employed_to_notinlf);
    calc_deep_rec(p_unemployed_to_employed), calc_deep_rec(p_unemployed_to_unemployed), ...
    calc_deep_rec(p_unemployed_to_notinlf);
    calc_deep_rec(p_notinlf_to_employed), calc_deep_rec(p_notinlf_to_unemployed), ...
    calc_deep_rec(p_notinlf_to_notinlf)
]

strong_recovery = [
    calc_strong_rec(p_employed_to_employed), calc_strong_rec(p_employed_to_unemployed), ...
    calc_strong_rec(p_employed_to_notinlf);
    calc_strong_rec(p_unemployed_to_employed), calc_strong_rec(p_unemployed_to_unemployed), ...
    calc_strong_rec(p_unemployed_to_notinlf);
    calc_strong_rec(p_notinlf_to_employed), calc_strong_rec(p_notinlf_to_unemployed), ...
    calc_strong_rec(p_notinlf_to_notinlf)
]

%% markov chain checks

% in order for markov chains to be calculated correctly, the sum of the
% rows must be 1.
if all(sum(long_run_smooth, 2), 1) && all(sum(deep_recession, 2), 1) && ...
    all(sum(strong_recovery, 2), 1)
    fprintf('All matrices row sum are 1. Transition matrix is valid.\n\n')
else
    fprintf('Error: row sum is not 1.\n')
    return;
end

% the matrix must be a square
if (size(long_run_smooth, 1) ~= size(long_run_smooth, 2))
   error('Transition matrix must be square.');
   return;
end

%% creating tables

% table for long run smooth
long_run_smooth_table = array2table(long_run_smooth, ...
    'VariableNames', {'To: Employed', 'To: Unemployed', 'To: Not in the Labor Force'}, ...
    'RowNames', {'From: Employed', 'From: Unemployed', 'From: Not in the Labor Force'});
disp('Long Run Smooth Transition Probabilities:');
disp(long_run_smooth_table);

% table for deep recession
deep_recession_table = array2table(deep_recession, ...
    'VariableNames', {'To: Employed', 'To: Unemployed', 'To: Not in the Labor Force'}, ...
    'RowNames', {'From: Employed', 'From: Unemployed', 'From: Not in the Labor Force'});
disp('Deep Recession Probabilities:');
disp(deep_recession_table);

% table for strong recovery
strong_recovery_table = array2table(strong_recovery, ...
    'VariableNames', {'To: Employed', 'To: Unemployed', 'To: Not in the Labor Force'}, ...
    'RowNames', {'From: Employed', 'From: Unemployed', 'From: Not in the Labor Force'});
disp('Strong Recovery Probabilities:');
disp(strong_recovery_table)

%% Markov chains for steps

% We assume that 1 step = 1 month
steps = 1000;

% call the markov function (see the last section) with the matrix, number
% of steps, and initial state values
final_prob_long_run_smooth = markov(long_run_smooth, steps, 0.14, 0.485, 0.375);
final_prob_deep_recession = markov(deep_recession, steps, 0.14, 0.485, 0.375);
final_prob_strong_rec = markov(strong_recovery, steps, 0.14, 0.485, 0.375);

% print the results
fprintf("\nLong Run Smooth Stable States:\n")
print_prob(final_prob_long_run_smooth, steps)
fprintf("\nDeep Recession Stable States:\n")
print_prob(final_prob_deep_recession, 1000)
fprintf("\nStrong Recovery Stable States:\n")
print_prob(final_prob_strong_rec, steps)

%% Transition Matrix Heat Maps (easier to compare probabilities from economic scenarios)

heatmap_graph(long_run_smooth, 'Long Run')
heatmap_graph(deep_recession, 'Deep Recession')
heatmap_graph(strong_recovery, 'Strong Recovery')

%% State transition diagram

stateLabels = {'Employment', 'Unemployment', 'Not in the Labor Force'};

% Plot the graph
transition_graph(digraph(long_run_smooth, stateLabels), 'Long Run')
transition_graph(digraph(deep_recession, stateLabels), 'Deep Recession')
transition_graph(digraph(strong_recovery, stateLabels), 'Strong Recovery')

%% Time Evolution Plot

time_evolution_plot(long_run_smooth, 1000, 0.14, 0.485, 0.375, 'Long Run')
time_evolution_plot(deep_recession, 1000, 0.14, 0.485, 0.375, 'Deep Recession')
time_evolution_plot(deep_recession, 1000, 0.14, 0.485, 0.375, 'Strong Recovery')

%% Function: Calculating the Proportion Matrix

function proportion_matrix = calculate_proportions(total_matrix, to_matrix,  dec_1999_total)
    % calculates the proportions month to month (ex: Jan. to Feb. etc)
    proportion_matrix = zeros(25, 12);
    for year = 1:25
        for month = 2:12
            proportion_matrix(year, month-1) = to_matrix(year, month) / total_matrix(year, month-1);
        end
    end
    % handling last column special case (Dec. to Jan.)
    proportion_matrix(1, 12) = to_matrix(1,1) / dec_1999_total;
    for year = 2:25 
        proportion_matrix(year, 12) = to_matrix(year, 1) / total_matrix(year-1, 12);
    end 
    % set the rest of the values in the matrix to 0
    proportion_matrix(25, 2:12) = 0;
end

%% Function: Calculating the averages

function average_value = calc_average(matrix)
    % calculates the average of the matrix (289 is hard coded because we
    % disregard the 0's)
    average_value = sum(matrix) / 289;
end

function deep_rec_value = calc_deep_rec(matrix)
    % calculates the deep recession values Dec - Jan 2008 to May - Jun 2009
    deep_rec_value = matrix(8, 12) + sum(matrix(9, :)) + sum(matrix(10, 1:5));
    deep_rec_value = deep_rec_value / 18;
end

function strong_rec_value = calc_strong_rec(matrix)
    % calculates the strong recession values May-Jun 2020 to Jan-Feb 2024
    strong_rec_value = matrix(25, 1) + sum(matrix(21, 5:12)) + sum(matrix(22:24, :), 'all');
    strong_rec_value  = strong_rec_value / 45;
end

%% Function: print

function print_prob = print_prob(matrix, steps)
    % print probabilities from the stable states
    fprintf("The probability to be Employed after %d steps is: %.4f\n", steps, matrix(1))
    fprintf("The probability to be Unemployed after %d steps is: %.4f\n", steps, matrix(2))
    fprintf("The probability to be Not in the Labor Force after %d steps is: %.4f\n", steps, matrix(3))
end

%% Markov Calculator

function final_state_vector = markov(matrix, steps, emp, unemp, nilf)
    % calculates the probability vector after a given amount of steps (assuming
    % 1 step = 1 month)
    
    % Determine the initial state vector (given by user input)
    initial_state = [emp; unemp; nilf];
        
    % Compute the final state probability after the given number of steps
    final_state_vector = initial_state' * matrix^steps;
end

%% Graph Functions

% creates a heatmap graph: indicates the probability from moving from one
% state to another (literally just easier to visualize)
function heatmap_graph = heatmap_graph(trans_matrix, type)
    stateLabels = {'Employment', 'Unemployment', 'Not in the Labor Force'};
    
    figure; % Creates a new figure
    heatmap(stateLabels, stateLabels, trans_matrix, ...
        'Title', strcat(type, ' Transition Probabilities'), ...
        'XLabel', 'To State', ...
        'YLabel', 'From State', ...
        'Colormap', cool, ... 
        'ColorbarVisible', 'on');
end

% state transition diagram: graphical representation of the markov chains
function transition_graph = transition_graph(digraph, type)
    figure;
    % edge thickness corresponds to transition probabilities
    weights = digraph.Edges.Weight * 10;
    plot(digraph, 'EdgeLabel',digraph.Edges.Weight, ...
        'LineWidth', weights, 'Layout','force');
    title(['State Transition Diagram for ', type]);
    xlabel('State');
    ylabel('Probability');
end

% time evolution plot: how the state distribution changes over time
function time_evolution_plot(trans_matrix, steps, emp, unemp, nilf, type)
    % Initialize the state vector as a row vector to match the 'markov' function
    initialState = [emp, unemp, nilf];
    
    % Create a matrix to store the evolution of the state distribution over time
    stateEvolution = zeros(steps + 1, 3);
    stateEvolution(1, :) = initialState;
    
    % Simulate the evolution of the state distribution over time
    for i = 2:steps + 1
        stateEvolution(i, :) = stateEvolution(i - 1, :) * trans_matrix;
    end
    
    % Plotting the evolution over time
    timeSteps = 0:steps;
    figure; % Creates a new figure window
    plot(timeSteps, stateEvolution(:, 1), 'b', 'LineWidth', 2); hold on;
    plot(timeSteps, stateEvolution(:, 2), 'r', 'LineWidth', 2);
    plot(timeSteps, stateEvolution(:, 3), 'g', 'LineWidth', 2);
    hold off;
    
    xlabel('Time Step');
    ylabel('Probability');
    title(['State Distribution Over Time for ', type]);
    legend({'Employment', 'Unemployment', 'Not in Labor Force'}, 'Location', 'best');
    grid on;
end


