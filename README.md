# Employment Data Analysis and Markov Chain Transitions

## Overview

This MATLAB script is designed to read employment data, calculate transition probabilities between different states (employment, unemployment, and not in the labor force), and analyze these transitions using Markov chains. The script includes data visualization tools to better understand the changes over time and under different economic scenarios.

## Requirements

- MATLAB environment
- `Employ_Data.mat` file containing the employment data matrix.

## Script Components

### 1. Initialization

```matlab
clc; 
clear;
close all;
```
Clears the command window, workspace, and closes all figure windows.

### 2. Data Loading

```matlab
load Employ_Data.mat
```
Loads the employment data from the `.mat` file.

### 3. Data Preparation

Defines column names (months) and row names (years) and converts the employment data matrix into a table for easy viewing.

```matlab
column_names = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', ...
    'Sep', 'Oct', 'Nov', 'Dec'};
row_names = {'2000', '2001', '2002', '2003', '2004', '2005', '2006', ...
    '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015',...
    '2016', '2017', '2018', '2019', '2020', '2021', '2022', '2023', '2024'};

employment_table = array2table(total_employment, 'VariableNames', ...
    column_names, 'RowNames', row_names);
disp(employment_table);
```

### 4. Proportion Matrices Calculation

Calculates the proportions of transitions from one state to another.

```matlab
dec_1999_employed = 134523;
dec_1999_unemployed = 5653; 
dec_1999_notinlf = 68655;

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

p_employed_to_notinlf(25, 2:12) = 0;
p_unemployed_to_notinlf(25, 2:12) = 0;
p_notinlf_to_notinlf(25, 2:12) = 0;
```

### 5. Transition Probabilities

Calculates the average transition probabilities and constructs transition matrices for different economic scenarios (long-run smooth, deep recession, strong recovery).

```matlab
avg_p_employed_to_employed = calc_average(p_employed_to_employed(:));
avg_p_unemployed_to_employed = calc_average(p_unemployed_to_employed(:));
avg_p_notinlf_to_employed = calc_average(p_notinlf_to_employed(:));

avg_p_employed_to_unemployed = calc_average(p_employed_to_unemployed(:));
avg_p_unemployed_to_unemployed = calc_average(p_unemployed_to_unemployed(:));
avg_p_notinlf_to_unemployed = calc_average(p_notinlf_to_unemployed(:));

avg_p_employed_to_notinlf = calc_average(p_employed_to_notinlf(:));
avg_p_unemployed_to_notinlf = calc_average(p_unemployed_to_notinlf(:));
avg_p_notinlf_to_notinlf = calc_average(p_notinlf_to_notinlf(:));

long_run_smooth = [
    avg_p_employed_to_employed, avg_p_employed_to_unemployed, avg_p_employed_to_notinlf;
    avg_p_unemployed_to_employed, avg_p_unemployed_to_unemployed, avg_p_unemployed_to_notinlf;
    avg_p_notinlf_to_employed, avg_p_notinlf_to_unemployed, avg_p_notinlf_to_notinlf
];

deep_recession = [
    calc_deep_rec(p_employed_to_employed), calc_deep_rec(p_employed_to_unemployed), ...
    calc_deep_rec(p_employed_to_notinlf);
    calc_deep_rec(p_unemployed_to_employed), calc_deep_rec(p_unemployed_to_unemployed), ...
    calc_deep_rec(p_unemployed_to_notinlf);
    calc_deep_rec(p_notinlf_to_employed), calc_deep_rec(p_notinlf_to_unemployed), ...
    calc_deep_rec(p_notinlf_to_notinlf)
];

strong_recovery = [
    calc_strong_rec(p_employed_to_employed), calc_strong_rec(p_employed_to_unemployed), ...
    calc_strong_rec(p_employed_to_notinlf);
    calc_strong_rec(p_unemployed_to_employed), calc_strong_rec(p_unemployed_to_unemployed), ...
    calc_strong_rec(p_unemployed_to_notinlf);
    calc_strong_rec(p_notinlf_to_employed), calc_strong_rec(p_notinlf_to_unemployed), ...
    calc_strong_rec(p_notinlf_to_notinlf)
];
```

### 6. Markov Chain Checks

Ensures that the transition matrices are valid Markov chains by checking if they are square matrices and if the sum of each row is 1.

```matlab
if all(sum(long_run_smooth, 2), 1) && all(sum(deep_recession, 2), 1) && ...
    all(sum(strong_recovery, 2), 1)
    fprintf('All matrices row sum are 1. Transition matrix is valid.\n\n')
else
    fprintf('Error: row sum is not 1.\n')
    return;
end

if (size(long_run_smooth, 1) ~= size(long_run_smooth, 2))
   error('Transition matrix must be square.');
   return;
end
```

### 7. Transition Matrix Tables

Creates and displays tables for the transition matrices.

```matlab
long_run_smooth_table = array2table(long_run_smooth, ...
    'VariableNames', {'To: Employed', 'To: Unemployed', 'To: Not in the Labor Force'}, ...
    'RowNames', {'From: Employed', 'From: Unemployed', 'From: Not in the Labor Force'});
disp('Long Run Smooth Transition Probabilities:');
disp(long_run_smooth_table);

deep_recession_table = array2table(deep_recession, ...
    'VariableNames', {'To: Employed', 'To: Unemployed', 'To: Not in the Labor Force'}, ...
    'RowNames', {'From: Employed', 'From: Unemployed', 'From: Not in the Labor Force'});
disp('Deep Recession Probabilities:');
disp(deep_recession_table);

strong_recovery_table = array2table(strong_recovery, ...
    'VariableNames', {'To: Employed', 'To: Unemployed', 'To: Not in the Labor Force'}, ...
    'RowNames', {'From: Employed', 'From: Unemployed', 'From: Not in the Labor Force'});
disp('Strong Recovery Probabilities:');
disp(strong_recovery_table);
```

### 8. Markov Chains for Steps

Simulates the Markov chains over a specified number of steps (months) and prints the stable state probabilities.

```matlab
steps = 1000;

final_prob_long_run_smooth = markov(long_run_smooth, steps, 0.14, 0.485, 0.375);
final_prob_deep_recession = markov(deep_recession, steps, 0.14, 0.485, 0.375);
final_prob_strong_rec = markov(strong_recovery, steps, 0.14, 0.485, 0.375);

fprintf("\nLong Run Smooth Stable States:\n")
print_prob(final_prob_long_run_smooth, steps)
fprintf("\nDeep Recession Stable States:\n")
print_prob(final_prob_deep_recession, 1000)
fprintf("\nStrong Recovery Stable States:\n")
print_prob(final_prob_strong_rec, steps)
```

### 9. Visualization

Generates heatmaps and state transition diagrams to visualize the transition probabilities.

```matlab
heatmap_graph(long_run_smooth, 'Long Run')
heatmap_graph(deep_recession, 'Deep Recession')
heatmap

_graph(strong_recovery, 'Strong Recovery')

stateLabels = {'Employment', 'Unemployment', 'Not in the Labor Force'};
transition_graph(digraph(long_run_smooth, stateLabels), 'Long Run')
transition_graph(digraph(deep_recession, stateLabels), 'Deep Recession')
transition_graph(digraph(strong_recovery, stateLabels), 'Strong Recovery')

time_evolution_plot(long_run_smooth, 1000, 0.14, 0.485, 0.375, 'Long Run')
time_evolution_plot(deep_recession, 1000, 0.14, 0.485, 0.375, 'Deep Recession')
time_evolution_plot(strong_recovery, 1000, 0.14, 0.485, 0.375, 'Strong Recovery')
```

### 10. Functions

Defines various helper functions used in the script for calculating proportions, averages, printing probabilities, simulating Markov chains, and creating visualizations.

```matlab
function proportion_matrix = calculate_proportions(total_matrix, to_matrix,  dec_1999_total)
    ...
end

function average_value = calc_average(matrix)
    ...
end

function deep_rec_value = calc_deep_rec(matrix)
    ...
end

function strong_rec_value = calc_strong_rec(matrix)
    ...
end

function print_prob = print_prob(matrix, steps)
    ...
end

function final_state_vector = markov(matrix, steps, emp, unemp, nilf)
    ...
end

function heatmap_graph = heatmap_graph(trans_matrix, type)
    ...
end

function transition_graph = transition_graph(digraph, type)
    ...
end

function time_evolution_plot(trans_matrix, steps, emp, unemp, nilf, type)
    ...
end
```

## Usage

1. Place `Employ_Data.mat` in the same directory as the script.
2. Run the script in MATLAB.
3. The script will display tables, transition matrices, and visualizations in the MATLAB command window and figure windows.

## Contact

For any questions or issues, please contact the script author.
