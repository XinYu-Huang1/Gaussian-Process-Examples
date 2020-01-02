% Script to generate random numbers drawn from a Rayleigh distribution.
% See Wikipedia article for definition of the Rayleigh
% probability density function and cumulative distribution function:
%     http://en.wikipedia.org/wiki/Rayleigh_distribution
% Uses the well known method of inverting the CDF of the desired probability
% function to get a function that can generate random numbers drawn
% from that desired distribution function,
% given as input numbers that have been drawn from a uniform distribution.
% Reference for that method:
% http://en.wikipedia.org/wiki/Inverse_transform_sampling
%

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
fontSize = 14;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
	cd(fileparts(which(mfilename)));
end

% Plot the Rayleigh cumulative distribution function
% over the range 0-20 with a sigma of 5.
x = 0:20;
sigma = 5;
rayleighCDF = 1 - exp(-x.^2 / (2*sigma^2));

subplot(3,1,1);
plot(x, rayleighCDF, 'LineWidth', 3);
caption = sprintf('Rayleigh CDF with sigma = %.2f', sigma);
title(caption, 'FontSize', fontSize);
set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full screen.
set(gcf,'name','Image Analysis Demo','numbertitle','off') 

% Ask user for a number of random numbers to generate.
userPrompt = 'Enter an integer number of random numbers to generate';
caNumberOfRandoms = inputdlg(userPrompt, 'Enter an integer',1,{'5000'});
if isempty(caNumberOfRandoms)
	% User clicked cancel
	return;
end
numberOfRandoms = round(str2num(cell2mat(caNumberOfRandoms)));
% Check for a valid integer.
if isempty(numberOfRandoms)
    % They didn't enter a number.  
    % They entered a character, symbols, or something else not allowed.
    numberOfRandoms = 5000;
    message = sprintf('I said it had to be an integer.\nI will use %d and continue.', numberOfRandoms);
    uiwait(warndlg(message));
end
    
% Generate numberOfRandoms uniformly distributed random numbers.
uniformlyDistributedRandomNumbers = rand(numberOfRandoms, 1);
subplot(3,2,3);
bar(uniformlyDistributedRandomNumbers, 'BarWidth', 1);
xlim([0 numberOfRandoms]);
caption = sprintf('%d Uniformly Distributed Numbers', numberOfRandoms);
title(caption, 'FontSize', fontSize);
xlabel('Element Number');
ylabel('Value');

%-----------------------------------------------------------------
% KEY PART, RIGHT HERE!!!!
% Invert the CDF of the Rayleigh function to get a function that can
% generate random numbers drawn from a Rayleigh distribution,
% given numbers drawn from a uniform distribution.
rayleighDistNumbers = sqrt(-log(1-uniformlyDistributedRandomNumbers)*(2*sigma^2));
%-----------------------------------------------------------------

% Plot the Rayleigh distributed numbers.
subplot(3,2,4);
bar(rayleighDistNumbers, 'BarWidth', 1);
xlim([0 numberOfRandoms]);
caption = sprintf('%d Rayleigh Distributed Numbers', numberOfRandoms);
title(caption, 'FontSize', fontSize);
xlabel('Element Number');
ylabel('Value');

% Get histogram of uniformly distributed numbers.
[countsU, binsU] = hist(uniformlyDistributedRandomNumbers, 50);
% Plot the uniformly distributed numbers.
subplot(3,2,5);
bar(binsU, countsU, 'BarWidth', 1);
caption = sprintf('Histogram of %d Uniformly Distributed Numbers', numberOfRandoms);
title(caption, 'FontSize', fontSize);
xlabel('Value');
ylabel('Count');

% Get histogram of Rayleigh distributed numbers.
% Observe that it's distribution is not flat like it is
% for the uniformly distributed numbers.
% It will take on the Rayleigh distribution shape.
[countsR, binsR] = hist(rayleighDistNumbers, 50);
subplot(3,2,6);
bar(binsR, countsR, 'BarWidth', 1);
caption = sprintf('Histogram of %d Rayleigh Distributed Numbers', numberOfRandoms);
title(caption, 'FontSize', fontSize);
xlabel('Value');
ylabel('Count');

message = sprintf('Done with demo.\nNote the difference in the histogram shapes\n(the bottom plots).');
uiwait(msgbox(message));

if numberOfRandoms > 300
	promptMessage = sprintf('Note the lists of random numbers (the middle plots) look biased higher\nbecause there are so many bars\nand the bars are so squished together.\n\nHere is what it looks like when you display fewer numbers');
	button = questdlg(promptMessage, 'Continue', 'Continue', 'Cancel', 'Continue');
	if strcmp(button, 'Cancel')
		return;
	end
	numberOfRandoms = 300;
	subplot(3,2,3);
	bar(uniformlyDistributedRandomNumbers(1 : numberOfRandoms), 'BarWidth', 1);
	xlim([0 numberOfRandoms]);
	caption = sprintf('%d Uniformly Distributed Numbers', numberOfRandoms);
	title(caption, 'FontSize', fontSize);

	subplot(3,2,4);
	bar(rayleighDistNumbers(1 : numberOfRandoms), 'BarWidth', 1);
	xlim([0 numberOfRandoms]);
	caption = sprintf('%d Rayleigh Distributed Numbers', numberOfRandoms);
	title(caption, 'FontSize', fontSize);
end
