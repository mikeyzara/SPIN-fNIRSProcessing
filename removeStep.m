function y = removeStep(intensity)
%% removeStep removes step-like noise as per Zhou et al. (2020):
% Step-like noise can be caused by sudden loss of contact between optodes and the skin, 
% or interposition of hair, during data collection. If not removed, the step-like noise 
% would skew fNIRS signals. To remove step-like artifacts in the data of each channel (y), 
% the derivative of y was first estimated as X. Any absolute values in X, i.e., abs(Xi), 
% which were two standard deviations or more above the mean of the absolute values of X, 
% were set as zeros. Response y (with step-like artifacts removed) was then recovered by 
% calculating the cumulative sum of the updated X.

% 1) Calculate the derivative of intensity (y) and store in X
X = diff(intensity,1,1);
% X = diff(intensity(:,1));


% 2) Identify values in abs(X) that are two or more standard deviations above abs(X)
A = abs(X) >= mean(abs(X),1) + 2*std(abs(X),0,1); % Logical vector identifying which points are 2 standard deviations or more above abs(X)
% A = abs(X) >= mean(abs(X)) + 2*std(abs(X)); 
B = ~A; % Flip the values of A. This logical vector will be used to set the identified points to 0

% 3) Set identified points to 0 
updateX = X.*B; %Dot product between X and B will set identified points to 0

% 4) Recover the response without step-like artifacts
y = cumsum([intensity(1,:); updateX],1);
end