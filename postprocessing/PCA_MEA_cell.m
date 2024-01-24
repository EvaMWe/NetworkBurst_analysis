%% PCA-ANALYSIS
% by Eva-Maria Weiss
% Version from 23.07.2020
% feature table: double array, columns contain features, rows contain
% samples; 

%


function [featureVector,PC, eigVal] = PCA_MEA_cell(featureTable)

%% STEP 1: get data
%variableNames = resultCell{1,2};

      
%% STEP 2:subtract the mean from each of the data dimensions
    
    meanVec = mean(featureTable,1);
    stdVec = std(featureTable,0,1);
    dataAdjust = (featureTable-meanVec)./stdVec;
    
% % %  %% STEP 3: Calculate Covariance matrix
% % %     C = cov(dataAdjust);
% % %     
% % % %% STEP 4: Calculate Eigenvalues and Eigenvectors
% % %     [rightEigVec, eigVal] = eig(C);
% % %     eigValVec = eigVal(eigVal ~= 0);
% % %     eigValVec_sort = sort(eigValVec,'descend'); 
% % %     idx_sort = arrayfun(@(a) find(eigValVec_sort == a),eigValVec);
% % %     featureVector = rightEigVec(:,idx_sort);
% *featureVector = coeff:right eigenvector, Each column  contains coefficients for one principal component. 
%          The columns are in the order of descending component variance, latent. 
% *PC = score: input x rotated to new basis of principle components;
% *eigVal = latent: eigenvalues of covariance matrix of x arragned in
% descending order
   [featureVector,PC, eigVal] = pca(dataAdjust); 
   

end
    
