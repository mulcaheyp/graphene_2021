clear all;clc
% load .mat files containing the features and traveling wave speeds
load('')
rec1 = formodelbuilding;
clear formodelbuilding
load('')
rec2 = formodelbuilding;
clear formodelbuilding
load('')
rec3 = formodelbuilding;
clear formodelbuilding

%% Train models
for z = 1:1000
    % Subselect from large training sets to make sure the training data is
    % enriched equally for all three records
    for1 = randi(size(rec1.ephys,1),size(rec2.ephys,1),1);
    for3 = randi(size(rec3.ephys,1),size(rec2.ephys,1),1);
    
    % get features and output concatenated
    featureinputE = [rec1.ephys(for1,:) ;rec2.ephys ;rec3.ephys(for3,:)];

    featureinputC = [rec1.caim(for1,:) ;rec2.caim ;rec3.caim(for3,:)];

    output = [rec1.output(for1,:) ;rec2.output ;rec3.output(for3,:)];

    features = [featureinputE featureinputC];
    
    % 5 fold cross validation
    crossval = 5; 
    nowithheld = floor(length(features)/crossval);
    forwitholding = [1:length(features)];

    for i = 1:crossval
        p = randperm(length(forwitholding),nowithheld);
        groups(i,:) = p;
        forwitholding(p) = [];
    end
    clear forwithholding
    
    % Train 5 bagged regression trees ensembles
    for i = 1:crossval
        % get training and testing set
        training = features;
        testing = training(groups(i,:),:);
        training(groups(i,:),:) = [];

        trainoutput = output;
        testingoutput = output(groups(i,:));
        trainoutput(groups(i,:)) = [];
        % train model
        mdl = TreeBagger(100,training,trainoutput,'Method','Regression','Surrogate','on',...
            'PredictorSelection','curvature','OOBPredictorImportance','on');
        % store importance of features
        imp(i,:) = mdl.OOBPermutedPredictorDeltaError;
        close all
        % predict and correlate with the testing data set
        preded = predict(mdl,testing);

        ccc = corr(preded,testingoutput);
        storecorrs(i) = ccc;
        storemodels{i} = mdl;    
    end

    clear training testing ccc preded mdl trainoutput testoutput
    clear groups
    
    % store the model that best predicts the testing set
    [~,idx] = max(storecorrs);
    gettheimps = imp(idx,:);
    clear imp
    bestmodel = storemodels{idx};
    storage_of_best_models{z} = bestmodel;
    storage_of_imps(z,:) = gettheimps;
    clear gettheimps
    clear bestmodel bestlevel idx storecorrs
end

%% Check performance across entire data set
% get features and the speeds
output = [rec1.output;rec2.output ;rec3.output];
featureinputE = [rec1.ephys ;rec2.ephys ;rec3.ephys];
featureinputC = [rec1.caim ;rec2.caim ;rec3.caim];
features = [featureinputE featureinputC];   

% predict the median speed and then correlate with experimentally observed
% distribution
for z = 1:1000
    prediction = predict(storage_of_best_models{z},features);
    rs = corr(prediction,output);
    rstore(z) = rs;   
end
[A,B] = sort(rstore);

% save the importances, correlations, and best performing models
forsaving.importances = storage_of_imps;
forsaving.corrs = rstore;
forsaving.bestperform = B(901:1000);
save('modelsinfo.mat','forsaving')
