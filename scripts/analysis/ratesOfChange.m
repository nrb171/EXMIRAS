% this is a script to process and analyze the results of an ensemble of EXMIRAS simulations.

% it is loosely broken up into three sections:
% 1.extracting the data from all of the results files. Preparing predictiors for predictor importance tests and neural network predictions.
% 2. plotting the results of the simulations
% 3. performing PCA on the results of the simulations to see if there are any correlations between the predictors and the results.
% 4. performing SVM and neural network predictions on the results of the simulations to see if there are any correlations between the predictors and the results.
% 5. plotting the results of the predictions.

% please note, this script has not been optimized, fully documented, or cleaned up. It is a work in progress and will be updated as needed.

% if you acquired this file from the APAR folder, please see https://github.com/nrb171/EXMIRAS for the latest version of this and other scripts.

files = dir('../../test/results/*.mat');
id = cell(1, length(files));
for ii = 1:length(files)
    file = files(ii);
    codes = strsplit(file.name, {'_','-','.'});
    id{ii} = strjoin(codes([3,5,6]), '-');
end
id = string(id);

vGrid = ["T", "qv","Zhh", "Zdr"];;
TGrid = 285:5:310;
ZhhGrid = 30:5:50;
hGrid = round(0.1:0.1:0.9,1);
qGrid = 0.1:0.1:0.9;
zdrGrid = 0.4:0.4:2;

shapeMap = NaN(numel(vGrid), numel(TGrid), numel(ZhhGrid), numel(hGrid), numel(zdrGrid),3);   
for vPlot = "T"
    disp(vPlot)

    fig1 = figure("Units","inches",Position=[0,0,3.25,3.25])
    hold on

    cmaps = {...
        interp1([1,9], [0.8,0,0.4; 0.8,0.8,0.4], [1:9]), ...
        interp1([1,9], [0,0.8,0.4; 0.8,0.8,0.4], [1:9]), ...
        interp1([1,9], [0.8,0.4,0; 0.8,0.4,0.8], [1:9]), ...
        interp1([1,9], [0.4,0.8,0; 0.4,0.8,0.8], [1:9]), ...
        interp1([1,9], [0,0.4,0.8; 0.8,0.4,0.8], [1:9]), ...
    };
    Results.hStart = [];
    Results.dBZstart = [];
    Results.tempStart = [];
    Results.zdrStart = [];
    Results.kurt = [];
    Results.skew = [];
    Results.strd = [];
    Results.coolRate = [];
    Results.totalCooling = [];
    for ii = 1:length(files)
        
        elementSplit = strsplit(id(ii),'-');
        ZhhId = str2num(elementSplit{3});
        elementIndex = find(ZhhId == [30:5:50]);
        color = cmaps{elementIndex}(1,:);
        cmaps{elementIndex} = circshift(cmaps{elementIndex}, -1);
        ExmirasRun = load(['../../test/', files(ii).name], 'ExmirasRun');
        ExmirasRun = ExmirasRun.ExmirasRun;
        sk = median(real(skewness(ExmirasRun.(vPlot))));
        kt = median(real(kurtosis(ExmirasRun.(vPlot)(~isnan(ExmirasRun.(vPlot))),1)));

        disp(ExmirasRun.ID)
        sd = median(std(ExmirasRun.(vPlot), 'omitmissing'), 'omitmissing'); 

        % save data into matrix ShapeMap
        shapeMap(vGrid == vPlot, TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4}),:) = [sd,sk, kt];

        Results.hStart(TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4})) = str2num(elementSplit{2})/100;
        Results.dBZstart(TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4})) = str2num(elementSplit{3});
        Results.tempStart(TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4})) = str2num(elementSplit{1});
        Results.kurt(TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4})) = kt;
        Results.skew(TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4})) = sk;
        Results.strd(TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4})) = sd;
        
        Results.totalCooling(TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4})) = median(min(ExmirasRun.T - ExmirasRun.T(1,:)));
        Results.coolRate(TGrid == str2num(elementSplit{1}), ZhhGrid == str2num(elementSplit{3}), hGrid == str2num(elementSplit{2})/100, zdrGrid == str2num(elementSplit{4})) = (median(min(ExmirasRun.T - ExmirasRun.T(1,:))))/max(ExmirasRun.Grids.tgrid, [], 'all')*3600;
        

        if elementSplit{1} == "290"
            disp('plotting')
            plot(linspace(min(ExmirasRun.Grids.tgrid(1:60:end))/60, max(ExmirasRun.Grids.tgrid(1:60:end))/60, numel(ExmirasRun.(vPlot)(1:120:end,1))), ...
                ExmirasRun.(vPlot)(1:120:end,1), ...
                'Color', color )
            xticks(0:15:180)
            xticklabels({'0', '15', '30', '45', '60', '75', '90', '105', '120', '135', '150', '165', '180'})
            xtickangle(45)
            %xlim([0, 3*3600])
            xlabel('Time [min]')
            ylabel(vPlot)
        end
    end

    if elementSplit{1} == "290"
        print2(fig1, sprintf('%s_290-K.png', vPlot), 'quality', '-r300')
        close all
    end

    % xq = TGrid;
    % yq = ZhhGrid;
    % zq = hGrid;
end

variablesToPlot = {'kurt', 'skew', 'strd', 'coolRate', 'totalCooling'};
fnames = fieldnames(Results);
for jj = 1:numel(fnames)
    Results.(fnames{jj})(Results.(fnames{jj}) == 0) = NaN
end
for jj = 1:length(variablesToPlot)

    for kk = 1:numel(zdrGrid)
        figure("Units","inches",Position=[0,0,3.25,3.25])
        yscale('linear')
        surf(Results.tempStart(:,:, 3,kk), ... 
            Results.dBZstart(:,:, 3,kk), ...
            Results.hStart(:,:,3,kk), ...
            Results.(variablesToPlot{jj})(:,:,3,kk), ...
            'EdgeColor', 'none') 
        view(15,20)
        xlabel('Temperature [K]')
        ylabel('Zhh [dBZ]')
        zlabel('relh [%]')
        zlim([0, 1])
        cb = colorbar;
        cb.Location = 'southoutside';
        colormap(jet(128))
        print2(gcf, sprintf('./%s_%d.png', variablesToPlot{jj}, kk), 'quality', '-r300')
        close all
    end
end

%% PCA
    vars = ["hStart", "dBZstart", "tempStart", "kurt", "skew", "strd" ];
    X = [];
    for v = vars
        X = [X,Results.(v)(:)];
    end
    X2 = [X, Results.coolRate(:), Results.totalCooling(:)];
    X2 = X2(~any(isnan(X2),2),:);
    nanMask = ~isnan(sum(X,2));
    X = X(nanMask,:);

    [coeff, score, latent, tsquared, explained, mu] = pca(X2);

    % X2 = [X, Results.coolRate(:), Results.totalCooling(:)];
    fig = figure("Units","inches",Position=[0,0,6.25,6.25])
    plotmatrix(X2)
    print2(fig, sprintf('./PCA.png'), 'quality', '-r300')

%% predictions
    % preprocessing
        vars = ["hStart", "dBZstart", "tempStart", "kurt", "strd" ];
        X = [];
        for v = vars
            X = [X,Results.(v)(:)];
        end 
        nanMask = ~isnan(sum(X,2));
        X = X(nanMask,:);

        shuffle = randperm(size(X,1));
        X = X(shuffle,:);
        YCoolingRate = Results.coolRate(nanMask);
        YCoolingRate = YCoolingRate(shuffle);
        YTotalCooling = Results.totalCooling(nanMask);
        YTotalCooling = YTotalCooling(shuffle);

        trainRatio = 0.8;

        trainingSubset = 1:round(trainRatio*size(X,1));
        testSubset = round(max(trainingSubset))+1:size(X,1);
        XVal = X(testSubset,:);
        YCoolingRateVal = YCoolingRate(testSubset);
        YTotalCoolingVal = YTotalCooling(testSubset);
        XTrain = X(trainingSubset,:);
        YCoolingRateTrain = YCoolingRate(trainingSubset);
        YTotalCoolingTrain = YTotalCooling(trainingSubset);

    % SVM
        nvSVM = {'Standardize',true,'KernelFunction','polynomial','KernelScale','auto'};
        mdlSVMCoolingRate = fitrsvm(XTrain, YCoolingRateTrain(:), nvSVM{:});

        % fprintf('Cooling Rate Convergence: %1.0f\n', mdlSVMCoolingRate.ConvergenceInfo.Converged);
        Error = YCoolingRateVal - mdlSVMCoolingRate.predict(XVal);
        fprintf('Cooling Rate Error: %1.2f\n', rms(Error));

        mdlSVMTotalCooling = fitrsvm(XTrain, YTotalCoolingTrain(:), nvSVM{:});
        % fprintf('Total Cooling Convergence: %1.0f\n', mdlSVMTotalCooling.ConvergenceInfo.Converged);
        Error = YTotalCoolingVal - mdlSVMTotalCooling.predict(XVal);
        fprintf('Total Cooling Error: %1.2f\n', rms(Error));
    
    % Neural net
        nvNN = {'Standardize',true, 'KFold',10};
        mdlNNCoolingRate = fitrnet(X, YCoolingRate(:));
        
        % fprintf('Cooling Rate Convergence: %1.0f\n', mdlNNCoolingRate.ConvergenceInfo.Converged);
        Error = YCoolingRateVal - mdlNNCoolingRate.predict(XVal);
        fprintf('Cooling Rate Error: %1.2f\n', rms(Error));

        means = mean(X,1);
        % colors = ["r", "g", "b", "c", "m"];
        

        mdlNNTotalCooling = fitrnet(XTrain, YTotalCoolingTrain(:));
        % fprintf('Total Cooling Convergence: %1.0f\n', mdlNNTotalCooling.ConvergenceInfo.Converged);
        Error = YTotalCoolingVal - mdlNNTotalCooling.predict(XVal);
        fprintf('Total Cooling Error: %1.2f\n', rms(Error));


    %% some plotting. This will need cleaned up later.
    % BLUE IS NN, RED IS SVM
        
        figure("Units","inches",Position=[0,0,3.25,3.25])
        % hold on
        error = [];
        for ii = 1:size(X,2)
            X2 = X;

            colsReplaceWithMean = (1:size(X,2) == ii)
            X2(:,colsReplaceWithMean) = repmat(means(colsReplaceWithMean), size(X2,1),1);

            predictionNN = mdlNNCoolingRate.predict(X2);
            predictionSVM = mdlSVMCoolingRate.predict(X2);
            errorNN(ii) = rms(predictionNN - YCoolingRate(:));
            errorSVM(ii) = rms(predictionSVM - YCoolingRate(:));

        end
        error = [errorNN; errorSVM];
        vars = ["hStart", "dBZstart", "tempStart", "kurt", "strd" ]
        ylabel('Predicted Cooling Rate Error')
        bar(vars,error)
        print2(gcf, './VI_CoolingRate.png', 'quality', '-r300')

        figure("Units","inches",Position=[0,0,3.25,3.25])
        % hold on
        error = [];
        for ii = 1:size(X,2)
            X2 = X;

            colsReplaceWithMean = (1:size(X,2) == ii)
            X2(:,colsReplaceWithMean) = repmat(means(colsReplaceWithMean), size(X2,1),1);

            predictionNN = mdlNNTotalCooling.predict(X2);
            predictionSVM = mdlSVMTotalCooling.predict(X2);
            errorNN(ii) = rms(predictionNN - YTotalCooling(:));
            errorSVM(ii) = rms(predictionSVM - YTotalCooling(:));

        end
        error = [errorNN; errorSVM];
        vars = ["hStart", "dBZstart", "tempStart", "kurt", "strd" ]
        ylabel('Predicted Total Cooling Error')
        bar(vars,error)
        print2(gcf, './VI_TotalCooling.png', 'quality', '-r300')

