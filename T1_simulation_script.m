clc, clear, close all

% Define the folder path 
folder = 'folder path';

% Get a list of files with the extension .DSC and T1
fileList = dir(fullfile(folder, '*T1*.DSC'));


% Create an empty table to store the fitting parameters
rowNames = {'Temp/K'; 'Field/mT'; 'T1mono/ns'; 'std1'; 'Rsq mono'; 'T1str/ns'; 'stdstr'; 'Streach';'Std_Streach';'Rsq streach'; 'T1bi1/ns'; 'std bi1'; 'T1bi2/ns';  'std bi2'; 'Rsq bi' };
parameterTable_T1 = table('RowNames', rowNames);

% Loop through each file
for i = 1:numel(fileList)
    % Get the current file name
    filename = fullfile(folder, fileList(i).name);

% Defining exponential functions to fit data
    linear  = @(a,t) a(1) + a(2)*t;                                                                    % linear function for baseline substraction
    exp1 = @(p,t) p(1) + p(2)*exp(-t/p(3));                                                         % simple monoexponential decay
    exp2 = @(p,t) p(1) + p(2)*exp(-(t/p(3)).^p(4));                                                 % stretched monoexponential decay
    exp3 = @(p,t) p(1) + p(2)*exp(-t/p(3)) + p(4)*exp(-t/p(5));                                     % biexponential decay


    [t,spc, par, FileName] = eprload(filename);

% define parameters of the experiment   
    exp_title = FileName;
% Extract numerical values before 'K' using regular expressions
% Field = regexp(exp_title, '(?<=\D)(\d+)(?=G)', 'match'); %Field
% Field = str2double(exp_title);
    Field=par.B0VL*1e3;
% Convert the extracted numbers to numeric format
    Temperature = regexp(exp_title, '(?<=\D)(\d+)(?=K)', 'match'); %Temperature
    Temperature = str2double(Temperature);

    spc = real(spc); % takes the real part
 
    if spc(1) > spc(10)  % if due to the phase change we decay rather then reversion recovery
       spc=-1*spc;
    end 
    spc = (spc - min(spc))/(max(spc) - min(spc)); % scaling and normalizing the function 
   
    a=(1/2.3);

    indices = find(spc > a ); % we estimate T1 value 
    t0=t(indices(1)); % we estimate T1 value 

% depe
    if Temperature > 20
       t1=t0;
       t2=t0;
    else
       t1=t0;
       t2=t0/5;
    end


    if spc(1) > spc(10)  % if due to the phase change we decay rather then reversion recovery
       spc=-1*spc;
    end 

% Fit the data

    l_bounds1 = [-inf,-inf,0];  %lower limit for parameters 
    guess1 =    [1 -10 t0];  % starting fitting point
    u_bounds1 = [inf,0,inf];    %upper limit for parameters

    l_bounds2 = [-inf,-inf,0,0];         %lower limit for parameters 
    guess2 =    [1 -10 t0 1.5];       % starting fitting point
    u_bounds2 = [inf,0,inf,3];           %upper limit for parameters

    l_bounds3 = [-inf,-inf,0,-inf,0];        %lower limit for parameters 
    guess3 =    [0.5 -1000 t1 -10 t2];   % starting fitting point
    u_bounds3 = [inf,0,inf,0,inf];           %upper limit for parameters



    options = optimoptions('lsqcurvefit');
    options.Algorithm = 'trust-region-reflective';
    options.MaxFunctionEvaluations = 5.0e+03;
    options.MaxIterations = 3000; % Adjust the maximum number of iterations
    options.Display = 'iter'; % Set the display option to 'iter' for iterative output
    options.StepTolerance = 1e-3; % Adjust the step size tolerance


    [t1fit1,resnorm1,residual1,exitflag1,output1,lambda1,J1]= lsqcurvefit(exp1, guess1, t, spc, l_bounds1,u_bounds1,options);
    [t1fit2,resnorm2,residual2,exitflag2,output2,lambda2,J2]= lsqcurvefit(exp2, guess2, t, spc, l_bounds2,u_bounds2,options);
    [t1fit3,resnorm3,residual3,exitflag3,output3,lambda3,J3]= lsqcurvefit(exp3, guess3, t, spc, l_bounds3,u_bounds3,options);
% Generate plot
    figure;
    sz=40;
    scatter(t, spc, sz,'MarkerEdgeColor', 'black', 'DisplayName', 'Experimental');
    hold on;
    x_fit = t;
    y_exp1_fit = exp1(t1fit1,x_fit);
    y_exp2_fit = exp2(t1fit2,x_fit);
    y_exp3_fit = exp3(t1fit3,x_fit);
%  plot(x_fit, y_fit, 'r', 'LineWidth', 2);
    hold on;
    plot(x_fit, y_exp1_fit, 'g', 'LineWidth', 2,  'DisplayName', 'Monoexponential');
    hold on;
    plot(x_fit, y_exp2_fit, "Color", [1 1 0.3], 'LineWidth', 2,  'DisplayName', 'Stretched Exponential');
    hold on;
    plot(x_fit, y_exp3_fit, 'Color', [1 0 0], 'LineWidth', 2, 'DisplayName', 'Biexponential');
    hold on;
%  plot(t,spc,'m','LineWidth',2);
    legend('Location', 'best')
    xlabel('Time / ns');
    ylabel('Normalised echo intensity / arb. unit.');
    filetitle = sprintf('Field_%dG_Temperature_%dK_T1_VO_Hf_experiment', Field, Temperature);
    title(filetitle,'Interpreter','None');

    figFilename = sprintf('Field_%dG_Temperature_%dK_T1_VO_Hf_fit.fig', Field, Temperature);
    savefig(figFilename);

% Save the figure in .png format
    pngFilename = sprintf('Field_%dG_Temperature_%dK_T1_VO_Hf_fit.png', Field, Temperature);
    saveas(gcf, pngFilename, 'png'); 


% Display the extracted numbers
    disp('Extracted field:');
    disp(Field);
    disp('Extracted temperature:');
    disp(Temperature);

%% Statistical error estimators standart deviation and R square

    stdv1 = sqrt(diag(inv(J1.'*J1)*var(residual1)));
    stdv2 = sqrt(diag(inv(J2.'*J2)*var(residual2)));
    stdv3 = sqrt(diag(inv(J3.'*J3)*var(residual3)));
    
    
    % Residuals
    
    SStot = sum((spc-mean(spc)).^2);                 % Total Sum-Of-Squares
    SSres1 = sum((spc-exp1(t1fit1,t)).^2);      % Residual Sum-Of-Squares
    Rsq1 = 1-SSres1/SStot;                         % R^2
    
    SSres2 = sum((spc-exp2(t1fit2,t)).^2);
    Rsq2 = 1-SSres2/SStot;
    
    SSres3 = sum((spc-exp3(t1fit3,t)).^2);
    Rsq3 = 1-SSres3/SStot;
  

%% write parameters to the table 
% Define the data
   %name = {'Temp'; 'Field'; 'T1mono'; 'std1'; 'Rsq mono'; 'T1str'; 'stdstr';'Streach';'Std_Streach';'Rsq streach'; 'T1bi1'; 'std bi1'; 'T1bi2';  'std bi2'; 'Rsq bi' };
   parameters = [Temperature; Field; t1fit1(3); stdv1(3); Rsq1; t1fit2(3); stdv2(3);t1fit2(4); stdv2(4); Rsq2; t1fit3(3);stdv3(3); t1fit3(5); stdv3(5); Rsq3];% % Create a table
   % 
% % append table
   parameterTable_T1=addvars(parameterTable_T1, parameters);
%%
% export data into txt file format 
% Construct the filename
   filename = sprintf('Field_%dG_Temperature_%dK_T1_VO_Hf_experiment.txt', Field, Temperature);
   filenamemono = sprintf('Field_%dG_Temperature_%dK_T1_VO_Hf_monoexp.txt', Field, Temperature);
   filenameStretchedExp = sprintf('Field_%dG_Temperature_%dK_T1_VO_Hf_StretchedExp.txt', Field, Temperature);
   filenamebiexp = sprintf('Field_%dG_Temperature_%dK_T1_VO_Hf_BiExp.txt', Field, Temperature);

% Export data to the file
   dlmwrite(filename, [t, spc], 'delimiter', '\t', 'precision', 6);
   dlmwrite(filenamemono, [x_fit, y_exp1_fit], 'delimiter', '\t', 'precision', 6);
   dlmwrite(filenameStretchedExp, [x_fit, y_exp2_fit], 'delimiter', '\t', 'precision', 6);
   dlmwrite(filenamebiexp, [x_fit, y_exp3_fit], 'delimiter', '\t', 'precision', 6);


end
disp(parameterTable_T1);
trasposeParam=rows2vars(parameterTable_T1);
writetable(trasposeParam,'T1_VO_Hf.xlsx')