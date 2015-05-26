%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: Kyle Booth
%Project: APS1022 - Financial Engineering II
%Date Due: Monday May 25, 2015 @ Noon
%Filename: Main.c (main file)
%Purpose: Main project file. Performs all analysis and produces results.

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0. PORTFOLIO BUILD, PHASE 0 - INPUT ASSETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%0. RAW DATA READ
%0.1 Data selection - Reads in data
sample = csvread('data\BCE.csv'); %used to determine size of input data
data = zeros(size(sample,1), size(sample,2)); %initialization of data matrix

%0.2 Data Read - Using dataRead function
data = dataRead(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. CLOSE PRICE MATRIX CREATION
%1.1 Parameters
period = size(sample,1)-1;
priceIndex = 5; %5 - adjusted close price, 4 - not adjusted close
closePrice = zeros(size(sample,1), size(data,3)); %initialization of closePrice matrix

%1.2 Matrix Creation
for i=1:size(data,3);    
    closePrice(:, i) = data(:,priceIndex,i);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. RATE OF RETURN CREATION
%2.1 Parameters
rateOfReturn = zeros(size(closePrice,1)-1, size(closePrice,2));

%2.2 Matrix Creation
for j=1:size(closePrice,2)
    for i=1:size(closePrice,1)-1
        rateOfReturn(i,j) = (closePrice(i,j) - closePrice(i+1,j))/(closePrice(i+1,j)); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3. MEAN RETURN CREATION
%3.1 Parameters
meanReturn = zeros(1,size(rateOfReturn,2));
meanIndex = 1;

%3.2 Matrix Creation
for j=1:size(rateOfReturn,2)
    for i=1:size(rateOfReturn,1)
        meanIndex = meanIndex*(1+rateOfReturn(i,j));
    end
    meanReturn(1,j) = (meanIndex^(1/(size(rateOfReturn,1)))) - 1;
    meanIndex = 1;
end
    
%3.3 Target Return Set Creation
targetReturn = zeros(1, ceil((max(meanReturn) - min(meanReturn))/(0.5/100)) + 1);
for i=1:size(targetReturn,2)
    if i == 1
        targetReturn(1,i) = min(meanReturn);
    elseif i == size(targetReturn,2)
        targetReturn(1,i) = max(meanReturn);
    else
    targetReturn(1,i) = targetReturn(1,i-1) + (0.5/100);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4. COVARIANCE CREATION & EIGENVALUE CHECKER
%4.1 Parameters
covarianceMatrix = zeros(size(data,3), size(data,3));
stdDevMatrix = zeros(1,size(rateOfReturn,2));
covarianceIndex = 0;

%4.2 Matrix Creation
for i=1:size(covarianceMatrix,1)
    for j=1:size(covarianceMatrix,2)
        for k=1:size(rateOfReturn,1)
            covarianceIndex = covarianceIndex + ((rateOfReturn(k,i) - meanReturn(1,i))*(rateOfReturn(k,j) - meanReturn(1,j)));
        end
        covarianceMatrix(i,j) = covarianceIndex/size(rateOfReturn,1);
        if (i == j)
            stdDevMatrix(1,j) = sqrt(covarianceMatrix(i,j));
        end
        covarianceIndex = 0;
    end
end

%4.3 Covariance Matrix Eigenvalue Check
lamda = eig(covarianceMatrix);
eigCheck = 0;

for i=1:size(lamda,1)
   if lamda(i,1) < 0
       eigCheck = eigCheck + 1;
   end
end
if eigCheck > 0 
    fprintf('WARNING: NEGATIVE EIGENVALUES DETECTED.\n\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5. BLACK-LITTERMAN METHOD FOR ADJUSTED MARKET EQUILIBRIUM RETURNS
%5.1 Parameters
totalCapitalization = 0;
capMatrix = zeros(1,size(rateOfReturn,2));
capWeightMatrix = zeros(1,size(rateOfReturn,2));
portfolioMean = 0;
portfolioVariance = 0;
riskAversion = 0;
riskFree = 2/100/12;
piReturn = zeros(1,size(rateOfReturn,2));
piReturnAlt = zeros(1,size(rateOfReturn,2));
    
%5.2 CapWeight Matrix Creation
for i=1:size(rateOfReturn,2)
    capMatrix(1,i) = data(1, 6, i);
    totalCapitalization = totalCapitalization + data(1,6,i);
end
for i=1:size(rateOfReturn,2)
    capWeightMatrix(1,i) = capMatrix(1,i)/totalCapitalization;
end
    
%5.3 Lambda (Risk Aversion Coefficient) Calculation
for i=1:size(rateOfReturn,2)
    portfolioMean = portfolioMean + meanReturn(1,i)*capWeightMatrix(1,i);
end
portfolioVariance = mtimes(mtimes(capWeightMatrix, covarianceMatrix), transpose(capWeightMatrix));
riskAversion = (portfolioMean - riskFree)/portfolioVariance
    
%5.4 Equilibrium Market Returns
piReturn = mtimes(mtimes(riskAversion, covarianceMatrix), transpose(capWeightMatrix));

%5.5 Custom View Implementation
%5.5.1 Parameters
numViews = 4;
qMatrix = zeros(numViews, 1);
pMatrix = zeros(numViews, size(rateOfReturn,2));
omega = zeros(numViews, numViews);
tau = 0.025;
%View 1: Bell Canada (BCE) will outperform Catamaran Corp (CCT) by 2.5%.
%View 2: TD Bank (TD) will outperform BMO Bank (BMO) by 3%.
%View 3: Valeant Pharmaceuticals (VRX) will underperform Power Corp. (POW) by 1%.
%View 4: Manulife Financial (MFC) returns will be 1%.
qMatrix = [0.0021; 0.0025; 0.001; 0.001];
pMatrix = [1, 0, -1, 0, 0, 0, 0, 0, 0, 0; 
           0, -1, 0, 0, 0, 0, 0, 0, 1, 0;
           0,  0, 0, 0, 0, 1, 0, 0, 0, -1;
           0,  0, 0, 0, 1, 0, 0, 0, 0, 0;];
omega = mtimes(mtimes(tau, pMatrix), mtimes(covarianceMatrix,transpose(pMatrix)));
piReturnAlt1 = inv(inv(mtimes(tau, covarianceMatrix)) + mtimes(mtimes(transpose(pMatrix),inv(omega)), pMatrix)); 
piReturnAlt2 = mtimes(inv(mtimes(tau, covarianceMatrix)), piReturn) + mtimes(mtimes(transpose(pMatrix),inv(omega)), qMatrix);
piReturnAlt = mtimes(piReturnAlt1, piReturnAlt2);
piReturnAlt = transpose(piReturnAlt);
piReturn = transpose(piReturn);
    
%5.6 Black-Litterman Target Return Set Creation
piTargetReturn = zeros(1, ceil((max(piReturnAlt) - min(piReturnAlt))/(0.25/100)) + 1);
for i=1:size(piTargetReturn,2)
    if i == 1
        piTargetReturn(1,i) = min(piReturnAlt);
    elseif i == size(piTargetReturn,2)
        piTargetReturn(1,i) = max(piReturnAlt);
    else
    piTargetReturn(1,i) = piTargetReturn(1,i-1) + (0.25/100);
    end
end

piTargetReturn2 = zeros(1, ceil((max(piReturn) - min(piReturn))/(0.25/100)) + 1);
for i=1:size(piTargetReturn2,2)
    if i == 1
        piTargetReturn2(1,i) = min(piReturn);
    elseif i == size(piTargetReturn2,2)
        piTargetReturn2(1,i) = max(piReturn);
    else
    piTargetReturn2(1,i) = piTargetReturn2(1,i-1) + (0.25/100);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%6. PORTFOLIO BUILD - MVO QUADRATIC OPTIMIZATION - NOMINAL MODEL
%6.1 Parameters

Q = covarianceMatrix;
c = zeros(size(covarianceMatrix(1,:)));
a = []; 
Aeq = -[piReturn; ones(size(piReturn(1,:)))]; %see 6.2 algorithm for Beq values
b = [];
ub = inf(size(covarianceMatrix(:,1)));
%lb = -inf(size(covarianceMatrix(:,1))); %IF SHORT SELLING ALLOWED
lb = zeros(size(covarianceMatrix(:,1))); 
xWeights = zeros(size(covarianceMatrix,1), size(piTargetReturn2,2));
xRisk = zeros(1, size(targetReturn,2));
portfolioRisk = zeros(1, size(piTargetReturn2,2));
portfolioReturn = zeros(1,size(piTargetReturn2,2));
portfolioScalar = 0;

%6.2 Optimization with Quad Prog
for i=1:size(piTargetReturn2,2)
    Beq = -[piTargetReturn2(1,i); ones(size(piTargetReturn2(1,i)))];
    [xWeights(:,i), xRisk(1,i)] = quadprog(Q, c, a, b, Aeq, Beq, lb, ub);
end

%6.3 Calculation of Portfolio Returns and Risk for Various R
for i=1:size(piTargetReturn2,2)
    for j=1:size(covarianceMatrix,1)
        portfolioScalar = portfolioScalar + (piReturn(1,j)*xWeights(j,i));
    end
    portfolioReturn(1,i) = portfolioScalar;
    portfolioRisk(1,i) = sqrt(2*xRisk(1,i));
    portfolioScalar = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%7. PORTFOLIO BUILD - MVO QUADRATIC OPTIMIZATION - BLACK-LITTERMAN MODEL
%7.1 Parameters

Q = covarianceMatrix;
c = zeros(size(covarianceMatrix(1,:)));
a = []; 
Aeq = -[piReturnAlt; ones(size(piReturnAlt(1,:)))]; %see 6.2 algorithm for Beq values
b = [];
ub = inf(size(covarianceMatrix(:,1)));
%lb = -inf(size(covarianceMatrix(:,1))); %IF SHORT SELLING ALLOWED
lb = zeros(size(covarianceMatrix(:,1))); 
pixWeights = zeros(size(covarianceMatrix,1), size(piTargetReturn,2));
pixRisk = zeros(1, size(piTargetReturn,2));
piPortfolioRisk = zeros(1, size(piTargetReturn,2));
piPortfolioReturn = zeros(1,size(piTargetReturn,2));
piPortfolioScalar = 0;

%7.2 Optimization with Quad Prog
for i=1:size(piTargetReturn,2)
    Beq = -[piTargetReturn(1,i); ones(size(piTargetReturn(1,i)))];
    [pixWeights(:,i), pixRisk(1,i)] = quadprog(Q, c, a, b, Aeq, Beq, lb, ub);
end

%7.3 Calculation of Portfolio Returns and Risk for Various R
for i=1:size(piTargetReturn,2)
    for j=1:size(covarianceMatrix,1)
        piPortfolioScalar = piPortfolioScalar + (piReturnAlt(1,j)*pixWeights(j,i));
    end
    piPortfolioReturn(1,i) = piPortfolioScalar;
    piPortfolioRisk(1,i) = sqrt(2*pixRisk(1,i));
    piPortfolioScalar = 0;
end

%7.4 Plot Results for Efficient Frontier
figure(1)
hold all
plot(portfolioRisk, portfolioReturn);
plot(piPortfolioRisk, piPortfolioReturn);
clear title;
title('APS1022 Project: Black-Litterman Portfolio Optimization');
xlabel('Volatility (Sigma)');
ylabel('Portfolio Return (Monthly, R)');