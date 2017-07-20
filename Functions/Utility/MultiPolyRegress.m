function reg = MultiPolyRegress(Data,R,PW,varargin)
%   Multi-variable polynomial regression analysis. A by-product of ongoing computational
%   materials science research at MINED@Gatech.(http://mined.gatech.edu/)
%
%   reg = MultiPolyRegress(Data,R,PW) performs multi-variable polynomial 
%   regression analysis on row stacked dimensional data matrix Data. Data is 
%   an m-by-n (m>n) matrix where m is the number of data points and n is the number of 
%   independent variables. R is the m-by-1 response column vector and PW is the degree
%   of the polynomial fit. 
%   
%   reg = MultiPolyRegress(...,PV) restricts individual dimensions of
%   Data to particular powers PV in the polynomial expansion. PV is an
%   m-by-1 vector. A PV of [2 1] would limit a 2-dimensional 2nd degree polynomial to 
%	the terms that have x^2, x and y, eliminating the terms with y^2.
%
%   reg = MultiPolyRegress(...,'legend') turns on the legend calculation. Legend
%	is a symbolic vector that contains the corresponding polynomial term for each 
%	coefficient in the fitted polynomial. Legend functionality REQUIRES SYMBOLIC TOOLBOX.
%
%   reg = MultiPolyRegress(...,'figure') adds a scatter plot for the fit. 
%
%   reg = MultiPolyRegress(...,'range') adjusts the normalization of
%   goodness of fit measures: mean of absolute error (mae) and standard deviation
%   of absolute error. i.e. by default, mae is defined mean(abs(y-yhat)./y),
%   however, when this switch is used, the definition is changed to 
%   mean(abs(y-yhat)./range(y)). It is useful when your y vector (R in the 
%   syntax of this code) contains values that are equal to or very close to
%   0.
%
%   reg is a struct with the following fields:
%          FitParameters: Section Header 
%            PowerMatrix: A matrix that describes the powers for each term
%                         of the polynomial fit. It is useful for
%                         evaluating any future points with the calculated
%                         fit. Refer to the "Compose" section in the code 
%                         on how to use it. 
%                 Scores: Is a diagnostic reference. Displays the raw value
%                         of individual polynomial terms for each data
%                         point, before multiplication with coefficients.
%                         In other words, it is the matrix X you would have
%                         input in to the Statistical Toolbox function
%                         "regress".                    
%           Coefficients: For the calculated fit.
%                 Legend: (OPTIONAL) Identity of the corresponding polynomial 
%                         terms for Coefficients or 'No Legend'.
%                   yhat: Estimated values by the fit.
%              Residuals: y-yhat or R-yhat in syntax ofthis code,
%          GoodnessOfFit: Section Header
%                RSquare: 1-SSE/TSE
%                    MAE: Normalized Mean of Absolute Error
%                 MAESTD: Standard Deviation of Absolute Error
%     LOOCVGoodnessOfFit: '-----------------'
%              CVRSquare: RSquare of LOOCV
%                  CVMAE: MAE of LOOCV
%               CVMAESTD: MAESTD of LOOCV
%
%   Author : Ahmet Cecen

    % Align Data
    if size(Data,2)>size(Data,1)
        Data=Data';
    end
    
    % Arrange Input Arguments
    PV = repmat(PW,[1,size(Data,2)]);
    LegendSwitch='legendoff';
    FigureSwitch='figureoff';
    NormalizationSwitch='1-to-1 (Default)';
    if nargin > 3
        for i=1:length(varargin)
            if strcmp(varargin{i},'legend') == 1
                LegendSwitch='legendon';
            end
            if strcmp(varargin{i},'figure') == 1
                FigureSwitch='figureon';
            end
            if strcmp(varargin{i},'range') == 1
                NormalizationSwitch='Range';
            end
            if isnumeric(varargin{i}) == 1
                PV=varargin{i};
            end
        end
    end
    
    % Function Parameters
    NData = size(Data,1); % Number of points
    NVars = size(Data,2); % Number of variables
    RowMultiB = '1';
    RowMultiC = '1';
    Lim = max(PV);      % Highest degree of polynomial
    
    % Initialize
    A=zeros(Lim^NVars,NVars);

    % Create Colums Corresponding to Mathematical Base
    for i=1:NVars
        A(:,i)=mod(floor((1:Lim^NVars)/Lim^(i-1)),Lim);
    end

    % Flip - Reduce - Augment
    A=fliplr(A); A=A(sum(A,2)<=Lim,:); Ab=diag(repmat(Lim,[1,NVars])); A=[A;Ab];

    % Degree Conditionals
    for i=1:NVars
        A=A(A(:,i)<=PV(i),:);
    end

    % Build Framework
	if strcmp(LegendSwitch,'legendon')==1
		B=sym(zeros(size(A,1),NVars));
		for i=1:NVars
			B(:,i)=sym(['x',num2str(i)]);
			RowMultiB=strcat(RowMultiB,['.*B(:,',num2str(i),')']);
			RowMultiC=strcat(RowMultiC,['.*C(:,',num2str(i),')']);
		end
		% Create a Legend for Coefficient Correspondence
		B=B.^A; Legend = eval(RowMultiB); %#ok<NASGU>
	else
		for i=1:NVars
			RowMultiC=strcat(RowMultiC,['.*C(:,',num2str(i),')']);
		end
		Legend='No Legend';
	end

    % Allocate
    NLegend = size(A,1);
    Scores = zeros(NData,NLegend);
    
    % Compose
    for i=1:NData
        current=repmat(Data(i,:),[NLegend,1]);
        C=current.^A; %#ok<NASGU>
        Scores(i,:) = eval(RowMultiC);
    end

	% Use  QR to avoid explicit inversion and check rank. 
    [QQ,RR,perm] = qr(Scores,0);

    p = sum(abs(diag(RR)) > size(Scores,2)*eps(RR(1)));
  
    if p < size(Scores,2)
        warning('Rank Deficiency within Polynomial Terms!');
        WARNING = 1;
        RR = RR(1:p,1:p);
        QQ = QQ(:,1:p);
        perm = perm(1:p);
    else
        WARNING = 0;
    end
    
	% Ordinary Least Squares Regression
    b = zeros(size(Scores,2),1);
	b(perm) = RR \ (QQ'*R);
	yhat = Scores*b;                     
    r = R-yhat;   
	
    % Goodness of Fit
    r2 = 1 - (norm(r)).^2/norm(R-mean(R))^2;
    if strcmp(NormalizationSwitch,'Range')==1
        mae = mean(abs(r./range(R)));
        maestd = std(abs(r./range(R))); 
    else
        mae = mean(abs(r./R));
        maestd = std(abs(r./R));
    end
    
	% Leave One Out Cross Validation
	H=QQ*QQ';
    rCV=r./(1-diag(H));

    % LOOCV Goodness of Fit
    CVr2 = 1 - (norm(rCV)).^2/norm(R-mean(R))^2; 
    if strcmp(NormalizationSwitch,'Range')==1
        CVmae = mean(abs(rCV./range(R)));
        CVmaestd = std(abs(rCV./range(R))); 
    else
        CVmae = mean(abs(rCV./R));
        CVmaestd = std(abs(rCV./R));
    end
    
    % Construct Output
    reg = struct('FitParameters','-----------------','PowerMatrix',A,'X',Scores, ...
        'Coefficients', b, 'Degree', PW, 'Legend', Legend, 'yhat', yhat, 'Residuals', r, ...
        'GoodnessOfFit','-----------------', 'RSquare', r2, 'MAE', mae, 'MAESTD', maestd, ...
        'Normalization',NormalizationSwitch,'LOOCVGoodnessOfFit','-----------------', 'CVRSquare', ...
        CVr2, 'CVMAE', CVmae, 'CVMAESTD', CVmaestd,'CVNormalization',NormalizationSwitch,'Data',Data,'y',R,'Warning',WARNING);
    
    % Generate model output function
    mfun = @(xi) sum(reg.Coefficients.*prod(repmat(xi,size(reg.PowerMatrix,1),1).^reg.PowerMatrix,2));
    reg.modelFun = mfun;
    
    % Get Goodness-of-Fit Measures, Confidence Interval Functions
    reg = fit_measures(reg);
    reg = conf_int(reg);
    
    % Optional Figure
    if strcmp(FigureSwitch,'figureon')==1
        figure1 = figure;
        axes1 = axes('Parent',figure1,'FontSize',16,'FontName','Times New Roman');
        line(yhat,R,'Parent',axes1,'Tag','Data','MarkerFaceColor',[1 0 0],...
            'MarkerEdgeColor',[1 0 0],...
            'Marker','o',...
            'LineStyle','none',...
            'Color',[0 0 1]);
        xlabel('yhat','FontSize',20,'FontName','Times New Roman');
        ylabel('y','FontSize',20,'FontName','Times New Roman');     
        title('Goodness of Fit Scatter Plot','FontSize',20,'FontName','Times New Roman');
        line([min([yhat,R]),max([yhat,R])],[min([yhat,R]),max([yhat,R])],'Parent',axes1,'Tag','Reference Ends','LineWidth',3,'color','black');
        axis tight
        axis square
        grid on
    end
end

function reg = fit_measures(reg)

% mu_e is number of repetitions of each data point
% prior is the prior probability of the model

mu_e = 1;
prior = 1;

n = size(reg.Residuals,1); reg.n = n;
p = size(reg.Coefficients,1); reg.p = p;
dof = n-p; reg.dof = dof;

reg.MSE = sum(reg.Residuals.^2)/n;
reg.model_var = sum(reg.Residuals.^2)/(n-p);
reg.AIC = n*(log(2*pi()*reg.MSE) + 1) + 2*p;
reg.PB = prior * 2^(-p/2) * reg.MSE^(-mu_e/2);

end

function reg = conf_int(reg)

% alpha = degree of confidence for t dist... default 0.05

alpha = 0.05;

% Now Compute Z, the model fit at each of the X and Y points

f_x = @(xi) prod(repmat(xi,size(reg.PowerMatrix,1),1).^reg.PowerMatrix,2);

reg.meanCIfun = @(x) tinv(1-alpha/2,reg.dof)*sqrt((f_x(x)'*inv(reg.X'*reg.X)*f_x(x))*reg.model_var);
reg.predCIfun = @(x) tinv(1-alpha/2,reg.dof)*sqrt((1+f_x(x)'*inv(reg.X'*reg.X)*f_x(x))*reg.model_var);
reg.CIoptpt = fminsearch(reg.meanCIfun,mean(reg.Data,1));
reg.CIoptval = reg.predCIfun(reg.CIoptpt);

end
    
    