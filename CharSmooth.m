function [Charcoal] = CharSmooth (Charcoal,Pretreatment,Smoothing,...
    Results)
% CharSmooth    Smooth charcoal record with 5 different methods.
%   [Charcoal] = CharSmooth (Charcoal,Pretreatment,Smoothing,Results)
% 
% Smooth interpolated charcoal series to estimate "background" charcoal 
% using 1 one of 5 methods:
% 1 -- Locally weighted scatter plot smooth using least squares linear 
%       polynomial fitting. (Matlab function in Curve Fitting Toolbox)
% 2 -- Lowess smoothing that is resistant to outliers. (Matlab function in
%       Curve Fitting Toolbox)
% 3 -- Moving average filter. (Matlab function in Curve Fitting Toolbox)
% 4 -- Running median. Each sample is assigned the median c.accI 
%       value within the smoothing window. This series is then smoothed
%       with a Lowess filter.
% 5 -- Running mode. Each sample is assigned the modal value from 
%       c.accI within the smoothing window. Within each window, 
%       charcoal.accI values are divided into 100 equally-spaced bins. This
%       series is then smoothed with a Lowess filter.

% Input variables:
% Charcoal -- output variable from CharPretreatment.m
% SmoothParams -- SmoothParams.method = 1, 2, 3, 4, or 5; this selectes the
%       method form the list above. SmoothParams.yr = years to smooth the
%       record over.
% 
% Output variables: this function adds one field to the input variable 
% Charocal:
% Charcoal.accIS: interpolated charcoal accumulation rate [# cm^-2 yr^-1],
%       smoothed according to input variables.

%% CREATE LOCAL VARIABLES
r = Pretreatment.yrInterp;
global plotData

%% SMOOTH CHAR SERIES TO ESTIMATE BACKGROUND CHAR
s = Smoothing.yr/r;  % number of data points over which to smooth the 
    % record. Value will be rounded when used in smooth function.
    
% Lowess
    charAccIS(:,1) = smooth(Charcoal.accI,s,'lowess'); % [# cm^-2 yr^-1] 
    
% Robust Lowess
    charAccIS(:,2) = smooth(Charcoal.accI,s,'rlowess'); % [# cm^-2 yr^-1] 

% Moving average
    charAccIS(:,3) = smooth(Charcoal.accI,s,'moving'); % [# cm^-2 yr^-1] 

% Running median
    for i = 1:length(Charcoal.accI)
        if i <= round(s/2)  % if 1/2 s twords start
            CHARi_t = Charcoal.accI(1:round(s)); % Charcoal.accI for year t
            charAccIS(i,4) = median(CHARi_t)';
        else
            if  i >= length(Charcoal.accI)-round(s) % if 1/2 s twords end
                CHARi_t = Charcoal.accI(length(Charcoal.accI)-...
                    round(s/2):end);
                charAccIS(i,4) = median(CHARi_t)';
            else    % else, you're in the middle of the record
                CHARi_t = Charcoal.accI(round(i-0.5*s):round(i+0.5*s));
                charAccIS(i,4) = median(CHARi_t)';
            end
        end
    end
        charAccIS(:,4) = smooth(charAccIS(:,4),s,'lowess'); 
            % [# cm^-2 yr^-1] smoothed with lowess filter 

% Running mode
    bin = 100;  % bins to divide Charcoal.accI into
    for i = 1:length(Charcoal.accI)
        if i <= round(s/2)  % if 1/2 s twords start
            CHARi_t = Charcoal.accI(1:round(s)); % Charcoal.accI for year t
            mode_bin = 0:range(CHARi_t)/bin:max(CHARi_t);
            [n,x] = hist(CHARi_t,bin);
            mode_in = x(find(n == max(n)));
            charAccIS(i,5) = median(mode_in);
        else
            if  i >= length(Charcoal.accI)-round(s)
                CHARi_t = Charcoal.accI(length(Charcoal.accI)-...
                    round(s/2):end);
                mode_bin = 0:range(CHARi_t)/bin:max(CHARi_t);
                [n,x] = hist(CHARi_t,bin);
                mode_in = x(find(n == max(n)));
                charAccIS(i,5) = median(mode_in);
            else
                CHARi_t = Charcoal.accI(round(i-0.5*s):round(i+0.5*s));
                mode_bin = 0:range(CHARi_t)/bin:max(CHARi_t);
                [n,x] = hist(Charcoal.accI(round(i-0.5*s):round(i+0.5*s)),...
                    bin);
                mode_in = x(find(n == max(n)));
                charAccIS(i,5) = median(mode_in);
            end
        end
    end
        charAccIS(:,5) = smooth(charAccIS(:,5),s,'lowess'); 
            % [# cm^-2 yr^-1] smoothed with lowess filter 
    % end running mode smoother

%%  PLOT RAW CHAR AND SMOOTHED SERIES, IF DESIRED
if plotData == 1 & Results.allFigures == 1
    figure (1);
    subplot(2,1,2)
    h = bar(Charcoal.ybpI,Charcoal.accI,1);
    set(h,'facecolor',[.5 .5 .5],'edgecolor',[.5 .5 .5])
    hold on
    plot(Charcoal.ybpI,charAccIS,'linewidth',1.5)
    ylim ([0, max(Charcoal.accI)]);
    %     ylim ([0, prctile(Charcoal.accI,99)]); 
    legend('C_i_n_t_e_r_p_o_l_a_t_e_d','Lowess','Robust Lowess',...
        'Moving Average','Moving Median','Moving Mode')
    xlabel ('time (cal. yr BP)')
    ylabel ('CHAR (# cm^-^2 yr^-^1)')
    set(gca,'xlim',[min(Charcoal.ybp)-100 max(Charcoal.ybp)],'box',...
        'off','tickdir','out','xdir','rev')
    title (['(b) C_i_n_t_e_r_p_o_l_a_t_e_d and different options for a ' num2str(Smoothing.yr) ' yr C_b_a_c_k_g_r_o_u_n_d']);
end

%% ADD charAccIS TO VARIABLE Charcoal%% CREATE OUTPUT VARIABLE Charcoal
Charcoal.accIS = charAccIS(:,Smoothing.method);