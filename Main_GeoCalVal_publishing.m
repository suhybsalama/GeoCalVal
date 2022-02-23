%% GEOCAL/VAL model of Salama et al 2012
% |*This code was created by Suhyb Salama, Department of Water Resources, ITC 
% Faculty, University of Twente, The Netherlands*|
% 
% |Salama, M. S.; et al. (2012): Technical Notes: Calibration and validation 
% of geophysical observation models. Biogeosciences 9,6: 2195-2201, DOI: 10.5194/bg-9-2195-2012|
%% *Way of operation:*
% *Input*
% 
% * Supply an input a file with a pair of X and Y;
% * The first column is the independent variables X; e.g. concentration/soil 
% moisture. The second column is the dependent variables Y; e.g. absorption/ backscattering;
% * Example of a linear system is Y=slope*X+ intercept;
% * The idea is to find "slope" and "intercept" using the calibration set and 
% estimate their accuracies using the validation set.
% 
% *Output:*
% 
% * Probability distribution functions of module coefficients (resulting from 
% the calibration) 
% * Probability distribution functions of the error metrics, resulting from 
% the validation 
%% Note for the water case the intercept should be close to zero
% *Start the timing*
%%
tic
clear all
%%
data_file='in_situ_Data_index_of_retrieved_ap.m';
insitu_data=load(data_file);
          
%% 
% *Concentrations*

          Xdata=insitu_data(:,4);
        
%% 
% *Absorption, ap*

          measured_ap=insitu_data(:,1);
      iop=load('derived_iop_CE.m');
      
%% 
% *Derived ap*

Ydata=iop(:,1);
%% 
% *Clean up the data, this Ydata>0.45 is an outlier in the retrievals*

index=find(Xdata>0&Ydata>0&Ydata<0.5&insitu_data(:,1)>0);
Xdata=Xdata(index);
Ydata=Ydata(index);
measured_ap=measured_ap(index);
%% 
% **
% 
% *Assuming sigma_x*

sigma_x=0.07*Xdata;
 
sigma_ap2=(measured_ap-Ydata).^2;
    
%% 
% *Get the length of the vector x data*

n=length(Xdata);

%% 
% *Minimum number in the data points in a set *

min_set=7;
 
%% 
% *Determin the number of combinations*

kk=min_set:n-min_set;
nr_of_combination=length(kk);


%% 
% *This will be used as a multiplier for the number of iterations in the 
% second loop*

limit=1e3;


%% 
% _*We will work with the indices of the data*_

 set_index=[1:n];  
ij=0;
container_Predicted_y=[];
for i=1:nr_of_combination;
%% 
% the kth combination
%%
i;
 k=kk(i);
%% 
% Limit the combinations to the possible ones use when data points are ~15
%%
  %nr_combination=nchoosek(n,k);
  %if nr_combination< limit;  m=nr_combination; end
%% 
% *This to define the number of sets to be generated from each combination.*
% 
% This number m is proportional to the number of actual combination we can 
% get how many sets 
%%
  nr_combi=nchoosek(n,k);
  m=round(log10(nr_combi))*10;
 
%% 
% * *Generate random seeds for the random number generator;*
% * *This is done to guarantee that the generated sets of one combination are 
% different. *

    m_seed=randi([1,1e8],1,m);
    store_seeds(i)={m_seed};
    
    
%% 
% *Number of combinations in each calibration set with a fixed number of 
% data points   * 

     for j=1:m;   
      ij=ij+1;
       
   set_map(ij)=k;
 
%% 
% *  For each set we generate a random stream form a random seed*

 m_streams = RandStream.create('mrg32k3a','NumStreams',1,'Seed',m_seed(j));
 
%% 
% * Sample K unique values form the set_index *

    cal_index=randsample(m_streams,set_index,k);
   
     Xcal=Xdata(cal_index); 
     Ycal=Ydata(cal_index); 
     
 val_index = setdiff(set_index,  cal_index) ;
 
 Xval= Xdata(val_index);
  Yval= Ydata(val_index);
%% Call the retival algorithm
%%
  sigma_Xcal(ij)=std(Xcal);
    mean_Xcal(ij)=mean(Xcal);
      sigma_Ycal(ij)=std(Ycal);
    mean_Ycal(ij)=mean(Ycal);
    
      sigma_Xval(ij)=std(Xval);
    mean_Xval(ij)=mean(Xval);
      sigma_Yval(ij)=std(Yval);
    mean_Yval(ij)=mean(Yval);
    
  [cal_slope(ij),cal_intercept(ij),c_R,std_CALslope(ij),std_CALintercept(ij)]=lsqfity(Xcal,Ycal);
   cal_R(ij)=c_R.^2;
  Predicted_Y=Xval.*cal_slope(ij)+cal_intercept(ij);
  Predicted_X= (Yval-cal_intercept(ij))./cal_slope(ij);
%% Call the validation algorithm
%%
    [val_slope(ij),val_intercept(ij),v_R]=lsqfitma(Xval,Predicted_X);
 val_R(ij)=v_R.^2;
    MAE(ij)=mean(abs(Xval-Predicted_X));
      RMSE(ij)=sqrt(mean((Xval-Predicted_X).^2));  
        
  RMSE_y(ij)=sqrt(mean((Predicted_Y-Yval).^2));
   MAE_y(ij)=mean(abs(Yval-Predicted_Y));
       %----------------------------------------------------------------%
%% Clear the memo...
%%
     clear  Xcal Xval Ycal Yval cal_index val_index   m_streams
     end
    m_sigma_Xcal(i)=std(sigma_Xcal);
    c(i)=mean(mean_Xcal);
      m_sigma_Ycal(i)=std(sigma_Ycal);
    m_mean_Ycal(i)=mean(mean_Ycal);
    
      m_sigma_Xval(i)=std(sigma_Xval);
    m_mean_Xval(i)=mean(mean_Xval);
      m_sigma_Yval(i)=std(sigma_Yval);
    m_mean_Yval(i)=mean(mean_Yval);
    
   clear m_seed
end
%% *Create figure*
%%
subplot(1,3,1);
  createhistograms(cal_slope,1000,'Slope for Chla, [cm^{2}.mg^{-1}]',[min(cal_slope) max(cal_slope)]);
    subplot(1,3,2);
  createhistograms(cal_intercept,500,'Intercept for Chla, [m^{-1}]',[min(cal_intercept) max(cal_intercept)]);
 subplot(1,3,3);
  createhistograms(MAE,4000,'MAE of derived Chla, [mg.m^{-3}]',[min(MAE) max(MAE)]);
%% 
% 
% 
% *  End the time*
%%
  toc