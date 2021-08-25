%%This routine 
% 1. extracts the calculated bridge features 
% 2. runs novelty detection
% 3. runs POE analysis 
%  to detect and assess damage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% developed by Sifat Sharmeen Muin on 07/24/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clear console
close all
clear all
clc

%% read data
Dir_base = 'C:\Users\Sifat Muin\Dropbox\Project_CSMIP_2019';    %Parent folder
addpath(strcat(Dir_base,'\H-MC_codes'));
out_data=xlsread('ExtractedFeatures_7LS_Avg.xlsx');             %database of features from recorded data
out_data=out_data(3:15,:);                                      %MLO data from row 3 to last

%% extract features for bridge structure
X=out_data(:,[3,4,6,7,8,9]); %all features
X_D0=X(out_data(:,10)==0,:); %undamaged condition DS0
X_D1=X(out_data(:,10)==1,:); %Minor damage condition DS1
X_D2=X(out_data(:,10)==2,:); %Moderate damage condition DS2
X_D3=X(out_data(:,10)==3,:); %extensive damage condition DS3
X_D4=X(out_data(:,10)==4,:); %extensive damage condition DS4
X_D5=X(out_data(:,10)==5,:); %extensive damage condition DS5
%X_D6=X(out_data(:,9)==6,:); %extensive damage condition DS3
%length of each damage state
N=size(X_D0,1);     %DS0
ND1=size(X_D1,1);   %DS1
ND2=size(X_D2,1);   %DS2
ND3=size(X_D3,1);   %DS3
ND4=size(X_D4,1);   %DS3
ND5=size(X_D5,1);   %DS3

%% Divide into tratining set and test set
N_tr=round(0.5*N);  %Number of training data (50% training and 50% testing data from undamaged condition)
N_test_D0=N-N_tr;   %Number of testing data

%randomly select test data
p=1:1:N;
p_test = randperm(N,N_test_D0); %randomized test indices
p_tr=setdiff(p,p_test);         %training indices
X_training=X_D0(p_tr,:);        %selected training samples


X_test=X;               %test al recorded data
N_test=size(X_test,1);  %length of test data

%extract output vector for test set
    tag_test=out_data(:,10);    %binary output for damage detection

%%  Novelty Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Training

    %   k fold cross validation k=5 
    N_v=round(N_tr/5);          %diving the training set into k strips
for j=1:5
    %train
    p_v=(N_v*(j-1)+1):1:N_v*j;  %validation strip indcies
    p_tr_s=1:1:N_tr;            %full training set indices
    p_t=setdiff(p_tr_s,p_v);    %reduced training strip indices
    X_tr=X_training(p_t,:);     %reduced training set
    X_v=X_training(p_v,:);      %validation set
    N_tr_s=length(p_t);         %length of training set
    for k=1:50
        alpha(1,j)=0.5;                 %initial critical p value
        beta(j)=chi2inv(1-alpha(k,j),4);%threshold Mahalobis distance
        MD2=mahal(X_tr,X_tr);           %square of mahalobis distance of the reduced training set
        label_nov_tr=MD2>beta(j);       %detected novellty of the reduced training set
        z=zeros(N_tr_s,1);              % output of the reduced training set; all zeros
            L(k,j)=sqrt(mse(label_nov_tr,z));   %loss function
            
            % Compare loss function and update alpha
                c=-1;
                if k>1
                    if L(k,j)>L(k-1,j)
                    c=1;
                    end
                end
            alpha(k+1,j)=alpha(k,j)-2*1/N_tr_s*L(k,j)*c;    %updated critical p value

    end
    %validate
    MD2_v=mahal(X_v,X_tr);                          %square of mahalobis distance of the validation set
    label_nov_v=MD2_v>(chi2inv(1-alpha(end,j),4));  %detected novellty of the validation set
    z_v=zeros(N_v,1);                               %output of the validation set; all zeros.
    
    Acc_M1_v(j)=length(find(label_nov_v==z_v))/N_v*100;     %accuracy for the jth split

end

Val_acc=mean(Acc_M1_v)      %mean validation accuracy 
[~,index]=(max(Acc_M1_v))   %split number that provides highest accuracy

%% Finalize model parameter

 Th=chi2inv(mean(alpha(end,:)),4);%critical value for 4 df chi squared distribution for P=alpha
  
%% % test model %%%%%%%%%%%%%%%%%%%%%%%
    % predict
    
        MD2_D=mahal(X_test,X_training);                 %Mahalonobis distance
        label_nov_test=zeros(size(X_test,1),1);         %create a vector for predicted outputs
    
        %loop to check with threshold distance
        for m=1:N_test
            if MD2_D(m)> Th
                 label_nov_test(m)=1; 
            end
        end
        
   % evaluate accuracy    
   Acc_M1=length(find(label_nov_test==tag))/N_test*100  %accuracy for novelty detection
   cMat(tag,label_nov_test)                             %confusion matrix for novelty detection
%% tsne visualization

    Y=tsne(X_test);
figure
    gscatter(Y(:,1),Y(:,2),tag,'gr')
    xlabel('tSNE_1'); ylabel('tSNE_2');
    title("Actual Damage States")
    legend("Undamaged","Damaged")
figure
    gscatter(Y(:,1),Y(:,2),label_nov_test,'gr');
    xlabel('tSNE_1'); ylabel('tSNE_2');
    title("Predicted Damage States with Novelty Detection")
    legend("Undamaged","Damaged")
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POE envelope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% extract data from the SDOF analysisresults
    outfile2='POE_Features-500MLO_update.xlsx'; %File with SDOF features
    out_data2=xlsread(outfile2);
    tag_a=out_data2(:,8);                       %output damage information 
    X_a=out_data2(:,[1,2,4,5,6,7]);             %selected features
    N_a=size(X_a,1);                            %%length of analytical data

%% statistical information of the analytical data
X_ad=X_a((tag_a>0),:);%all damage data
mu_ad=mean(X_ad);
s_ad=cov(X_ad);
 
 X_ad0=X_a((tag_a==0),:);%undammage data
 mu_ad0=mean(X_ad0);
 s_ad0=cov(X_ad0);
 
 X_ad1=X_a((tag_a==1),:);%minor damage data
 mu_ad1=mean(X_ad1);
 s_ad1=cov(X_ad1);
 
 X_ad2=X_a((tag_a==2),:);%moderate damage data
 mu_ad2=mean(X_ad2);
 s_ad2=cov(X_ad2);
 
 X_ad3=X_a((tag_a==3),:);%extensive damage data
 mu_ad3=mean(X_ad3);
 s_ad3=cov(X_ad3);
 

%% train the POE model 

I=tag_a>0;              %indicator vactor for damaged events
%k fold validation k=5 
    p_a=1:1:N_a;        %all indices of training set
    Na_v=round(N_a/5);  %dividing the training set into 5 splits

    for j=1:5
        p_a_v=(Na_v*(j-1)+1):1:Na_v*j;      %validation strip indcies
        p_a_tr=setdiff(p_a,p_a_v);          %reduced  training set indices
        Xa_tr=X_a(p_a_tr,:);                % reduced training set for jth fold
        Xa_v=X_a(p_a_v,:);                  %validation set for jth fold
        Na_tr=length(p_a_tr);               %length of the training set
        for k=1:25
            a(1,j)=0.2;             %initial threshold probability for undamaged events
            b(1,j)=.5;              % initial scaling factor for the bandwidth vector
            h=std(X_a)*b(k);        %bandwidth vecor for kernel density
            I_tr=I(p_a_tr);         %indicator vector for the reduced training set
            
            %train
                for i=1:Na_tr
                    dist_tr=(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2);         %distance calc
                    if sum(exp(-dist_tr))==0                        %Probabilty of Exceedance calculation
                        POE_tr(i)=1;
                    else
                        POE_tr(i)=sum(I_tr.*(exp(-dist_tr)))/sum((exp(-dist_tr)));%POE estimate
                    end
                end
             %training set test   
                    label_POE_tr=zeros(size(Xa_tr,1),1);    %create predicted output vector
                    
                        for m=1:Na_tr
                            if POE_tr(m)> a(k,j)    %check with threshold probability
                                label_POE_tr(m)=1;  %predict output
                            end
                        end
                     %calculated loss and update threshold
                        L_POE(k,j)=sqrt(mse(label_POE_tr,I_tr));%loss function
                         c=-1;
                          if k>1
                            if L_POE(k,j)>L_POE(k-1,j)
                                     c=1;
                            end
                          end
                          
                            if abs(L_POE(k,j)) < 0.2
                                a(k+1,j)=a(k,j);                %updated threshold probability for undamaged events
                                b(k+1,j)=b(k,j);                %updated scaling factor for the bandwidth vector
                            else
                                a(k+1,j)=a(k,j)-0.2*L_POE(k,j)*c;   %updated threshold probability for undamaged events
                                b(k+1,j)=b(k,j)+0.2*L_POE(k,j)*c;   %updated scaling factor for the bandwidth vector
                            end
        end
    %validate
            I_v=I(p_a_v);%indicator vector for validation set
             %validation set test
            label_POE_v=zeros(Na_v,1);          %create predicted output vector
                 for m=1:Na_v
                     POE_v(m)=sum(I_tr.*(exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2)))); %Probabilty of Exceedance calculation
                        if POE_v(m)>a(end,j)        %check with threshold probability
                            label_POE_v(m)=1;       %predict output
                        end
                 end
        Acc_M2_v(j)=length(find(label_POE_v==I_v))/Na_v*100; %accuracy for the jth split
        %%error calculation 
         L_tr(j)=L_POE(k,j);                %trainng loss
         L_v(j)=sqrt(mse(label_POE_v,I_v)); %validation loss
    end
clear y

Val_acc=mean(Acc_M2_v)%average validationaccuracy
[~,index]=(max(Acc_M2_v))

%% finalize the model parameters
Th_POE=mean(a(end,:));%threshold probability value
h=std(X_a)*mean(b(end,:));%bandwitdth matrix

for i=1:N_a
    POE(i)=sum(I.*(exp(-((X_a-X_a(i,:)).^2/(2*h.^2))))/sum(exp(-(X_a-X_a(i,:)).^2/(2*h.^2)))); %final POE valuesof the SDOF set
end  

%% test damage detection
    %calculate POE of damage detection for test set
        for j=1:N_test
             dist_test=(X_a-X_test(j,:)).^2/(2*h.^2);
                if sum(exp(-dist_test))==0
                    POE_test(j)=1;
                else
                    POE_test(j)=sum(I.*(exp(-dist_test)))/sum((exp(-dist_test)));%POE values for the test set
                end
        end
    label_POE_test=zeros(size(X_test,1),1); %create output vector

    for m=1:N_test
        if label_nov_test(m)==1         % for the test samples that were detected as novelty
            if POE_test(m)> Th_POE      %check with the threashold POE
                 label_POE_test(m)=1; 
            else
                label_POE_test(m)=0;
            end
        else
                label_POE_test(m)=0;
        end
    end
%% Results
    Acc_M2=length(find(label_POE_test==tag))/N_test*100 %accuracy achieved by the platform; printed on console
    cMat(tag,label_POE_test)                            %confusion matrix
   %tSNE visualization 
   figure
        gscatter(Y(:,1),Y(:,2),label_POE_test,'gr');
        xlabel('tSNE_1'); ylabel('tSNE_2');
        title("Predicted Damage States with H-MC")
        legend("Undamaged","Damaged")
%% train damage assessment model
I_d1=tag_a==1;
I_d2=tag_a==2;
I_d3=tag_a==3;
I_d4=tag_a==4;
I_d5=tag_a==5;


      for n=1:5
        p_a_v=(Na_v*(n-1)+1):1:Na_v*n;%validation strip indcies
        p_a_tr=setdiff(p_a,p_a_v);%reduced  training set indices
        Xa_tr=X_a(p_a_tr,:);% reduced training set for jth fold
        Xa_v=X_a(p_a_v,:);%validation set for jth fold
        Na_tr=length(p_a_tr);%length of the training set
        for k=1:25
            a_DA(1,n)=.3;
            b_DA(1,n)=.1;% initial scaling factor for the bandwidth vector
            h=std(X_a)*b_DA(k,n);%bandwidth vecor for kernel density
            I_tr_1=I_d1(p_a_tr);%indicator vector for DS1
            I_tr_2=I_d2(p_a_tr);%indicator vector for DS2
            I_tr_3=I_d3(p_a_tr);%indicator vector for DS3
            I_tr_4=I_d4(p_a_tr);%indicator vector for DS3
            I_tr_5=I_d5(p_a_tr);%indicator vector for DS3

            tag_tr=tag_a(p_a_tr);
            %train
                for i=1:Na_tr
                    POE_1(i)=sum(I_tr_1.*(exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))));%POE estimate
                    POE_2(i)=sum(I_tr_2.*(exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))));%POE estimate
                    POE_3(i)=sum(I_tr_3.*(exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))));%POE estimate
                    POE_4(i)=sum(I_tr_4.*(exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))));%POE estimate
                    POE_5(i)=sum(I_tr_5.*(exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_tr(i,:)).^2/(2*h.^2))));%POE estimate
                     
                
                end
             %training set test   
            label_POE_DA_tr=zeros(size(Xa_tr,1),1);
            
                for m=1:Na_tr
                    if (POE_1(m)+POE_2(m)+POE_3(m)+POE_4(m)+POE_5(m))>a_DA(k,n)
                    [~,i3]=max( [POE_1(m),POE_2(m),POE_3(m),POE_4(m),POE_5(m)]);
                    label_POE_DA_tr(m)=i3;
                        
                    else
                        label_POE_DA_tr(m)=0;
                    end
                end
                
            L_POE_DA(k,n)=sqrt(mse(label_POE_DA_tr,tag_tr));%loss function
            c=1;
                if k>1
                    if L_POE_DA(k-1,n)>=L_POE_DA(k,n)
                        c=-1;
                    end
                end
                
                if L_POE_DA(k,n)< 0.2
                    a_DA(k+1,n)=a_DA(k,n);
                    b_DA(k+1,n)=b_DA(k,n);
                else
                    
                a_DA(k+1,n)=a_DA(k,n)-0.2*L_POE_DA(k,n)*c;%updated threshold probability for undamaged events
                b_DA(k+1,n)=b_DA(k,n)+0.005*L_POE_DA(k,n)*c;%updated scaling factor for the bandwidth vector
                end
                
                end
%validate
    tag_v=tag_a(p_a_v);
    %validation set test
    label_POE_DA_v=zeros(Na_v,1);
        for m=1:Na_v
            
            POE_v_1(m)=sum(I_tr_1.*(exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))));
            POE_v_2(m)=sum(I_tr_2.*(exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))));
            POE_v_3(m)=sum(I_tr_3.*(exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))));
            POE_v_4(m)=sum(I_tr_4.*(exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))));
            POE_v_5(m)=sum(I_tr_5.*(exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))))/sum((exp(-(Xa_tr-Xa_v(m,:)).^2/(2*h.^2))));
             
            if (POE_1(m)+POE_2(m)+POE_3(m)+POE_4(m)+POE_5(m))>a_DA(end,n)
                [~,i3]=max([POE_v_1(m),POE_v_2(m),POE_v_3(m),POE_v_4(m),POE_v_5(m)]);
                     label_POE_DA_v(m)=i3;              
                else
                    label_POE_DA_v(m)=0;    
                end
        end
      
  Acc_M2_DA_v(n)=length(find(label_POE_DA_v==tag_v))/Na_v*100; %accuracy for the jth split
  %error calculation 
  L_tr(n)=L_POE_DA(end,n);%trainng loss
  L_v(n)=sqrt(mse(label_POE_DA_v,tag_v));%validation loss
      end
      Val_acc=mean(Acc_M2_DA_v)
      [~,index]=(max(Acc_M2_DA_v))
    %% finalize the model parameters
h=std(X_a)*mean(b_DA(end,:));

   %% test damage assessment

 label_POE_DA_test=zeros(size(X_test,1),1);
 I_d1=tag_a==1;
 I_d2=tag_a==2;
 I_d3=tag_a==3;
  
 for j=1:N_test
POE_test_D1(j)=sum(I_d1.*(exp(-(X_a-X_test(j,:)).^2/(2*h.^2))))/sum((exp(-(X_a-X_test(j,:)).^2/(2*h.^2))));%POE value for DS1
POE_test_D2(j)=sum(I_d2.*(exp(-(X_a-X_test(j,:)).^2/(2*h.^2))))/sum((exp(-(X_a-X_test(j,:)).^2/(2*h.^2))));%POE value for DS2
POE_test_D3(j)=sum(I_d3.*(exp(-(X_a-X_test(j,:)).^2/(2*h.^2))))/sum((exp(-(X_a-X_test(j,:)).^2/(2*h.^2))));%POE value for DS3
POE_test_D4(j)=sum(I_d4.*(exp(-(X_a-X_test(j,:)).^2/(2*h.^2))))/sum((exp(-(X_a-X_test(j,:)).^2/(2*h.^2))));%POE value for DS3
POE_test_D5(j)=sum(I_d5.*(exp(-(X_a-X_test(j,:)).^2/(2*h.^2))))/sum((exp(-(X_a-X_test(j,:)).^2/(2*h.^2))));%POE value for DS3
  
 end
 % assign damage state to highest probability value
for m=1:N_test
    if label_POE_test(m)==1
    [~,i2]=max([POE_test_D1(m),POE_test_D2(m),POE_test_D3(m),POE_test_D4(m),POE_test_D5(m)]);
     label_POE_DA_test(m)=i2;
    else
        label_POE_DA_test(m)=0;
    end
 end
Acc_M2_D=length(find(label_POE_DA_test==tag_test))/N_test*100 %accuracy for damage assessment
cMat6(tag_test,label_POE_DA_test)%confusion matrix for damage assessment

%% write the actual and predicted results (optional)

%C=[X_test tag_test label_POE_DA_test];
%writematrix(C,strcat('Bridge_testset_DamageDetection',num2str(N_test),'.xlsx'));
