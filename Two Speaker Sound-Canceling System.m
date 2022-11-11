clear
clearvars
clc
%h.Visible = 'off'



f = 1000; %frequency

T = 1/f; %period





t=0;

    %this is sound power at sound source
    speaker1pos = [0,0];
    speaker2pos = [1,0];



    

    x_0 = str2double(inputdlg('Enter x axis listener position','X position',[1 1]))
    y_0 = str2double(inputdlg('Enter y axis listener position','Y position',[1 1]))

    observerpos = [x_0,y_0];
    
    
    di1 = pdist([speaker1pos;observerpos],'euclidean'); %diagonal distance from speaker 1 to obs
    di2 = pdist([speaker2pos;observerpos],'euclidean'); %diagonal distance from speaker 2 to obs

    
    t1 = (di1/343); %time to arrival of speaker 1 in ms
    t2 = (di2/343); %time to arrival of speaker 2 in ms

    newphase = ((t2-t1)*2)/(2*pi());


    phi1 = (360/(1/T))*(t1); %phase shift of speaker 1
    phi2 = (360/(1/T))*(t2); %phase shift of speaker 2


    
    y1 = sin(2*pi*f*t-phi1);
    y2 = -sin(2*pi*f*t-phi2);
 


    L_w1 = 10*log10(y1/(10^-12)); %sound power level of speaker 1
    L_w2 = 10*log10(y2/(10^-12)); %sound power level of speaker 2



    L_p_y1 = L_w1 - 10*log10(di1) - 8; %sound pressure level of speaker 1
    L_p_y2 = L_w2 - 10*log10(di2) - 8; %sound pressure level of speaker 2



    L_p_y1max = max(L_p_y1);  %max SPL of speaker 1
    L_p_y2max = max(L_p_y2); %max SPL of speaker 2

 
    L_p_tot = 10*log10((10.^(L_p_y1/10))+(10.^(L_p_y2/10)));
  
    
    
    
    
    
    
    
    

    
    stepsize = 0.1;

    
    i=0;
    j=0;
    x_sample = 0:stepsize:100;
    y_sample = 0:stepsize:100;
    
    G_Lp_matrix = zeros(1,1);
    
    magfactor = 50;
    for i = 1:50
        for j = 1:50
            G_di1 = pdist([speaker1pos;[i/magfactor,j/magfactor]], 'euclidean');
            G_di2 = pdist([speaker2pos;[i/magfactor,j/magfactor]],'euclidean');
 


            G_Lp_y1 = L_w1 - 10*log10(G_di1) - 8;
            G_Lp_y2 = L_w2 - 10*log10(G_di2) - 8;


      
           G_Lp_matrix(i,j) = real(10*log10(10.^(G_Lp_y1/10)+10.^(G_Lp_y2/10)));
        end
    end

    surf(G_Lp_matrix)












