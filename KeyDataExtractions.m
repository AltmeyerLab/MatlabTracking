%% Figure out the number of cells that are there and assign first column
time_fr = input('Enter the number of frames: ');
c_time = input('Enter the column with the time: ');
c_x = input('Enter the column with the x-coordinate: ');
c_y = input('Enter the column with the y-coordinate: ');
TimeIntv = input('Enter the time-interval between frames in hours (e.g. 0.25 for 15min): ');
TimeTreat = input('Enter the frame number prior to addition of drug: ');
Col53BP1foc = input('Enter the column with the 53BP1 foci counts (16): ');
ColSumInt53BP1 = input('Enter the column with the SumInt of 53BP1 foci (18): ');
ColIntStart = input('Enter the number of the first of the 7 consecutive columns with post-staining data (23): ');
n = 1;
keytable = zeros(50,size(track,2));
for i = 1:size(track,1)
    if i == 1 || track(i,1) ~= track(i-1,1)
        keytable(n,1) = track(i,1);
        n = n + 1;
    end
end
%% Generate matrix for each cell track
TrackMat = zeros(time_fr,size(track,2),size(keytable,1));
for i = 1:size(keytable,1)
    for m = 1:time_fr
        for s = 1:size(track,1)
            if track(s,c_time) == m && track(s,1) == keytable(i,1)
                TrackMat(m,:,i) = track(s,:);
            end
        end
    end
end
fprintf('Track matrix generated\n');
%% Lineage assignment (Cols 2,3,4)
% Code for tagging mother cells with ID
for l = 1:size(keytable,1)
    keytable(l,2) = ceil((keytable(l,1))/8);
end
% Code for tagging first generation daughter cells with ID
for l = 1:size(keytable,1)
    if rem(keytable(l,1),2) == 1
        keytable(l,3) = 1;
    else
        keytable(l,3) = 2;
    end
end
% Code for tagging second generation daughter cells
for l = 1:size(keytable,1)
    if rem(keytable(l,1),8)<=4 && rem(keytable(l,1),8)~=0
        keytable(l,4) = 1;
    else
        keytable(l,4) = 2;
    end
end
fprintf('Lineages assigned in keytable columns 2-4\n');
%% Extraction of division times for mother cells (Col 5)
for i = 1:size(keytable,1)
    for m = 1:time_fr
        if rem(keytable(i,1),8) == 2 && TrackMat(m,c_time,i) ~= 0
            keytable(i,5) = TrackMat(m,c_time,i);
            keytable(i-1,5) = TrackMat(m,c_time,i);
            break
        end
    end
end
fprintf('Mother cell division times assigned keytable column 5\n');
%% Extraction of division times for daughter 1 cells (Col 6)
for i = 1:size(keytable,1)
    for m = 1:time_fr
        if (rem(keytable(i,1),8) == 3 || rem(keytable(i,1),8) == 4) && TrackMat(m,c_time,i) ~= 0
            keytable(i,6) = TrackMat(m,c_time,i);
            if keytable(i,1) == keytable(i-2,1)+2
                keytable(i-2,6) = TrackMat(m,c_time,i);
            elseif keytable(i,1) == keytable(i-1,1)+2
                keytable(i-1,6) = TrackMat(m,c_time,i);
            end
            break
        end
    end
end
%% Extraction of division times for daughter 2 cells (Col 7)
for i = 1:size(keytable,1)
    for m = 1:time_fr
        if i>4 && (rem(keytable(i,1),8) >= 5 || rem(keytable(i,1),8) == 0) && TrackMat(m,c_time,i) ~= 0
            keytable(i,7) = TrackMat(m,c_time,i);
            if keytable(i,1) == keytable(i-4,1)+4
                keytable(i-4,7) = TrackMat(m,c_time,i);
            elseif keytable(i,1) == keytable(i-3,1)+4
                keytable(i-3,7) = TrackMat(m,c_time,i);
            elseif keytable(i,1) == keytable(i-2,1)+4
                keytable(i-2,7) = TrackMat(m,c_time,i);
            elseif keytable(i,1) == keytable(i-1,1)+4
                keytable(i-1,7) = TrackMat(m,c_time,i);
            end
            break
        end
    end
end
%% Calculation of division times for first and second generations (Col 8+9)
for i = 1:size(keytable,1)
    if (rem(keytable(i,1),8) == 1 || rem(keytable(i,1),8) == 2) && keytable(i,6) ~= 0 
        keytable(i,8) = (keytable(i,6)-keytable(i,5))*TimeIntv;
    end
end
for i = 1:size(keytable,1)
    if rem(keytable(i,1),8) <= 4 && rem(keytable(i,1),8) ~= 0 && keytable(i,7) ~= 0 
        keytable(i,9) = (keytable(i,7)-keytable(i,6))*TimeIntv;
    end
end
%% Determination of time between treatment and division and foci and sumint53BP1 prior to NCS and foci and sumint53BP1 after NCS (Col 10 + 11 + 12 + 13 + 14)
for i = 1:size(keytable,1)
    for m = 1:time_fr
        if rem(keytable(i,1),8) == 1 || rem(keytable(i,1),8) == 2
            keytable(i,10) = TimeTreat - keytable(i,5);
            keytable(i,11) = max(TrackMat(TimeTreat-5:TimeTreat,Col53BP1foc,i));
            keytable(i,12) = max(TrackMat(TimeTreat-5:TimeTreat,ColSumInt53BP1,i));
            keytable(i,13) = max(TrackMat(TimeTreat+1:TimeTreat+8,Col53BP1foc,i));
            keytable(i,14) = max(TrackMat(TimeTreat+1:TimeTreat+8,ColSumInt53BP1,i));
        end
    end
end
%% Addition of post-staining information (Col 15-26)
PoststainingCol = input('Enter number of columns with post-staining data: ');
for i = 1:size(keytable,1)
    keytable(i,15:15+PoststainingCol-1) = TrackMat(time_fr,ColIntStart:ColIntStart+PoststainingCol-1,i);
end
%% Determination of lack of division (Col 27) or disappearance (Col 28) prior to post-staining images
for i = 1:size(keytable,1)
    if keytable(i,6) == 0
        keytable(i,15+PoststainingCol) = 1;
    elseif TrackMat(time_fr,1,i) == 0
        keytable(i,15+PoststainingCol+1) = 1;
    end
end
%% Complement columns in which data should be copied between daughter cells
for i = 1:size(keytable,1)
    for j = 1:size(keytable,1)
        if keytable(i,2) == keytable(j,2)
            keytable(j,[5 10]) = keytable(i,[5 10]);
        end
        if keytable(i,2) == keytable(j,2) && keytable(i,3) == keytable(j,3)
            keytable(j,[8 11 12 13 14]) = keytable(i,[8 11 12 13 14]);
        end
    end
end
%% Generation and saving of final table with data
VarNames(1:40) = {'Track ID','Mother ID','Daughter 1 ID','Daughter 2 ID','Mother Division frame'...
    ,'Daughter 1 Division frame','Daughter 2 Division frame','Daughter 1 division time','Daughter 2 division time','Time between NCS and mother division'...
    ,'53BP1 foci pre NCS','SumInt53BP1 pre NCS','Max foci post NCS','SumInt53BP1 post NCS','gH2AX foci PostStaining'...
    ,'SumIntgH2AX','Total Intensity DAPI','Mean Intensity DAPI','Total Intensity gH2AX','Mean Intensity gH2AX'...
    ,'Total Intensity p21','Mean Intensity p21','Total Intensity p53','Mean Intensity p53','Total Intensity pRb'...
    ,'Mean Intensity pRb','No_division','Disappear_by_end','Empty1','Empty2','Empty3','Empty4','Empty5','Empty6','Empty7','Empty8','Empty9','Empty10','Empty11','Empty12'};
KeyTable_final = array2table(keytable);
KeyTable_final.Properties.VariableNames = VarNames;
%%
writetable(KeyTable_final,'KeyTable_Untreated.xls');