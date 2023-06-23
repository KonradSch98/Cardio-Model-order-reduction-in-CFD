%%% this script writes all the data used for calculation in one .bin file
%%% which allows muh faster loading process compared to reading every 
%%% FVM file for all the time steps

%time steps
steps = 250;
%vector size of one time step
vec_size = 5136636;
cd 'D:\Eigene Dokumente\Uni\BA\BA_cardioMOR\Testdaten\Dynamisch\Unsteady_5_cycles_newflow'
folder_content = dir();
X_comp = zeros(vec_size,steps);
tic
k=4;
for i = 1:steps   
    %open and read file
    fileid = folder_content(k).name;
    data = readmatrix(fileid);
    snap = reshape(data(:,7:9), [],1);
    %store data in snapshot matrices
    X_comp(:,i) = snap;
    k=k+1;
end
toc

%% write whole matrix to memory efficient bin file
tic
id=fopen('5cn_all.bin','w');
fwrite(id,X_comp,'double');
fclose(id);
toc

%% reread bin file for correctness check
tic
id=fopen('5cn_all.bin');
XX_bin = fread(id,[vec_size steps],'double');
fclose(id);
toc
