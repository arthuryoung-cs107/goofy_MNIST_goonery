clear
close all

fig_pos_old = [10, 455, 350, 350; 360, 455, 350, 350; 710, 455, 350, 350; 1060, 455, 350, 350; 10, 30, 350, 350; 360, 30, 350, 350; 710, 30, 350, 350; 1060, 30, 350, 350];
fig_pos = [1 700 173 130; 175 700 173 130; 349 700 173 130; 523 700 173 130; 697 700 173 130; 871 700 173 130; 1045 700 173 130; 1219 700 173 130];
fig_pos1 = [1 498 173 130; 175 498 173 130; 349 498 173 130; 523 498 173 130; 697 498 173 130; 871 498 173 130; 1045 498 173 130; 1219 498 173 130];
fig_pos2 = [1 300 173 130; 175 300 173 130; 349 300 173 130; 523 300 173 130; 697 300 173 130; 871 300 173 130; 1045 300 173 130; 1219 300 173 130];
fig_pos3 = [1 100 173 130; 175 100 173 130; 349 100 173 130; 523 100 173 130; 697 100 173 130; 871 100 173 130; 1045 100 173 130; 1219 100 173 130];


figure1 = figure('Name', 'view1', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'view2', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'view3', 'Renderer', 'painters', 'Position', fig_pos(3,:));
figure4 = figure('Name', 'view4', 'Renderer', 'painters', 'Position', fig_pos(4,:));
figure5 = figure('Name', 'view5', 'Renderer', 'painters', 'Position', fig_pos(5,:));
figure6 = figure('Name', 'view6', 'Renderer', 'painters', 'Position', fig_pos(6,:));
figure7 = figure('Name', 'view7', 'Renderer', 'painters', 'Position', fig_pos(7,:));
figure8 = figure('Name', 'view8', 'Renderer', 'painters', 'Position', fig_pos(8,:));

figure9 = figure('Name', 'view9', 'Renderer', 'painters', 'Position', fig_pos1(1,:));
figure10 = figure('Name', 'view10', 'Renderer', 'painters', 'Position', fig_pos1(2,:));
figure11 = figure('Name', 'view11', 'Renderer', 'painters', 'Position', fig_pos1(3,:));
figure12 = figure('Name', 'view12', 'Renderer', 'painters', 'Position', fig_pos1(4,:));
figure13 = figure('Name', 'view13', 'Renderer', 'painters', 'Position', fig_pos1(5,:));
figure14 = figure('Name', 'view14', 'Renderer', 'painters', 'Position', fig_pos1(6,:));
figure15 = figure('Name', 'view15', 'Renderer', 'painters', 'Position', fig_pos1(7,:));
figure16 = figure('Name', 'view16', 'Renderer', 'painters', 'Position', fig_pos1(8,:));

figure17 = figure('Name', 'view9', 'Renderer', 'painters', 'Position', fig_pos2(1,:));
figure18 = figure('Name', 'view10', 'Renderer', 'painters', 'Position', fig_pos2(2,:));
figure19 = figure('Name', 'view11', 'Renderer', 'painters', 'Position', fig_pos2(3,:));
figure20 = figure('Name', 'view12', 'Renderer', 'painters', 'Position', fig_pos2(4,:));
figure21 = figure('Name', 'view13', 'Renderer', 'painters', 'Position', fig_pos2(5,:));
figure22 = figure('Name', 'view14', 'Renderer', 'painters', 'Position', fig_pos2(6,:));
figure23 = figure('Name', 'view15', 'Renderer', 'painters', 'Position', fig_pos2(7,:));
figure24 = figure('Name', 'view16', 'Renderer', 'painters', 'Position', fig_pos2(8,:));

figure25 = figure('Name', 'view9', 'Renderer', 'painters', 'Position', fig_pos3(1,:));
figure26 = figure('Name', 'view10', 'Renderer', 'painters', 'Position', fig_pos3(2,:));
figure27 = figure('Name', 'view11', 'Renderer', 'painters', 'Position', fig_pos3(3,:));
figure28 = figure('Name', 'view12', 'Renderer', 'painters', 'Position', fig_pos3(4,:));
figure29 = figure('Name', 'view13', 'Renderer', 'painters', 'Position', fig_pos3(5,:));
figure30 = figure('Name', 'view14', 'Renderer', 'painters', 'Position', fig_pos3(6,:));
figure31 = figure('Name', 'view15', 'Renderer', 'painters', 'Position', fig_pos3(7,:));
figure32 = figure('Name', 'view16', 'Renderer', 'painters', 'Position', fig_pos3(8,:));

mnist_it = 0;
mnist_in_raw = readmatrix(['../MNist_csv_raw/MNIST' num2str(mnist_it) '.csv']);
mnist_in_N = readmatrix(['../MNist_csv_50pcentMasked/MNIST' num2str(mnist_it) '_N.csv']);
mnist_in_RR = readmatrix(['../MNist_csv_50pcentMasked_RR/MNIST' num2str(mnist_it) '_RR.csv']);
mnist_in_CRR = readmatrix(['../MNist_csv_CRR/MNIST' num2str(mnist_it) '_CRR.csv']);


X_raw = mnist_in_raw(2:size(mnist_in_raw, 1),2:size(mnist_in_raw, 2) );
X_N = mnist_in_N;
X_RR = mnist_in_RR;
X_CRR = (mnist_in_CRR(2:size(mnist_in_CRR, 1),2:size(mnist_in_CRR, 2) ));


% 
figure(figure1.Number)
imshow(X_raw)

figure(figure2.Number)
imshow(X_N)

figure(figure3.Number)
imshow(X_RR)

figure(figure4.Number)
imshow(X_CRR)

% X_all = 
X_CRR_1 = zeros([8, size(mnist_in_RR)]);
X_CRR_2 = zeros([8, size(mnist_in_RR)]);
X_CRR_3= zeros([8, size(mnist_in_RR)]);

mu_count = readmatrix('../C_code/MNist_interim/interim_spec_file.dat');

for i = 1:mu_count
    raw1 = readmatrix(['../C_code/MNist_interim/MNIST' num2str(mnist_it) '_CRR_1_' num2str(i) '.csv']);
    raw2 = readmatrix(['../C_code/MNist_interim/MNIST' num2str(mnist_it) '_CRR_2_' num2str(i) '.csv']);
    raw3 = readmatrix(['../C_code/MNist_interim/MNIST' num2str(mnist_it) '_CRR_3_' num2str(i) '.csv']);
    
    X_CRR_1(i, :, :) = (raw1(2:size(mnist_in_CRR, 1),2:size(mnist_in_CRR, 2) ));
    X_CRR_2(i, :, :) = (raw2(2:size(mnist_in_CRR, 1),2:size(mnist_in_CRR, 2) ));
    X_CRR_3(i, :, :) = (raw3(2:size(mnist_in_CRR, 1),2:size(mnist_in_CRR, 2) ));

end

figure_it = figure8.Number;

for i = 1:mu_count
    figure_it = figure_it + 1;
    figure(figure_it)
    imshow(reshape(X_CRR_1(i, :, :), size(mnist_in_RR) ));
end

for i = 1:mu_count
    figure_it = figure_it + 1;
    figure(figure_it)
    imshow(reshape(X_CRR_2(i, :, :), size(mnist_in_RR) ));
end

for i = 1:mu_count
    figure_it = figure_it + 1;
    figure(figure_it)
    imshow(reshape(X_CRR_3(i, :, :), size(mnist_in_RR) ));
end
