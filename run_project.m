%%
close all;
clearvars; clc;

%% INPUT ARGS
disp('Load input args...');
tic;

input.dir = dir(fullfile(pwd, 'dataset_cipic\standard_hrir_database\'));
input.names = {input.dir(4:48).name}.';
input.Fs = 44100;
input.N = 256;
input.bases = [1,5,10,20,30,50,75,100,128];

varb.fnScale = -input.N/2 : input.N/2-1;
varb.fScale = (input.Fs/input.N) * varb.fnScale;
toc;

%% LOAD HRTF'S, CALCULATE GRAND MEAN
disp('Load HRTFs...')
tic;
query = horzcat(linspace(0,0,5).', (0:45:180).');
varb.azQuery = query(:,1);
varb.elQuery = query(:,2);
varb.sum = 0;

varb.subj2 = hrtf_load(input.names{12}, query, 1); % KEMAR large pinnae visualization
varb.subj2.handles(2).XLim = [3000,22050];
varb.subj2.handles(3).XLim = [3000,22050];

for sub = 1:numel(input.names)
%for sub = 1:3
    subj(sub) = hrtf_load(input.names{sub}, query, 0);
    varb.sum = varb.sum + squeeze(sum(subj(sub).DTF_l, [1,2])) + squeeze(sum(subj(sub).DTF_r, [1,2]));
end

subj(1).grandMean = varb.sum ./ (2 * numel(subj(1).DTF_l(:,:,1)) * numel(subj));

for sub = 1:numel(subj)
    subj(sub).DTF2_l = subj(sub).DTF_l ./ reshape(subj(1).grandMean, 1,1,numel(subj(1).grandMean));
    subj(sub).DTF2_r = subj(sub).DTF_r ./ reshape(subj(1).grandMean, 1,1,numel(subj(1).grandMean));
end

toc;

%% REDUCE

varb.fields = {'OnL','OnR','ITD','hrir_l','hrir_r'};
subj = rmfield(subj, varb.fields);


%% LOAD ANTHROPOMETRIC DATA
%{
anthro = anthro_load(0);
%}

%% CALCULATE COVARIANCE MATRICES AND GRAND SUMS
% each cell of cov_f represents an angle. each NxN matrix in the cell
% represents the covariance between 20LOG10 frequency magnitudes across subjects.

disp('Covariance and correlation...')
tic;

co.covar_f = cell(size(subj(1).DTF_l(:,:,1))); % '_f' means covariance between frequency bins
varb.covarSum = zeros(input.N/2, input.N/2);
varb.correSum = zeros(input.N/2, input.N/2);

for i=1:numel(co.covar_f(:,1))
    for j=1:numel(co.covar_f(1,:))	% for each angle...
        varb.A = zeros(2*numel(subj), input.N/2);
        for k = 1:numel(subj)           % populate rows of A with each subject's 20log10(DTF), L+R
            varb.A(k,:) = 20*log10(abs(squeeze(subj(k).DTF2_l( i, j, 1:(input.N/2) ))));
            varb.A(k+numel(subj),:) = 20*log10(abs(squeeze(subj(k).DTF2_r( i, j, 1:(input.N/2) ))));
        end
        co.covar_f{i,j} = cov(varb.A);  % and find the covariance between frequencies across subjects
        co.corre_f{i,j} = corrcov(co.covar_f{i,j}); 
        varb.covarSum = varb.covarSum + co.covar_f{i,j};
        varb.correSum = varb.correSum + co.corre_f{i,j};
    end
end
co.covarMean_f = varb.covarSum / numel(co.covar_f);
co.correMean_f = varb.correSum / numel(co.corre_f);

[eigen.vec, eigen.val] = eig(co.covarMean_f);
%[eigen.vec, eigen.val] = eig(co.correMean_f);

% more sanity check
%{
a = eigen.vec * eigen.val;
b = co.covarMean_f * eigen.vec;
c=a-b;
%}
toc;

%% sanity check covar/corr figure

figure;

subplot(2,2,1);
surf(varb.fScale(129:end), varb.fScale(129:end), co.covar_f{12,15}, 'EdgeAlpha', 0);
view(2); colorbar;
xlim([varb.fScale(130), varb.fScale(end)]);
ylim([varb.fScale(130), varb.fScale(end)]);
ax = gca; ax.XScale = 'log'; ax.YScale = 'log'; ax.TickDir = 'out';
title('Covariance sample');

subplot(2,2,2);
surf(varb.fScale(129:end), varb.fScale(129:end), co.covarMean_f, 'EdgeAlpha', 0);
view(2); colorbar;
xlim([varb.fScale(130), varb.fScale(end)]);
ylim([varb.fScale(130), varb.fScale(end)]);
ax = gca; ax.XScale = 'log'; ax.YScale = 'log'; ax.TickDir = 'both';
title('Covariance mean');

subplot(2,2,3);
surf(varb.fScale(129:end), varb.fScale(129:end), co.corre_f{12,15}, 'EdgeAlpha', 0);
view(2); colorbar;
xlim([varb.fScale(130), varb.fScale(end)]);
ylim([varb.fScale(130), varb.fScale(end)]);
xlabel('Frequency (Hz)');
ax = gca; ax.XScale = 'log'; ax.YScale = 'log'; ax.TickDir = 'out';
title('Correlation sample');

subplot(2,2,4);
surf(varb.fScale(129:end), varb.fScale(129:end), co.correMean_f, 'EdgeAlpha', 0);
view(2); colorbar;
xlim([varb.fScale(130), varb.fScale(end)]);
ylim([varb.fScale(130), varb.fScale(end)]);
xlabel('Frequency (Hz)');
ax = gca; ax.XScale = 'log'; ax.YScale = 'log'; ax.TickDir = 'out';
title('Correlation mean');

%% FIND STRONGEST EIGENVALUES AND DECOMPOSE

disp('Decomposing...');
tic;

for b = 1:numel(input.bases)
    eigen.vecProj = eigen.vec( :, (end-input.bases(b)+1):end );
    eigen.valProj = eigen.val( (end-input.bases(b)+1):end, (end-input.bases(b)+1):end ); % build reduced covar matrix

    for sub = 1:numel(subj) 
        
        subj(sub).DTF2Weights_l(b).data = cell( numel(subj(sub).DTF2_l(:,1,1)), numel(subj(sub).DTF2_l(1,:,1)), input.bases(b) );
        subj(sub).DTF2Weights_r(b).data = cell( numel(subj(sub).DTF2_r(:,1,1)), numel(subj(sub).DTF2_r(1,:,1)), input.bases(b) );
        subj(sub).DTF2synth_l(b).data = cell( numel(subj(sub).DTF2_r(:,1,1)), numel(subj(sub).DTF2_r(1,:,1)), input.N/2 );
        
        subj(sub).DTF2Weights_l(b).numBases = input.bases(b);
        subj(sub).DTF2Weights_r(b).numBases = input.bases(b);
        
        %{
        subj(sub).DTF2Weights_l(b) = zeros( numel(subj(sub).DTF2_l(:,1,1)), numel(subj(sub).DTF2_l(1,:,1)), input.bases(b) );
        subj(sub).DTF2Weights_r(b) = zeros( numel(subj(sub).DTF2_r(:,1,1)), numel(subj(sub).DTF2_r(1,:,1)), input.bases(b) );
        %}
        for i=1:numel(subj(sub).DTF2_l(:,1,1))
            for j=1:numel(subj(sub).DTF2_l(1,:,1))
                varb.weights_l = eigen.vecProj.' * 20*log10(abs(squeeze(subj(sub).DTF2_l(i,j,1:128))));
                varb.weights_r = eigen.vecProj.' * 20*log10(abs(squeeze(subj(sub).DTF2_r(i,j,1:128))));
                varb.synth_l = eigen.vecProj * varb.weights_l;
                varb.synth_r = eigen.vecProj * varb.weights_r;
                subj(sub).DTF2Weights_l(b).data{i,j} = varb.weights_l;
                subj(sub).DTF2Weights_r(b).data{i,j} = varb.weights_r;
                subj(sub).DTF2synth_l(b).data{i,j} = varb.synth_l;
                subj(sub).DTF2synth_r(b).data{i,j} = varb.synth_r;
                %{
                for k=1:(input.bases(b))
                    subj(sub).DTF2Weights_l(i,j,k) = varb.weights_l(k);
                    subj(sub).DTF2Weights_r(i,j,k) = varb.weights_r(k);
                end
                %}
                %{
                for k=1:(input.N/2)
                    subj(sub).DTF2synth_l(b).data(i,j,k) = varb.synth_l(k);
                    subj(sub).DTF2synth_r(b).data(i,j,k) = varb.synth_r(k);
                end
                %}
            end
        end
    end
end
toc;
%% sanity check figure

figure;
i=1;

subplot(3,2,1);
hold on;
varb.H_l = abs(fftshift(squeeze(subj(i).HRTF_l(12,5,:))));
plot(varb.fScale, 20*log10(varb.H_l), 'k', 'LineWidth', 1.5);
varb.H_l = abs(fftshift(squeeze(subj(i).meanHRTF_l)));
plot(varb.fScale, 20*log10(varb.H_l), 'k:', 'LineWidth', 1.5);
hold off;
varb.ax = gca;
varb.ax.XScale = 'log';
xlim([100,input.Fs/2]);
grid on;
legend('L. HRTF', 'L. mean', 'Location', 'southwest');
title('HRTF');

subplot(3,2,3);
hold on;
varb.H_l = abs(fftshift(squeeze(subj(i).DTF_l(12,5,:))));
plot(varb.fScale, 20*log10(varb.H_l), 'k', 'LineWidth', 1.5);
varb.H_l = abs(fftshift(squeeze(subj(1).grandMean)));
plot(varb.fScale, 20*log10(varb.H_l), 'k:', 'LineWidth', 1.5);
hold off;
varb.ax = gca;
varb.ax.XScale = 'log';
xlim([100,input.Fs/2]);
grid on;
legend('L. DTF', 'LR Grand mean', 'Location', 'southwest');
title('DTF = (HRTF - Lmean)');

subplot(3,2,5);
hold on;
varb.H_l = abs(fftshift(squeeze(subj(i).DTF2_l(12,5,:))));
plot(varb.fScale, 20*log10(varb.H_l), 'k', 'LineWidth', 1.5, 'DisplayName', 'L. DTF2');
%{
for b = 1:numel(input.bases)
    varb.H_l = subj(i).DTF2synth_l(b).data{12,5};
    varb.dispName = strcat(num2str(input.bases(b)), ' bases');
    plot(varb.fScale(129:end), varb.H_l, 'DisplayName', varb.dispName);
end
%}
hold off;
varb.ax = gca;
varb.ax.XScale = 'log';
xlim([100,input.Fs/2]);
grid on;
xlabel('Frequency (Hz)');
legend('Location', 'southwest');
title('DTF2 = (DTF - Grand Mean)');

subplot(3,2,2);
hold on;
varb.H_r = abs(fftshift(squeeze(subj(i).HRTF_r(12,5,:))));
plot(varb.fScale, 20*log10(varb.H_r), 'k', 'LineWidth', 1.5);
varb.H_r = abs(fftshift(squeeze(subj(i).meanHRTF_r)));
plot(varb.fScale, 20*log10(varb.H_r), 'k:', 'LineWidth', 1.5);
hold off;
varb.ax = gca;
varb.ax.XScale = 'log';
xlim([100,input.Fs/2]);
grid on;
legend('R. HRTF', 'R. mean', 'Location', 'southwest');
title('HRTF');

subplot(3,2,4);
hold on;
varb.H_r = abs(fftshift(squeeze(subj(i).DTF_r(12,5,:))));
plot(varb.fScale, 20*log10(varb.H_r), 'k', 'LineWidth', 1.5);
varb.H_r = abs(fftshift(squeeze(subj(1).grandMean)));
plot(varb.fScale, 20*log10(varb.H_r), 'k:', 'LineWidth', 1.5);
hold off;
varb.ax = gca;
varb.ax.XScale = 'log';
xlim([100,input.Fs/2]);
grid on;
legend('R. DTF', 'LR Grand mean', 'Location', 'southwest');
title('DTF = (HRTF - Lmean)');

subplot(3,2,6);
hold on;
varb.H_r = abs(fftshift(squeeze(subj(i).DTF2_r(12,5,:))));
plot(varb.fScale, 20*log10(varb.H_r), 'k', 'LineWidth', 1.5, 'DisplayName', 'R. DTF2');
%{
for b = 1:numel(input.bases)
    varb.H_r = subj(i).DTF2synth_r(b).data{12,5};
    varb.dispName = strcat(num2str(input.bases(b)), ' bases');
    plot(varb.fScale(129:end), varb.H_r, 'DisplayName', varb.dispName);
end
%}
hold off;
varb.ax = gca;
varb.ax.XScale = 'log';
xlim([100,input.Fs/2]);
grid on;
xlabel('Frequency (Hz)');
legend('Location', 'southwest');
title('DTF2 = (DTF - Grand Mean)');

%% HRTF2 = ADD MEANS BACK INTO DTF2
disp('Reconstructing HRTFs...');
tic;
for b = 1:numel(input.bases)
    for sub = 1:numel(subj)
        subj(sub).HRTF2_l(b).data = zeros(size(subj(sub).HRTF_l(:,:,1), input.N/2)); % initialize container
        subj(sub).HRTF2_r(b).data = zeros(size(subj(sub).HRTF_r(:,:,1), input.N/2));
        subj(sub).MSE_l(b).data = zeros(size(subj(sub).HRTF_l(:,:,1)));
        subj(sub).MSE_r(b).data = zeros(size(subj(sub).HRTF_r(:,:,1)));
        subj(sub).L2_l(b).data = zeros(size(subj(sub).HRTF_l(:,:,1)));
        subj(sub).L2_r(b).data = zeros(size(subj(sub).HRTF_r(:,:,1)));
        for i = 1:numel(subj(sub).HRTF_l(:,1,1))
            for j = 1:numel(subj(sub).HRTF_l(1,:,1))
                varb.backin_l = subj(sub).DTF2synth_l(b).data{i,j} + 20*log10(abs(subj(1).grandMean(1:128))) + 20*log10(abs(subj(sub).meanHRTF_l(1:128)));
                varb.backin_r = subj(sub).DTF2synth_r(b).data{i,j} + 20*log10(abs(subj(1).grandMean(1:128))) + 20*log10(abs(subj(sub).meanHRTF_r(1:128)));
                varb.SE_l = sum((varb.backin_l - squeeze(20*log10(abs(subj(sub).HRTF_l(i,j,1:128))))).^2);
                varb.SE_r = sum((varb.backin_r - squeeze(20*log10(abs(subj(sub).HRTF_r(i,j,1:128))))).^2);
                varb.MSE_l = (1/128) * varb.SE_l; 
                varb.MSE_r = (1/128) * varb.SE_r;
                varb.L2_l = sqrt(varb.SE_l);
                varb.L2_r = sqrt(varb.SE_r);
                subj(sub).MSE_l(b).data(i,j) = varb.MSE_l; % MSE of dB
                subj(sub).MSE_r(b).data(i,j) = varb.MSE_r; % MSE of dB
                subj(sub).L2_l(b).data(i,j) = varb.L2_l; % L2 Loss of dB
                subj(sub).L2_r(b).data(i,j) = varb.L2_r; % L2 Loss of dB
                for k = 1:(input.N/2)
                    subj(sub).HRTF2_l(b).data(i,j,k) = varb.backin_l(k);
                    subj(sub).HRTF2_r(b).data(i,j,k) = varb.backin_r(k);
                end
            end
        end
    end
end
toc;

%% SHOW PRINCIPAL COMPONENTS

figure;
input.show = 5;

for b = 1:input.show
    subplot(5,1,b);
    plot(varb.fScale(129:end), eigen.vecProj(:,end-b+1), 'k', 'LineWidth', 1.5);
    varb.ax = gca;
    varb.ax.XScale = 'log';
    xlim([100,input.Fs/2]);
    grid on;
end

subplot(5,1,1); title('Components');
subplot(5,1,5); xlabel('Frequency (Hz)');

%%
figure;

i=1;
subplot(2,2,1);
hold on;
varb.H_l = abs(fftshift(squeeze(subj(i).HRTF_l(12,5,:))));
plot(varb.fScale, 20*log10(varb.H_l), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Original HRTF');
%for b = 1:numel(input.bases)
for b = [1,3,5,6,8]
    varb.H_l = squeeze(subj(i).HRTF2_l(b).data(12,5,:));
    varb.dispName = strcat(num2str(input.bases(b)), ' bases');
    plot(varb.fScale(129:end), varb.H_l, 'DisplayName', varb.dispName);
end
hold off;
varb.ax = gca;
varb.ax.XScale = 'log';
xlim([100,input.Fs/2]);
grid on;
xlabel('Frequency (Hz)');
legend('Location', 'southwest');
title('HRTF2 = (DTF2 + Grand Mean + Lmean)');

subplot(2,2,2);
hold on;
varb.H_r = abs(fftshift(squeeze(subj(i).HRTF_r(12,5,:))));
plot(varb.fScale, 20*log10(varb.H_r), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Original HRTF');
%for b = 1:numel(input.bases)
for b = [1,3,5,6,8]
    varb.H_r = squeeze(subj(i).HRTF2_r(b).data(12,5,:));
    varb.dispName = strcat(num2str(input.bases(b)), ' bases');
    plot(varb.fScale(129:end), varb.H_r, 'DisplayName', varb.dispName);
end
hold off;
varb.ax = gca;
varb.ax.XScale = 'log';
xlim([100,input.Fs/2]);
grid on;
xlabel('Frequency (Hz)');
legend('Location', 'southwest');
title('HRTF2 = (DTF2 + Grand Mean + Rmean)');

subplot(2,2,3);
varb.holdMSE = zeros(numel(input.bases), 1);
varb.holdL2 = zeros(numel(input.bases), 1);
for b = 1:numel(input.bases)
    varb.holdMSE(b) = subj(i).MSE_l(b).data(12,5);
    varb.holdL2(b) = subj(i).L2_l(b).data(12,5);
end
hold on;
plot(input.bases, varb.holdMSE, 'k-*', 'DisplayName', 'MSE', 'LineWidth', 1.5);
plot(input.bases, varb.holdL2, 'k-x', 'DisplayName', 'L2 Norm', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('Number of Bases');
legend('Location', 'northeast');
title('Left Loss');

subplot(2,2,4);
varb.holdMSE = zeros(numel(input.bases), 1);
varb.holdL2 = zeros(numel(input.bases), 1);
for b = 1:numel(input.bases)
    varb.holdMSE(b) = subj(i).MSE_r(b).data(12,5);
    varb.holdL2(b) = subj(i).L2_r(b).data(12,5);
end
hold on;
plot(input.bases, varb.holdMSE, 'k-*', 'DisplayName', 'MSE', 'LineWidth', 1.5);
plot(input.bases, varb.holdL2, 'k-x', 'DisplayName', 'L2 Norm', 'LineWidth', 1.5);
hold off;
grid on;
xlabel('Number of Bases');
legend('Location', 'northeast');
title('Right Loss');
