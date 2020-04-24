function subj = hrtf_load(in_subj, in_query, showFigure) % showFigure: 0=no figure; 1=layout1; 2=layout2
%% LOAD FILES

    % FOR TESTING
    %{
    in_subj = 'subject_003';
    in_query = [[0,0]; [0,45]; [0,90]]; % [[az1, el1]; [az2, el2]; ...]
    showFigure = 1;
    %}


    input.directory = fullfile(pwd, 'dataset_cipic\standard_hrir_database\');
    input.subject = in_subj;
    input.query = in_query;
    input.N = 256;
    input.Fs = 44100;
    
    if isempty(in_query)
       input.query = [0,0];
    end
    
    subj = load(fullfile(input.directory, input.subject, '\hrir_final.mat'));
    varb.azQuery = input.query(:,1);
    varb.elQuery = input.query(:,2);

%% CALCULATE IID ( sum(H(t).^2) ) AND HRTF
    subj.IID_l = zeros(size(subj.hrir_l(:,:,1)));
    subj.IID_r = zeros(size(subj.hrir_r(:,:,1)));
    subj.HRTF_l = complex(zeros([size(subj.hrir_l(:,:,1)), input.N]));
    subj.HRTF_r = complex(zeros([size(subj.hrir_r(:,:,1)), input.N]));
    varb.fft_l = complex(zeros(input.N, 1));
    varb.fft_r = complex(zeros(input.N, 1));

    subj.IID_l = sum(subj.hrir_l.^2, 3);
    subj.IID_r = sum(subj.hrir_r.^2, 3);

    subj.HRTF_l = fft(subj.hrir_l, input.N, 3);
    subj.HRTF_r = fft(subj.hrir_r, input.N, 3);

% OLD    
%{    
    for i = 1:numel(subj.hrir_l(:,1,1))
        for j = 1:numel(subj.hrir_l(i,:,1))
            varb.fft_l = fft(squeeze(subj.hrir_l(i,j,:)),input.N);
            varb.fft_r = fft(squeeze(subj.hrir_r(i,j,:)),input.N);
            for k = 1:input.N
                subj.HRTF_l(i,j,k) = varb.fft_l(k);
                subj.HRTF_r(i,j,k) = varb.fft_r(k);
            end
        end
    end
%}
%% CALCULATE PER-SUBJECT PER-EAR MEANS
    subj.meanHRTF_l = complex(zeros(input.N,1));
    subj.meanHRTF_r = complex(zeros(input.N,1));
    subj.DTF_l = complex(zeros([size(subj.HRTF_l(:,:,1)), input.N]));
    subj.DTF_r = complex(zeros([size(subj.HRTF_r(:,:,1)), input.N]));
    subj.meanHRTF = complex(zeros(input.N,1));
    
    for i=1:numel(subj.HRTF_l(1,1,:)) % for each fft bin
%         varb.sum_l = sum( subj.HRTF_l(:,:,i) , [1,2]); % unweighted by IID
%         varb.sum_r = sum( subj.HRTF_r(:,:,i) , [1,2]);

        varb.sum_l = sum( subj.HRTF_l(:,:,i) .* subj.IID_l , [1,2]); % weighted by IID
        varb.sum_r = sum( subj.HRTF_r(:,:,i) .* subj.IID_r , [1,2]);

        subj.meanHRTF_l(i) = varb.sum_l / numel(subj.HRTF_l(:,:,i));
        subj.meanHRTF_r(i) = varb.sum_r / numel(subj.HRTF_r(:,:,i));
        
        subj.DTF_l(:,:,i) = subj.HRTF_l(:,:,i) ./ subj.meanHRTF_l(i);
        subj.DTF_r(:,:,i) = subj.HRTF_r(:,:,i) ./ subj.meanHRTF_r(i);
        
%         subj.meanHRTF(i) = (subj.meanHRTF_l(i)+subj.meanHRTF_r(i)) / 2; % in case u need it
    end

%% PREP VARB
    varb.xScale = 1:25;
    varb.yScale = 1:50;
    varb.azScale = [-80, -65, -55, -45:5:45, 55, 65, 80];
    varb.elScale = -45 + 5.625*(0:49);
    varb.fnScale = -input.N/2 : input.N/2-1;
    varb.fScale = (input.Fs/input.N) * varb.fnScale;
    [varb.coordx, varb.coordy] = meshgrid(varb.azScale, varb.elScale);
    for i = 1:numel(varb.azQuery)
        [varb.xError(i), varb.xQuery(i)] = min(abs(varb.azQuery(i) - varb.azScale));
        [varb.yError(i), varb.yQuery(i)] = min(abs(varb.elQuery(i) - varb.elScale));
    end
    varb.azQueryActual = varb.azScale(varb.xQuery);
    varb.elQueryActual = varb.elScale(varb.yQuery);
    varb.legend = vertcat(varb.azQueryActual, varb.elQueryActual).';
    
    %% PLOT FIGURES
    if (showFigure)
        %%
        figure;

        subj.handles(1) = subplot(2,3,1);
        hold on;
        for i=1:numel(varb.xQuery)
           plot(varb.azScale, (subj.OnR(:,varb.yQuery(i)) - subj.OnL(:,varb.yQuery(i))), 'LineWidth', 1.5);
        end
        hold off;
        grid on;
        xlim([min(varb.azScale), max(varb.azScale)]);
        varb.leg = legend(num2str(varb.legend));
        title(varb.leg, strcat(input.subject, ' (az el)'));
        varb.leg.Title.Visible = 1;
        varb.leg.Title.Interpreter = 'none';
        xlabel('Azimuth (\circ)'); ylabel('ITD, L to R (ms)');
        title('ITD Slice');

        %%
        subj.handles(2) = subplot(4,3,2);
        if (showFigure == 1) % if 1, show DTF
            hold on;
            for i=1:numel(varb.xQuery)
                varb.H_l = abs(fftshift(squeeze(subj.DTF_l(varb.xQuery(i),varb.yQuery(i),:))));
                plot(varb.fScale, 20*log10(varb.H_l), 'LineWidth', 1.5);
            end      
            hold off;
            varb.ax = gca;
            varb.ax.XScale = 'log';
            xlim([100,input.Fs/2]);
            ylabel('dB Amplitude (?/bin)');
            grid on;
            title('(HRTF - Avg.), Left');
        end
        if (showFigure == 2) % if 2, show HRIR
            varb.max = 0;
            hold on;
            for i=1:numel(varb.xQuery)
                varb.ir = squeeze(subj.hrir_l(varb.xQuery(i),varb.yQuery(i),:));
                plot(1:200, varb.ir, 'LineWidth', 1.5);
                if max(abs(varb.ir)) > varb.max
                    varb.max = max(abs(varb.ir));
                end
            end
            hold off;
            ylim([-varb.max, varb.max]);
            ylabel('Amplitude (?)');
            grid on;
            title('HRIR, Left');
        end

        subj.handles(3) = subplot(4,3,3);
        if (showFigure == 1) % if 1, show DTF
            hold on;
            for i=1:numel(varb.xQuery)
                %varb.H_r = abs(fftshift(squeeze(subj.HRTF_r(varb.xQuery(i),varb.yQuery(i),:)) ./ subj.meanHRTF_r));
                varb.H_r = abs(fftshift(squeeze(subj.DTF_r(varb.xQuery(i),varb.yQuery(i),:))));
                plot(varb.fScale, 20*log10(varb.H_r), 'LineWidth', 1.5);
            end      
            hold off;
            varb.ax = gca;
            varb.ax.XScale = 'log';
            xlim([100,input.Fs/2]);
            ylabel('dB Amplitude (?/bin)');
            grid on;
            title('(HRTF - Avg.), Right');
        end
        if (showFigure == 2) % if 2, show HRIR
            varb.max = 0;
            hold on;
            for i=1:numel(varb.xQuery)
                varb.ir = squeeze(subj.hrir_r(varb.xQuery(i),varb.yQuery(i),:));
                plot(1:200, varb.ir, 'linewidth', 1.5);
                if max(abs(varb.ir)) > varb.max
                    varb.max = max(abs(varb.ir));
                end
            end            
            hold off;
            ylim([-varb.max, varb.max]);
            ylabel('Amplitude (?)');
            grid on;
            title('HRIR, Right');
        end

        %%        
        subj.handles(4) = subplot(4,3,5);
        hold on;

        varb.H_l = abs(fftshift(squeeze(subj.meanHRTF_l))); % plot mean
        plot(varb.fScale, 20*log10(varb.H_l), 'k--', 'LineWidth', 1.5);    
        varb.H_l = abs(fftshift(squeeze(subj.meanHRTF))); % plot mean
        plot(varb.fScale, 20*log10(varb.H_l), 'k:', 'LineWidth', 1.5);    
        for i=1:numel(varb.xQuery)
            varb.H_l = abs(fftshift(squeeze(subj.HRTF_l(varb.xQuery(i),varb.yQuery(i),:))));
            plot(varb.fScale, 20*log10(varb.H_l), 'LineWidth', 1.5);
        end    

        hold off;
        varb.ax = gca;
        varb.ax.XScale = 'log';
        xlim([100,input.Fs/2]);
        xlabel('Frequency'); ylabel('dB Amplitude (?/bin)');
        grid on;
        legend('L avg', 'LR avg', 'Location', 'southwest');
        title('HRTF, Left');

        %%
        subj.handles(5) = subplot(4,3,6);
        hold on;
        
        varb.H_r = abs(fftshift(squeeze(subj.meanHRTF_r)));
        plot(varb.fScale, 20*log10(varb.H_r), 'k--', 'LineWidth', 1.5);        
        varb.H_r = abs(fftshift(squeeze(subj.meanHRTF)));
        plot(varb.fScale, 20*log10(varb.H_r), 'k:', 'LineWidth', 1.5);        
        for i=1:numel(varb.xQuery)
            varb.H_r = abs(fftshift(squeeze(subj.HRTF_r(varb.xQuery(i),varb.yQuery(i),:))));
            plot(varb.fScale, 20*log10(varb.H_r), 'LineWidth', 1.5);
        end
        
        hold off;
        varb.ax = gca;
        varb.ax.XScale = 'log';
        xlim([100,input.Fs/2]);
        xlabel('Frequency'); ylabel('dB Amplitude (?/bin)');
        grid on;
        legend('R avg', 'LR avg.', 'Location', 'southwest');
        title('HRTF, Right');

        %%
        subj.handles(6) = subplot(2,3,4);
        surf(varb.coordx.', varb.coordy.', subj.OnR-subj.OnL, 'EdgeAlpha', 0);
        hold on;
        for i=1:numel(varb.xQuery)
            stem3(varb.azScale(varb.xQuery(i)), varb.elScale(varb.yQuery(i)), max(max(subj.ITD)), 'o', 'LineWidth', 1.5);
        end
        hold off;
        xlim([min(varb.azScale), max(varb.azScale)]); ylim([min(varb.elScale), max(varb.elScale)]);
        view(2);
        colormap(gca, gray(30));
        varb.c = colorbar;
        varb.c.Label.String = 'ITD, L to R (ms)';
        xlabel('Azimuth (\circ)'); ylabel('Elevation (\circ)');
        title('ITD');
%%
        subj.handles(7) = subplot(2,3,5);
        surf(varb.coordx.', varb.coordy.', 20*log10(subj.IID_l), 'EdgeAlpha', 0);
        hold on;
        view(2)
        for i=1:numel(varb.xQuery)
            stem3(varb.azScale(varb.xQuery(i)), varb.elScale(varb.yQuery(i)), max(max(subj.ITD)), 'o', 'LineWidth', 1.5);
        end
        hold off;
        xlim([min(varb.azScale), max(varb.azScale)]); ylim([min(varb.elScale), max(varb.elScale)]);
        colormap(gca, gray(30));
        varb.c = colorbar;
        varb.c.Label.String = 'dB Energy (J)';
        xlabel('Azimuth (\circ)');
        ylabel('Elevation (\circ)');
        title('IID, Left');

        subj.handles(8) = subplot(2,3,6);
        surf(varb.coordx.', varb.coordy.', 20*log10(subj.IID_r), 'EdgeAlpha', 0);
        hold on;
        view(2);
        for i=1:numel(varb.xQuery)
            stem3(varb.azScale(varb.xQuery(i)), varb.elScale(varb.yQuery(i)), max(max(subj.ITD)), 'o', 'LineWidth', 1.5);
        end
        hold off;
        xlim([min(varb.azScale), max(varb.azScale)]); ylim([min(varb.elScale), max(varb.elScale)]);
        colormap(gca, gray(30));
        varb.c = colorbar;
        varb.c.Label.String = 'dB Energy (J)';
        title('IID, Right');
        xlabel('Azimuth (\circ)');
        ylabel('Elevation (\circ)');
    end
end