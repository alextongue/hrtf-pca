function out = anthro_load(showFigure)
    %% PLOT THETA AND D, L-R SCATTERS
    anthro = load('D:\Desktop\FA19\ECE225A\DB_CIPIC\anthropometry\anthro.mat');
    if (showFigure==1)
        figure;

        anthro.handles(1) = subplot(2,5,1);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.theta(i,1), anthro.theta(i,3), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.theta(i,1), anthro.theta(i,3), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('\theta_1 Left');
        ylabel('\theta_1 Right');
        axis square;
        title('Pinna Rotation Angle');

        anthro.handles(2) = subplot(2,5,6);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.theta(i,2), anthro.theta(i,4), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.theta(i,2), anthro.theta(i,4), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('\theta_2 Left');
        ylabel('\theta_2 Right');
        axis square;
        title('Pinna Flare Angle');

        anthro.handles(3) = subplot(2,5,2);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.D(i,1), anthro.D(i,9), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.D(i,1), anthro.D(i,9), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('D_1 Left');
        ylabel('D_1 Right');
        axis square;
        title('Cavum Concha Height');

        anthro.handles(4) = subplot(2,5,3);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.D(i,2), anthro.D(i,10), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.D(i,2), anthro.D(i,10), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('D_2 Left');
        ylabel('D_2 Right');
        axis square;
        title('Cymba Concha Height');

        anthro.handles(5) = subplot(2,5,4);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.D(i,3), anthro.D(i,11), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.D(i,3), anthro.D(i,11), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('D_3 Left');
        ylabel('D_3 Right');
        axis square;
        title('Cavum Concha Width');

        anthro.handles(6) = subplot(2,5,5);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.D(i,4), anthro.D(i,12), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.D(i,4), anthro.D(i,12), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('D_4 Left');
        ylabel('D_4 Right');
        axis square;
        title('Fossa Height');

        anthro.handles(7) = subplot(2,5,7);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.D(i,5), anthro.D(i,13), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.D(i,5), anthro.D(i,13), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('D_5 Left');
        ylabel('D_5 Right');
        axis square;
        title('Pinna Height');

        anthro.handles(8) = subplot(2,5,8);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.D(i,6), anthro.D(i,14), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.D(i,6), anthro.D(i,14), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('D_6 Left');
        ylabel('D_6 Right');
        axis square;
        title('Pinna Width');

        anthro.handles(9) = subplot(2,5,9);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.D(i,7), anthro.D(i,15), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.D(i,7), anthro.D(i,15), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('D_7 Left');
        ylabel('D_7 Right');
        axis square;
        title('Intertragal Incisure Width');

        anthro.handles(10) = subplot(2,5,10);
        hold on;
        for i=1:numel(anthro.sex)
            switch anthro.sex(i)
                case 'M'
                    scatter(anthro.D(i,8), anthro.D(i,16), 'ro', 'LineWidth', 1.5);
                case 'F'
                    scatter(anthro.D(i,8), anthro.D(i,16), 'bo', 'LineWidth', 1.5);
            end
        end
        hold off;
        grid on;
        xlabel('D_8 Left');
        ylabel('D_8 Right');
        axis square;
        title('Cavum Concha Depth');
    end
    out = anthro;
    
    %% PLOT X'S
end
