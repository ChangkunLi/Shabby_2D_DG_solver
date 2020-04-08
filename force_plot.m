clear;

citymax = 0;
pmax = 0;

for icity = 0:citymax
    Data = cell(pmax+1,1);
    for p = 0:pmax
        filename = ['city' num2str(icity) '_p' num2str(p) '_Force.dat'];
        f = fopen(filename);
        force = fscanf(f,'%f');
        fclose(f);
        total = size(force,1);
        force = reshape(force,[8,total/8]);
        force = force';
        Data{p+1} = force;
    end

    for i = 1:4
        figure;  % Fx

        for p = 0:pmax
            Nt = size(Data{p+1},1);
            t = linspace(0,2,Nt);
            plot(t, Data{p+1}(:,2*i-1),'LineWidth',3,'DisplayName',['p = ' num2str(p)]);
            hold on;
        end

        legend;
        xlabel('T');
        ylabel('F_x (N)');
        title(['city' num2str(icity) ', B' num2str(i)]);
        ax = gca;
        ax.FontSize = 20;

        figure;  % Fy

        for p = 0:pmax
            Nt = size(Data{p+1},1);
            t = linspace(0,2,Nt);
            plot(t, Data{p+1}(:,2*i),'LineWidth',3,'DisplayName',['p = ' num2str(p)]);
            hold on;
        end    

        legend;
        xlabel('T');
        ylabel('F_y (N)');
        title(['city' num2str(icity) ', B' num2str(i)]);
        ax = gca;
        ax.FontSize = 20;    
    end
end

