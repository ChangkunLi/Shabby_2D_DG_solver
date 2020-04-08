clear;

screen = get(0,'ScreenSize');
W_g = screen(3); H_g = screen(4);

for p = 1:1
    for j = 0:5
        filename = ['city0_p' num2str(p) '_State_' num2str(j) '.dat'];
        f = fopen(filename);
        Data = fscanf(f,'%f');
        fclose(f);
        total = size(Data,1);
        Data = reshape(Data,[9, total/9]);
        X = Data([1,4,7],:);
        Y = Data([2,5,8],:);
        H = Data([3,6,9],:);
        
        w_g = 0.5*W_g;              % width of the window
        h_g = 0.5*H_g;              % hight of the window
        
        figure('Color',[1 1 1],'Position',[0,0,w_g,h_g]);
        %patch(X,Y,H); 
        patch(X,Y,H,'EdgeColor','none');    
        axis equal;
        axis tight;
        xlabel('X');
        ylabel('Y');
        title(['city0, p = ' num2str(p) ', t = ' num2str(j*0.05)]);
        ax = gca;
        ax.FontSize = 20;
        colorbar;
        colormap jet;
        if j == 0
            caxis([0.7 1.3]);
        else
            caxis([0.9 1.1]);
        end
    end
end