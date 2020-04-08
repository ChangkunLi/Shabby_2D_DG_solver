clear;

%% Full State boundary
filename = 'fullstate_test_p0_Res.dat';
f = fopen(filename);
Residual = fscanf(f,'%f');
fclose(f);

figure;
semilogy(abs(Residual),'LineWidth',3);
xlabel('Number of iteration');
ylabel('Residual');
title('Full State boundary condition, p=0');
ax = gca;
ax.FontSize = 20;

filename = 'fullstate_test_p1_Res.dat';
f = fopen(filename);
Residual = fscanf(f,'%f');
fclose(f);

figure;
semilogy(abs(Residual),'LineWidth',3);
xlabel('Number of iteration');
ylabel('Residual');
title('Full State boundary condition, p=1');
ax = gca;
ax.FontSize = 20;

filename = 'fullstate_test_p2_Res.dat';
f = fopen(filename);
Residual = fscanf(f,'%f');
fclose(f);

figure;
semilogy(abs(Residual),'LineWidth',3);
xlabel('Number of iteration');
ylabel('Residual');
title('Full State boundary condition, p=2');
ax = gca;
ax.FontSize = 20;


%% Wall boundary

filename = 'wall_test_p0_Res.dat';
f = fopen(filename);
Residual = fscanf(f,'%f');
fclose(f);

figure;
semilogy(abs(Residual),'LineWidth',3);
xlabel('Number of iteration');
ylabel('Residual');
title('Wall boundary condition, p=0');
ax = gca;
ax.FontSize = 20;

filename = 'wall_test_p1_Res.dat';
f = fopen(filename);
Residual = fscanf(f,'%f');
fclose(f);

figure;
semilogy(abs(Residual),'LineWidth',3);
xlabel('Number of iteration');
ylabel('Residual');
title('Wall boundary condition, p=1');
ax = gca;
ax.FontSize = 20;

filename = 'wall_test_p2_Res.dat';
f = fopen(filename);
Residual = fscanf(f,'%f');
fclose(f);

figure;
semilogy(abs(Residual),'LineWidth',3);
xlabel('Number of iteration');
ylabel('Residual');
title('Wall boundary condition, p=2');
ax = gca;
ax.FontSize = 20;