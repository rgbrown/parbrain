close all;
clear all;
clc;

[S, fn, fc] = data_reader();
[t, Y] = fn();
Yall{1} = Y;
T(1) = t;
k = 2;
while ~isempty(Y)
    [t, Y] = fn();
   
    Yall{k} = Y;
    if ~isempty(t)
        T(k) = t;
    end
    k = k + 1;
end
Yall(end) = [];
fc()
%T(end) = [];
%% Get out the data (radius etc.) from one block  - which one? x{idx}(1) = 1st
getdata = @(idx) cellfun(@(x) x{idx}(1), Yall);
i_radius = getdata(1);
R_k      = getdata(2);
N_Na_k   = getdata(3);
N_K_k    = getdata(4);
N_HCO3_k = getdata(5);
N_Cl_k   = getdata(6);
N_Na_s   = getdata(7);
N_K_s    = getdata(8);
N_HCO3_s = getdata(9);
K_p      = getdata(10);                              
w_k      = getdata(11);
      
      
ca_i     = getdata(12);
ca_sr_i  = getdata(13);
v_i      = getdata(14);
w_i      = getdata(15);
ip3_i    = getdata(16);

K_i      = getdata(17);
ca_j     = getdata(18);
ca_er_j  = getdata(19);
v_j      = getdata(20);
ip3_j    = getdata(21);
Mp       = getdata(22);
AMp      = getdata(23);
AM       = getdata(24);

figure
subplot(5,2,1)
plot(T, ca_i)
ylabel('Ca_i')
title('SMC')

subplot(5,2,2)
plot(T, ca_j)
ylabel('Ca_i')
title('EC')

subplot(5,2,3)
plot(T, ip3_i)
ylabel('IP3_i')

subplot(5,2,4)
plot(T, ip3_j)
ylabel('IP3_j')

subplot(5,2,5)
plot(T, ca_sr_i)
ylabel('s_i')

subplot(5,2,6)
plot(T, ca_er_j)
ylabel('s_j')

subplot(5,2,7)
plot(T, v_i)
ylabel('v_i')

subplot(5,2,8)
plot(T, v_j)
ylabel('v_j')

subplot(5,2,9)
plot(T, w_i)
ylabel('w_i')

subplot(5,2,10)
plot(T, i_radius*1e6)
ylabel('Radius (micro m)')
% hier muss man ranzoomen!




% [S, fn, fc] = data_reader();
% [t, Y] = fn();
% idx = 4;
% %%
% h = imagesc(Y{idx}, [0 10]);
% 
% 
% %%
% while 1
%     tprev = t;
%     [t, Y] = fn();
%     if isempty(Y)
%         break
%     end
%     set(h, 'CData', Y{idx});
%     pause(t - tprev);
% end
% fc();