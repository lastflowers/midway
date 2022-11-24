% The MIT License
% 
% Copyright (c) 2022 Jaehyung Jung
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.


clc; clear; close all;

%% Set path and load necessary files
path_dataset = 'path_to_dataset';
path_toolbox = 'path_to_ycb_video_toolbox';
path_pose = 'cosypose\0022\';
addpath(path_toolbox);
addpath('functions');

mat_file = dir([path_dataset, '*.mat']);
img_file = dir([path_dataset, '*-color.png']);
obj_model{1,1} = load([path_toolbox, 'models\005_tomato_soup_can.mat']);
obj_model{1,2} = load([path_toolbox, 'models\007_tuna_fish_can.mat']);
obj_model{1,3} = load([path_toolbox, 'models\009_gelatin_box.mat']);
obj_model{1,4} = load([path_toolbox, 'models\025_mug.mat']);

%% Load CosyPose estimates
% Exclude brick(21): unstable estimates
meta_0 = load(fullfile(mat_file(1,1).folder, mat_file(1,1).name));
meta_0.cls_indexes(5,:) = [];
meta_0.poses(:,:,5) = [];

% parse pose measurements per an object
No = size(meta_0.cls_indexes,1);
npy_file = cell(1,No);
npy_idx = cell(1,No);
for j = 1:No
    idx = meta_0.cls_indexes(j,1);
    if idx < 10
        npy_name = ['obj_00000', num2str(idx), '*'];
    else
        npy_name = ['obj_0000', num2str(idx), '*'];
    end
    npy_file{1,j} = dir([path_pose, npy_name]);
    N_j = size(npy_file{1,j},1);
    time_idx = zeros(1,N_j);
    for jj = 1:N_j
        time_idx(1,jj) = str2double(npy_file{1,j}(jj,1).name(end-9:end-4));
    end
    tmp.index = idx;
    tmp.time_index = time_idx;
    npy_idx{1,j} = tmp;
end


%% Set constants
opt = globals();  % global variables from the YCB-V toolbox
isPlot = 0;
d2r = pi/180;
r2d = 1/d2r;
N = size(mat_file,1);
w = 0.7;
a = 1;
b = [w, 1-w];
draw_color = [0, 128, 0; 128, 0, 0; 0, 0, 128; 160, 160, 0; 0, 128, 128]/255;

% estimator dimensions
Nc = 6;
Nd = Nc + No*Nc;
qid = 1:3;
pid = 4:6;
oqid = 1:3;
opid = 4:6;


%% Filter initialization
cov = eye(Nd);
cov(pid, pid) = (0.01)^2 * eye(3);
cov(qid ,qid) = (0.1*d2r)^2 * eye(3);
Q = eye(6);
Q(pid,pid) = (8.66e-3)^2 * eye(3);
Q(qid,qid) = (0.5196*d2r)^2* eye(3);

std_psi(1,1) = d2r*4;
std_psi(1,2) = d2r*8;
R = zeros(6,6,2);
R(qid,qid,1) = (std_psi(1,1))^2*eye(3);
R(pid,pid,1) = (5e-3)^2 * eye(3);
R(:,:,2) = R(:,:,1);
R(qid(3),qid(3),2) = (std_psi(1,2))^2;

Tc0_g = [meta_0.rotation_translation_matrix; 0, 0, 0, 1];
Tg_ci = inv(Tc0_g);
Tg_o = zeros(4,4,1,No);  % row, col, hypothesis, object
for j = 1:No
    Tc0_oj = double(readNPY(fullfile(npy_file{1,j}(1,1).folder, npy_file{1,j}(1,1).name)));
    Tg_o(:,:,1,j) = Tg_ci * Tc0_oj;
    cov(6*j+1:6*j+3, 6*j+1:6*j+3) = (4*d2r)^2 * eye(3);
    cov(6*j+4:6*j+6, 6*j+4:6*j+6) = (0.1)^2 * eye(3);
end


%% allocate memory
est_p = zeros(3,N);
est_eul = zeros(3,N);
est_po = zeros(3,N,No);
est_eulo = zeros(3,N,No);
est_rel_p = zeros(3,N,No);
est_rel_eul = zeros(3,N,No);
est_p(:,1) = Tg_ci(1:3,4);
est_eul(:,1) = dcm2euler(Tg_ci(1:3,1:3));

err_p = zeros(3,N);
err_q = zeros(3,N);
err_po = zeros(3,N,No);
err_qo = zeros(3,N,No);
err_rel_p = zeros(3,N,No);
err_rel_q = zeros(3,N,No);
time_elapsed = zeros(1,N);

std_p = zeros(3,N);
std_q = zeros(3,N);
std_po = zeros(3,N,No);
std_qo = zeros(3,N,No);
std_p(:,1) = sqrt(diag(cov(pid,pid)));
std_q(:,1) = sqrt(diag(cov(qid,qid)));

true_p = zeros(3,N);
true_eul = zeros(3,N);
true_po = zeros(3,N,No);
true_eulo = zeros(3,N,No);
true_p(:,1) = Tg_ci(1:3,4);
true_eul(:,1) = dcm2euler(Tg_ci(1:3,1:3));

for j = 1:No
    true_Tco = [meta_0.poses(:,:,j); 0, 0, 0, 1];
    true_Tgo = Tg_ci*true_Tco;
    oqid_j = 6*j+1:6*j+3;
    opid_j = 6*j+4:6*j+6;

    est_po(:,1,j) = Tg_o(1:3,4,1,j);
    est_eulo(:,1,j) = dcm2euler(Tg_o(1:3,1:3,1,j));
    true_po(:,1,j) = true_Tgo(1:3,4);
    true_eulo(:,1,j) = dcm2euler(true_Tgo(1:3,1:3));

    std_qo(:,1,j) = sqrt(diag(cov(oqid_j,oqid_j)));
    std_po(:,1,j) = sqrt(diag(cov(opid_j,opid_j)));

    err_Tgo = -Log_se3(Tg_o(:,:,1,j) / true_Tgo);
    err_qo(:,1,j) = err_Tgo(1:3,1);
    err_po(:,1,j) = err_Tgo(4:6,1);

    est_Tco = Tg_ci\Tg_o(:,:,1,j);
    err_rel = Log_se3(true_Tco\est_Tco);
    err_rel_q(:,1,j) = err_rel(1:3,1);
    err_rel_p(:,1,j) = err_rel(4:6,1);
end


%% Main loop
if isPlot
    % set figure
    h0 = figure(1);
    s0 = subplot(1,2,1);
    s1 = subplot(1,2,2);
end
for i = 2:N
    fprintf('[%d / %d] image ...\n', i, N);
    meta_i = load(fullfile(mat_file(i,1).folder, mat_file(i,1).name));

    % read pose measurements
    img_i = imread(fullfile(img_file(i,1).folder, img_file(i,1).name));
    Y_i = zeros(4,4,No);
    meas_flag = false(1,No);
    for j = 1:No
        idx_i = find(i == npy_idx{1,j}.time_index);
        if ~isempty(idx_i)
           Y_i(:,:,j) = double(readNPY(fullfile(npy_file{1,j}(idx_i,1).folder, npy_file{1,j}(idx_i,1).name)));
           meas_flag(1,j) = true;             
        end
    end

    % Predicton for every hypothesis
    Phi = eye(Nd);
    for n = 1:size(a,2)
        G = Adj_se3(Tg_ci(:,:,n));
        cov(:,:,n) = Phi*cov(:,:,n)*Phi';
        cov(1:6,1:6,n) = cov(1:6,1:6,n) + G*Q*G';
    end

    % filter update for every hypothesis
    Na = size(a,2);
    Nb = size(b,2);
    Tgci_upd = zeros(4,4,Na*Nb);
    Tgo_upd = zeros(4,4,Na*Nb,No);
    cov_upd = zeros(Nd,Nd,Na*Nb);
    ua = zeros(1,Na*Nb);
    cnt = 1;
    for n = 1:Na % per prior hypothesis
        Tgci_n = Tg_ci(:,:,n);
        Tinv = inv(Tgci_n);
        Ad_Tinv = Adj_se3(Tinv);
        cov_n = cov(:,:,n);

        for m = 1:Nb % per likelihood hypothesis
            N_vis = sum(meas_flag);
            H = zeros(6*N_vis,Nd);
            R_m = zeros(6*N_vis);
            innov = zeros(6*N_vis,1);
            cnt_j = 1;
            
            for j = 1:No % per object
                if meas_flag(1,j)  % if measured
                    Tgo_nj = Tg_o(:,:,n,j);
                    Ad_Yg = Adj_se3(Tgo_nj);
                    ridx = 6*(cnt_j-1)+1:6*cnt_j;
                    cidx = 6*j+1:6*(j+1);
                    H(ridx,1:6) = -Ad_Tinv;
                    H(ridx,cidx) = Ad_Tinv;

                    dY = (Tgci_n\Tgo_nj)/Y_i(:,:,j);
                    innov(ridx,1) = Log_se3(dY);

                    R_m(ridx,ridx) = Ad_Tinv*Ad_Yg*R(:,:,m)*Ad_Yg'*Ad_Tinv';
                    cnt_j = cnt_j + 1;
                end
            end

            S_nm = H*cov_n*H' + R_m;
            K = cov_n*H'/S_nm;
            I_KH = eye(Nd,Nd) - K*H;
            cov_upd(:,:,cnt) = I_KH*cov_n*I_KH' + K*R_m*K';
            delx = K*innov;

            Tgci_upd(:,:,cnt) = Exp_se3(-delx(1:6))*Tgci_n;
            for j = 1:No
                ridx = 6*j+1:6*j+6;
                Tgo_upd(:,:,cnt,j) = Exp_se3(-delx(ridx,1))*Tg_o(:,:,n,j);
            end

            b_m = b(1,m);
            gam_nm = (2*pi)^(-3)*det(S_nm)^(-0.5) * exp(-0.5*innov'*(S_nm\innov));
            if gam_nm == 0
                gam_nm = 1e-10;
            end
            ua(1,cnt) = a(1,n) * b_m * gam_nm;

            cnt = cnt +1;
        end
    end
    Tg_ci = Tgci_upd;
    Tg_o = Tgo_upd;
    cov = cov_upd;
    a = ua/sum(ua);

    % prune small weights
    prun_idx = a < 5e-6;
    a(:,prun_idx) = [];
    Tg_ci(:,:,prun_idx) = [];
    Tg_o(:,:,prun_idx,:) = [];
    cov(:,:,prun_idx) = [];
    a = a/sum(a);

    % Midway-merge
    tmp_a = a;
    tmp_Tgci = Tg_ci;
    tmp_Tgo = Tg_o;
    tmp_cov = cov; 
    if size(a,2) >= 16
        while size(tmp_a,2) > 12
            % Merge the largest and smallest weights
            [~, i_idx] = max(tmp_a);
            [~, j_idx] = min(tmp_a);  

            uw_i = tmp_a(1,i_idx);
            uw_j = tmp_a(1,j_idx);
            Tc_i = tmp_Tgci(:,:,i_idx);
            Tc_j = tmp_Tgci(:,:,j_idx);
            To_i = squeeze(tmp_Tgo(:,:,i_idx,:));
            To_j = squeeze(tmp_Tgo(:,:,j_idx,:));
            cov_i = tmp_cov(:,:,i_idx);
            cov_j = tmp_cov(:,:,j_idx);

            [w_ij, Tc_ij, To_ij, cov_ij] = midwayMerge(uw_i, Tc_i, To_i, cov_i, uw_j, Tc_j, To_j, cov_j);

            tmp_a(:,[i_idx, j_idx]) = [];
            tmp_Tgci(:,:,[i_idx, j_idx]) = [];
            tmp_Tgo(:,:,[i_idx, j_idx],:) = [];
            tmp_cov(:,:,[i_idx, j_idx]) = [];

            N_h = size(tmp_Tgci,3);
            tmp_a(:,N_h+1) = w_ij;
            tmp_Tgci(:,:,N_h+1) = Tc_ij;
            for j = 1:No
                tmp_Tgo(:,:,N_h+1,j) = To_ij(:,:,j);
            end
            tmp_cov(:,:,N_h+1) = cov_ij;
        end
    end
    a = tmp_a;
    Tg_ci = tmp_Tgci;
    Tg_o = tmp_Tgo;
    cov = tmp_cov;

    % Summarize estimates
    if size(a,2) > 1
        [~,idx0] = max(a);
        tmp_a = a;
        tmp_a(:,idx0) = [];
        [smax_val,~] = max(tmp_a);
        idx1 = find(a == smax_val);
        idx1 = idx1(1,1);
        [~, Tc_hat, To_hat, cov_hat] = midwayMerge(a(1,idx0), Tg_ci(:,:,idx0), squeeze(Tg_o(:,:,idx0,:)), cov(:,:,idx0), ...
            a(1,idx1), Tg_ci(:,:,idx1), squeeze(Tg_o(:,:,idx1,:)), cov(:,:,idx1));
    else
        Tc_hat = Tg_ci(:,:,1);
        To_hat = Tg_o(:,:,1,:);
        cov_hat = cov(:,:,1);
    end
    
    % Draw
    if isPlot
        set(0, 'CurrentFigure', h0); cla;
        set(gcf, 'CurrentAxes', s0); cla;
        imshow(img_i); hold on;
        set(gcf, 'CurrentAxes', s1); cla;
        imshow(img_i); hold on;   
        for j = 1:No
            color_j = draw_color(j,:);          
            obj_j = obj_model{1,j};
            est_Tcoj = Tc_hat\To_hat(:,:,j);
            x2d_meas = project(obj_j.obj.v', meta_i.intrinsic_matrix, Y_i(1:3,:,j));
            x2d = project(obj_j.obj.v', meta_i.intrinsic_matrix, est_Tcoj(1:3,:));
            
            set(0, 'CurrentFigure', h0);
            set(gcf, 'CurrentAxes', s0);
            patch('vertices', x2d_meas, 'faces', obj_j.obj.f3', ...
                'FaceColor', color_j, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
            set(gcf,'color','w');
            set(gca,'color','w');
            title('cosypose single view');
            
            set(gcf, 'CurrentAxes', s1);
            patch('vertices', x2d, 'faces', obj_j.obj.f3', ...
                'FaceColor', color_j, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
            set(gcf,'color','w');
            set(gca,'color','w');
            title('GM-IEKF with midway-merge');
        end
        pause(0.005);
    end

    % save estimates
    est_p(:,i) = Tc_hat(1:3,4);
    est_eul(:,i) = dcm2euler(Tc_hat(1:3,1:3));
    std_p(:,i) = sqrt(diag(cov_hat(pid,pid)));
    std_q(:,i) = sqrt(diag(cov_hat(qid,qid)));

    % save true (only for reference)
    true_Tcg = [meta_i.rotation_translation_matrix; 0, 0, 0, 1];
    true_Tgc = inv(true_Tcg);
    true_p(:,i) = true_Tgc(1:3,4);
    true_eul(:,i) = dcm2euler(true_Tgc(1:3,1:3));

    % save error
    err_Tgc = -Log_se3(Tc_hat / true_Tgc);
    err_q(:,i) = err_Tgc(1:3,1);
    err_p(:,i) = err_Tgc(4:6,1);
    for j = 1:No
        ridx = 6*j+1:6*j+6;

        est_po(:,i,j) = To_hat(1:3,4,j);
        est_eulo(:,i,j) = dcm2euler(To_hat(1:3,1:3,j));
        std_qo(:,i,j) = sqrt(diag(cov_hat(ridx(1:3),ridx(1:3))));
        std_po(:,i,j) = sqrt(diag(cov_hat(ridx(4:6),ridx(4:6))));

        match_j = find(meta_i.cls_indexes == npy_idx{1,j}.index);
        true_Tco = [meta_i.poses(:,:,match_j); 0, 0, 0, 1];
        true_Tgo = true_Tcg\true_Tco;
        true_po(:,i,j) = true_Tgo(1:3,4);
        true_eulo(:,i,j) = dcm2euler(true_Tgo(1:3,1:3));

        err_Tgo = -Log_se3(To_hat(:,:,j) / true_Tgo);
        err_qo(:,i,j) = err_Tgo(1:3,1);
        err_po(:,i,j) = err_Tgo(4:6,1);

        est_Tco = Tc_hat\To_hat(:,:,j);
        err_rel = Log_se3(true_Tco\est_Tco);
        err_rel_q(:,i,j) = err_rel(1:3,1);
        err_rel_p(:,i,j) = err_rel(4:6,1);
    end
end

%% Report estimation error
rmse_q = sqrt(sum(sqrt(sum(err_q.^2,1)).^2/N));
rmse_p = sqrt(sum(sqrt(sum(err_p.^2,1)).^2/N));
rmse_qo = sqrt(sum(sqrt(sum(err_qo.^2,1)).^2/N));
rmse_po = sqrt(sum(sqrt(sum(err_po.^2,1)).^2/N));
rmse_rel_q = sqrt(sum(sqrt(sum(err_rel_q.^2,1)).^2/N));
rmse_rel_p = sqrt(sum(sqrt(sum(err_rel_p.^2,1)).^2/N));

fprintf('Robot RMSE: %.2f [cm]  %.2f [deg] \n', 100*rmse_p, r2d*rmse_q);
for j = 1:No
    rmse_po_j = mean(sqrt(sum(rmse_po(:,:,j).^2,1)));
    rmse_qo_j = mean(sqrt(sum(rmse_qo(:,:,j).^2,1)));
    rmse_relp_j = mean(sqrt(sum(rmse_rel_p(:,:,j).^2,1)));
    rmse_relq_j = mean(sqrt(sum(rmse_rel_q(:,:,j).^2,1)));

    fprintf('%d Object RMSE: %.2f [cm]  %.2f [deg]  Rel RMSE: %.2f [cm]  %.2f [deg]\n', ...
        meta_0.cls_indexes(j,1), 100*rmse_po_j, r2d*rmse_qo_j, 100*rmse_relp_j, r2d*rmse_relq_j);
end 


%% Plot: ground-truth vs. estimates
figure;
for i = 1:3
    subplot(3,2,2*i-1);
    plot(r2d*true_eul(i,:), 'k'); hold on;
    plot(r2d*est_eul(i,:), 'r'); grid on;
    if i == 1
        ylabel('\phi [deg]');
        title('Robot euler angles');
        legend('ground-truth', 'gm-iekf (midway)');
    elseif i == 2
        ylabel('\theta [deg]');
    elseif i == 3
        ylabel('\psi [deg]');
        xlabel('image index');
    end
end
for i = 1:3
    subplot(3,2,2*i);
    plot(true_p(i,:), 'k'); hold on;
    plot(est_p(i,:), 'r'); grid on;
    if i == 1
        title('Robot position');
        ylabel('p^g_x [m]');
    elseif i == 2
        ylabel('p^g_y [m]');
    elseif i == 3
        ylabel('p^g_z [m]');
        xlabel('image index');
    end
end

for j = 1:No
    obj_idx = meta_i.cls_indexes(j,1);
    figure;
    for i = 1:3
        subplot(3,2,2*i-1);
        plot(r2d*true_eulo(i,:,j), 'k'); hold on;
        plot(r2d*est_eulo(i,:,j), 'r'); grid on;
        if i == 1
            title([num2str(obj_idx), ' object euler angles']);
            ylabel('\phi_o [deg]');
            legend('ground-truth', 'gm-iekf (midway)');
        elseif i == 2
            ylabel('\theta_o [deg]');
        elseif i == 3
            ylabel('\psi_o [deg]');
            xlabel('image index');
        end
    end
    for i = 1:3
        subplot(3,2,2*i);
        plot(true_po(i,:,j), 'k'); hold on;
        plot(est_po(i,:,j), 'r'); grid on;
        if i == 1
            title([num2str(obj_idx), ' object position']);
            ylabel('p^g_{ox} [m]');
        elseif i == 2
            ylabel('p^g_{oy} [m]');
        elseif i == 3
            ylabel('p^g_{oz} [m]');
            xlabel('image index');
        end
    end
end