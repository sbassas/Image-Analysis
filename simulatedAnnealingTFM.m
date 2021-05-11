%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated Annealing TFM, May 10 2021, Séraphin Bassas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE: The purpose of this script is to converge on the optimal 
% traction force parameter that will, when it is fed into the Peñas et al. 
% script TFM_model.m, reach a target average displacement value obtained in 
% experimental TFM. By comparing the displacement values obtained from
% running TFM on images obtained experimentally and from the TFM_model.m
% simulation, one can validate the results obtained from experimental TFM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREDIT: The TFM_model.m script can be found in the following paper: A. 
% Jorge-Peñas, A. Muñoz-Barrutia, E.M. de-Juan-Pardo & C. Ortiz-
% de-Solorzano (2014): Validation tool for traction force microscopy, 
% Computer Methods in Biomechanics and Biomedical Engineering,
% DOI: 10.1080/10255842.2014.903934 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE: This script uses Simulated Annealing to optimize the output 
% obtained from running the TFM_model.m simulation for different traction 
% force input parameter settings, which is represented by the variable 
% called m_force in the TFM_model.m script. In this algorithm, each 
% completed simulation, which is comprised of a set of one force input 
% parameter, several output images, and one ouput displacement matrix, 
% represents a candidate solution to the optimization problem. But, to 
% generate a candidate, one must input a force into the TFM_simulator 
% function in this script, adapated from the TFM_model.m script, which 
% yields simulated images and computes displacements. Further, neighbours 
% of a candidate solution are chosen deterministically (± 1e-5 or ±1e-6) 
% and stochastically (± rand*1e-5). And finally, the evaluation function 
% that gives a score to each candidate solution is simply the product of 
% the difference of zero value pixels between the null force image and the 
% distorted image and the absolute value of the difference between the 
% target displacement and the computed displacement of a candidate 
% solution. The output of this script, the optimized force parameter, is
% 1x4 variable with the force parameter at (1,1), the force application
% radius at (1,2), the black pixel count in distorted image at (1,3), and
% the averaged displacement at (1,4).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISCLAIMER: This script optimizes the force input parameters for the
% script from the article cited above. This script assumes that images
% obtained using the TFM_model.m script best reproduce fluorescent bead 
% matrix deformations the way the parameters (excluding force) are preset
% in the TFM_simulation function. However, for any external use, these 
% parameters (fr, cell_l, cell_w, fr_a, etc.) as they appear in the 
% TFM_simulation function, should be tuned to the performance of your own
% experimental setup for best results. This script also assumes that the
% values that the force parameter can take are bounded above and below, as
% forces too weak produce unpercievable deformation, and forces to great
% distort beads in the simulated images beyond reasonable thresholds.

%% clear
clc
clear all

%% user defined parameters

base_im = imread('after.tiff'); %null displacement image
m_force_lower = 1e-5; %lower bound on m_force parameter
m_force_higher = 1e-4; %upper bound on m_force parameter
fr= [8]*1e-6; %List of force application gaussian area radius in micron (focal adhesion)

iterations = 9; % #iterations before termination, this was the number I typically found at which the algorithm converged on a solution
disp_target = 6.5146e-07; % taget displacement value to converge on

%counting the number of black pixels (value 0) in null force image
init_black_px = size(find(~base_im),1);

%% converging
latestCandidateSolution = TFM_simulator(base_im, (m_force_lower+ m_force_higher)/2, fr, init_black_px,0); % where the candidate solution from the last iteration is stored
lastScore = computeEvaluationFunctionScore(latestCandidateSolution{1,3}, latestCandidateSolution{1,4}, disp_target);
scoreLastIteration = latestCandidateSolution;
tempScore = 0; 
bestScore = lastScore; % where the best overall score is stored
bestSolution = latestCandidateSolution; % where the best overall candidate solution is stored
i=1; % current iteration of the algorithm, incremented after each iteration
probability = 0.8; % initial probability with which a "bad solution" will be chosen to be your candidate solution
% temperature = 0.8; % temperature parameter for a more complex probability variable 

% while loop may be chosen to terminate after a certain amount of
% iterations or whenever the latestCandidateSolution variable stops 
% changing between iterations
while i < iterations 
    
    % storing all the neighbouring force input parameter values computed
    % using neighbourhood function
    neighbourhood = findNeighbours(m_force_lower, m_force_higher, latestCandidateSolution{1,1});
    % storing score of the candidate solution under consideration
    scoreLastIteration = computeEvaluationFunctionScore(latestCandidateSolution{1,3}, latestCandidateSolution{1,4}, disp_target);
    
    %for each neigbour check if its evaluation score is better
    for j = 1:size(neighbourhood,2)
        
        %if the force is NaN then skip to next neighbour
        if ~isnan( neighbourhood(1,j) )
            
            %perform simulation with new force parameter
            tempRun = TFM_simulator( base_im, neighbourhood(1,j), fr, init_black_px,i );
            tempScore = computeEvaluationFunctionScore(tempRun{1,3}, tempRun{1,4}, disp_target);
            % updating the probability of choosing a "bad solution" according to the Boltzmann distribution
            % this bugged the program though
%             probability = exp( - (scoreLastIteration - tempScore) / temperature);
            
            % checking if this neighbouring candidate scores better than
            % the best overall score
            if tempScore < bestScore
                
                bestScore = tempScore;
                bestSolution = tempRun;
            end
            
            % checking if this neighbouring candidate scores better than
            % the candidate solution under consideration at this iteration
            % without breaking out of the loop, to consider all neighbours
            if tempScore < lastScore
                
                latestCandidateSolution = tempRun;
                lastScore = tempScore;
            
            % checking if this neighbouring candidate solution should be 
            % accepted as a "bad solution" with a certain probability in 
            % which case we break out of the while loop
            else
                x=rand;
                if x < probability
                    latestCandidateSolution = tempRun;
                    lastScore = tempScore;
                    break
                end
            end
            
        end
    end
    %updating the probability of choosing a bad solution
    probability = probability*0.6;
    i=i+1;
end

% At this point the algorithm is finished and the best candidate solution
% can be found stored in the variable 'bestSolution' or saved using the
% following set of commands.

% newdir=['Simulated annealing/Target_disp',bestSolution];
% mkdir(newdir)
% save([newdir,'/OptimizedForceValue.mat'],'bestSolution','-v7.3')            

disp('Simulation finished')


% **************************** FUNCTIONS BELOW ****************************
% *************************************************************************
%% evaluation function of a simulation
function score = computeEvaluationFunctionScore(pixels, proximity, target)
score = pixels* pixels * abs(target - proximity);
end


%% function to determine neighbors of a TFM_simulation candidate solution
function returnVal = findNeighbours(lower_force, upper_force, m_force)

neighbours = zeros(1,6);
neighbours(1,1) = m_force - 1e-5;
neighbours(1,2) = m_force - 1e-6;
neighbours(1,3) = m_force + 1e-6;
neighbours(1,4) = m_force + 1e-5;
neighbours(1,5) = m_force - rand*1e-5;
neighbours(1,6) = m_force + rand*1e-5;
% for each force value in the 'neighbour' variable, check if it fits within
% the upper and lower bounds of permissible force values, otherwise discard
% it
for i=1:size(neighbours,2)
    if  neighbours(1,i) < lower_force || neighbours(1,i) > upper_force
        neighbours(1,i) = NaN;
    end
end
returnVal = neighbours;
end

%% converted 3D TFM simulator to function so that it takes input parameters
%%changed the file location where the images are saved and removed some of
%%the data from original TFM_model.m script
function returnData=TFM_simulator(base_im_frame, m_force, fr, init_black_px,iter)
str = sprintf('%d',iter);
%===========================================
%Simulation parameters
%===========================================

%Virtual cell parameters
%Cell dimension and position parameters
cell_l = 75e-6;%cell length
cell_w = 18e-6;%cell width
%Cell main axis orientation (*****LOOPED OVER LIST*****)
cell_angle = [40];

%Simulated force time profile (*****CHOOSE ONE OR SEVERAL IN CELL ARRAY*****)
F_shape = {'Triangle'};%Choose one of: 'Triangle';'Square';'Gaussian';'Custom';
n_peaks = 1;% number of peaks
peak_freq = 1;%frequency of peaks in Hz
peak_dr = 0.5;%Duty ratio of force: peak force/contraction cycle duration
F_ar = 1;%Force application guassian radius dimension aspect ratio

%Material parameters
E = 10e4;%Yongs modulus
nu = 0.49;%Poisson ratio

%Video parameters
dimx =size(base_im_frame,2);%Image dimension in pixels (*****MUST MATCH THAT OF THE MODEL IMAGE*****)
dimy = size(base_im_frame,1);
fps = 3; %frame per seconds
mu_pix = 0.1894;%micron per pixel

%===========================================

%parameters conversion
A=(1+nu)/(pi*E);
D = mu_pix*1e-6;%conversion factor
Apx=D^2;%px area in micron
cell_x = dimx/2*D;%
cell_y = dimy/2*D;
dur_vid = n_peaks/peak_freq;% duration in seconds
dur_peak = 1/peak_freq;% duration in seconds
n_frame = dur_vid*fps;%number of frame per video
n_frame_peak = floor(fps*dur_peak);%number of frames for one contraction cycle
end_peaks = n_frame_peak;
for n = n_peaks+1
    end_peaks = [end_peaks, n_frame_peak*n_peaks];
end
n_frame_n0peak = floor(fps*dur_peak*peak_dr);%number of frames during non-zero active contraction
time_vec = linspace(0,dur_vid,n_frame);%time vector
time_vec_peak = linspace(0,dur_peak,n_frame_peak);%time vector for one contraction cycle
time_vec_n0peak = linspace(0,dur_peak*peak_dr,n_frame_n0peak);%time vector

% initialize some variables for speed
Trx = zeros(dimy,dimx, n_frame);
Try = zeros(dimy,dimx,n_frame);
fu = zeros(dimy,dimx,n_frame);
fv = zeros(dimy,dimx,n_frame);

F_tot = zeros(1,n_frame);
Fx_tot = zeros(1,n_frame);
Fy_tot = zeros(1,n_frame);

%Meashgrid
[X,Y] = meshgrid(1:dimx,1:dimy);
xdata(:,:,1) = X*D;
xdata(:,:,2) = Y*D;

%===========================================
%loop over Force time profiles
for m = 1:size(F_shape,2)
    F_type = F_shape{m};
    
    %Loop over the list of force application radius
    for n = 1:size(fr,2)
        
        %Force application guassian radius dimension aspect ratio
        fr_x = fr(n)*F_ar;
        fr_y = fr(n);
        %Norm of the force application area radiuses
        norm = 1/(2*pi*(fr_x^2+fr_y^2));
        
        %%Square force up and down in time
        sq_peak = [zeros(1,floor((n_frame_peak-n_frame_n0peak)/2)),...
            m_force*ones(1,n_frame_n0peak),...
            zeros(1,n_frame_peak-floor((n_frame_peak-n_frame_n0peak)/2)-n_frame_n0peak)];
        %%Triangle force up-down
        tri_peak = [zeros(1,floor((n_frame_peak-n_frame_n0peak)/2)),...
            (n_frame_n0peak/fps-abs(linspace(0,n_frame_n0peak/fps,n_frame_n0peak)-linspace(n_frame_n0peak/fps,0,n_frame_n0peak)))*m_force/(n_frame_n0peak/fps),...
            zeros(1,floor((n_frame_peak-n_frame_n0peak)/2))];
        %%Gaussian force peak in time
        gau_peak = m_force*exp(-((time_vec_peak-time_vec_n0peak(floor(n_frame_peak/2)))/(peak_dr/2)).^2);
        gau_peak(1) = 0; gau_peak(end) = 0;
        f_sh_vecs = {sq_peak,tri_peak,gau_peak};
        
        switch F_type
            
            case 'Gaussian'
                force_vec = f_sh_vecs{1};
                for i = 1:n_peaks-1
                    force_vec = [force_vec, f_sh_vecs{1}];
                end
            case 'Square'
                force_vec = f_sh_vecs{2};
                for i = 1:n_peaks-1
                    force_vec = [force_vec, f_sh_vecs{2}];
                end
            case 'Triangle'
                force_vec = f_sh_vecs{3};
                for i = 1:n_peaks-1
                    force_vec = [force_vec, f_sh_vecs{3}];
                end
            case 'Custom'
                force_vec = f_sh_vecs{1};
                for i = 1:n_peaks-1
                    force_vec = [force_vec, f_sh_vecs{i+1}];
                end
        end
        
        
        %         figure(1)
        %         plot(time_vec, force_vec) %force time profile
        %plot(linspace(0,peak_dr,length(force)), force) %force time profile
        
        %mask parameter
        phi = linspace(0,2*pi,180);
        cosphi = cos(phi);
        sinphi = sin(phi);
        
        %prepare movie figures
        
        % Plot displacement
        zlim = m_force;
        %         fig1 = figure(1);%'Name','Displacement magnitude');
        %         surf(xdata(:,:,1),xdata(:,:,2),sqrt(fu(:,:,1).^2+fv(:,:,1).^2),sqrt(fu(:,:,1).^2+fv(:,:,1).^2),'EdgeColor','none')
        %hold on
        %quiver3(xdata(:,:,1),xdata(:,:,2),-ones(dim,dim)*zlim,u,v,zeros(dim,dim),2,'k');
        %         axis tight manual
        %         ax = fig1.CurrentAxes;
        %         ax.XLim = [0 dimx*D];
        %         ax.YLim = [0 dimy*D];
        %         ax.ZLim = [-zlim zlim];
        %         ax.CLim = [0 zlim];
        %         ax.NextPlot = 'replaceChildren';
        
        %Plot tractions
        zlim2 = 2.5e2;
        %         fig2 = figure(2);%'Name','Traction magnitude');
        %         surf(xdata(:,:,1),xdata(:,:,2),Trx(:,:,1),'EdgeColor','none')
        %hold on
        %quiver3(xdata(:,:,1),xdata(:,:,2),-ones(dim,dim)*zlim,u,v,zeros(dim,dim),2,'k');
        %         axis tight manual
        %         ax2 = fig2.CurrentAxes;
        %         ax2.XLim = [0 dimx*D];
        %         ax2.YLim = [0 dimy*D];
        %         ax2.ZLim = [-zlim2 zlim2];
        %         ax2.CLim = [0 zlim2];
        %         ax2.NextPlot = 'replaceChildren';
        
        
        %loop over cell orientation angle
        for k = 1:size(cell_angle,2)
            
            
            % initialize some variables for speed
            u_fft = zeros(dimy,dimx);
            v_fft = zeros(dimy,dimx);
            Trx_save = zeros(dimy,dimx,n_frame);
            Try_save = zeros(dimy,dimx,n_frame);
            fu_save = zeros(dimy,dimx,n_frame);
            fv_save = zeros(dimy,dimx,n_frame);
            x = zeros(dimy,dimx,n_frame);
            y = zeros(dimy,dimx,n_frame);
            
            Mov1 = struct('cdata',cell(1,n_frame),'colormap',cell(1,n_frame));
            Mov2 = struct('cdata',cell(1,n_frame),'colormap',cell(1,n_frame));
            U = zeros(1,n_frame);
            av_disp = zeros(1,n_frame);
            disp2_prev = zeros(dimy,dimx,2);
            prev_frame = base_im_frame;
            
            
            
            %loop over force magnitude
            for l = 1:length(force_vec)
                
                %cell orientation angle
                cell_theta = cell_angle(k);
                
                
                %calculate ellipse mask
                % Calculate cell ellipse
                R = [ cosd(cell_theta)   sind(cell_theta)
                    -sind(cell_theta)   cosd(cell_theta)];
                xy = [cell_l/2*cosphi; cell_w/2*sinphi];
                xy = R*xy;
                x_el = xy(1,:) + cell_x;
                y_el = xy(2,:) + cell_y;
                mask=poly2mask(x_el./D,y_el./D,dimy,dimx);
                
                
                %calculate blue ellipse
                %calculate based on predicted force propagation
                fu_pred=D;
                %                 d=m_force*(1+nu)/(fu_pred*E*pi)*1e-6/D;
                an=0.75*sqrt(3)*cell_l/2;
                bn=1/0.75*sqrt(3)*cell_w/2;
                xyn = [an*cosphi; bn*sinphi];
                xyn = R*xyn;
                xn = xyn(1,:) + cell_x;
                yn = xyn(2,:) + cell_y;
                %generate blue ellipse mask
                mask_ellipse=poly2mask(xn/D,yn/D,dimy,dimx);
                
                %                 figure('Name','Mask')
                %                 imagesc(mask)
                
                %force application position
                fx_p_1 = cell_x-cell_l/2*cosd(-cell_theta);
                fy_p_1 = cell_y-cell_l/2*sind(-cell_theta);
                fx_p_2 = cell_x+cell_l/2*cosd(-cell_theta);
                fy_p_2 = cell_y+cell_l/2*sind(-cell_theta);
                
                %Distribute force as a 2D Gaussian
                para = [norm,fx_p_1,fr_x,fy_p_1,fr_y,cell_theta];
                %Rotate Gaussian axis
                xdatarot(:,:,1)= xdata(:,:,1)*cos(para(6)) - xdata(:,:,2)*sin(para(6));
                xdatarot(:,:,2)= xdata(:,:,1)*sin(para(6)) + xdata(:,:,2)*cos(para(6));
                x0rot = para(2)*cos(para(6)) - para(4)*sin(para(6));
                y0rot = para(2)*sin(para(6)) + para(4)*cos(para(6));
                %Generate normalized Gaussian stress distribution 1st end of the cell
                T = para(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*para(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*para(5)^2) )    );
                
                %Repeat for 2nd end of the cell
                para = [-norm,fx_p_2,fr_x,fy_p_2,fr_y,cell_theta];
                
                xdatarot(:,:,1)= xdata(:,:,1)*cos(para(6)) - xdata(:,:,2)*sin(para(6));
                xdatarot(:,:,2)= xdata(:,:,1)*sin(para(6)) + xdata(:,:,2)*cos(para(6));
                x0rot = para(2)*cos(para(6)) - para(4)*sin(para(6));
                y0rot = para(2)*sin(para(6)) + para(4)*cos(para(6));
                
                %Sum two stress fields
                T = T+(  para(1)*exp(   -((xdatarot(:,:,1)-x0rot).^2/(2*para(3)^2) + (xdatarot(:,:,2)-y0rot).^2/(2*para(5)^2) )    ));
                
                %                 %Plot simulated stress field
                %                 figure('Name','Simulated stress field')
                %                 surf(T,'EdgeColor','none')
                
                %Project traction stress magnitude to x and y axis
                Trx = T*force_vec(l)*cosd(-cell_theta);
                Try = T*force_vec(l)*sind(-cell_theta);
                
                %Integrate force in the cellmask
                Fx = Apx*Trx.*mask_ellipse;
                Fy = Apx*Try.*mask_ellipse;
                Trt = sqrt( Trx.^2 + Try.^2 );
                Fr = Apx*Trt;%;.*mask;
                F_tot(l)=sum(sum(Fr));
                Fx_tot(l)=sum(sum(Fx));
                Fy_tot(l)=sum(sum(Fy));
                
                %Remove mean before fft2
                Tx1_0=(Trx)-nanmean(nanmean(Trx)).*ones(size(Trx,1),size(Trx,2));
                Ty1_0=(Try)-nanmean(nanmean(Try)).*ones(size(Try,1),size(Try,2));
                
                %Calculate the Fourier transform of the traction
                Tx_k = fft2(Tx1_0);
                Ty_k = fft2(Ty1_0);
                
                % Generate Fourier space coordinates vectors
                Nx=size(Tx_k,2);
                Ny=size(Ty_k,1);
                dkx = 1/(Nx*D);
                dky = 1/(Ny*D);
                kx = [0:fix(Nx/2)-1,-fix(Nx/2):-1]*dkx*2*pi;
                ky = [0:fix(Ny/2)-1,-fix(Ny/2):-1]*dky*2*pi;
                
                %Operate forward stress to displacement calculation in the Fourier
                %space
                ii=0;
                %loop over image pixels
                for i=ky(1:end)
                    ii=ii+1;
                    jj=0;
                    
                    for j=kx(1:end)
                        jj=jj+1;
                        
                        kn=sqrt(i^2+j^2);
                        Tnx=Tx_k(ii,jj);
                        Tny=Ty_k(ii,jj);
                        Tn=[Tnx;Tny];
                        
                        %                         r = sqrt(i^2+j^2);
                        
                        % Implement the Green tensor in the Fourier space
                        G=A*2*pi/(kn^3)*[(1-nu)*kn^2+nu*i^2,-nu*i*j; -nu*i*j,(1-nu)*kn^2+nu*j^2];
                        
                        %Take the product of the matrix in the Fourier domain
                        dn = G*Tn;
                        
                        u_fft(ii,jj)=dn(1);
                        v_fft(ii,jj)=dn(2);
                        
                        %Make sure no drift displacement is generated by nulling
                        %the zero-th order displacement in Fourier space
                        u_fft(1,1)=0;
                        v_fft(1,1)=0;
                        
                        
                    end
                end
                
                %Take the inverse Fourier transform of the displacement field
                fu(:,:,l) = real(ifft2(u_fft));
                fv(:,:,l) = real(ifft2(v_fft));
                % Apx=(conversion*1e-6)^2*(x1(1,2)-x1(1,1))*(y1(2,1)-y1(1,1));
                % Fx=Apx*Trx;
                
                
                
                %Plot x-displacement
                %                 fig1 = figure(fig1);
                %                 surf(ax,xdata(:,:,1),xdata(:,:,2),sqrt(fu(:,:,l).^2+fv(:,:,l).^2),sqrt(fu(:,:,l).^2+fv(:,:,l).^2),'EdgeColor','none')
                %hold on
                %quiver3(ax,xdata(:,:,1),xdata(:,:,2),-ones(dim,dim)*zlim,u,v,zeros(dim,dim),2);
                %                 drawnow
                %                 Mov1(l) = getframe(ax);
                %                 ax.NextPlot = 'replaceChildren';
                
                %Plot x-traction
                %                 fig2 = figure(fig2);
                %                 surf(ax2,xdata(:,:,1),xdata(:,:,2),sqrt(Trx.^2+Try.^2),sqrt(Trx.^2+Try.^2),'EdgeColor','none')
                %hold on
                %quiver3(ax,xdata(:,:,1),xdata(:,:,2),-ones(dim,dim)*zlim,u,v,zeros(dim,dim),2);
                %                 drawnow
                %                 Mov2(l) = getframe(ax2);
                %                 ax2.NextPlot = 'replaceChildren';
                
                %                 newdir=['Simulated annealing/',F_type,'_signal/fr',num2str(fr_x),'_iter_',str];
                %                 mkdir(newdir)
                
                %Warp beads image and write movie of deforming image
                % Warp base image
                disp2(:,:,1)=fu(:,:,l);
                disp2(:,:,2)=fv(:,:,l);
                beads_d = imwarp(base_im_frame,-disp2./D);
                %                 imwrite(beads_d,[newdir,'/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta),'.tif'],'WriteMode','append')
                
                
                %%Calculate Strain Energy
                Trx_vec = Trx(:);
                Try_vec = Try(:);
                fu_f = fu(:,:,l);
                fv_f = fv(:,:,l);
                u_vec = fu_f(:);
                v_vec = fv_f(:);
                dx = D;
                dy = D;
                U(l) = .5*sum((Trx_vec.*u_vec+Try_vec.*v_vec)*dx*dy);
                %%Calculate mean displacement
                av_disp(l) = mean(sqrt(fu_f(:).^2+fv_f(:).^2));
                averaged_avdisp = mean(av_disp);
                
                if l == 1+ floor(fps/2)
                    peak_black_px = size(find(~beads_d),1);
                end
                
            end
            
            %%Save data for Matlab GUI calculation
            %             newdir=['Simulated annealing/',F_type,'_signal','/fr',num2str(fr_x),'_iter_',str,'/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta)];
            %             mkdir(newdir)
            %             mkdir([newdir,'/Datasets/Average Displacement'])
            %             save([newdir,'/Datasets/Average Displacement','/sim_av_disp.mat'],'av_disp','-v7.3')
            %             save([newdir,'/Datasets/Average Displacement','/avg_avDisp.mat'],'averaged_avdisp','-v7.3')
            %             save([newdir,'/Datasets/Average Displacement','/evaluation_criteria.mat'],'eval_func_data','-v7.3')
            %             mkdir([newdir,'/Mask'])
            %             save([newdir,'/Mask/',F_type(1:3),'_',num2str(m_force),'N_rad_',num2str(fr_x),'m_phi_',num2str(cell_theta),'.mat'],'mask','-v7.3')
            
        end
    end
    
end
%Create return variable with data for computation of evaluation function
returnData = { m_force, fr, abs(init_black_px - peak_black_px), averaged_avdisp};
end