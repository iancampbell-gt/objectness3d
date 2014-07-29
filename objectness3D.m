function [Iout,whatScale,Voutx1,Vouty1,Voutz1,Voutx2,Vouty2,Voutz2,Voutx3,Vouty3,Voutz3,Lam1,Lam2,Lam3]=objectness3D(I,M,options)
%
% Modification of FrangiFilter3D from the Matlab file exchange, originally
% by D. Kroon. Computes objectness in 3D based upon the algorithm of Luca
% Antiga ("Generalizing vesselness with respect to dimensionality and
% shape", The Insight Journal, July 2007, http://hdl.handle.net/1926/576).
%
% INPUTS:
%       I: [IxJxK] matrix representing image stack to filter
%       M: integer dimensional order of shapes to segment (M > 3)
%           0 = blobs
%           1 = tubes
%           2 = plates
%           3 = hyper-plates (not supported in this function)
%       options : Struct with input options,
%           .ScaleRange : The range of sigmas used, default [1 8]
%           .ScaleRatio : Step size between sigmas, default 2
%           .Alpha : Frangi vesselness constant, threshold on Lambda2/Lambda3
%					   determines if its a line (vessel) or a plane like structure
%					   default .5;
%           .Beta  : Frangi vesselness constant, which determines the deviation
%					   from a blob like structure, default .5;
%           .Gamma     : Frangi vesselness constant which gives
%					   the threshold between eigenvalues of noise and 
%					   vessel structure. A rule of thumb is dividing the 
%					   the grey values of the vessels by 4 - 6, default 500;
%           .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.
%           .verbose : Show debug information, default true
%
% OUTPUTS:
%       Iout: [IxJxK] matrix representing filtered image stack
%       whatScale: [IxJxK] matrix of the maximum sigma range used to detect
%           object at that scale
%       VoutxN, VoutyN, VoutzN: [IxJxK] normalized eigenvectors numbered corresponding
%           to eigenvalue N (sorted smallest to largest). Note that x,y,z are kind of a
%           misnomer as they are really i, j, and k. Matlab orderes data in [row, col]
%           which is really [y, x] so for plotting purposes you want to plot
%           Vouty1,Voutx1,Voutz1
%       Lam1, Lam2, Lam3: [IxJxK] eigenvalues of each voxel, sorted lowest to largest
%           from 1-3
%
% PREREQUISITES:
%       need to have a functioning version of eig3volume using:
%           %   mex eig3volume.c
%
% BASED UPON:
%       http://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
%
% USAGE:
%       % run the filter
%       options = struct('ScaleRange', [0.1 2], 'ScaleRatio', 0.2, 'Alpha', 0.05, 'Beta', 1, 'Gamma', 250, 'verbose', false,'BlackWhite',false);
%       [Iout,whatScale,Voutx1,Vouty1,Voutz1,Voutx2,Vouty2,Voutz2,Voutx3,Vouty3,Voutz3,Lam1,Lam2,Lam3]=objectness3D(im,2,options);
%       % plot the output
%       n = 10; % let's just look at a random slice in the middle of the image stack
%       figure;imagesc(Iout(:,:,n));colormap gray;axis equal tight
%       figure;sl=slice(Iout,[],[],n); set(sl, 'EdgeColor','none');set(sl,'FaceAlpha',0.5); colormap gray
%       hold on;
%       [x,y,z] = meshgrid(1:size(im,1),1:size(im,2),1:size(im,3));
%       % note that Voutx1 and Vouty1 have swapped places compared to what you might logically expect
%       quiver3(x(:,:,n),y(:,:,n),z(:,:,n),whatScale(:,:,n).*Vouty1(:,:,n),whatScale(:,:,n).*Voutx1(:,:,n),whatScale(:,:,n).*Voutz1(:,:,n));
%       quiver3(x(:,:,n),y(:,:,n),z(:,:,n),whatScale(:,:,n).*Vouty2(:,:,n),whatScale(:,:,n).*Voutx2(:,:,n),whatScale(:,:,n).*Voutz2(:,:,n));
%       hold off;
%
% VERSION HISTORY:
%		FrangiFilter3D: Originally by Dirk-Jan Kroon, posted on Matlab file exchange
% 		objectness3D: Modified version by I. Campbell, 16 July 2013



% Constants vesselness function

defaultoptions = struct('ScaleRange', [0.1 2], 'ScaleRatio', 0.2, 'Alpha', 0.05, 'Beta', 1, 'Gamma', 250, 'verbose',false,'BlackWhite',false);

% Process inputs
if(~exist('options','var')), 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(options))), 
        warning('objectness:unknownoption','unknown options found');
    end
end

% Use single or double for calculations
if(~isa(I,'double')), I=single(I); end

sigmas=options.ScaleRange(1):options.ScaleRatio:options.ScaleRange(2);
sigmas = sort(sigmas, 'ascend');

% Frangi filter for all sigmas
for i = 1:length(sigmas),
    % Show progress
    if(options.verbose)
        disp(['Current Filter Sigma: ' num2str(sigmas(i)) ]);
    end
    
    % Calculate 3D hessian
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigmas(i));

    if(sigmas(i)>0)
        % Correct for scaling
        c=(sigmas(i)^2);
        Dxx = c*Dxx; Dxy = c*Dxy;
        Dxz = c*Dxz; Dyy = c*Dyy;
        Dyz = c*Dyz; Dzz = c*Dzz;
    end
    
    % Calculate eigen values
        [Lambda{1},Lambda{2},Lambda{3},Vx1,Vy1,Vz1,Vx2,Vy2,Vz2,Vx3,Vy3,Vz3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
%		% the following lines replace the preceeding if you can't get eig3volume to compile as a MEX file, but are very slow
%         Lambda{1}=zeros(size(Dxx));
%         Lambda{2}=Lambda{1};
%         Lambda{3}=Lambda{1};
%         Vx1=Lambda{1};Vy1=Lambda{2};Vz1=Lambda{3};
%         Vx2=Lambda{1};Vy2=Lambda{2};Vz2=Lambda{3};
%         Vx3=Lambda{1};Vy3=Lambda{2};Vz3=Lambda{3};
%         for l=1:size(Dxx,2)
%             for k=1:size(Dxx,3)
%                 for j=1:size(Dxx,1)
%                     [V,D]=eig([Dxx(j,k,l), Dxy(j,k,l), Dxz(j,k,l);Dxy(j,k,l),Dyy(j,k,l),Dyz(j,k,l);Dxz(j,k,l),Dyz(j,k,l),Dzz(j,k,l)]);
%                     [~,IX] = sort([abs(D(1,1)),abs(D(2,2)),abs(D(3,3))]);
%                     Lambda{1}(j,k,l)=D(IX(1),IX(1));
%                     Lambda{2}(j,k,l)=D(IX(2),IX(2));
%                     Lambda{3}(j,k,l)=D(IX(3),IX(3));
%                     Vx1(j,k,l)=V(1,IX(1));
%                     Vy1(j,k,l)=V(2,IX(1));
%                     Vz1(j,k,l)=V(3,IX(1));
%                     Vx2(j,k,l)=V(1,IX(2));
%                     Vy2(j,k,l)=V(2,IX(2));
%                     Vz2(j,k,l)=V(3,IX(2));
%                     Vx3(j,k,l)=V(1,IX(3));
%                     Vy3(j,k,l)=V(2,IX(3));
%                     Vz3(j,k,l)=V(3,IX(3));
%                 end;
%             end;
%         end;
    
    % Free memory
    clear Dxx Dyy  Dzz Dxy  Dxz Dyz;

    % Calculate absolute values of eigen values
    LambdaAbs{1}=abs(Lambda{1});
    LambdaAbs{2}=abs(Lambda{2});
    LambdaAbs{3}=abs(Lambda{3});

    % The Vesselness Features
    if M == 2
        Ra = inf;
    elseif M == 1
        Ra = LambdaAbs{M+1}./LambdaAbs{3}.^(1/(2-M));
    elseif M == 0
        Ra = LambdaAbs{M+1}./((LambdaAbs{2}.^(1/(2-M))) .* (LambdaAbs{3}.^(1/(2-M))));
    else
        error('objectness:M','M out of bounds');
    end
    
    if M == 0
        Rb = 0;
    elseif M == 1
        Rb = LambdaAbs{M} ./ ((LambdaAbs{2}.^(1/(3-M))) .* (LambdaAbs{3}.^(1/(3-M))));
    elseif M == 2
        Rb = LambdaAbs{M} ./ (LambdaAbs{3});
    else
        error('objectness:M','M out of bounds');
    end

    % Second order structureness. S = sqrt(sum(L^2[i])) met i =< D
    S = sqrt(LambdaAbs{1}.^2+LambdaAbs{2}.^2+LambdaAbs{3}.^2);
    A = 2*options.Alpha^2; B = 2*options.Beta^2;  C = 2*options.Gamma^2;
	
    % Free memory
    clear LambdaAbs1 LambdaAbs2 LambdaAbs3

    %Compute Vesselness function
    expRa = (1-exp(-(Ra.^2./A)));
    expRb = exp(-(Rb.^2./B));
    expS  = (1-exp(-S.^2./C));

    % Free memory
    clear S A B C Ra Rb

    %Compute Vesselness function
    Voxel_data = expRa .* expRb .* expS;
        
    % Free memory
    clear expRa expRb expRc;
    
    for j = M+1:3
        if(options.BlackWhite)
            Voxel_data(Lambda{j} < 0)=0;
        else
            Voxel_data(Lambda{j} > 0)=0;
        end
    end
        
%     % Remove NaN values
%     Voxel_data(~isfinite(Voxel_data))=0;
    if nnz(~isfinite(Voxel_data(:))) > 0
        display('NaN or Inf data found!');
    end

    % Add result of this scale to output
    if(i==1)
        Iout=Voxel_data;
        whatScale = ones(size(I),class(Iout));
        Voutx1=Vx1; Vouty1=Vy1; Voutz1=Vz1;
        Voutx2=Vx2; Vouty2=Vy2; Voutz2=Vz2;
        Voutx3=Vx3; Vouty3=Vy3; Voutz3=Vz3;
        Lam1 = Lambda{1}; Lam2 = Lambda{2}; Lam3 = Lambda{3};
    else
        whatScale(Voxel_data>Iout)=i;
        Voutx1(Voxel_data>Iout)=Vx1(Voxel_data>Iout);
        Vouty1(Voxel_data>Iout)=Vy1(Voxel_data>Iout);
        Voutz1(Voxel_data>Iout)=Vz1(Voxel_data>Iout);
        Voutx2(Voxel_data>Iout)=Vx2(Voxel_data>Iout);
        Vouty2(Voxel_data>Iout)=Vy2(Voxel_data>Iout);
        Voutz2(Voxel_data>Iout)=Vz2(Voxel_data>Iout);
        Voutx3(Voxel_data>Iout)=Vx3(Voxel_data>Iout);
        Vouty3(Voxel_data>Iout)=Vy3(Voxel_data>Iout);
        Voutz3(Voxel_data>Iout)=Vz3(Voxel_data>Iout);
        Lam1(Voxel_data>Iout)=Lambda{1}(Voxel_data>Iout);
        Lam2(Voxel_data>Iout)=Lambda{2}(Voxel_data>Iout);
        Lam3(Voxel_data>Iout)=Lambda{3}(Voxel_data>Iout);
        
        % Keep maximum filter response
        Iout=max(Iout,Voxel_data);
    end
end

