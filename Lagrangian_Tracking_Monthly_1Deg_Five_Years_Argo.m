%% #########################################################################
%% Calculate Four-Year Lagrangian Trajectories based on Annual Mean Subduction Script
%  #########################################################################
%  2023/12/08, ZHI LI, UNSW SYDNEY CMSI, zhi.li4@unsw.edu.au
%  Monthly Mean, 2005-2021
%  Vertically Averaged over 700 m
%  Scripps RG-Argo Data
%  Geo Velocity are 
%     1. Ref to NM2000m
%     2. Ref to YoMaHa 1000m
%  #########################################################################

%  #########################################################################
%  2024/03/05, ZHI LI, UNSW SYDNEY CMSI, zhi.li4@unsw.edu.au
%  Releasing Water Particles from September
%  Then Tracking Paricles Monthly, 2005-2021 Monthly Climatology Geo Vel 
%  0-500 m Averaged
%  Ref to YoMaHa 1000m
%  #########################################################################

%% ########################################################################
%% CALCULATING THE MONTHLY LAGRANGIAN TRAJECTORIES
% #########################################################################
%  V2:   2025/02/20: Correct the dx dy and trajectory estimates 
% #########################################################################
%  V2.2: 2025/02/22: Added a checking and correction scheme to check if 
%  boundaries are correctly located around the new locations
%  2024/03/05, ZHI LI, UNSW SYDNEY CMSI, zhi.li4@unsw.edu.au
% #########################################################################
clc; clear

% SET UP ##################################################################
% Spatial resolution of the data
resolu = 1;  
% Bounadries of the trajectory's movement
limit_south = -64.5;    
limit_north = - 4.5;
limit_west  =   0.5;
limit_east  = 399.5; % Extend further east to 400E/40E

% Distance of each degree longitude/latitude
load('lat_RGArgo.mat')
lat_RGArgo=lat_RGArgo(1:61);
load('lon_RGArgo.mat')  
lon0=lon_RGArgo;
lon0(361:400,1)=(360.5:1:399.5)'; clear lon_RGArgo
lon_RGArgo=lon0;
clear lon0

[initial_x,initial_y]=meshgrid(lon_RGArgo, lat_RGArgo);
initial_x=initial_x'; initial_y=initial_y';

[~,dx,dy]=function_Cgrid_Area_Distance(lon_RGArgo,lat_RGArgo);
% #########################################################################


% #########################################################################
% Release water particles from moon=#  
for release_mon = 9
    clc
    disp('Water particles releasing from Sep......')  
    
    for year = 20052021
        disp(['year # ',num2str(year)]) 
        tic

        %% Ocean Geostrophic Current Velocities in Clim-Monthly Mean, RG-Argo, 2005-2021, 0-500 m Averaged
        % #########################################################################
        % U V YoMaHa@1000m
        disp('    loading YoMaHa@1000m geo_u v ......')
        cd('/Users/z5195509/Documents/Data/Argo_RG18/Geo_Vel/Ref#YoMaHa1000_Annual')
        load(['uv_geo_TW_YoHoMa_1000m_RG_20052021.mat'],'u_geo','v_geo','lon_RGArgo','lat_RGArgo','depth10');
        u_geo0(:,:,:,:,1)=u_geo; clear u_geo
        v_geo0(:,:,:,:,1)=v_geo; clear v_geo
        % #########################################################################
        
        % #########################################################################
%         % UV from NM-2000m
%         disp('    loading NM@2000m geo_u v ......')
%         cd('/Users/z5195509/Documents/Data/Argo_RG18/Geo_Vel/Ref#NM2000')
%         load(['uv_geo_TW_NM_2000m_RG18_20052021.mat'],'u_geo','v_geo');
%         u_geo0(:,:,:,:,2)=u_geo; clear u_geo
%         v_geo0(:,:,:,:,2)=v_geo; clear v_geo
        % #########################################################################
        u_geo=nanmean(u_geo0,5); clear u_geo0
        v_geo=nanmean(v_geo0,5); clear v_geo0
        
        u0(:,:,1:12)=squeeze(nanmean(u_geo(:,1:61,1:51,:),3)); 
        v0(:,:,1:12)=squeeze(nanmean(v_geo(:,1:61,1:51,:),3)); 
        u0(:,:,13:24)=squeeze(nanmean(u_geo(:,1:61,1:51,:),3)); 
        v0(:,:,13:24)=squeeze(nanmean(v_geo(:,1:61,1:51,:),3)); 
        u0(:,:,25:48)=u0(:,:,1:24);
        v0(:,:,25:48)=v0(:,:,1:24);
        u0(:,:,49:72)=u0(:,:,1:24);
        v0(:,:,49:72)=v0(:,:,1:24);
        clear u_geo v_geo month depth10 lon_RGArgo

        lat_RGArgo=lat_RGArgo(1:61,1);
        lon_RGArgo(:,1)=  (0.5:399.5)';
        u(  1: 20,:,:) =u0(341:360,:,:); %extended to 0-20E
        v(  1: 20,:,:) =v0(341:360,:,:); %extended to 0-20E
        u( 21:380,:,:) =u0(  1:360,:,:); %20E-380E
        v( 21:380,:,:) =v0(  1:360,:,:); %20E-380E
        u(381:400,:,:) =u0(  1: 20,:,:); %extended to 380-400E
        v(381:400,:,:) =v0(  1: 20,:,:); %extended to 380-400E
        clear u0 v0

        % 60 months scince the release month
        u=u(:,:,release_mon:release_mon+59);
        v=v(:,:,release_mon:release_mon+59);
        

        %% Pre-allocate Positions of Particles during Month=1-60
        [clonx,claty]=size(u(:,:,1));
        positions_x    =nan(clonx,claty,60); % lon,lat,month
        positions_y    =positions_x;
        positions_left =positions_x; 
        positions_right=positions_x; 
        positions_down =positions_x; 
        positions_up   =positions_x;

        % Ocean Current Velocities (Corrected from RK4 Metric)
        particle_u=positions_x; 
        particle_v=positions_x; 


        %% Initial Positions of Particles on Month=1
        disp('   Month #1')
        positions_x(:,:,1)=initial_x(:,:); % Initial x-coordinates
        positions_y(:,:,1)=initial_y(:,:); % Initial y-coordinates

        positions_left =positions_x;
        positions_right=positions_x;
        positions_down =positions_y;
        positions_up   =positions_y;


        %% Initial Ocean Current Velocity of Particles on Month=1
        particle_u(:,:,1)=u(:,:,1);
        particle_v(:,:,1)=v(:,:,1);

        disp('   particle_u(isnan(particle_u))=NaN')
        disp('   particle_v(isnan(particle_v))=NaN')
        particle_u(isnan(particle_u))=NaN;
        particle_v(isnan(particle_v))=NaN;
        u(isnan(u))=NaN;
        v(isnan(v))=NaN;


        %% LAGRANGIAN PARTICLE TRACKING FROM MONTH=2 TO MONTH=60
        for month = 1 : 59
            disp('  ')
            disp(['   Month #',num2str(month+1),' Year #',num2str(year)])    
            % Ocean Current Veocity at Current and Previous Months
            % month=month+1, from mon=2
            particle_v_curMON(:,:)=v(:,:,month+1); 
            particle_u_curMON(:,:)=u(:,:,month+1);

            % Using Runge-Kutta methods to get a better monthly mean velocity     
            disp(['    RK4 Velocity: Pariticles Move to Mon#',num2str(month+1),' Yr#',num2str(year)])
            % RK4 (month=month, from mon=1)
            particle_v_prvMON(:,:)=v(:,:,month);
            particle_u_prvMON(:,:)=u(:,:,month);
            
            %seconds_per_week=7*24*3600; % Seven days a week
            seconds_per_mon=30.5*24*3600;

                %% ############################################################
                %% Below Is To Derive The Particle Velocity Using RK-4 Method
                %  particle_u(:,:,1:23);
                %  particle_v(:,:,1:23);
                %% K1 (t0,x0,y0) 
                disp('        K1...')
                RK_vel_u1=particle_u(:,:,month);
                RK_vel_v1=particle_v(:,:,month);

                %% ############################################################
                %% K2 (t0+h/2,x0+h/2*RK_vel_u1,y0+h/2*RK_vel_v1)
                  disp('        K2...')
                  % new horizontal location of every grids:(lon1,lat1)*(u1,v1)-->(lon2,lat2) 
                  distance_x_RK   = RK_vel_u1.*seconds_per_mon./2; % moving distance for half month: velocity*time (u1,v1)*month/2=(dlon1,dlat1)
                  distance_y_RK   = RK_vel_v1.*seconds_per_mon./2;   
                  % Compute new positions (lon2, lat2)
                  positions_x_rk4 = distance_x_RK./dx*resolu + positions_x(:,:,month); % distance/grid distance = numeber of moving grid
                  positions_y_rk4 = distance_y_RK./dy*resolu + positions_y(:,:,month); % (lon1,lat1)+(dlon1,dlat1)=(lon2,lat2): new location

                  % Four neighbors around the new location
                  if month==1
                     positions_left_rk4  = positions_x(:,:,1) + floor(distance_x_RK./dx)*resolu; % West  boundary of the new location
                     positions_right_rk4 = positions_left_rk4 + resolu;                          % East  boundary
                     positions_down_rk4  = positions_y(:,:,1) + floor(distance_y_RK./dy)*resolu; % South boundary
                     positions_up_rk4    = positions_down_rk4 + resolu;                          % North boundary
                  else
                     positions_left_rk4  = positions_left(:,:,month) + floor(distance_x_RK./dx + (positions_x(:,:,month)-positions_left(:,:,month))./resolu)*resolu;  
                     positions_right_rk4 = positions_left_rk4        + resolu;
                     positions_down_rk4  = positions_down(:,:,month) + floor(distance_y_RK./dy + (positions_y(:,:,month)-positions_down(:,:,month))./resolu)*resolu;   
                     positions_up_rk4    = positions_down_rk4        + resolu; 
                  end   
                  clear distance_x_RK distance_y_RK
                      
                  % Interpolate mid-month UV to current location=K2
                  RK_vel_v_mid(:,:,1)=particle_v_prvMON; 
                  RK_vel_v_mid(:,:,2)=particle_v_curMON; 
                  RK_vel_v_mid=nanmean(RK_vel_v_mid,3);   % Mid-Month UV
                  RK_vel_u_mid(:,:,1)=particle_u_prvMON; 
                  RK_vel_u_mid(:,:,2)=particle_u_curMON; 
                  RK_vel_u_mid=nanmean(RK_vel_u_mid,3);   % Mid-Month UV
                  for i=1:clonx
                      for j=1:claty
                          x00=positions_left_rk4(i,j); x01=x00; % left grid
                          x10=positions_right_rk4(i,j);x11=x10; % right grid
                          y00=positions_down_rk4(i,j); y10=y00; % south grid
                          y01=positions_up_rk4(i,j);   y11=y01; % north grid
                          % check isnan=0
                          if isnan(positions_x_rk4(i,j))==0
                             % check within boundary, if yes, do the interpolation
                             if x00>limit_west && x10<=limit_east && y00>limit_south && y01<=limit_north % winthin boundary
                                % 2D Interpolation
                                number_nan=isnan(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu));
                                if number_nan==0 % no NaN
                                   % U
                                   RK_vel_u_mid1=(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)-...
                                                  RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)).*(positions_x_rk4(i,j)-x00)./resolu+...
                                                  RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu);
                                   RK_vel_u_mid1=squeeze(RK_vel_u_mid1);

                                   RK_vel_u_mid2=(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)-...
                                                  RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)).*(positions_x_rk4(i,j)-x01)./resolu+...
                                                  RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu);
                                   RK_vel_u_mid2=squeeze(RK_vel_u_mid2);   

                                   RK_vel_u_mid3=(RK_vel_u_mid2-RK_vel_u_mid1).*(positions_y_rk4(i,j)-y00)./resolu+RK_vel_u_mid1;
                                   clear RK_vel_u_mid1 RK_vel_u_mid2

                                   % V
                                   RK_vel_v_mid1=(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)-...
                                                  RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)).*(positions_x_rk4(i,j)-x00)./resolu+...
                                                  RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu);
                                   RK_vel_v_mid1=squeeze(RK_vel_v_mid1);

                                   RK_vel_v_mid2=(RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)-...
                                                  RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)).*(positions_x_rk4(i,j)-x01)./resolu+...
                                                  RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu);
                                   RK_vel_v_mid2=squeeze(RK_vel_v_mid2); 

                                   RK_vel_v_mid3=(RK_vel_v_mid2-RK_vel_v_mid1).*(positions_y_rk4(i,j)-y00)./resolu+RK_vel_v_mid1; %
                                   clear RK_vel_v_mid1 RK_vel_v_mid2

                                elseif number_nan>=1 && number_nan<=3 % 1-3 NaN
                                   % U
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)))=0;
                                   RK_vel_u_mid3(1,1)=(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu))./(4-number_nan);
                                   % V
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)))=0;
                                   RK_vel_v_mid3(1,1)=(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu))./(4-number_nan);    
                                elseif number_nan==4  %4/all NaN 
                                   RK_vel_u_mid3(1,1)=NaN;
                                   RK_vel_v_mid3(1,1)=NaN;
                                end % 2D interpolation end
                             else % out of boundary
                                RK_vel_u_mid3(1,1)=NaN;
                                RK_vel_v_mid3(1,1)=NaN;
                             end % boundary loop end

                             particle_ud(i,j)=RK_vel_u_mid3;
                             particle_vd(i,j)=RK_vel_v_mid3;
                             clear RK_vel_u_mid3 RK_vel_v_mid3

                          else % isnan==1
                             particle_ud(i,j)=NaN;
                             particle_vd(i,j)=NaN; 
                          end % isnan loop end
                      end %lat
                      clear x00 x01 x10 x11 y00 y01 y10 y11 
                  end % lon
                  clear i j
                  RK_vel_u2=particle_ud;
                  RK_vel_v2=particle_vd;
                  clear particle_ud particle_vd RK_vel_u_mid RK_vel_v_mid positions_x_rk4 positions_y_rk4
                  clear positions_left_rk4 positions_right_rk4 positions_down_rk4 positions_up_rk4


                %% ############################################################
                %% K3 (t0+h/2,x0+h/2*RK_vel_u2,y0+h/2*RK_vel_v2)
                  disp('        K3...')
                  % new horizontal location of every grids:(lon1,lat1)*(u1,v1)-->(lon2,lat2) 
                  distance_x_RK   = RK_vel_u2.*seconds_per_mon./2;   % moving distance for one month: velocity*time (u1,v1)*month=(dlon1,dlat1)
                  distance_y_RK   = RK_vel_v2.*seconds_per_mon./2;   % velocity on density class
                  positions_x_rk4 = distance_x_RK./dx*resolu + positions_x(:,:,month); % distance/grid distance = numeber of moving grid
                  positions_y_rk4 = distance_y_RK./dy*resolu + positions_y(:,:,month); % (lon1,lat1)+(dlon1,dlat1)=(lon2,lat2): new location

                  % Four neighbors around the new location
                  if month==1
                     positions_left_rk4  = positions_x(:,:,1) + floor(distance_x_RK./dx)*resolu; % West  boundary of the new location
                     positions_right_rk4 = positions_left_rk4 + resolu;                          % East  boundary
                     positions_down_rk4  = positions_y(:,:,1) + floor(distance_y_RK./dy)*resolu; % South boundary
                     positions_up_rk4    = positions_down_rk4 + resolu;                          % North boundary
                  else
                     positions_left_rk4  = positions_left(:,:,month) + floor(distance_x_RK./dx + (positions_x(:,:,month)-positions_left(:,:,month))./resolu)*resolu;  
                     positions_right_rk4 = positions_left_rk4        + resolu;
                     positions_down_rk4  = positions_down(:,:,month) + floor(distance_y_RK./dy + (positions_y(:,:,month)-positions_down(:,:,month))./resolu)*resolu;   
                     positions_up_rk4    = positions_down_rk4        + resolu; 
                  end   
                  clear distance_x_RK distance_y_RK
                  
                  % Interpolate mid-month UV to current location=K3
                  RK_vel_v_mid(:,:,:,1)=particle_v_prvMON; 
                  RK_vel_v_mid(:,:,:,2)=particle_v_curMON; 
                  RK_vel_v_mid=nanmean(RK_vel_v_mid,3);   % Mid-Month UV
                  RK_vel_u_mid(:,:,:,1)=particle_u_prvMON; 
                  RK_vel_u_mid(:,:,:,2)=particle_u_curMON; 
                  RK_vel_u_mid=nanmean(RK_vel_u_mid,4);   % Mid-Month UV
                  for i=1:clonx
                      for j=1:claty
                          x00=positions_left_rk4(i,j); x01=x00; % left grid
                          x10=positions_right_rk4(i,j);x11=x10; % right grid
                          y00=positions_down_rk4(i,j); y10=y00; % south grid
                          y01=positions_up_rk4(i,j);   y11=y01; % north grid
                          % check isnan=0
                          if isnan(positions_x_rk4(i,j))==0
                             % check within boundary, if yes, do the interpolation
                             if x00>limit_west && x10<=limit_east && y00>limit_south && y01<=limit_north % winthin boundary
                                % 2D Interpolation
                                number_nan=isnan(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu));
                                if number_nan==0 % no NaN
                                   % U
                                   RK_vel_u_mid1=(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)-...
                                                  RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)).*(positions_x_rk4(i,j)-x00)./resolu+...
                                                  RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu);
                                   RK_vel_u_mid1=squeeze(RK_vel_u_mid1);

                                   RK_vel_u_mid2=(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)-...
                                                  RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)).*(positions_x_rk4(i,j)-x01)./resolu+...
                                                  RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu);
                                   RK_vel_u_mid2=squeeze(RK_vel_u_mid2);   

                                   RK_vel_u_mid3=(RK_vel_u_mid2-RK_vel_u_mid1).*(positions_y_rk4(i,j)-y00)./resolu+RK_vel_u_mid1;
                                   clear RK_vel_u_mid1 RK_vel_u_mid2

                                   % V
                                   RK_vel_v_mid1=(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)-...
                                                  RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)).*(positions_x_rk4(i,j)-x00)./resolu+...
                                                  RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu);
                                   RK_vel_v_mid1=squeeze(RK_vel_v_mid1);

                                   RK_vel_v_mid2=(RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)-...
                                                  RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)).*(positions_x_rk4(i,j)-x01)./resolu+...
                                                  RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu);
                                   RK_vel_v_mid2=squeeze(RK_vel_v_mid2); 

                                   RK_vel_v_mid3=(RK_vel_v_mid2-RK_vel_v_mid1).*(positions_y_rk4(i,j)-y00)./resolu+RK_vel_v_mid1; %
                                   clear RK_vel_v_mid1 RK_vel_v_mid2

                                elseif number_nan>=1 && number_nan<=3 % 1-3 NaN
                                   % U
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)))=0;
                                   RK_vel_u_mid3(1,1)=(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu))./(4-number_nan);
                                   % V
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)))=0;
                                   RK_vel_v_mid3(1,1)=(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu))./(4-number_nan);  
                                elseif number_nan==4  %4/all NaN 
                                   RK_vel_u_mid3(1,1)=NaN;
                                   RK_vel_v_mid3(1,1)=NaN;
                                end % 2D interpolation end
                             else % out of boundary
                                RK_vel_u_mid3(1,1)=NaN;
                                RK_vel_v_mid3(1,1)=NaN;
                             end % boundary loop end

                             particle_ud(i,j)=RK_vel_u_mid3;
                             particle_vd(i,j)=RK_vel_v_mid3;
                             clear RK_vel_u_mid3 RK_vel_v_mid3

                          else % isnan==1
                             particle_ud(i,j)=NaN;
                             particle_vd(i,j)=NaN; 
                          end % isnan loop end
                      end %lat
                      clear x00 x01 x10 x11 y00 y01 y10 y11 
                  end % lon
                  RK_vel_u3=particle_ud;
                  RK_vel_v3=particle_vd;
                  clear particle_ud particle_vd RK_vel_u_mid RK_vel_v_mid positions_x_rk4 positions_y_rk4   
                  clear positions_left_rk4 positions_right_rk4 positions_down_rk4 positions_up_rk4


                %% ############################################################
                %% K4 (t0+h,x0+h*RK_vel_u3,y0+h*RK_vel_v3)
                  disp('        K4...')
                  % new horizontal location of every grids:(lon1,lat1)*(u1,v1)-->(lon2,lat2) 
                  distance_x_RK   = RK_vel_u3.*seconds_per_mon;   % moving distance for one month: velocity*time (u1,v1)*month=(dlon1,dlat1)
                  distance_y_RK   = RK_vel_v3.*seconds_per_mon;   % velocity on density class
                  positions_x_rk4 = distance_x_RK./dx*resolu + positions_x(:,:,month); % distance/grid distance = numeber of moving grid
                  positions_y_rk4 = distance_y_RK./dy*resolu + positions_y(:,:,month); % (lon1,lat1)+(dlon1,dlat1)=(lon2,lat2): new location

                  % Four neighbors around the new location
                  if month==1
                     positions_left_rk4  = positions_x(:,:,1) + floor(distance_x_RK./dx)*resolu; % West  boundary of the new location
                     positions_right_rk4 = positions_left_rk4 + resolu;                          % East  boundary
                     positions_down_rk4  = positions_y(:,:,1) + floor(distance_y_RK./dy)*resolu; % South boundary
                     positions_up_rk4    = positions_down_rk4 + resolu;                          % North boundary
                  else
                     positions_left_rk4  = positions_left(:,:,month) + floor(distance_x_RK./dx + (positions_x(:,:,month)-positions_left(:,:,month))./resolu)*resolu;  
                     positions_right_rk4 = positions_left_rk4        + resolu;
                     positions_down_rk4  = positions_down(:,:,month) + floor(distance_y_RK./dy + (positions_y(:,:,month)-positions_down(:,:,month))./resolu)*resolu;   
                     positions_up_rk4    = positions_down_rk4        + resolu; 
                  end   
                  clear distance_x_RK distance_y_RK
                  
                  % Interpolate current month UV to current location=K4
                  RK_vel_v_mid(:,:,:)=particle_v_curMON;
                  RK_vel_u_mid(:,:,:)=particle_u_curMON;
                  for i=1:clonx
                      for j=1:claty
                          x00=positions_left_rk4(i,j); x01=x00; % left grid
                          x10=positions_right_rk4(i,j);x11=x10; % right grid
                          y00=positions_down_rk4(i,j); y10=y00; % south grid
                          y01=positions_up_rk4(i,j);   y11=y01; % north grid
                          % check isnan=0
                          if isnan(positions_x_rk4(i,j))==0                   
                             % check within boundary, if yes, do the interpolation
                             if x00>limit_west && x10<=limit_east && y00>limit_south && y01<=limit_north % winthin boundary
                                % 2D Interpolation
                                number_nan=isnan(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu))+...
                                           isnan(RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu));
                                if number_nan==0 % no NaN
                                   % U
                                   RK_vel_u_mid1=(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)-...
                                                  RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)).*(positions_x_rk4(i,j)-x00)./resolu+...
                                                  RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu);
                                   RK_vel_u_mid1=squeeze(RK_vel_u_mid1);

                                   RK_vel_u_mid2=(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)-...
                                                  RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)).*(positions_x_rk4(i,j)-x01)./resolu+...
                                                  RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu);
                                   RK_vel_u_mid2=squeeze(RK_vel_u_mid2);   

                                   RK_vel_u_mid3=(RK_vel_u_mid2-RK_vel_u_mid1).*(positions_y_rk4(i,j)-y00)./resolu+RK_vel_u_mid1;
                                   clear RK_vel_u_mid1 RK_vel_u_mid2

                                   % V
                                   RK_vel_v_mid1=(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)-...
                                                  RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)).*(positions_x_rk4(i,j)-x00)./resolu+...
                                                  RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu);
                                   RK_vel_v_mid1=squeeze(RK_vel_v_mid1);

                                   RK_vel_v_mid2=(RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)-...
                                                  RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)).*(positions_x_rk4(i,j)-x01)./resolu+...
                                                  RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu);
                                   RK_vel_v_mid2=squeeze(RK_vel_v_mid2); 

                                   RK_vel_v_mid3=(RK_vel_v_mid2-RK_vel_v_mid1).*(positions_y_rk4(i,j)-y00)./resolu+RK_vel_v_mid1; %
                                   clear RK_vel_v_mid1 RK_vel_v_mid2

                                elseif number_nan>=1 && number_nan<=3 % 1-3 NaN
                                   % U
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)))=0;
                                   RK_vel_u_mid(isnan(RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)))=0;
                                   RK_vel_u_mid3(1,1)=(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)+...
                                                       RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu))./(4-number_nan);
                                   % V
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)))=0;
                                   RK_vel_v_mid(isnan(RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)))=0;
                                   RK_vel_v_mid3(1,1)=(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)+...
                                                       RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu))./(4-number_nan);     
                                elseif number_nan==4  %4/all NaN 
                                   RK_vel_u_mid3(1,1)=NaN;
                                   RK_vel_v_mid3(1,1)=NaN;
                                end % 2D interpolation end
                             else % out of boundary
                                RK_vel_u_mid3(1,1)=NaN;
                                RK_vel_v_mid3(1,1)=NaN;
                             end % boundary loop end

                             particle_ud(i,j)=RK_vel_u_mid3;
                             particle_vd(i,j)=RK_vel_v_mid3;
                             clear RK_vel_u_mid3 RK_vel_v_mid3

                          else % isnan==1
                             particle_ud(i,j)=NaN;
                             particle_vd(i,j)=NaN; 
                          end % isnan loop end
                      end %lat
                      clear x00 x01 x10 x11 y00 y01 y10 y11 
                  end % lon
                  RK_vel_u4=particle_ud;
                  RK_vel_v4=particle_vd;
                  clear particle_ud particle_vd RK_vel_u_mid RK_vel_v_mid positions_x_rk4 positions_y_rk4   
                  clear positions_left_rk4 positions_right_rk4 positions_down_rk4 positions_up_rk4

                %% Particle Velocity Using RK-4 Method For Tracking Month:Month+1: U_rk4 = 1/6*(K1 + 2*K2 + 2*K3 + K4)
                  disp('        U_RK...')
                  Ku(:,:,1)=RK_vel_u1;
                  Ku(:,:,2)=RK_vel_u2;
                  Ku(:,:,3)=RK_vel_u2;
                  Ku(:,:,4)=RK_vel_u3;
                  Ku(:,:,5)=RK_vel_u3;
                  Ku(:,:,6)=RK_vel_u4;

                  Kv(:,:,1)=RK_vel_v1;
                  Kv(:,:,2)=RK_vel_v2;
                  Kv(:,:,3)=RK_vel_v2;
                  Kv(:,:,4)=RK_vel_v3;
                  Kv(:,:,5)=RK_vel_v3;
                  Kv(:,:,6)=RK_vel_v4;

                  particle_u(:,:,month)=nanmean(Ku,3);
                  particle_v(:,:,month)=nanmean(Kv,3);          
        %           particle_u(:,:,month)=(RK_vel_u1 + 2.*RK_vel_u2 + 2.*RK_vel_u3 + RK_vel_u4)./6;
        %           particle_v(:,:,month)=(RK_vel_v1 + 2.*RK_vel_v2 + 2.*RK_vel_v3 + RK_vel_v4)./6;
                  clear RK_vel_*1 RK_vel_*2 RK_vel_*3 RK_vel_*4 

                  clear particle_u_prvMON particle_v_prvMON 
                  particle_u(isnan(particle_u))=NaN;
                  particle_v(isnan(particle_v))=NaN;


              %% Moving Distance and New Position of Particles
                disp(['    New Positions Month# ',num2str(month+1)])
                % new horizontal location of every grids:(lon1,lat1)*(u1,v1)-->(lon2,lat2) 
                distance_x_RK            = particle_u(:,:,month).*seconds_per_mon;% moving distance for one month: velocity*time (u1,v1)*month=(dlon1,dlat1)
                distance_y_RK            = particle_v(:,:,month).*seconds_per_mon;% velocity on density class
                positions_x(:,:,month+1) = distance_x_RK./dx*resolu + positions_x(:,:,month); % distance/grid distance = numeber of moving grid
                positions_y(:,:,month+1) = distance_y_RK./dy*resolu + positions_y(:,:,month); % (lon1,lat1)+(dlon1,dlat1)=(lon2,lat2): new location

                % Four neighbors around the new location
                if month==1
                   positions_left(:,:,month+1)  = positions_x(:,:,1)          + floor(distance_x_RK./dx)*resolu;  % left hand of the new location
                   positions_right(:,:,month+1) = positions_left(:,:,month+1) + resolu;                           % right hand
                   positions_down(:,:,month+1)  = positions_y(:,:,1)          + floor(distance_y_RK./dy)*resolu;  % south side
                   positions_up(:,:,month+1)    = positions_down(:,:,month+1) + resolu;                           % north side
                else
                   positions_left(:,:,month+1)  = positions_left(:,:,month)   + floor(distance_x_RK./dx+(positions_x(:,:,month)-positions_left(:,:,month))./resolu)*resolu;  
                   positions_right(:,:,month+1) = positions_left(:,:,month+1) + resolu;
                   positions_down(:,:,month+1)  = positions_down(:,:,month)   + floor(distance_y_RK./dy+(positions_y(:,:,month)-positions_down(:,:,month))./resolu)*resolu;   
                   positions_up(:,:,month+1)    = positions_down(:,:,month+1) + resolu; 
                end   
                clear distance_x_RK distance_y_RK
                  
                % Checking if the boundaries are located correctly
                  for i = 1:size(positions_x,1)
                      for j = 1:size(positions_x,2)
                          if positions_x(i,j,month+1) - positions_left(i,j,month+1) >= resolu && positions_x(i,j,month+1) - positions_left(i,j,month+1) < 2*resolu
                             disp(['          X-boundary relocated at i#',num2str(i),' j#',num2str(j)])
                             positions_left(i,j,month+1)  = positions_left(i,j,month+1)  + resolu;
                             positions_right(i,j,month+1) = positions_right(i,j,month+1) + resolu;
                             
                          elseif positions_x(i,j,month+1) - positions_left(i,j,month+1) >= 2*resolu 
                             disp(['          Warning'])
                             disp(['          X-boundary relocated at i#',num2str(i),' j#',num2str(j)])
                             positions_left(i,j,month+1)  = positions_left(i,j,month+1)  + floor((positions_x(i,j,month+1) - positions_left(i,j,month+1))./resolu).*resolu;
                             positions_right(i,j,month+1) = positions_right(i,j,month+1) + floor((positions_x(i,j,month+1) - positions_left(i,j,month+1))./resolu).*resolu;
                          end
                          
                          if positions_y(i,j,month+1) - positions_down(i,j,month+1) >= resolu  && positions_y(i,j,month+1) - positions_down(i,j,month+1) < 2*resolu
                             disp(['          Y-boundary relocated at i#',num2str(i),' j#',num2str(j)])
                             positions_down(i,j,month+1)  = positions_down(i,j,month+1)  + resolu;
                             positions_up(i,j,month+1)    = positions_up(i,j,month+1)    + resolu;
                             
                          elseif positions_y(i,j,month+1) - positions_down(i,j,month+1) >= 2*resolu 
                             disp(['          Warning'])
                             disp(['          Y-boundary relocated at i#',num2str(i),' j#',num2str(j)])
                             positions_down(i,j,month+1)  = positions_down(i,j,month+1)  + floor((positions_y(i,j,month+1) - positions_down(i,j,month+1))./resolu).*resolu;
                             positions_up(i,j,month+1)    = positions_up(i,j,month+1)    + floor((positions_y(i,j,month+1) - positions_down(i,j,month+1))./resolu).*resolu;
                          end
                      end
                  end
                  
                % Ocean Current Velocity at New Positions, Month=month+1
                RK_vel_v_mid(:,:)=particle_v_curMON;
                RK_vel_u_mid(:,:)=particle_u_curMON;
                for i=1:clonx
                    for j=1:claty
                        x00=positions_left(i,j,month+1); x01=x00; % left grid
                        x10=positions_right(i,j,month+1);x11=x10; % right grid
                        y00=positions_down(i,j,month+1); y10=y00; % south grid
                        y01=positions_up(i,j,month+1);   y11=y01; % north grid

                        % check isnan=0
                        if isnan(positions_x(i,j,month+1))==0

                           % check within boundary, if yes, do the interpolation
                           if x00>limit_west && x10<=limit_east && y00>limit_south && y01<=limit_north % winthin boundary

                              % 2D Interpolation
                              number_nan=isnan(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu))+...
                                         isnan(RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu))+...
                                         isnan(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu))+...
                                         isnan(RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu));
                              if number_nan==0 % no NaN
                                 % U
                                 RK_vel_u_mid1=(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)-...
                                                RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)).*(positions_x(i,j,month+1)-x00)./resolu+...
                                                RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu);
                                 RK_vel_u_mid1=squeeze(RK_vel_u_mid1);

                                 RK_vel_u_mid2=(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)-...
                                                RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)).*(positions_x(i,j,month+1)-x01)./resolu+...
                                                RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu);
                                 RK_vel_u_mid2=squeeze(RK_vel_u_mid2); 

                                 RK_vel_u_mid3(1,1)=(RK_vel_u_mid2-RK_vel_u_mid1).*(positions_y(i,j,month+1)-y00)./resolu+RK_vel_u_mid1;
                                 clear RK_vel_u_mid1 RK_vel_u_mid2

                                 % V
                                 RK_vel_v_mid1=(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)-...
                                                RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)).*(positions_x(i,j,month+1)-x00)./resolu+...
                                                RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu);
                                 RK_vel_v_mid1=squeeze(RK_vel_v_mid1);

                                 RK_vel_v_mid2=(RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)-...
                                                RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)).*(positions_x(i,j,month+1)-x01)./resolu+...
                                                RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu);
                                 RK_vel_v_mid2=squeeze(RK_vel_v_mid2);  

                                 RK_vel_v_mid3(1,1)=(RK_vel_v_mid2-RK_vel_v_mid1).*(positions_y(i,j,month+1)-y00)./resolu+RK_vel_v_mid1; %
                                 clear RK_vel_v_mid1 RK_vel_v_mid2

                              elseif number_nan>=1 && number_nan<=3 % 1-3 NaN
                                 % U
                                 RK_vel_u_mid(isnan(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)))=0;
                                 RK_vel_u_mid(isnan(RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)))=0;
                                 RK_vel_u_mid(isnan(RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)))=0;
                                 RK_vel_u_mid(isnan(RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)))=0;
                                 RK_vel_u_mid3(1,1)=(RK_vel_u_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)+...
                                                     RK_vel_u_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)+...
                                                     RK_vel_u_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)+...
                                                     RK_vel_u_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu))./(4-number_nan);
                                 % V
                                 RK_vel_v_mid(isnan(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)))=0;
                                 RK_vel_v_mid(isnan(RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)))=0;
                                 RK_vel_v_mid(isnan(RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)))=0;
                                 RK_vel_v_mid(isnan(RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu)))=0;
                                 RK_vel_v_mid3(1,1)=(RK_vel_v_mid((x10-limit_west)/resolu,(y10-limit_south)/resolu)+...
                                                     RK_vel_v_mid((x00-limit_west)/resolu,(y00-limit_south)/resolu)+...
                                                     RK_vel_v_mid((x11-limit_west)/resolu,(y11-limit_south)/resolu)+...
                                                     RK_vel_v_mid((x01-limit_west)/resolu,(y01-limit_south)/resolu))./(4-number_nan);      
                              elseif number_nan==4  %4/all NaN 
                                 % Velocity
                                 RK_vel_u_mid3(1,1)=NaN;
                                 RK_vel_v_mid3(1,1)=NaN;
                              end % 2D interpolation end
                           else % out of boundary
                              RK_vel_u_mid3(1,1)=NaN;
                              RK_vel_v_mid3(1,1)=NaN;
                           end % boundary loop end

                           particle_ud(i,j)=RK_vel_u_mid3;
                           particle_vd(i,j)=RK_vel_v_mid3; 
                           clear RK_vel_u_mid3 RK_vel_v_mid3

                        else % isnan==1
                           particle_ud(i,j)=NaN;
                           particle_vd(i,j)=NaN; 
                        end % isnan loop end
                    end %lat
                    clear x00 x01 x10 x11 y00 y01 y10 y11 

                end % lon
                particle_u(:,:,month+1)=particle_ud;
                particle_v(:,:,month+1)=particle_vd;
                clear particle_ud particle_vd
                clear RK_vel_u_mid RK_vel_v_mid

                particle_u(isnan(particle_u))=NaN;
                particle_v(isnan(particle_v))=NaN;
                clear particle_u_curMON particle_v_curMON
                % toc

        end % Month Loop End


        %% Saving One Year Data
        %  #########################################################################
        cd('/Users/z5195509/Documents/5_Oceanic_Teleconnections_ITF_ACC/6_2_Lagrangian_Trajectory/')
        disp('  Saving 5-Year Lagrangian Trajectories from YoMaHa@1000m...')
        save(['cal_Lagrangian_Trajectory_Four_Year_RG18_',num2str(year),'_V11_R1_Monthly_Tracking_V22_500m_YoMaHa1000.mat'],...
              'positions_x','positions_y','particle_u','particle_v','lon_RGArgo','lat_RGArgo')

        clear positions_x positions_y particle_u particle_v
        clear u v 
        disp('  Finishing Calculation For 1-Year Lagrangian Trajectories...')
        % #########################################################################

    end % Year Loop
end % Month Loop
toc
% #########################################################################
% #########################################################################
