%%%% THIS CODE IDENTIFIES DIFFERENT REGIONS OF THE PHASE DIAGRAM BASED ON
%%%% THE COMPARISON OF DPHI-VS-INITIAL ANGLE WITH THAT FOR A CIRCULAR CASE,
%%%% WHICH HAS DPHI =0 FOR ALL ANGLES. IT CAN BE RUN AS SEPARATE CODE BY
%%%% SUPPLYING THE VALUE OF ecc, OR CAN BE CALLED IN THE CODE 'Plot_phase_diagram.m'
%%%% Also modified to take into account emergence of an oblique angle. If
%%%% we had the functional form of dphi, we could use findroot to determine
%%%% how many roots there are in the interval [0,pi/2]. Since we do not, we
%%%% use the criteria on the sign of dphi in this interval as well as the
%%%% sign of dphi near 0 and pi/2 to count the number of roots. There is no
%%%% tolerance set on how close dphi can be to 0, in order to identify the
%%%% stability boundaries correctly. Thus even a really small negative value of
%%%% dphi will be allowed to determine the nature of an FP. 

Lm_vec = 0.1:0.1:2;
ls_vec = 0.1:0.1:1;
ecc =0.2;
plt_phi = 0; 
z_st_pi2_re= zeros(size(Lm_vec,2),size(ls_vec,2));        %%%% zero stable pi/2 reachable
pi2_st_0_re = zeros(size(Lm_vec,2),size(ls_vec,2));       %%%% pi/2 reachable & stable, 0 reachable
only_z_reach = zeros(size(Lm_vec,2),size(ls_vec,2));      %%%% only zero reachable and stable
only_pi2_reach = zeros(size(Lm_vec,2),size(ls_vec,2));    %%%  only pi/2 reachable
both_not = zeros(size(Lm_vec,2),size(ls_vec,2));          %%%% both not reachable
both_yes = zeros(size(Lm_vec,2),size(ls_vec,2));          %%%% both FPs stable
only_oblique = zeros(size(Lm_vec,2),size(ls_vec,2)); %    %%% Only oblique angle stable: 0 either unreachble or unstable
ob_plus_0 = zeros(size(Lm_vec,2),size(ls_vec,2));   
%%%% Oblique angle and 0 stable. This corresponds to the case when 
%%%% there are two roots of dphi ==0 in the interval [0,pi/2], excluding
%%%% the end points. 

iLm=7; ils = 9;
dir  = 'Phase_bds_data';
for iLm = 1:size(Lm_vec,2)  
for ils = 1:size(ls_vec,2)
    [iLm ils]
     Lm= Lm_vec(iLm); 
      
     len_sp = ls_vec(ils);
    % Lm = 0.7; ls = 0.9;
        
     load(fullfile(dir,(['b=',num2str(1),'ecc=',num2str(ecc)...
         ,'lc=',num2str(len_sp),'Lm=',num2str(Lm),'.mat'])));
       
     figure(1)
     plot(rad2deg(phi_sp(:,1)),rad2deg(phi_sp(:,end)));  hold on;
     %%%% to compare the angular displacement with the no
     %%%% displacement case (circular)
     plot(rad2deg(phi_sp(:,1)),rad2deg(phi_sp(:,1))); hold on
     drawnow
     %figure;
     if plt_phi==1
        figure;
        for iphi=1:size(phi_initial,2)
        plot(phi_sp(iphi,:),'linewidth',1.5); drawnow; hold on; 
        end
     end
        
     Npts = size(phi_initial,2);         
     %%%% get the angular displacement
     dphi = phi_sp(:,end) - phi_sp(:,1); 
        
     %%%% define following to determine between the cases i) no
     %%%% displacement, ii) all angles converging to 0 and iii) all
     %%%% angles converging to pi/2, iv) Only pi/2 reachable and stable
     %%%% v) Only 0 reachable (this should never happen) and vi) both
     %%%% stable. 
        
     Iv = find(dphi(2:end-1)==0); 
     Ivneg =  find(dphi(2:end-1)<= 0);        
                                  
     if (size(Iv,1) == Npts-2) 
     %%% case(i) : no angle moves
        cse = 'No FP reachable' 
        both_not(iLm,ils) =1; 
            
     elseif (size(Iv,1) ~= Npts-2 && size(Ivneg,1) <1)
     %%%% All angles reach periphery and for all of them displacement is
     %%%% >0, i.e. they move towards pi/2, which is a stable FP.
     %%%% so displacement strictly positive. Both 0 and pi/2 reachable but pi/2 stable
        cse = 'both accessible but only pi/2 stable'
        pi2_st_0_reach(iLm,ils) = 1;   
            
     elseif (size(Iv,1) ~= Npts-2 && size(Ivneg,1) >=1)
     %%%% angular displacement is negative or 0 for some angles, 4 subcases
     %%%%  a) all negative : size(Ivneg) == Npts -2, 0 stable
     %%%%  b) some negative, others 0 : pi/2 not reachable, this case
     %%%%     should not be possible
     %%%%  c) some 0, others positive :  0 not reachable
     %%%%  d) some negative, some positive :  both stable and reachable
        
         if size(Ivneg,1) == Npts-2
         %%% all points converge to 0
             cse ='0 stable'
             z_st_pi2_re(iLm,ils) =1; 
         else
             %Ites = diff(Ivneg); 
             ind = Ivneg(end);   %%%% define the last point which converged to 0, excluding 0
             ind_trans = 1+ind;  %%%% raise ind by one to include 0
             %phi_ast(iLm,ils) = rad2deg(phi_initial(ind_trans));  %%%% transition point

             if ind_trans == size(phi_initial,2)-1  
             %%% checks the case when the last point befor pi/2 has negative displacement
             %%% if true: this case subdivides into 3 cases i) all angles have
             %%% negative displacement :  alrady covered in cse =
             %%% '0_stable', ii) angle near pi/2 move away from pi/2 and
             %%% angles near 0 move away from 0, so both unstable and an
             %%% oblique angle stable, iii) angle near 0 move towards 0,
             %%% and angles near pi/2 move away from pi/2, giving a stable
             %%% oblique angle and 0. This corresponds to a new root
             %%% emerging. 
             %%% check for dphi(2) 
                 if dphi(2) >=0    %%% 0 unreachable/unstable
                     cse = 'Only oblique stable'
                     only_oblique(iLm,ils) = 1;
                 else    %%% 0 stable
                     cse = 'Oblique and 0 stable'
                     ob_plus_0(iLm,ils)=1;
                 end         
             elseif ind_trans ~= size(phi_initial,2)-1
                 
                 %%%% angles that do not move near 0 or pi/2
                 v0 = find(dphi(1:ind_trans)==0); 
                 vpihalf = find(dphi(ind_trans+1:end)==0); 

                if size(v0,1) == ind_trans
                    %%% angles below phi_ast do not move
                    cse = 'Only pi/2 reachable & stable'
                    only_pi2_reach(iLm,ils) = 1; 
                elseif (size(vpihalf,1) == Npts-ind_trans )
                    %%% angles above phi_ast do not move
                    cse =  'only 0 reachable & stable' %%% possible when  lc/l >1  
                    only_z_reach(iLm,ils) = 1;
                elseif (size(find(dphi(2:ind_trans)<0),2) ~= 0 && size((dphi(ind_trans:end)>0),2)~=0)
                    %%% angles near 0 go to 0, angles near pi/2 go to pi/2
                    cse = 'both reachable and stable'               
                    both_yes(iLm,ils) = 1; 
                
                end
            
            end
        end
     end       
end
end