Lm_vec = 0.1:0.1:2;
ls_vec = 0.1:0.1:1;
 
plt_phi = 0; 
z_stable = zeros(size(Lm_vec,2),size(ls_vec,2));    %%%% only zero stable
pi2_stable = zeros(size(Lm_vec,2),size(ls_vec,2));   %%%% only pi/2 reachable
z_mat = zeros(size(Lm_vec,2),size(ls_vec,2));    %%%% only zero reachable
pi2_mat = zeros(size(Lm_vec,2),size(ls_vec,2));  %%% only pi/2 reachable
both_not = zeros(size(Lm_vec,2),size(ls_vec,2));   %%%% both not reachable
both_yes = zeros(size(Lm_vec,2),size(ls_vec,2));   %%%% both FPs reachable

iLm=9; ils = 6;
dir  = 'Phase_bds_data';
for iLm = 1:size(Lm_vec,2)  
for ils = 1:size(ls_vec,2)
    [iLm ils]
        Lm= Lm_vec(iLm); 
        len_sp = ls_vec(ils);
        
        load(fullfile(dir,...
        (['b=',num2str(1),'ecc=',num2str(ecc),'lc=',num2str(len_sp),'Lm=',num2str(Lm),'.mat'])));
        figure(1)
        plot(rad2deg(phi_sp(:,1)),rad2deg(phi_sp(:,end)));  hold on;
        %%%% to compare the total angular displacement with the no
        %%%% displacement case (circular)
        plot(rad2deg(phi_sp(:,1)),rad2deg(phi_sp(:,1))); hold on
        drawnow
        
        if plt_phi==1
        figure;
        for iphi=1:45
        plot(phi_sp(iphi,:),'linewidth',1.5); drawnow; hold on; 
        end
        end
        
        Npts = size(phi_initial,2);         
        %%%% get the total angular displacement
        dphi = phi_sp(:,end) - phi_sp(:,1); 
        
        %%%% define following to determine between the cases i) no
        %%%% displacement, ii) all angles converging to 0 and iii) all
        %%%% angles converging to pi/2, iv) Only pi/2 reachable and stable
        %%%% v) Only 0 reachable (this should never happen) and vi) both
        %%%% stable. 
        
        Iv = find(dphi(2:end-1)==0); 
        Ivneg =  find(dphi(2:end-1)<= 0);        
                                  
        if (size(Iv,1) == Npts-2) 
        %%%% case(i) : no angle moves
            cse = 'No FP reachable' 
            both_not(iLm,ils) =1; 
            
        elseif (size(Iv,1) ~= Npts-2 && size(Ivneg,1) <1)
        %%%% some angles reach periphery and all of them end up near pi/2,
        %%%% so displacement strictly positive. Both 0 and pi/2 reachable but pi/2 stable
            cse = 'pi/2 stable'
            pi2_stable(iLm,ils) = 1;   
            
        elseif (size(Iv,1) ~= Npts-2 && size(Ivneg,1) >=1)
        %%%% angular displacement is negative or 0 for some angles, 4
        %%%% subcases
        %%%%  a) all negative : size(Ivneg) == Npts -2, 0 stable
        %%%%  b) some negative, others 0 : pi/2 not reachable, this case
        %%%%     should not be possible
        %%%%  c) some 0, others positive :  0 not reachable
        %%%%  d) some negative, some positive :  both stable and reachable
        
        if size(Ivneg,1) == Npts-2
            %%% all points converge to 0
            cse ='0 stable'
            z_stable(iLm,ils) =1; 
        else
            %Ites = diff(Ivneg); 
            ind = Ivneg(end);
            ind_trans = 1+ind; 
            phi_ast(iLm,ils) = rad2deg(phi_initial(ind_trans)); 
            
            v0 = find(dphi(1:ind_trans)==0); 
            vpi2 = find(dphi(ind_trans+1:end)==0); 
            
            if size(v0,1) == ind_trans
                %%% angles below phi_ast do not move
                cse = 'Only pi/2 reachable & stable'
                pi2_mat(iLm,ils) = 1; 
            elseif size(vpi2,1) == Npts-ind_trans
                %%% angles above phi_ast do not move
                cse =  'Only 0 reachable & stable'   %%%% this should not be possible. 
                z_mat(iLm,ils) = 1;
            else %(dphi(2:ind_trans) <=0 && dphi(ind_trans+1:end) >= 0)
                %%% both angles reachable and stable
                cse = 'both reachable and stable'
                both_yes(iLm,ils) = 1; 
                
            end
            
        end
        end
         
end
end