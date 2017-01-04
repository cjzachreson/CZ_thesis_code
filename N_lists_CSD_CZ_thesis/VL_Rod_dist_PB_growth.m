
function [ d_ij ] = VL_Rod_dist_PB_growth( dXP, dYP, theta_i, theta_j, l_i, l_j)
%Vega-Lago algorithm (in 2D) for shortest distance between two rods

% Vega and Lago algorithm for calculating shortest distance
            % between two rods see http://www.sklogwiki.org/SklogWiki/
            % index.php/Source_code_for_the_minimum_distance_between_two_rods
            % this version is in 2D
        
        
            
        Rij = [dXP, dYP];% vector connecting centroids of the two rods 
        Ui = [cos(theta_i), sin(theta_i)];% unit vector describing orientation of rod i
        Uj = [cos(theta_j), sin(theta_j)];% unit vector describing orientation of rod j
        XLiD2 = l_i/2; % half the length of rod i
        XLjD2 = l_j/2;% half the length of rod j
        % Ro2 = square of the shortest distance between the two rods
        
        % ********
        
        %Step 1 
        Rij2 = Rij(1)^2 + Rij(2)^2; % euclidean distance squared 
        RijEUi = Rij * Ui';
        RijEUj = Rij * Uj';
        UiEUj = Ui * Uj';
        CC = 1 - UiEUj^2;
        
        brk = false;
        
       while brk == false
        
        % check whether rods are parallel
        if CC < exp(-6)
            if sqrt(Rij2) >= (XLiD2 + XLjD2)
                
                % rods are parallel but in end-on orientation
                
                XLanda = XLiD2 * sign(RijEUi);

                XMu = XLjD2 * -sign(RijEUj);
   
                if (abs(XMu) > XLjD2)
                    XMu = XLjD2 * sign(XMu);
                end
                
              Ro2 = Rij2 + XLanda^2 + XMu^2 - 2 * XLanda * XMu * UiEUj ...
              + 2 * XMu * RijEUj - 2 * XLanda * RijEUi;
  
          brk = true;
          pt = 1;
          break
             
                %********
                
            else % rods are parallel and in side-on orientation

                %first asign to points at tips of opposite particle 
                
                XLanda1 = RijEUi + XLjD2 * sign(UiEUj);%good
                XLanda2 = RijEUi - XLjD2 * sign(UiEUj);%good
                
                XMu1 = -RijEUj + XLiD2 * sign(UiEUj);
                XMu2 = -RijEUj - XLiD2 * sign(UiEUj);               
                
                %then decide if they need to come in to extremities
                XLanda1(abs(XLanda1) > XLiD2) = XLiD2 * sign(XLanda1); 
                XLanda2(abs(XLanda2) > XLiD2) = XLiD2 * sign(XLanda2);
                
                XMu1(abs(XMu1) > XLjD2) = XLjD2 * sign(XMu1);
                XMu2(abs(XMu2) > XLjD2) = XLjD2 * sign(XMu2);  
                
                XLanda = (XLanda1 + XLanda2) / 2;
                XMu = (XMu1 + XMu2) / 2;
                
                
                 
            if RijEUi == 0 && sqrt(Rij2) <= (XLiD2 + XLjD2) && l_i == l_j
                     
                XLanda = 0;
                XMu = 0;
                
            end
                 
              Ro2 = Rij2 + XLanda^2 + XMu^2 - 2 * XLanda * XMu * UiEUj ...
              + 2 * XMu * RijEUj - 2 * XLanda * RijEUi;
          
           
                          %********
                %goto(160)rods are parallel go to the end
     
          brk = true;
          
 
          pt = 2;
          break
         
                %********

            end
        end % end of parallel check
        
 

        % Step 1
        % Evauate XLanda' and XMu'
        
        XLanda = (RijEUi - UiEUj * RijEUj) / CC;
        XMu = (-RijEUj + UiEUj * RijEUi) / CC;

        % Step 2
 
        if ((abs(XLanda) <= XLiD2) && (abs(XMu) <= XLjD2)) 
            
            %********
            %goto(160) last step
             Ro2 = Rij2 + XLanda^2 + XMu^2 - 2 * XLanda * XMu * UiEUj + ...
                 2 * XMu * RijEUj - 2 * XLanda * RijEUi;
              
             brk = true;
             pt = 3;
             break
        
            %********
        end
        
        % End Step 2
        
        AuxIi = abs(XLanda) - XLiD2;
        AuxIj = abs(XMu) - XLjD2;
        
        % Step 3

        if AuxIi > AuxIj
            
            XLanda = XLiD2 * sign(XLanda);
           
            % Step 4 
                
            XMu = XLanda * UiEUj - RijEUj;
   
        if (abs(XMu) > XLjD2)
            XMu = XLjD2 * sign(XMu);
        end
                % Back to Step 3
            else
                XMu = XLjD2 * sign(XMu);
                % back to Step 4
                XLanda = XMu * UiEUj + RijEUi;
                
        if (abs(XLanda) > XLiD2) 
            XLanda = XLiD2 * sign(XLanda);
        end
        
        end
        
              
                    % Step 5 
       Ro2 = Rij2 + XLanda^2 + XMu^2 - 2 * XLanda * XMu * UiEUj ...
              + 2 * XMu * RijEUj - 2 * XLanda * RijEUi;
          brk = true;
          pt = 4;
          
       end

        d_ij = sqrt(Ro2);
        d_ij = abs(d_ij);
        
        
        
        if imag(d_ij) ~= 0
            X_i = X_i
            Y_i = Y_i
            X_j = X_j
            X_i = X_i
            pt = pt
            pause
        end