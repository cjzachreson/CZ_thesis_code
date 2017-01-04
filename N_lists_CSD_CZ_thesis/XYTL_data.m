classdef XYTL_data
    %contains the necessary data and functions for plotting an instance of the simulation
    %the root property is a string that defines the location of the data in
    %the directory, 
    
    properties 
        directory % location of the data
        raw % the raw, unsorted data as read from the text file
        s % size, the number of frames to be plotted
        X % X coordinates
        Y % Y coordinats
        theta % orientations
        l % lengths of each rod-shaped particle
        L % linear side length of the space
        Vx % velocity in x dimension
        Vy %velocity in y dimension
        V_mag % velocity magnitude
        Vx_dir % x component of velocity direction
        Vy_dir % y component of velocity direction
    end
    
    methods
        
        function this = XYTL_data(dir, L)
            this.L = L;
            this.directory = dir;
            this.raw = dlmread(this.directory);
            this.s = nnz(this.raw(1,:))/4;
            this.X = this.raw(:, 1:4:end-3);
            this.Y = this.raw(:, 2:4:end-2);
            this.theta = this.raw(:, 3:4:end-1);
            this.l = this.raw(:, 4:4:end);
            this.Vx = diff(this.X, 1, 2);
                this.Vx = min(abs(this.Vx), 2*L - abs(this.Vx)) .* sign(this.Vx);
            this.Vy = diff(this.Y, 1, 2);
                this.Vy = min(abs(this.Vy), 2*L - abs(this.Vy)) .* sign(this.Vy);
            
            this.V_mag = sqrt(this.Vx.^2 + this.Vy.^2);
            this.Vx_dir = this.Vx ./ this.V_mag;
            this.Vy_dir = this.Vy ./ this.V_mag;
        end
        
        function Quiver(this, t, pos, V)
            
            if pos
            
             quiver(this.X(:,t) - 1/2*cos(this.theta(:,t)) ...
             .* (this.l(:,t)), this.Y(:,t) - 1/2*sin(this.theta(:,t)).* (this.l(:,t)),...
             cos(this.theta(:,t)) .* (this.l(:,t)), sin(this.theta(:,t)) .* (this.l(:,t)),...
             'autoscale', 'off','ShowArrowHead', 'off')
             axis([-this.L, this.L, -this.L, this.L])
             axis square
             set(gca,'XTickLabel','');
             set(gca,'YTickLabel','');
            
            end
        
            if (V && (t < this.s))
           
                hold on
                quiver(this.X(:,t), this.Y(:,t), ...
                    this.Vx(:,t).* 3, this.Vy(:,t) .* 3,...
                'autoscale', 'off','ShowArrowHead', 'on')
                        axis([-this.L, this.L, -this.L, this.L])
                axis square
                set(gca,'XTickLabel','');
                set(gca,'YTickLabel','');
               
                hold off
            
            end
        end
            
        function Quiver_selected(this, t, indices, color, pos, V)
            
            if pos
            
             quiver(this.X(indices,t) - 1/2*cos(this.theta(indices,t)) ...
             .* (this.l(indices,t)), this.Y(indices,t) - 1/2*sin(this.theta(indices,t)).* (this.l(indices,t)),...
             cos(this.theta(indices,t)) .* (this.l(indices,t)), sin(this.theta(indices,t)) .* (this.l(indices,t)),...
             'autoscale', 'off','ShowArrowHead', 'off', 'Color', [color(1), color(2), color(3)], 'LineWidth', 3*1/(this.L/80))
             axis([-this.L, this.L, -this.L, this.L])
             axis square
             set(gca,'XTickLabel','');
             set(gca,'YTickLabel','');
             
            end
        
            if (V && (t < this.s))
           
              
                quiver(this.X(indices,t), this.Y(indices,t), ...
                    this.Vx(indices,t), this.Vy(indices,t),...
                'autoscale', 'off','ShowArrowHead', 'on','Color', [color(1), color(2), color(3)])
                        axis([-this.L, this.L, -this.L, this.L])
                axis square
                set(gca,'XTickLabel','');
                set(gca,'YTickLabel','');
               
               
            
            end
            
            
        end
        
        
            
            
            
            
    end
    
end

