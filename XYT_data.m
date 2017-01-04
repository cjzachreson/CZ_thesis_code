classdef XYT_data
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
        L % linear side length of the space
%         Vx % velocity in x dimension
%         Vy %velocity in y dimension
%         V_mag % velocity magnitude
%         Vx_dir % x component of velocity direction
%         Vy_dir % y component of velocity direction
    end
    
    methods
        
        function this = XYT_data(dir, L)
            this.L = L;
            this.directory = dir;
            this.raw = dlmread(this.directory);
            this.s = nnz(this.raw(1,:))/3;
            
            this.X = this.raw(:, 1:3:end-2);
            this.Y = this.raw(:, 2:3:end-1);
            this.theta = this.raw(:, 3:3:end);
            this.raw = [];
%             this.Vx = diff(this.X, 1, 2);
%                 this.Vx = min(abs(this.Vx), 2*L - abs(this.Vx)) .* sign(this.Vx);
%             this.Vy = diff(this.Y, 1, 2);
%                 this.Vy = min(abs(this.Vy), 2*L - abs(this.Vy)) .* sign(this.Vy);
%             
%             this.V_mag = sqrt(this.Vx.^2 + this.Vy.^2);
%             this.Vx_dir = this.Vx ./ this.V_mag;
%             this.Vy_dir = this.Vy ./ this.V_mag;
        end
        
        function Quiver(this, t, pos)%, V)
            
            if pos
            
             quiver(this.X(:,t) - 1/2*cos(this.theta(:,t)) ...
             , this.Y(:,t) - 1/2*sin(this.theta(:,t)),...
             cos(this.theta(:,t)), sin(this.theta(:,t)),...
             'autoscale', 'off','ShowArrowHead', 'on')
             axis([-this.L, this.L, -this.L, this.L])
             axis square
             set(gca,'XTickLabel','');
             set(gca,'YTickLabel','');
            
            end
        
%             if (V && (t < this.s))
%            
%                 hold on
%                 quiver(this.X(:,t), this.Y(:,t), ...
%                     this.Vx(:,t).* 3, this.Vy(:,t) .* 3,...
%                 'autoscale', 'off','ShowArrowHead', 'on')
%                         axis([-this.L, this.L, -this.L, this.L])
%                 axis square
%                 set(gca,'XTickLabel','');
%                 set(gca,'YTickLabel','');
%                
%                 hold off
%             
%             end
            
            
        end
        
        
            
            
            
            
    end
    
end

