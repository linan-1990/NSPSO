% Description: This code shows an example of Non-dominated Sorting 
% -------------Multi-objective Particle Swarm Optimization Algorithm. 
% -------------The objective function is defined in fun.m, which contains
% -------------2 objective functions.
% Author: Nan LI (Nick)
% Company: RMIT Uniersity, Australia
% Email: s3468780@student.rmit.edu.au

clear all;
clc;

format long;

% Data = xlsread('Response Surface.xlsx');
% global Ti Vi CO2 PMV Energy;
% Ti = [17;18;19;20;21];
% Vi = [0.1;0.2;0.3;0.4;0.5];
% PMV1 = Data(1:5,3);
% PMV2 = Data(6:10,3);
% PMV3 = Data(11:15,3);
% PMV4 = Data(16:20,3);
% PMV5 = Data(21:25,3);
% PMV = [PMV1,PMV2,PMV3,PMV4,PMV5];
% 
% Energy1 = Data(1:5,4);
% Energy2 = Data(6:10,4);
% Energy3 = Data(11:15,4);
% Energy4 = Data(16:20,4);
% Energy5 = Data(21:25,4);
% Energy = [Energy1,Energy2,Energy3,Energy4,Energy5];

parameters = [0.5,1]; % the parameter vector of PSO
swarm_scale = 100; % the number of the swarm's elements
swarm_dimension = 2;
x_search_range = [-100,100];
y_search_range = [-100,100]; % the search range of y
x_range_min = x_search_range(1);
x_range_max = x_search_range(2);
y_range_min = y_search_range(1);
y_range_max = y_search_range(2);
vmax = [0.25*(x_range_max-x_range_min),0.25*(y_range_max-y_range_min)];

max_iter = 500; % the max iterations' number for first process

current_iter = 1; % the number of iterations processed

px = x_range_min + (x_range_max-x_range_min).*rand(swarm_scale,1); % Initial swarm.x
py = y_range_min + (y_range_max-y_range_min).*rand(swarm_scale,1); % Initial swarm.y
pos = [px,py]; % Initial swarm
v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
fitness1 = zeros(swarm_scale,1);
fitness2 = zeros(swarm_scale,1);
Front_cnt = 0;

while current_iter<max_iter
%     if(current_iter>5000)
%         parameters = [0.3,1];
%     end
    for l = 1:swarm_scale
        [fitness1(l),fitness2(l)] = fun(pos(l,1),pos(l,2)); % Calculate Fitness Value for every swarm
    end
    
    if current_iter==1 || current_iter==50
        pBest = pos; 
        pBest_Value = [fitness1,fitness2];
        [new_fitness1,index1] = sort(pBest_Value(:,1)); % Sort the objective function1's values for the swarm
        gBest = pos(index1(1),:); 
        gBest_Value = [fitness1(index1(1)),fitness2(index1(1))];
    else
        % start non-dominated sorting
        if mod(current_iter,2)==1
            [new_fitness1,index1] = sort(pBest_Value(:,1)); % Sort the objective function1's values for the swarm
            % [new_fitness2,index2] = sort(fitness2); % Sort the objective function2's values for the swarm
            Front_cnt = 1; % Initial the total amount of front_plane
            Element_cnt = ones(swarm_scale,1); % Initial the element amount in each front_plane
            flag = zeros(swarm_scale,1); % this flag indicates if this point has been clarified
            temp_flag = 0;
            sort_list = zeros(swarm_scale,swarm_scale);
            sort_list(1,1) = index1(1);
            for i = 1:swarm_scale
                % update pBest 
                if fitness1(i)<pBest_Value(i,1) && fitness2(i)<pBest_Value(i,2)
                    pBest(i) = pos(i); 
                    pBest_Value(i,:) = [fitness1(i),fitness2(i)];
                end
                if temp_flag == 1
                    temp_flag = 0;
                    Front_cnt = Front_cnt + 1; % Add a new front_plane
                    sort_list(Front_cnt,Element_cnt(Front_cnt)) = index1(i);
                end

                for j = 1:swarm_scale-i
                    if flag(i+j)==0
                        if pBest_Value((index1(i)),2)>pBest_Value((index1(i+j)),2) % non-dominate
                            Element_cnt(Front_cnt) = Element_cnt(Front_cnt) + 1;
                            sort_list(Front_cnt,Element_cnt(Front_cnt)) = index1(i+j);
                            flag(i+j) = Front_cnt; % set the flag to 1 showing that this point has been clarified
                            break;
                        else
                            temp_flag = 1;
                        end
                    end
                end  
            end
        else
            [new_fitness1,index1] = sort(pBest_Value(:,2)); % Sort the objective function1's values for the swarm
            % [new_fitness2,index2] = sort(fitness2); % Sort the objective function2's values for the swarm
            Front_cnt = 1; % Initial the total amount of front_plane
            Element_cnt = ones(swarm_scale,1); % Initial the element amount in each front_plane
            flag = zeros(swarm_scale,1); % this flag indicates if this point has been clarified
            temp_flag = 0;
            sort_list = zeros(swarm_scale,swarm_scale);
            sort_list(1,1) = index1(1);
            for i = 1:swarm_scale
                % update pBest 
                if fitness1(i)<pBest_Value(i,1) && fitness2(i)<pBest_Value(i,2)
                    pBest(i) = pos(i); 
                    pBest_Value(i,:) = [fitness1(i),fitness2(i)];
                end
                if temp_flag == 1
                    temp_flag = 0;
                    Front_cnt = Front_cnt + 1; % Add a new front_plane
                    sort_list(Front_cnt,Element_cnt(Front_cnt)) = index1(i);
                end

                for j = 1:swarm_scale-i
                    if flag(i+j)==0
                        if pBest_Value((index1(i)),1)>pBest_Value((index1(i+j)),1) % non-dominate
                            Element_cnt(Front_cnt) = Element_cnt(Front_cnt) + 1;
                            sort_list(Front_cnt,Element_cnt(Front_cnt)) = index1(i+j);
                            flag(i+j) = Front_cnt; % set the flag to 1 showing that this point has been clarified
                            break;
                        else
                            temp_flag = 1;
                        end
                    end
                end  
            end
        end
        % end non-dominated sorting
        
        temp = randi([1 Element_cnt(1)+Element_cnt(2)+Element_cnt(3)]);
        if temp<=Element_cnt(1)
            gBest = pBest(sort_list(1,temp),:); % identify the gBest point
            gBest_Value = pBest_Value(sort_list(1,temp),:);
        elseif temp<=(Element_cnt(1)+Element_cnt(2))
            gBest = pBest(sort_list(2,temp-Element_cnt(1)),:); % identify the gBest point
            gBest_Value = pBest_Value(sort_list(2,temp-Element_cnt(1)),:);
        else
            gBest = pBest(sort_list(3,temp-Element_cnt(1)-Element_cnt(2)),:); % identify the gBest point
            gBest_Value = pBest_Value(sort_list(3,temp-Element_cnt(1)-Element_cnt(2)),:);
        end
    end
    v = parameters(1).*v + parameters(2)*rand().*(pBest-pos) + parameters(2)*rand().*(repmat(gBest,swarm_scale,1)-pos); % Velocity
    
    for i=1:swarm_scale
        for j=1:swarm_dimension
            if v(i,j)>vmax(j)
               v(i,j)=vmax(j);
            elseif v(i,j)<-vmax(j)
                   v(i,j)=-vmax(j);
            end
        end
    end
%     if current_iter<50
%         for i=1:swarm_scale
%             v(i,2)=2*v(i,1)*rand();
%         end
%     end
    
    pos = pos+v; % Calculate the new Postion Vector
%     for p = 1:swarm_scale
%         for q = 1:swarm_dimension
%             if pos(p,q)>x_range_max
%                 pos(p,q)=x_range_max;
%             elseif pos(p,q)<x_range_min
%                 pos(p,q)=x_range_min;
%             end
%         end
%     end
    current_iter = current_iter + 1;
    
    if Front_cnt==1
        break;
    end
    
    if current_iter<50
        
        y_search_range = [-max(pos(:,1)),max(pos(:,1))]; % the search range of y
        
        y_range_min = y_search_range(1);
        y_range_max = y_search_range(2);
            
        pos(:,2)=y_range_min + (y_range_max-y_range_min).*rand(swarm_scale,1);
    end
    
    if current_iter==50
        x_search_range = [-10,10];
        y_search_range = [-2,2]; % the search range of y
        x_range_min = x_search_range(1);
        x_range_max = x_search_range(2);
        y_range_min = y_search_range(1);
        y_range_max = y_search_range(2);
        vmax = [0.25*(x_range_max-x_range_min),0.25*(y_range_max-y_range_min)];

        px = x_range_min + (x_range_max-x_range_min).*rand(swarm_scale,1); % Initial swarm.x
        py = y_range_min + (y_range_max-y_range_min).*rand(swarm_scale,1); % Initial swarm.y
        pos = [px,py]; 
        v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
    end
    
    subplot(1,2,1);
    plot(pBest_Value(:,1),pBest_Value(:,2),'b.','markersize',25); % plot the population
%             set(gca,'XLim',[0 10]);
%             set(gca,'YLim',[0 10]);
    set(gca,'FontSize',20);
    subplot(1,2,2);
%             plot(pBest(:,1),pBest(:,2),'g.'); % plot the population
    plot(pos(:,1),pos(:,2),'g.','markersize',25); % plot the population
    set(gca,'FontSize',20);
%         set(gca,'XLim',[-2 2]);
%         set(gca,'YLim',[-2 2]);
    pause(1/current_iter);
    
end
% plot(pBest_Value(:,1),pBest_Value(:,2),'b.'); % plot the population
disp('Done!');

% % start diversity control---Algorithm 1
% [sort_pBest,Index2] = sort(pBest_Value(:,1));
% pos = pBest;
% v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
% current_iter2 = 0;
% while current_iter2<max_iter2
%     for i = 2:swarm_scale-1
%         v(Index2(i)) = (pos(Index2(i+1))+pos(Index2(i-1)))/2-pos(Index2(i));  
%     end
%     pos = pos + v;
%     for l = 2:swarm_scale-1
%         [fitness1(l),fitness2(l)] = fun(pos(l,1),pos(l,2),pos(l,3));  % Calculate Fitness Value for every swarm except the pole two
%         if fitness1(l)>pBest_Value(l,1) && fitness2(i)>pBest_Value(l,2)
%             %
%         else
%             pBest(l) = pos(l); 
%             pBest_Value(l,:) = [fitness1(l),fitness2(l)];
%         end
%     end
%     [sort_pBest,Index2] = sort(pBest_Value(:,1));
%     pos = pBest;
%     v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
%     current_iter2 = current_iter2 +1;
%     
%     plot(pBest_Value(:,1),pBest_Value(:,2),'b.'); % plot the population
%     pause(0.01);
% end
% plot(pBest_Value(:,1),pBest_Value(:,2),'c.'); % plot the population
% disp('Done!')
% % end diversity control---Algorithm 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start diversity control---Algorithm 2
% [sort_pBest,Index2] = sort(pBest_Value(:,1));
% pos = pBest;
% v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
% distance = zeros(swarm_scale,swarm_dimension);
% current_iter2 = 0;
% while current_iter2<max_iter2
%     distance(1) = (pBest_Value(Index2(1),1)-pBest_Value(Index2(2),1))^2+(pBest_Value(Index2(1),2)-pBest_Value(Index2(2),2))^2;
%     for i = 2:swarm_scale-1
%         % calculate the distance
%         distance(i) = (pBest_Value(Index2(i),1)-pBest_Value(Index2(i+1),1))^2+(pBest_Value(Index2(i),2)-pBest_Value(Index2(i+1),2))^2;
%         if distance(i)>distance(i-1)
%             v(Index2(i)) = (pos(Index2(i+1))-pos(Index2(i)))/2;
%         else
%             v(Index2(i)) = (pos(Index2(i-1))-pos(Index2(i)))/2;
%         end    
%     end
%     pos = pos + v;
%     for l = 2:swarm_scale-1
%         [fitness1(l),fitness2(l)] = fun(pos(l,1),pos(l,2),pos(l,3));  % Calculate Fitness Value for every swarm except the pole two
%         if fitness1(l)>pBest_Value(l,1) && fitness2(i)>pBest_Value(l,2)
%             %
%         else
%             pBest(l) = pos(l); 
%             pBest_Value(l,:) = [fitness1(l),fitness2(l)];
%         end
%     end
%     [sort_pBest,Index2] = sort(pBest_Value(:,1));
%     pos = pBest;
%     v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
%     distance = zeros(swarm_scale,1);
%     current_iter2 = current_iter2 +1;
%     
%     plot(pBest_Value(:,1),pBest_Value(:,2),'b.'); % plot the population
%     pause(0.01);
% end
% plot(pBest_Value(:,1),pBest_Value(:,2),'c.'); % plot the population
% disp('Done!')
% end diversity control---Algorithm 2