% Description: This code shows an example of Non-dominated Sorting 
% -------------Multi-objective Particle Swarm Optimization Algorithm. 
% -------------The objective function is defined in fun.m, which contains
% -------------2 objective functions.
% Author: Nan LI (Nick)
% Company: RMIT Uniersity, Australia
% Email: s3468780@student.rmit.edu.au

clear all;
clc;

parameters = [0.7,2]; % the parameter vector of PSO
swarm_scale = 50; % the number of the swarm's elements
swarm_dimension = 1;
x_search_range = [-1000,1000]; % the search range of x
max_iter = 1000; % the max iterations' number for first process
max_iter2 = 500; % the max iterations' number for diversity control

current_iter = 1; % the number of iterations processed

x_range_min = x_search_range(1);
x_range_max = x_search_range(2);

px = x_range_min + (x_range_max-x_range_min).*rand(swarm_scale,swarm_dimension); % Initial swarm.x
pos = px; % Initial swarm
v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
fitness1 = zeros(swarm_scale,1);
fitness2 = zeros(swarm_scale,1);
% selection = zeros(swarm_scale,1);
Front_cnt = 0;

% while current_iter<max_iter % the process of iteration
%     for l = 1:swarm_scale
%         [fitness1(l),fitness2(l)] = fun(pos(l)); % Calculate Fitness Value for every swarm
% %         plot(fitness1(:,1),fitness2(:,1),'b.'); % plot the population
%         if current_iter==1  % Initial the pBest Point
%             pBest(l) = pos(l); 
%             pBest_Value(l,:) = [fitness1(l),fitness2(l)]; 
%         elseif Selection_fun(fitness1(l),fitness2(l))<selection(l)
%             pBest(l) = pos(l); 
%             pBest_Value(l,:) = [fitness1(l),fitness2(l)]; 
%         end
%         selection(l) = Selection_fun(fitness1(l),fitness2(l));
%     end
%     
%     [sort_selection,index] = sort(selection);
%     temp = randi([1 fix(swarm_scale/2)]);
%     gBest = pos(index(temp*2)); % identify the gBest point
%     gBest_Value = [fitness1(index(1)),fitness2(index(1))];
%     v = parameters(1).*v + parameters(2)*rand().*(pBest-pos) + parameters(2)*rand().*(repmat(gBest,swarm_scale,1)-pos); % Velocity
%     pos = pos+v; % Calculate the new Position Vector
%     current_iter = current_iter + 1;
% end
% 
% plot(fitness1(:,1),fitness2(:,1),'b.'); % plot the population

while current_iter<max_iter
    for l = 1:swarm_scale
        [fitness1(l),fitness2(l)] = feval('fun',pos(l,:)); % Calculate Fitness Value for every swarm
    end
    
    if current_iter==1
        pBest = pos; 
        pBest_Value = [fitness1,fitness2];
        [new_fitness1,index1] = sort(pBest_Value(:,1)); % Sort the objective function1's values for the swarm
        gBest = pos(index1(1),:); 
        gBest_Value = [fitness1(index1(1)),fitness2(index1(1))];
    else
        % start non-dominated sorting
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
        % end non-dominated sorting
        
%         % start diversity control
%         if Element_cnt(1)>swarm_scale/3
%             distance = zeros(Element_cnt(1),1);
%             for k = 1:Element_cnt(1)-1
%                 distance(k) = (pBest_Value(sort_list(1,k),1)-pBest_Value(sort_list(1,k+1),1))^2+(pBest_Value(sort_list(1,k),2)-pBest_Value(sort_list(1,k+1),2))^2;   
%             end
%             distance(k+1)=distance(k); % for the last element, the distance equals to that of the privious one
%             [new_distance,Index]=sort(distance,'descend');
%             temp = randi([1 fix(Element_cnt(1)/50)+1]);
%             gBest = pBest(sort_list(1,Index(temp))); % identify the gBest point from the top 2% of the list
%             gBest_Value = pBest_Value(sort_list(1,Index(temp)));
%         else
%             temp = randi([1 Element_cnt(1)]);
%             gBest = pBest(sort_list(1,temp)); % identify the gBest point
%             gBest_Value = pBest_Value(sort_list(1,temp));
%         end
%         % end diversity control
            
        % no diversity control
        temp = randi([1 Element_cnt(1)]);
        gBest = pBest(sort_list(1,temp),:); % identify the gBest point
        gBest_Value = pBest_Value(sort_list(1,temp),:);
    end
    v = parameters(1).*v + parameters(2)*rand().*(pBest-pos) + parameters(2)*rand().*(repmat(gBest,swarm_scale,1)-pos); % Velocity
    pos = pos+v; % Calculate the new Position Vector
    current_iter = current_iter + 1;
    
    if Front_cnt==1
        break;
    end
    
    plot(pBest_Value(:,1),pBest_Value(:,2),'b.'); % plot the population
    pause(0.01);
end

% start diversity control---Algorithm 1
[sort_pBest,Index2] = sort(pBest_Value(:,1));
pos = pBest;
v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
current_iter2 = 0;
while current_iter2<max_iter2
    for i = 2:swarm_scale-1
        v(Index2(i)) = (pos(Index2(i+1))+pos(Index2(i-1)))/2-pos(Index2(i));  
    end
    pos = pos + v;
    for l = 2:swarm_scale-1
        [fitness1(l),fitness2(l)] = feval('fun',pos(l,:));  % Calculate Fitness Value for every swarm except the pole two
        if fitness1(l)>pBest_Value(l,1) && fitness2(i)>pBest_Value(l,2)
            %
        else
            pBest(l) = pos(l); 
            pBest_Value(l,:) = [fitness1(l),fitness2(l)];
        end
    end
    [sort_pBest,Index2] = sort(pBest_Value(:,1));
    pos = pBest;
    v = zeros(swarm_scale,swarm_dimension); % Initial swarm's velocity
    current_iter2 = current_iter2 +1;
    
    plot(pBest_Value(:,1),pBest_Value(:,2),'b.'); % plot the population
    pause(0.01);
end
plot(pBest_Value(:,1),pBest_Value(:,2),'b.'); % plot the population
disp('Done!')
% end diversity control---Algorithm 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % start diversity control---Algorithm 2
% [sort_pBest,Index2] = sort(pBest_Value(:,1));
% pos = pBest;
% v = zeros(swarm_scale,1); % Initial swarm's velocity
% distance = zeros(swarm_scale,swarm_dimension);
% current_iter2 = 0;
% while current_iter2<max_iter2
%     distance(1) = (pBest_Value(Index2(1),1)-pBest_Value(Index2(2),1))^2+(pBest_Value(Index2(1),2)-pBest_Value(Index2(2),2))^2;
%     for i = 2:swarm_scale-1
%         % calculate the distance
%         distance(i) = (pBest_Value(Index2(i),1)-pBest_Value(Index2(i+1),1))^2+(pBest_Value(Index2(i),2)-pBest_Value(Index2(i+1),2))^2;
%         if distance(i)>distance(i-1)
%             v(Index2(i)) = (pos(Index2(i+1))-pos(Index2(i)))./2;
%         else
%             v(Index2(i)) = (pos(Index2(i-1))-pos(Index2(i)))./2;
%         end    
%     end
%     pos = pos + v;
%     for l = 2:swarm_scale-1
%         [fitness1(l),fitness2(l)] = feval('fun',pos(l,:));  % Calculate Fitness Value for every swarm except the pole two
%         if fitness1(l)>pBest_Value(l,1) && fitness2(i)>pBest_Value(l,2)
%             %
%         else
%             pBest(l) = pos(l); 
%             pBest_Value(l,:) = [fitness1(l),fitness2(l)];
%         end
%     end
%     [sort_pBest,Index2] = sort(pBest_Value(:,1));
%     pos = pBest;
%     v = zeros(swarm_scale,1); % Initial swarm's velocity
%     distance = zeros(swarm_scale,1);
%     current_iter2 = current_iter2 +1;
%     
%     plot(pBest_Value(:,1),pBest_Value(:,2),'b.'); % plot the population
%     pause(0.01);
% end
% plot(pBest_Value(:,1),pBest_Value(:,2),'c.'); % plot the population
% disp('Done!')
% % end diversity control---Algorithm 2