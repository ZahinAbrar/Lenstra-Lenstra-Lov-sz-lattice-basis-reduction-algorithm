% original LLL Algorithm %
function [Basis]=LLL_modified(Basis)
[x,y] = meshgrid(-40:40,-40:40);
d=reshape([x,y],[],2);
Lattice = Basis*reshape([x,y],[],2)';
figure
scatter(Lattice(1,:),Lattice(2,:),'fill');
hold on
x_points= [0 Basis(1,1) Basis(1,1)+Basis(1,2) Basis(1,2) ];
y_points= [0 Basis(2,1) Basis(2,1)+Basis(2,2) Basis(2,2)];
patch(x_points,y_points,rand(1,3))
xlim([-10 40]);
ylim([-10 20]);
hold off


for i = 1:(size(Basis,2))
   b_star(:,i)=Basis(:,i);
   mu(i,i)=1;
    for j = 1: i-1
        mu(i,j)= (b_star(:,j)'*Basis(:,i))/n(1,j);
        b_star(:,i)= b_star(:,i)-mu(i,j)*b_star(:,j);
    end
    n(1,i)= b_star(:,i)'*b_star(:,i);
end

k =2 ;
  while(k<=size(Basis,2))
     [mu,Basis]=Reduce(k,k-1,mu,Basis);
  
        if n(1,k) < (3/4-(mu(k,k-1))^2)*n(1,k-1);
           mu_scalar= mu(k,k-1);
          B= n(1,k) + mu_scalar^2 * n(1,k-1);
          mu(k,k-1)= mu_scalar * n(1,k-1)/B;
          n(k) = n(1,k-1)*n(1,k) /B;
          n(1,k-1)= B;
          
          
          a=Basis(:,k-1);
          b=Basis(:,k);
          Basis(:,k-1)=b;
          Basis(:,k)=a;
          [x,y] = meshgrid(-40:40,-40:40);
          figure
        d=reshape([x,y],[],2);
        Lattice = Basis*reshape([x,y],[],2)';
        hold on
        scatter(Lattice(1,:),Lattice(2,:),'fill');
          x_points= [0 Basis(1,1) Basis(1,1)+Basis(1,2) Basis(1,2) ];
           y_points= [0 Basis(2,1) Basis(2,1)+Basis(2,2) Basis(2,2)];
           patch(x_points,y_points,rand(1,3))
           xlim([-10 40]);
            ylim([-10 20]);
          hold off
          for j = 1:k-2
                c=mu(k-1,j); 
                d=mu(k,j);
                mu(k,j)=c;
                mu(k-1,j)=d;
                
          end 
        
          for i = k+1 : size(Basis,2)
              t = mu(i,k);
              mu(i,k)= mu(i,k-1)-(mu_scalar*t);
              mu(i,k-1)= t + mu(k,k-1)* mu(i,k);
          end
                
        if k>2
        k = k-1;
        end
        continue
      end
  
        for l = k-2:-1:1
           [Basis,mu]=Reduce(k,k-1,mu,Basis); 
        end
            
        
           k = k+1
           [x,y] = meshgrid(-40:40,-40:40);
            d=reshape([x,y],[],2);
            Lattice = Basis*reshape([x,y],[],2)';
            figure()
            scatter(Lattice(1,:),Lattice(2,:),'fill');
            hold on
                x_points= [0 Basis(1,1) Basis(1,1)+Basis(1,2) Basis(1,2) ];
             y_points= [0 Basis(2,1) Basis(2,1)+Basis(2,2) Basis(2,2)];
             patch(x_points,y_points,rand(1,3))
           xlim([-10 40]);
            ylim([-10 20]);
          hold off
         
  end

                                                                                                          
                                
                              
                    
            