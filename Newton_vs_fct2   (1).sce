clear;  clf;
// Définition de la fonction de Rosenbrock
function [f]= fct1(x,y)
       f=(1-x)^2+100*(y-x^2)^2  
endfunction

//Définition du gradiant
function j = grad(x)
       j(1)=400*x(1)^3-400*x(1)*x(2)+2*x(1)-2;
       j(2)=200*(x(2)-x(1)^2);
endfunction

//Définition de la hessienne 
function [H]= hessi(x)
       H(1,1)=1200*x(1)^2-400*x(2)+2 ;     H(1,2)=-400*x(1);
       H(2,1)=-400*x(1);                               H(2,2)=200;
endfunction
    
// Représentation des lignes des niveaux pour la fonction de Rosenbrock
k=linspace(-3,3);
y=linspace(-3,3);    
z=feval(k,y,fct1);
xset("fpf"," ");
subplot(121)
plot3d(k,y,z)
subplot(122)
contour(k,y,z,80)

// La méthode de newton 
function [sol]=newton(x0,grad,hessi)
       N=10^6;
       eps=10^-4;
       xx=x0;
       i=0;
       tic();
       while(i<N)
              i=i+1;
              H = hessi(xx);
              xn = xx - inv(H)*grad(xx);
              plot(xn(1),xn(2),'g.');
              printf("iteration %d",i);
              disp(xn); printf("\n");
              if(norm(grad(xn))<eps) then  // solution trouvée 
                     t=toc();
                     printf("la solution est"); disp(xn);
                     printf("atteinte apres %d iterations \n",i);
                     printf("et apres %f secondes",t);
                     sol=resume(xn); 
              end; 
              xx=xn;
       end;
       sol=xn
       printf('pas de convergence apres %d iterations \n',i);
       abort;    // L'exécution s'arrêtera ici
endfunction

    
//intialisation
x=[-1,1.5]';
[sol]=newton(x);
