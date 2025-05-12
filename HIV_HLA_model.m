%This Matlab function implements the evolutionary epidemiological HIV HLA
%model explored in Herbert et al 2025

%The outputs of the function have been given explanatory names

%The inputs required by the model are as follows (see manuscript for further explanation):

%mu1=the background mortality rate of the population;
%beta1 = the HIV transmission parameter
%mu2_UU,mu2_II,mu2_FF,mu2_UI,mu2_FI,mu2_UF = genotype specfic HIV mortality
%rates in the absence of treatment
%transmission_prob_UU,transmission_prob_II,transmission_prob_FF,transmission_prob_UI,transmission_prob_FI,transmission_prob_UF = genotype specific vertical transmission probabilities for HIV
%age1=rate of aging from the child class to the adult class
%birth_rate= birth rate
%treatment=rate at which individuals are put onto ART
%time=length of simulation (unit depends on whatever units have been used
%for the rates parameters of the model
%init_pop = initial state of the population

function [T,Y,total_infected_not_treated,total_infected_and_treated,proportions_alleles,total_population,total_adult_pop,total_child_pop]=HIV_HLA_model(mu1,beta1,mu2_UU,mu2_II,mu2_FF,mu2_UI,mu2_FI,mu2_UF,transmission_prob_UU,transmission_prob_II,transmission_prob_FF,transmission_prob_UI,transmission_prob_FI,transmission_prob_UF,age1,birth_rate,treatment,time,init_pop)

pars=[mu1,beta1,mu2_UU,mu2_II,mu2_FF,mu2_UI,mu2_FI,mu2_UF,transmission_prob_UU,transmission_prob_II,transmission_prob_FF,transmission_prob_UI,transmission_prob_FI,transmission_prob_UF,age1,birth_rate,treatment];

[T,Y]= ode45(@equations, 0:1:time,init_pop,odeset('nonnegative',1:24,'AbsTol', 1e-14), pars);  %uses the solver ode45 to numerically solve the system of equations found in "equations"

total_population=sum(Y,2);

total_child_pop=sum(Y(:,1:4:21),2);

total_adult_pop=sum(Y,2)-total_child_pop;


total_infected_not_treated=sum(Y(:,3:4:23),2);
total_infected_and_treated=sum(Y(:,4:4:24),2);

totalUU=sum(Y(:,1:4),2);
totalII=sum(Y(:,5:8),2);
totalFF=sum(Y(:,9:12),2);
totalUI=sum(Y(:,13:16),2);
totalFI=sum(Y(:,17:20),2);
totalUF=sum(Y(:,21:24),2);

propU=(2*totalUU+totalUI+totalUF)./(2*total_population);
propI=(2*totalII+totalUI+totalFI)./(2*total_population);
propF=(2*totalFF+totalFI+totalUF)./(2*total_population);

proportions_alleles=[propU,propI,propF];


end


%this is the function "equations", which is formatted in a particular way
%so that it can be easily used by the solver in the other function, and so
%that it can use the parameters included in the other function as an input
%argument.
function dydt = equations(~,y,pars)

%the following lines of code convert each element of "pars" back to their
%individual names, so that you can recognise them more easily. 
mu1=pars(1);
beta1=pars(2);
mu2_UU=pars(3);
mu2_II=pars(4);
mu2_FF=pars(5);
mu2_UI=pars(6);
mu2_FI=pars(7);
mu2_UF=pars(8);
transmission_prob_UU=pars(9);
transmission_prob_II=pars(10);
transmission_prob_FF=pars(11);
transmission_prob_UI=pars(12);
transmission_prob_FI=pars(13);
transmission_prob_UF=pars(14);
age1=pars(15);
birth_rate=pars(16);
treatment=pars(17);

mu2=[mu2_UU;mu2_II;mu2_FF;mu2_UI;mu2_FI;mu2_UF];


total_adult_pop=0;
adults_by_geno=zeros(6,1);

for geno=1:1:6

marker=(geno-1)*4;

total_adult_pop=total_adult_pop+y(marker+2)+y(marker+3)+y(marker+4);

adults_by_geno(geno)=y(marker+2)+y(marker+3)+y(marker+4);

end

%working out force of infection
lambda=beta1*(sum(y(3:4:23))./total_adult_pop);

%genotype contribution matrices
UU_matrix=[1,0,0,0.5,0,0.5; 0,0,0,0,0,0; 0,0,0,0,0,0;     0.5,0,0,0.25,0,0.25; 0,0,0,0,0,0; 0.5,0,0,0.25,0,0.25];

II_matrix=[0,0,0,0,0,0; 0,1,0,0.5,0.5,0; 0,0,0,0,0,0;     0,0.5,0,0.25,0.25,0; 0,0.5,0,0.25,0.25,0; 0,0,0,0,0,0; ];

FF_matrix=[0,0,0,0,0,0; 0,0,0,0,0,0; 0,0,1,0,0.5,0.5;     0,0,0,0,0,0; 0,0,0.5,0,0.25,0.25; 0,0,0.5,0,0.25,0.25];

UI_matrix=[0,1,0,0.5,0.5,0; 1,0,0,0.5,0,0.5; 0,0,0,0,0,0;      0.5,0.5,0,0.5,0.25,0.25; 0.5,0,0,0.25,0,0.25; 0,0.5,0,0.25,0.25,0];

FI_matrix=[0,0,0,0,0,0;  0,0,1,0,0.5,0.5; 0,1,0,0.5,0.5,0;     0,0,0.5,0,0.25,0.25; 0,0.5,0.5,0.25,0.5,0.25; 0,0.5,0,0.25,0.25,0];

UF_matrix=[0,0,1,0,0.5,0.5; 0,0,0,0,0,0; 1,0,0,0.5,0,0.5;      0,0,0.5,0,0.25,0.25; 0.5,0,0,0.25,0,0.25;  0.5,0,0.5,0.25,0.25,0.5];

%working out assumed reproductive contributions of men and women

UU_contribF=adults_by_geno(1)-y(3)*transmission_prob_UU;
II_contribF=adults_by_geno(2)-y(7)*transmission_prob_II;
FF_contribF=adults_by_geno(3)-y(11)*transmission_prob_FF;
UI_contribF=adults_by_geno(4)-y(15)*transmission_prob_UI;
FI_contribF=adults_by_geno(5)-y(19)*transmission_prob_FI;
UF_contribF=adults_by_geno(6)-y(23)*transmission_prob_UF;

UU_prop_F=UU_contribF/(UU_contribF+II_contribF+FF_contribF+UI_contribF+FI_contribF+UF_contribF);
II_prop_F=II_contribF/(UU_contribF+II_contribF+FF_contribF+UI_contribF+FI_contribF+UF_contribF);
FF_prop_F=FF_contribF/(UU_contribF+II_contribF+FF_contribF+UI_contribF+FI_contribF+UF_contribF);
UI_prop_F=UI_contribF/(UU_contribF+II_contribF+FF_contribF+UI_contribF+FI_contribF+UF_contribF);
FI_prop_F=FI_contribF/(UU_contribF+II_contribF+FF_contribF+UI_contribF+FI_contribF+UF_contribF);
UF_prop_F=UF_contribF/(UU_contribF+II_contribF+FF_contribF+UI_contribF+FI_contribF+UF_contribF);


genotypes_prop_F=[UU_prop_F;II_prop_F;FF_prop_F;UI_prop_F;FI_prop_F;UF_prop_F];

genotypes_prop_M=adults_by_geno./total_adult_pop;


%working out a matrix of possible female x male genotype pairings...
%Female genotype proportions go down the rows, male genotype proportions go
%across the columns... Thus, each row of the matrix represents a single
%female genotype and how that genotype might pair with all the different
%male genotypes
Reprod_matrix=genotypes_prop_F*genotypes_prop_M';

%check_sum_reprod_matrix=sum(sum(Reprod_matrix))

%now we apportion the number of births we want in each time step to each
%genotype, according to our pre defined matrices which tell us the
%proportions of genotypes that should come from each type of pairing of two
%adult genotypes

totbirths=birth_rate*sum(y);

birth(1)=sum(sum(Reprod_matrix.*UU_matrix))*totbirths;
birth(2)=sum(sum(Reprod_matrix.*II_matrix))*totbirths;
birth(3)=sum(sum(Reprod_matrix.*FF_matrix))*totbirths;
birth(4)=sum(sum(Reprod_matrix.*UI_matrix))*totbirths;
birth(5)=sum(sum(Reprod_matrix.*FI_matrix))*totbirths;
birth(6)=sum(sum(Reprod_matrix.*UF_matrix))*totbirths;

dydt=[];

for geno=1:1:6

marker=(geno-1)*4;

dydt=[dydt;

birth(geno)-(age1+mu1)*y(marker+1);
age1*y(marker+1)-(lambda+mu1)*y(marker+2);
lambda*y(marker+2)-(mu1+mu2(geno)+treatment)*y(marker+3);
treatment*y(marker+3)-(mu1)*y(marker+4);

];

end

end





