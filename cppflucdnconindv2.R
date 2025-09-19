#fluctuation

# load Packages 
library(TMB);

#cpp 
cppflucdnconindv2 <- '
 //Fluctuation instead of the model in Schnute and Hilborn (1993, CJFAS);
 //The SH model has annual population sizes continue to increase. 
 //Fix N0 and K; //Otherwise the numerical optimization fails;  
 //Under the same conditions (with the constant N0 and K), we must compare the models
 //As of 18 Feb 2025
 
 #include <TMB.hpp> 

 // pass missing values
 template<class Type>
 bool isNA(Type x){
    return R_IsNA(asDouble(x));
 }

 // dlnorm
 template<class Type>
 Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0) {
 Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
 if(give_log) 
      return logres; 
 else 
      return exp(logres);
 };

 // square
 template<class Type>
 Type square(Type x){
    return pow(x,2.0); 
 };

 //objective function
 template<class Type>
 Type objective_function<Type>::operator() () {
 
 //data
 DATA_MATRIX(matIt);   //It: multi_fisheries CPUE indices as a matrix at time t 
                       //Under the same conditions (with the constant P0 and K), we must compare the models
 DATA_SCALAR(N0);     
 DATA_SCALAR(K);       //
 DATA_VECTOR(perrs); //standard deviation of the process errors  
 
 std::cout<<"matIt: "<<matIt<<std::endl;

 //parameters
 PARAMETER(logr); 

 //derived quantities
 Type r=exp(logr);
 
 int nt=matIt.rows();    //the number of time steps;
 int ncpue=matIt.cols(); //the number of cpue series;
 //int n=It.size();
 std::cout<<"nt and ncpue: "<<nt<<" "<<ncpue<<std::endl;
 
 //identity matrix
 matrix<Type> matrho(ncpue,ncpue);
 vector<Type> matdiag(ncpue);
 matrho.fill(Type(0));
 matdiag.fill(Type(1.0));
 matrho.diagonal()=matdiag;
 std::cout<<"matrho: "<<matrho<<std::endl;
 
 //process
 vector<Type> Nt(nt+1);   //the logistic pop growth 
 Nt(0)=N0; 
 for(int i=0;i<nt;i++)  
    Nt(i+1)=(Nt(i)+Nt(i)*r*(1-Nt(i)/K))*exp(perrs(i)); //the discrete logistic growth dynamics
    
 matrix<Type> matlogIt = log(matIt.array());  //transformation to logIt; 
 std::cout<<"matlogIt: "<<matlogIt<<std::endl; 
 
 matrix<Type> matlogpredIt(nt,ncpue);
 
 vector<Type> q(ncpue);
 matrix<Type> residuals_logItvslogNt(nt,ncpue);
 for(int j=0;j<ncpue;j++) {
   for(int i=0;i<nt;i++) 
      residuals_logItvslogNt(i,j)=log(matIt(i,j))-log(Nt(i));
      
   q(j)=exp( residuals_logItvslogNt.col(j).sum()/nt);       // the same as the form of MLE of q  
   //q(j)=exp( sum(residuals_logItvslogNt.col(j))/nt);       // the same as the form of MLE of q  
 };
     
 for(int j=0;j<ncpue;j++)
    for(int i=0;i<nt;i++)   
        matlogpredIt(i,j)=log(q(j))+log(Nt(i));
 
 std::cout<<"matlogpredIt: "<<matlogpredIt<<std::endl; 
        
 matrix<Type> matpredIt = exp(matlogpredIt.array()); // instead of //matrix<Type> matpredIt(nt,ncpue);     
                                                     //predicted It;matpredIt(i,j)=exp(matlogpredIt(i,j));
 
 matrix<Type> residualSquared(nt,ncpue);     //(logIt-log(q*Nt))^2;

 vector<Type> vecsig2_o(ncpue); 
 
 for(int j=0;j<ncpue;j++) {
   for(int i=0;i<nt;i++)   
      residualSquared(i,j)=square(matlogIt(i,j)-log(q(j)*Nt(i) ) ); 
      
   vecsig2_o(j)= residualSquared.col(j).sum()/nt;   //matrix.col(0) = 1st column of matrix
 };
 
 matrix<Type> Sigma(ncpue,ncpue); //variance-covariance matrix for the Its, 
                                  //which follow a multivariate normal distribution
 for(int j=0;j<ncpue;j++) {
   for(int i=0;i<ncpue;i++) 
     Sigma(j,i)=sqrt(vecsig2_o(j))*sqrt(vecsig2_o(i))*matrho(j,i);
 };
 std::cout<<"Sigma: "<<std::endl;
 std::cout<<Sigma<<std::endl; 
 
 Type nll=0.0;   //objective function as the negative logarithm; //initialization as 0.0;
 
 //=======================================================
 //observation equation
 vector<Type> residuals_logIt(ncpue);
 
 //For the multivariate normal distributions
 using namespace density;
 MVNORM_t<Type> neg_log_dmvnorm(Sigma);
 
 for(int i=0;i<nt;i++) {
      residuals_logIt = (matlogIt.row(i) - matlogpredIt.row(i));
      std::cout << "residuals_logIt: "  << residuals_logIt << std::endl;  
      nll += neg_log_dmvnorm(residuals_logIt);  //return the negative log density unlike the other density function
 };
 
 //report
 REPORT(nll);
 REPORT(matIt);
 REPORT(matpredIt);
 REPORT(q);
 REPORT(matrho);
 REPORT(vecsig2_o);
 REPORT(Sigma);
 REPORT(perrs);
 
 ADREPORT(r);
 ADREPORT(q);
 ADREPORT(Sigma);
 ADREPORT(Nt);
 
 return nll;
 
 }; 
'

#generate cpp file 
write(cppflucdnconindv2, file="cppflucdnconindv2.cpp"); #together

#compile
compile("cppflucdnconindv2.cpp");
dyn.load(dynlib("cppflucdnconindv2"));
