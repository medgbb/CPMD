#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <armadillo>
#include <iomanip>

using namespace arma;
using namespace std;



double overlapIntegral(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);
double overlapIntegral_derivative(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);

double kineticIntegral(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);

double nuclearAttarctionIntegral(const int Z, const uint p, const uint q,
                                 const uint Rp, const uint Rq, vec alpha, mat R);

double electronInteractionIntegral(const int p, const int r, const int q, const int s,
                                   const int  Rp, const int  Rr, const int  Rq, const int  Rs,
                                   vec alpha,mat R, mat S);

double nuclearRepulsion(const mat R);

double errorFunction(double arg);

double calculateEnergy(double ****Q, const mat h, const vec C);
vec normalize(vec C, mat S);







double errorFunction_derivative(double arg);

double kineticIntegral_derivative(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);
double nuclearAttarctionIntegral_derivative(const int Z, const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R);
double electronInteractionIntegral_derivative(const int p, const int r, const int q, const int s,
                                              const int  Rp, const int  Rr, const int  Rq, const int  Rs,
                                              vec alpha, mat R);

int main()
{

    /*-----------------------------------------------------------------------------------------------------------*/


    //System configuration:
    uint nBasisFunc = 4;
    uint nNuclei    = 2;
    uint Z          = 1;

    uint e_nSteps  = 50;
    uint n_nSteps  = 1.0;

    double e_dt    = 0.1;
    double n_dt    = 4.3;
    double e_gamma = 1.0;
    double n_gamma = 5.0;

    double mu = 4.0;
    double M  = 1000*mu;

    double Eg;
    double a, b, c, lambda;


    vec alpha = zeros(nBasisFunc);
    alpha(0) = 13.00773;
    alpha(1) = 1.962079;
    alpha(2) = 0.444529;
    alpha(3) = 0.1219492;

    /*-----------------------------------------------------------------------------------------------------------*/
    //Initilize:
    uint nOrbitals = nBasisFunc * nNuclei;

    mat h = zeros(nOrbitals,nOrbitals);
    mat G = zeros(nOrbitals,nOrbitals);
    mat S = zeros(nOrbitals,nOrbitals);
    mat F = zeros(nOrbitals,nOrbitals);


    mat dh = zeros(nOrbitals,nOrbitals);
    mat dG = zeros(nOrbitals,nOrbitals);
    mat dS = zeros(nOrbitals,nOrbitals);
    mat dF = zeros(nOrbitals,nOrbitals);

    mat R      = zeros(nNuclei,3);
    mat Rplus  = zeros(nNuclei,3);
    mat Rminus = zeros(nNuclei,3);

    double q = 0.125;
    vec C      = ones(nOrbitals)*q;
    vec Cminus = ones(nOrbitals)*q;
    vec Cplus  = ones(nOrbitals)*q;

    vec2 lambdaVec;

    double ****Q;
    Q = new double***[nOrbitals];
    for (uint i = 0; i < nOrbitals; ++i) {
        Q[i] = new double**[nOrbitals];

        for (uint j = 0; j < nOrbitals; ++j){
            Q[i][j] = new double*[nOrbitals];

            for (uint k = 0; k < nOrbitals; ++k){
                Q[i][j][k] = new double[nOrbitals];
            }
        }
    }


    //One nuclei at -0.5ex and the other at 0.5ex
    Rminus(0,0) = -1.0;
    Rminus(1,0) =  1.0;

    R(0,0) = -0.5;
    R(1,0) =  0.5;

    /*-----------------------------------------------------------------------------------------------------------*/

    for(uint nStep = 0; nStep < n_nSteps; nStep++){


        //Set up the h and S matrix:
        for(uint Rp = 0; Rp < R.n_rows; Rp++){
            for(uint Rq = 0; Rq < R.n_rows; Rq++){

                for(uint p=0; p <alpha.size(); p++){
                    for(uint q=0; q <alpha.size(); q++){

                        S(p+Rp*4,q+Rq*4) = overlapIntegral(p,q,Rp,Rq,alpha,R);
                        h(p+Rp*4,q+Rq*4) = kineticIntegral(p,q,Rp,Rq,alpha,R)
                                + nuclearAttarctionIntegral(Z,p,q,Rp,Rq,alpha,R);
                    }
                }
            }

        }

        /*-----------------------------------------------------------------------------------------------------------*/
        //Set up the Q array:
        for(uint Rp = 0; Rp < R.n_rows; Rp++){
            for(uint Rr = 0; Rr < R.n_rows; Rr++){
                for(uint Rq = 0; Rq < R.n_rows; Rq++){
                    for(uint Rs = 0; Rs < R.n_rows; Rs++){

                        for(uint p=0; p <alpha.size(); p++){
                            for(uint r=0; r <alpha.size(); r++){
                                for(uint q=0; q <alpha.size(); q++){
                                    for(uint s=0; s <alpha.size(); s++){

                                        Q[p+Rp*4][r+Rr*4][q+Rq*4][s+Rs*4] = electronInteractionIntegral(p,r,q,s,Rp,Rr,Rq,Rs,alpha,R,S);
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }


        /*-----------------------------------------------------------------------------------------------------------*/

        //Loop over time
        for(uint step=0; step <= e_nSteps; step++){

            //Zero out elements
            G = zeros(nOrbitals,nOrbitals);
            Eg = 0.0;

            //Set up the G matrix:
            for(uint p=0; p < nOrbitals; p++){
                for(uint q=0; q < nOrbitals; q++){
                    for(uint r=0; r < nOrbitals; r++){
                        for(uint s=0; s < nOrbitals; s++){
                            G(p,q) += Q[p][r][q][s]*C(r)*C(s);
                        }
                    }
                }
            }

            /*-----------------------------------------------------------------------------------------------------------*/
            //Set up the F matrix
            F = h + G;
            /*-----------------------------------------------------------------------------------------------------------*/
            //Calculate Ctilde:
            Cplus = (2*C - (1-e_gamma*e_dt*0.5)*Cminus - e_dt*e_dt* F * C)/(1+e_gamma*e_dt*0.5);

            //Calculate lambda:
            a = dot(C, S*S*S*C);
            b = -2*dot(S*Cplus,S*C);
            c = dot(Cplus,S*Cplus)-1.0;

            lambdaVec(0)  = (-b + sqrt(b*b - 4*a*c)) /(2*a);
            lambdaVec(1)  = (-b - sqrt(b*b - 4*a*c)) /(2*a);

            if(lambdaVec(1) < 0 && lambdaVec(0) < 0 ){
                cerr << "negative roots!!" <<endl;
                cerr << lambdaVec <<endl;
                exit (EXIT_FAILURE);

            }

            if(lambdaVec(0) < lambdaVec(1) && lambdaVec(0) > 0.0 ){
                lambda = lambdaVec(0);
            }else{
                lambda = lambdaVec(1) ;
            }

            //Calculate C(t+h):
            Cplus -= lambda*S*C;

            //update C
            Cminus = C;
            C = Cplus;

            /*-----------------------------------------------------------------------------------------------------------*/
            //Calculate energy:
            Eg = calculateEnergy(Q,h,C);
            /*-----------------------------------------------------------------------------------------------------------*/


            cout.precision(8);
            cout <<"Energy: " << Eg <<" step: " << step << endl;
        }// Endof time loop electrons



        //Set up the dh and dS matrix:
        for(uint Rp = 0; Rp < R.n_rows; Rp++){
            for(uint Rq = 0; Rq < R.n_rows; Rq++){

                for(uint p=0; p <alpha.size(); p++){
                    for(uint q=0; q <alpha.size(); q++){

                        dS(p+Rp*4,q+Rq*4) = overlapIntegral_derivative(p,q,Rp,Rq,alpha,R);
                        dh(p+Rp*4,q+Rq*4) = kineticIntegral_derivative(p,q,Rp,Rq,alpha,R)
                                + nuclearAttarctionIntegral_derivative(Z,p,q,Rp,Rq,alpha,R);
                    }
                }
            }

        }

        /*-----------------------------------------------------------------------------------------------------------*/
        //Set up the Q array:
        for(uint Rp = 0; Rp < R.n_rows; Rp++){
            for(uint Rr = 0; Rr < R.n_rows; Rr++){
                for(uint Rq = 0; Rq < R.n_rows; Rq++){
                    for(uint Rs = 0; Rs < R.n_rows; Rs++){

                        for(uint p=0; p <alpha.size(); p++){
                            for(uint r=0; r <alpha.size(); r++){
                                for(uint q=0; q <alpha.size(); q++){
                                    for(uint s=0; s <alpha.size(); s++){

                                        Q[p+Rp*4][r+Rr*4][q+Rq*4][s+Rs*4] = electronInteractionIntegral(p,r,q,s,Rp,Rr,Rq,Rs,alpha,R,S);
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }

        dG = zeros(nOrbitals,nOrbitals);

        //Set up the G matrix:
        for(uint p=0; p < nOrbitals; p++){
            for(uint q=0; q < nOrbitals; q++){
                for(uint r=0; r < nOrbitals; r++){
                    for(uint s=0; s < nOrbitals; s++){
                        dG(p,q) += Q[p][r][q][s]*C(r)*C(s);
                    }
                }
            }
        }

        /*-----------------------------------------------------------------------------------------------------------*/
        //Set up the F matrix
        dF = dh + dG;
        /*-----------------------------------------------------------------------------------------------------------*/

        Rplus = 2*R - (1-n_gamma*n_dt*0.5)*Rminus - n_dt;


    }// End of time loop nuclei


    /*-----------------------------------------------------------------------------------------------------------*/

    //De-Allocate memory to prevent memory leak
    for (uint i = 0; i < nOrbitals; ++i) {
        for (uint j = 0; j < nOrbitals; ++j){
            for (uint k = 0; k < nOrbitals; ++k){
                delete [] Q[i][j][k];
            }
            delete [] Q[i][j];
        }
        delete [] Q[i];
    }
    delete [] Q;

    return 0;
}
/*################################################################################################################*/
/*End of main function
##################################################################################################################*/


double calculateEnergy(double ****Q, const mat h, const vec C){
    double Eg = 0.0;


    //One-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint q=0; q < C.n_elem; q++){
            Eg += C(p)*C(q)*h(p,q);
        }
    }
    Eg = 2*Eg;

    //Two-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint r=0; r < C.n_elem; r++){
            for(uint q=0; q< C.n_elem; q++){
                for(uint s=0; s < C.n_elem; s++){
                    Eg +=Q[p][r][q][s]*C(p)*C(q)*C(r)*C(s);
                }
            }
        }
    }



    //Nuclear repulsion term
    //    Eg +=nuclearRepulsion(R);

    return Eg;

}

/*-----------------------------------------------------------------------------------------------------------*/
vec normalize(vec C, mat S){
    double normFactor= 0.0;

    for(uint i= 0; i < C.n_elem; i++){
        for(uint j= 0; j < C.n_elem; j++){
            normFactor += C(i)*S(i,j)*C(j);
        }
    }
    return C/sqrt(normFactor);
}
/*-----------------------------------------------------------------------------------------------------------*/
double errorFunction(double t){

    if (t < 1.0E-6){
        return 1.0;
    }

    else{
        t = sqrt(t);
        double f = 1.0/t * erf(t)*sqrt(acos(-1))/2.0;
        return f;
    }

}

/*-----------------------------------------------------------------------------------------------------------*/
double overlapIntegral(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R){

    double factor = 1.0/(alpha(p)+ alpha(q));
    double Rpq = dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-alpha(p)*alpha(q)*factor*Rpq);
    double overlap = pow(acos(-1)*factor,3.0/2.0)*expTerm;

    return overlap;

}
/*-----------------------------------------------------------------------------------------------------------*/
double kineticIntegral(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R){

    double factor  = (alpha(p)*alpha(q))/(alpha(p)+ alpha(q));
    double Rpq = dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-factor*Rpq);
    double kin     = 0.5*factor*(6-4*factor*Rpq)*pow(acos(-1)/(alpha(p)+ alpha(q)),3.0/2.0)*expTerm;

    return kin;
}
/*-----------------------------------------------------------------------------------------------------------*/
double nuclearAttarctionIntegral(const int Z, const uint p, const uint q,
                                 const uint Rp, const uint Rq, vec alpha, mat R){

    double factor = 1.0/(alpha(p)+ alpha(q));
    double Rpq =dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq));
    double expTerm = exp(-alpha(p)*alpha(q)*factor*Rpq);

    rowvec Rmc = (alpha(p)*R.row(Rp) + alpha(q)*R.row(Rq))*factor;

    double F0p = errorFunction(1.0/factor*dot(Rmc-R.row(0),Rmc-R.row(0)));
    double F0q = errorFunction(1.0/factor*dot(Rmc-R.row(1),Rmc-R.row(1)));

    double nucAtt = -2*Z*factor*acos(-1)*expTerm*(F0p+F0q);

    return nucAtt;

}
/*-----------------------------------------------------------------------------------------------------------*/
double electronInteractionIntegral(const int p, const int r, const int q, const int s,
                                   const int Rp, const int Rr, const int Rq, const int Rs,
                                   vec alpha, mat R, mat S){


    double A = alpha[p] + alpha[q];
    double B = alpha[r] + alpha[s];

    rowvec Ra = (alpha[p]*R.row(Rp) + alpha[q]*R.row(Rq))/A;
    rowvec Rb = (alpha[r]*R.row(Rr) + alpha[s]*R.row(Rs))/B;


    double t = (A*B/(A + B))*dot(Ra-Rb,Ra-Rb);

    double arg = 2*sqrt(A*B/(acos(-1)*(A+B)))*errorFunction(t)*S(p+Rp*4,q+Rq*4)*S(s+Rs*4,r+Rr*4);

    return arg;


}
/*-----------------------------------------------------------------------------------------------------------*/
double nuclearRepulsion(const mat R){


    return 1/sqrt(dot(R.row(0) - R.row(1),R.row(0)- R.row(1)));
}




/*-----------------------------------------------------------------------------------------------------------*/
double errorFunction_derivative(double t){

    if (t < 1.0E-6){
        return 1.0/3.0;
    }

    else{
        double f = errorFunction(t);
        t = sqrt(t);
        return (exp(-t)-f)/(2*t) ;
    }

}


//*-----------------------------------------------------------------------------------------------------------*/
double overlapIntegral_derivative(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R){

    double factor  = (alpha(p)*alpha(q))/(alpha(p)+ alpha(q));
    double X     = sqrt(dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq)));
    double Spq     = overlapIntegral(p, q, Rp,Rq,alpha,R); // Use precalculated???
    double overlap = -2*factor*X*Spq;

    return overlap;

}
/*-----------------------------------------------------------------------------------------------------------*/
double kineticIntegral_derivative(const uint p, const uint q, const uint Rp, const uint Rq, vec alpha, mat R){

    double factor  = (alpha(p)*alpha(q))/(alpha(p)+ alpha(q));
    double X       = sqrt(dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq)));
    double dSdX    = overlapIntegral_derivative(p, q, Rp, Rq,alpha,R);
    double Spq     = overlapIntegral(p, q, Rp, Rq,alpha,R); // Use precalculated???
    double kin     = -4*factor*factor*X*Spq + (3*factor - 2*factor*factor*X*X)*dSdX;

    return kin;
}



/*-----------------------------------------------------------------------------------------------------------*/
double nuclearAttarctionIntegral_derivative(const int Z, const uint p, const uint q,
                                            const uint Rp, const uint Rq, vec alpha, mat R){

    double pq    = (alpha(p)+ alpha(q));
    double X     = sqrt(dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq)));
    double theta = 2*sqrt(pq/acos(-1));
    double Spq   = overlapIntegral(p, q, Rp,Rq,alpha,R);
    double nucAtt;


    if(X == 0){
        double t = pq*X;
        nucAtt   = 2*Z*theta*Spq*X*pq*errorFunction_derivative(t);
    }
    else{
        double dSdX    = overlapIntegral_derivative(p, q, Rp, Rq, alpha, R);
        rowvec P = (p*R.row(Rp)+ q*R.row(Rq))/pq;

        double F0p = errorFunction(X*X*alpha[p]*alpha[p]/pq);
        double F0q = errorFunction(X*X*alpha[q]*alpha[q]/pq);

        double dF0p = errorFunction_derivative(X*X*alpha[p]*alpha[p]/pq);
        double dF0q = errorFunction_derivative(X*X*alpha[q]*alpha[q]/pq);

        nucAtt = Z*theta*dSdX*(F0p + F0q) + 2*Z*(theta/pq)*Spq *(dF0p*alpha[p]*alpha[p] + dF0q*alpha[q]*alpha[q])*X;

    }


    return nucAtt;
}

/*-----------------------------------------------------------------------------------------------------------*/
double electronInteractionIntegral_derivative(const int p, const int r, const int q, const int s,
                                              const int Rp, const int Rr, const int Rq, const int Rs,
                                              vec alpha, mat R){

    double X     = sqrt(dot(R.row(Rp)-R.row(Rq),R.row(Rp)-R.row(Rq)));
    double A = alpha[p] + alpha[q];
    double B = alpha[r] + alpha[s];

    rowvec P = (p*R.row(Rp)+ q*R.row(Rq))/A;
    rowvec Q = (r*R.row(Rr)+ s*R.row(Rs))/B;

    double rho = 2*sqrt((A*B)/(acos(-1)*(A+B)));
    double t     = dot(P-Q,P-Q)*(A*B)/(A+B);


    double dSpqdX    = overlapIntegral_derivative(p, q, Rp, Rq, alpha, R);
    double dSrsdX    = overlapIntegral_derivative(r, s, Rr, Rs, alpha, R);

    double Spq = overlapIntegral(p, q, Rp, Rq, alpha, R);
    double Srs = overlapIntegral(r, s, Rr, Rs, alpha, R);


    double F0  = errorFunction(t);
    double dF0 = errorFunction_derivative(t);

    double arg = rho*dSpqdX*Srs*F0
            + rho*dSrsdX*Spq*F0
            + rho*Spq*Srs*dF0
            * (A*B)/(A+B) * 2*dot(P-Q,P-Q)/X;

    return arg;


}
