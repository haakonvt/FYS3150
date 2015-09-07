#include "DiffSolv.h"
#include "lib.h"

DiffusionSolvers::DiffusionSolvers(int Nx_arg, int Nt_arg, double dx_arg, double dt_arg, int cutoff_arg){
    // Constructor initializes different constants:
    pi = 4*atan(1);
    Nx = Nx_arg;
    Nt = Nt_arg;
    dx = dx_arg;
    dt = dt_arg;
    cutoff = cutoff_arg;
    dtxx  = dt_arg/(dx_arg*dx_arg);
    dtxx2 = 1/(1+4*dtxx);
}

void DiffusionSolvers::FTCS_2D(mat u){
    // This code ONLY work with equal spatial steplengths: dx = dy
    mat u_prev  = u;

    /* The time-evolving solution can be saved as a "3D matrix", cube in armadillo.
    // The 2D solution at timestep k is then: sol.slice(k)
    cube sol = cube(Nx,Nx,Nt+1);
    sol.slice(0) = u;  // Initial condition */

    // Triple loop over time and [x,y]-grid:
    for (t=1; t<=Nt; t++){
        for (i=1; i<=Nx-2; i++){
            for (j=1; j<=Nx-2; j++){
                //cout << "i= " << i << ", j= " << j << endl;
                u(i,j) = (1-4*dtxx)*u_prev(i,j) + dtxx*(u_prev(i+1,j) + u_prev(i-1,j) + u_prev(i,j+1) + u_prev(i,j-1));
                //cout << "uprev(i,j) = " << u_prev(i,j) << ", u(i,j) = " <<u(i,j)<<endl;
            }
        }
        //sol.slice(t) = u; // Save current solution
        u_prev       = u;   // Update u_previous for next loop
    }
    // Save the time-evolving solution to file
    //sol.save("ExplForwEuler2D.dat",raw_ascii);

    // Only save 'last' solution (typically steady state solution for large t)
    u.save("ExplForwEuler2D.dat",raw_ascii);
    //cout <<u<<endl;

}




void DiffusionSolvers::BTCS_2D(mat u){
    // This code ONLY work with equal spatial steplengths: dx = dy and _known_ bcs

    double gamma = 1 + 4*dtxx;
    int lim   = Nx-2;
    int limsq = lim*lim;

    /*TEST SOLVE(A,b)
    u = randu(Nx,Nx);
    cout << u << endl;*/

    // Write u(x,y,t=const.) as a vector. Then u(i,j) = u_vec(lim*(i-1)+j-1)
    mat u_inner_matrix =  u.submat(1, 1, lim, lim);
    vec u_vec = vectorise(u_inner_matrix.t());  // 'Vectorise' concatenates columns, so a transpose is needed

    // SAVE SOLUTION WITH TIME EVOLUTION
    //mat sol    = zeros(limsq,Nt);    // Solution matrix with columns being u(x,y) in "vector-form"
    //sol.col(0) = u_vec;

    // Setup of matrix A in Au=b
    mat fivebandmat = zeros(limsq,limsq);
    fivebandmat.diag(0)    =  gamma*ones(limsq);
    fivebandmat.diag(1)    = -dtxx*ones(limsq-1);
    fivebandmat.diag(-1)   = -dtxx*ones(limsq-1);
    fivebandmat.diag(lim)  = -dtxx*ones(limsq-lim);
    fivebandmat.diag(-lim) = -dtxx*ones(limsq-lim);

    j=1;
    for (i=1; i<=limsq-1; i++){
        if (j == lim){
            fivebandmat.diag( 1)(i-1) = 0;
            fivebandmat.diag(-1)(i-1) = 0;
            j = 0;
        }
        j++;
    }

    for(t=1; t<Nt; t++){

        // Setup of RHS, b in Au=b
        vec b       = zeros(limsq);
        int b_index = 0;
        for (i=1; i<=lim; i++){
            for (j=1; j<=lim; j++){

                // Only known values (boundary cond.) can be put into RHS, b:
                if(i==1){
                    b(b_index) += dtxx*u(i-1,j);
                }
                if(j==1){
                    b(b_index) += dtxx*u(i,j-1);
                }
                if(i==lim){
                    b(b_index) += dtxx*u(i+1,j);
                }
                if(j==lim){
                    b(b_index) += dtxx*u(i,j+1);
                }

                b(b_index) += u_vec(lim*(i-1) + j-1); //This translates to element u_with_boundaries(i,j)
                b_index++;

                /* CHECK THAT MAT AND VEC GIVE SAME ANSWER
            cout << "mat: " << u_inner_matrix(i-1,j-1) << endl;
            cout << "vec: " << u_vec(lim*(i-1) + j-1) << endl;
            */

            }
        }

        // Solve Au=b, where u is written as a vector.
        // See report for details.
        u_vec = solve(fivebandmat,b);

        // NB: Solution only saves inner matrix (written as vector)
        // sol.col(t) = u_vec;

        // b.save("b.dat",raw_ascii);
        // u_vec = solve(fivebandmat,b);
        // cout << "new u:"<<endl;
        // cout << u_vec<<endl;
    }

    // sol.save("btcs_2d.dat",raw_ascii);

    // SAVE JUST SOLUTION AT LAST TIMESTEP
    u_vec.save("BTCS2D.dat",raw_ascii);

}




vec DiffusionSolvers::TriDiag(vec aa, vec ab, vec ac, vec b, int Nx){

    vec v = zeros(Nx); // Solution vector

    //aa(0) = 0;
    //ac(Nx-1) = 0;

    // Forward substitution
    ac(0) = ac(0)/ab(0);
    b(0)  = b(0)/ab(0);

    for (i=1; i<=Nx-2; i++){
        ac(i) = ac(i) / (ab(i) - aa(i)*ac(i-1));
        b(i)  = (b(i) - aa(i)*b(i-1)) / (ab(i) - aa(i)*ac(i-1));
    }

    b(Nx-1) = (b(Nx-1)-aa(Nx-1)*b(Nx-2)) / (ab(Nx-1)-aa(Nx-1)*ac(Nx-2));

    // Backwards substitution gives the solutions
    v(Nx-1) = b(Nx-1);

    for (i=Nx-1; i>=1; i--){
        v(i-1) = b(i-1) - ac(i-1)*v(i);
    }
    return v;

}


void DiffusionSolvers::ExplForwEuler(vec &v){
    vec v_next = zeros(Nx);
    mat sol    = zeros(Nx,Nt);
    vec us     = 1-linspace(0,1,Nx);

    sol.col(0) = v+us; // Save real solution u = v + us

    double dtxx = dt/(dx*dx);

    for(j=1; j<Nt; j++){
        for(i=1; i<Nx-1; i++){
            v_next(i) = dtxx*v(i-1) + (1-2*dtxx)*v(i) + dtxx*v(i+1);
        }
        v = v_next;

        sol.col(j) = v+us;
    }
    sol.save("ExplForwEuler.dat",raw_ascii);
}


void DiffusionSolvers::ImplBackEuler(vec &v){
    vec a = -dtxx*ones(Nx);
    vec b = (1 + 2*dtxx)*ones(Nx);
    vec c = a;

    mat sol    = zeros(Nx,Nt);
    vec us     = 1-linspace(0,1,Nx);
    sol.col(0) = v+us; // Save real solution u = v + us

    for(j=1; j<Nt; j++){
        v = DiffusionSolvers::TriDiag(a,b,c,v,Nx);

        sol.col(j) = v+us;
    }

    sol.save("ImplBackEuler.dat",raw_ascii);
}


void DiffusionSolvers::CrankNicolson(vec &v){
    double d1 = dtxx;
    double d2 = 2-2*d1;
    double d3 = 2+2*d1;

    vec a = -d1*ones(Nx);
    vec b = d3*ones(Nx);
    vec c = a;

    vec vv     = v;
    mat sol    = zeros(Nx,Nt);
    vec us     = 1-linspace(0,1,Nx);
    sol.col(0) = v+us; // Save real solution u = v + us

    for (j=1; j<Nt; j++){
        for (i=1; i<Nx-1; i++){
            vv(i) = d1*v(i-1) + d2*v(i) + d1*v(i+1);

        }
        v = DiffusionSolvers::TriDiag(a, b, c, vv, Nx);
        sol.col(j) = v+us;
    }
    sol.save("CrankNicolson.dat",raw_ascii);

}


void DiffusionSolvers::AnalyticSolution(vec u, vec x, vec t){
    mat sol = mat(Nx,Nt);
    double pi = 4*atan(1.0);

    for(i=0; i<Nt; i++){
        u = 1-x;
        for(n=1; n<cutoff; n++){
            u += -2*sin(n*pi*x)*exp(-n*n*pi*pi*t(i))/(pi*n);
        }
        sol.col(i) = u;
    }

    sol.save("u_analytic.dat", raw_ascii);

}


void DiffusionSolvers::UniStep(vec x){

    // x is a (column) vector with the coordinates of all the particles

    double step   = sqrt(2*dt);
    double x0_lim = 1e-15;
    double rand;
    long seed   = -1;
    int N;
    int new_particles;

    for (t=0; t<=Nt; t++){
        N = x.n_rows; // Total number of particles will deviate
        new_particles = 0;

        for (i=0; i<N; i++){
            rand = ran0(&seed); // Movement right or left will be choosen randomly with equal probability

            if (x(i) < x0_lim){  // Always maintaining the same number of particles at x=0
                new_particles++;
            }
            if (rand>0.5){      // 50/50: Particle number 'i' jumps to the right
                x(i) += step;}


            else {              // If not right, then particle 'i' have to jump left
                x(i) -= step;
            }

            if (abs(x(i)) < x0_lim) {    // Particle 'i' has returned 'home'. Must remove the addition of a new one.
                new_particles--;
            }
        }

        j = 0;
        while (j<N){ // Loop over all current particles
                     // i.e. if a particle is removed, the next one inherits its index
            if (x(j) > 1){
                x.shed_row(j); // Remove particles that have 'made it across'
                N--;        // One less particle in total
            }

            else {
                if (x(j) < 0){
                    x.shed_row(j); // Remove particles that have moved to negative x-values
                    N--;          // One less particle in total
                }

                else {
                    j++;           // Move on to next index if no particles are to be 'deleted'
                }
            }
        }
        // Add new particles all with coordinate x=0
        x.insert_rows(0,new_particles); // fills with zeros by default
        }
    // Save final distribution of particles
    x.save("UniRandWalk.dat", raw_ascii);
}




void DiffusionSolvers::GaussStep(vec x, double bin0){

    // x is a (column) vector with the coordinates of all the particles

    double step;
    double randR, randT;
    long seed   = -1;
    int N;
    int new_particles;

    for (t=0; t<=Nt; t++){
        N = x.n_rows; // Total number of particles will deviate
        new_particles = 0;

        for (i=0; i<N; i++){
            randR = sqrt(-log(1-ran0(&seed)));
            randT = 2*pi*ran0(&seed);

            // The step length will now be chosen randomly from a gaussian distribution:
            step       = sqrt(2*dt)*randR*cos(randT);
            //cout << "particle is located at: " << x(i) << endl;

            if (x(i) < bin0){  // Always maintaining the same number of particles in the first bin
                new_particles++;
                //cout << "added a particle" << endl;
            }

            x(i) += step; // Step is a number between (-inf,inf), actually between +-6, see report

            if (x(i) < bin0 && x(i) > 0){   // If the particle stayed in or returned to bin0
                new_particles--;            // ..I need to 'add one less' new particle.
            }
        }

        j = 0;
        while (j<N){ // Loop over all current particles i.e. if a particle is removed,
                     // the next one inherits its index so j don't have to advance one step
            if (x(j) > 1){
                x.shed_row(j); // Remove particles that have 'made it across' the cleft
                N--;           // One less particle in total (affects the while-loop)
            }

            else {
                if (x(j) < 0){
                    x.shed_row(j); // Remove particles that have moved to negative x-values
                    N--;           // One less particle in total
                }

                else {
                    j++;           // Move on to next index, if the particle are not to be 'deleted'
                }
            }
        }
        // Add new particles all with coordinate x=0
        x.insert_rows(0,new_particles); // fills with zeros by default

        /*cout << "New particles: " << new_particles << endl;
        cout << "Total number of particles: " << x.n_rows << endl;
        _sleep(200);*/

        }
    // Save final distribution of particles
    x.save("GaussRandWalk.dat", raw_ascii);

}
