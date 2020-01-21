//
// Created by giniunaite on 14/01/20.
// this code will be used for CiL and CoA using Lennard Jones potential. Langevin model will be used
//



#include "Aboria.h"
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXd, MatrixXd


VectorXi proportions(int n_seed) {

    bool domain_growth = true; // if false change length_x to 1100, true 300
    bool CiLonly = false;

    int length_x;
    if (domain_growth == false) {
        length_x = 1100;
    } else {
        length_x = 300;// length in x direction of the chemoattractant matrix
    }

    int length_y = 120;
    double Lt_old = length_x;
    int real_length_y = 120;
    const double final_time = 1080; // 18hrs number of timesteps, 1min - 1timestep, from 6h tp 24hours. 1440 if 24hrs
    double t = 0.0; // initialise time
    double dt = 0.05; // time step
    double dx = 1; // maybe 1/100
    int counter = 0; // to count simulations
    const size_t N = 5; // initial number of cells
    double sigma = 2.0;
    double mean = 0.0; // mean movement in x direction
    double cell_radius = 7.5;//0.5; // radius of a cell
    const double diameter =
            2 * cell_radius; // diameter of a cell

    double D =20.0; // diffusion coefficient for Brownian motion

    // CiL and CoA paramters, Lennard Jones
    double searchCIL = 1.0*diameter; // CiL occurs when the centres of two cells are within this distance
    vdouble2 sum_forces; // some of all the forces exerted on one cell, always set to zero for a new cell
    vdouble2 bound_sum_forces; // some of all the forces exerted on one cell, always set to zero for a new cell
    double npow = 1.0; // Lennard-Jones type model powers, attractive
    double mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive

    //parameters for Lennard-Jones model
    double sigma_max = 15.0;//15.0; // cell diamter ( don't know why I wrote this before: maximum distance at which the cells can have an effect on each other)
    double search_param = 50.0; //radius to search for cell-cell interactions
    double f0 = 1.0;//1 before strength of the force
    double eps_ij = 10.0; // depth of potential well
    vdouble2 force_ij; // store the value of force
    double dzdr = 0.0; // Lennard-Jones potential initialise
    vdouble2 change; // difference between the positions of two cells
    double distance = 0.0; // distance between two cells
    vdouble2 random_vector;


    double furthestCell;

// cell variable
    ABORIA_VARIABLE(radius, double, "radius");
    ABORIA_VARIABLE(direction, vdouble2, "direction");// stores the direction a particle moved
    typedef Particles<std::tuple<radius, direction>, 2> particle_type; // 2 stands for dimension

// will use stored value of the position of a particle
    typedef particle_type::position position;


/*
 * Domain growth section
 * */
// domain growth parameters, domain_length

    double L_inf = 867.6;
    double L0 = 300;
    double t_s = 12.77;
    double alpha = 0.288;
    double k_0 = 291.2;

    VectorXd Gamma = VectorXd::Zero(length_x);
    VectorXd Gamma_old = VectorXd::Zero(length_x);

    for (int i = 0; i < length_x; i++) {
        Gamma(i) = i;
        Gamma_old(i) = Gamma(i);
    }

    // initialise the number of particles
    particle_type particles(N);

    // initialise random number generator for particles entering the domain, appearing at the start in x and uniformly in y
    std::default_random_engine gen;
    // domain width 120
    std::uniform_real_distribution<double> uniform(cell_radius, length_y - 1 - cell_radius);// when length is 120

    //for double width (240)
    //std::uniform_real_distribution<double> uniform(cell_radius+ double(real_length_y)/2, length_y - double(real_length_y)/2 -1 - cell_radius);

    // for domain width 60
    //std::uniform_real_distribution<double> uniform(cell_radius, length_y - 1 - cell_radius);// when length is 120

    /*
* compact initialisation of particles
*/

    for (int i = 0; i < N; ++i) {


        get<radius>(particles[i]) = cell_radius;

        // for domain width 120

            get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
                                                                0.5 * double(length_y - 1) /
                                                                double(N)); // x=radius, uniformly in


//            get<position>(particles[i+5]) = vdouble2(3*cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y
//
//
//            get<position>(particles[i+10]) = vdouble2(5*cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y
//
//
//            get<position>(particles[i+15]) = vdouble2(7*cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                            0.5 * double(length_y - 1) /
//                                                            double(N)); // x=radius, uniformly in y
//
//            get<position>(particles[i+20]) = vdouble2(9*cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                                 0.5 * double(length_y - 1) /
//                                                                 double(N)); // x=radius, uniformly in y






////             // for domain width 240
//            get<position>(particles[i]) = vdouble2(cell_radius,
//                                                   real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                   0.5 * double(real_length_y - 1) /
//                                                   double(N)); // x=radius, uniformly in y


//            get<position>(particles[i+5]) = vdouble2(3*cell_radius, real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                     0.5 * double(real_length_y - 1) /
//                                                     double(N)); // x=radius, uniformly in y
//
//
//            get<position>(particles[i+10]) = vdouble2(5*cell_radius,
//                                                      real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                      0.5 * double(real_length_y - 1) /
//                                                      double(N)); // x=radius, uniformly in y
//
//
//            get<position>(particles[i+15]) = vdouble2(7*cell_radius,
//                                                      real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                      0.5 * double(real_length_y - 1) /
//                                                      double(N)); // x=radius, uniformly in y
//
//            get<position>(particles[i+20]) = vdouble2(9*cell_radius,
//                                                      real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                      0.5 * double(real_length_y - 1) /
//                                                      double(N)); // x=radius, uniformly in y




//        // for domain width 60
//           get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                            0.5 * double(length_y - 1) /
//                                                            double(N)); // x=radius, uniformly in y

    }


    particles.update_positions();


    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));


    vtkWriteGrid("particles", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed( n_seed); // choose different seeds to obtain different random numbers
    //std::uniform_real_distribution<double> uniformpi(-M_PI/3.5,M_PI/3.5 +  M_PI); // 0 to up, M_PI/2 move straigth M_PI downt

    //std::normal_distribution<double> normal(M_PI/2,sigma); // normal distribution for filopodia

    // for Lennard-Jones
    std::normal_distribution<double> normalX(mean,1.0); // mean 1 variance 1
    std::normal_distribution<double> normalY(0.0,1.0); // mean 0 variance 1

    int countfalse = 0;
    int counttrue = 0;

    //for each timestep
    while (t < final_time) {
        //while (furthestCell < 1000.0) {


        //      insert new cells
        //if (particles.size()<25) {
        //if (counter % 100 == 0){
        bool free_position = true;
        particle_type::value_type f;

        get<position>(f) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y

        /*
         * loop over all neighbouring leaders within "dem_diameter" distance
         */
        for (auto tpl = euclidean_search(particles.get_query(), get<position>(f), diameter); tpl != false; ++tpl) {

            vdouble2 diffx = tpl.dx();

            if (diffx.norm() < diameter) {
                free_position = false;
                break;
            }
        }

        // all the cells are of the same type


        if (free_position) {

            particles.push_back(f);
        }


        particles.update_positions();
        //}
        // end of insert new cells

        t = t + dt;

        counter = counter + 1;
/*
      * Domain growth
      * */

        /*
    * Domain growth from here
   * */

        if (domain_growth == true) {
            //cout << "here" << endl;
            /*
             * Piecewise constant // all linear, for presentation
             * */

            Gamma(length_x - 1) = (L_inf * exp(alpha * (24.0 / final_time * t - t_s)) /
                                   (L_inf / L0 + exp(alpha * (24.0 / final_time * t - t_s)) - 1)) +
                                  k_0;

            for (int i = 0; i < length_x - 1; i++) {
                Gamma(i) = (double(i) / (length_x - 1)) * Gamma(length_x - 1);
            }



            /*
             * Domain growth
             * */


// comment this if domain does not grow
            /// update positions uniformly based on the domain growth

            vdouble2 x; // use variable x for the position of cells
            int pos;


            for (int i = 0; i < particles.size(); i++) {

                x = get<position>(particles[i]);
                // since I do not know how to do it for general case, I will do it for my specific

//            if (x[0] > Gamma(length_x - 1)) { // these are very extreme cases
//                get<position>(particles)[i] += vdouble2(Gamma(length_x - 1) - Gamma_old(length_x - 1), 0);
//            } else {

                get<position>(particles)[i] *= vdouble2((Gamma(length_x - 1)) / (Gamma_old(length_x - 1)),
                                                        1); // update position based on changes in Gamma
                //}

            }
            //cout << "Gamma " << Gamma(length_x -1 ) << endl;
// comment this if domain does not grow, up to here
            Gamma_old = Gamma;

        }
        /*
        * Domain growth to here
       * */

        /*
      * Update the position of particles
      * */


        //  create a random list of cell ids
        int check_rep = 0; // check for repetitions, 0 no rep, 1 rep


//        std::default_random_engine gen2;
//        gen2.seed(t * n_seed); // different seeds
//        std::uniform_real_distribution<double> uniform_particles(0, particles.size()); // can only move forward

//        VectorXi particle_id = VectorXi::Zero(particles.size());
//
//        for (int i = 0; i < particles.size(); i++) {
//
//            check_rep = 1; // set to 1 to enter the while loop
//            while (check_rep == 1) {
//                check_rep = 0; // it is initially zero and then will be changed to 1 if it is equivalent to others
//                particle_id(i) = uniform_particles(gen2);
//
//
//                for (int j = 0; j < i; j++) {
//                    if (particle_id(i) == particle_id(j)) { check_rep = 1; }
//                }
//            }
//        }

        MatrixXd positions = MatrixXd::Zero(particles.size(), 2);

        // position update
        for (int j = 0; j < particles.size(); j++) {
            vdouble2 x; // use variable x for the position of cells
//

            //new version

//            bool free_position = false; // check if the neighbouring position is free
//
//        int count =0;
//
//
//
//            while (free_position == false && count < 5){
//                count = count + 1;
//                x = get<position>(particles[particle_id(j)]);
//
//                // update the position of a cell based on the deterministic and random forces exerted on it
//
//
//                // deterministic force, from all the cells within the threshold distance
//
//
//                sum_forces = vdouble2(0, 0);
//                bound_sum_forces = vdouble2(0, 0);
//                force_ij = vdouble2(0, 0);
//
//                for (auto k = euclidean_search(particles.get_query(), x, search_param); k != false; ++k) {
//
//                    // make sure it is not the same cell
//                    if (get<id>(*k) != get<id>(particles[particle_id(j)])) {
//                        change = get<position>(particles[particle_id(j)]) - get<position>(*k);
//                        distance = change.norm();
//
//                        dzdr = npow * eps_ij*(2* pow((sigma_max / distance),mpow));// - pow((sigma_max / distance),npow)); // Lennard-Jones potential, comment the part with npow if I want to ignore attraction
//
//                        //dzdr = (npow *  pow((sigma_max / distance),npow));// - 2*mpow * pow(sigma_max/distance,mpow)); //WRONG
//
//                        force_ij = f0 * dzdr * change / distance;
//                    }
//
//                    sum_forces += force_ij;
//
//                }
//
//                //random force
//
//                random_vector[0] = normalX(gen1);
//                random_vector[1] = normalY(gen1);
//
//                x = get<position>(particles[particle_id(j)]) + dt * (sum_forces+ random_vector);// + bound_sum_forces);
//
//
//                // boundary forces before the position is updated
//
//                if(x[0]<cell_radius){
//                    bound_sum_forces += (cell_radius - x[0])/x[0] * vdouble2(x[0],0);
//                }
//                if(x[1]<cell_radius){
//                    bound_sum_forces += (cell_radius - x[1])/x[1] * vdouble2(0,x[1]);
//
//                }
//                if(x[1]> length_y - cell_radius -1){
//                    bound_sum_forces += (length_y - cell_radius - 1 - x[1])/x[1] * vdouble2(0,x[1]);
//
//                }
//                if(x[0]> Gamma(length_x-1) - cell_radius){
//                    bound_sum_forces += (Gamma(length_x-1) - cell_radius - x[0])/x[0] * vdouble2(x[0],0);
//
//                }
//
//                x = x + (bound_sum_forces);
//
//                free_position = true;
//                // if this loop is entered, it means that there is another cell where I want to move
//                for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {
//
//                    if (get<id>(*k) !=
//                        get<id>(particles[particle_id(j)])) { // check if it is not the same particle
//                        free_position = false;
//                        break;
//                    }
//                }
//
//                // update the position if the place they want to move to is free and not out of bounds
//
//                if (free_position == true){
//                    get<position>(particles[particle_id(j)]) = x;
//                    counttrue = counttrue +1;
//                }
//                //get<position>(particles[particle_id(j)]) = x;
//
//                //cout << direction.norm() << endl; // this is how much a cell moved
//                if (count >1){
//                    //cout << count << endl;
//
//                }
//                if (count == 5){
//                    countfalse = countfalse +1;
//                }
//
//            }

// end of new version




            // old version from here
            x = get<position>(particles[j]);

            // update the position of a cell based on the deterministic and random forces exerted on it


            // deterministic force, from all the cells within the threshold distance


            sum_forces = vdouble2(0, 0);
            bound_sum_forces = vdouble2(0, 0);
            force_ij = vdouble2(0, 0);

            for (auto k = euclidean_search(particles.get_query(), x, search_param); k != false; ++k) {

                // make sure it is not the same cell
                if (get<id>(*k) != get<id>(particles[j])) {
                    change = get<position>(particles[j]) - get<position>(*k);
                    distance = change.norm();

//                    if(get<id>(*k) < N && particle_id(j) >= N){
//                        eps_ij =5;
//                    }
//                    else{
//                        eps_ij =1;
//                    }
                    if (CiLonly == true) {
                        dzdr = npow * eps_ij * (2 * pow((sigma_max / distance), mpow));
                    } else {
                        dzdr = npow * eps_ij * (2 * pow((sigma_max / distance), mpow) - pow((sigma_max / distance),
                                                                                            npow)); // Lennard-Jones potential, comment the part with npow if I want to ignore attraction
                    }


                    force_ij = f0 * dzdr * change / distance;
                }

                sum_forces += force_ij;

            }

            //random force, tryining to make the first part of this biased.
//        if (particle_id(j) < N){
//            random_vector[0] = normalX(gen1)+0.5;
//
//        }
//        else{
//            random_vector[0] = normalX(gen1);
//        }


            random_vector[0] = normalX(gen1);
            random_vector[1] = normalY(gen1);
//
//            cout << "random " << random_vector/random_vector.norm() << endl;
//            cout << "sum forces " << sum_forces << endl;



            x = get<position>(particles[j]) + dt * (sum_forces) + sqrt(2.0 * dt * D) *
                                                                  (random_vector);// + bound_sum_forces); I could have random_vector/random_vector.norm()


            // boundary forces before the position is updated

            if (x[0] < cell_radius) {
                bound_sum_forces += (cell_radius - x[0]) / x[0] * vdouble2(x[0], 0);
            }
            if (x[1] < cell_radius) {
                bound_sum_forces += (cell_radius - x[1]) / x[1] * vdouble2(0, x[1]);

            }
            if (x[1] > length_y - cell_radius - 1) {
                bound_sum_forces += (length_y - cell_radius - 1 - x[1]) / x[1] * vdouble2(0, x[1]);

            }
            if (x[0] > Gamma(length_x - 1) - cell_radius) {
                bound_sum_forces += (Gamma(length_x - 1) - cell_radius - x[0]) / x[0] * vdouble2(x[0], 0);

            }

            x = x + (bound_sum_forces);


            bool free_position = true; // check if the neighbouring position is free

            // if this loop is entered, it means that there is another cell where I want to move
            for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                if (get<id>(*k) !=
                    get<id>(particles[j])) { // check if it is not the same particle
                    free_position = false;
                    break;
                }
            }

            // update the position if the place they want to move to is free and not out of bounds


//            if (freeposition && x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
//                (x[1]) < length_y - 1 - cell_radius){

            if (free_position == true) {
                counttrue = counttrue + 1;
            }
            if (free_position == false) {
                countfalse = countfalse + 1;
            }


            vdouble2 direction;
            direction = (sum_forces + random_vector);


            positions(j, 0) = x[0];
            positions(j, 1) = x[1];


            //if (free_position){
            //get<position>(particles[j]) = x;
            //cout << direction.norm() << endl; // this is how much a cell moved

            //}

            // old version up to here


        }


        for (int ipos = 0; ipos < particles.size(); ipos++) {

            get<position>(particles)[ipos] = vdouble2(positions(ipos, 0), positions(ipos, 1));
        }


        particles.update_positions(); // not sure if needed here

        if (counter % 100 == 0) {
              //if (t > 1078){
            //if (furthestCell > 980) { // for the one to see how long it takes till they reach the end


//                // save at every time step
//            #ifdef HAVE_VTK
//                vtkWriteGrid("TRIALD20eps5MoreCellsCiLParticles", t, particles.get_grid(true));
//            #endif
//

            //}
//

            // for the one to see how long it takes till they reach the end
//        // postion of five cells at the front
//        vector<double> xArray;
//        // positions of cells most at the front
//        for (int j = 0; j < particles.size(); j++) {
//            vdouble2 xi = get<position>(particles[j]);
//            xArray.push_back(xi[0]);
//        }
//        sort(xArray.begin(), xArray.end(),greater<double>()); // sort descending order
//
//        furthestCell = xArray[0];
//        //cout << furthestCell << endl;


        }
        //cout << "Final t " << t << endl;
    }

//    /*
// * return the density of cells in domain_partition parts of the domain
// */
    const int domain_partition = int(round(Gamma(length_x -1 ) / double(55))); // number of intervals of 50 \mu m



    VectorXi proportions = VectorXi::Zero(domain_partition); // integer with number of cells in particular part

    double one_part = Gamma(length_x -1 ) / double(domain_partition);


    for (int i = 0; i < domain_partition; i++) {

        for (int j = 0; j < particles.size(); j++) {
            vdouble2 x = get<position>(particles[j]);
            if (i * one_part < x[0] && x[0] < (i + 1) * one_part) {
                proportions(i) += 1;
            }
        }

    }
    //position of leaders
//    for (int i = 0; i<N; i++){
//        vdouble2 xposi = get<position>(particles[i]);
//        //cout << "position leaders " << get<position>(particles[i]) << endl;
//        //cout << xposi[0] << endl;
//    }

    // postion of five cells at the front
    vector<double> xArray;
    // positions of cells most at the front
    for (int j = 0; j < particles.size(); j++) {
        vdouble2 xi = get<position>(particles[j]);
        xArray.push_back(xi[0]);
    }
    sort(xArray.begin(), xArray.end(),greater<double>()); // sort descending order
    for (int i = 0; i<N; i++){

        cout << xArray[i] << endl;
    }
//
//    cout << "true " << counttrue << endl;
//    cout << "false " << countfalse << endl;

    return proportions;

}



// parameter analysis
int main() {

    const int number_parameters = 1; // parameter range
    const int sim_num = 20;

    VectorXi vector_check_length = proportions(2); //just to know what the length is
//
    int num_parts = vector_check_length.size(); // number of parts that I partition my domain
//    cout << "length " << vector_check_length.size() << endl;
    //   int num_parts = 20; // for 1800 timesteps
    MatrixXf sum_of_all = MatrixXf::Zero(num_parts, number_parameters); // sum of the values over all simulations

//     n would correspond to different seeds
////     parallel programming
#pragma omp parallel for
    for (int n = 1; n < sim_num; n++) {


        //initialise the matrix to store the values
        MatrixXi numbers = MatrixXi::Zero(num_parts, number_parameters);


        numbers.block(0, 0, num_parts, 1) = proportions( n);

        // This is what I am using for MATLAB
        ofstream output2("sepdataD20eps10m2CoACiLGrowingDomain" + to_string(n) + ".csv");

        for (int i = 0; i < numbers.rows(); i++) {

            for (int j = 0; j < numbers.cols(); j++) {

              output2 << numbers(i, j) << ", ";

                sum_of_all(i, j) += numbers(i, j);

            }
            output2 << "\n" << endl;
        }

    }
//    /*
//    * will store everything in one matrix, the entries will be summed over all simulations
//    */

    //ofstream output3("FIRST075_twice_speed_later.csv");
    ofstream output3("DATAD20eps10m2CoACiLGrowingDomain.csv"); // at the end GrowingDomain


    for (int i = 0; i < num_parts; i++) {

        for (int j = 0; j < number_parameters; j++) {

            output3 << sum_of_all(i, j) << ", ";

        }
        output3 << "\n" << endl;
    }


}