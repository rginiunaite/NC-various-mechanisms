//
// Created by giniunaite on 03/04/19.
//

//
// Very simple implementation of Lennard-Jones model
//




#include "Aboria.h"
#include <random>
#include <map>
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
//int  main(){


// parameters

// model parameters
// model parameters

/*
 * Initial length of the domain is 300 \mu m, I will create a matrix of chemoattractant which will be of length 300/10= 30.
 * All the length parameters will be scaled consistently, i.e. /10.
 */
    bool first_part_grows = true; //false if final part grows faster
    double space_grid_controller = 100.0;

    double domain_length = 3.42; //this variable is for the actual domain length, since it will be increasing
    double Lt_old = domain_length;
    int length_x =
            int(domain_length*space_grid_controller); // length in x direction of the chemoattractant matrix
    double initial_domain_length = domain_length;
    int real_length_y = 120;
    const int length_y = int(1.2 * space_grid_controller); // length in y direction of the chemoattractant matrix
    const double final_time = 54; // number of timesteps, 1min - 1timestep, from 6h tp 24hours.
    double final_length = 1014;//real 1014

// parameters for the dynamics of chemoattractant concentration

    double D = 0.0002;//0.00001;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0.0; // initialise time
    double dt = 0.01; // time step
    double dt_init = dt;
    int number_time = int(1 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1.0; // space step in x direction, double to be consistent with other types
    double dy = 1.0; // space step in y direction
    //double kai = 0.00001;//0;//0.1 // to 1 /h production rate of chemoattractant
//    double y_init = 120.0;

    // reaction rate
    double k_reac = 0;//0.00001;//0.1;//0.105;//0.03;//.205; // reaction term


    // cell parameters

    double cell_radius = 7.5;//0.5; // radius of a cell
    const double diameter =
            2 * cell_radius; // diameter of a cell
    const size_t N = 5; // initial number of cells
    double l_filo_y = 27.5;//2; // sensing radius, filopodia + cell radius
    double l_filo_x = 27.5; // sensing radius, it will have to be rescaled when domain grows
    double l_filo_x_in = l_filo_x; // this value is used for rescaling when domain grows based on initial value
    double l_filo_max = 45; // this is the length when two cells which were previously in a chain become dettached
    //double diff_conc = 0.1; // sensing threshold, i.e. how much concentration has to be bigger, so that the cell moves in that direction
    int freq_growth = 1; // determines how frequently domain grows (actually not relevant because it will go every timestep)
    int insertion_freq = 1; // determines how frequently new cells are inserted, regulates the density of population
    double speed_l = 70.0/60.0;//10 * 0.05;// 0.05;//1;//0.05; // speed of a leader cell
    double increase_fol_speed = 1.3;
    double speed_f = increase_fol_speed * speed_l;//0.05;//0.1;//0.08; // speed of a follower cell
    double dettach_prob = 0.5; // probability that a follower cell which is on trail looses the trail
    double chemo_leader = 0.9; //0.5; // phenotypic switching happens when the concentration of chemoattractant is higher than this (presentation video 0.95), no phenotypic switching
    double eps = 1; // for phenotypic switching, the distance has to be that much higher
    const int filo_number = 1; // number of filopodia sent
    int same_dir = 0; // number of steps in the same direction +1, because if 0, then only one step in the same direction
    bool random_pers = true; // persistent movement also when the cell moves randomly
    int count_dir = 0; // this is to count the number of times the cell moved the same direction, up to same_dir for each cell
    double lam = 1.0;//72/(100)/10; // to 1000 /h chemoattractant internalisation
    int value = 0; // value of the Gamma(value), where Gamma is close to a cell center

    double magn = 2.0; // magnitude of CiL effect
    double searchCIL = 1.5*diameter; // CiL occurs when the centres of two cells are within this distance
    vdouble2 sum_forces; // some of all the forces exerted on one cell, always set to zero for a new cell
    vdouble2 bound_sum_forces; // some of all the forces exerted on one cell, always set to zero for a new cell
    double npow = 2.0; // Lennard-Jones type model powers
    double mpow = 1.0*npow; //  Lennard-Jones type model powers

    // distance to the track parameters
    double dist_thres = 1;
    int closest_time;
    int leader_track;

    //int n_seed = 0;
    double diff_conc = 0.05;


    // domain growth parameters, sode in domain growth, domain_length

//    double L_inf = 867.6;
//    double L0 = 300;
//    double t_s = 12.77;
//    double alpha = 0.288;
//    double k_0 = 291.2;


    //parameters for Lennard-Jones model
    double sigma_max = 15.0; // maximum distance at which the cells can have an effect on each other
    double search_param = 50.0; //radius to search for cell-cell interactions
    double f0 = 1.0;// strength of the force
    double eps_ij = 10.0; // depth of potential well
    vdouble2 force_ij; // store the value of force
    double dzdr = 0.0; // Lennard-Jones potential
    vdouble2 change; // difference between the positions of two cells
    double distance = 0.0; // distance between two cells
    vdouble2 random_vector;


    
    /*
      * strain rate
      * */


    //piecewise constant, two parts
    // 1 part n_faster times faster than the other part
    double n_faster = 2.0;

    double thetasmall = 1.0; // first thetasmall is growing
    int theta1 = int(thetasmall * length_x);

    double alpha1;
    double alpha2;

    //check whether first part grows
    double thetasmalltemp = thetasmall;
    if (first_part_grows ==true){
        double xvar = final_length/ (n_faster * double(length_x)*thetasmalltemp + double(length_x)*(1-thetasmalltemp)); // solve: 2 *xvar * length_x * thetasmall + x * length(1-thetasmall) = final_length

        double ratio1 = n_faster * double(length_x)*thetasmalltemp * xvar / (double(length_x)*thetasmalltemp);

        alpha1 = log (ratio1)/final_time;

        double ratio2 =  double(length_x)*(1- thetasmalltemp) * xvar / (double(length_x)*(1- thetasmalltemp));

        alpha2 = log (ratio2)/final_time;
    }

    else if(first_part_grows == false){
        double xvar = final_length/ (double(length_x)*thetasmall + n_faster*double(length_x)*(1-thetasmall)); // solve: 2 *xvar * length_x * thetasmall + x * length(1-thetasmall) = final_length

        double ratio1 = double(length_x)*thetasmall * xvar / (double(length_x)*thetasmall);

        alpha1 = log (ratio1)/final_time;

        double ratio2 =  double(length_x)*(1- thetasmall) * xvar*n_faster / (double(length_x)*(1- thetasmall));

        alpha2 = log (ratio2)/final_time;
    }




    VectorXd strain = VectorXd::Zero(length_x);


    // constant
//    for (int i = 0; i < length_x; i++) {
//        strain(i) = alpha;// 0;//linear_par * double(theta1) / double(space_grid_controller);//0;
//    }


    // first part it is linear growth
    for (int i = 0; i < theta1; i++) {
        strain(i) = alpha1;//linear_par * double(theta1) /
        // double(space_grid_controller);//linear_par * (double(i) / double(space_grid_controller));
    }


    // second part is constant
    for (int i = theta1; i < length_x; i++) {
        strain(i) = alpha2;// 0.002;//0.5 * linear_par * double(theta1) / double(space_grid_controller); // constant to where it was
        //strain(i,j) = linear_par*theta1/(theta1- (theta2-1))*(i-(theta2-1)); // linearly decreasing, if I do this do not forget to change Gamma
    }



    // growth function

    // I will mainly use its first derivative with respect to x

    VectorXd Gamma_x = VectorXd::Zero(length_x);
    VectorXd Gamma = VectorXd::Zero(length_x);
    VectorXd Gamma_t = VectorXd::Zero(length_x);
    VectorXd Gamma_old = VectorXd::Zero(length_x);

    for (int i = 0; i < length_x; i++) {
        Gamma_x(i) = exp(0 * strain(i));
        //Gamma(i) = Gamma_initial * i / space_grid_controller;
        Gamma(i) = i;
        Gamma_old(i) = Gamma(i);
    }

    double Gamma_initial = Gamma(length_x - 1);


    // for total length
    double Lt = 0;
    double Ltdot = 0;

    /*
    * initialise a matrix that stores values of concentration of chemoattractant
    */

    MatrixXd chemo = MatrixXd::Zero(length_x, length_y);
    MatrixXd chemo_new = MatrixXd::Zero(length_x, length_y);

//    // uniform initial conditions
//    for (int i = 0; i < length_x; i++) {
//        for (int j = 0; j < length_y; j++) {
//            chemo(i, j) = 1; // uniform concentration initially
//            chemo_new(i, j) = 1; // this is for later updates
//        }
//    }e

    // non uniform initial conditions
    double beta = 1.0; // up to here the initial chemo concentration is C_0
    double C0 = 1.0; // initially non-zero, afterwards zero
    double n = 10.0; // for 1 - 0.5 cos(n \pi x)

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 1;//cos(
            // M_PI * i / space_grid_controller);//1;//C0 - 0.5 * cos( M_PI * i/space_grid_controller * n);
        }
    }

    // initialise internalisation matrix
    MatrixXd intern = MatrixXd::Zero(length_x, length_y);


    // four columns for x, y, z, u (z is necessary for paraview)

    // form a matrix which would store x,y,z,u

    MatrixXd chemo_3col(length_x * length_y, 4), chemo_3col_ind(length_x * length_y,
                                                                2); // need for because that is how paraview accepts data, third dimension is just zeros

    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k < length_x * length_y) {
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                chemo_3col_ind(k, 0) = i;
                chemo_3col_ind(k, 1) = j;
                chemo_3col(k, 2) = 0;
                k += 1;
            }
        }
    }

    // save the x coordinates, scaling only based on the grid
    for (int i = 0; i < length_x * length_y; i++) {
//        chemo_3col(i, 0) = Gamma_initial * chemo_3col_ind(i, 0) / double(space_grid_controller);
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);// / double(space_grid_controller);

    }



    // u column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
    }


    // y coordinates, 1D so nothing changes
    for (int i = 0; i < length_x * length_y; i++) {
        // chemo_3col(i, 1) = y_init * chemo_3col_ind(i, 1) / double(space_grid_controller);
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);// / double(space_grid_controller);
    }


    int counter = 0;






    /*
     * 2D domain with a few randomly placed particles
     */

    /*
     * initial cells of fixed radius
     */

    //ABORIA_VARIABLE(velocity,vdouble2,"velocity")
    ABORIA_VARIABLE(radius, double, "radius");
    ABORIA_VARIABLE(direction, vdouble2, "direction");// stores the direction a particle moved
    ABORIA_VARIABLE(persistence_extent, int,
                    "persistence extent");// stores whether cell moves only one step in current direction or in a process of moving persistently
    ABORIA_VARIABLE(same_dir_step, int,
                    "same dir step");// the number which stores how many steps in the same direction are made.
    ABORIA_VARIABLE(attached_to_id, int, "attached_to_id");
    ABORIA_VARIABLE(type, int, "type");// 0 if a cell is a leader, 1 if follower
    ABORIA_VARIABLE(chain_type, int, "chain_type"); // leaders form different chain types
    ABORIA_VARIABLE(chain, int,"chain"); // stores whether a follower is part of the chain or no, 0 if it is not part of
    ABORIA_VARIABLE(scaling, int, "scaling ") // stores the value that the cell position is scaled down to initial
    // the chain and then increasing integer if it is. If it is attached to a leader, it is 1, and then increasing order.
    // stores the distance to the closest neighbour, if less than thresold
    typedef Particles<std::tuple<radius, type, attached_to_id, direction, chain, chain_type, persistence_extent, same_dir_step,scaling>, 2> particle_type; // 2 stands for dimension

    // will use stored value of the position of a particle
    typedef particle_type::position position;

    // initialise the number of particles
    particle_type particles(N);

    // initialise random number generator for particles entering the domain, appearing at the start in x and uniformly in y
    std::default_random_engine gen;
    std::uniform_real_distribution<double> uniform(cell_radius , length_y - 1 - cell_radius);// when length is 120

    //for double length (240)
    //std::uniform_real_distribution<double> uniform(cell_radius+ double(real_length_y)/2, length_y - double(real_length_y)/2 -1 - cell_radius);
    // for 180 length
    //std::uniform_real_distribution<double> uniform(cell_radius+ double(real_length_y)/4, length_y - double(real_length_y)/4 -1 - cell_radius);



    /*
     * compact initialisation of particles
     */

    for (int i = 0; i < N; ++i) {


        get<radius>(particles[i]) = cell_radius;
        get<type>(particles[i]) = 0; // initially all cells are leaders

        //get<position>(p) = vdouble2(cell_radius,(i+1)*diameter); // x=2, uniformly in y

        // one line of cells
        get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=radius, uniformly in y


        // two lines of cells
//        if (i < N/2){
//            get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N/2) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N/2)); // x=radius, uniformly in y
//        }
//        else{
//            get<position>(particles[i]) = vdouble2(3*cell_radius, (i-N/2 + 1) * double(length_y - 1) / double(N/2) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N/2)); // x=radius, uniformly in y
//        }


        // for attraction
        // get<position>(particles[i]) = vdouble2(cell_radius, 50 +i*30); // x=radius, uniformly in y


        // for double length (240)
//        get<position>(particles[i]) = vdouble2(cell_radius, real_length_y/2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                            0.5 * double(real_length_y - 1) /
//                                                            double(N)); // x=radius, uniformly in y

        // for length 180
//        get<position>(particles[i]) = vdouble2(cell_radius, real_length_y/4 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                            0.25 * double(real_length_y - 1) /
//                                                            double(N)); // x=radius, uniformly in y


        get<persistence_extent>(particles[i]) = 0;
        get<same_dir_step>(particles[i]) = 0;

    }




    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(-100,-100),5* vdouble2(length_x, length_y), vbool2(false, false));

    // save particles before they move

    vtkWriteGrid("CILParticles", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(t * n_seed); // choose different seeds to obtain different random numbers
    //std::uniform_real_distribution<double> uniformpi(-M_PI/3.5,M_PI/3.5 +  M_PI); // 0 to up, M_PI/2 move straigth M_PI downt

    //std::normal_distribution<double> normal(M_PI/2,2); // mean 0 variance 1

    // uniformly random
    //std::uniform_real_distribution<double> uniformpi(0,2*M_PI); // 0 to up, M_PI/2 move straigth M_PI downt

    // for Lennard-Jones
    std::normal_distribution<double> normalX(0.5,2); // mean 1 variance 1
    std::normal_distribution<double> normalY(0,2); // mean 0 variance 1




    //for each timestep
    while (t < final_time) {



//      insert new cells


        bool free_position = false;
        particle_type::value_type f;
        //get<radius>(f) = cell_radius;


        get<position>(f) = vdouble2( cell_radius, uniform(gen)); // x=2, uniformly in y
        free_position = true;
        /*
         * loop over all neighbouring leaders within "dem_diameter" distance
         */
        for (auto tpl = euclidean_search(particles.get_query(), get<position>(f), 1.5*diameter); tpl != false; ++tpl) {

            vdouble2 diffx = tpl.dx();

            if (diffx.norm() < diameter) {
                free_position = false;
                break;
            }
        }

        // all the cells are of the same type
        get<type>(f) = 0;


        if (free_position) {
            get<chain>(f) = 0;
            get<chain_type>(f) = -1;
            get<attached_to_id>(f) = -1;
            particles.push_back(f);
        }


        particles.update_positions();


        t = t + dt;

        counter = counter + 1;



        /*
         * Domain growth from here
        * */


        for (int i = 0; i < length_x; i++) {
            Gamma_x(i) = exp(t * strain(i));
        }


        Gamma(0) = 0; // this is assumption, since I cannot calculate it
        //cout << "Gamma(0) " << 0 << " value " << Gamma(0) << endl;

        for (int i = 1; i < length_x;i++) {

            Gamma(i) = Gamma_x(i) * dx + Gamma(i - 1);

        }


        for (int i = 0; i < length_x; ++i) {

            Gamma_t(i) = (Gamma(i) - Gamma_old(i)) / dt;

        }


//
        /// update positions uniformly based on the domain growth

        vdouble2 x; // use variable x for the position of cells
        double x0=0;
        int pos;

        for (int i = 0; i < particles.size(); i++) {

            x = get<position>(particles[i]);

            // since I do not know how to do it for general case, I will do it for my specific

            // analytical for two different regions
//            if (x[0] > Gamma_old(theta1 - 1)) {
//
//                x0 = (x[0] - Gamma_old(theta1 - 1)) *
//                        ((Gamma(length_x - 1) - Gamma(theta1 - 1)) /
//                     (Gamma_old(length_x - 1) - Gamma_old(theta1 - 1))) + Gamma(theta1 - 1);
//
//                get<position>(particles)[i] = vdouble2(x0,x[1]);
//
//            } else {
//                get<position>(particles)[i] *= vdouble2((Gamma(theta1 - 1)) / (Gamma_old(theta1 - 1)),
//                                                        1); // update position based on changes in Gamma
//            }
            int j = 0;
            while (x[0]> Gamma_old(j)){
                value = j;
                j = j+1;
                //cout << "value " << value << endl;

            }


            get<scaling>(particles)[i] = value;

//            x[0] = x[0] + Gamma(value)-Gamma_old(value);


            get<position>(particles)[i] += vdouble2(Gamma(value)-Gamma_old(value),0);



        }


        Gamma_old = Gamma;
//
//
//
//
//
//
//
//
//
        // save the chemoattractant concentration with properly rescaled coordinates

        int counting_first = 0;
        int counting_final = 0;

        for (int a = 0; a < length_x; a++) {
            counting_first = length_y * a;
            counting_final = counting_first + length_y;
            for (int k = counting_first; k < counting_final; k++) {
                chemo_3col(k, 0) = Gamma(a);
            }
        }


        // u column
        for (int i = 0; i < length_x * length_y; i++) {
            chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
        }



        /*
 * Comment up to here without domain growth
 * */




        /*
         * Update the position of particles
         * */


        //  create a random list of cell ids
        int check_rep = 0; // check for repetitions, 0 no rep, 1 rep


        std::default_random_engine gen2;
        gen2.seed(t * n_seed); // different seeds
        std::uniform_real_distribution<double> uniform_particles(0, particles.size()); // can only move forward

        VectorXi particle_id = VectorXi::Zero(particles.size());

        for (int i = 0; i < particles.size(); i++) {

            check_rep = 1; // set to 1 to enter the while loop
            while (check_rep == 1) {
                check_rep = 0; // it is initially zero and then will be changed to 1 if it is equivalent to others
                particle_id(i) = uniform_particles(gen2);


                for (int j = 0; j < i; j++) {
                    if (particle_id(i) == particle_id(j)) { check_rep = 1; }
                }
            }
        }



        for (int j = 0; j < particles.size(); j++) {



            vdouble2 x; // use variable x for the position of cells

            x = get<position>(particles[particle_id(j)]);
            // update the position of a cell based on the deterministic and random forces exerted on it




            // deterministic force, from all the cells within the threshold distance


            sum_forces = vdouble2(0, 0);
            bound_sum_forces = vdouble2(0, 0);
            force_ij = vdouble2(0, 0);

            for (auto k = euclidean_search(particles.get_query(), x, search_param); k != false; ++k) {

                // make sure it is not the same cell
                if (get<id>(*k) != get<id>(particles[particle_id(j)])) {
                    change = get<position>(particles[particle_id(j)]) - get<position>(*k);
                    distance = change.norm();

                    dzdr = npow * eps_ij*(2* pow((sigma_max / distance),mpow) - pow((sigma_max / distance),npow)); // Lennard-Jones potential

                    force_ij = f0 * dzdr * change / distance;
//                    cout << "change " << change << endl;
//                    cout << "particle id " << get<id>(particles[particle_id(j)]) << endl;
//                    cout << "distance " << distance << endl;
//                    cout << "specific force " << force_ij << endl;
                }

                sum_forces += force_ij;

            }

            if(t>120 & t <150){
                cout << "sum forces " << sum_forces << endl;

            }

            //random force

            random_vector[0] = normalX(gen1);
            random_vector[1] = normalY(gen1);

            x = get<position>(particles[particle_id(j)]) + dt * (sum_forces+ random_vector);// + bound_sum_forces);


            // boundary forces before the position is updated

            if(x[0]<cell_radius){
                bound_sum_forces += (cell_radius - x[0])/x[0] * vdouble2(x[0],0);
            }
            if(x[1]<cell_radius){
                bound_sum_forces += (cell_radius - x[1])/x[1] * vdouble2(0,x[1]);

            }
            if(x[1]> length_y - cell_radius -1){
                bound_sum_forces += (length_y - cell_radius - 1 - x[1])/x[1] * vdouble2(0,x[1]);

            }

            if(t>120 & t <150){
                cout << "bound sum forces " << bound_sum_forces << endl;

            }

            x = x + (bound_sum_forces);



            //if (x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
            //    (x[1]) < length_y - 1 - cell_radius) {
            get<position>(particles[particle_id(j)]) = x;

            //}

            particles.update_positions(); // not sure if needed here

        }

        //particles.update_positions(); // not sure if needed here


        if (counter % 100 == 0) {


            // save at every time step
#ifdef HAVE_VTK
            vtkWriteGrid("BIASEDNOVELCILParticles", t, particles.get_grid(true));
#endif




            ofstream output("BIASEDNOVELCILMatrix" + to_string(int(t)) + ".csv");


            output << "x, y, z, u" << "\n" << endl;


            for (int i = 0; i < length_x * length_y; i++) {
                for (int j = 0; j < 4; j++) {
                    output << chemo_3col(i, j) << ", ";
                }
                output << "\n" << endl;
            }



        }



    }


    //    /*
// * return the density of cells in domain_partition parts of the domain
// */
    const int domain_partition = round(Gamma(length_x-1) / double(55));; // number of intervals of 50 \mu m



    VectorXi proportions = VectorXi::Zero(domain_partition); // integer with number of cells in particular part

    double one_part = Gamma(length_x-1) / double(domain_partition);


    for (int i = 0; i < domain_partition; i++) {

        for (int j = 0; j < particles.size(); j++) {
            vdouble2 x = get<position>(particles[j]);
            if (i * one_part < x[0] && x[0] < (i + 1) * one_part) {
                proportions(i) += 1;
            }
        }

    }

    return proportions;

}


// parameter analysis
int main() {

    const int number_parameters = 1; // parameter range
    const int sim_num = 1;

//    VectorXi vector_check_length = proportions(2); //just to know what the length is
//
//    int num_parts = vector_check_length.size(); // number of parts that I partition my domain
//    cout << "length " << vector_check_length.size() << endl;
    int num_parts = 20; // for 1800 timesteps
    MatrixXf sum_of_all = MatrixXf::Zero(num_parts, number_parameters); // sum of the values over all simulations

    // n would correspond to different seeds
    // parallel programming
#pragma omp parallel for
    for (int n = 0; n < sim_num; n++) {


        //initialise the matrix to store the values
        MatrixXi numbers = MatrixXi::Zero(num_parts, number_parameters);

        //#pragma omp parallel for
        //        for (int i = 0; i < number_parameters; i++) {

        // for (int j = 0; j < 1; j++) {

        numbers.block(0, 0, num_parts, 1) = proportions( n);

        //}
        // }


//        // This is what I am using for MATLAB
//        ofstream output2("growing-width180-normal2var" + to_string(n) + ".csv");
//
//        for (int i = 0; i < numbers.rows(); i++) {
//
//            for (int j = 0; j < numbers.cols(); j++) {
//
//                output2 << numbers(i, j) << ", ";
//
//              //  sum_of_all(i, j) += numbers(i, j);
//
//            }
//            output2 << "\n" << endl;
//        }

    }
    /*
    * will store everything in one matrix, the entries will be summed over all simulations
    */

//    //ofstream output3("FIRST075_twice_speed_later.csv");
//    ofstream output3("width180-uni.csv");
//
//
//    for (int i = 0; i < num_parts; i++) {
//
//        for (int j = 0; j < number_parameters; j++) {
//
//            output3 << sum_of_all(i, j) << ", ";
//
//        }
//        output3 << "\n" << endl;
//    }


}



