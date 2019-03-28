//
// Created by giniunaite on 21/02/19.
//


/*
 * All the cells move randomly, no chemoattractant, in addition, there is CiL, cells not only can not overlap, they also repel each other.
 * */



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

    double space_grid_controller = 100.0;

    double domain_length = 3.0; //this variable is for the actual domain length, since it will be increasing
    double Lt_old = domain_length;
    int length_x =
            int(domain_length) * int(space_grid_controller); // length in x direction of the chemoattractant matrix
    double initial_domain_length = domain_length;
    int real_length_y = 120;
    const int length_y = int(1.2 * space_grid_controller); // length in y direction of the chemoattractant matrix
    const double final_time = 72; // number of timesteps, 1min - 1timestep, from 6h tp 24hours.

// parameters for the dynamics of chemoattractant concentration

    double D = 0.0002;//0.00001;//0.05; // to 10^5 \nu m^2/h diffusion coefficient
    double t = 0.0; // initialise time
    double dt = 0.05; // time step
    double dt_init = dt;
    int number_time = int(1 / dt_init); // how many timesteps in 1min, which is the actual simulation timestep
    double dx = 1.0 / space_grid_controller; // space step in x direction, double to be consistent with other types
    double dy = 1.0 / space_grid_controller; // space step in y direction
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

    double magn = 2.0; // magnitude of CiL effect
    double searchCIL = 1.5*diameter; // CiL occurs when the centres of two cells are within this distance

    // distance to the track parameters
    double dist_thres = 1;
    int closest_time;
    int leader_track;


    //int n_seed = 0;
    double diff_conc = 0.05;


    // domain growth parameters, sode in domain growth, domain_length

    double L_inf = 867.6;
    double L0 = 300;
    double t_s = 12.77;
    double alpha = 0.288;
    double k_0 = 291.2;




    /*
    * strain rate
    * */

    //piecewise constant, two parts
    double thetasmall = 1.00; // first thetasmall is growing
    int theta1 = int(thetasmall * length_x);





    VectorXd strain = VectorXd::Zero(length_x);




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
    // the chain and then increasing integer if it is. If it is attached to a leader, it is 1, and then increasing order.
    // stores the distance to the closest neighbour, if less than thresold
    typedef Particles<std::tuple<radius, type, attached_to_id, direction, chain, chain_type, persistence_extent, same_dir_step>, 2> particle_type; // 2 stands for dimension

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

        get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=radius, uniformly in y



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
    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));

    // save particles before they move

    vtkWriteGrid("particles", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(t * n_seed); // choose different seeds to obtain different random numbers
    //std::uniform_real_distribution<double> uniformpi(-M_PI/3.5,M_PI/3.5 +  M_PI); // 0 to up, M_PI/2 move straigth M_PI downt

    //std::normal_distribution<double> normal(M_PI/2,2); // mean 0 variance 1

    // uniformly random
    std::uniform_real_distribution<double> uniformpi(0,2*M_PI); // 0 to up, M_PI/2 move straigth M_PI downt






    //for each timestep
    while (t < final_time) {



//      insert new cells


        bool free_position = false;
        particle_type::value_type f;
        //get<radius>(f) = cell_radius;


        get<position>(f) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y
        free_position = true;
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


        /*
         * Piecewise constant // all linear, for presentation
         * */

        Gamma(theta1-1) = (L_inf * exp(alpha * (24.0/double(final_time)*t - t_s)) / (L_inf / L0 + exp(alpha * (24.0/double(final_time)*t - t_s)) - 1)) +
                          k_0;

        for (int i = 0; i < theta1-1; i++) {
            Gamma(i) =  double(i)/(theta1-1) * Gamma(theta1-1);
        }





        /*
         * this time I will find Gamma_x from Gamma, before it was other way around, assume theta is the last one
         */


        for (int i = 0; i < theta1 - 1; i++) {
//            Gamma(i) =
//                    Gamma_initial * (double(i) / double(space_grid_controller)) * (Gamma_x(i)); // linearly increasing
            Gamma_x(i) = (Gamma(i + 1) - Gamma(i)) / dx; // linearly increasing

        }

        Gamma_x(theta1 - 1) = Gamma_x(theta1 - 2);


        for (int i = 0; i < theta1; i++) {
//            Gamma(i) =
//                    Gamma_initial * (double(i) / double(space_grid_controller)) * (Gamma_x(i)); // linearly increasing
            strain(i) = log(Gamma_x(i)) / t; // linearly increasing

        }



        /*
         * Domain growth
         * */

        Lt = thetasmall + thetasmall * exp(alpha * t);
        Lt = Gamma(length_x - 1);


        //   Ltdot = alpha * thetasmall * exp(alpha * t); // piecewise, exact

        Ltdot = (Lt - Lt_old) / dt; //could be used for both, especially should be used if derivative is unknown

        Lt_old = Lt;


        // I need Gamma_t for cos verification as well

        for (int i = 0; i < length_x; ++i) {

            Gamma_t(i) = (Gamma(i) - Gamma_old(i)) / dt;

        }



        /// update positions uniformly based on the domain growth

        vdouble2 x; // use variable x for the position of cells
        int pos;

        for (int i = 0; i < particles.size(); i++) {

            x = get<position>(particles[i]);
            // since I do not know how to do it for general case, I will do it for my specific

            if (x[0] > Gamma(theta1 - 1)) {
                get<position>(particles)[i] += vdouble2(Gamma(theta1 - 1) - Gamma_old(theta1 - 1), 0);
            } else {
                get<position>(particles)[i] *= vdouble2((Gamma(theta1 - 1)) / (Gamma_old(theta1 - 1)),
                                                        1); // update position based on changes in Gamma
            }

        }


        Gamma_old = Gamma;









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

            vdouble2 xnbh; // position of the neighbour


            double x_in;

            vdouble2 change; // vector which will incorporate change due to CiL for one cell
            vdouble2 changenbh; // vector which will incorporate change due to CiL for the neighbouring cell


            /*
             * CiL, check if there are close neighbours and all of the bounce of
             * */


            for (auto k = euclidean_search(particles.get_query(), x, searchCIL); k != false; ++k) {

                change = get<position>(particles[particle_id(j)]) - get<position>(*k);
                changenbh = get<position>(*k) - get<position>(particles[particle_id(j)]);

                xnbh =  get<position>(*k) + magn*changenbh/changenbh.norm();

                x += magn * change/change.norm();

                bool free_position = true; // check if the neighbouring position is free

                bool free_position_nbh = true; // check if the neighbouring position is free


                //for the main cell
                // if this loop is entered, it means that there is another cell where I want to move
                for (auto p = euclidean_search(particles.get_query(), x, diameter); p != false; ++p) {

                    if (get<id>(*p) !=
                        get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                        free_position = false;
                    }
                }


                // update the position if the place they want to move to is free and not out of bounds
                if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
                    (x[1]) < length_y - 1 - cell_radius) {
                    get<position>(particles)[particle_id(j)] +=
                            magn * change/change.norm(); // update if nothing is in the next position

                }

                // for the neighbour

                for (auto j = euclidean_search(particles.get_query(), xnbh, diameter); j != false; ++j) {

                    if (get<id>(*j) !=
                        get<id>(*k)) { // check if it is not the same particle
                        free_position_nbh = false;
                    }
                }


                // update the position if the place they want to move to is free and not out of bounds
                if (free_position_nbh && xnbh[0] > cell_radius && xnbh[0] < Gamma(length_x - 1) && (xnbh[1]) > cell_radius &&
                    (xnbh[1]) < length_y - 1 - cell_radius) {
                    get<position>(*k) +=
                            magn* changenbh/changenbh.norm(); // update if nothing is in the next position

                }



            }
            particles.update_positions(); // not sure if needed here


            /*
            * end of CIL
             * */


            x = get<position>(particles[particle_id(j)]); // make sure I still have it's real position


            // create an array to store random directions
            array<double, filo_number + 1> random_angle;

            // choose the number of angles where the filopodia is sent
            for (int k = 0; k < filo_number + 1; k++) {

                // uniform

                double random_angle_tem = uniformpi(gen1);

                // either one can choose any angle, even though it would lead to a movement outside the domain
                random_angle_tem = uniformpi(gen1);
                random_angle[k] = random_angle_tem;// either one can choose any angle, even though it would lead to a movement outside the domain

//                // or the filopodia can only be sent inside the domain

//                while (((x[0] + sin(random_angle_tem) * l_filo_y) < cell_radius ||
//                        ((x[0] + sin(random_angle_tem) * l_filo_y)) >
//                                Gamma(length_x - 1) || (x[1] + cos(random_angle_tem) * l_filo_y) < cell_radius ||
//                        (x[1] + cos(random_angle_tem) * l_filo_y) > length_y - cell_radius)) {
//                    random_angle_tem = uniformpi(gen1);
//                }
//                random_angle[k] = random_angle_tem;


                // normal

//
//                double random_angle_tem = normal(gen1);
//
//                // either one can choose any angle, even though it would lead to a movement outside the domain
//                //random_angle_tem = normal(gen1);
//
//                // this is for truncated normal distribution
//                do{
//                    random_angle_tem = normal(gen1);
//                }while( random_angle_tem< - M_PI/2.0 || random_angle_tem>(3.0*M_PI)/2.0);
//
//                random_angle[k] = random_angle_tem;// either one can choose any angle, even though it would lead to a movement outside the domain
//


//                // or the filopodia can only be sent inside the domain

//                while (((x[0] + sin(random_angle_tem) * l_filo_y) < cell_radius ||
//                        ((x[0] + sin(random_angle_tem) * l_filo_y)) >
//                                Gamma(length_x - 1) || (x[1] + cos(random_angle_tem) * l_filo_y) < cell_radius ||
//                        (x[1] + cos(random_angle_tem) * l_filo_y) > length_y - cell_radius)) {
//                    do{
//                        random_angle_tem = normal(gen1);
//                    }while( random_angle_tem< - M_PI/2.0 || random_angle_tem>(3.0*M_PI)/2.0);
//                }
//                random_angle[k] = random_angle_tem;







            }


            x += speed_l * vdouble2(sin(random_angle[filo_number]), cos(random_angle[filo_number]));


            if (x[0] < Gamma(theta1 - 1)) {
                x_in = x[0] * (theta1) / Gamma(theta1 - 1);
                l_filo_x = l_filo_x_in * theta1 / Gamma(theta1 - 1);
            } else {
                x_in = x[0] - (Gamma(theta1 - 1) - theta1);
            }


            bool free_position = true; // check if the neighbouring position is free

            // if this loop is entered, it means that there is another cell where I want to move
            for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

                if (get<id>(*k) !=
                    get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                    free_position = false;
                }
            }


            // update the position if the place they want to move to is free and not out of bounds
            if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
                (x[1]) < length_y - 1 - cell_radius) {


                get<position>(particles)[particle_id(j)] +=
                        speed_l * vdouble2(sin(random_angle[filo_number]),
                                           cos(random_angle[filo_number])); // update if nothing is in the next position
                get<direction>(particles)[particle_id(j)] =
                        speed_l * vdouble2(sin(random_angle[filo_number]),
                                           cos(random_angle[filo_number]));


            }
        }


        particles.update_positions(); // not sure if needed here



        if (counter % 20 == 0) {

            // save at every time step
#ifdef HAVE_VTK
            vtkWriteGrid("CILParticles", t, particles.get_grid(true));
#endif




//            ofstream output("CILMatrix" + to_string(int(t)) + ".csv");
//
//
//            output << "x, y, z, u" << "\n" << endl;
//
//
//            for (int i = 0; i < length_x * length_y; i++) {
//                for (int j = 0; j < 4; j++) {
//                    output << chemo_3col(i, j) << ", ";
//                }
//                output << "\n" << endl;
//            }


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

        //for (int j = 0; j < 1; j++) {

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



