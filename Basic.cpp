//
// Created by Rasa on 08/01/2020.
// this is a simple random walk model with non-growing domain
//  SAMPLE INSIDE has to be commented if instead of only sampling inside, I would assume that a cell does not move if it attempts to move outside the domain
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

    // Initialise variables

    bool domain_growth = false; // if false change length_x to 1100, true 300
    int length_x;
    if (domain_growth == false){
        length_x = 1100;
    }
    else{
        length_x = 300;// length in x direction of the chemoattractant matrix, starts from 300
    }

    int length_y = 120;
    double Lt_old = length_x;
    int real_length_y = 120;
    const double final_time = 1080; // 18hrs number of timesteps, 1min - 1timestep, from 6h tp 24hours. 1440 if 24hrs
    double t = 0.0; // initialise time
    double dt = 1.0; // time step
    double dx = 1; // maybe 1/100
    int counter = 0; // to count simulations
    double speed_l = 0.8;//10 * 0.05;// 0.05;//1;//0.05; // speed of a leader cell
    const size_t N = 5; // initial number of cells
    double sigma = 2.0;
    double mean =0.0;

    vdouble2 random_vector; // for cell position updates

    // cell paramters
    double cell_radius = 7.5;//0.5; // radius of a cell
    const double diameter =
            2 * cell_radius; // diameter of a cell

    // cell variable
    ABORIA_VARIABLE(radius, double, "radius");
    typedef Particles<std::tuple<radius>, 2> particle_type; // 2 stands for dimension

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
    

    // growth function

    // I will mainly use its first derivative with respect to x

    VectorXd Gamma_x = VectorXd::Zero(length_x);
    VectorXd Gamma = VectorXd::Zero(length_x);
    VectorXd Gamma_t = VectorXd::Zero(length_x);
    VectorXd Gamma_old = VectorXd::Zero(length_x);

    for (int i = 0; i < length_x; i++) {
        Gamma(i) = i;
        Gamma_old(i) = Gamma(i);
    }



    // for total length
    double Lt = 0;
    double Ltdot = 0;



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
                                                            double(N)); // x=radius, uniformly in y

////             // for domain width 240
//            get<position>(particles[i]) = vdouble2(cell_radius,
//                                                   real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                   0.5 * double(real_length_y - 1) /
//                                                   double(N)); // x=radius, uniformly in y

            // for domain width 60
//           get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                            0.5 * double(length_y - 1) /
//                                                            double(N)); // x=radius, uniformly in y

    }

    // initialise neighbourhood search, note that the domain will grow in x direction, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(0, 0), 5 * vdouble2(length_x, length_y), vbool2(false, false));


    vtkWriteGrid("particles", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed(t * n_seed); // choose different seeds to obtain different random numbers
    //std::uniform_real_distribution<double> uniformpi(-M_PI/3.5,M_PI/3.5 +  M_PI); // 0 to up, M_PI/2 move straigth M_PI downt

    std::normal_distribution<double> normal(M_PI/2,sigma); // normal distribution for filopodia

    // to have the same as Lennard-Jones
//    std::normal_distribution<double> normalX(mean,1.0); // mean 1 variance 1
//    std::normal_distribution<double> normalY(0.0,1.0); // mean 0 variance 1


    //for each timestep
    while (t < final_time) {



        //      insert new cells

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


    // position update
    for (int j = 0; j < particles.size(); j++) {


        vdouble2 x; // use variable x for the position of cells
        x = get<position>(particles[particle_id(j)]);
        vdouble2 x_temp; // temp to check whether out of bounds

        double x_in;

  //       this is what I had with an angle
        // create an array to store random directions
        double random_angle;

            // normal


            double random_angle_tem = normal(gen1);

            // either one can choose any angle, even though it would lead to a movement outside the domain
            //random_angle_tem = normal(gen1);

            // this is for truncated normal distribution, assume that it has to be inside the domain, sample inside the domain
            int i = 0;
            do{
                i = i+1;
                random_angle_tem = normal(gen1);
                x_temp = x + speed_l * vdouble2(sin(random_angle_tem), cos(random_angle_tem));

            }while( random_angle_tem< - M_PI/2.0 || random_angle_tem>(3.0*M_PI)/2.0 || x_temp[0] < cell_radius || x_temp[0] > Gamma(length_x - 1)-cell_radius || (x_temp[1]) < cell_radius ||
                                                                                       (x_temp[1]) > length_y - 1 - cell_radius); // SAMPLE INSIDE, for REJECT I commented the second part
            random_angle = random_angle_tem;// either one can choose any angle, even though it would lead to a movement outside the domain


        x += speed_l * vdouble2(sin(random_angle), cos(random_angle));

// to make consistent with Lennard-Jones

////        random_vector[0] = normalX(gen1);
////        random_vector[1] = normalY(gen1);
////
////        random_vector = random_vector/random_vector.norm();
//        // this is for truncated normal distribution, assume that it has to be inside the domain, sample inside the domain
//            int i = 0;
//            do{
//                i = i+1;
//                random_vector[0] = normalX(gen1);
//                random_vector[1] = normalY(gen1);
//
//                random_vector = random_vector/random_vector.norm();
//                x_temp = x + speed_l * random_vector;
//
//            }while( x_temp[0] < cell_radius || x_temp[0] > Gamma(length_x - 1)-cell_radius || (x_temp[1]) < cell_radius ||
//                                                                                       (x_temp[1]) > length_y - 1 - cell_radius); // SAMPLE INSIDE, for REJECT I commented the second part
//
//
//        x += speed_l * random_vector;

        //end of the same as Lennard Jones


        bool free_position = true; // check if the neighbouring position is free

        // if this loop is entered, it means that there is another cell where I want to move
        for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {

            if (get<id>(*k) !=
                get<id>(particles[particle_id(j)])) { // check if it is not the same particle
                free_position = false;
            }
        }


        // update the position if the place they want to move to is free and not out of bounds
        if (free_position && x[0] > cell_radius && x[0] < Gamma(length_x - 1) -cell_radius && (x[1]) > cell_radius &&
            (x[1]) < length_y - 1 - cell_radius) {

// when I look at the angle
            get<position>(particles)[particle_id(j)] +=
                    speed_l * vdouble2(sin(random_angle),
                                       cos(random_angle)); // update if nothing is in the next position
//           get<position>(particles)[particle_id(j)] += speed_l * random_vector; // update if nothing is in the next position, same as Lennard Jones


        }
    }



//            if (counter % 20 == 0) {
//
//                // save at every time step
//            #ifdef HAVE_VTK
//                vtkWriteGrid("GrowingParticlesSigma10width120grows", t, particles.get_grid(true));
//            #endif
//
//
//            }
//////


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
    for (int i = 0; i<N; i++){
        vdouble2 xposi = get<position>(particles[i]);
        //cout << "position leaders " << get<position>(particles[i]) << endl;
        cout << xposi[0] << endl;
    }



    return proportions;







    }


// parameter analysis
int main() {

    const int number_parameters = 1; // parameter range
    const int sim_num = 1;

    VectorXi vector_check_length = proportions(2); //just to know what the length is
//
    int num_parts = vector_check_length.size(); // number of parts that I partition my domain
//    cout << "length " << vector_check_length.size() << endl;
 //   int num_parts = 20; // for 1800 timesteps
    MatrixXf sum_of_all = MatrixXf::Zero(num_parts, number_parameters); // sum of the values over all simulations

//     n would correspond to different seeds
//     parallel programming
#pragma omp parallel for
    for (int n = 1; n < sim_num; n++) {


        //initialise the matrix to store the values
        MatrixXi numbers = MatrixXi::Zero(num_parts, number_parameters);


        numbers.block(0, 0, num_parts, 1) = proportions( n);

//        // This is what I am using for MATLAB
//        ofstream output2("sepdata18hrs-growing-width120-normal1p0var" + to_string(n) + ".csv");
//
//        for (int i = 0; i < numbers.rows(); i++) {
//
//            for (int j = 0; j < numbers.cols(); j++) {
//
//              output2 << numbers(i, j) << ", ";
//
//                sum_of_all(i, j) += numbers(i, j);
//
//            }
//            output2 << "\n" << endl;
//        }
//
    }
    /*
    * will store everything in one matrix, the entries will be summed over all simulations
    */

//    //ofstream output3("FIRST075_twice_speed_later.csv");
//    ofstream output3("DATA18hrsgrowing-width120-normal1p0var.csv");
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