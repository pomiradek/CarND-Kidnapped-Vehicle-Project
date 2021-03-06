/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if (!initialized()) {
        num_particles = 50; // better results than 20, almost same as 100
        default_random_engine gen;
        double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
        std_x = std[0];
        std_y = std[1];
        std_theta = std[2];

        // This line creates a normal (Gaussian) distribution for x
        normal_distribution<double> dist_x(x, std_x);

        // Create normal distributions for y and theta
        normal_distribution<double> dist_y(y, std_y);
        normal_distribution<double> dist_theta(theta, std_theta);

        for (int i = 0; i < num_particles; ++i) {
            weights.push_back(1.0);

            Particle particle;
            particle.id = i;
            particle.x = dist_x(gen);
            particle.y = dist_y(gen);
            particle.theta = dist_theta(gen);
            particle.weight = 1.0;

            particles.push_back(particle);

        }
	    is_initialized = true;
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
    std_x = std_pos[0];
    std_y = std_pos[1];
    std_theta = std_pos[2];

    // This creates a normal (Gaussian) distribution centered around 0
    normal_distribution<double> dist_x(0, std_x);
    normal_distribution<double> dist_y(0, std_y);
    normal_distribution<double> dist_theta(0, std_theta);

    for (int i = 0; i < num_particles; ++i) {

        double pred_x, pred_y, pred_theta;

        if (fabs(yaw_rate) < 0.0001) {
            pred_x = particles[i].x + velocity * cos(particles[i].theta) * delta_t;
            pred_y = particles[i].y + velocity * sin(particles[i].theta) * delta_t;
            pred_theta = particles[i].theta;
        } else {
            double vel_yaw = velocity/yaw_rate;
            pred_x = particles[i].x + vel_yaw * (sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta));
            pred_y = particles[i].y + vel_yaw * (cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t));
            pred_theta = particles[i].theta+yaw_rate*delta_t;
        }

        particles[i].x = pred_x + dist_x(gen);
        particles[i].y = pred_y + dist_y(gen);
        particles[i].theta = pred_theta + dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	double distance;
    for (auto &observation : observations) {
        distance = 9999999999999999.0;
        for (auto &i : predicted) {
            double new_distance = dist(observation.x , observation.y, i.x, i.y);
            if (new_distance < distance) {
                 distance = new_distance;
                observation.id = i.id;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    double exponent;
    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];

    double total_weight = 0.0;
    double gauss_norm = (1.0/(2.0 * M_PI * sig_x * sig_y));
    double sig_x_2 = 2.0 * sig_x * sig_x;
    double sig_y_2 = 2.0 * sig_y * sig_y;

    // TODO add some defensive logic to handle the situation where none of the landmarks are inside the sensor range
    // even if it doesn't occur in the simulation. At the moment, if this happens the weight will left with the initial value of 1.

    for (unsigned int j = 0; j < num_particles; ++j) {

        // Map to MAP coordinates system
        vector<LandmarkObs> observations_in_map_cor_sys;
	    for (unsigned int i = 0; i < observations.size(); i++) {
            LandmarkObs observation_in_map_cor_sys{};
            observation_in_map_cor_sys.id = i;
            // transform to map x coordinate
            observation_in_map_cor_sys.x = particles[j].x + (cos(particles[j].theta) * observations[i].x) - (sin(particles[j].theta) * observations[i].y);
            // transform to map y coordinate
            observation_in_map_cor_sys.y = particles[j].y + (sin(particles[j].theta) * observations[i].x) + (cos(particles[j].theta) * observations[i].y);
            observations_in_map_cor_sys.push_back(observation_in_map_cor_sys);
        }

        // Filter only landmarks in range of radar
        vector<LandmarkObs> valid_landmarks;
        for (auto k : map_landmarks.landmark_list) {
            if (dist(particles[j].x, particles[j].y, k.x_f, k.y_f) <= sensor_range) {
                LandmarkObs valid_landmark{};
                valid_landmark.id = k.id_i;
                valid_landmark.x = k.x_f;
                valid_landmark.y = k.y_f;
                valid_landmarks.push_back(valid_landmark);
            }
        }

        dataAssociation(valid_landmarks, observations_in_map_cor_sys);

        particles[j].weight = 1.0;

        for (auto &observations_in_map_cor_sy : observations_in_map_cor_sys) {

            for (auto &valid_landmark : valid_landmarks) {
                if (observations_in_map_cor_sy.id == valid_landmark.id) {
                    exponent = (pow(observations_in_map_cor_sy.x - valid_landmark.x, 2))/sig_x_2 + (pow(
                            observations_in_map_cor_sy.y - valid_landmark.y, 2))/sig_y_2;
                    particles[j].weight *= (gauss_norm * exp(-exponent));
                }
            }
        }
        total_weight += particles[j].weight;
	}

	// Normalize weight to 0-1
	for (unsigned int i = 0; i < particles.size(); i++) {
        particles[i].weight /= total_weight;
        weights[i] = particles[i].weight;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> resampled_particles;
    vector<double> resampled_weights;

    double beta = 0.0;
    double maxw = 2.0 * *max_element(weights.begin(), weights.end());

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<int> dindex(0, num_particles-1);
    uniform_real_distribution<double> dweight(0.0, maxw);

    /*
     * discrete_distribution<int> discrete_dist(weights.begin(), weights.end());

    vector<Particle> resampled_particles(num_particles);

    for (int i = 0; i < num_particles; ++i) {
        int index = discrete_dist(gen_);
        resampled_particles.at(i) = particles.at(index);
    }
    particles = resampled_particles;
     */

    int index = dindex(gen);

    for (int i = 0; i < num_particles; i++) {

        beta += dweight(gen);

        while (beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        resampled_particles.push_back(particles[index]);
        resampled_weights.push_back(particles[index].weight);
    }
    particles = resampled_particles;
    weights = resampled_weights;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
