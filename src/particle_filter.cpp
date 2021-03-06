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
   
    if (is_initialized)
        return;
    num_particles=100;

    std::default_random_engine gen;

    std::normal_distribution<double> N_x(x, std[0]);
    std::normal_distribution<double> N_y(y, std[1]);
    std::normal_distribution<double> N_theta(theta, std[2]);

    for (int i=0; i<num_particles; i++){
        Particle particle;
        particle.id=i;
        particle.x=N_x(gen);
        particle.y=N_y(gen);
        particle.theta=N_theta(gen);
        particle.weight=1;
        particles.push_back(particle);
        weights.push_back(1);

    }
    is_initialized=true;
   

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
   
    std::default_random_engine gen;
    for (int i=0; i<num_particles; i++){
        double new_x;
        double new_y;
        double new_theta;

        if (yaw_rate == 0){
            new_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
            new_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
            new_theta = particles[i].theta;
        } else {
            new_x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t)-sin(particles[i].theta));
            new_y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
            new_theta = particles[i].theta + yaw_rate * delta_t;
        }
        std::normal_distribution<double> N_x(new_x, std_pos[0]);
        std::normal_distribution<double> N_y(new_y, std_pos[1]);
        std::normal_distribution<double> N_theta(new_theta, std_pos[2]);

        particles[i].x = N_x(gen);
        particles[i].y = N_y(gen);
        particles[i].theta = N_theta(gen);
    }

   
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    for (int i = 0; i< observations.size(); i++){
        double closest_dist = std::numeric_limits<double>::max();
        double closest_id;
        for (int k = 0; k <predicted.size(); k++){
            double dist = pow(observations[i].x - predicted[k].y, 2) + pow(observations[i].y - predicted[k].y, 2);
            if (closest_dist > dist){
                closest_id = k;
                closest_dist = dist;
            }
        }
        observations[i].id = closest_id;
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
   
    double sensor_range_sq = pow(sensor_range, 2);
    for (int p = 0; p< particles.size(); p++){
        vector<int> associations;
        vector<double> sense_x;
        vector<double> sense_y;

        vector<LandmarkObs> trans_observations;
        vector<LandmarkObs> close_landmarks;
        LandmarkObs obs;
        const double cost = cos(particles[p].theta);
        const double sint = sin(particles[p].theta);
        for (int i=0; i<observations.size(); i++){
       
            LandmarkObs trans_obs;
            obs = observations[i];

            // vehicle -> map coordinates
            trans_obs.x = particles[p].x + obs.x*cost-obs.y*sint;
            trans_obs.y = particles[p].y + obs.x*sint+obs.y*cost;
            trans_observations.push_back(trans_obs);
        }

        for (int i = 0; i < map_landmarks.landmark_list.size(); i++){
            LandmarkObs obj;
            obj.x = map_landmarks.landmark_list[i].x_f;
            obj.y = map_landmarks.landmark_list[i].y_f;
            if (pow(particles[p].x - obj.x, 2) + pow(particles[p].y - obj.y, 2) < sensor_range_sq){
                obj.id = map_landmarks.landmark_list[i].id_i;
                close_landmarks.push_back(obj);


            }
        }

        dataAssociation(close_landmarks, trans_observations);

        particles[p].weight = 1;
        for (int i=0; i<trans_observations.size(); i++){
                 double x_obs = trans_observations[i].x;
                 double y_obs = trans_observations[i].y;
                 double association = trans_observations[i].id;
        
                 double mx = close_landmarks[association].x;
                 double my = close_landmarks[association].y;

                 sense_x.push_back(x_obs);
                 sense_y.push_back(y_obs);
                 associations.push_back(close_landmarks[association].id);
                 

                 // multiply weight by pertinence of obseration, set minimum weight for observations with value below certain threshold
                 const double sig_x = std_landmark[0];
                 const double sig_y = std_landmark[1];
                 const double mpsig_xy = 1/(2*M_PI)*sig_x*sig_y;
                 const double sig_xx2 = 2*sig_x*sig_x;
                 const double sig_yy2 = 2*sig_y*sig_y;
       
                 double multiplier = mpsig_xy*exp(-(pow(x_obs - mx, 2)/sig_xx2 + (pow(y_obs - my,2)/sig_yy2)));
                 // multiply weight by pertinence of obseration, set minimum weight for observations with value below certain threshold
                 const double minw = pow (10, -5);
                 // penalize more observations with 0 weight
                 const double mtw = pow (10, -5);
                 if (multiplier > 0){
                     if (multiplier > minw)
                         particles[p].weight*=multiplier;
                     else
                         particles[p].weight*=minw;}
                 else
                     particles[p].weight*=mtw;
             }
        particles[p]=SetAssociations(particles[p], associations, sense_x, sense_y);
        weights[p] = particles[p].weight;
   
   }
   
    
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   
    std::default_random_engine gen;
    std::discrete_distribution<int> distribution(weights.begin(), weights.end());

    vector<Particle> resample_particles;
    for (int i = 0; i<num_particles; i++)
        resample_particles.push_back(particles[distribution(gen)]);
    particles = resample_particles;
   
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
   

    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
   
    return(particle);
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
