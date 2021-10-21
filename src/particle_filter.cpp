/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 80;  // TODO: Set the number of particles

  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i=0; i<num_particles; i++){
    Particle particle;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    this->particles.push_back(particle);
    this->weights.push_back(1.0);

  }
  this->is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;

  for (int i=0; i< this->num_particles; i++){
    if (yaw_rate == 0){
      this->particles[i].x += velocity * delta_t * cos(this->particles[i].theta);
      this->particles[i].y += velocity * delta_t * sin(this->particles[i].theta);
    }else{
      this->particles[i].x += (velocity / yaw_rate) * (sin(this->particles[i].theta + delta_t * yaw_rate) - sin(particles[i].theta));
      this->particles[i].y += (velocity / yaw_rate) * (cos(this->particles[i].theta) - cos(this->particles[i].theta + delta_t * yaw_rate));
      this->particles[i].theta += yaw_rate * delta_t;
    }
    normal_distribution<double> dist_x(this->particles[i].x, std_pos[0]);
    normal_distribution<double> dist_y(this->particles[i].y, std_pos[1]);
    normal_distribution<double> dist_theta(this->particles[i].theta, std_pos[2]);

    this->particles[i].x = dist_x(gen);
    this->particles[i].y = dist_y(gen);
    this->particles[i].theta = dist_theta(gen);

  }


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
  for (size_t i=0; i<observations.size(); i++){
    vector<double> landmark_dist;
    for (size_t j=0; j<predicted.size(); j++){
      double d = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      landmark_dist.push_back(d);
    }
    int minElementIndex = std::min_element(landmark_dist.begin(),landmark_dist.end()) - landmark_dist.begin();
    observations[i].id = minElementIndex;
  }



}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */


  vector <LandmarkObs> observations_map;
  for (std::size_t p=0; p<this->particles.size(); p++){
    // convert observations from local (car) coordinates to global (map) coordinates
    for (std::size_t o=0; o<observations.size(); o++){
      LandmarkObs obs_map;
      obs_map.id = observations[o].id;
      obs_map.x = this->particles[p].x + cos(this->particles[p].theta) * observations[o].x - sin(this->particles[p].theta) * observations[p].y;
      obs_map.y = this->particles[p].y + sin(this->particles[p].theta) * observations[o].x + cos(this->particles[p].theta) * observations[p].y;
      observations_map.push_back(obs_map);
      }

    // get the map landmarks that are in the sensor range
    vector <LandmarkObs> range_landmarks;
    for (std::size_t i=0; i<map_landmarks.landmark_list.size(); i++){
      if (dist(this->particles[p].x, this->particles[p].y, map_landmarks.landmark_list[i].x_f, map_landmarks.landmark_list[i].y_f) < sensor_range ){
        LandmarkObs landmark;

        landmark.id = map_landmarks.landmark_list[i].id_i;
        landmark.x = map_landmarks.landmark_list[i].x_f;
        landmark.y = map_landmarks.landmark_list[i].y_f;

        range_landmarks.push_back(landmark);
      }
    }

    //compute coordiates of nearest landmarks to the observations
    this->dataAssociation(range_landmarks, observations_map);
    particles[p].weight = 1;
    for(std::size_t i=0; i<observations_map.size(); i++){
      double sig_x = std_landmark[0];
      double sig_y = std_landmark[1];
      double mu_x = range_landmarks[observations_map[i].id].x;
      double mu_y = range_landmarks[observations_map[i].id].y;
      double weight_i = multiv_prob(sig_x, sig_y, observations_map[i].x, observations_map[i].y, mu_x, mu_y);
      this->particles[p].weight *= weight_i;
    }
    this->weights[p] = this->particles[p].weight;
  }
  //normalize the particle weights
  auto sum_of_weights = std::accumulate(this->weights.begin(), this->weights.end(), 0.0);
  for(size_t i=0; i<this->weights.size(); i++){
    this->weights.at(i) = weights[i]/sum_of_weights;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> sampled_particles;

  std::default_random_engine gen;
  std::discrete_distribution<> d(this->weights.begin(), this->weights.end());

  for(std:: size_t i=0; i<particles.size(); i++){
    sampled_particles.push_back(this->particles.at(d(gen)));
  }
  this->particles = sampled_particles;


}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
