## Particle-Filter
This repository contains the necessary code in C++ to localize the position of a self-driving car. The Particle filter algorithm uses the data that the sensors(RADAR, LIDAR) mounted on the car receive to localize it within required accuracy.

The required algorithm and classes used to process the data("particle_filter.h", "particle_filter.cpp") can be found in the src folder.

"Particle_Filter_in_action" is a video file that shows the particle filter in action. The green lines are the observations made by Lidars and Radars and the blue circle represents the estimated position of the car and its heading. The error in the position and yaw estimation can also be found in the video.  
