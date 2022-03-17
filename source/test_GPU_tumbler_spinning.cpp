// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2019 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Luning Fang
// =============================================================================
// drum spinning
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>
#include <ctime>

#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/core/ChMathematics.h"

#include "chrono_gpu/ChGpuData.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuJsonParser.h"
#include "chrono_gpu/utils/ChGpuVisualization.h"

#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_thirdparty/cxxopts/ChCLI.h"


using namespace chrono;
using namespace chrono::gpu;


// powder group
// drum height and diameter about 3 inches
// first figure out how drum height influence the result
double drum_diameter = 7.62;

// forward declaration of specs
bool GetProblemSpecs(int argc,
                     char** argv,
                     double& drum_height,
                     double& particle_radius,
                     double& mu_s,
                     double& mu_r,
                     double& step_size,
                     std::string& output_dir,
                     double& drum_speed);
// ============================================

void tokenizeCSVLine(std::ifstream& istream, std::vector<float>& data) {
    std::string line;
    std::getline(istream, line);  // load in current line
    std::stringstream lineStream(line);
    std::string cell;

    // iterate over cells
    while (std::getline(lineStream, cell, ',')) {
        data.push_back(std::stof(cell));
    }
}

// load sphere positions from a checkpoint file
void loadParticlePosition(std::string infile, std::vector<ChVector<float>> &sphere_positions) {
    // file stream to load in
    std::ifstream ptFile(infile);
    std::string tmp_line;
    std::getline(ptFile, tmp_line);  // skip first header line

    while (ptFile.good()) {
        std::vector<float> line_data;
        tokenizeCSVLine(ptFile, line_data);

        if (line_data.size() != 0){
            ChVector<float> curr_pos(line_data.at(0), line_data.at(1),line_data.at(2));
            sphere_positions.push_back(curr_pos);
        }
    }
}

std::vector<ChVector<float>> generateBoundaryLayer(double cylinder_radius,
                                                   double sphere_diameter,
                                                   double cylinder_width) {
    int numLayers = round(cylinder_width / sphere_diameter);
    int particles_per_layer = (int)(std::ceil(2 * CH_C_PI * cylinder_radius / sphere_diameter));

    double theta_interval = 2 * CH_C_PI / (double)particles_per_layer;
    std::vector<ChVector<float>> boundary_layer;

    double particle_x, particle_y, particle_z;
    double theta;
    for (int iz = 0; iz < numLayers; iz++) {
        particle_z = -cylinder_width / 2.0f + sphere_diameter * 0.5f + iz * sphere_diameter;

        for (int j = 0; j < particles_per_layer; j++) {
            theta = j * theta_interval;
            particle_x = cos(theta) * (cylinder_radius - sphere_diameter / 2.0f);
            particle_y = sin(theta) * (cylinder_radius - sphere_diameter / 2.0f);

            boundary_layer.push_back(ChVector<float>(particle_x, particle_y, particle_z));
        }
    }
    std::cout << "Created " << boundary_layer.size() << "particles for boundary" << std::endl;
    return boundary_layer;
}

// boundary layer position update
void updateBoundaryLayerPosition(double cyl_radius,
                                 std::vector<int> boundary_idx,
                                 double omega,
                                 ChSystemGpu& gpu_sys,
                                 double dt) {
    int index;
    ChVector<float> position;
    double vx, vy;
    double sin_theta, cos_theta;
    ChVector<double> position_curr;
    for (int i = 0; i < boundary_idx.size(); i++) {
        index = boundary_idx.at(i);
        position = gpu_sys.GetParticlePosition(index);
        sin_theta = position.y() / cyl_radius;
        cos_theta = position.x() / cyl_radius;
        vx = -omega * cyl_radius * sin_theta;
        vy = omega * cyl_radius * cos_theta;

        position_curr.x() = position.x() + dt * vx;
        position_curr.y() = position.y() + dt * vy;
        position_curr.z() = position.z();

        gpu_sys.SetParticlePosition(index, position_curr);
    }
}

// boundary layer velocity update
void updateBoundaryLayerVelocity(double cyl_radius, std::vector<int> boundary_idx, double omega, ChSystemGpu& gpu_sys) {
    int index;
    ChVector<float> position;
    ChVector<float> velocity;double drum_height = 7.62;

    double sin_theta, cos_theta;
    double vx, vy;
    for (int i = 0; i < boundary_idx.size(); i++) {
        index = boundary_idx.at(i);
        position = gpu_sys.GetParticlePosition(index);

        sin_theta = position.y() / cyl_radius;
        cos_theta = position.x() / cyl_radius;

        vx = -omega * cyl_radius * sin_theta;
        vy = omega * cyl_radius * cos_theta;

        velocity.x() = -omega * cyl_radius * sin_theta;
        velocity.y() = omega * cyl_radius * cos_theta;
        velocity.z() = 0.0f;

        gpu_sys.SetParticleVelocity(index, velocity);
    }
}

// find boundary layer index
std::vector<int> findBoundaryLayerIndex(ChSystemGpu& gpu_sys, double cyl_radius, int nb) {
    std::vector<int> boundary_idx;
    for (int i = 0; i < nb; i++) {
        ChVector<float> pos = gpu_sys.GetParticlePosition(i);
        double dist = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());
        if (std::abs(dist - cyl_radius) < 1e-5) {
            boundary_idx.push_back(i);
        }
    }
    return boundary_idx;
}

int main(int argc, char* argv[]) {

    // Parse command line arguments
    double sphere_radius = 0.0265;
    double drum_height = 0.5;
    double mu_s = 0.16;
    double mu_r = 0.09;
    double step_size = 5e-6;
    std::string output_dir = "test";
    double drum_omega_rpm = 10;

    if (!GetProblemSpecs(argc, argv, drum_height, sphere_radius, mu_s, mu_r, step_size, output_dir, drum_omega_rpm)){
        std::cout << "Error parsing parameters. "  << std::endl;
        return 1;
    }

    // caltech test
    // double drum_diameter = 30;
    // double drum_height = 1.6;

    double sphere_density = 2.5;

    double box_X = drum_diameter;
    double box_Y = drum_diameter;
    double box_Z = drum_height * 1.2;

    // double mu_s2s = 0.16;
    // double mu_s2w = 0.45;

    double mu_s2s = mu_s;
    double mu_s2w = 3*mu_s;

    double rolling_fr_s2s = mu_r;
    double rolling_fr_s2w = mu_r;

    double drum_omega = drum_omega_rpm * 2.0f * CH_C_PI / 60.0f;  // drumming spinning rad/s

    double time_end = 120.0f/drum_omega_rpm;
    ChSystemGpu gpu_sys(sphere_radius, sphere_density, ChVector<float>(box_X, box_Y, box_Z));

    double grav_X = 0.0;
    double grav_Y = -981.0;
    double grav_Z = 0.0;

    double cor_p = 0.97;
    double cor_w = 0.82;
    double youngs_modulus = 7e8;  // 70Mpa = 7e7Pa = 7e8 g/(cms^2)
    double poisson_ratio = 0.24;

    gpu_sys.UseMaterialBasedModel(true);

    gpu_sys.SetYoungModulus_SPH(youngs_modulus);
    gpu_sys.SetYoungModulus_WALL(youngs_modulus);

    gpu_sys.SetRestitution_SPH(cor_p);
    gpu_sys.SetRestitution_WALL(cor_w);

    gpu_sys.SetPoissonRatio_SPH(poisson_ratio);
    gpu_sys.SetPoissonRatio_WALL(poisson_ratio);

    gpu_sys.SetGravitationalAcceleration(ChVector<float>(grav_X, grav_Y, grav_Z));

    gpu_sys.SetFrictionMode(chrono::gpu::CHGPU_FRICTION_MODE::MULTI_STEP);

    gpu_sys.SetStaticFrictionCoeff_SPH2SPH(mu_s2s);
    gpu_sys.SetStaticFrictionCoeff_SPH2WALL(mu_s2w);

    gpu_sys.SetRollingMode(CHGPU_ROLLING_MODE::SCHWARTZ);
    gpu_sys.SetRollingCoeff_SPH2SPH(rolling_fr_s2s);
    gpu_sys.SetRollingCoeff_SPH2WALL(rolling_fr_s2w);


    // check to see if output directory exist
    std::string out_dir = GetChronoOutputPath() + output_dir;
    filesystem::path output_folder(out_dir);
    if (!output_folder.exists()){
        std::cout << "ERROR: output folder " << output_dir << " does not exist. " << std::endl; 
        return 1;
    }

    std::string subfolder = out_dir + "/omega_" + std::to_string((int)drum_omega_rpm) + "_rpm";
    filesystem::create_directory(filesystem::path(subfolder));
    std::cout << "Created folder " << subfolder.c_str() << std::endl;


    gpu_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    gpu_sys.SetFixedStepSize(step_size);
    gpu_sys.SetBDFixed(true);
    gpu_sys.SetVerbosity(CHGPU_VERBOSITY::INFO);

    // create front plate at location drum_height/2
    ChVector<float> front_plate_pos(0.0, 0.0, drum_height / 2.0);
    ChVector<float> front_plate_normal(0.0, 0.0, -1.0f);
    size_t front_plate_id = gpu_sys.CreateBCPlane(front_plate_pos, front_plate_normal, true);

    // bottom plate at location -H/2, for recording force
    ChVector<float> back_plate_pos(0.0, 0.0, -drum_height / 2.0);
    ChVector<float> back_plate_normal(0.0, 0.0, 1.0f);
    size_t back_plate_id = gpu_sys.CreateBCPlane(back_plate_pos, back_plate_normal, true);

    ChVector<double> spin_omega(0.0f, 0.0f, drum_omega);
    gpu_sys.SetBCPlaneRotation(front_plate_id, ChVector<double>(0.0f, 0.0f, 0.0f), spin_omega);
    gpu_sys.SetBCPlaneRotation(back_plate_id, ChVector<double>(0.0f, 0.0f, 0.0f), spin_omega);

    std::vector<ChVector<float>> filler_points;

    std::string initial_points_filename = out_dir + "/settled.csv";
    loadParticlePosition(initial_points_filename, filler_points);
    std::cout << "loaded " << filler_points.size() << " particles for filling" << std::endl;
    if (filler_points.size() == 0){
        std::cout << "ERROR: No particle loaded." << std::endl;
        return 1;
    }


    std::vector<ChVector<float>> boundary_points =
        generateBoundaryLayer(drum_diameter / 2.0f, sphere_radius * 2.0f, drum_height);

    std::vector<bool> fixed_points;
    std::vector<ChVector<float>> all_points;
    all_points.insert(all_points.end(), boundary_points.begin(), boundary_points.end());
    all_points.insert(all_points.end(), filler_points.begin(), filler_points.end());

    fixed_points.insert(fixed_points.end(), boundary_points.size(), true);
    fixed_points.insert(fixed_points.end(), filler_points.size(), false);

    gpu_sys.SetParticles(all_points);
    gpu_sys.SetParticleFixed(fixed_points);

    gpu_sys.Initialize();

    std::vector<int> boundary_idx =
        findBoundaryLayerIndex(gpu_sys, drum_diameter / 2.0f - sphere_radius, all_points.size());

    float curr_time = 0.0f;
    int currframe = 0;

    char filename[100];
    ChVector<float> velocity;

    // write initial position
    sprintf(filename, "%s/step%06d.csv", subfolder.c_str(), currframe);
    gpu_sys.WriteParticleFile(std::string(filename));

    // int steps_per_frame = 50000;
    //    int steps_per_frame = 2000; // how often to output positions
    int steps_per_frame = 2500;  // how often to output positions

    double frame_step = 1e-4;  // how often to update the boundary
    std::cout << "boundary particles updated at every " << frame_step << " seconds" << std::endl;

    


    // write specs before simulation 


    if (!GetProblemSpecs(argc, argv, drum_height, sphere_radius, mu_s, mu_r, step_size, output_dir, drum_omega_rpm)){
        std::cout << "Error parsing parameters. "  << std::endl;
        return 1;
    }

    // Write file with stats for the settling phase
    std::ofstream outf;
    outf.open(subfolder + "/params.info", std::ios::out);
    outf << "Number particles:             " << filler_points.size() << endl;
    outf << "Particle radius (cm):         " << sphere_radius << endl;
    outf << "Drum depth (cm):              " << drum_height << endl;
    outf << "Sliding friction coefficient: " << mu_s << endl;
    outf << "Rolling friction coefficient: " << mu_r << endl;
    outf << "Drum rotation speed (rpm):    " << drum_omega_rpm << endl;
    outf << "Simulation step size:         " << step_size << endl;

    clock_t cpu_start_time = std::clock();
    while (curr_time < time_end) {
        gpu_sys.AdvanceSimulation(frame_step);

        updateBoundaryLayerPosition(drum_diameter / 2.0f - sphere_radius, boundary_idx, drum_omega, gpu_sys,
                                    frame_step);

        updateBoundaryLayerVelocity(drum_diameter / 2.0f - sphere_radius, boundary_idx, drum_omega, gpu_sys);

        if (currframe % steps_per_frame == 0) {

            // only write particle info during the last half revolution at 0.1 sec interval
            if (curr_time > 0.75 * time_end){
                sprintf(filename, "%s/step%06d.csv", subfolder.c_str(), (int)(currframe / steps_per_frame) + 1);
                gpu_sys.WriteParticleFile(std::string(filename));
            }

            std::cout << curr_time << std::endl;
        }

        curr_time += frame_step;
        currframe++;
    }
    clock_t end_time = std::clock();
    double computation_time = (end_time - cpu_start_time)/CLOCKS_PER_SEC;
    outf << "Simulation end time:          " << computation_time << endl;


    return 0;
}


bool GetProblemSpecs(int argc,
                     char** argv,
                     double& drum_height,
                     double& particle_radius,
                     double& mu_s,
                     double& mu_r,
                     double& step_size,
                     std::string& output_dir,
                     double& drum_speed){
    ChCLI cli(argv[0], "Tumbler test - Spinning phase");
    cli.AddOption<double>("Experiment", "drum_height", "height of the tumbler [cm]", std::to_string(drum_height));
    cli.AddOption<double>("Experiment", "particle_radius", "particle radius [cm]", std::to_string(particle_radius));
    cli.AddOption<double>("Simulation", "mu_s", "sliding friction coefficient", std::to_string(mu_s));
    cli.AddOption<double>("Simulation", "mu_r", "rolling friction coefficient", std::to_string(mu_r));
    cli.AddOption<double>("Simulation", "step_size", "simulation step size ", std::to_string(step_size));
    cli.AddOption<std::string>("Output", "output_directory", "output directory name", output_dir);
    cli.AddOption<double>("Experiment", "drum_speed", "drum rotation speed [rpm]", std::to_string(particle_radius));

    if (!cli.Parse(argc, argv)) {
        return false;
    }

    drum_height = cli.GetAsType<double>("drum_height");
    particle_radius = cli.GetAsType<double>("particle_radius");
    mu_s = cli.GetAsType<double>("mu_s");
    mu_r = cli.GetAsType<double>("mu_r");
    step_size = cli.GetAsType<double>("step_size");
    output_dir = cli.GetAsType<std::string>("output_directory");
    drum_speed = cli.GetAsType<double>("drum_speed");

    return true;

}