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
// used to initialize things
// =============================================================================

#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>

#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"
#include "chrono/assets/ChTriangleMeshShape.h"
#include "chrono/core/ChMathematics.h"

#include "chrono_gpu/ChGpuData.h"
#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/utils/ChGpuVisualization.h"

#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_thirdparty/cxxopts/ChCLI.h"


using namespace chrono;
using namespace chrono::gpu;

// unit cgs  (caltech experiment)
// double drum_diameter = 30;
// double drum_height = 1.6;

// cecily test in paper
// double drum_diameter = 6;
// double drum_height = 0.5;

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
                     std::string& output_dir);
// ============================================



std::vector<ChVector<double>> generateBoundaryLayer(double cylinder_radius,
                                                   double sphere_diameter,
                                                   double cylinder_width) {
    int numLayers = round(cylinder_width / sphere_diameter);
    int particles_per_layer = (int)(std::ceil(2 * CH_C_PI * cylinder_radius / sphere_diameter));

    double theta_interval = 2 * CH_C_PI / (double)particles_per_layer;
    std::vector<ChVector<double>> boundary_layer;

    double particle_x, particle_y, particle_z;
    double theta;
    for (int iz = 0; iz < numLayers; iz++) {
        particle_z = -cylinder_width / 2.0f + sphere_diameter * 0.5f + iz * sphere_diameter;

        for (int j = 0; j < particles_per_layer; j++) {
            theta = (double)j * theta_interval;
            particle_x = cos(theta) * (cylinder_radius - sphere_diameter / 2.0f);
            particle_y = sin(theta) * (cylinder_radius - sphere_diameter / 2.0f);

            boundary_layer.push_back(ChVector<double>(particle_x, particle_y, particle_z));
        }
    }
    std::cout << "Created " << boundary_layer.size() << " particles for boundary" << std::endl;
    return boundary_layer;
}

int main(int argc, char* argv[]) {

    // Parse command line arguments
    double sphere_radius = 0.0265;
    double drum_height = 0.5;
    double mu_s = 0.16;
    double mu_r = 0.09;
    double step_size = 5e-6;
    std::string output_dir = "";

    if (!GetProblemSpecs(argc, argv, drum_height, sphere_radius, mu_s, mu_r, step_size, output_dir)){
        std::cout << "Error parsing parameters. "  << std::endl;
        return 1;
    }


    double sphere_density = 2.5;
    double box_X = drum_diameter;
    double box_Y = drum_diameter;
    double box_Z = drum_height;

    double mu_s2s = mu_s;
    double mu_s2w = mu_s;
    double rolling_fr_s2s = mu_r;
    double rolling_fr_s2w = mu_r;

    double time_end = 0.6f;

    double frame_step = 0.01;
    std::cout << "radius: " << sphere_radius << std::endl;

    ChSystemGpu gpu_sys(sphere_radius, sphere_density, ChVector<float>(box_X, box_Y, box_Z));

    double grav_X = 0.0;
    double grav_Y = -981.0;
    double grav_Z = 0.0;

    bool use_material_based_model = true;
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


    // create directory
    std::string out_dir = GetChronoOutputPath() + output_dir;
    filesystem::create_directory(filesystem::path(out_dir));
    std::cout << "Successfully creating directory " << out_dir << std::endl;

    gpu_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    gpu_sys.SetFixedStepSize(step_size);
    gpu_sys.SetBDFixed(true);

    // sampling particles
    // modify sampler distance to change filler ratio
    double filler_ratio = 2.3;
    utils::HCPSampler<float> sampler(filler_ratio * sphere_radius);
    // fill out particles
    std::vector<ChVector<float>> filler_points;
    const float filler_bottom = -box_Z / 2.0 + sphere_radius;
    const float filler_radius = box_X / 2.f - sphere_radius * 2.0f;
    const float filler_top = box_Z / 2.0 - sphere_radius;

    ChVector<float> filler_center(0.0, 0.0, filler_bottom);
    while (filler_center.z() < filler_top) {
        auto points = sampler.SampleCylinderZ(filler_center, filler_radius, 0);
        filler_points.insert(filler_points.end(), points.begin(), points.end());
        filler_center.z() += filler_ratio * sphere_radius;
    }

    std::vector<ChVector<double>> boundary_points =
        generateBoundaryLayer(drum_diameter / 2.0f, sphere_radius * 2.0f, drum_height);
    std::vector<bool> fixed_points;

    std::vector<ChVector<float>> all_points;
    all_points.insert(all_points.end(), boundary_points.begin(), boundary_points.end());
    all_points.insert(all_points.end(), filler_points.begin(), filler_points.end());

    std::vector<ChVector<float>> initial_velocity;
    initial_velocity.insert(initial_velocity.end(), boundary_points.size(), ChVector<float>(0.0, 0.0, 0.0));
    // randomize initial velocity to e
    double velo_mag = 0.05;
    for (int i = 0; i < filler_points.size(); i++) {
        ChVector<float> rand_velo(ChRandom() - 0.5, ChRandom() - 0.5, ChRandom() - 0.5);
        rand_velo.Normalize();
        initial_velocity.push_back(rand_velo * velo_mag);
    }

    fixed_points.insert(fixed_points.end(), boundary_points.size(), true);
    fixed_points.insert(fixed_points.end(), filler_points.size(), false);

    gpu_sys.SetParticles(all_points, initial_velocity);
    gpu_sys.SetParticleFixed(fixed_points);

    gpu_sys.Initialize();

    float curr_time = 0.0f;

    int currframe = 0;
    char filename[100];
    ChVector<float> velocity;
    while (curr_time < time_end) {

//        suppress writing settling data to speed up
        // sprintf(filename, "%s/step%04d.csv", out_dir.c_str(), currframe);
        // gpu_sys.WriteParticleFile(std::string(filename));
        // currframe++;
        std::cout << "settling phase: "  << curr_time << " sec" << std::endl;

        gpu_sys.AdvanceSimulation(frame_step);

        curr_time += frame_step;
    }


    char initial_points_filename[300];
    sprintf(initial_points_filename, "%s/settled.csv", out_dir.c_str());
    std::cout << "write settled particle file " << initial_points_filename << std::endl;

    std::ofstream outstream(initial_points_filename, std::ios::out);
    outstream << "x, y, z" << std::endl;
    ChVector<float> pos;
    // write position up to half of the drum
    for (int i = 0; i < all_points.size(); i++) {
        if (gpu_sys.IsFixed(i) == false) {
            pos = gpu_sys.GetParticlePosition(i);
            if (pos.y() <= 0.0f) {
                outstream << std::setprecision(7) << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
            }
        }
    }

    outstream.close();

    return 0;
}



bool GetProblemSpecs(int argc,
                     char** argv,
                     double& drum_height,
                     double& particle_radius,
                     double& mu_s,
                     double& mu_r,
                     double& step_size,
                     std::string& output_dir){

    ChCLI cli(argv[0], "Tumbler test - Settling phase");
    cli.AddOption<double>("Experiment", "drum_height", "height of the tumbler [cm]", std::to_string(drum_height));
    cli.AddOption<double>("Experiment", "particle_radius", "particle radius [cm]", std::to_string(particle_radius));
    cli.AddOption<double>("Simulation", "mu_s", "sliding friction coefficient", std::to_string(mu_s));
    cli.AddOption<double>("Simulation", "mu_r", "rolling friction coefficient", std::to_string(mu_r));
    cli.AddOption<double>("Simulation", "step_size", "simulation step size ", std::to_string(step_size));
    cli.AddOption<std::string>("Output", "output_directory", "output directory name", output_dir);

    if (!cli.Parse(argc, argv)) {
        return false;
    }

    drum_height = cli.GetAsType<double>("drum_height");
    particle_radius = cli.GetAsType<double>("particle_radius");
    mu_s = cli.GetAsType<double>("mu_s");
    mu_r = cli.GetAsType<double>("mu_r");
    step_size = cli.GetAsType<double>("step_size");
    output_dir = cli.GetAsType<std::string>("output_directory");

    return true;
}