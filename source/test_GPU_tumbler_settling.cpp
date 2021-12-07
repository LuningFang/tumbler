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
double drum_height = 7.62;


std::vector<ChVector<double>> generateBoundaryLayer(double cylinder_radius,
                                                   double sphere_diameter,
                                                   double cylinder_width) {
    int numLayers = (int)(cylinder_width / sphere_diameter);
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

    if (argc == 2){
        drum_height = atof(argv[1]);
    }

    if (argc != 1 && argc != 2){
        printf("incorrect input arguments. \n run ./test_GPU_tumbler_settling or pass the depth of the drum (cm).\n");
    }

    double sphere_radius = 0.0265;
    double sphere_density = 2.5;
    double box_X = drum_diameter;
    double box_Y = drum_diameter;
    double box_Z = drum_height;

    double mu_s2s = 0.16;
    double mu_s2w = 0.45;
    double rolling_fr_s2s = 0.09;
    double rolling_fr_s2w = 0.09;

    double step_size = 1e-6;
    //    double step_size = 5e-6;

    double time_end = 0.6f;
    //    double time_end = 0.01f;

    double frame_step = 0.05;
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


    // gpu_sys.SetCohesionRatio(params.cohesion_ratio);
    // gpu_sys.SetAdhesionRatio_SPH2MESH(params.adhesion_ratio_s2m);
    // gpu_sys.SetAdhesionRatio_SPH2WALL(params.adhesion_ratio_s2w);

    gpu_sys.SetGravitationalAcceleration(ChVector<float>(grav_X, grav_Y, grav_Z));

    gpu_sys.SetFrictionMode(chrono::gpu::CHGPU_FRICTION_MODE::MULTI_STEP);

    gpu_sys.SetStaticFrictionCoeff_SPH2SPH(mu_s2s);
    gpu_sys.SetStaticFrictionCoeff_SPH2WALL(mu_s2w);

    gpu_sys.SetRollingMode(CHGPU_ROLLING_MODE::SCHWARTZ);
    gpu_sys.SetRollingCoeff_SPH2SPH(rolling_fr_s2s);
    gpu_sys.SetRollingCoeff_SPH2WALL(rolling_fr_s2w);

    std::string out_dir = GetChronoOutputPath() + "tumbler_settling/";
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
        gpu_sys.AdvanceSimulation(frame_step);

        sprintf(filename, "%s/step%04d.csv", out_dir.c_str(), currframe);
        gpu_sys.WriteParticleFile(std::string(filename));

        curr_time += frame_step;
        std::cout << curr_time << std::endl;

        currframe++;
    }

    std::string initial_points_filename = "tumbler_initial_positions.csv";

    std::cout << "write file " << initial_points_filename << std::endl;

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
