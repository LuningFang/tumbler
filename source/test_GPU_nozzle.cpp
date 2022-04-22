// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2022 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Luning
// =============================================================================
// Chrono::Gpu simulation of material flowing out of a nozzle
// =============================================================================

#include <iostream>
#include <string>
#include <cmath>
#include "chrono/core/ChGlobal.h"
#include "chrono/utils/ChUtilsSamplers.h"

#include "chrono_gpu/physics/ChSystemGpu.h"
#include "chrono_gpu/ChGpuData.h"
#include "chrono/core/ChMathematics.h"

#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_thirdparty/cxxopts/ChCLI.h"

using namespace chrono;
using namespace chrono::gpu;


// forward declaration of specs
bool GetProblemSpecs(int argc,
                     char** argv,
                     double& nozzle_diameter,
                     double& nozzle_angle,
                     double& mu_s,
                     double& mu_r,
                     std::string& output_dir);

void setupContactModel(ChSystemGpu& gran_sys){
    float kn = 1e7;
    float kt = 1e7;
    float gn = 5e4;
    float gt = 5e4;

    // normal force model
    gran_sys.SetKn_SPH2SPH(kn);
    gran_sys.SetKn_SPH2WALL(kn);
    gran_sys.SetGn_SPH2SPH(gn);
    gran_sys.SetGn_SPH2WALL(gn);

    // tangential force model 
    gran_sys.SetKt_SPH2SPH(kt);
    gran_sys.SetKt_SPH2WALL(kt);
    gran_sys.SetGt_SPH2SPH(gt);
    gran_sys.SetGt_SPH2WALL(gt);
}

int main(int argc, char* argv[]) {
    gpu::SetDataPath(std::string(PROJECTS_DATA_DIR) + "gpu/");

    // units in cm and gram
    double aperture_diameter = 0.1f;
    double nozzle_angle_deg = 60;
    float sphere_radius = 0.005;
    float sphere_density = 2.5;
    float grav_Z = -980;
    std::string output_dir = "nozzle";
    float box_X = 2;
    float box_Y = 2;
    float box_Z = 2;
    float step_size = 1e-6;
    float time_end = 3;
    double mu_s = 0.15;
    double mu_r = 0.2;


    // parse command line arguments
    if (!GetProblemSpecs(argc, argv, aperture_diameter, nozzle_angle_deg, mu_s, mu_r, output_dir)){
        std::cout << "Error parsing parameters. "  << std::endl;
        return 1;
    }


    std::cout << "nozzle diameter: " << aperture_diameter << " cm" << std::endl;
    std::cout << "nozzle angle   : " << nozzle_angle_deg << " deg" << std::endl;

    double nozzle_angle_rad = (nozzle_angle_deg/180.0) * CH_C_PI;

    // Setup simulation
    ChSystemGpu gran_sys(sphere_radius, sphere_density,
                         ChVector<float>(box_X, box_Y, box_Z));

    setupContactModel(gran_sys);

    gran_sys.SetFrictionMode(CHGPU_FRICTION_MODE::MULTI_STEP);    
    gran_sys.SetStaticFrictionCoeff_SPH2SPH(mu_s);
    gran_sys.SetStaticFrictionCoeff_SPH2WALL(3*mu_s);
    gran_sys.SetRollingMode(CHGPU_ROLLING_MODE::SCHWARTZ);
    gran_sys.SetRollingCoeff_SPH2SPH(mu_r);
    gran_sys.SetRollingCoeff_SPH2WALL(mu_r);

    gran_sys.SetGravitationalAcceleration(ChVector<float>(0, 0, grav_Z));
    gran_sys.SetBDFixed(true);

    ChVector<float> center_pt(0, 0, 0 - box_Z/3);

    float cone_slope = tan(nozzle_angle_rad/2.0f);
    float cone_offset = aperture_diameter / 2.f / cone_slope;
    float hmax = box_Z;
    float hmin = center_pt.z() + cone_offset;


    // Fill box with bodies
    std::vector<ChVector<float>> body_points;

    // padding in sampler
    float fill_epsilon = 2.02f;
    // padding at top of fill
    float fill_gap = fill_epsilon*sphere_radius;

    chrono::utils::PDSampler<float> sampler(fill_epsilon * sphere_radius);


    // width we want to fill to
    float fill_width = aperture_diameter/2 - sphere_radius;
    // height that makes this width above the cone
    float fill_height = fill_width;

    // fill to top
    float fill_bottom = hmin + fill_gap;
    // float fill_top = box_Z/6;
    float fill_top = 0;

    printf("width is %f, bot is %f, top is %f, height is %f\n", fill_width, fill_bottom, fill_top, fill_height);
    // fill box, layer by layer
    ChVector<> center(0, 0, fill_bottom);
    // shift up for bottom of box
    center.z() += fill_gap;

    while (center.z() < fill_top) {
        float diff = center.z() - center_pt.z();
        fill_width = diff * cone_slope - sphere_radius;
        auto points = sampler.SampleCylinderZ(center, fill_width, 0);
        body_points.insert(body_points.end(), points.begin(), points.end());
        center.z() += fill_epsilon * sphere_radius;
    }

    gran_sys.SetParticles(body_points);

    float sphere_mass =
        (4.f / 3.f) * sphere_density * sphere_radius * sphere_radius * sphere_radius;
    printf("%d spheres with total mass %f \n", body_points.size(), body_points.size() * sphere_mass);

    // set time integrator
    gran_sys.SetTimeIntegrator(CHGPU_TIME_INTEGRATOR::CENTERED_DIFFERENCE);
    gran_sys.SetFixedStepSize(step_size);

    // set friction mode
    filesystem::create_directory(filesystem::path(output_dir));


    // Finalize settings and initialize for runtime
    gran_sys.CreateBCConeZ(center_pt, cone_slope, hmax, hmin, false, false);

    float cyl_rad = fill_width + 8;
    printf("top of cone is at %f, cone tip is at %f, botoom of cone is %f\n", hmax, center_pt.z(), hmin);

    ChVector<float> plane_center(0, 0, hmin);
    ChVector<float> plane_normal(0, 0, 1);

    printf("center is %f, %f, %f, plane center is is %f, %f, %f\n", center_pt.x(), center_pt.y(), center_pt.z(),
           plane_center.x(), plane_center.y(), plane_center.z());
    size_t cone_plane_bc_id = gran_sys.CreateBCPlane(plane_center, plane_normal, false);

    // put a plane at the bottom of the box to count forces
    ChVector<float> box_bottom(0, 0, -box_Z / 2.f + 2.f);
    gran_sys.Initialize();

    // number of times to capture force data per second
    int captures_per_second = 200;
    // number of times to capture force before we capture a frame
    int captures_per_frame = 4;

    // assume we run for at least one frame
    float frame_step = 1.0f / captures_per_second;
    float curr_time = 0;
    int currcapture = 0;
    int currframe = 0;

    std::cout << "capture step is " << frame_step << std::endl;

    float t_remove_plane = .5;
    bool plane_active = false;

    char filename[100];
    // Run settling experiments
    while (curr_time < time_end) {
        if (!plane_active && curr_time > t_remove_plane) {
            gran_sys.DisableBCbyID(cone_plane_bc_id);
        }
        gran_sys.AdvanceSimulation(frame_step);
        curr_time += frame_step;

        // if this frame is a render frame
        if (currcapture % captures_per_frame == 0) {
            printf("rendering frame %u, curr time %f\n", currframe, curr_time);
            sprintf(filename, "%s/step%06d.csv", output_dir.c_str(), currframe++);
            gran_sys.WriteParticleFile(std::string(filename));
        }
        currcapture++;
    }

    return 0;
}


bool GetProblemSpecs(int argc,
                     char** argv,
                     double& nozzle_diameter,
                     double& nozzle_angle,
                     double& mu_s,
                     double& mu_r,
                     std::string& output_dir){

    ChCLI cli(argv[0], "Nozzle test - toy model");
    cli.AddOption<double>("Experiment", "nozzle_diameter", "diameter of nozzle [cm]", std::to_string(nozzle_diameter));
    cli.AddOption<double>("Experiment", "nozzle_angle", "angle of nozzle [deg]", std::to_string(nozzle_angle));
    cli.AddOption<double>("Simulation", "mu_s", "sliding friction coefficient", std::to_string(mu_s));
    cli.AddOption<double>("Simulation", "mu_r", "rolling friction coefficient", std::to_string(mu_r));
    cli.AddOption<std::string>("Output", "output_directory", "output directory name", output_dir);

    if (!cli.Parse(argc, argv)) {
        return false;
    }

    nozzle_diameter = cli.GetAsType<double>("nozzle_diameter");
    nozzle_angle = cli.GetAsType<double>("nozzle_angle");
    mu_s = cli.GetAsType<double>("mu_s");
    mu_r = cli.GetAsType<double>("mu_r");
    output_dir = cli.GetAsType<std::string>("output_directory");

    return true;

}