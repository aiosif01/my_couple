#ifndef CELLS_H_
#define CELLS_H_

#include "biodynamo.h"
#include "precice_adapter.h"
#include "my_cell.h" // Make sure this is included

namespace bdm {

inline int Simulate(int argc, const char** argv) {
  Simulation simulation(argc, argv);
  auto* rm = simulation.GetResourceManager();

  // --- Create agents ---
  Log::Info("Simulate", "Creating agents...");
  
  // Define OpenFOAM domain bounds for verification
  const double domain_min_x = 0.0;
  const double domain_max_x = 1.0;
  const double domain_min_y = 0.0;
  const double domain_max_y = 1.0;
  const double domain_min_z = 0.0;
  const double domain_max_z = 1.0;
  
  // OpenFOAM mesh details from blockMeshDict (50×50×50 cells)
  const int openfoam_cells_per_dim = 50;
  const double openfoam_cell_size = 1.0 / openfoam_cells_per_dim; // 0.02 units
  
  Log::Info("Simulate", "OpenFOAM domain: ", domain_min_x, "-", domain_max_x, " × ",
            domain_min_y, "-", domain_max_y, " × ", domain_min_z, "-", domain_max_z);
  Log::Info("Simulate", "OpenFOAM mesh: ", openfoam_cells_per_dim, "×", openfoam_cells_per_dim, 
            "×", openfoam_cells_per_dim, " cells (cell size: ", openfoam_cell_size, ")");
  
  // BioDynaMo agent parameters
  const int agents_per_dim = 10; // 10×10×10 = 1000 agents
  double cell_diameter = 0.01; // Cell diameter (smaller than OpenFOAM cell size)
  
  // We'll still define an initial temperature as a fallback value
  // But we'll read it from OpenFOAM when possible
  double initial_temp = 300.0; // Default initial temperature if preCICE fails
  
  // Calculate step size to evenly distribute BioDynaMo agents across the domain
  int step = openfoam_cells_per_dim / agents_per_dim; // Step = 5 for 10 agents per dimension
  
  Log::Info("Simulate", "Creating BioDynaMo agents in a ", agents_per_dim, "×", 
            agents_per_dim, "×", agents_per_dim, " grid (", 
            agents_per_dim*agents_per_dim*agents_per_dim, " total agents)");
  
  int cell_count = 0;
  
  // Create BioDynaMo agents at the centers of selected OpenFOAM cells
  for (int i = 0; i < agents_per_dim; i++) {
    for (int j = 0; j < agents_per_dim; j++) {
      for (int k = 0; k < agents_per_dim; k++) {
        // Calculate OpenFOAM cell indices (evenly distributed)
        int of_cell_i = i * step + step/2; // Add half step to get to middle cells
        int of_cell_j = j * step + step/2;
        int of_cell_k = k * step + step/2;
        
        // Calculate the exact center position of this OpenFOAM cell
        double pos_x = (of_cell_i + 0.5) * openfoam_cell_size;
        double pos_y = (of_cell_j + 0.5) * openfoam_cell_size;
        double pos_z = (of_cell_k + 0.5) * openfoam_cell_size;
        
        // Create a BioDynaMo agent (cell) at the exact center of this OpenFOAM cell
        MyCell* cell = new MyCell({pos_x, pos_y, pos_z});
        cell->SetDiameter(cell_diameter);
        cell->SetTemperature(initial_temp); // Initially set to default value, will be updated with OpenFOAM data
        rm->AddAgent(cell);
        
        // Log cell creation (limit output)
        if (cell_count % 100 == 0 || cell_count < 5 || cell_count > 995) {
          Log::Info("Simulate", "Created agent ", cell_count, 
                   " at position (", pos_x, ", ", pos_y, ", ", pos_z, ")",
                   " - center of OpenFOAM cell [", of_cell_i, ", ", of_cell_j, ", ", of_cell_k, "]");
        }
        
        cell_count++;
      }
    }
  }
  
  Log::Info("Simulate", "Created ", cell_count, " agents at centers of OpenFOAM volume cells");
  
  // Verify cells were created
  uint64_t total_cells = rm->GetNumAgents();
  Log::Info("Simulate", "Total agents in resource manager: ", total_cells);
  
  // Verify all cells are within domain
  int cells_in_domain = 0;
  int cells_outside_domain = 0;
  
  rm->ForEachAgent([&](Agent* agent) {
    if (auto* cell = dynamic_cast<MyCell*>(agent)) {
      const auto& pos = cell->GetPosition();
      bool in_domain = (pos[0] >= domain_min_x && pos[0] <= domain_max_x &&
                       pos[1] >= domain_min_y && pos[1] <= domain_max_y &&
                       pos[2] >= domain_min_z && pos[2] <= domain_max_z);
      
      if (in_domain) {
        cells_in_domain++;
      } else {
        cells_outside_domain++;
        Log::Warning("Simulate", "DEBUG: After creation, cell ID ", cell->GetUid(), 
                    " position (", pos[0], ",", pos[1], ",", pos[2], 
                    ") is OUTSIDE OpenFOAM domain bounds!");
      }
    }
  });
  
  Log::Info("Simulate", "DEBUG: Cell position verification - ", cells_in_domain, 
           " cells INSIDE domain, ", cells_outside_domain, " cells OUTSIDE domain");
  
  if (total_cells == 0) {
    Log::Error("Simulate", "ERROR: No cells were created! Stopping simulation.");
    return 1;
  }

  // --- Set up preCICE coupling ---
  // CRITICAL: Follow this exact order to avoid mesh modification errors
  
  // 1. Create preCICE adapter
  Log::Info("Simulate", "Creating preCICE adapter...");
  PreciceAdapter adapter("../precice-config.xml", "cells");
  
  // 2. Register agent positions with preCICE mesh BEFORE initialization
  Log::Info("Simulate", "Registering agent positions with preCICE mesh...");
  adapter.UpdateMesh(simulation);
  
  // 3. Initialize preCICE AFTER mesh registration is complete
  Log::Info("Simulate", "Initializing preCICE connection...");
  adapter.Initialize();
  
  // --- EXPLICIT AGENT TEMPERATURE INITIALIZATION FROM OPENFOAM ---
  Log::Info("Simulate", "Initializing agent temperatures from OpenFOAM data...");
  std::vector<double> initial_temperatures;
  adapter.ReadTemperature(initial_temperatures);
  
  if (!initial_temperatures.empty()) {
    std::cout << "INITIALIZATION: Successfully received " << initial_temperatures.size() 
              << " initial temperature values from OpenFOAM" << std::endl;
    
    // Apply initial temperature values to cells
    size_t idx = 0;
    double min_init_temp = std::numeric_limits<double>::max();
    double max_init_temp = std::numeric_limits<double>::lowest();
    double sum_init_temp = 0.0;
    int cells_initialized = 0;
    
    rm->ForEachAgent([&](Agent* agent) {
      if (auto* cell = dynamic_cast<MyCell*>(agent)) {
        if (idx < initial_temperatures.size()) {
          double temp = initial_temperatures[idx];
          double default_temp = cell->GetTemperature();
          cell->SetTemperature(temp);
          
          // Log the first few cells for verification
          if (idx < 5 || idx > initial_temperatures.size() - 5) {
            const auto& pos = cell->GetPosition();
            std::cout << "INITIALIZATION: Cell at (" << pos[0] << ", " << pos[1] << ", " << pos[2]
                     << ") initialized with temperature " << temp 
                     << " (was default: " << default_temp << ")" << std::endl;
          }
          
          // Track temperature stats
          min_init_temp = std::min(min_init_temp, temp);
          max_init_temp = std::max(max_init_temp, temp);
          sum_init_temp += temp;
          cells_initialized++;
          
          idx++;
        } else {
          if (idx == initial_temperatures.size()) { // Print only once to avoid log spam
            std::cout << "INITIALIZATION: More cells than temperature values available" << std::endl;
          }
          idx++;
        }
      }
    });
    
    // Log initialization statistics
    if (cells_initialized > 0) {
      double avg_temp = sum_init_temp / cells_initialized;
      std::cout << "INITIALIZATION: Successfully initialized " << cells_initialized 
                << " agent temperatures from OpenFOAM data" << std::endl;
      std::cout << "INITIALIZATION: Temperature stats - Min: " << min_init_temp 
                << ", Max: " << max_init_temp << ", Avg: " << avg_temp << std::endl;
    } else {
      std::cout << "INITIALIZATION: WARNING - No cells were initialized with temperature data" << std::endl;
    }
  } else {
    std::cout << "INITIALIZATION: Failed to receive initial temperature data from OpenFOAM. "
              << "Using default temperature value (" << initial_temp << ")" << std::endl;
  }
  
  // --- Run simulation ---
  double dt = adapter.GetMaxTimeStep();
  Log::Info("Simulate", "Starting simulation with dt = ", dt);
  
  int timestep = 0;
  while (adapter.IsCouplingOngoing()) {
    timestep++;
    Log::Info("Simulate", "Starting timestep ", timestep);
    
    // Read temperature data from preCICE
    std::vector<double> temperatures;
    adapter.ReadTemperature(temperatures);
    
    // Use direct console output to ensure visibility regardless of log level
    if (!temperatures.empty()) {
      std::cout << "MAIN LOOP: Successfully received " << temperatures.size() << " temperature values in timestep " << timestep << std::endl;
      
      // Apply temperature values to cells
      size_t idx = 0;
      double min_applied_temp = std::numeric_limits<double>::max();
      double max_applied_temp = std::numeric_limits<double>::lowest();
      double sum_applied_temp = 0.0;
      int cells_updated = 0;
      
      rm->ForEachAgent([&](Agent* agent) {
        if (auto* cell = dynamic_cast<MyCell*>(agent)) {
          if (idx < temperatures.size()) {
            double temp = temperatures[idx];
            double old_temp = cell->GetTemperature();
            cell->SetTemperature(temp);
            
            // Only print the first cell's temperature for each timestep to avoid excessive output
            if (idx == 0) {
              const auto& pos = cell->GetPosition();
              std::cout << "MAIN LOOP: First cell at (" << pos[0] << ", " << pos[1] << ", " << pos[2]
                       << ") temperature changed from " << old_temp << " to " << temp
                       << " (delta: " << temp - old_temp << ")" << std::endl;
            }
            
            // Track temperature stats
            min_applied_temp = std::min(min_applied_temp, temp);
            max_applied_temp = std::max(max_applied_temp, temp);
            sum_applied_temp += temp;
            cells_updated++;
            
            idx++;
          } else {
            if (idx == temperatures.size()) { // Only print once to avoid excessive output
              std::cout << "MAIN LOOP: More cells than temperature values available" << std::endl;
            }
            idx++;
          }
        }
      });
      
      // Log temperature application statistics with direct console output
      if (cells_updated > 0) {
        double avg_temp = sum_applied_temp / cells_updated;
        std::cout << "MAIN LOOP: Temperature stats for timestep " << timestep
                 << " - Min: " << min_applied_temp << ", Max: " << max_applied_temp
                 << ", Avg: " << avg_temp << ", Cells updated: " << cells_updated << std::endl;
      } else {
        std::cout << "MAIN LOOP: WARNING - No cells were updated with temperature data in timestep " << timestep << std::endl;
      }
    } else {
      std::cout << "MAIN LOOP: No temperature data received in timestep " << timestep << std::endl;
    }
    
    // Run one simulation step
    simulation.GetScheduler()->Simulate(1);
    
    // Advance preCICE
    Log::Info("Simulate", "Advancing preCICE with dt = ", dt);
    adapter.Advance(dt);
    
    double new_dt = adapter.GetMaxTimeStep();
    Log::Info("Simulate", "Timestep ", timestep, " completed. New dt = ", new_dt);
    dt = new_dt;
  }
  
  Log::Info("Simulate", "Finalizing preCICE...");
  adapter.Finalize();
  
  Log::Info("Simulate", "Simulation completed successfully after ", timestep, " timesteps");
  return 0;
}

}  // namespace bdm

#endif  // CELLS_H_