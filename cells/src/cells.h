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
  
  Log::Info("Simulate", "DEBUG: OpenFOAM domain bounds: X[", domain_min_x, ",", domain_max_x, 
           "], Y[", domain_min_y, ",", domain_max_y, "], Z[", domain_min_z, ",", domain_max_z, "]");
  
  // Create cells in a grid pattern
  const int num_cells = 10;
  for (int i = 0; i < num_cells; i++) {
    double x = 0.1 + (i % 3) * 0.3; // Spread across x (3 columns)
    double y = 0.1 + (i / 3) * 0.3; // Spread across y (4 rows)
    double z = 0.5;                 // All at same z height
    
    // Verify position is within domain bounds
    bool in_domain = (x >= domain_min_x && x <= domain_max_x &&
                     y >= domain_min_y && y <= domain_max_y &&
                     z >= domain_min_z && z <= domain_max_z);
    
    if (!in_domain) {
      Log::Warning("Simulate", "DEBUG: Cell position (", x, ",", y, ",", z, 
                  ") is OUTSIDE OpenFOAM domain bounds!");
    }
    
    // Create cell and add to simulation
    MyCell* cell = new MyCell({x, y, z});
    cell->SetDiameter(0.01);
    cell->SetTemperature(300.0); // Initial temperature
    rm->AddAgent(cell);
    
    Log::Info("Simulate", "Created cell ", i, " at position (", x, ", ", y, ", ", z, 
             "), domain check: ", (in_domain ? "INSIDE" : "OUTSIDE"));
  }
  
  // Verify cells were created
  uint64_t total_cells = rm->GetNumAgents();
  Log::Info("Simulate", "Created ", total_cells, " cells");
  
  // Added verification: Check all cell positions again after creation
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