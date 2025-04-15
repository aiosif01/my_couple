#ifndef CELLS_H_
#define CELLS_H_

#include "biodynamo.h"
#include "precice_adapter.h"
#include "my_cell.h"
#include <algorithm>
#include <map>

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
  // We'll create agents exactly at OpenFOAM cell centers for better mapping
  const int spacing = 5; // Use every 5th OpenFOAM cell
  double cell_diameter = openfoam_cell_size * 0.5; // Cell diameter smaller than OF cell size
  
  // We'll still define an initial temperature as a fallback value
  // But we'll read it from OpenFOAM when possible
  double initial_temp = 300.0; // Default initial temperature if preCICE fails
  
  Log::Info("Simulate", "Creating BioDynaMo agents at OpenFOAM cell centers");
  
  int cell_count = 0;
  std::map<int, MyCell*> of_cell_to_agent_map; // Maps OF cell indices to BDM agents
  
  // Create BioDynaMo agents at the centers of selected OpenFOAM cells
  for (int i = 0; i < openfoam_cells_per_dim; i += spacing) {
    for (int j = 0; j < openfoam_cells_per_dim; j += spacing) {
      for (int k = 0; k < openfoam_cells_per_dim; k += spacing) {
        // Calculate the exact center position of this OpenFOAM cell
        double pos_x = (i + 0.5) * openfoam_cell_size;
        double pos_y = (j + 0.5) * openfoam_cell_size;
        double pos_z = (k + 0.5) * openfoam_cell_size;
        
        // Create a BioDynaMo agent (cell) at the exact center of this OpenFOAM cell
        MyCell* cell = new MyCell({pos_x, pos_y, pos_z});
        cell->SetDiameter(cell_diameter);
        cell->SetTemperature(initial_temp); // Initially set to default value, will be updated with OpenFOAM data
        
        // Calculate the linear index of this OpenFOAM cell
        int of_cell_index = i + j * openfoam_cells_per_dim + k * openfoam_cells_per_dim * openfoam_cells_per_dim;
        
        // Store in our mapping
        of_cell_to_agent_map[of_cell_index] = cell;
        
        // Add the agent to the simulation
        rm->AddAgent(cell);
        
        // Log cell creation (limit output)
        if (cell_count < 5 || cell_count % 100 == 0) {
          Log::Info("Simulate", "Created agent ", cell_count, 
                   " at position (", pos_x, ", ", pos_y, ", ", pos_z, ")",
                   " - OpenFOAM cell index: ", of_cell_index);
        }
        
        cell_count++;
      }
    }
  }
  
  Log::Info("Simulate", "Created ", cell_count, " agents at centers of OpenFOAM volume cells");
  
  // Add special test cells at cardinal directions for debugging
  Log::Info("Simulate", "------ Adding Cardinal Direction Test Cells ------");
  const std::vector<std::pair<std::string, Real3>> testPoints = {
    {"Origin", {0.1, 0.1, 0.1}},
    {"X-axis", {0.9, 0.1, 0.1}},
    {"Y-axis", {0.1, 0.9, 0.1}},
    {"Z-axis", {0.1, 0.1, 0.9}},
    {"Center", {0.5, 0.5, 0.5}}
  };

  for (const auto& [name, pos] : testPoints) {
    MyCell* cell = new MyCell(pos);
    cell->SetDiameter(cell_diameter);
    cell->SetTemperature(350.0); // Higher temperature to distinguish them
    rm->AddAgent(cell);
    
    Log::Info("Simulate", "Added test cell '", name, "' at position (", 
              pos[0], ",", pos[1], ",", pos[2], ")");
  }
  
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
        Log::Warning("Simulate", "Cell at position (", pos[0], ",", pos[1], ",", pos[2], 
                    ") is OUTSIDE OpenFOAM domain bounds!");
      }
    }
  });
  
  Log::Info("Simulate", "Cell position verification - ", cells_in_domain, 
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
  // Check if we'll need to write initial data (safe to call before initialization)
  bool will_need_initial_data = adapter.WillRequireInitialData();
  Log::Info("Simulate", "preCICE adapter ", (will_need_initial_data ? "will" : "will not"), " require initial data");

  // Initialize preCICE
  adapter.Initialize();

  // *** REMOVE THIS CODE BLOCK THAT'S CAUSING THE ERROR ***
  // if (adapter.RequiresInitialData()) {
  //   Log::Info("Simulate", "Writing initial data to satisfy preCICE requirements...");
  //   // This is a placeholder - implement WriteData method if needed
  //   // adapter.WriteData();
  // }

  // --- EXPLICIT AGENT TEMPERATURE INITIALIZATION FROM OPENFOAM ---
  Log::Info("Simulate", "Reading initial temperature data from OpenFOAM...");
  std::vector<double> initial_temperatures;
  adapter.ReadTemperature(initial_temperatures);

  // Log temperature data status
  if (!initial_temperatures.empty()) {
    Log::Info("Simulate", "Successfully received ", initial_temperatures.size(), 
              " initial temperature values from OpenFOAM (range: ", 
              *std::min_element(initial_temperatures.begin(), initial_temperatures.end()), " to ",
              *std::max_element(initial_temperatures.begin(), initial_temperatures.end()), ")");
  } else {
    Log::Warning("Simulate", "No initial temperature data received from OpenFOAM!");
  }
  
  if (!initial_temperatures.empty()) {
    Log::Info("Simulate", "Successfully received ", initial_temperatures.size(), 
              " initial temperature values from OpenFOAM");
    
    // Apply initial temperature values to cells - using agent UID to OpenFOAM cell mapping
    const auto& cell_agent_map = adapter.GetCellAgentMap();
    
    if (!cell_agent_map.empty()) {
      Log::Info("Simulate", "Using cell-agent mapping with ", cell_agent_map.size(), " entries");
      
      // Create map from AgentUid to temperature
      std::map<AgentUid, double> agent_temperatures;
      
      for (size_t i = 0; i < cell_agent_map.size() && i < initial_temperatures.size(); i++) {
        agent_temperatures[cell_agent_map[i].first] = initial_temperatures[i];
      }
      
      // Track temperature stats
      double min_init_temp = std::numeric_limits<double>::max();
      double max_init_temp = std::numeric_limits<double>::lowest();
      double sum_init_temp = 0.0;
      int cells_initialized = 0;
      
      // Apply temperatures to agents
      rm->ForEachAgent([&](Agent* agent) {
        if (auto* cell = dynamic_cast<MyCell*>(agent)) {
          auto agent_uid = cell->GetUid();
          
          // Check if this agent has a temperature mapping
          auto it = agent_temperatures.find(agent_uid);
          if (it != agent_temperatures.end()) {
            double temp = it->second;
            cell->SetTemperature(temp);
            
            // Track temperature stats
            min_init_temp = std::min(min_init_temp, temp);
            max_init_temp = std::max(max_init_temp, temp);
            sum_init_temp += temp;
            cells_initialized++;
            
            // Set color based on temperature (blue to red)
            double norm_temp = (temp - 300.0) / 150.0; // Normalize to 0-1 range
            Double3 color = {
              std::min(1.0, std::max(0.0, norm_temp)),     // Red
              0.0,                                         // Green
              std::min(1.0, std::max(0.0, 1.0 - norm_temp)) // Blue
            };
            cell->SetCellColor(color);
            
            // Log a few examples
            if (cells_initialized < 5 || cells_initialized % 200 == 0) {
              const auto& pos = cell->GetPosition();
              Log::Info("Simulate", "Cell at (", pos[0], ",", pos[1], ",", pos[2],
                      ") initialized with temperature ", temp);
            }
          }
        }
      });
      
      // Log initialization statistics
      if (cells_initialized > 0) {
        double avg_temp = sum_init_temp / cells_initialized;
        Log::Info("Simulate", "Successfully initialized ", cells_initialized, 
                  " agent temperatures from OpenFOAM data");
        Log::Info("Simulate", "Temperature stats - Min: ", min_init_temp, 
                  ", Max: ", max_init_temp, ", Avg: ", avg_temp);
      } else {
        Log::Warning("Simulate", "No cells were initialized with temperature data");
      }
    } else {
      Log::Warning("Simulate", "Cell-agent mapping is empty!");
    }
  } else {
    Log::Warning("Simulate", "Failed to receive initial temperature data from OpenFOAM. ",
                "Using default temperature value (", initial_temp, ")");
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
      std::cout << "TIMESTEP " << timestep << ": Received " << temperatures.size() 
                << " temperature values" << std::endl;
      
      // Track temperature ranges
      double min_temp = std::numeric_limits<double>::max();
      double max_temp = std::numeric_limits<double>::lowest();
      
      for (size_t i = 0; i < temperatures.size(); i++) {
        min_temp = std::min(min_temp, temperatures[i]);
        max_temp = std::max(max_temp, temperatures[i]);
      }
      
      std::cout << "TIMESTEP " << timestep << ": Temperature range [" 
                << min_temp << ", " << max_temp << "]" << std::endl;
      
      // Check if values are changing between timesteps
      static double prev_max_temp = 0;
      static double prev_min_temp = 0;
      
      if (timestep > 1) {
        std::cout << "TIMESTEP " << timestep << ": Temperature change: Min delta = " 
                  << (min_temp - prev_min_temp) << ", Max delta = "
                  << (max_temp - prev_max_temp) << std::endl;
      }
      
      prev_min_temp = min_temp;
      prev_max_temp = max_temp;
      
      // Apply temperature values to cells - using agent UID to OpenFOAM cell mapping
      const auto& cell_agent_map = adapter.GetCellAgentMap();
      
      if (!cell_agent_map.empty()) {
        // Create map from AgentUid to temperature
        std::map<AgentUid, double> agent_temperatures;
        
        for (size_t i = 0; i < cell_agent_map.size() && i < temperatures.size(); i++) {
          agent_temperatures[cell_agent_map[i].first] = temperatures[i];
        }
        
        // Track temperature stats for agents
        double min_applied_temp = std::numeric_limits<double>::max();
        double max_applied_temp = std::numeric_limits<double>::lowest();
        double sum_applied_temp = 0.0;
        int cells_updated = 0;
        
        // Apply temperatures to agents
        rm->ForEachAgent([&](Agent* agent) {
          if (auto* cell = dynamic_cast<MyCell*>(agent)) {
            auto agent_uid = cell->GetUid();
            
            // Check if this agent has a temperature mapping
            auto it = agent_temperatures.find(agent_uid);
            if (it != agent_temperatures.end()) {
              double temp = it->second;
              double old_temp = cell->GetTemperature();
              cell->SetTemperature(temp);
              
              // Track temperature stats
              min_applied_temp = std::min(min_applied_temp, temp);
              max_applied_temp = std::max(max_applied_temp, temp);
              sum_applied_temp += temp;
              cells_updated++;
              
              // Set color based on temperature (blue to red)
              double norm_temp = (temp - 300.0) / 150.0; // Normalize to 0-1 range
              Double3 color = {
                std::min(1.0, std::max(0.0, norm_temp)),     // Red
                0.0,                                         // Green
                std::min(1.0, std::max(0.0, 1.0 - norm_temp)) // Blue
              };
              cell->SetCellColor(color);
              
              // Only log a small sample of cells
              if (cells_updated == 1 || cells_updated % 200 == 0) {
                const auto& pos = cell->GetPosition();
                std::cout << "TIMESTEP " << timestep << ": Cell at (" 
                          << pos[0] << ", " << pos[1] << ", " << pos[2]
                          << ") temperature changed from " << old_temp 
                          << " to " << temp << " (delta: " 
                          << temp - old_temp << ")" << std::endl;
              }
            }
          }
        });
        
        // Log temperature application statistics
        if (cells_updated > 0) {
          double avg_temp = sum_applied_temp / cells_updated;
          std::cout << "TIMESTEP " << timestep << ": Temperature stats - Min: " 
                    << min_applied_temp << ", Max: " << max_applied_temp
                    << ", Avg: " << avg_temp << ", Cells updated: " 
                    << cells_updated << std::endl;
        } else {
          std::cout << "TIMESTEP " << timestep << ": WARNING - No cells were updated with temperature data" << std::endl;
        }
      } else {
        std::cout << "TIMESTEP " << timestep << ": WARNING - Cell-agent mapping is empty!" << std::endl;
      }
    } else {
      std::cout << "TIMESTEP " << timestep << ": No temperature data received" << std::endl;
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