#ifndef PRECICE_ADAPTER_H_
#define PRECICE_ADAPTER_H_

#include "biodynamo.h"
#include "precice/precice.hpp" // Includes span definition
#include <vector>
#include <string>

#include "my_cell.h" // Include the full definition of MyCell

namespace bdm {

class PreciceAdapter {
 public:
  PreciceAdapter(const std::string& config_file, const std::string& participant_name)
      : interface_(participant_name, config_file, 0, 1),
        mesh_name_("CellMesh"),
        temperature_data_name_("T") {
    Log::Info("PreciceAdapter", "Adapter created for participant: ", participant_name);
  }

  void Initialize() {
    Log::Info("PreciceAdapter", "Initializing preCICE interface...");
    interface_.initialize();
    Log::Info("PreciceAdapter", "preCICE interface initialized successfully.");
  }

  void UpdateMesh(Simulation& simulation) {
    // We can't check if preCICE is initialized in this version of the API
    // Instead, we'll use a flag that we maintain ourselves
    static bool meshAlreadySet = false;
    
    if (meshAlreadySet) {
      // Skip mesh modification after first setup
      Log::Info("PreciceAdapter", "UpdateMesh: Mesh already registered with preCICE, skipping modification");
      return;
    }
    
    // Continue with mesh setup
    Log::Info("PreciceAdapter", "UpdateMesh: Setting up mesh vertices...");
    auto* rm = simulation.GetResourceManager();
    positions_.clear();
    vertex_ids_.clear();

    // Count MyCell agents for logging
    int cellCount = 0;
    rm->ForEachAgent([&](Agent* agent) {
      if (dynamic_cast<MyCell*>(agent)) {
        cellCount++;
      }
    });
    
    Log::Info("PreciceAdapter", "UpdateMesh: Found ", cellCount, " MyCell agents");

    // Collect positions of all MyCell agents
    rm->ForEachAgent([&](Agent* agent) {
      if (auto* my_cell = dynamic_cast<MyCell*>(agent)) {
        const auto& pos = my_cell->GetPosition();
        // Store x, y, z coordinates for each cell
        positions_.push_back(pos[0]);
        positions_.push_back(pos[1]);
        positions_.push_back(pos[2]);
        
        Log::Info("PreciceAdapter", "UpdateMesh: Added cell at position (", 
                 pos[0], ", ", pos[1], ", ", pos[2], ")");
      }
    });

    // Calculate number of vertices from positions array (3 coords per vertex)
    size_t num_vertices = positions_.size() / 3;

    if (num_vertices == 0) {
      Log::Warning("PreciceAdapter", "UpdateMesh: No MyCell agents found.");
      return; // Exit early if no vertices
    }

    // Resize vertex IDs array to match number of vertices
    vertex_ids_.resize(num_vertices);
    
    // Register the mesh vertices with preCICE
    Log::Info("PreciceAdapter", "UpdateMesh: Registering ", num_vertices, " vertices with preCICE");
    
    interface_.setMeshVertices(
        mesh_name_,
        precice::span<const double>(positions_.data(), positions_.size()), 
        precice::span<int>(vertex_ids_.data(), vertex_ids_.size()));
    
    Log::Info("PreciceAdapter", "UpdateMesh: Successfully registered ", num_vertices, " vertices with preCICE");
    
    // Mark that mesh has been set
    meshAlreadySet = true;
  }

  void ReadTemperature(std::vector<double>& temperatures) {
     size_t num_vertices = vertex_ids_.size();

     if (num_vertices == 0) {
         // Use std::cout directly for critical debugging - will be visible regardless of log level
         std::cout << "CRITICAL DEBUG: No vertices registered with preCICE, cannot read temperature data" << std::endl;
         temperatures.clear();
         return;
     }

     temperatures.resize(num_vertices); // Ensure output buffer is correctly sized
     
     // Use direct console output for critical debugging
     std::cout << "CRITICAL DEBUG: Attempting to read temperature data for " << num_vertices << " vertices" << std::endl;

     // *** Use the 5-argument span-based readData ***
     double relative_read_time = 0.0; // Usually 0.0 for current time step

     try {
         interface_.readData(
             mesh_name_,                                                      // 1. Mesh Name
             temperature_data_name_,                                          // 2. Data Name
             precice::span<const int>(vertex_ids_.data(), vertex_ids_.size()), // 3. Vertex IDs
             relative_read_time,                                              // 4. Relative Read Time
             precice::span<double>(temperatures.data(), temperatures.size())); // 5. Output Data

         // Enhanced debugging: Print temperature values with direct console output
         if (!temperatures.empty()) {
             std::cout << "CRITICAL DEBUG: Successfully read " << temperatures.size() << " temperature values" << std::endl;
             std::cout << "CRITICAL DEBUG: First few temperature values: ";
             
             // Print first few temperatures (max 5) with detailed info
             int max_to_print = std::min(static_cast<size_t>(5), temperatures.size());
             for (int i = 0; i < max_to_print; i++) {
                 std::cout << temperatures[i];
                 if (i < max_to_print - 1) std::cout << ", ";
             }
             std::cout << std::endl;
             
             // Also calculate min/max/avg for easy verification
             if (!temperatures.empty()) {
                 double min_temp = *std::min_element(temperatures.begin(), temperatures.end());
                 double max_temp = *std::max_element(temperatures.begin(), temperatures.end());
                 double sum = std::accumulate(temperatures.begin(), temperatures.end(), 0.0);
                 double avg_temp = sum / temperatures.size();
                 
                 std::cout << "CRITICAL DEBUG: Temperature stats - Min: " << min_temp 
                           << ", Max: " << max_temp << ", Avg: " << avg_temp << std::endl;
             }
         } else {
             std::cout << "CRITICAL DEBUG: Temperature data array is empty after readData!" << std::endl;
         }
     } catch (const std::exception& e) {
         std::cout << "CRITICAL DEBUG: Exception reading temperature data: " << e.what() << std::endl;
         temperatures.clear();
     }
  }

  // --- NO CHANGES BELOW THIS LINE ---
  void Advance(double dt) {
    interface_.advance(dt);
  }

  bool IsCouplingOngoing() {
    bool ongoing = interface_.isCouplingOngoing();
    return ongoing;
  }

  double GetMaxTimeStep() {
      double dt = interface_.getMaxTimeStepSize();
      return dt;
  }

  void Finalize() {
    Log::Info("PreciceAdapter", "Finalizing preCICE interface...");
    try {
        interface_.finalize();
        Log::Info("PreciceAdapter", "preCICE interface finalized.");
    } catch (const std::exception& e) {
        Log::Error("PreciceAdapter", "Exception during preCICE finalize: ", e.what());
    }
  }

 private:
  precice::Participant interface_;
  std::string mesh_name_;
  std::string temperature_data_name_;
  std::vector<double> positions_;
  std::vector<int> vertex_ids_;
};

} // namespace bdm

#endif // PRECICE_ADAPTER_H_