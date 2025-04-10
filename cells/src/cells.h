#ifndef CELLS_H_
#define CELLS_H_

#include "biodynamo.h"
#include "precice_adapter.h"

namespace bdm {

inline int Simulate(int argc, const char** argv) {
  Simulation simulation(argc, argv);

  PreciceAdapter adapter("../precice-config.xml", "cells");
  
  // Define domain bounds matching OpenFOAM's 1x1x1 cube
  double min_bound = 0.0;
  double max_bound = 1.0;

  auto construct = [](const Real3& position) {
    Cell* cell = new Cell(position);
    cell->SetDiameter(0.01);
    return cell;
  };

  // Create agents within the domain
  ModelInitializer::CreateAgentsRandom(min_bound, max_bound, 10, construct);

  // Update mesh with initial agent positions BEFORE initialization
  adapter.UpdateMesh(simulation);

  // Initialize preCICE
  adapter.Initialize();

  double dt = 0.01;

  while (adapter.IsCouplingOngoing()) {
    // Update mesh with current positions (if agents move)
    // adapter.UpdateMesh(simulation);

    // Read temperatures from preCICE
    std::vector<double> temperatures;
    adapter.ReadTemperature(temperatures);

    // Access agents and print temperatures
    auto* rm = simulation.GetResourceManager();
    size_t temp_idx = 0;
    for (size_t i = 0; i < rm->GetNumAgents(); ++i) {
      auto* agent = rm->GetAgent(AgentHandle(i));
      if (auto* cell = dynamic_cast<Cell*>(agent)) {
        if (temp_idx < temperatures.size()) {
          double temp = temperatures[temp_idx++];
          std::cout << "Cell at " << cell->GetPosition() << " has T = " << temp << "\n";
        }
      }
    }

    // Simulate one step (may update agent positions)
    simulation.GetScheduler()->Simulate(1);

    // Advance preCICE
    adapter.Advance(dt);
  }

  adapter.Finalize();
  std::cout << "Simulation completed successfully!\n";
  return 0;
}

}  // namespace bdm

#endif  // CELLS_H_