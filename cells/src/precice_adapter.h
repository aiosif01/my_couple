#ifndef PRECICE_ADAPTER_H
#define PRECICE_ADAPTER_H

#include "biodynamo.h"
#include "precice/precice.hpp"

namespace bdm {

class PreciceAdapter {
 public:
  PreciceAdapter(const std::string& config_file, const std::string& participant_name)
      : interface_(participant_name, config_file, 0, 1),
        mesh_name_("CellMesh") {}

  void Initialize() {
    interface_.initialize();
  }

  void UpdateMesh(Simulation& simulation) {
    auto* rm = simulation.GetResourceManager();
    positions_.clear();

    for (size_t i = 0; i < rm->GetNumAgents(); ++i) {
      auto* agent = rm->GetAgent(AgentHandle(i));
      if (auto* cell = dynamic_cast<Cell*>(agent)) {
        const auto& pos = cell->GetPosition();
        positions_.insert(positions_.end(), {pos[0], pos[1], pos[2]});
      }
    }

    size_t num_vertices = positions_.size() / 3;
    vertex_ids_.resize(num_vertices);

    interface_.setMeshVertices(
        mesh_name_,
        precice::span<const double>(positions_.data(), positions_.size()),
        precice::span<int>(vertex_ids_.data(), vertex_ids_.size()));
  }

  void ReadTemperature(std::vector<double>& temperatures) {
    temperatures.resize(vertex_ids_.size());

    interface_.readData(
        mesh_name_,
        "T",
        precice::span<const int>(vertex_ids_.data(), vertex_ids_.size()),
        0.0,  // relativeReadTime
        precice::span<double>(temperatures.data(), temperatures.size()));
  }

  void Advance(double dt) {
    interface_.advance(dt);
  }

  bool IsCouplingOngoing() {
    return interface_.isCouplingOngoing();
  }

  void Finalize() {
    interface_.finalize();
  }

 private:
  precice::Participant interface_;
  std::string mesh_name_;
  std::vector<double> positions_;
  std::vector<int> vertex_ids_;
};

}  // namespace bdm

#endif  // PRECICE_ADAPTER_H