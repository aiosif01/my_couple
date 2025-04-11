#ifndef MY_CELL_H_
#define MY_CELL_H_

#include "biodynamo.h"

namespace bdm {

// Define a custom cell type MyCell that inherits from bdm::Cell
class MyCell : public Cell {
  BDM_AGENT_HEADER(MyCell, Cell, 1); // Macro needed for BioDynaMo's agent system

 private: // Keep internal data private
  double temperature_ = 0.0; // Initialize temperature, default to 0 or an expected initial value
  Double3 cell_color_ = {0.0, 0.0, 1.0}; // Default blue color (RGB), needed for visualization

 public:
  MyCell() : Base() {}
  explicit MyCell(const Real3& position) : Base(position) {}

  // Method to set the temperature (e.g., from preCICE data)
  void SetTemperature(double temp) { temperature_ = temp; }

  // Method to get the stored temperature (e.g., for use in behaviors)
  double GetTemperature() const { return temperature_; }

  // Methods for cell coloring
  void SetCellColor(const Double3& color) { cell_color_ = color; }
  const Double3& GetCellColor() const { return cell_color_; }

  // Add any other custom properties or behaviors specific to MyCell later.
};

} // namespace bdm

#endif // MY_CELL_H_