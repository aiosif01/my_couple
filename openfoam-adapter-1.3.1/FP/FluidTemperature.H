#ifndef FLUIDTEMPERATURE_H
#define FLUIDTEMPERATURE_H

#include "CouplingDataUser.H" // Base class
#include "fvMesh.H"
#include "volFields.H" // For volScalarField

namespace preciceAdapter
{
namespace FP
{
class FluidTemperature : public CouplingDataUser
{
private:
    // Pointer to the OpenFOAM Temperature field
    const Foam::volScalarField* T_;

public:
    // Constructor: Takes mesh and the name of the temperature field
    FluidTemperature(const Foam::fvMesh& mesh, const std::string nameT);

    // Destructor (usually empty unless allocating memory here)
    ~FluidTemperature() override = default;

    // --- Implement virtual methods from CouplingDataUser ---

    // Write temperature data FROM OpenFOAM TO buffer
    std::size_t write(double* dataBuffer, bool meshConnectivity, const unsigned int dim) override;

    // Read method (NO-OP as OpenFOAM writes this data)
    void read(double* dataBuffer, const unsigned int dim) override;

    // Check if volumeCenters location is supported
    bool isLocationTypeSupported(const bool meshConnectivity) const override;

    // Get the data name string ("FluidTemperature")
    std::string getDataName() const override;

    // Initialize (likely not needed)
    // void initialize() override;
};
} // namespace FP
} // namespace preciceAdapter
#endif // FLUIDTEMPERATURE_H

