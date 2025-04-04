#ifndef PARTICLEPOSITION_H
#define PARTICLEPOSITION_H

#include "CouplingDataUser.H" // Base class
#include "fvMesh.H"
#include "volFields.H"      // For volVectorField if using temporary storage

namespace preciceAdapter
{
namespace FP
{

class ParticlePosition : public CouplingDataUser
{
private:
    const Foam::fvMesh& mesh_;
    // Optional: Pointer to store positions if needed persistently between calls
    // Foam::volVectorField* positionField_;
    // const std::string positionFieldName_ = "preciceParticlePositions"; // Example name

public:
    // Constructor
    ParticlePosition(const Foam::fvMesh& mesh);

    // Destructor
    ~ParticlePosition() override = default;

    // --- Implement virtual methods from CouplingDataUser ---

    // Read positions FROM the buffer (sent by BioDynaMo)
    void read(double* dataBuffer, const unsigned int dim) override;

    // Write positions TO the buffer (Does nothing here as OF reads positions)
    std::size_t write(double* dataBuffer, bool meshConnectivity, const unsigned int dim) override;

    // Check if volumeCenters location is supported
    bool isLocationTypeSupported(const bool meshConnectivity) const override;

    // Get the data name string ("ParticlePosition")
    std::string getDataName() const override;

    // Optional: Initialize temporary storage field if needed
    // void initialize() override;
};

} // namespace FP
} // namespace preciceAdapter

#endif // PARTICLEPOSITION_H

