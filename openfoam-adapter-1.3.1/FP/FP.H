#ifndef FP_H
#define FP_H

#include "fvMesh.H"
#include "dictionary.H"
#include "Interface.H" // Adjust path if necessary

namespace preciceAdapter
{
class Interface; // Forward declaration

namespace FP
{
class FluidParticle
{
private:
    const Foam::fvMesh& mesh_;
    // Add any FP-specific config members here if needed
    std::string nameT_ = "T"; // Default name for Temperature field

public:
    FluidParticle(const Foam::fvMesh& mesh);
    ~FluidParticle();
    bool configure(const Foam::dictionary& dict);
    bool addWriters(std::string dataName, Interface* interface);
    bool addReaders(std::string dataName, Interface* interface);
};
} // namespace FP
} // namespace preciceAdapter
#endif // FP_H

