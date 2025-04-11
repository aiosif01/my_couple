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
    
    // Configuration parameters
    std::string nameT_ = "T"; // Default name for Temperature field
    
    // Status tracking for better debugging
    bool isConfigured_ = false;
    bool tempFieldExists_ = false;

public:
    FluidParticle(const Foam::fvMesh& mesh);
    ~FluidParticle();
    
    // Configuration method with enhanced validation
    bool configure(const Foam::dictionary& dict);
    
    // Writer/reader setup with additional validation
    bool addWriters(std::string dataName, Interface* interface);
    bool addReaders(std::string dataName, Interface* interface);
    
    // Validation helpers
    bool validateTemperatureField() const;
    
    // Getters for diagnostics
    const std::string& getTemperatureFieldName() const { return nameT_; }
    bool isTemperatureFieldValid() const { return tempFieldExists_; }
};
} // namespace FP
} // namespace preciceAdapter
#endif // FP_H