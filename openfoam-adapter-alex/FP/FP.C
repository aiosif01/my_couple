#include "FP.H"
#include "Interface.H"
#include "Utilities.H"

// Include headers for FP data handlers
#include "FluidTemperature.H"
#include "ParticlePosition.H"

using namespace Foam;

namespace preciceAdapter
{
namespace FP
{
FluidParticle::FluidParticle(const Foam::fvMesh& mesh) 
: mesh_(mesh),
  isConfigured_(false),
  tempFieldExists_(false)
{
    DEBUG(adapterInfo("FP Module: Constructed", "debug"));
}

FluidParticle::~FluidParticle()
{
    DEBUG(adapterInfo("FP Module: Destroyed", "debug"));
}

bool FluidParticle::validateTemperatureField() const
{
    bool exists = mesh_.foundObject<volScalarField>(nameT_);
    
    if (exists) {
        // Field exists, let's check its properties
        const volScalarField& field = mesh_.lookupObject<volScalarField>(nameT_);
        
        DEBUG(adapterInfo("FP Module: Temperature field validation - " + nameT_ + 
                         " exists with " + std::to_string(field.internalField().size()) + 
                         " internal cells", "debug"));
        
        // Sample values for debugging if available
        if (field.internalField().size() > 0) {
            DEBUG(adapterInfo("FP Module: Temperature field sample value: " + 
                             std::to_string(field.internalField()[0]), "debug"));
        }
        
        return true;
    } else {
        DEBUG(adapterInfo("FP Module: Temperature field '" + nameT_ + "' NOT found in mesh registry", "debug"));
        return false;
    }
}

bool FluidParticle::configure(const Foam::dictionary& dict)
{
    DEBUG(adapterInfo("FP Module: Beginning configuration...", "debug"));
    
    // Find the FP sub-dictionary
    const dictionary* FPDict = dict.findDict("FP");
    if (FPDict) {
        DEBUG(adapterInfo("FP Module: Found 'FP' sub-dictionary in preciceDict", "debug"));
        
        // Read temperature field name
        if (FPDict->found("nameT")) {
            nameT_ = FPDict->get<word>("nameT");
            DEBUG(adapterInfo("FP Module: Using user-specified Temperature field name: '" + nameT_ + "'", "debug"));
        } else {
            DEBUG(adapterInfo("FP Module: Using default Temperature field name: '" + nameT_ + "'", "debug"));
        }
        
        // Read other FP-specific parameters if needed
        // ...
    } else {
        DEBUG(adapterInfo("FP Module: No 'FP' sub-dictionary found in preciceDict. Using defaults.", "debug"));
    }
    
    // Validate temperature field existence
    tempFieldExists_ = validateTemperatureField();
    
    if (!tempFieldExists_) {
        adapterInfo("FP Module: WARNING - Temperature field '" + nameT_ + 
                   "' not found. Temperature coupling might not work correctly.", "warning");
    }
    
    isConfigured_ = true;
    DEBUG(adapterInfo("FP Module: Configuration completed. Status: " + 
                     std::string(tempFieldExists_ ? "OK" : "WARNINGS"), "debug"));
    
    return true; // Continue even with warnings
}

bool FluidParticle::addWriters(std::string dataName, Interface* interface)
{
    DEBUG(adapterInfo("FP Module: Checking writer for data: '" + dataName + "'", "debug"));
    
    if (!isConfigured_) {
        adapterInfo("FP Module: WARNING - Module not configured before adding writers", "warning");
    }
    
    bool found = false;

    if (dataName == "T") // Temperature data that OpenFOAM sends
    {
        // Check if the field exists before creating the writer
        if (mesh_.foundObject<volScalarField>(nameT_)) {
            DEBUG(adapterInfo("FP Module: Temperature field '" + nameT_ + "' found, adding writer", "debug"));
            
            interface->addCouplingDataWriter(dataName, new FluidTemperature(mesh_, nameT_));
            DEBUG(adapterInfo("FP Module: Successfully added writer for Temperature field '" + nameT_ + "'", "debug"));
            found = true;
        } else {
            adapterInfo("FP Module: ERROR - Cannot add Temperature writer because field '" + 
                       nameT_ + "' does not exist", "error");
            
            // Additional debug info to help troubleshoot
            DEBUG(adapterInfo("FP Module: Available fields in registry:", "debug"));
            // This would need to get registry objects, but we'll leave as is for simplicity
        }
    }
    // Add other data writers here with 'else if' blocks
    
    return found;
}

bool FluidParticle::addReaders(std::string dataName, Interface* interface)
{
    DEBUG(adapterInfo("FP Module: Checking reader for data: '" + dataName + "'", "debug"));
    
    if (!isConfigured_) {
        adapterInfo("FP Module: WARNING - Module not configured before adding readers", "warning");
    }
    
    bool found = false;

    if (dataName == "ParticlePosition") // Example - data that OpenFOAM receives
    {
        DEBUG(adapterInfo("FP Module: Adding reader for ParticlePosition", "debug"));
        interface->addCouplingDataReader(dataName, new ParticlePosition(mesh_));
        DEBUG(adapterInfo("FP Module: Successfully added reader for ParticlePosition", "debug"));
        found = true;
    }
    // Add other data readers here with 'else if' blocks
    
    return found;
}

} // namespace FP
} // namespace preciceAdapter
