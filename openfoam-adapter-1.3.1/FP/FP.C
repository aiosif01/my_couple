#include "FP.H" // Adjust path
#include "Interface.H" // Adjust path
#include "Utilities.H" // For logging

// Include headers for YOUR specific FP data handlers
#include "FluidTemperature.H" // Example
#include "ParticlePosition.H" // Example

using namespace Foam;

namespace preciceAdapter
{
namespace FP
{
FluidParticle::FluidParticle(const Foam::fvMesh& mesh) : mesh_(mesh)
{
    DEBUG(adapterInfo("FP Module: Constructed.", "debug"));
}

FluidParticle::~FluidParticle()
{
    DEBUG(adapterInfo("FP Module: Destroyed.", "debug"));
}

bool FluidParticle::configure(const Foam::dictionary& dict)
{
    DEBUG(adapterInfo("FP Module: Configuring...", "debug"));
    const dictionary* FPDict = dict.findDict("FP");
    if (FPDict) {
         DEBUG(adapterInfo("FP Module: Reading FP sub-dictionary.", "debug"));
         // Example: Allow overriding the default Temperature field name
         nameT_ = FPDict->lookupOrDefault<word>("nameT", "T");
         DEBUG(adapterInfo("FP Module: Using Temperature field name: " + nameT_, "debug"));
         // Read other FP-specific parameters if needed
    } else {
         DEBUG(adapterInfo("FP Module: No 'FP' sub-dictionary found in preciceDict. Using defaults.", "debug"));
    }
    return true; // Return false on error
}

bool FluidParticle::addWriters(std::string dataName, Interface* interface)
{
    // Called by Adapter::configure to add data *written* by OpenFOAM
    DEBUG(adapterInfo("FP Module: Checking writer for data: " + dataName, "debug"));
    bool found = false;

    if (dataName == "T") // Data OF sends
    {
        interface->addCouplingDataWriter(dataName, new FluidTemperature(mesh_, nameT_));
        DEBUG(adapterInfo("FP Module: Added writer for Temperature.", "debug"));
        found = true;
    }
    // Add 'else if' blocks for other data OpenFOAM writes

    return found;
}

bool FluidParticle::addReaders(std::string dataName, Interface* interface)
{
    // Called by Adapter::configure to add data *read* by OpenFOAM
    DEBUG(adapterInfo("FP Module: Checking reader for data: " + dataName, "debug"));
    bool found = false;

    // if (dataName == "ParticlePosition") // Data OF receives
    // {
    //      interface->addCouplingDataReader(dataName, new ParticlePosition(mesh_));
    //      DEBUG(adapterInfo("FP Module: Added reader for ParticlePosition.", "debug"));
    //      found = true;
    // }
    // Add 'else if' blocks for other data OpenFOAM reads

    return found;
}
} // namespace FP
} // namespace preciceAdapter

