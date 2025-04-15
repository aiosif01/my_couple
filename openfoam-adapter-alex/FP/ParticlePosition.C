#include "ParticlePosition.H"
#include "Utilities.H"
#include "volFields.H"
#include "fvMesh.H"
#include "cellSet.H"

using namespace Foam;

namespace preciceAdapter
{
namespace FP
{

// Constructor
ParticlePosition::ParticlePosition(const Foam::fvMesh& mesh)
: mesh_(mesh)
{
    dataType_ = vector; // Positions are vectors
    DEBUG(adapterInfo("ParticlePosition: Constructed.", "debug"));
}

// Read implementation: Copy positions from buffer into OpenFOAM (or use directly)
void ParticlePosition::read(double* dataBuffer, const unsigned int dim)
{
    DEBUG(adapterInfo("ParticlePosition: Reading data...", "debug"));
    if (dataType_ != vector || dim != 3) {
        adapterInfo("ParticlePosition::read - ERROR: Expecting 3D vector data.", "error");
        return;
    }

    int bufferIndex = 0;
    int vectorCount = 0; // Keep track of how many vectors we read

    // Logic for volumeCenters coupling
    if (this->locationType_ == LocationType::volumeCenters)
    {
        if (cellSetNames_.empty()) // Entire volume
        {
            // IMPORTANT: Need a reliable way to know how many vectors (positions) to read.
            // Reading based on mesh_.nCells() is likely WRONG. The number of coupling
            // points comes from the preCICE setup (number of vertices in the mesh defined
            // in Interface::configureMesh). This value isn't directly available here.
            // For now, we assume the buffer is correctly sized by preCICE.
            
            // Placeholder: Read just one vector for demonstration if count is unknown
            // Replace this logic with a proper loop when count is known.
            if (true) { // Replace 'true' with condition checking if buffer has data
                 Foam::vector pos;
                 pos.x() = dataBuffer[bufferIndex++];
                 pos.y() = dataBuffer[bufferIndex++];
                 pos.z() = dataBuffer[bufferIndex++];
                 vectorCount++;
                 DEBUG(adapterInfo("  Read Pos (point 0): (" + 
                                 std::to_string(pos.x()) + ", " + 
                                 std::to_string(pos.y()) + ", " + 
                                 std::to_string(pos.z()) + ")", "debug"));
            }
            adapterInfo("ParticlePosition::read - Loop for entire volume needs correct coupling point count!", "error");
        }
        else // Specific cellSets
        {
            DEBUG(adapterInfo("ParticlePosition: Reading positions for specific cellSet(s).", "debug"));
            for (const auto& cellSetName : cellSetNames_)
            {
                try {
                    // Create cell set - this will throw if it doesn't exist
                    cellSet region(mesh_, cellSetName);
                    
                    // If cell set is empty, skip
                    if (region.empty()) {
                        DEBUG(adapterInfo("ParticlePosition: CellSet '" + cellSetName + "' is empty. Skipping.", "debug"));
                        continue;
                    }
                    
                    const labelList& cells = region.toc();
                    DEBUG(adapterInfo("ParticlePosition: CellSet '" + cellSetName + 
                                     "' contains " + std::to_string(cells.size()) + " cells.", "debug"));
                    
                    // Loop through the cells using index-based loop
                    for (label i = 0; i < cells.size(); i++)
                    {
                        // Safety check for buffer access
                        if (bufferIndex + 2 >= dim * mesh_.nCells()) {
                            adapterInfo("ParticlePosition::read - ERROR: Buffer index would exceed buffer size.", "error");
                            return;
                        }
                        
                        Foam::vector pos;
                        pos.x() = dataBuffer[bufferIndex++];
                        pos.y() = dataBuffer[bufferIndex++];
                        pos.z() = dataBuffer[bufferIndex++];
                        vectorCount++;

                        DEBUG(adapterInfo("  Read Pos #" + std::to_string(i) + " for cell " + 
                                       std::to_string(cells[i]) + " in set '" + cellSetName + "': (" + 
                                       std::to_string(pos.x()) + ", " + 
                                       std::to_string(pos.y()) + ", " + 
                                       std::to_string(pos.z()) + ")", "debug"));
                    }
                }
                catch (const Foam::error& e) {
                    adapterInfo("ParticlePosition::read - ERROR: CellSet '" + cellSetName + 
                               "' not found: " + e.message(), "error");
                }
            }
        }
    } // End volumeCenters logic

    // --- Handle Boundary Patches (if applicable, likely not for this scenario) ---
    if (!patchIDs_.empty()) {
        DEBUG(adapterInfo("ParticlePosition: Reading positions for boundary patches (if any).", "debug"));
        adapterInfo("ParticlePosition::read - Boundary patch handling for positions not implemented.", "warning");
    }

    DEBUG(adapterInfo("ParticlePosition: Finished reading data. Read " + std::to_string(vectorCount) + 
                    " vectors (" + std::to_string(bufferIndex) + " doubles).", "debug"));
}

// Write implementation (No-Op for this class)
std::size_t ParticlePosition::write(double* dataBuffer, bool meshConnectivity, const unsigned int dim)
{
    // OpenFOAM reads ParticlePosition, does not write it.
    DEBUG(adapterInfo("ParticlePosition: Write called (NO-OP).", "debug"));
    return 0; // Return 0 elements written
}

// Support check for location type
bool ParticlePosition::isLocationTypeSupported(const bool meshConnectivity) const
{
    // Designed for volume coupling
    return locationType_ == LocationType::volumeCenters;
}

// Get the data name
std::string ParticlePosition::getDataName() const
{
    // MUST match preciceDict and FP::addReaders
    return "ParticlePosition";
}

} // namespace FP
} // namespace preciceAdapter