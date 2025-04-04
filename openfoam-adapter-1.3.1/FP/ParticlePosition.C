#include "ParticlePosition.H" // Use correct header name if different
#include "Utilities.H"
#include "volFields.H"   // For volVectorField (if using temp storage)
#include "fvMesh.H"      // Include fvMesh.H explicitly
#include "cellSet.H"     // For handling cellSets

using namespace Foam; // Use Foam namespace

namespace preciceAdapter
{
namespace FP
{

// Constructor
ParticlePosition::ParticlePosition(const Foam::fvMesh& mesh)
: mesh_(mesh) // Store mesh reference if needed by methods
// , positionField_(nullptr) // Initialize if using temporary field
{
    dataType_ = vector; // Positions are vectors
    DEBUG(adapterInfo("ParticlePosition: Constructed.", "debug"));
}

// Optional: Initialize temporary field if using one
// void ParticlePosition::initialize() { ... create and register field ... }

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
            // A truly robust solution might query preCICE or pass the expected size.

            // Example: Loop based on an estimated/assumed count (replace with correct size!)
            // label nCouplingPoints = get_expected_coupling_point_count(); // Need a way to get this!
            // for (label i = 0; i < nCouplingPoints; ++i)

            // Alternative (less safe): Loop until buffer seems exhausted? Requires buffer size.

            // --- Minimal implementation using a placeholder ---
            // Let's assume for now we *don't* know the exact count and just log
            // This part needs refinement based on how Interface reports size.
            DEBUG(adapterInfo("ParticlePosition: Reading positions for entire volume (actual count depends on preCICE mesh definition).", "debug"));
            // Placeholder: Read just one vector for demonstration if count is unknown
            // Replace this logic with a proper loop when count is known.
            if (true) { // Replace 'true' with condition checking if buffer has data
                 Foam::vector pos;
                 pos.x() = dataBuffer[bufferIndex++];
                 pos.y() = dataBuffer[bufferIndex++];
                 pos.z() = dataBuffer[bufferIndex++];
                 vectorCount++;
                 DEBUG(adapterInfo("  Read Pos (point 0): " + pos.str(), "debug"));
            }
             adapterInfo("ParticlePosition::read - Loop for entire volume needs correct coupling point count!", "error"); // Highlight this limitation
        }
        else // Specific cellSets
        {
            DEBUG(adapterInfo("ParticlePosition: Reading positions for specific cellSet(s).", "debug"));
            for (const auto& cellSetName : cellSetNames_)
            {
                cellSet region(mesh_, cellSetName);
                const labelList& cells = region.toc();
                // Loop through the cells associated with this coupling point definition
                for (const label currentCell : cells) // Now using 'currentCell'
                {
                    Foam::vector pos;
                    pos.x() = dataBuffer[bufferIndex++];
                    pos.y() = dataBuffer[bufferIndex++];
                    pos.z() = dataBuffer[bufferIndex++];
                    vectorCount++;

                    // Use 'currentCell' in the debug message to resolve the warning
                    DEBUG(adapterInfo("  Read Pos for coupling point corresponding to cell " + std::to_string(currentCell) + " in set '" + cellSetName + "': " + pos.str(), "debug"));

                    // If you later decide to store positions, you'd use currentCell here:
                    // if(positionField_ && currentCell >= 0 && currentCell < positionField_->internalField().size()) {
                    //     positionField_->internalFieldRef()[currentCell] = pos;
                    // }
                }
            }
        }
        // if (positionField_) positionField_->correctBoundaryConditions();
    } // End volumeCenters logic

    // --- Handle Boundary Patches (if applicable, likely not for this scenario) ---
    if (!patchIDs_.empty()) {
         DEBUG(adapterInfo("ParticlePosition: Reading positions for boundary patches (if any).", "debug"));
        // ... (Loop through patches and faces, reading from buffer) ...
        // Remember to increment bufferIndex and vectorCount appropriately.
        adapterInfo("ParticlePosition::read - Boundary patch handling for positions not fully implemented.", "warning");
    }

    DEBUG(adapterInfo("ParticlePosition: Finished reading data. Read " + std::to_string(vectorCount) + " vectors (" + std::to_string(bufferIndex) + " doubles).", "debug"));
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

