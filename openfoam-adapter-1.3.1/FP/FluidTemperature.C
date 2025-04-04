#include "FluidTemperature.H" // Use correct header name if different
#include "Utilities.H"
#include "volFields.H"
#include "fvMesh.H"      // Include fvMesh.H explicitly if needed beyond volFields.H
#include "cellSet.H"     // For handling cellSets

using namespace Foam; // Use Foam namespace

namespace preciceAdapter
{
namespace FP
{

// Constructor
FluidTemperature::FluidTemperature(
    const Foam::fvMesh& mesh, // Use Foam::fvMesh explicitly
    const std::string nameT)
: T_(nullptr) // Initialize pointer
{
    dataType_ = scalar; // Temperature is scalar
    DEBUG(adapterInfo("FluidTemperature: Constructing...", "debug"));

    // Lookup the temperature field in the mesh object registry
    if (mesh.foundObject<volScalarField>(nameT))
    {
        T_ = &mesh.lookupObject<volScalarField>(nameT); // Get address of looked-up object
        DEBUG(adapterInfo("FluidTemperature: Found temperature field '" + nameT + "'.", "debug"));
    }
    else
    {
        // If field not found, it's a fatal error for writing
        // Use the adapterInfo utility which handles error reporting
        adapterInfo("FluidTemperature: ERROR - Could not find volScalarField '" + nameT + "'.", "error");
        // adapterInfo throws and potentially exits, so T_ check later might be redundant
        // but doesn't hurt.
    }
}

// Write implementation: Copy T from OpenFOAM to buffer
std::size_t FluidTemperature::write(double* dataBuffer, bool meshConnectivity, const unsigned int dim)
{
    DEBUG(adapterInfo("FluidTemperature: Writing data...", "debug"));
    if (!T_) {
        adapterInfo("FluidTemperature::write - ERROR: Temperature field pointer is null. Cannot write.", "error");
        return 0; // Return 0 elements written on error
    }
     // Check data type (optional sanity check)
    if (dataType_ != scalar) {
         adapterInfo("FluidTemperature::write - ERROR: Expecting scalar data.", "error");
         return 0;
    }


    int bufferIndex = 0;

    // --- Handle Volume Coupling Data ---
    if (this->locationType_ == LocationType::volumeCenters)
    {
        const volScalarField::Internal& field = T_->internalField();

        if (cellSetNames_.empty()) // Couple entire internal field
        {
            DEBUG(adapterInfo("FluidTemperature: Writing temperature for entire volume.", "debug"));
            // Use indexed loop
            for (label cellI = 0; cellI < field.size(); ++cellI)
            {
                dataBuffer[bufferIndex++] = field[cellI]; // Terminate statement
            }
        }
        else // Couple only specified cellSets
        {
            DEBUG(adapterInfo("FluidTemperature: Writing temperature for specific cellSet(s).", "debug"));
            for (const auto& cellSetName : cellSetNames_)
            {
                 // Access mesh via the field T_
                cellSet region(T_->mesh(), cellSetName);
                const labelList& cells = region.toc();
                for (const label currentCell : cells)
                {
                    // Ensure cell index is valid before accessing
                    if (currentCell >= 0 && currentCell < field.size()) { // Check lower bound too
                        dataBuffer[bufferIndex++] = field[currentCell]; // Terminate statement
                    } else {
                        adapterInfo("FluidTemperature::write - WARNING: Invalid cell index "
                                    + std::to_string(currentCell) + " from cellSet '"
                                    + cellSetName + "'. Skipping.", "warning");
                    }
                }
            }
        }
    }

    // --- Handle Boundary Patches (if specified in preciceDict) ---
    if (!patchIDs_.empty())
    {
        DEBUG(adapterInfo("FluidTemperature: Writing temperature for boundary patches.", "debug"));
        for (int patchID : patchIDs_)
        {
            // Get reference to the boundary field for this patch
            const fvPatchScalarField& patchField = T_->boundaryField()[patchID];
            // Use indexed loop for boundary field
            for (label faceI = 0; faceI < patchField.size(); ++faceI)
            {
                dataBuffer[bufferIndex++] = patchField[faceI]; // Terminate statement
            }
        }
    }

    DEBUG(adapterInfo("FluidTemperature: Finished writing " + std::to_string(bufferIndex) + " values.", "debug"));
    return bufferIndex; // Return number of elements written
}

// Read implementation (No-Op for this class)
void FluidTemperature::read(double* dataBuffer, const unsigned int dim)
{
    // OpenFOAM writes FluidTemperature, does not read it.
    DEBUG(adapterInfo("FluidTemperature: Read called (NO-OP).", "debug"));
    // No return statement needed for void function
}

// Support check for location type
bool FluidTemperature::isLocationTypeSupported(const bool meshConnectivity) const
{
    // Allow volume coupling and face center coupling
    return (locationType_ == LocationType::volumeCenters || locationType_ == LocationType::faceCenters);
}

// Get the data name
std::string FluidTemperature::getDataName() const
{
    // Ensure this matches the name used in preciceDict and FP::addWriters
    return "Temperature"; // Or just "FluidTemperature" if preferred
}

// Initialize method (likely not needed for reading a standard field)
// void FluidTemperature::initialize() { }

} // namespace FP
} // namespace preciceAdapter

