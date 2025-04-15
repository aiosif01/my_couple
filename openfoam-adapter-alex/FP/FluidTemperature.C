#include "FluidTemperature.H"
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
FluidTemperature::FluidTemperature(
    const Foam::fvMesh& mesh,
    const std::string nameT)
: T_(nullptr),
  mesh_(mesh),
  fieldName_(nameT)
{
    dataType_ = scalar; // Temperature is scalar
    Info << "FP DEBUG: FluidTemperature constructor called for field '" << nameT << "'" << endl;

    // Lookup the temperature field in the mesh object registry
    if (mesh.foundObject<volScalarField>(nameT))
    {
        T_ = &mesh.lookupObject<volScalarField>(nameT);
        Info << "FP DEBUG: Successfully found temperature field '" << nameT << "'" << endl;
        
        // Log some field information for debugging
        Info << "FP DEBUG: Field dimensions: " << T_->dimensions() << endl;
        Info << "FP DEBUG: Field size: " << T_->internalField().size() << " cells" << endl;
        
        // Log first few values if available
        if (T_->internalField().size() > 0) {
            Info << "FP DEBUG: First few temperature values: ";
            int numToPrint = std::min(5, static_cast<int>(T_->internalField().size()));
            for (int i = 0; i < numToPrint; i++) {
                Info << T_->internalField()[i] << " ";
            }
            Info << endl;
        }
    }
    else
    {
        Info << "FP DEBUG: ERROR - Could not find volScalarField '" << nameT << "'" << endl;
        adapterInfo("FluidTemperature: ERROR - Could not find volScalarField '" + nameT + 
                   "'. Temperature coupling will not work correctly.", "error");
    }
}

// Helper method to validate the field
bool FluidTemperature::validateField() const 
{
    if (!T_) {
        Info << "FP DEBUG: Temperature field pointer is null" << endl;
        adapterInfo("FluidTemperature: Temperature field pointer is null.", "error");
        return false;
    }
    
    if (T_->internalField().size() == 0) {
        Info << "FP DEBUG: Temperature field has zero size" << endl;
        adapterInfo("FluidTemperature: Temperature field has zero size.", "warning");
        // Continue anyway as this might be valid in some cases
    }
    
    return true;
}

// Initialize method - perform additional validation
void FluidTemperature::initialize()
{
    Info << "FP DEBUG: FluidTemperature::initialize called" << endl;
    
    // Validate temperature field and configuration
    if (!validateField()) {
        Info << "FP DEBUG: Field validation failed during initialization" << endl;
        adapterInfo("FluidTemperature: Field validation failed during initialization.", "warning");
    }
    
    // Log coupling configuration
    Info << "FP DEBUG: Coupling configuration - location type: " << 
         (locationType_ == LocationType::faceCenters ? "faceCenters" : 
          locationType_ == LocationType::volumeCenters ? "volumeCenters" : "other") << endl;
    
    if (!patchIDs_.empty()) {
        Info << "FP DEBUG: Coupled patches: ";
        for (size_t i = 0; i < patchIDs_.size(); i++) {
            if (i > 0) Info << ", ";
            Info << patchIDs_[i];
        }
        Info << endl;
    }
    
    if (!cellSetNames_.empty()) {
        Info << "FP DEBUG: Coupled cell sets: ";
        for (size_t i = 0; i < cellSetNames_.size(); i++) {
            if (i > 0) Info << ", ";
            Info << cellSetNames_[i];
        }
        Info << endl;
        
        // Verify each cell set
        for (const auto& name : cellSetNames_) {
            try {
                cellSet cs(mesh_, name);
                Info << "FP DEBUG: CellSet '" << name << "' exists with " << cs.size() << " cells" << endl;
            } catch (const Foam::error& e) {
                Info << "FP DEBUG: ERROR accessing cellSet '" << name << "': " << e.message() << endl;
            } catch (...) {
                Info << "FP DEBUG: ERROR accessing cellSet '" << name << "': unknown error" << endl;
            }
        }
    } else {
        Info << "FP DEBUG: No cell sets specified for coupling" << endl;
    }
}

// Enhanced write implementation with better error handling and debugging
std::size_t FluidTemperature::write(double* dataBuffer, bool meshConnectivity, const unsigned int dim)
{
    Info << "FP DEBUG: FluidTemperature::write starting for field '" << fieldName_ << "'" << endl;
    Info << "FP DEBUG: T_ pointer is " << (T_ ? "valid" : "NULL") << endl;
    if (T_) {
        Info << "FP DEBUG: T field has " << T_->internalField().size() << " internal cells" << endl;
        
        // Print first few values if available
        if (T_->internalField().size() > 0) {
            Info << "FP DEBUG: First few T values: ";
            for (int i = 0; i < std::min(5, static_cast<int>(T_->internalField().size())); i++) {
                Info << T_->internalField()[i] << " ";
            }
            Info << endl;
        }
    }

    Info << "FP DEBUG: Location type is " << 
        (locationType_ == LocationType::volumeCenters ? "volumeCenters" :
         locationType_ == LocationType::faceCenters ? "faceCenters" : "other") << endl;

    Info << "FP DEBUG: Number of cellSets: " << cellSetNames_.size() << endl;
    for (const auto& name : cellSetNames_) {
        Info << "FP DEBUG: CellSet name: " << name << endl;
        
        try {
            cellSet cs(T_->mesh(), name);
            Info << "FP DEBUG: CellSet '" << name << "' exists with " << cs.size() << " cells" << endl;
        } catch (const Foam::error& e) {
            Info << "FP DEBUG: CellSet '" << name << "' error: " << e.message() << endl;
        } catch (...) {
            Info << "FP DEBUG: CellSet '" << name << "' does not exist or cannot be accessed" << endl;
        }
    }
    
    // Validate temperature field
    if (!validateField()) {
        Info << "FP DEBUG: ERROR - Temperature field invalid. Cannot write." << endl;
        adapterInfo("FluidTemperature::write - ERROR: Temperature field invalid. Cannot write.", "error");
        return 0;
    }
    
    // Check data type
    if (dataType_ != scalar) {
        Info << "FP DEBUG: ERROR - Expecting scalar data but field is not scalar." << endl;
        adapterInfo("FluidTemperature::write - ERROR: Expecting scalar data but field is not scalar.", "error");
        return 0;
    }

    int bufferIndex = 0;
    int valuesWritten = 0;

    // --- Handle Volume Coupling Data ---
    if (this->locationType_ == LocationType::volumeCenters)
    {
        const volScalarField::Internal& field = T_->internalField();
        
        if (cellSetNames_.empty()) // Couple entire internal field
        {
            Info << "FP DEBUG: Writing temperature for entire volume (" << 
                field.size() << " cells)" << endl;
            
            // Safety check
            if (field.size() == 0) {
                Info << "FP DEBUG: WARNING - Temperature field has zero internal cells." << endl;
                adapterInfo("FluidTemperature::write - WARNING: Temperature field has zero internal cells.", "warning");
            }
            
            // Use indexed loop for safety
            for (label cellI = 0; cellI < field.size(); ++cellI)
            {
                dataBuffer[bufferIndex++] = field[cellI];
                valuesWritten++;
            }
        }
        else // Couple only specified cellSets
        {
            Info << "FP DEBUG: Writing temperature for specific cellSet(s)" << endl;
            for (const auto& cellSetName : cellSetNames_)
            {
                Info << "FP DEBUG: Processing cellSet '" << cellSetName << "'" << endl;
                
                // Try to create the cellSet to check if it exists
                try {
                    cellSet region(T_->mesh(), cellSetName);
                    
                    if (region.empty()) {
                        Info << "FP DEBUG: WARNING - CellSet '" << cellSetName << 
                            "' exists but is empty. Skipping." << endl;
                        adapterInfo("FluidTemperature::write - WARNING: CellSet '" + cellSetName + 
                                   "' exists but is empty. Skipping.", "warning");
                        continue;
                    }
                    
                    Info << "FP DEBUG: CellSet '" << cellSetName << "' has " << 
                        region.size() << " cells" << endl;
                    
                    const labelList& cells = region.toc();
                    
                    // Log some cell IDs for debugging
                    if (cells.size() > 0) {
                        Info << "FP DEBUG: First few cell IDs: ";
                        for (int i = 0; i < std::min(5, static_cast<int>(cells.size())); i++) {
                            Info << cells[i] << " ";
                        }
                        Info << endl;
                    }
                    
                    for (const label currentCell : cells)
                    {
                        // Ensure cell index is valid before accessing
                        if (currentCell >= 0 && currentCell < field.size()) {
                            dataBuffer[bufferIndex++] = field[currentCell];
                            valuesWritten++;
                        } else {
                            Info << "FP DEBUG: WARNING - Invalid cell index " << 
                                currentCell << " from cellSet '" << cellSetName << "'. Skipping." << endl;
                            adapterInfo("FluidTemperature::write - WARNING: Invalid cell index "
                                        + std::to_string(currentCell) + " from cellSet '"
                                        + cellSetName + "'. Skipping.", "warning");
                        }
                    }
                } catch (const Foam::error& e) {
                    Info << "FP DEBUG: ERROR - CellSet '" << cellSetName << 
                        "' access error: " << e.message() << endl;
                    adapterInfo("FluidTemperature::write - WARNING: CellSet '" + cellSetName + 
                               "' access error: " + e.message(), "warning");
                    continue;
                }
            }
        }
    }

    // --- Handle Boundary Patches ---
    if (!patchIDs_.empty())
    {
        Info << "FP DEBUG: Writing temperature for boundary patches" << endl;
        for (int patchID : patchIDs_)
        {
            // Safety check for patch validity
            if (patchID < 0 || patchID >= T_->boundaryField().size()) {
                Info << "FP DEBUG: WARNING - Invalid patchID " << 
                    patchID << ". Skipping." << endl;
                continue;
            }
            
            // Access patch field
            const scalarField& TPatch = T_->boundaryField()[patchID];
            
            Info << "FP DEBUG: Processing patch ID " << patchID << 
                " with " << TPatch.size() << " values" << endl;
            
            // Copy values to buffer
            forAll(TPatch, i)
            {
                dataBuffer[bufferIndex++] = TPatch[i];
                valuesWritten++;
            }
        }
    }
    
    // Log statistics about the write operation
    lastWriteCount_ = valuesWritten;
    Info << "FP DEBUG: FluidTemperature::write completed with " << 
        valuesWritten << " temperature values written" << endl;
    
    return bufferIndex;
}

// Read implementation (No-Op for this class)
void FluidTemperature::read(double* dataBuffer, const unsigned int dim)
{
    // OpenFOAM writes FluidTemperature, does not read it.
    Info << "FP DEBUG: FluidTemperature::read called (NO-OP - this is write-only)" << endl;
}

// Support check for location type
bool FluidTemperature::isLocationTypeSupported(const bool meshConnectivity) const
{
    // For temperature, we support volume coupling and face center coupling
    bool isSupported = (locationType_ == LocationType::volumeCenters || 
                       locationType_ == LocationType::faceCenters);
    
    if (!isSupported) {
        Info << "FP DEBUG: Location type check failed. Current type is not supported." << endl;
    }
    
    return isSupported;
}

// Get the data name
std::string FluidTemperature::getDataName() const
{
    return "T"; // Match what's used in FP::addWriters
}

} // namespace FP
} // namespace preciceAdapter
