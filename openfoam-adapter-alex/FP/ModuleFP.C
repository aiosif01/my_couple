// FP/ModuleFP.C - Central include file for the FP module

#include "FP.C"                // Includes the FluidParticle class implementation
#include "FP.H"                // Include the header first for proper dependencies
#include "FluidTemperature.C"  // Temperature writer implementation
#include "FluidTemperature.H"  // Include header first
#include "ParticlePosition.C"  // Position reader implementation
#include "ParticlePosition.H"  // Include header first
// Add #include for any other .C files in the FP/ directory

// The include order matters - headers should come before implementations
// to ensure proper symbol resolution