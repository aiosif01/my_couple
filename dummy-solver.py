#!/usr/bin/env python3

import precice
import numpy as np
import time

# Configure preCICE
participant_name = "dummy"
config_file_name = "cavity/precice-config.xml"
solver_process_index = 0
solver_process_size = 1

# Create preCICE participant
participant = precice.Participant(participant_name, config_file_name, solver_process_index, solver_process_size)

# Get mesh ID
mesh_name = "Cavity-Mesh"
mesh_dimensions = participant.get_mesh_dimensions(mesh_name)

# Define the number of vertices you expect from the cavity mesh
# This is just a placeholder - the actual vertices will come from cavity
# For testing, we won't worry about the exact number

# Initialize the participant
print("Initializing dummy solver")
participant.initialize()

# Coupling time loop
while participant.is_coupling_ongoing():
    # Get max time step size from preCICE
    dt = participant.get_max_time_step_size()
    
    # Read data (if available)
    if participant.requires_reading_checkpoint():
        print("Dummy: Reading checkpoint")
    
    # Read force data from cavity
    vertex_ids = participant.get_mesh_vertex_ids_from_positions(mesh_name, [])
    
    if len(vertex_ids) > 0:
        force = np.zeros((len(vertex_ids), mesh_dimensions))
        participant.read_data(mesh_name, "Force", vertex_ids, dt, force)
        print(f"Dummy: Read forces from Cavity: {force.shape}, max={np.max(force)}, min={np.min(force)}")
        
        # Simple processing - just return some stress based on the force
        # In a real solver, this would be more sophisticated
        stress = np.zeros_like(force)
        for i in range(len(force)):
            stress[i] = force[i] * 0.5  # Simple dummy relation
        
        # Write the stress data to be sent to cavity
        participant.write_data(mesh_name, "Stress", vertex_ids, stress)
        print(f"Dummy: Wrote stress to Cavity: {stress.shape}")
    else:
        print("Dummy: No vertices found in mesh")
    
    # Advance the coupling
    print(f"Dummy: Advancing with dt={dt}")
    participant.advance(dt)
    
    if participant.requires_writing_checkpoint():
        print("Dummy: Writing checkpoint")

# Finalize preCICE
participant.finalize()
print("Dummy solver completed")