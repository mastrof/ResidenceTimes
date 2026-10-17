export HaltUntilTumble

#==
Surface-contact rule for bacteria–phytoplankton encounters.
Collision detection itself is delegated to MicrobeAgents' `Obstacles` interface;
only the response to a contact is defined here.
==#

"""
    HaltUntilTumble()

Contact rule: the microbe is pinned at the surface (rest of the step discarded)
and its speed is set to zero, heading preserved. It resumes motion at its next
motile-state switch, since `reorient_step!` resamples the speed (`update_speed!`).

Unlike MicrobeAgents' `Halt`, the microbe cannot leave the surface by rotational
diffusion alone.
"""
struct HaltUntilTumble end

function MicrobeAgents.contact!(::HaltUntilTumble, microbe, c::SurfaceContact, model)
    microbe.speed = zero(microbe.speed)
    return zero(c.displacement)
end
