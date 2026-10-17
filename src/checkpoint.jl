export save_checkpoint, load_checkpoint, isidentical

"""
    save_checkpoint(filename, model)
Serialize the full `model` (agents, rng, clock, properties) to `filename`.
Stepping the loaded model reproduces the original run bit-for-bit,
provided the same Julia and package versions are used to read it back.
"""
function save_checkpoint(filename, model)
    mkpath(dirname(filename))
    serialize(filename, model)
    return filename
end

"""
    load_checkpoint(filename)
Load a model saved with `save_checkpoint`.
"""
load_checkpoint(filename) = deserialize(filename)

"""
    isidentical(a, b)
Bitwise recursive comparison of two objects (e.g. two models),
field by field; floats are compared with `isequal`.
"""
isidentical(a, b) = a === b || (typeof(a) == typeof(b) && _isidentical(a, b))
function _isidentical(a, b)
    if a isa Union{Number,AbstractString,Symbol}
        isequal(a, b)
    elseif a isa AbstractArray
        size(a) == size(b) && all(eachindex(a)) do i
            isassigned(a, i) == isassigned(b, i) &&
                (!isassigned(a, i) || isidentical(a[i], b[i]))
        end
    elseif a isa AbstractDict
        issetequal(keys(a), keys(b)) && all(isidentical(a[k], b[k]) for k in keys(a))
    elseif fieldcount(typeof(a)) == 0
        isequal(a, b)
    else
        all(isidentical(getfield(a, k), getfield(b, k)) for k in 1:fieldcount(typeof(a)))
    end
end
