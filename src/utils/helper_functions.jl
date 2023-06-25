export vec_from_span

# Helper function to add filler points to start/end span
function vec_from_span(span; num_points=100)
    if length(span) == 1
        vec_out = collect(range( 0, span[1], num_points))

    elseif length(span) == 2
        vec_out = collect(range( span[1], span[2], num_points))
        
    else
        return span
    end
end