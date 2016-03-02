#include "copy_polyhedron.h"

template <class Poly_B, class Poly_A>
void poly_copy(Poly_B& poly_b, const Poly_A& poly_a)
{
        poly_b.clear();
        Copy_polyhedron_to<Poly_A, Poly_B> modifier(poly_a);
        poly_b.delegate(modifier);
} 
