#include "fpotencia_libs.h"

namespace fPotencia {
#ifndef COMMON_F
#define	COMMON_F

    /*
     * Performs the Kron reduction of an impedance matrix
     * Zabc: Primitive matrix
     * p: number of poles (that are not neutral)
     * m: number of neutrals
     * 
     * p+m=number of conductors
     */
    cx_mat Kron_cx(cx_mat Zprimitive, uint p, uint m);

    /*
     * Performs the Kron reduction of an impedance matrix
     * Zabc: Primitive matrix
     * p: number of poles (that are not neutral)
     * m: number of neutrals
     * 
     * p+m=number of conductors
     */
    mat Kron(mat Zprimitive, uint p, uint m);

    /*
     * Converts the ABC impedances to sequence (012) Impedances    
     */
    cx_mat abc2seq(cx_mat Z);

    /*
     * This function adds a 3x3 matrix M to the three phase admittance matrix Y
     * at the buses location i, j
     */
    void addMatToY(cx_mat &Y, cx_mat3 M, int i, int j);
    
    
#endif
}