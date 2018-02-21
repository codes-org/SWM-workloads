#include "cubiclattice.h"

Triplet
CUBICLAT::id2tripletZ(Triplet in_partition, Idz in_idz)
    {
        //in_partition : the full count of blocks in the (a|b|c)-axis direction
        //in_id :  an offset of the current block in the lattice
        //          0 <= in_id < in_partition.abc()
        // Returns the tuplet (x,y,z) of the block specified by in_id
        //      0<= x < in_partition.a
        //      0<= y < in_partition.b
        //      0<= z < in_partition.c
        // The returned triplet will in the Z-ordering.
        Triplet t;
        const Idz ab = in_partition.ab();
        t.c = in_idz / ab;
        in_idz -= t.c * ab;
        t.b = in_idz / in_partition.a;
        t.a = in_idz - t.b * in_partition.a;
        return t;
    }

Idz
CUBICLAT::tripletZ2id(Triplet in_partition, Triplet in_locationZ)
    {
        //in_partition : the full count of blocks in the (a|b|c)-axis direction
        //in_location : the location (x,y,z) of the current block within the lattice
        //              ordered in the Z-ordering.
        // Returns the id associate with in_location
        Idz idz =0;
        idz = in_locationZ.a + in_locationZ.b * in_partition.a + in_locationZ.c * in_partition.ab();
        return idz;
    }
