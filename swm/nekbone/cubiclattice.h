#ifndef CUBIC_LATTICE_H
#define CUBIC_LATTICE_H

#include "triplet.h"

namespace CUBICLAT
{
    Triplet id2tripletZ(Triplet in_partition, Idz in_idz);
    Idz tripletZ2id(Triplet in_partition, Triplet in_locationZ);
}//namespace CUBICLAT
#endif //CUBIC_LATTICE_H
