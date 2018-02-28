#ifndef TRIPLET_H
#define TRIPLET_H
 
typedef long Coord_t;
typedef long Idz;  //An accumulation of Coord_t is an Idz.                                                                                                                                          
 
class Triplet  //Keep it POD                                                                                                                                                                         
{
    public:
            Idz a,b,c;
             
                bool isinLattice(Triplet in_lattice)
                            {
                                            if(a<0 || a>=in_lattice.a) return false;
                                                        if(b<0 || b>=in_lattice.b) return false;
                                                                    if(c<0 || c>=in_lattice.c) return false;
                                                                                return true;
                                                                                        }
                 
                    inline Idz abc() const { return (a*b*c); }
                        inline Idz ab()  const { return (a*b); }
                         
                            void operator+=(Triplet in)        { a+=in.a; b+=in.b; c+=in.c; }
                                void operator*=(Triplet in_scale)  { a*=in_scale.a; b*=in_scale.b; c*=in_scale.c; }
                                 
                                    Triplet():a(0),b(0),c(0) {}
                                        Triplet(Coord_t in_a, Coord_t in_b, Coord_t in_c): a(in_a),b(in_b),c(in_c) {}
};
 
 
#endif // TRIPLET_H   
