#ifndef _NDINDEX_HPP
#define _NDINDEX_HPP

#include <assert.h>

template <int N>
class RowMajorIndexer {
	public:
	RowMajorIndexer()
    {
    }

	void index_to_tuple(const int shape[], int index, int tuple[]) const
    {
        assert(index >= 0);

        int tmp = index;

        int count = 1;
        for (int i=N-1; i>=0; i--) {
            tuple[i] = tmp % shape[i];
            tmp /= shape[i];
            count *= shape[i];
        }

        // Safety check
        assert(index < count);
    }

	void tuple_to_index(const int shape[], const int tuple[], int * index) const
    {
        int tmp = 0;
        for (int i=0; i<N; i++) {
            tmp *= shape[i];
            tmp += tuple[i];
        }
        *index = tmp;
    }

};


#endif
