#ifndef __BOOST_PTREE_ARRAY_TO_STD_VECTOR__H__
#define __BOOST_PTREE_ARRAY_TO_STD_VECTOR__H__
 
#include <assert.h>
#include <iostream>

template <typename T> std::vector<T> boost_ptree_array_to_std_vector(boost::property_tree::ptree const& pt, boost::property_tree::ptree::key_type const& key, std::vector<T> def, bool disallow_empty_arrays=true)
{
         std::vector<T> r;
         for (auto& item : pt.get_child(key))
         {
             cout << item.second.get_value<T>() << endl;
            r.push_back(item.second.get_value<T>());
         }
         if(r.empty()) {
                      std::cerr << "\nERROR: boost_ptree_array_to_std_vector matched a tag (" << key << ") that is not a valid array!\n" << std::endl;
                      assert(r.empty() == false);
                     }
               return r;
}
#endif
