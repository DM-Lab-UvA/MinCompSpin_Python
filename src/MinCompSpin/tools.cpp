#include <string>
#include <map>
#include <algorithm> // for std::reverse()

/******************************************************************************/
/***************************   Constant variables   ***************************/
/******************************************************************************/
const __uint32_t un32 = 1;  // we only need the first bit here -->> will be used to check the lowest bit of other unsigned integer (with arbitrarily larger number of bits)
//const unsigned int n_max = 128;  // for bitset

/******************************************************************************/
/*******************   Convert Integer to Binary string   *********************/
/******************************************************************************/
std::string int_to_bstring(__int128_t bool_nb, unsigned int n)
{
    std::string s;
    do
    {
        s.push_back( ((bool_nb & un32)?'1':'0') );
    } while(bool_nb >>= 1);

    reverse(s.begin(), s.end());
    s = (std::string(n - s.length(), '0')).append(s);

    return s;
}

// Using bitset: Very bad performance (~ double the time)
/*std::string int_to_bstring_Bitset(__int128_t bool_nb)     // Very bad performance (~ double the time)
{
    std::bitset<n_max> hi{ static_cast<unsigned long long>(bool_nb >> 64) },
            lo{ static_cast<unsigned long long>(bool_nb) },
            bits{ (hi << 64) | lo };
    return bits.to_string();
}*/

/******************************************************************************/
/****************   Count number of set bits of an integer  *******************/
/******************************************************************************/
unsigned int Bitset_count(__int128_t bool_nb)      // Using trick; by far the fastest option in all cases (faster than using bitset library or checking bits one by one)
{
    unsigned int count;
    for (count=0; bool_nb; count++)
        bool_nb &= bool_nb - 1;
    return count;
}

/******************************************************************************/
/*****************************    Check MCM   *********************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.

std::pair<bool, unsigned int> check_partition(std::map<unsigned int, __int128_t> Partition)
{
  std::map<unsigned int, __int128_t>::iterator Part;
  __int128_t sum = 0;
  unsigned int rank = 0; 

  for (Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    sum |= (*Part).second;
    rank += Bitset_count((*Part).second); 
  }

  return std::make_pair((Bitset_count(sum) == rank), rank); 
}


