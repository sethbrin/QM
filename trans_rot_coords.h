/**
 * some helper functions
 *
 */

#ifndef TRANS_ROT_COORDS_H_
#define TRANS_ROT_COORDS_H_

#include <vector>
#include <assert>

namespace math
{
  template <typename T>
  class Coord
  {
  public:
    Coord():
      values_(std::vector<T>()),
      nDim_(0)
    {}

    Coord(vector<T> values):
      values_(values),
      nDim_(values.size())
    {}

    vector<T>& getCoords()
    {
      return values_;
    }

    T getCoords(int dim)
    {
      assert(dim < nDim_);

      return values_[dim];
    }

    int getDimension()
    {
      return nDim_;
    }

  private:
    std::vector<T> values_;
    int nDim_;
  }

  class helper
  {
  public:
    static Coord<T>& translate(const Coord<T>& vec, const Coord<T>& dVec)
    {
      assert(vec.getDimension() == dVec.getDimension());
      Coord<T> res(vec);

      for (int i=0; i<dVec.getDimension(); i++)
      {
        res[i] += dVec[i];
      }

      return res;
    }
  }
}
#endif
