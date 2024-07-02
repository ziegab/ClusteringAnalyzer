#ifndef TilesConstantsEB_h
#define TilesConstantsEB_h

namespace util {
static constexpr int32_t ceil(float num) {
  return (static_cast<float>(static_cast<int32_t>(num)) == num)
             ? static_cast<int32_t>(num)
             : static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
}
}  // namespace util

struct TilesConstantsEB {
  static constexpr int tileSize =5;
  static constexpr int minDim1 = -85;
  static constexpr int maxDim1 = 85;
  static constexpr int minDim2 = 0;
  static constexpr int maxDim2 = 360;
  static constexpr int nColumns = util::ceil((maxDim1 - minDim1) / tileSize);
  static constexpr int nRows = util::ceil((maxDim2 - minDim2) / tileSize);
  static constexpr float invDim1BinSize = nColumns / (maxDim1 - minDim1);
  static constexpr float invDim2BinSize = nRows / (maxDim2 - minDim2);
  static constexpr int nTiles = nColumns * nRows;
  //static constexpr float showerSigma = 0.05f; // in unit of xtals
};


#endif  // TilesConstantsEB_h